"""Build a count-level scoring harness from X-Atlas/Orion HCT116.

## Why this is the harness the project has been missing

Everything with a count-level number attached to it so far was measured on `h1-fast`: 40
perturbations of **H1 hESC**, scored against the **2025** panel. Both are wrong for this task
in ways already known to matter:

- **Wrong biology.** H1 is embryonic; the 2026 contexts are somatic lines from different
  tissues. A +20% transfer cosine on H1 produced a *negative* `pds` move on A/B/C (exp 010).
- **Wrong gene population.** The 2025 panel is common-essential enriched; the 2026 panel was
  built to *exclude* common-essential genes (0/2,058 overlap with Replogle K562-essential).
  So `h1-fast` measures the DE metrics on exactly the gene class the competition removed.

Orion HCT116 fixes both. It is a **somatic** line (colorectal), it carries **290 of the 300
VCC panel targets**, and the cell-level counts are on disk. Holding it out and predicting from
K562-gwps + Orion-HEK293T is the same somatic-to-somatic, VCC-panel transfer the competition
asks for -- and unlike LODO it produces all six metrics from real emitted counts.

## Design

- Perturbations: VCC-panel targets present in **all three** donor lines, so the prediction is
  defined the same way it is in production, keeping the `--n-perts` with the most cells.
- 400 cells per perturbation and 18,400 control cells, matching the competition exactly.
- Gene axis: the official 18,533 in `gene_names.csv` order. Orion genes outside it are
  dropped, genes it lacks are zero.
- No depth thinning: Orion runs at ~16k UMI/cell against the reference's 20k, so thinning
  could only make the gap worse. The realised median is recorded in the metadata.

## The sampling bug this script shipped until 2026-08-27, and why it mattered

The original loop filled both pools **greedily, first-come-first-served**: controls until
`--n-control` was reached, perturbation cells until each hit `--cells-per-pert`, then broke
out of the batch loop early. Controls are ~4.85% of every batch, so 18,400 of them were
satisfied inside the first ~20 batches, while perturbation cells -- 1-3 per batch per target
-- had to accumulate across all 109. Batch depth varies enormously (median non-targeting UMI
24,850 in batches 1-20 against 16,589 in 21-60 and 18,007 in 61-109), so the two pools ended
up drawn from different parts of the run.

The result was a harness whose **control cells are 38% deeper than its perturbation cells**
(median UMI 16,798 against 12,178, ratio 0.725) in **132/132 perturbations**. The Wilcoxon is
a rank test and CPM normalisation does not fix the per-cell zero fraction, so uniformly
shallower perturbed cells read as uniformly *down*: the reference's DE calls came out **90.8%
negative**, against ~0.50 on h1-test, Replogle K562 and RPE1. In the raw parquets the pools
are balanced within batch at ratio 0.9896, so the builder manufactured the whole effect.

That harness is the instrument experiments 031, 061, 065, 090, 091, 092 and 093 were tuned
on, and it is why experiment 093 measured `+0.0202` avg6 and the board returned `-0.00291`.
Full diagnosis: `docs/HARNESS-DEPTH-CONFOUND.md`.

`--sampling reservoir` (the default) fixes it: every pool is filled by uniform reservoir
sampling over **all** batches, so both are drawn from the same batch distribution and the
depth mismatch cannot arise. `--sampling greedy` reproduces the original file bit-for-bit
for anyone who needs to re-read an old result.

    uv run python scripts/build_orion_harness.py --line HCT116 --n-perts 80
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp

CONTROL_PATTERNS = (r"^non.?targeting", r"^ntc", r"^control", r"^safe.?harbou?r", r"^neg")
CONTROL_LABEL = "non-targeting"
PERT_COL = "target_gene"


def is_control(name: str) -> bool:
    n = str(name).lower()
    return any(re.match(p, n) for p in CONTROL_PATTERNS)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--line", default="HCT116")
    ap.add_argument("--data", default="data/external/orion/data")
    ap.add_argument("--meta", default="data/external/orion/metadata/gene_metadata.parquet")
    ap.add_argument("--donors", nargs="+",
                    default=["data/donors/k562_gwps", "data/donors/orion_hek293t"])
    ap.add_argument("--n-perts", type=int, default=80)
    ap.add_argument("--cells-per-pert", type=int, default=400)
    ap.add_argument("--n-control", type=int, default=18400)
    ap.add_argument("--out", default="data/harness/orion-hct116.h5ad")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--thin-controls", action="store_true",
        help=(
            "Binomially thin control counts so their depth distribution matches the "
            "perturbation pool's, quantile for quantile. Use with --sampling greedy: it "
            "removes the depth confound without the power loss that uniform reservoir "
            "sampling costs on Orion's heterogeneous batches."
        ),
    )
    ap.add_argument(
        "--sampling", default="reservoir", choices=["reservoir", "greedy"],
        help=(
            "How both pools are filled. 'reservoir' samples uniformly over every batch, so "
            "controls and perturbation cells come from the same batch distribution and "
            "cannot end up at different sequencing depths. 'greedy' is the original "
            "first-come-first-served fill with the early break, which produced a 0.725 "
            "pert/control depth ratio and a 90.8 percent negative reference DE set; keep it only "
            "to reproduce the old file."
        ),
    )
    args = ap.parse_args()

    spec_genes = list(pd.read_csv("data/challenge/gene_names.csv")["gene_name"])
    tindex = {g: i for i, g in enumerate(spec_genes)}
    gm = pd.read_parquet(args.meta).sort_values("gene_token_id")
    n_tok = int(gm["gene_token_id"].max()) + 1
    tok_name = np.array([""] * n_tok, dtype=object)
    tok_name[gm["gene_token_id"].to_numpy()] = gm["gene_name"].to_numpy()
    # token -> official column, -1 where the symbol is not on the official axis
    tok2col = np.full(n_tok, -1, dtype=np.int64)
    for t, nm in enumerate(tok_name):
        c = tindex.get(nm)
        if c is not None:
            tok2col[t] = c
    print(f"gene axis: {n_tok:,} Orion tokens -> {int((tok2col >= 0).sum()):,} of "
          f"{len(spec_genes):,} official genes")

    # Candidate perturbations: on the VCC panel and carried by every donor line, so the
    # harness asks the production question rather than an easier one.
    panel = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
    have = []
    for d in args.donors:
        have.append(set((Path(d) / "perts.txt").read_text().split("\n")))
    candidates = {g for g in panel if all(g in h for h in have)}
    print(f"{len(candidates)} panel targets carried by every donor")

    files = sorted(Path(args.data).glob(f"{args.line}_Batch*.parquet"),
                   key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))
    print(f"{len(files)} batches\n")

    rng = np.random.default_rng(args.seed)
    kept: dict[str, list] = {}
    n_kept: dict[str, int] = {}
    ctrl: list = []
    n_ctrl = 0
    # Reservoir bookkeeping: how many cells of each class we have *seen*, which is what
    # makes the retained sample uniform over the whole run rather than over its prefix.
    seen_ctrl = 0
    seen_pert: dict[str, int] = {}
    reservoir = args.sampling == "reservoir"

    def to_official(idx_list, val_list):
        idx = np.asarray(idx_list, dtype=np.int64)
        col = tok2col[idx]
        m = col >= 0
        return col[m].astype(np.int32), np.asarray(val_list, dtype=np.float32)[m]

    for i, f in enumerate(files, 1):
        t = pq.read_table(f, columns=["gene_token_id", "gene_expression", "gene_target",
                                      "pass_guide_filter"]).to_pandas()
        t = t[t["pass_guide_filter"] == 1]
        if t.empty:
            continue
        tgt = t["gene_target"].to_numpy()
        gi = t["gene_token_id"].to_numpy()
        gv = t["gene_expression"].to_numpy()
        for j in range(len(t)):
            name = tgt[j]
            if is_control(name):
                seen_ctrl += 1
                if n_ctrl < args.n_control:
                    ctrl.append(to_official(gi[j], gv[j]))
                    n_ctrl += 1
                elif reservoir:
                    k = int(rng.integers(0, seen_ctrl))
                    if k < args.n_control:
                        ctrl[k] = to_official(gi[j], gv[j])
                continue
            if name not in candidates:
                continue
            seen_pert[name] = seen_pert.get(name, 0) + 1
            buf = kept.setdefault(name, [])
            if len(buf) < args.cells_per_pert:
                buf.append(to_official(gi[j], gv[j]))
                n_kept[name] = len(buf)
            elif reservoir:
                k = int(rng.integers(0, seen_pert[name]))
                if k < args.cells_per_pert:
                    buf[k] = to_official(gi[j], gv[j])
        full = sum(1 for v in n_kept.values() if v >= args.cells_per_pert)
        if i % 20 == 0 or i == len(files):
            print(f"  {i}/{len(files)} | {full} perturbations full | "
                  f"{n_ctrl:,}/{args.n_control:,} controls", flush=True)
        # A uniform reservoir has to see every batch; stopping early is exactly the bug.
        if not reservoir and full >= args.n_perts and n_ctrl >= args.n_control:
            print("  enough cells collected, stopping early", flush=True)
            break

    # Rank by how many cells the run actually offered, not by how many we retained --
    # under reservoir sampling every full perturbation retains exactly `cells_per_pert`,
    # so `n_kept` is a constant and cannot order anything.
    names = sorted([g for g, v in n_kept.items() if v >= args.cells_per_pert],
                   key=lambda g: (-seen_pert.get(g, 0), g))[:args.n_perts]
    if not names or n_ctrl < args.n_control:
        raise SystemExit(f"insufficient: {len(names)} full perturbations, {n_ctrl} controls")
    print(f"\nkeeping {len(names)} perturbations x {args.cells_per_pert} cells "
          f"+ {n_ctrl:,} controls")

    rows, obs = [], []
    for g in names:
        rows.extend(kept[g][:args.cells_per_pert])
        obs.extend([g] * args.cells_per_pert)
    rows.extend(ctrl)

    if args.thin_controls:
        # Quantile-match the control depth distribution to the perturbation pool's by
        # binomial thinning of raw counts.
        #
        # `--sampling reservoir` removes the depth gap by drawing both pools from the same
        # batches, but Orion's 109 batches differ enough that uniformly-drawn controls are
        # heterogeneous, which widens the reference distribution and costs power: measured
        # n_real per perturbation falls from 326 to 12. Thinning takes the other route --
        # keep the tight, single-region control pool the greedy fill produces and remove
        # only the depth difference, which is the thing the rank test actually reacts to.
        #
        # Binomial thinning is the right operator because it is exactly how a shallower
        # sequencing run would have produced these cells: each UMI survives independently.
        n_pert_rows = len(names) * args.cells_per_pert
        ptot = np.array([v.sum() for _, v in rows[:n_pert_rows]], dtype=np.float64)
        ctot = np.array([v.sum() for _, v in rows[n_pert_rows:]], dtype=np.float64)
        q = (np.argsort(np.argsort(ctot)) + 0.5) / ctot.size
        target = np.quantile(ptot, q)
        keep_p = np.clip(target / np.maximum(ctot, 1.0), 0.0, 1.0)
        thinned = 0
        for i, pr in enumerate(keep_p):
            if pr >= 1.0:
                continue
            c, v = rows[n_pert_rows + i]
            nv = rng.binomial(v.astype(np.int64), pr).astype(np.float32)
            m = nv > 0
            rows[n_pert_rows + i] = (c[m], nv[m])
            thinned += 1
        after = np.array([v.sum() for _, v in rows[n_pert_rows:]], dtype=np.float64)
        print(f"thinned {thinned:,}/{ctot.size:,} control cells: "
              f"median UMI {np.median(ctot):,.0f} -> {np.median(after):,.0f} "
              f"against perturbation median {np.median(ptot):,.0f}")
    obs.extend([CONTROL_LABEL] * n_ctrl)

    indptr = np.zeros(len(rows) + 1, dtype=np.int64)
    np.cumsum([len(r[0]) for r in rows], out=indptr[1:])
    X = sp.csr_matrix((np.concatenate([r[1] for r in rows]),
                       np.concatenate([r[0] for r in rows]), indptr),
                      shape=(len(rows), len(spec_genes)))
    libs = np.asarray(X.sum(1)).ravel()
    print(f"matrix {X.shape}, {X.nnz:,} nonzeros, median library {np.median(libs):,.0f} UMI")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=X, obs=pd.DataFrame({PERT_COL: obs},
                                     index=[f"c{i}" for i in range(len(obs))]),
               var=pd.DataFrame(index=spec_genes)).write_h5ad(out, compression="lzf")
    Path(str(out) + ".meta.json").write_text(json.dumps(dict(
        source=f"X-Atlas/Orion {args.line}", n_perturbations=len(names),
        cells_per_pert=args.cells_per_pert, n_control=n_ctrl,
        median_umi=float(np.median(libs)), nnz=int(X.nnz),
        genes=len(spec_genes), perts=names), indent=2))
    print(f"wrote {out} ({out.stat().st_size/1e9:.2f} GB)")


if __name__ == "__main__":
    main()
