"""Build a donor effect table from X-Atlas/Orion (Xaira).

Orion is genome-wide FiCS Perturb-seq in HCT116 (colorectal) and HEK293T (kidney
epithelial): 8M cells, all 18,903 protein-coding genes targeted, ~16k UMI/cell. It matters
because Replogle K562-gwps is our only donor covering the VCC panel (272/300), so every
prediction currently leans on one cell line -- and Nadig 2025 measured cross-line effect
correlations from 0.40 (K562/RPE1) to 0.75 (HepG2/RPE1) depending on tissue similarity,
so a single donor is a single point of failure. Orion also carries 25 of the 28 panel
targets K562-gwps misses.

Storage format is one row per cell with parallel sparse lists (`gene_token_id`,
`gene_expression`), so group sums accumulate directly from the lists -- no matrix is ever
built, and memory stays flat across all 109 batches.

    uv run python scripts/build_orion_donor.py --line HCT116
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp

TS = 5e4
CONTROL_PATTERNS = (r"^non.?targeting", r"^ntc", r"^control", r"^safe.?harbou?r", r"^neg")


def is_control(name: str) -> bool:
    n = str(name).lower()
    return any(re.match(p, n) for p in CONTROL_PATTERNS)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--line", default="HCT116")
    ap.add_argument("--data", default="data/external/orion/data")
    ap.add_argument("--meta", default="data/external/orion/metadata/gene_metadata.parquet")
    ap.add_argument("--out", default=None)
    ap.add_argument("--min-cells", type=int, default=20)
    args = ap.parse_args()

    out = Path(args.out or f"data/donors/orion_{args.line.lower()}")
    gm = pd.read_parquet(args.meta).sort_values("gene_token_id")
    n_genes = int(gm["gene_token_id"].max()) + 1
    genes = np.array([""] * n_genes, dtype=object)
    genes[gm["gene_token_id"].to_numpy()] = gm["gene_name"].to_numpy()
    print(f"gene axis: {n_genes:,} tokens, {int((genes != '').sum()):,} named")

    files = sorted(Path(args.data).glob(f"{args.line}_Batch*.parquet"),
                   key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))
    if not files:
        raise SystemExit(f"no {args.line} parquet files in {args.data}")
    print(f"{len(files)} batches\n")

    # Accumulate per-target sums batch by batch. Each batch becomes one CSR matrix built
    # straight from the sparse lists, then a sparse indicator matmul does the group sum --
    # np.add.at over ~10^10 (cell, gene) entries would take hours.
    sums: dict[str, np.ndarray] = {}
    counts: dict[str, int] = defaultdict(int)
    ctrl_sum = np.zeros(n_genes, dtype=np.float64)
    ctrl_n = 0

    for i, f in enumerate(files, 1):
        t = pq.read_table(f, columns=["gene_token_id", "gene_expression",
                                      "gene_target", "pass_guide_filter"]).to_pandas()
        t = t[t["pass_guide_filter"] == 1]
        if t.empty:
            continue
        lengths = t["gene_token_id"].map(len).to_numpy()
        indptr = np.zeros(len(t) + 1, dtype=np.int64)
        np.cumsum(lengths, out=indptr[1:])
        idx = np.concatenate(t["gene_token_id"].to_numpy()).astype(np.int32)
        val = np.concatenate(t["gene_expression"].to_numpy()).astype(np.float64)
        X = sp.csr_matrix((val, idx, indptr), shape=(len(t), n_genes))

        tgt = t["gene_target"].to_numpy()
        ctrl_mask = np.fromiter((is_control(x) for x in tgt), bool, len(tgt))
        if ctrl_mask.any():
            ctrl_sum += np.asarray(X[ctrl_mask].sum(0)).ravel()
            ctrl_n += int(ctrl_mask.sum())

        pert = tgt[~ctrl_mask]
        if len(pert):
            Xp = X[~ctrl_mask]
            uniq, codes = np.unique(pert, return_inverse=True)
            ind = sp.csr_matrix(
                (np.ones(len(codes), np.float64), (codes, np.arange(len(codes)))),
                shape=(len(uniq), Xp.shape[0]))
            block = np.asarray((ind @ Xp).todense())
            for k, name in enumerate(uniq):
                if name in sums:
                    sums[name] += block[k]
                else:
                    sums[name] = block[k].copy()
            for name, c in zip(*np.unique(pert, return_counts=True)):
                counts[name] += int(c)
        if i % 10 == 0 or i == len(files):
            print(f"  {i}/{len(files)} batches | {len(sums):,} targets | "
                  f"{ctrl_n:,} control cells", flush=True)

    if ctrl_n == 0:
        raise SystemExit(f"no control cells found; gene_target examples: "
                         f"{list(sums)[:10]}")

    b_ctrl = np.log1p(TS * ctrl_sum / ctrl_sum.sum()).astype(np.float32)
    names = [g for g, c in counts.items() if c >= args.min_cells and sums[g].sum() > 0]
    names.sort()
    delta = np.vstack([
        np.log1p(TS * sums[g] / sums[g].sum()).astype(np.float32) - b_ctrl for g in names
    ])

    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "delta.npy", delta)
    np.save(out / "b_ctrl.npy", b_ctrl)
    (out / "perts.txt").write_text("\n".join(names))
    (out / "genes.txt").write_text("\n".join(map(str, genes)))
    cells = [counts[g] for g in names]
    norms = np.linalg.norm(delta, axis=1)
    (out / "meta.json").write_text(json.dumps(dict(
        source=f"X-Atlas/Orion {args.line}", target_sum=TS,
        n_perturbations=len(names), n_genes=n_genes, n_control_cells=ctrl_n,
        min_cells=args.min_cells, dropped=len(counts) - len(names),
        cells_per_pert=dict(min=int(min(cells)), median=int(np.median(cells)),
                            max=int(max(cells))),
        effect_norm=dict(median=float(np.median(norms)), min=float(norms.min()),
                         max=float(norms.max())),
    ), indent=2))
    print(f"\nwrote {out}: delta {delta.shape}, {ctrl_n:,} control cells")
    print(f"  cells/pert median {np.median(cells):.0f}, effect norm median {np.median(norms):.2f}")

    vcc = set(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
    hit = vcc & set(names)
    print(f"  VCC panel coverage: {len(hit)}/300 = {100*len(hit)/300:.0f}%")


if __name__ == "__main__":
    main()
