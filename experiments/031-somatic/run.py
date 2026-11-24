"""031 — everything, scored at count level on a SOMATIC, VCC-panel reference.

Every count-level number in this project so far came from `h1-fast`: 40 perturbations of H1
hESC against the 2025 panel. Two things are wrong with that and both are documented failures:
H1 is embryonic where the 2026 contexts are somatic (a +20% transfer cosine on H1 produced a
*negative* `pds` move on A/B/C), and the 2025 panel is common-essential enriched — exactly the
gene class the 2026 panel was built to exclude.

`data/harness/orion-hct116.h5ad` is 132 real VCC-panel targets in HCT116 (colorectal, somatic)
with 18,400 control cells. **HCT116 is one of our donors, so it is dropped from the donor set
here** — predictions come from K562-gwps + Orion-HEK293T only. That is leave-one-donor-out at
count level, on the right biology and the right genes.

Scored on **raw** metric values read straight from `agg_results.csv` and rescaled with the
official 2026 anchors, per the translation rule in the log — no local replicate bundle needed,
which is what made `h1-test` a six-hour dead end.

Three open questions, one grid:

  gamma    does whitening's +0.070 LODO gain survive on somatic VCC-panel data?
  q        magnitude sharpening (exp 030): large `pds` gain, `nmae` cost estimated only by a
           proxy whose alpha calibration the count-level run contradicts
  spread   the per-gene variance channel (exp 027): +0.028 scaled `reach` on h1-fast
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, emit_counts, TS   # noqa: E402
from vcc2026.multi_donor import MultiDonorModel              # noqa: E402
from vcc2026.whiten import Whitener                          # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"
# Official 2026 anchors, midpoints of the A/B/C ranges (AGENTS.md section 4).
ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}


def official(raw: dict) -> dict:
    out = {}
    for k, (b, r) in ANCHORS.items():
        if k in raw:
            v = (raw[k] - b) / (r - b)
            out[k] = max(0.0, min(1.0, v)) if k.startswith("expr_mse") else v
    out["avg6"] = float(np.mean(list(out.values())))
    return out


def confidence(model, perts, weights):
    G = len(model.target_genes)
    out = np.zeros((len(perts), G))
    for i, p in enumerate(perts):
        have = [n for n in model.names if p in model.pert_index[n]]
        if not have:
            continue
        V = np.zeros((len(have), G))
        for j, n in enumerate(have):
            V[j, model.tcols[n]] = model.s[n][model.pert_index[n][p]][model.dcols[n]]
        agree = np.all(np.sign(V) == np.sign(V[0]), axis=0) & (V[0] != 0)
        out[i] = np.exp(np.log(np.abs(V) + 1e-12).mean(0)) * agree
    return out


def build(real, model, wcols, W, alpha, gamma, q, spread, tight_k, min_cpm, out,
          n_emit=400, seed=0, basal_beta=0.0):
    a = ad.read_h5ad(real)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    labels = a.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(a.var_names)
    cm = labels == CONTROL
    Xc = X[cm].tocsr()
    Pc = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())
    libs = np.asarray(Xc.sum(1)).ravel()
    perts = sorted(set(labels) - {CONTROL})
    w = {n: 1.0 for n in model.names}

    E = np.vstack([model.effect(p, 1.0, 0.0, weights=w) for p in perts]).astype(np.float64)
    if gamma > 0:
        before_op = float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
        E[:, wcols] = W.apply(E[:, wcols], gamma, 0.0, renorm=False)
        E[:, wcols] *= before_op / (float(np.median(np.linalg.norm(E[:, wcols], axis=1)))
                                    or 1.0)
    if q != 1.0:
        med = float(np.median(np.linalg.norm(E, axis=1))) or 1.0
        E = np.sign(E) * np.abs(E) ** q
        E *= med / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
    # alpha is applied last so it always means "emitted magnitude", whatever precedes it
    E *= alpha * model.effect_scale / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)

    # tight_k = 0 with spread > 0 means a GLOBAL spread: every gene gets the same fraction
    # of Poisson noise. That is the arm that buys the `expr_mse` jackknife credit, which a
    # minimum-variance emitter forfeits entirely (Chat log, 2026-08-22).
    if basal_beta:
        # Per-gene modulation by the TARGET's expression level. Experiment 017 measured this
        # for `pds` only (+0.003 at beta=0.25, harmful by 1.0). The untested question is
        # `nmae`: it is scored over the reference's DE genes, which are predominantly the
        # well-expressed ones, and whitening systematically down-weights exactly those.
        e = np.expm1(b_ctrl)
        cpm_t = 1e6 * e / max(e.sum(), 1e-9)
        med = float(np.median(np.linalg.norm(E, axis=1))) or 1.0
        E = E * (np.maximum(cpm_t, 1e-3) ** basal_beta)[None, :]
        E *= med / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)

    C = confidence(model, perts, w) if (spread > 0 and tight_k) else None
    global_spread = spread if (spread > 0 and not tight_k) else 0.0
    rng = np.random.default_rng(seed)
    data, indices, indptr, obs = [], [], [0], []
    for i, p in enumerate(perts):
        prof = model.profile_from_effect(E[i], b_ctrl, min_cpm=min_cpm,
                                         gate_unexpressed=True)
        pick = rng.choice(len(libs), n_emit, replace=len(libs) < n_emit)
        if C is not None:
            k = min(tight_k, int((C[i] > 0).sum()))
            tight = np.zeros(len(genes), bool)
            if k:
                tight[np.argpartition(-C[i], k - 1)[:k]] = True
            d, idx, ip = emit_counts(prof, libs[pick], 0.0, rng, mode="quasi",
                                     spread=np.where(tight, 0.0, spread))
        elif global_spread > 0:
            d, idx, ip = emit_counts(prof, libs[pick], 0.0, rng, mode="quasi",
                                     spread=np.full(len(genes), global_spread))
        else:
            d, idx, ip = emit_counts(prof, libs[pick], 0.0, rng, mode="even")
        data.append(d)
        indices.append(idx)
        indptr.extend(list(ip[1:] + indptr[-1]))
        obs += [p] * n_emit
    data.append(Xc.data)
    indices.append(Xc.indices)
    indptr.extend(list(Xc.indptr[1:] + indptr[-1]))
    obs += [CONTROL] * Xc.shape[0]

    M = sp.csr_matrix((np.concatenate(data).astype(np.float32),
                       np.concatenate(indices).astype(np.int32),
                       np.asarray(indptr, dtype=np.int64)), shape=(len(obs), len(genes)))
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=M, obs=pd.DataFrame({PERT_COL: obs},
                                     index=[f"c{i}" for i in range(len(obs))]),
               var=pd.DataFrame(index=genes)).write_h5ad(out, compression="lzf")
    return out


def score_one(t):
    tag, pred, real, wd = t
    run = Path(wd) / "run"
    r = subprocess.run(["uv", "run", "cell-eval2", "run", "-ap", pred, "-ar", real,
                        "-o", str(run), "--preset", "vcc2026", "--pert-col", PERT_COL,
                        "--control", CONTROL], capture_output=True, text=True)
    agg = run / "agg_results.csv"
    if r.returncode != 0 or not agg.exists():
        return tag, {"error": ("run: " + (r.stderr or r.stdout))[-900:]}
    raw = pd.read_csv(agg).set_index("statistic")
    d = {m: float(raw.loc["mean", m]) for m in raw.columns if m != "statistic"}
    return tag, {"raw": d, "official": official(d)}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real", default="data/harness/orion-hct116.h5ad")
    ap.add_argument("--donors", nargs="+",
                    default=["k562_gwps=data/donors/k562_gwps",
                             "orion_hek293t=data/donors/orion_hek293t"])
    ap.add_argument("--outdir", default="experiments/031-somatic/out")
    ap.add_argument("--configs", default="experiments/031-somatic/grid.json")
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--whiten-axis", default="inter", choices=["inter", "union"],
                    help="genes the operator covers: every donor (inter) or any donor (union)")
    ap.add_argument(
        "--exclude-panel",
        default="auto",
        choices=["auto", "real", "challenge"],
        help=(
            "target-gene coordinates excluded before dense whitening. 'auto' uses "
            "the full challenge panel when the real targets are its subset, else "
            "the held-out real panel."
        ),
    )
    ap.add_argument(
        "--basis-prune",
        default="none",
        choices=["none", "essential", "random"],
        help=(
            "drop rows from the whitening basis before fitting. 'essential' "
            "drops perturbations naming a common-essential gene (experiment 090 "
            "showed adding such rows hurts monotonically, -0.0004 per 1000 rows); "
            "'random' drops the same COUNT at random, the size control that "
            "separates basis composition from basis size."
        ),
    )
    ap.add_argument("--essential-list", default="data/donors/k562_essential/perts.txt")
    ap.add_argument("--prune-seed", type=int, default=0)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    real_backed = ad.read_h5ad(args.real, backed="r")
    genes = np.asarray(real_backed.var_names)
    real_targets = sorted(
        set(real_backed.obs[PERT_COL].astype(str)) - {CONTROL}
    )
    real_backed.file.close()
    donors = {}
    for item in args.donors:
        name, path = item.split("=", 1)
        donors[name] = Donor.load(Path(path))
    model = MultiDonorModel(donors, genes, gene_norm=False, eb=False)
    print(
        f"donors {list(donors)} | held-out real={args.real} "
        f"targets={len(real_targets)}"
    )

    if args.whiten_axis == "union":
        shared = set()
        for n in model.names:
            shared |= set(model.target_genes[model.tcols[n]])
    else:
        shared = set(model.target_genes[model.tcols[model.names[0]]])
        for n in model.names[1:]:
            shared &= set(model.target_genes[model.tcols[n]])
    # The whitening matrix is dense: a target gene coordinate can be mixed
    # into every scored coordinate. The old H1 runs excluded the wrong 300
    # genes, leaving 80/85 H1 self-target coordinates in the operator. For the
    # HCT116 harness, however, its 132 targets are a subset of the 2026 panel
    # and production excludes all 300, so auto deliberately keeps that superset.
    challenge_targets = set(
        pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].astype(str)
    )
    real_target_set = set(real_targets)
    if args.exclude_panel == "challenge":
        excluded_targets = challenge_targets
    elif args.exclude_panel == "real":
        excluded_targets = real_target_set
    else:
        excluded_targets = (
            challenge_targets
            if real_target_set <= challenge_targets
            else real_target_set
        )
    shared -= excluded_targets
    ti = {g: i for i, g in enumerate(model.target_genes)}
    wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
    cache = outdir / f"whitener_{args.whiten_axis}.npz"
    cache_meta = cache.with_suffix(".json")
    expected_meta = {
        "axis": args.whiten_axis,
        "donors": list(model.names),
        "n_cols": int(len(wcols)),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(excluded_targets),
        "basis_prune": args.basis_prune,
        "prune_seed": args.prune_seed if args.basis_prune == "random" else None,
    }
    if (
        cache.exists()
        and cache_meta.exists()
        and json.loads(cache_meta.read_text()) == expected_meta
    ):
        W = Whitener.load(cache)
    else:
        essential = set()
        if args.basis_prune != "none":
            essential = {
                line.strip()
                for line in Path(args.essential_list).read_text().splitlines()
                if line.strip()
            }
        prune_rng = np.random.default_rng(args.prune_seed)
        rows = []
        dropped = 0
        for n in model.names:
            S = np.zeros((model.s[n].shape[0], len(model.target_genes)), np.float32)
            S[:, model.tcols[n]] = model.s[n][:, model.dcols[n]]
            block = model._renorm[n] * S[:, wcols].astype(np.float64)
            del S
            if args.basis_prune != "none":
                names = [None] * block.shape[0]
                for pert, idx in model.pert_index[n].items():
                    if 0 <= idx < len(names):
                        names[idx] = pert
                is_ess = np.array([nm in essential for nm in names], dtype=bool)
                if args.basis_prune == "random":
                    # same count, drawn uniformly: isolates composition from size
                    keep = np.ones(block.shape[0], dtype=bool)
                    k = int(is_ess.sum())
                    if k:
                        keep[prune_rng.choice(block.shape[0], k, replace=False)] = False
                else:
                    keep = ~is_ess
                dropped += int((~keep).sum())
                block = block[keep]
            rows.append(block)
        B = np.vstack(rows)
        del rows
        print(
            f"fitting whitener {B.shape} on {len(wcols)} genes "
            f"(prune={args.basis_prune}, dropped={dropped})",
            flush=True,
        )
        W = Whitener.fit(B)
        W.save(cache)
        cache_meta.write_text(json.dumps(expected_meta, indent=2))
        del B

    grid = json.loads(Path(args.configs).read_text())
    jobs = []
    for cfg in grid:
        tag = "_".join(f"{k}{v}" for k, v in cfg.items())
        wd = outdir / "runs" / tag
        wd.mkdir(parents=True, exist_ok=True)
        pred = wd / "pred.h5ad"
        if not pred.exists():
            print(f"  build {tag}", flush=True)
            build(args.real, model, wcols, W, cfg.get("a", 1.08), cfg.get("g", 0.85),
                  cfg.get("q", 1.0), cfg.get("spread", 0.0), cfg.get("k", 0),
                  cfg.get("cpm", 5.0), pred, seed=cfg.get("seed", 0),
                  basal_beta=cfg.get("b", 0.0))
        jobs.append((tag, str(pred), args.real, str(wd)))

    print(f"\nscoring {len(jobs)} configs, {args.jobs} workers\n", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(score_one, j): j[0] for j in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            tag, res = fut.result()
            res["tag"] = tag
            rows.append(res)
            if "error" in res:
                print(f"[{n}/{len(jobs)}] {tag}: {res['error'][:200]}", flush=True)
            else:
                o, r = res["official"], res["raw"]
                print(f"[{n}/{len(jobs)}] {tag}: avg6={o['avg6']:.4f} "
                      f"pds={r['pds_cosine']:.4f}({o['pds_cosine']:.3f}) "
                      f"nmae={r['de_wilcoxon_lfc_nmae']:.4f}"
                      f"({o['de_wilcoxon_lfc_nmae']:.3f}) "
                      f"reach={r['de_wilcoxon_direction_reach_raw']:.4f}"
                      f"({o['de_wilcoxon_direction_reach_raw']:.3f}) "
                      f"fid={r['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
                      f"jac={r['de_wilcoxon_sig_jaccard']:.4f}", flush=True)
            (outdir / "results.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
