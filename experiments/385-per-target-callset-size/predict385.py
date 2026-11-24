#!/usr/bin/env python3
"""385 step 2: is n_real,p predictable from build-time observables?

Reads (read-only): 323 work/reference.npz + reference.json (n_real cross-check),
effects.npz::cap15, stage.npz::uniform (pre-whitening donor-consensus norm),
360 work/CRISPRGeneEffect_subset.csv. Writes: work/predict385.json + stdout table.
No emission, no scorer, <1 GiB.
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
WORK = HERE / "work"

TAUS = (0.05, 0.1, 0.2, 0.5)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    z = np.load(BASE / "work/reference.npz", allow_pickle=False)
    genes_r = z["genes"].astype(str)
    targets = z["targets"].astype(str)
    keep = z["keep"]
    q = z["q"]
    uni = genes_r[keep]
    n_real = np.empty(len(targets), dtype=int)
    for i, t in enumerate(targets):
        sig = q[i] < 0.05
        sig &= uni != t  # scorer excludes the perturbation's own target gene
        n_real[i] = int(sig.sum())
    # cross-check against the recorded n_real table
    rec = json.loads((BASE / "work/reference.json").read_text())
    check = np.array([rec["n_real_by_target"][t] for t in targets])
    assert np.array_equal(n_real, check), "n_real mismatch vs reference.json"
    print(f"n_real: n=85 median={np.median(n_real):.0f} mean={n_real.mean():.0f} "
          f"min={n_real.min()} max={n_real.max()} zeros={(n_real == 0).sum()}")

    ze = np.load(BASE / "work/effects.npz", allow_pickle=False)
    genes_e = ze["genes"].astype(str)
    assert list(ze["targets"].astype(str)) == list(targets)
    cap15 = ze["cap15"]
    zs = np.load(BASE / "work/stage.npz", allow_pickle=False)
    assert list(zs["targets"].astype(str)) == list(targets)
    pos = {g: i for i, g in enumerate(genes_e)}
    shared = np.array([pos[g] for g in genes_r if g in pos])
    uni_cols = np.array([pos[g] for g in uni if g in pos])
    print(f"axis map: {len(shared)}/{len(genes_r)} ref genes in effect axis; "
          f"{len(uni_cols)}/{len(uni)} universe genes mapped")

    obs = {}
    obs["cap15_L2_full"] = np.linalg.norm(cap15, axis=1)
    obs["cap15_L2_universe"] = np.linalg.norm(cap15[:, uni_cols], axis=1)
    obs["prewhite_L2_full"] = np.linalg.norm(zs["uniform"], axis=1)
    for tau in TAUS:
        obs[f"cap15_n_gt_{tau}"] = (np.abs(cap15[:, uni_cols]) > tau).sum(1).astype(float)

    # DepMap essentiality of the target: Jurkat row + 17-model mean (HeLa absent, H1 absent)
    with open(ROOT / "experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv") as f:
        r = csv.reader(f)
        hdr = next(r)
        rows = list(r)
    sym_of = [c.split(" (")[0] for c in hdr[1:]]
    scol = {s: j for j, s in enumerate(sym_of)}
    models = [x[0] for x in rows]
    mat = np.array([[float(v) if v not in ("", "NA") else np.nan for v in x[1:]] for x in rows])
    jur = mat[models.index("ACH-000995")]
    dep_jur, dep_mean, dep_n = [], [], 0
    for t in targets:
        if t in scol:
            dep_n += 1
            dep_jur.append(abs(jur[scol[t]]))
            dep_mean.append(float(np.nanmean(np.abs(mat[:, scol[t]]))))
        else:
            dep_jur.append(np.nan)
            dep_mean.append(np.nan)
    obs["depmap_absJurkat"] = np.array(dep_jur)
    obs["depmap_absMean17"] = np.array(dep_mean)
    print(f"depmap coverage: {dep_n}/85 H1 targets in CRISPRGeneEffect columns")

    y = np.log(n_real + 1)
    print(f"\n{'observable':>22} {'spearman':>9} {'p':>9} | LOO Pearson Spearman residSD")
    out = {"n": len(targets), "predictors": {}}
    best = 0.0
    for name, x in obs.items():
        ok = np.isfinite(x)
        rho, p = spearmanr(x[ok], n_real[ok])
        # LOO OLS of log(n+1) on log(x+eps), valid rows only
        xx = np.log(x[ok] + 1e-12)
        yy = y[ok]
        pred = np.empty_like(yy)
        A = np.column_stack([np.ones_like(xx), xx])
        for i in range(len(xx)):
            m = np.ones(len(xx), bool)
            m[i] = False
            beta, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
            pred[i] = beta @ A[i]
        loo_p = float(np.corrcoef(pred, yy)[0, 1])
        loo_s, _ = spearmanr(pred, yy)
        resid = float(np.std(yy - pred, ddof=1))
        base = float(np.std(yy, ddof=1))
        print(f"{name:>22} {rho:+.4f} {p:.2e} | {loo_p:+.4f} {float(loo_s):+.4f} {resid:.3f} (const {base:.3f})")
        out["predictors"][name] = {"spearman": float(rho), "p": float(p),
                                   "loo_pearson": loo_p, "loo_spearman": float(loo_s),
                                   "loo_resid_sd": resid, "const_sd": base,
                                   "n_valid": int(ok.sum())}
        best = max(best, abs(float(rho)))
    # 3-predictor LOO combo on rows valid everywhere
    combo_names = ["cap15_L2_universe", "cap15_n_gt_0.1", "depmap_absJurkat"]
    ok = np.all([np.isfinite(obs[n]) for n in combo_names], axis=0)
    X = np.column_stack([np.log(obs[n][ok] + 1e-12) for n in combo_names])
    yy = y[ok]
    A = np.column_stack([np.ones(len(X)), X])
    pred = np.empty(len(X))
    for i in range(len(X)):
        m = np.ones(len(X), bool)
        m[i] = False
        beta, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
        pred[i] = beta @ A[i]
    loo_p = float(np.corrcoef(pred, yy)[0, 1])
    loo_s, _ = spearmanr(pred, yy)
    resid = float(np.std(yy - pred, ddof=1))
    print(f"{'LOO-3combo':>22} {'':>9} {'':>9} | {loo_p:+.4f} {float(loo_s):+.4f} {resid:.3f} (n={ok.sum()})")
    out["combo"] = {"names": combo_names, "loo_pearson": loo_p,
                    "loo_spearman": float(loo_s), "loo_resid_sd": resid, "n_valid": int(ok.sum())}
    best = max(best, abs(float(loo_s)))
    out["best_abs_spearman"] = best
    out["kill"] = bool(best < 0.4)
    (WORK / "predict385.json").write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(f"\nbest |Spearman| = {best:.4f} -> {'KILL (stop, no pricing)' if best < 0.4 else 'CLEARS 0.4 gate'}")


if __name__ == "__main__":
    main()
