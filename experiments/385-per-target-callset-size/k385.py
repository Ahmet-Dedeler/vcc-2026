#!/usr/bin/env python3
"""385 step 3a: build per-target K_p rules from the LOO DepMap predictor.

Predictor: log(n_real+1) ~ OLS on log(depmap_absMean17) (the gate-clearing
observable, Spearman +0.51), fitted LEAVE-ONE-OUT so every target's n̂ is
out-of-sample. K_p = clip(round(c * n̂), 400, 4474), c in {1, 2}.
CallSet(n_head_p, n_cpm_p, 1.9), n_head_p = min(400, K_p) — the protected
400-|H| head stays inside whenever K_p >= 400 (always, given the floor).
Shuffled control: same K multiset (c=1) dealt via fixed-seed derangement.

Writes work/k385.json {targets, n_real, n_hat_loo, K_c1, K_c2, K_shuf}.
"""
import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
WORK = HERE / "work"

KMIN, KMAX, HEAD = 400, 4474, 400


def derangement(n, seed=385):
    rng = np.random.default_rng(seed)
    while True:
        p = rng.permutation(n)
        if bool(np.all(p != np.arange(n))):
            return p


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    z = np.load(BASE / "work/reference.npz", allow_pickle=False)
    genes_r = z["genes"].astype(str)
    targets = z["targets"].astype(str)
    uni = genes_r[z["keep"]]
    q = z["q"]
    n_real = np.array([( (q[i] < 0.05) & (uni != t) ).sum() for i, t in enumerate(targets)])

    with open(ROOT / "experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv") as f:
        r = csv.reader(f)
        hdr = next(r)
        rows = list(r)
    scol = {c.split(" (")[0]: j for j, c in enumerate(hdr[1:])}
    models = [x[0] for x in rows]
    assert "ACH-000995" in models  # Jurkat sanity anchor; H1 verified absent (see RESULTS)
    mat = np.array([[float(v) if v not in ("", "NA") else np.nan for v in x[1:]] for x in rows])
    x = np.array([np.nanmean(np.abs(mat[:, scol[t]])) if t in scol else np.nan for t in targets])
    missing = [t for t in targets if t not in scol]
    print("median-imputed (absent from DepMap columns):", missing)
    x[~np.isfinite(x)] = np.nanmedian(x)  # build-time-legitimate: median of the 83 valid

    # LOO OLS in log space
    lx, y = np.log(x + 1e-12), np.log(n_real + 1)
    A = np.column_stack([np.ones_like(lx), lx])
    pred = np.empty_like(y)
    for i in range(len(y)):
        m = np.ones(len(y), bool)
        m[i] = False
        beta, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
        pred[i] = beta @ A[i]
    n_hat = np.expm1(pred)
    n_hat = np.clip(n_hat, 0, None)

    out = {"targets": targets.tolist(), "n_real": n_real.tolist(),
           "n_hat_loo": [float(v) for v in n_hat]}
    for c in (1, 2):
        K = np.clip(np.round(c * n_hat).astype(int), KMIN, KMAX)
        out[f"K_c{c}"] = K.tolist()
        print(f"c={c}: K min={K.min()} med={np.median(K):.0f} mean={K.mean():.0f} max={K.max()} "
              f"at_floor={(K == KMIN).sum()} at_cap={(K == KMAX).sum()}")
    K1 = np.array(out["K_c1"])
    out["K_shuf"] = K1[derangement(len(K1))].tolist()
    assert sorted(out["K_shuf"]) == sorted(out["K_c1"]), "shuffled must preserve multiset"
    print(f"n_real: min={n_real.min()} med={np.median(n_real):.0f} mean={n_real.mean():.0f} max={n_real.max()}")
    print(f"n_hat_loo: min={n_hat.min():.1f} med={np.median(n_hat):.1f} mean={n_hat.mean():.1f} max={n_hat.max():.1f}")
    (WORK / "k385.json").write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
