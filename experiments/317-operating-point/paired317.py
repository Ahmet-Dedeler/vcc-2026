"""317 — paired per-TARGET intervals for the sweep.

EXPERIMENTAL UNIT: one perturbation (target gene) on h1-test, n = 85, PAIRED between an arm
and `base` (identical effect matrix, identical library draw, one knob moved). Intervals are
percentile bootstraps over the 85 targets, 10,000 resamples, seed 0.

Five of the six members are reported per target by `cell-eval2 run` (`results.csv`).
`expr_mse_unbiased_capped_norm` is a single aggregate over the panel and has NO per-target
value, so it gets no interval — only the aggregate delta.
"""
import csv
from pathlib import Path
import numpy as np

PER_TARGET = {"pds": "pds_cosine", "nmae": "de_wilcoxon_lfc_nmae",
              "fid": "de_wilcoxon_direction_fidelity_yield_raw",
              "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}


def per_target(path):
    out = {m: {} for m in PER_TARGET}
    for r in csv.DictReader(Path(path).open()):
        for m, c in PER_TARGET.items():
            if r["metric"] == c:
                out[m][r["perturbation"]] = float(r["value"])
    return out


def paired(a_path, b_path, n=10000, seed=0):
    """b - a, per member, paired over the targets both report."""
    A, B = per_target(a_path), per_target(b_path)
    rng = np.random.default_rng(seed)
    res = {}
    for m in PER_TARGET:
        keys = sorted(set(A[m]) & set(B[m]))
        dv = np.array([B[m][k] - A[m][k] for k in keys])
        # `reach` is NaN for a target the reference gives no confident DE genes; the scorer
        # drops those from its own mean, so the paired test drops them too.
        dv = dv[np.isfinite(dv)]
        if dv.size == 0:
            continue
        bs = np.array([rng.choice(dv, dv.size, replace=True).mean() for _ in range(n)])
        res[m] = {"n": int(dv.size), "mean": float(dv.mean()),
                  "lo": float(np.percentile(bs, 2.5)), "hi": float(np.percentile(bs, 97.5)),
                  "n_up": int((dv > 0).sum()), "n_dn": int((dv < 0).sum())}
    return res
