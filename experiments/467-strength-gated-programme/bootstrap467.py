"""467: target-bootstrap intervals for the H1 claims (added 2026-09-26 after a peer review).

The seed sd (0.0002) measures emission re-rolls only. For a claim about a new panel, the uncertainty that matters is target
sampling: resample h1-test's 85 targets with replacement (4,000 draws per seed), recompute each member as the mean over
the resampled targets (nmae over its eligible ones), and combine with 467's spans. Also: the arms without the four
reach-loss targets (CNN3, COL6A1, ZNF664, AASS).
usage: uv run python experiments/467-strength-gated-programme/bootstrap467.py   -> work/bootstrap467.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
F = "per-target-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.csv"
M = {"fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
     "nmae": "de_wilcoxon_lfc_nmae", "jac": "de_wilcoxon_sig_jaccard"}
SPANS = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37}   # pds and mse deltas are exactly 0 in every arm
FOUR = ["CNN3", "COL6A1", "ZNF664", "AASS"]


def load(d):
    df = pd.concat([pd.read_csv(d / F.format(s=s)).assign(seed=s) for s in (0, 1, 2)])
    return df.pivot_table(index=["seed", "perturbation"], columns="metric", values="value")


def ov(dd, idx):
    return sum(np.nanmean(dd[M[k]].values[idx]) / SPANS[k] for k in SPANS) / 6


base = load(ROOT / "experiments/435-zs-on-coarse/work/h1")
arms = {"ungated": load(ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull"),
        "gated33": load(HERE / "work/h1_gated33")}
rng = np.random.default_rng(467)
out = {}
for seed in (0, 1, 2):
    B = base.loc[seed]
    D = {k: (a.loc[seed].loc[B.index] - B) for k, a in arms.items()}
    D["gate_benefit"] = arms["gated33"].loc[seed].loc[B.index] - arms["ungated"].loc[seed].loc[B.index]
    n = len(B); idx0 = np.arange(n); keep = np.flatnonzero(~B.index.isin(FOUR))
    res = {}
    for k, dd in D.items():
        boots = np.array([ov(dd, rng.integers(0, n, n)) for _ in range(4000)])
        res[k] = {"point": round(ov(dd, idx0), 5), "ci95": [round(float(np.quantile(boots, .025)), 5), round(float(np.quantile(boots, .975)), 5)],
                  "p_le_0": round(float((boots <= 0).mean()), 3), "without_the_four": round(ov(dd, keep), 5)}
    out[f"seed{seed}"] = res
    print(f"seed {seed}: " + "; ".join(f"{k} {v['point']:+.4f} [{v['ci95'][0]:+.4f}, {v['ci95'][1]:+.4f}] P<=0 {v['p_le_0']:.3f} (w/o 4: {v['without_the_four']:+.4f})" for k, v in res.items()))
(HERE / "work/bootstrap467.json").write_text(json.dumps(out, indent=1))
