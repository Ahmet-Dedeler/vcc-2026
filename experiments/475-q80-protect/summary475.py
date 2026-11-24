"""475 summary: fix vs 435's production-emitter run at the same seed (== `orig` by the KA), per member, on all 85
h1-test targets and on the 67 non-essential ones, with a target bootstrap (3,000 draws per seed) of the overall-equiv."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
F = "per-target-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.csv"
M = {"fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
     "nmae": "de_wilcoxon_lfc_nmae", "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine"}
SP = {"fid": 0.30, "reach": 0.90, "nmae": -0.60, "jac": 0.37, "pds": 0.45}
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}


def load(d):
    return pd.concat([pd.read_csv(d / F.format(s=s)).assign(seed=s) for s in (0, 1, 2) if (d / F.format(s=s)).exists()]) \
             .pivot_table(index=["seed", "perturbation"], columns="metric", values="value")


base, fix = load(ROOT / "experiments/435-zs-on-coarse/work/h1"), load(HERE / "work/h1_fix")
rng = np.random.default_rng(475)
ov = lambda d, idx: sum(np.nanmean(d[M[k]].values[idx]) / SP[k] for k in SP) / 6
out = {}
for seed in sorted(set(fix.index.get_level_values(0))):
    B = base.loc[seed]; D = fix.loc[seed].loc[B.index] - B
    for sub in ("all", "nonessential"):
        idx = np.arange(len(B)) if sub == "all" else np.flatnonzero(~B.index.isin(list(ESS)))
        mem = {k: round(float(np.nanmean(D[v].values[idx])), 5) for k, v in M.items()}
        bo = np.array([ov(D, idx[rng.integers(0, len(idx), len(idx))]) for _ in range(3000)])
        out.setdefault(sub, {})[f"seed{seed}"] = {**mem, "overall_equiv": round(ov(D, idx), 5),
                                                  "ci95": [round(float(np.quantile(bo, .025)), 5), round(float(np.quantile(bo, .975)), 5)],
                                                  "p_le_0": round(float((bo <= 0).mean()), 3)}
        print(sub, seed, out[sub][f"seed{seed}"], flush=True)
(HERE / "work/summary475.json").write_text(json.dumps(out, indent=1))
