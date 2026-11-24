"""470 summary: raw member deltas vs 435's shipped CLEAN4 run (same seed), with 467's gated33 and ungated rows, on all 85
h1-test targets and on the 67 non-essential ones (the 2026 panel is essential-free; FINDINGS instruments, 467)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]; HERE = Path(__file__).resolve().parent
F = "per-target-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.csv"
M = {"fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
     "nmae": "de_wilcoxon_lfc_nmae", "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm"}
spans = {"pds": 0.45, "mse": -0.955, "nmae": -0.60, "fid": 0.30, "reach": 0.90, "jac": 0.37}   # 467's H1-like spans
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}


def load(d):
    parts = [pd.read_csv(d / F.format(s=s)).assign(seed=s) for s in (0, 1, 2) if (d / F.format(s=s)).exists()]
    if not parts:
        return None
    return pd.concat(parts).pivot_table(index=["perturbation", "seed"], columns="metric", values="value")


base = load(ROOT / "experiments/435-zs-on-coarse/work/h1")
variants = [("467 ungated arc", ROOT / "experiments/462-broad16-lab-signature/work/h1_arcfull"),
            ("467 gated33 (= wc4)", ROOT / "experiments/467-strength-gated-programme/work/h1_gated33")] + \
           [(d.name[3:], d) for d in sorted((HERE / "work").glob("h1_*"))]
out = {}
for name, d in variants:
    a = load(d)
    if a is None:
        continue
    dd = (a - base).dropna(how="all")
    seeds = sorted(set(dd.index.get_level_values(1)))
    res = {"seeds": len(seeds)}
    for sub in ("all", "nonessential"):
        x = dd if sub == "all" else dd[~np.isin(dd.index.get_level_values(0), list(ESS))]
        per = {k: float(x[v].groupby(level=1).mean().mean()) for k, v in M.items() if v in x.columns}
        per["overall_equiv"] = float(sum(per.get(k, 0.0) / spans[k] for k in spans) / 6)
        res[sub] = {k: round(v, 5) for k, v in per.items()}
    out[name] = res
    print(f"{name:22s} all {res['all']}  |  non-essential {res['nonessential']}  ({res['seeds']} seeds)")
(HERE / "work/summary470.json").write_text(json.dumps(out, indent=1))
