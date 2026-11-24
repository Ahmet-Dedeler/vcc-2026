"""450: read the H1 emitter runs against the shipped construction (435's coarse-MDE base, same seeds), paired by seed.
Overall in board units uses 435's H1->board slopes (ka435.SLOPE): sum(slope * d_raw) / 6.
usage: uv run python experiments/450-same-line-programme/read450.py"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
TAG = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
BASE = ROOT / "experiments/435-zs-on-coarse/work/h1"
ARMS = {
    "prior: Pan programme (all targets)": "h1_pan",
    "prior: Pan programme, CLEAN4 fill": "h1_pan_c4fill",
    "rows w0.5, CLEAN4 prior": "h1b_w0.5_clean4",
    "rows w1, CLEAN4 prior": "h1b_w1_clean4",
    "rows w2, CLEAN4 prior": "h1b_w2_clean4",
    "rows w1, Pan prior": "h1b_w1_pan",
    "rows w1 centred, Pan prior": "h1b_w1_pan_c",
    "rows w1 + per-target prior, Pan default": "h1c_w1_pan",
    "rows w1 centred + per-target prior, Pan default": "h1c_w1_pan_c",
    "rows w1 + per-target prior, CLEAN4 default": "h1c_w1_clean4",
    "rows w4, CLEAN4 prior": "h1b_w4_clean4",
    "rows w8, CLEAN4 prior": "h1b_w8_clean4",
    "rows w32 (~replace), CLEAN4 prior": "h1b_w32_clean4",
    "rows w2 + per-target prior, Pan default": "h1c_w2_pan",
    "rows w2 + per-target prior, CLEAN4 default": "h1c_w2_clean4",
    "rows w4 + per-target prior, CLEAN4 default": "h1c_w4_clean4",
    "NOISY bank x1 (signs 0.66): w0.5 + pt, CLEAN4": "h1c_w0.5_clean4_n1s0",
    "NOISY bank x1 (signs 0.66): w1 + pt, CLEAN4": "h1c_w1_clean4_n1s0",
    "NOISY bank x1 (signs 0.66): w2 + pt, CLEAN4": "h1c_w2_clean4_n1s0",
    "NOISY bank x2 (signs 0.61): w0.5 + pt, CLEAN4": "h1c_w0.5_clean4_n2s0",
    "NOISY bank x2 (signs 0.61): w1 + pt, CLEAN4": "h1c_w1_clean4_n2s0",
    "NOISY bank x2 (signs 0.61): w2 + pt, CLEAN4": "h1c_w2_clean4_n2s0",
    "FLIP 0.2 bank (signs ~0.66): w0.5 + pt": "h1c_w0.5_clean4_f0.2s0",
    "FLIP 0.2 bank (signs ~0.66): w1 + pt": "h1c_w1_clean4_f0.2s0",
    "FLIP 0.2 bank (signs ~0.66): w2 + pt": "h1c_w2_clean4_f0.2s0",
    "FLIP 0.3 bank (signs ~0.61): w0.5 + pt": "h1c_w0.5_clean4_f0.3s0",
    "FLIP 0.3 bank (signs ~0.61): w1 + pt": "h1c_w1_clean4_f0.3s0",
    "FLIP 0.3 bank (signs ~0.61): w2 + pt": "h1c_w2_clean4_f0.3s0",
}


def raw(p):
    return json.loads(p.read_text())["raw"]


print(f"{'arm':50s} seeds " + " ".join(f"{m:>9s}" for m in SLOPE) + "   overall(board units)")
for name, d in ARMS.items():
    rows = []
    for s in (0, 1, 2):
        f = HERE / "work" / d / f"metrics-{TAG.format(s=s)}.json"
        b = BASE / f"metrics-{TAG.format(s=s)}.json"
        if f.exists() and b.exists():
            a, c = raw(f), raw(b)
            rows.append({m: a[m] - c[m] for m in SLOPE})
    if not rows:
        continue
    mean = {m: float(np.mean([r[m] for r in rows])) for m in SLOPE}
    ov = [sum(SLOPE[m] * r[m] for m in SLOPE) / 6 for r in rows]
    sd = f" ± {np.std(ov, ddof=1):.4f}" if len(ov) > 1 else ""
    print(f"{name:50s} {len(rows):5d} " + " ".join(f"{mean[m]:+9.4f}" for m in SLOPE) + f"   {np.mean(ov):+.4f}{sd}")
