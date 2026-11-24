#!/usr/bin/env python3
"""614 descriptive row: the gate's G1 / G2 pairs on H1 thinned to K562's depth (H1@11k), beside 611's H1@20k numbers.
Raw member deltas at f 1 from each run's metrics json (evaluate323.score), mean of the seeds present."""
import json
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
def raw(d, tag):
    return json.loads((W / d / f"metrics-{tag}.json").read_text())["raw"]
rows = {"G1 zs": [], "G2 MDE": []}
for s in (0, 1, 2):
    pre = f"base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
    try:
        a, b = raw("g1t", pre), raw("g1t", pre + "_dnk250zs"); rows["G1 zs"].append({m: b[m] - a[m] for m in a})
    except FileNotFoundError:
        pass
    try:
        a, b = raw("g2t_nomde", f"base_s{s}_bps_prodem"), raw("g2t_mde", pre); rows["G2 MDE"].append({m: b[m] - a[m] for m in a})
    except FileNotFoundError:
        pass
out = {}
for k, v in rows.items():
    out[k] = {"seeds": len(v), **{m: float(np.mean([r[m] for r in v])) if v else float("nan") for m in ("nmae", "fid", "jac", "reach", "pds", "mse")}}
    print(f"H1@11k {k} ({len(v)} seeds): " + "  ".join(f"{m} {out[k][m]:+.5f}" for m in ("nmae", "fid", "jac", "reach")))
print("H1@20k (611, f 1): zs nmae −0.00183; coarse MDE nmae −0.00966, reach −0.00549")
(W / "h1t614.json").write_text(json.dumps(out, indent=1))
