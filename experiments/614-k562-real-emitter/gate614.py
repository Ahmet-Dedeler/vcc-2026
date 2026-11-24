#!/usr/bin/env python3
"""614 trust-gate verdict, mechanical (TRUSTGATE.md, cf5f202f). Raw members are evaluate323.score() output (cell-eval2's
DE member code) at f 1, read from each run's metrics json; mean over seeds 0-2.
  G1 zs:   Δnmae = zs − pre-zs (same run) in [-0.0065, -0.00073]
  G2 MDE:  Δnmae = coarse MDE − no MDE   in [-0.0504, -0.0056]
  G3 call set: Δfid < 0 and Δjac > 0     (cp − the zs base; if G3's runs are absent, G3 is PENDING)"""
import json
import sys
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
def raw(d, tag):
    return json.loads((W / d / f"metrics-{tag}.json").read_text())["raw"]
G1, G2, G3 = [], [], []
for s in (0, 1, 2):
    pre = f"base_s{s}_zr-5_prodem_mde0.1cap0.1min0.3coarse"   # run496 tag order: _zr before _prodem (fixed 16:43Z, pre-result)
    try:
        a, b = raw("g1", pre), raw("g1", pre + "_dnk250zs")
        G1.append({m: b[m] - a[m] for m in a})
    except FileNotFoundError:
        pass
    try:
        a, b = raw("g2nomde", f"base_s{s}_prodem"), raw("g2mde", f"base_s{s}_prodem_mde0.1cap0.1min0.3coarse")
        G2.append({m: b[m] - a[m] for m in a})
    except FileNotFoundError:
        pass
    try:
        a, b = raw("g1", pre + "_dnk250zs"), raw("g3cp", pre + "_dnk250zs")
        G3.append({m: b[m] - a[m] for m in a})
    except FileNotFoundError:
        pass
def mean(rows, m):
    return float(np.mean([r[m] for r in rows])) if rows else float("nan")
out = {}
for name, rows in (("G1", G1), ("G2", G2), ("G3", G3)):
    out[name] = {"seeds": len(rows), **{m: mean(rows, m) for m in ("pds", "mse", "nmae", "fid", "reach", "jac")}}
    print(f"{name} ({len(rows)} seeds): " + "  ".join(f"{m} {out[name][m]:+.5f}" for m in ("nmae", "fid", "jac", "reach", "pds", "mse")))
g1 = len(G1) == 3 and -0.0065 <= out["G1"]["nmae"] <= -0.00073
g2 = len(G2) == 3 and -0.0504 <= out["G2"]["nmae"] <= -0.0056
g3 = None if len(G3) < 3 else (out["G3"]["fid"] < 0 and out["G3"]["jac"] > 0)
print(f"G1 {'PASS' if g1 else 'FAIL'}; G2 {'PASS' if g2 else 'FAIL'}; G3 {'PENDING' if g3 is None else ('PASS' if g3 else 'FAIL')}")
verdict = "PENDING (G3)" if (g3 is None and g1 and g2) else ("PASS" if (g1 and g2 and g3) else "FAIL")
print(f"614 TRUST GATE: {verdict}")
(W / "gate614.json").write_text(json.dumps({"verdict": verdict, **out}, indent=1))
