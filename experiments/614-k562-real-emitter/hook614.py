#!/usr/bin/env python3
"""614 hook forms vs PRED-614-levers.md (registered first): raw member deltas at f 1 on the zs stats, mean of seeds 0-2.
  c4h − no hook = hk_c4h − G1's zs;  idall3 − c4h = hk_id − hk_c4h.  Read as the gate verdict allows (gate614.json)."""
import json
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
def raw(d, s):
    return json.loads((W / d / f"metrics-base_s{s}_zr-5_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.json").read_text())["raw"]
rows = {"c4h − no hook": [], "idall3 − c4h": [], "idall3 − no hook": []}
for s in (0, 1, 2):
    try:
        b, c, i = raw("g1", s), raw("hk_c4h", s), raw("hk_id", s)
    except FileNotFoundError:
        continue
    rows["c4h − no hook"].append({m: c[m] - b[m] for m in b})
    rows["idall3 − c4h"].append({m: i[m] - c[m] for m in b})
    rows["idall3 − no hook"].append({m: i[m] - b[m] for m in b})
v = json.loads((W / "gate614.json").read_text())["verdict"] if (W / "gate614.json").exists() else "absent"
print(f"614 gate verdict: {v}")
out = {"gate": v}
for k, r in rows.items():
    out[k] = {"seeds": len(r), **{m: float(np.mean([x[m] for x in r])) if r else float("nan") for m in ("nmae", "fid", "jac", "reach")}}
    print(f"{k:18s} ({len(r)} seeds): " + "  ".join(f"{m} {out[k][m]:+.5f}" for m in ("nmae", "fid", "jac", "reach")))
p1 = out["idall3 − c4h"]["nmae"] > 0
p2 = out["c4h − no hook"]["nmae"] > 0
print(f"registered: idall3 − c4h nmae > 0 at f 1: {'MATCH' if p1 else 'MISS'}; c4h − no hook nmae > 0 at f 1: {'MATCH' if p2 else 'MISS'}")
(W / "hook614.json").write_text(json.dumps(out, indent=1))
