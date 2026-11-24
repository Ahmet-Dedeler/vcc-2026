#!/usr/bin/env python3
"""612 A3 verdict, mechanical (KGID prereg amendment A3, 2c24e322). The real emitter's idkg3 − c4h raw nmae (score612.json,
3 seeds) must lie within ±0.0006 of 613d's h1 KGID − c4h at every regime:
  f 1 in [-0.0026, -0.0014];  f 0.5 in [-0.0009, +0.0003];  f 0.35 in [-0.0005, +0.0007].
Written before idkg3's 3-seed score exists. Prints PASS/FAIL per regime and the verdict."""
import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
BANDS = {"1.0": (-0.0026, -0.0014), "0.5": (-0.0009, 0.0003), "0.35": (-0.0005, 0.0007)}
TARGET = {"1.0": -0.00202, "0.5": -0.00034, "0.35": 0.00006}
d = json.loads((HERE / "work/score612.json").read_text())
r = d.get("idkg3-c4h")
if r is None or any(r[k].get("seeds", 0) < 3 for k in BANDS):
    sys.exit("idkg3-c4h is not scored on 3 seeds yet")
ok = True
for k, (lo, hi) in BANDS.items():
    v = r[k]["nmae"]
    p = lo <= v <= hi
    ok &= p
    print(f"f {k:4s}: idkg3 − c4h nmae {v:+.5f} (per-seed {[round(x, 5) for x in r[k]['per_seed']['nmae']]}) vs 613d {TARGET[k]:+.5f}, band [{lo:+.4f}, {hi:+.4f}]: {'PASS' if p else 'FAIL'}")
print(f"A3 VERDICT: {'PASS: 613d becomes KGID’s held-out instrument (void if 496a reads < −0.0005 on 10-05)' if ok else 'FAIL: A1’s FAIL branch stays (TIE → plain; KGID needs WIN)'}")
