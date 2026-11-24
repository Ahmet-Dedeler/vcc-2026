#!/usr/bin/env python3
"""612 idall3 verdict, mechanical (612/PREDICTION-idall3.md, 633f3a39). The real emitter's idall3 − c4h raw nmae (score612.json,
3 seeds) must lie within ±0.0006 of 613d's h1 ungated-v3 − c4h at every regime:
  f 1 in [-0.0022, -0.0010];  f 0.5 in [-0.0029, -0.0017];  f 0.35 in [-0.0029, -0.0017].
Written before idall3's 3-seed score exists. Prints PASS/FAIL per regime and the verdict."""
import json
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
BANDS = {"1.0": (-0.0022, -0.0010), "0.5": (-0.0029, -0.0017), "0.35": (-0.0029, -0.0017)}
TARGET = {"1.0": -0.00160, "0.5": -0.00232, "0.35": -0.00225}
d = json.loads((HERE / "work/score612.json").read_text())
r = d.get("idall3-c4h")
if r is None or any(r[k].get("seeds", 0) < 3 for k in BANDS):
    sys.exit("idall3-c4h is not scored on 3 seeds yet")
ok = True
for k, (lo, hi) in BANDS.items():
    v = r[k]["nmae"]
    p = lo <= v <= hi
    ok &= p
    print(f"f {k:4s}: idall3 − c4h nmae {v:+.5f} (per-seed {[round(x, 5) for x in r[k]['per_seed']['nmae']]}) vs 613d {TARGET[k]:+.5f}, band [{lo:+.4f}, {hi:+.4f}]: {'PASS' if p else 'FAIL'}")
print(f"idall3 VERDICT: {'PASS: idall3 replaces KGID in the queue (coordinator rule); finalize its board prereg before 10-05' if ok else 'FAIL: the hook family stays as is (KGID per A3)'}")
