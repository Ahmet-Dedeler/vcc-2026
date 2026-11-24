"""426 known answer: run426.py `base S --bulk ps` (mde off) must equal 416's base_sS_bps bit-for-bit (all stats arrays)
and its members exactly. usage: uv run python ka426.py [SEED]"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
a = np.load(HERE / f"work/h1/stats-base_s{seed}_bps.npz")
b = np.load(HERE.parent / f"416-de-carrier-menu/work/h1/stats-base_s{seed}_bps.npz")
ok = True
for k in ("targets", "sums", "jk", "p", "q", "lfc"):
    eq = np.array_equal(a[k], b[k], equal_nan=(k != "targets"))
    ok &= eq
    print(f"  {k:6s} bit-equal: {eq}")
ma = json.loads((HERE / f"work/h1/metrics-base_s{seed}_bps.json").read_text())["raw"]
mb = json.loads((HERE.parent / f"416-de-carrier-menu/work/h1/metrics-base_s{seed}_bps.json").read_text())["raw"]
for m in mb:
    # members to 1e-12: the scorer's polars aggregation can differ in the last ulp between identical inputs (419)
    print(f"  {m:6s} {ma[m]!r} vs 416 {mb[m]!r}  {'==' if abs(ma[m] - mb[m]) <= 1e-12 else 'DIFF'}")
    ok &= abs(ma[m] - mb[m]) <= 1e-12
out = {"seed": seed, "pass": bool(ok)}
(HERE / f"work/known-answer-base_s{seed}_bps.json").write_text(json.dumps(out, indent=1))
print("KA PASS" if ok else "KA FAIL")
sys.exit(0 if ok else 1)
