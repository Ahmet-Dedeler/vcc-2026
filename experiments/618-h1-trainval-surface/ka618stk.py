#!/usr/bin/env python3
"""618 stacks KA: the plain dump run (stk_base, --dump 126) must reproduce the gate's G1 zs arm (g1) bit for bit: same
flags, and --dump does not change stats (612's chain612e KA). Both zs and pre-zs stats, seeds present in both."""
import sys
from pathlib import Path
import numpy as np
W = Path(__file__).resolve().parent / "work"
TAG = "stats-base_s{s}_bps_zr-5_prodem_mde0.1cap0.1min0.3coarse{z}.npz"
ok, n = True, 0
for s in (0, 1, 2):
    for z in ("", "_dnk250zs"):
        a, b = W / "g1" / TAG.format(s=s, z=z), W / "stk_base" / TAG.format(s=s, z=z)
        if not (a.exists() and b.exists()):
            continue
        A, B = np.load(a), np.load(b)
        same = all(np.array_equal(A[k], B[k], equal_nan=np.issubdtype(A[k].dtype, np.floating)) for k in A.files if k in B.files)
        n += 1
        ok &= same
        print(f"s{s}{z or ' (pre-zs)'}: stk_base == g1: {'PASS' if same else 'FAIL'}")
print(f"KA stk_base == g1 over {n} stats files: {'PASS' if ok and n else 'FAIL'}")
sys.exit(0 if ok and n else 1)
