#!/usr/bin/env python3
"""661 KA: flip618.py h1 none s0 (612's exact c4kg flags) == 612's h1_c4kg s0 stats, bit for bit (q, p, lfc, sums, jk),
pre-zs and zs. Exit 0 iff both match."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ok, n = True, 0
for z in ("", "_dnk250zs"):
    t = f"stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse{z}.npz"
    a, b = ROOT / "experiments/612-key-gated-hook/work/h1_c4kg" / t, HERE / "work/h1_octnone" / t
    if not (a.exists() and b.exists()):
        print(f"missing {t}: 612 {a.exists()} 661 {b.exists()}")
        ok = False
        continue
    A, B = np.load(a), np.load(b)
    same = all(np.array_equal(A[k], B[k], equal_nan=np.issubdtype(A[k].dtype, np.floating)) for k in ("q", "p", "lfc", "sums", "jk"))
    n += 1
    ok &= same
    print(f"{t}: flip618 h1 none == 612 c4kg: {'PASS' if same else 'FAIL'}")
print(f"KA661 ({n} files): {'PASS' if ok and n == 2 else 'FAIL'}")
sys.exit(0 if ok and n == 2 else 1)
