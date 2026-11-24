"""521 — who carries the recurring 'bottom 10 % by strength -> 10 calls' row, and where the null targets sit in the ranking.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/diag521.py
"""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
import lib521 as L
from policy521 import Table, per_target_contrib, SGN
T = Table(0)
P = np.load(L.WORK / "pred521_h1.npz")
tg = [str(t) for t in P["targets"]]
z = P["S4__S3_p_vip_thr4"]
n1 = T.nreal["1.0"]
o = np.argsort(z, kind="stable")
print("bottom 17 by S4 strength rank (bottom 10% = first 8):")
res = {}
for reg in ("1.0", "0.5", "0.35", "cmp"):
    R0 = T.rows[f"{reg}|base"]
    Kv = [None] * 85
    for i in o[:8]:
        Kv[i] = 10
    c = per_target_contrib(T.assemble(reg, "headcpm", "b", Kv), R0)
    res[reg] = c
print(f"{'rank':>4} {'target':>9} {'S4':>7} {'n_real':>6} {'ess':>4} " + " ".join(f"{'c@' + r:>9}" for r in res))
for j, i in enumerate(o[:17]):
    print(f"{j:4d} {tg[i]:>9} {z[i]:+.3f} {n1[i]:6d} {str(bool(P['ess'][i])):>5} " + " ".join(f"{res[r][i]:+9.5f}" for r in res))
print("sum over the 8:", {r: float(res[r][o[:8]].sum()) for r in res})
nulls = np.flatnonzero(n1 == 0)
rk = np.empty(85, int); rk[o] = np.arange(85)
print("null targets (f 1):", [(tg[i], int(rk[i])) for i in nulls], "(rank 0 = weakest by strength)")
lt3 = np.flatnonzero(n1 <= 1)
print("n_real <= 1 targets:", [(tg[i], int(n1[i]), int(rk[i])) for i in lt3])
