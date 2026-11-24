"""513c: who occupies the head of our reach walk on H1? Block genes (top 2,200 by control CPM: the called-for-every-target
block, whose sign is the shrunk CLEAN4 programme except the q80-protected strongest) vs the rest. Accuracy of each class
by walk position, split by pool size.
    uv run python experiments/513-reach-walk-anatomy/block513.py
"""
from __future__ import annotations

import numpy as np

from walk513 import ALPHA, load, order_keys

d = load(); G = d["genes"]
ctrl = d["ctrl"]
rank = np.empty(len(ctrl), int); rank[np.argsort(-ctrl)] = np.arange(len(ctrl))
block = rank < 2200
print(f"genes on axis {len(G)}; block genes {block.sum()}; ctrl CPM at block edge {np.sort(ctrl)[::-1][2199]:.1f}")
bins = [(0, 10), (10, 30), (30, 100), (100, 10**6)]
stats = {}
for k, t in enumerate(d["targets"]):
    own = G == t
    pool = np.nonzero((d["rq"][k] < ALPHA) & ~own)[0]
    if len(pool) == 0:
        continue
    nr = len(pool)
    band = "<20" if nr < 20 else "20-100" if nr < 100 else ">=100"
    order = order_keys(d["oq"][k], d["op"][k], d["ol"][k], G)
    inpool = np.zeros(len(G), bool); inpool[pool] = True
    pos = 0
    for g in order:
        if not inpool[g]:
            continue
        t_ = d["rl"][k][g]
        if not np.isfinite(t_) or t_ == 0:
            continue
        ok = 1 if (d["ol"][k][g] != 0 and np.sign(d["ol"][k][g]) == np.sign(t_)) else 0
        b = next(i for i, (lo, hi) in enumerate(bins) if lo <= pos < hi)
        key = (band, "block" if block[g] else "rest", b)
        s = stats.setdefault(key, [0, 0]); s[0] += ok; s[1] += 1
        called = d["oq"][k][g] < ALPHA
        key2 = (band, "called" if called else "uncalled")
        s = stats.setdefault(key2, [0, 0]); s[0] += ok; s[1] += 1
        pos += 1
for band in ["<20", "20-100", ">=100"]:
    print(f"== pool {band}")
    for cls in ["block", "rest"]:
        print("   " + cls.ljust(6) + "  ".join(
            f"pos[{lo},{hi if hi < 10**6 else 'end'}) n={stats.get((band, cls, i), [0, 0])[1]:5d} acc={stats.get((band, cls, i), [0, 1])[0] / max(stats.get((band, cls, i), [0, 1])[1], 1):.3f}"
            for i, (lo, hi) in enumerate(bins)))
    for cls in ["called", "uncalled"]:
        s = stats.get((band, cls), [0, 0])
        print(f"   {cls:8s} n={s[1]:6d} acc={s[0] / max(s[1], 1):.3f}")
