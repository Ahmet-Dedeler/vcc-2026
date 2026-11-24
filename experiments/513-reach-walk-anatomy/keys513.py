"""513b: sign accuracy by rank position inside the reach walk, and reach under alternative confidence orders.
Orders compared (all computable from our own emitted stats, i.e. deliverable only if the emitter can realise them):
  prod   the scorer's order (sig first, q, p, |lfc|)             -- known answer 0.26926
  lfc    |our lfc| descending (sig first)                          -- magnitude as confidence
  lfcall |our lfc| descending, ignoring the sig-first rule
  oracle truth |lfc| order with OUR signs (upper bound on ordering alone; not deliverable)
    uv run python experiments/513-reach-walk-anatomy/keys513.py
"""
from __future__ import annotations

import numpy as np

from walk513 import ALPHA, P0, load, order_keys


def kstar(walk):
    if len(walk) == 0:
        return 0
    cum = np.cumsum(walk) / np.arange(1, len(walk) + 1)
    ok = np.nonzero(cum >= P0)[0]
    return int(ok[-1] + 1) if len(ok) else 0


def walk_for(order, pool, rl, ol):
    inpool = np.zeros(len(rl), bool); inpool[pool] = True
    w = []
    for g in order:
        if not inpool[g]:
            continue
        t = rl[g]
        if not np.isfinite(t) or t == 0:
            continue
        w.append(1 if (ol[g] != 0 and np.sign(ol[g]) == np.sign(t)) else 0)
    return np.asarray(w, int)


def main():
    d = load(); G = d["genes"]
    orders = {}
    res = {k: [] for k in ("prod", "lfc", "lfcall", "oracle")}
    pos_acc = {k: np.zeros(6) for k in ("prod", "lfc")}; pos_n = {k: np.zeros(6) for k in ("prod", "lfc")}
    bins = [(0, 10), (10, 25), (25, 50), (50, 100), (100, 250), (250, 10**6)]
    nrs = []
    for k, t in enumerate(d["targets"]):
        own = G == t
        pool = np.nonzero((d["rq"][k] < ALPHA) & ~own)[0]
        if len(pool) == 0:
            continue
        nrs.append(len(pool))
        ol, rl = d["ol"][k], d["rl"][k]
        sig = d["oq"][k] < ALPHA
        o = {"prod": order_keys(d["oq"][k], d["op"][k], ol, G),
             "lfc": np.lexsort((G, -np.abs(np.nan_to_num(ol)), ~sig)),
             "lfcall": np.lexsort((G, -np.abs(np.nan_to_num(ol)))),
             "oracle": np.lexsort((G, -np.abs(np.nan_to_num(rl))))}
        for name, od in o.items():
            w = walk_for(od, pool, rl, ol)
            res[name].append(kstar(w) / len(pool))
            if name in pos_acc and len(pool) >= 100:
                for i, (lo, hi) in enumerate(bins):
                    seg = w[lo:hi]
                    pos_acc[name][i] += seg.sum(); pos_n[name][i] += len(seg)
    nrs = np.array(nrs)
    for name, v in res.items():
        v = np.array(v)
        print(f"{name:7s} reach {v.mean():.4f}  | n_real<20 {v[nrs<20].mean():.3f}  20-100 {v[(nrs>=20)&(nrs<100)].mean():.3f}"
              f"  100-1000 {v[(nrs>=100)&(nrs<1000)].mean():.3f}  >=1000 {v[nrs>=1000].mean():.3f}")
    print("sign accuracy by walk position, targets with n_real >= 100:")
    for name in pos_acc:
        print(f"  {name:5s} " + "  ".join(f"[{lo},{hi if hi<10**6 else 'end'}) {pos_acc[name][i]/max(pos_n[name][i],1):.3f}"
                                      for i, (lo, hi) in enumerate(bins)))


if __name__ == "__main__":
    main()
