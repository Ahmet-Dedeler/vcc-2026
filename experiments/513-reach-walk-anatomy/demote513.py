"""513d: analytic price of demoting the high-CPM block in our reach order (H1, production emitter stats).
Every variant keeps our means, call set and signs (so fid, jac, nmae, pds, mse are unchanged by construction); only the
order in which pool genes are walked changes. Variants (stable within each class, by our own p order):
  prod       the scorer's order on today's p-values (known answer 0.26926)
  demote     called non-block, then called block, then uncalled (non-block and block by p)
  demote2    called non-block, then uncalled non-block, then block (called, then uncalled)
  blocklast  every non-block gene (called first) before every block gene
  nonblock   drop block genes from the walk entirely (NOT deliverable: shows the ceiling of the class split)
Block = top-B by control CPM on the 10,779-gene axis; B scanned.
    uv run python experiments/513-reach-walk-anatomy/demote513.py
"""
from __future__ import annotations

import json

import numpy as np

from walk513 import ALPHA, load, order_keys
from keys513 import kstar, walk_for

d = load(); G = d["genes"]; ctrl = d["ctrl"]
rank = np.empty(len(ctrl), int); rank[np.argsort(-ctrl)] = np.arange(len(ctrl))
out = {}
for B in (1000, 2200, 3000, 4000):
    block = rank < B
    res = {k: [] for k in ("prod", "demote", "demote2", "blocklast", "nonblock")}
    nrs = []
    for k, t in enumerate(d["targets"]):
        own = G == t
        pool = np.nonzero((d["rq"][k] < ALPHA) & ~own)[0]
        if len(pool) == 0:
            continue
        nrs.append(len(pool))
        ol, rl = d["ol"][k], d["rl"][k]
        prod = order_keys(d["oq"][k], d["op"][k], ol, G)
        called = d["oq"][k] < ALPHA
        pos = np.empty(len(G), int); pos[prod] = np.arange(len(G))
        cls_demote = np.where(called & ~block, 0, np.where(called & block, 1, 2))
        cls_demote2 = np.where(called & ~block, 0, np.where(~called & ~block, 1, np.where(called & block, 2, 3)))
        cls_last = np.where(~block, np.where(called, 0, 1), np.where(called, 2, 3))
        orders = {"prod": prod,
                  "demote": np.lexsort((pos, cls_demote)),
                  "demote2": np.lexsort((pos, cls_demote2)),
                  "blocklast": np.lexsort((pos, cls_last))}
        for name, od in orders.items():
            res[name].append(kstar(walk_for(od, pool, rl, ol)) / len(pool))
        nb_pool = pool[~block[pool]]
        w = walk_for(prod, nb_pool, rl, ol)
        res["nonblock"].append(kstar(w) / len(pool))
    nrs = np.array(nrs)
    row = {}
    for name, v in res.items():
        v = np.array(v)
        row[name] = dict(reach=float(v.mean()), lt20=float(v[nrs < 20].mean()),
                         mid=float(v[(nrs >= 20) & (nrs < 100)].mean()), ge100=float(v[nrs >= 100].mean()))
    out[B] = row
    print(f"B={B}")
    for name, r in row.items():
        print(f"   {name:9s} reach {r['reach']:.4f} (Δ {r['reach'] - row['prod']['reach']:+.4f})   <20 {r['lt20']:.3f}  20-100 {r['mid']:.3f}  >=100 {r['ge100']:.3f}")
json.dump(out, open("work/demote513.json", "w"), indent=1)
