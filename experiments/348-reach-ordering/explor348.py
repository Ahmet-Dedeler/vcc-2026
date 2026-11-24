#!/usr/bin/env python3
"""348 exploratory (POST-HOC, hypothesis-generating only): sig-first+|design|,
down-first+|design|. Cannot overturn the pre-registered kill verdict."""
from __future__ import annotations
import sys
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reach348 as R348

d = R348.load()
n, g = len(d["targets"]), len(d["rgenes"])
qrm, lrm = R348.build_long(d)
pp, qp, lp = d["pp"], d["qp"], d["lp"]
ref_sig = qrm < R348.ALPHA
own = np.equal.outer(d["targets"], d["rgenes"])
pool = ref_sig & ~own
n_conf = pool.sum(1)
adjud = np.isfinite(lrm) & (lrm != 0)
commit = np.isfinite(lp) & (lp != 0)
match = adjud & commit & (np.sign(lp) == np.sign(lrm))
feat = np.tile(d["rgenes"], (n, 1))
D = R348.design_matrix(d)
ad_m = np.where(np.isfinite(np.abs(D)), -np.abs(D), np.inf)
sig = (np.isfinite(qp) & (qp < R348.ALPHA))
lp_neg = lp < 0


def lex(*keys):
    idx = np.empty((n, g), dtype=int)
    for i in range(n):
        idx[i] = np.lexsort((feat[i],) + tuple(k[i] for k in keys))
    return idx


o7d = lex(ad_m, ~sig)      # sig-first, then |design|
o7e = lex(ad_m, ~lp_neg)   # down-first, then |design|
o4 = lex(ad_m)
o1 = lex(np.where(np.isfinite(np.abs(lp)), -np.abs(lp), np.inf),
         np.where(np.isfinite(pp), pp, np.inf),
         np.where(np.isfinite(qp), qp, np.inf), ~sig)
for name, idx in (("O1_current", o1), ("O4_abs_design", o4),
                  ("O7d_sig_design", o7d), ("O7e_down_design", o7e)):
    Rr, _ = R348.walk_reach(idx, match & pool, adjud & pool, n_conf)
    print(f"{name}: mean reach = {np.nanmean(Rr):.6f}", flush=True)
