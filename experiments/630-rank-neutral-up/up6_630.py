"""630 v6 PROBE: library drift as the compensation (driver-only; module 7e6b3e52 untouched). NOT a shipping candidate yet.

The capacity bound is the compensation currency: today each count a recipient receives is paid with a sub-gate count
moved back to the donor ballast row (libraries exact), and untested genes are ~0.55 % of a carrier's library. This probe
asks what happens if the recipient's library is simply allowed to GROW (and the donor ballast row's to shrink):
  * per-gene block totals stay exact (counts only move between rows), so pds is unchanged; mse's jackknife may move;
  * every other gene's CPM in a recipient falls by L/(L+k) (~0.5 % for k ~ 100 on a 20k library), so its log value
    drops by ~0.005 and may cross a few controls: a small U change spread over MANY genes;
  * the ballast row's library shrinks by a tiny fraction of ~995k.
`up_pass6` = up_pass2 by source transformation with (1) an unlimited sub_left and (2) no sub-gate payment, (3) the
"a library moved" assert removed. Gene j's own U accounting inside the pass uses the ORIGINAL libraries (approximate);
the driver re-scores the real block, and the call-flip guard and the six members are computed on the real counts.
`all_gene_dauc(B0, B1, cc)` reports every kept gene's |ΔAUC| between two blocks on their own libraries (the coordinator's
bar is |ΔAUC| <= 1e-4 per gene), so the probe shows how much of the bar library drift spends.
Variant grammar (rnup630): <v3 variant>_L.
"""
from __future__ import annotations

import inspect

import numpy as np

import src_rank_up as M
from src_rank_up import ZS


def _build():
    src = inspect.getsource(M.up_pass2)
    edits = [
        ("def up_pass2(", "def up_pass6("),
        ("    sub_left = X[np.ix_(car, sidx)].sum(1).astype(np.int64)\n",
         "    sub_left = np.full(n_car, np.int64(1) << 40)                 # v6 probe: unlimited (library drift)\n"),
        ("    for rc, bb, t in sub_moves:\n        give = ZS._take(X[rc, sidx], t)\n        X[rc, sidx] -= give\n        X[bb, sidx] += give\n",
         "    pass                                                          # v6 probe: no sub-gate payment\n"),
        ("    if not np.array_equal(X.sum(1), Lr):\n        raise AssertionError(\"a library moved\")\n",
         "    pass                                                          # v6 probe: libraries drift by design\n"),
    ]
    for old, new in edits:
        if src.count(old) != 1:
            raise AssertionError(f"up6: anchor not found exactly once: {old[:60]!r}")
        src = src.replace(old, new)
    # the final per-gene U re-check uses the ORIGINAL scales s; with drifting libraries it is not exact -> drop it
    a = "    for j in np.flatnonzero(np.any(X[:, gidx] != B[:, gidx], axis=0)):\n"
    i = src.index(a)
    k = src.index("    au = np.abs(dU[dU != 0])")
    src = src[:i] + "    # v6 probe: U re-check skipped (libraries drift); the driver measures every gene's real |dAUC|\n" + src[k:]
    ns = dict(M.__dict__)
    exec(compile(src, "<up_pass6>", "exec"), ns)
    return ns["up_pass6"]


up_pass6 = _build()


def all_gene_dauc(B0, B1, cc):
    """|ΔAUC| of every kept gene between blocks B0 and B1, each on its OWN libraries (the scorer's view)."""
    keep = np.asarray(cc.keep, dtype=bool)
    gidx = np.flatnonzero(keep)
    s0, s1 = np.asarray(B0).sum(1) / 1e6, np.asarray(B1).sum(1) / 1e6
    n, m = B0.shape[0], cc.cs.shape[1]
    out = np.zeros(gidx.size)
    X0, X1 = np.asarray(B0)[:, gidx], np.asarray(B1)[:, gidx]
    ch = np.flatnonzero(np.any(X0 != X1, axis=0) | True)            # every gene: drifted libraries touch all of them
    for j in ch:
        d = int(ZS.ranks_of(X1[:, j], s1, cc.cs[j]).sum()) - int(ZS.ranks_of(X0[:, j], s0, cc.cs[j]).sum())
        out[j] = abs(d) / (2.0 * n * m)
    return out
