"""280 Q1b -- the donor test AT THE HEADLINE CELL SCALE (~240 cells/side), which §3.1's 2x2 cannot reach.

Added analysis (not in the preregistration verbatim; the preregistration's §3.2 asks for a within-donor
comparator "at the same per-side cell count" and a donor has only ~130 cells per target, so the comparator is
built by re-partitioning the SAME eight (donor, guide) pseudobulks in three ways at the same per-side size):

  DONOR split    side A = both guides of 2 donors      | side B = both guides of the other 2   (donors differ)
  GUIDE split    side A = guide-1 of ALL FOUR donors   | side B = guide-2 of all four          (donors identical)
  BOTH split     side A = guide-1 of 2 donors          | side B = guide-2 of the other 2       (donors + guides)

All three have 4 pseudobulk rows and ~the same cells per side, and every side is thinned to a common
per-target counting budget.  median(GUIDE) - median(DONOR) > 0  <=>  the donor is a real unit at this scale.
"""
from __future__ import annotations

import itertools
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from lib280 import CONDS, DORDER, SEED, Bank, boot_median, boot_median_diff, cosines, residualise, thin  # noqa: E402

P = lambda *a: print(*a, flush=True)  # noqa: E731
KEEP = "--all" not in sys.argv
THIN = "--thin" in sys.argv
TAG = ("keepDE" if KEEP else "allguides") + ("_thin" if THIN else "")
bank = Bank(".")
SC, CC = bank.scored, bank.cols_cd4
o = bank.obs
tg = o[o.gtype == "targeting"]
if KEEP:
    tg = tg[tg.keep_for_DE]
IDX = {(c, d, t): g.sort_values("guide") for (c, d, t), g in tg.groupby(["cond", "D", "pert"])}
ALLT = sorted(set(tg.pert))
rows = []

for cond in CONDS:
    ok = [t for t in ALLT if all((cond, d, t) in IDX and len(IDX[(cond, d, t)]) >= 2 for d in DORDER)]
    P(f"\n=== {cond}: {len(ok)} targets with 2 guides in ALL FOUR donors ===")
    cell = {t: {d: [(int(IDX[(cond, d, t)].i.values[r]), float(IDX[(cond, d, t)].total.values[r]),
                     float(IDX[(cond, d, t)].n_cells.values[r])) for r in (0, 1)] for d in DORDER} for t in ok}
    # common per-side counting budget: 2 x the smallest single (donor, guide) total (a side always has 2 rows)
    budget = {t: 2.0 * min(v[1] for d in DORDER for v in cell[t][d]) for t in ok}
    rng = np.random.default_rng(SEED)

    def side(spec, thin_to=THIN):
        """spec: list of (donor, guide-rank).  Returns residual rows (n_ok, n_scored) and cells/target."""
        donors = sorted({d for d, _ in spec})
        b = bank.ctrl_profile(cond, donors)[CC]
        M = np.zeros((len(ok), len(CC)))
        nc = np.zeros(len(ok))
        for k, t in enumerate(ok):
            S = np.zeros(bank.X.shape[1])
            for d, r in spec:
                i, _, c = cell[t][d][r]
                S = S + bank.counts([i])
                nc[k] += c
            if thin_to:
                S = thin(S, budget[t], rng)
            M[k] = bank.profile(S)[CC] - b
        R = residualise(M)
        R[:, ~SC] = 0.0
        return R[:, SC], nc

    res = {}
    # DONOR split: 3 disjoint pairings
    acc, cA = [], []
    for pr in (("D1", "D2"), ("D1", "D3"), ("D1", "D4")):
        oth = tuple(d for d in DORDER if d not in pr)
        A, ca = side([(d, r) for d in pr for r in (0, 1)])
        B, cb = side([(d, r) for d in oth for r in (0, 1)])
        acc.append(cosines(A, B)); cA.append((ca + cb) / 2)
    res["DONOR"] = (np.nanmean(np.vstack(acc), 0), np.mean(np.vstack(cA), 0))
    # GUIDE split: donors identical on both sides
    A, ca = side([(d, 0) for d in DORDER])
    B, cb = side([(d, 1) for d in DORDER])
    res["GUIDE"] = (cosines(A, B), (ca + cb) / 2)
    # BOTH split: donors and guides differ (averaged over the 3 pairings x 2 guide assignments)
    acc, cA = [], []
    for pr in (("D1", "D2"), ("D1", "D3"), ("D1", "D4")):
        oth = tuple(d for d in DORDER if d not in pr)
        for r in (0, 1):
            A, ca = side([(d, r) for d in pr])
            B, cb = side([(d, 1 - r) for d in oth])
            acc.append(cosines(A, B)); cA.append((ca + cb) / 2)
    res["BOTH"] = (np.nanmean(np.vstack(acc), 0), np.mean(np.vstack(cA), 0))

    for k in ("DONOR", "GUIDE", "BOTH"):
        c, nc = res[k]
        m, lo, hi, n = boot_median(c)
        P(f"  {k:<6} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  mean {np.nanmean(c):+.4f}"
          f"  <=0 {np.nanmean(c<=0):.3f}  cells/side median {np.median(nc):.0f}")
        rows.append(dict(cond=cond, split=k, n=n, median=m, lo=lo, hi=hi, mean=float(np.nanmean(c)),
                         frac_le0=float(np.nanmean(c <= 0)), cells_per_side=float(np.median(nc))))
    for lab, a, b_ in (("GUIDE - DONOR  (donor is a real unit if > 0)", "GUIDE", "DONOR"),
                       ("GUIDE - BOTH", "GUIDE", "BOTH"), ("DONOR - BOTH", "DONOR", "BOTH")):
        d, lo, hi, n = boot_median_diff(res[a][0], res[b_][0])
        star = "   <-- PRIMARY (matched scale)" if a == "GUIDE" and b_ == "DONOR" and cond == "Stim48hr" else ""
        P(f"    {lab:<45} {d:+.4f} [{lo:+.4f},{hi:+.4f}]  n={n}{star}")
        rows.append(dict(cond=cond, split=lab, n=n, median=d, lo=lo, hi=hi, mean=np.nan,
                         frac_le0=np.nan, cells_per_side=np.nan))
    np.savez(f"q1b_{cond}_{TAG}.npz", targets=np.asarray(ok),
             **{k: res[k][0] for k in res}, **{f"cells_{k}": res[k][1] for k in res})

pd.DataFrame(rows).to_csv(f"q1b_matched_scale_{TAG}.csv", index=False)
