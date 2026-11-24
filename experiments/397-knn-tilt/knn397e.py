#!/usr/bin/env python3
"""397e -- HEAD flips where the six-bank neighbourhood (K=25, shrink 400) disagrees with the delivered sign at confidence >= tau,
on top of 399's block signs. Reach reads the head first, so this is the risky member. H1 analytic."""
import numpy as np
exec(open("~/Code/virtual-cell-challenge-2026/experiments/397-knn-tilt/knn397d.py").read().split("for m in (50.0, 400.0):")[0])
m = 400.0; PK = np.full((nT, nG), np.nan)
for i in range(nT):
    if not has[i]: continue
    t = targets[i]; v = GWn[gwp[t]]; s_ = GWn @ v; s_[gwp[t]] = -np.inf; nb = np.argsort(-s_)[:K]; P, nn = votes(nb, FAM)
    PK[i] = np.where(np.isfinite(P) & np.isfinite(P_marg6), (nn * np.nan_to_num(P) + m * np.nan_to_num(P_marg6)) / (nn + m), np.where(np.isfinite(P), P, P_marg6))
conf = np.abs(PK - 0.5); sgk = np.where(np.isfinite(PK), np.where(PK > 0.5, -1.0, 1.0), sg6[None, :])
SG399 = np.where(cb, sgk, SG0); truth = np.sign(lrm)
M399 = metrics(SG399, np.abs(E0))
def contrast(name, SG):
    M = metrics(SG, np.abs(E0)); row = {}
    for k in M:
        dv = M[k] - M399[k]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (4000, len(okk))); bs = dv[okk][D].mean(1)
        row[k] = f"{np.nanmean(dv):+.5f} [{np.quantile(bs,.025):+.5f}, {np.quantile(bs,.975):+.5f}]"
    print(name, row, flush=True)
hc = head & sig & adjud & ~own & np.isfinite(PK); dis = hc & (sgk != SG0)
for tau in (.1, .15, .2, .3):
    mk = dis & (conf >= tau)
    print(f"head disagree conf>={tau}: n={mk.sum():5d} ({mk.sum()/hc.sum():.3%} of head calls) delivered right {(SG0 == truth)[mk].mean():.3f} -> kNN right {(sgk == truth)[mk].mean():.3f}")
    contrast(f"  flip them", np.where(mk, sgk, SG399))
# and the rank position of those genes in the delivered order (do they sit at the reach front?)
rank_cur = np.argsort(cur, axis=1)
mk = dis & (conf >= .2); print("median delivered rank of flipped head genes:", int(np.median(rank_cur[mk])), " (head = first 400 by design)")
