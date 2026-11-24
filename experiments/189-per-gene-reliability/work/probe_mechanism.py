"""Why does a prior that is overwhelming on the donor surface (z = +41.6) buy nothing on
h1-test?  Two candidates: (a) it is already inside |H|, (b) it measures the REFERENCE's
measurement precision, which our prediction does not control.

(a) is tested by correlating theta_g with |H|'s own gene-level marginal, and by asking whether
theta orders sign accuracy WITHIN a fixed |H| decile.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from harness189 import Surface189, check
from common189 import spearman

T = Surface189(); check(T)
absH = np.abs(T.H)
m = T.prior_covered
gH = absH.mean(0)                       # the gene's mean |H| across targets
gHpool = np.array([absH[T.pool[:,g], g].mean() if T.pool[:,g].any() else np.nan
                   for g in range(len(T.genes))])
print(f"\nspearman(theta, gene-mean |H|)                {spearman(T.theta[m], gH[m]):+.4f}")
ok = m & np.isfinite(gHpool)
print(f"spearman(theta, gene-mean |H| inside the pool) {spearman(T.theta[ok], gHpool[ok]):+.4f}")
print(f"spearman(theta, expression)                   {spearman(T.theta[m], T.b_ctrl[m]):+.4f}")
print(f"spearman(theta, reference |real_lfc| marginal) ", end="")
rl = np.array([np.abs(T.S.real_lfc[T.pool[:,g], g]).mean() if T.pool[:,g].any() else np.nan
               for g in range(len(T.genes))])
ok2 = m & np.isfinite(rl)
print(f"{spearman(T.theta[ok2], rl[ok2]):+.4f}")

# does theta order sign accuracy WITHIN a fixed |H| decile? (pooled over the DE pool)
print(f"\nsign accuracy by theta tercile, WITHIN each |H| decile of the DE pool:")
rows = []
for p in range(len(T.perts)):
    idx = np.flatnonzero(T.pool[p])
    if idx.size < 30: continue
    o = idx[np.argsort(-absH[p, idx])]
    e = (np.arange(11)*o.size)//10
    for d in range(10):
        blk = o[e[d]:e[d+1]]
        if blk.size < 3: continue
        th = T.theta[blk]
        q = np.argsort(th)
        t = blk.size//3
        if t == 0: continue
        rows.append((d, T.correct[p, blk[q[:t]]].mean(), T.correct[p, blk[q[-t:]]].mean(),
                     t))
rows = np.array([(a,b,c,d) for a,b,c,d in rows])
print(f"{'|H| decile':>11}{'low-theta acc':>15}{'high-theta acc':>16}{'diff':>9}{'n blocks':>10}")
tot_lo = tot_hi = tot_n = 0
for d in range(10):
    s = rows[rows[:,0]==d]
    if not len(s): continue
    lo = np.average(s[:,1], weights=s[:,3]); hi = np.average(s[:,2], weights=s[:,3])
    tot_lo += (s[:,1]*s[:,3]).sum(); tot_hi += (s[:,2]*s[:,3]).sum(); tot_n += s[:,3].sum()
    print(f"{d+1:>11}{lo:>15.4f}{hi:>16.4f}{hi-lo:>+9.4f}{len(s):>10}")
print(f"{'ALL':>11}{tot_lo/tot_n:>15.4f}{tot_hi/tot_n:>16.4f}{(tot_hi-tot_lo)/tot_n:>+9.4f}")
