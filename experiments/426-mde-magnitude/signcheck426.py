"""426: did the MDE arm change delivered signs or the called set / order (the members reach reads)?"""
import sys
import numpy as np
b = np.load("work/h1/stats-base_s0_bps.npz"); m = np.load(f"work/h1/stats-base_s0_bps_{sys.argv[1]}.npz")
Lb, Lm = b["lfc"], m["lfc"]; qb, qm = b["q"], m["q"]; pb, pm = b["p"], m["p"]
fin = np.isfinite(Lb) & np.isfinite(Lm)
flip = fin & (np.sign(Lb) != np.sign(Lm))
print(f"sign flips: {flip.sum()} of {fin.sum()} ({flip.mean():.5f}); among our calls (base q<.05): {(flip & (qb < .05)).sum()}")
cb, cm = qb < .05, qm < .05
print(f"calls base {cb.sum()} mde {cm.sum()}  lost {(cb & ~cm).sum()} gained {(cm & ~cb).sum()}")
# order among base calls: spearman of -log p within each target's calls
from scipy.stats import spearmanr
rs = []
for k in range(Lb.shape[0]):
    c = cb[k] & cm[k]
    if c.sum() > 20:
        rs.append(spearmanr(np.log(pb[k, c] + 1e-300), np.log(pm[k, c] + 1e-300)).correlation)
print(f"within-target order of shared calls, spearman(log p): median {np.median(rs):.4f}  min {np.min(rs):.4f}")
ab = np.abs(np.nan_to_num(Lb)); am = np.abs(np.nan_to_num(Lm))
print(f"|lfc| ratio mde/base median over all genes {np.median(am[ab > 0] / ab[ab > 0]):.3f}; mean |lfc| base {ab.mean():.4f} mde {am.mean():.4f}")
