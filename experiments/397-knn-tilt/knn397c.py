#!/usr/bin/env python3
"""397c -- the neighbourhood signs on NON-CALLED genes (reach walks the uncalled tail; nmae reads the gate): replace the delivered
sign of uncalled, non-head genes where the kNN confidence |P-0.5| >= tau; 399's signs on the called block. H1 analytic."""
import numpy as np
exec(open("~/Code/virtual-cell-challenge-2026/experiments/397-knn-tilt/knn397.py").read().split("# ---- kNN tilts")[0].replace('print("KNOWN ANSWER', 'pass  # print("KNOWN ANSWER'))
has = np.array([t in gwp for t in targets]); K = 25; m = 50.0
PK = np.full((nT, nG), np.nan)
for i, t in enumerate(targets):
    if not has[i]: continue
    v = GWn[gwp[t]]; s_ = GWn @ v; s_[gwp[t]] = -np.inf; nb = np.argsort(-s_)[:K]
    P, nn = votes(nb, ["k562_essential", "nadig_jurkat"])
    PK[i] = np.where(np.isfinite(P) & np.isfinite(P_marg_3), (nn * np.nan_to_num(P) + m * np.nan_to_num(P_marg_3)) / (nn + m), np.where(np.isfinite(P), P, P_marg_3))
conf = np.abs(PK - 0.5); sgk = np.where(np.isfinite(PK), np.where(PK > 0.5, -1.0, 1.0), sg4[None, :])
SGK = np.where(cb, sgk, SG0); truth = np.sign(lrm)
unc = ~sig & ~head & adjud & ~own & np.isfinite(PK)
for lo, hi in ((0, .1), (.1, .2), (.2, .3), (.3, .5)):
    mk = unc & (conf >= lo) & (conf < hi)
    print(f"UNCALLED conf[{lo},{hi}): n={mk.sum():7d}  kNN precision={(sgk == truth)[mk].mean():.3f}  delivered={(SG0 == truth)[mk].mean():.3f}  in ref gate={pool[mk].mean():.3f}  gate-weighted kNN={(sgk == truth)[mk & pool].mean() if (mk & pool).any() else float('nan'):.3f} delivered={(SG0 == truth)[mk & pool].mean() if (mk & pool).any() else float('nan'):.3f}")
M399 = metrics(SGK, np.abs(E0))
def contrast(name, SG):
    M = metrics(SG, np.abs(E0)); row = {}
    for k in M:
        dv = M[k] - M399[k]; okk = np.flatnonzero(np.isfinite(dv)); D = rng.integers(0, len(okk), (4000, len(okk))); bs = dv[okk][D].mean(1)
        row[k] = f"{np.nanmean(dv):+.5f} [{np.quantile(bs,.025):+.5f}, {np.quantile(bs,.975):+.5f}]"
    print(name, row)
for tau in (.1, .2, .3):
    SG = np.where(unc & (conf >= tau), sgk, SGK); contrast(f"uncalled kNN sign where conf>={tau} (vs 399 signs)", SG)
