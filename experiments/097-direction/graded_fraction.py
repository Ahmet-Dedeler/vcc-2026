"""Does the GRADED FRACTION matter?  The largest unquantified thing between 097 and the builder.

The per-gene alpha grades the genes with a nonzero donor-agreement value and damps the rest to one
shared midrank.  That fraction is 20.8 % / 20.6 % on the two harnesses.  On the real 300-target
panel it is 10.6 % under all-donor unanimity and 75.6 % under the mean-over-pairs form (miller,
`shrink_coverage.log`) -- because unanimity is an INTERSECTION over pairs and pairwise is a UNION,
so the two extrapolate from the measured 21 % in OPPOSITE directions.  Neither reproduces it, and
a two-donor harness cannot tell them apart (with two donors the two definitions are the same
function).

So the question is not which key: it is whether the multiplier is sensitive to the graded fraction
at all.  Emulated here by thresholding the agreement vector at a quantile, which moves the graded
fraction from 5 % to 95 % on a fixed key and a fixed surface.
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import rankdata
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor"); sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor
DON = ("k562_gwps", "orion_hek293t")
for nm, real, pred in [("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
    f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
    f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")]:
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    dh, d, w = t["dh"], t["d"], t["w"]
    tg = np.array([str(x) for x in pp[t["keep"]]]); P, G = dh.shape
    Ms, hv = {}, {}
    for nd in DON:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in rg]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        M = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, xx in enumerate(tg):
            if xx in pi: M[r, ok] = s[pi[xx]][cols[ok]]; h[r] = True
        Ms[nd] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nd] = h; del dn, s
    A_, B_ = Ms[DON[0]], Ms[DON[1]]
    ag = np.sqrt(np.abs(A_)*np.abs(B_)) * ((np.sign(A_)==np.sign(B_)) & (A_!=0) & (B_!=0))
    # a CONTINUOUS surrogate for the same ranking, so the graded fraction can be dialled:
    # geometric magnitude with the sign gate applied as a large penalty rather than as a zero.
    geo = np.sqrt(np.abs(A_)*np.abs(B_))
    sgn = ((np.sign(A_) == np.sign(B_)) & (A_ != 0) & (B_ != 0)).astype(float)
    def kap(DH):
        X = float((np.einsum("ij,ij->i", DH, d)*w).sum())
        B = float(_across_pert_budget(DH, t["cols"], w))
        return X/np.sqrt(B*t["sumD"])
    k0 = kap(dh)
    nz = float(np.mean(ag > 0))
    print(f"\n##### {nm}: the real key grades {100*nz:.1f} % of genes")
    print(f"  {'graded fraction f':<22}" + "".join(f"{'eta='+str(e):>10}" for e in (0.5,1.0,2.0)))
    for f in (0.05, 0.106, 0.21, 0.40, 0.60, 0.756, 0.95):
        # keep the top f of genes per row by (sign-gated geo, then geo), grade THEM by midrank,
        # damp everything else to the midrank of the damped block -- the same shape the real key has
        row = []
        for eta in (0.5, 1.0, 2.0):
            al = np.empty_like(geo)
            for i in range(P):
                key = geo[i]*sgn[i]
                thr = np.quantile(key, 1-f)
                sel = key > max(thr, 0.0)
                r = np.empty(G)
                if sel.sum() == 0:
                    r[:] = 1.0
                else:
                    r[~sel] = (0.5*(~sel).sum())/G
                    r[sel] = ((~sel).sum() + rankdata(key[sel], method="average"))/G
                al[i] = r
            row.append(kap((al**eta)*dh)/k0)
        print(f"  f={f:<20.3f}" + "".join(f"{v:>10.3f}" for v in row))
