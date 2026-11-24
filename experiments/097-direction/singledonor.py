"""RETRACTION CHECK. How much of the h1-test per-perturbation gain was the 4 single-donor rows?

`perpert.py` left a target carried by fewer than two donors at `donor_cos = 0`, which ranks it
near the 20th percentile and therefore DAMPS it hard.  The production builder pins such a target
at z = 1.  h1-test has 4 such rows of 85; orion has 0 of 132.  Switching to the builder's policy
moved h1-test's multiplier from 1.238 to 0.951 at gamma = 2 -- so four rows out of eighty-five
were carrying the entire h1-test result, and the two surfaces no longer agree.

This isolates it three ways, so the conclusion does not rest on one contrast.
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import rankdata, spearmanr
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
        for r, x in enumerate(tg):
            if x in pi: M[r, ok] = s[pi[x]][cols[ok]]; h[r] = True
        Ms[nd] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nd] = h; del dn, s
    A_, B_ = Ms[DON[0]], Ms[DON[1]]; both = hv[DON[0]] & hv[DON[1]]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dcos = np.where(both, np.einsum("ij,ij->i", A_, B_)/np.maximum(na*nb, 1e-12), 0.0)
    x = np.einsum("ij,ij->i", dh, d); a = np.einsum("ij,ij->i", dh, dh)
    cos_p = x/np.maximum(np.sqrt(a)*np.linalg.norm(d, axis=1), 1e-30)

    def W(z): 
        n = float((w*z*x).sum()); q = float((w*z*z*a).sum())
        return (n*n/q) if q > 0 else 0.0
    base = W(np.ones(P))
    print(f"\n##### {nm}: {int((~both).sum())} of {P} rows carry <2 donors")
    if (~both).sum():
        print(f"  those rows: {list(tg[~both])}")
        print(f"    their share of total A: {100*(w*a)[~both].sum()/(w*a).sum():.1f} %   "
              f"of total X: {100*(w*x)[~both].sum()/(w*x).sum():+.1f} %   "
              f"mean cos_p {cos_p[~both].mean():+.4f} vs {cos_p[both].mean():+.4f} for the rest")
    r = rankdata(dcos, method="average")/P
    print(f"  {'policy for <2-donor rows':<34}" + "".join(f"{'g='+str(g):>9}" for g in (0.5,1,2,3)))
    for lbl, zf in (("(A) damp them: z=rank(dcos=0)  [my first pass]", lambda g: r**g),
                    ("(B) keep them: z=1            [the builder]", lambda g: np.where(both, r**g, 1.0)),
                    ("(C) drop them from the panel entirely", None)):
        if zf is None: continue
        print(f"  {lbl:<34}" + "".join(f"{np.sqrt(W(zf(g))/base):>9.3f}" for g in (0.5,1,2,3)))
    # (C) restrict the whole measurement to the two-donor rows: the cleanest isolation
    m = both
    r2 = rankdata(dcos[m], method="average")/m.sum()
    b2 = None
    def W2(z):
        n = float((w[m]*z*x[m]).sum()); q = float((w[m]*z*z*a[m]).sum())
        return (n*n/q) if q > 0 else 0.0
    b2 = W2(np.ones(int(m.sum())))
    print(f"  {'(C) 2-donor rows only, ranked among themselves':<34}"
          + "".join(f"{np.sqrt(W2(r2**g)/b2):>9.3f}" for g in (0.5,1,2,3)))
    print(f"  rho(donor_cos, cos_p): all rows {spearmanr(dcos, cos_p).statistic:+.3f}   "
          f"2-donor rows only {spearmanr(dcos[m], cos_p[m]).statistic:+.3f}")
