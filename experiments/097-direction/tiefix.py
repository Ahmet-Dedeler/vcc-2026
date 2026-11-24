"""Rule-3 check on my own per-gene key: ordinal ranks vs MIDRANKS.

`np.argsort(np.argsort(v))` breaks ties by GENE INDEX.  The agreement vector
sqrt(|A||B|)*1[sign agree] is exactly zero wherever either donor is silent (k562_gwps covers only
7,679 of 18,533 panel genes) or the signs disagree, so a large majority of every row is tied at
zero and my alpha assigned them an arbitrary ramp ordered by position in the gene axis.  A row
with FEWER THAN TWO donors is all-zero and got a pure index ramp.

The `miller` builder found the same defect in the production `agreement_shrink_vector`.  Measured
here on my own numbers rather than assumed to be small.
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
        dn = Donor.load(f"{ROOT}/data/donors/{nd}")
        m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in rg]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        M = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(tg):
            if x in pi:
                M[r, ok] = s[pi[x]][cols[ok]]; h[r] = True
        Ms[nd] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nd] = h
        del dn, s
    A_, B_ = Ms[DON[0]], Ms[DON[1]]
    both = hv[DON[0]] & hv[DON[1]]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dcos = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na*nb, 1e-12), 0.0)
    ag = np.sqrt(np.abs(A_)*np.abs(B_)) * ((np.sign(A_)==np.sign(B_)) & (A_!=0) & (B_!=0))
    ord_r = np.empty_like(ag); mid_r = np.empty_like(ag)
    for i in range(P):
        ord_r[i] = (np.argsort(np.argsort(ag[i])).astype(float) + 0.5) / G
        mid_r[i] = rankdata(ag[i], method="average") / G
        if np.unique(ag[i]).size < 2:            # nothing to rank: no damping at all
            mid_r[i] = 1.0
    zc = (rankdata(dcos, method="average") - 0.5) / P
    nzero = np.array([(ag[i] == 0).mean() for i in range(P)])
    print(f"\n##### {nm}: rows with <2 donors {int((~both).sum())}/{P}; "
          f"share of each row TIED AT ZERO: median {100*np.median(nzero):.1f} %")

    def kap(DH):
        X = float((np.einsum("ij,ij->i", DH, d)*w).sum())
        B = float(_across_pert_budget(DH, t["cols"], w))
        A = float((np.einsum("ij,ij->i", DH, DH)*w).sum())
        return X/np.sqrt(B*t["sumD"]), (A-B)/B
    k0, e0 = kap(dh)
    print(f"  {'arm':<34}{'ordinal (as measured)':>23}{'MIDRANK (correct)':>20}{'eps mid':>10}")
    for lbl, mk in (("per-gene ^1", 1.0), ("per-gene ^2", 2.0)):
        a1, _ = kap((ord_r**mk)*dh); a2, e2 = kap((mid_r**mk)*dh)
        print(f"  {lbl:<34}{a1/k0:>22.3f}x{a2/k0:>19.3f}x{e2:>+10.5f}")
    for lbl, g, mk in (("STACK pert^2 x gene^1", 2.0, 1.0), ("STACK pert^2 x gene^2", 2.0, 2.0)):
        z = (zc**g)[:, None]
        a1, _ = kap(z*(ord_r**mk)*dh); a2, e2 = kap(z*(mid_r**mk)*dh)
        print(f"  {lbl:<34}{a1/k0:>22.3f}x{a2/k0:>19.3f}x{e2:>+10.5f}")
    a, e = kap((zc**2)[:, None]*dh)
    print(f"  {'per-pert donor_cos^2 (unaffected)':<34}{a/k0:>22.3f}x{'-':>20}{e:>+10.5f}")
