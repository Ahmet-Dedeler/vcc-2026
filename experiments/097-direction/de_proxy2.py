"""The DE-side proxy, re-run with the CORRECTED keys (midrank alpha, floored z, damped
single-donor rows).  Per-gene damping never flips a sign, so `fid` precision can only move
through membership; this measures the sign accuracy of the top-K genes by |dh| per perturbation,
which is the ordering-driven part rule 1 says transfers."""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import rankdata
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor"); sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
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
    dh, d = t["dh"], t["d"]
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
    A_, B_ = Ms[DON[0]], Ms[DON[1]]; both = hv[DON[0]] & hv[DON[1]]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dcos = np.where(both, np.einsum("ij,ij->i", A_, B_)/np.maximum(na*nb, 1e-12), 0.0)
    ag = np.sqrt(np.abs(A_)*np.abs(B_)) * ((np.sign(A_)==np.sign(B_)) & (A_!=0) & (B_!=0))
    al = np.empty_like(ag)
    for i in range(P):
        al[i] = rankdata(ag[i], method="average")/G
        if np.unique(ag[i]).size < 2: al[i] = 1.0
    z = np.maximum((rankdata(dcos, method="average")-0.5)/P, 0.3) ** 2
    print(f"\n  {nm}: {P} perturbations")
    print(f"  {'arm':<26}" + "".join(f"{'K='+str(K):>10}" for K in (20, 50, 100, 300, 1000)))
    for lbl, W in (("base", np.ones_like(dh)), ("pert^2 floor.3", z[:, None]*np.ones_like(dh)),
                   ("gene^1 midrank", al), ("STACK p2f3 x g1", z[:, None]*al)):
        DH = W*dh; row = []
        for K in (20, 50, 100, 300, 1000):
            acc = []
            for i in range(P):
                idx = np.argpartition(-np.abs(DH[i]), K)[:K]
                m = d[i][idx] != 0
                if m.sum(): acc.append(float(np.mean(np.sign(DH[i][idx][m]) == np.sign(d[i][idx][m]))))
            row.append(np.mean(acc))
        print(f"  {lbl:<26}" + "".join(f"{v:>10.4f}" for v in row))
