"""Permutation null for the per-GENE half, shuffling ACROSS GENES.

validate.py's null shuffled the per-gene key across ROWS, which tests row assignment.  The
per-gene arm's claim is that the GENE ORDERING carries information, so the null has to break that:
permute the per-gene agreement values within each perturbation row and re-measure.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor"); sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor
from scipy.stats import rankdata
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
        del dn, s
    A_, B_ = Ms[DON[0]], Ms[DON[1]]
    ag = np.sqrt(np.abs(A_)*np.abs(B_)) * ((np.sign(A_)==np.sign(B_)) & (A_!=0) & (B_!=0))
    al = np.empty_like(ag)
    for i in range(P):
        al[i] = rankdata(ag[i], method="average") / G      # MIDRANKS -- see tiefix.py
        if np.unique(ag[i]).size < 2:
            al[i] = 1.0

    def kap(DH):
        X = float((np.einsum("ij,ij->i", DH, d) * w).sum())
        B = float(_across_pert_budget(DH, t["cols"], w))
        return X / np.sqrt(B * t["sumD"])
    k0 = kap(dh)
    for eta in (1.0, 2.0):
        obs = kap((al ** eta) * dh) / k0
        rng = np.random.default_rng(3); null = np.empty(200)
        for s_ in range(200):
            perm = np.empty_like(al)
            for i in range(P):
                perm[i] = al[i][rng.permutation(G)]
            null[s_] = kap((perm ** eta) * dh) / k0
        print(f"  {nm:<10} eta={eta:<4} observed {obs:.3f}  null p50 {np.percentile(null,50):.3f} "
              f"p95 {np.percentile(null,95):.3f}  p99 {np.percentile(null,99):.3f}  "
              f"p={float((null>=obs).mean()):.4f}")
