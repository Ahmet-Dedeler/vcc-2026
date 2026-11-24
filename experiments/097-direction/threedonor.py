"""Does the key survive being built the way PRODUCTION would build it -- from THREE donors?

The harness model is a two-donor LODO transfer, so everything measured so far uses a two-donor
key.  Production has three donor slots and would use the MEAN PAIRWISE cosine.  A third donor can
be added to the KEY without touching the prediction, which is the right directional test.

LEAK RULES: h1-test may use orion_hct116 (it is not the h1 harness); orion-hct116-thinned may NOT.
For orion the only leak-free third donor with any coverage is h1_2025 -- coverage is printed, and
if it is near zero the arm is reported dead rather than run.
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
CFG = [("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        ["k562_gwps", "orion_hek293t"], ["orion_hct116", "rpe1", "nadig_jurkat", "k562_essential"]),
       ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        ["k562_gwps", "orion_hek293t"], ["h1_2025", "rpe1", "nadig_jurkat", "k562_essential"])]
for nm, real, pred, base_d, extra_d in CFG:
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    dh, d, w = t["dh"], t["d"], t["w"]
    tg = np.array([str(x) for x in pp[t["keep"]]]); P, G = dh.shape
    x = np.einsum("ij,ij->i", dh, d); a = np.einsum("ij,ij->i", dh, dh)
    cos_p = x/np.maximum(np.sqrt(a)*np.linalg.norm(d, axis=1), 1e-30)
    def kap(DH):
        X = float((np.einsum("ij,ij->i", DH, d)*w).sum())
        B = float(_across_pert_budget(DH, t["cols"], w))
        return X/np.sqrt(B*t["sumD"])
    k0 = kap(dh)
    M, H, COV = {}, {}, {}
    for nd in base_d + extra_d:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in rg]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        Mx = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, xx in enumerate(tg):
            if xx in pi: Mx[r, ok] = s[pi[xx]][cols[ok]]; h[r] = True
        M[nd] = Mx / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        H[nd] = h; COV[nd] = int(h.sum()); del dn, s
    print(f"\n##### {nm}: {P} targets.  donor coverage: "
          + ", ".join(f"{k} {COV[k]}" for k in COV))
    def key(dons):
        tot, cnt = np.zeros(P), np.zeros(P)
        for i in range(len(dons)):
            for j in range(i+1, len(dons)):
                A_, B_ = M[dons[i]], M[dons[j]]; msk = H[dons[i]] & H[dons[j]]
                na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
                c = np.einsum("ij,ij->i", A_, B_)/np.maximum(na*nb, 1e-12)
                tot += np.where(msk, c, 0.0); cnt += msk
        return np.where(cnt > 0, tot/np.maximum(cnt, 1), 0.0), cnt
    for lbl, dons in [("2-donor (as measured)", base_d)] + \
                     [(f"3-donor +{e}", base_d + [e]) for e in extra_d if COV[e] >= 0.5*P]:
        kk, cnt = key(dons)
        z = np.maximum((rankdata(kk, method="average")-0.5)/P, 0.3)
        r = spearmanr(kk, cos_p).statistic
        row = [kap((z**g)[:, None]*dh)/k0 for g in (1.0, 2.0, 3.0)]
        print(f"  {lbl:<26} rho={r:+.3f}  pairs/target median {np.median(cnt):.0f}  "
              f"M(g=1,2,3) = " + " ".join(f"{v:.3f}" for v in row))
    for e in extra_d:
        if COV[e] < 0.5*P:
            print(f"  {'+'+e:<26} DEAD: covers {COV[e]}/{P} targets")
