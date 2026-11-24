"""Unanimity vs PAIRWISE per-gene agreement, measured directly on the one surface with 3 donors.

miller measured the graded fraction on the real panel: 3-donor targets read 0.106 under all-donor
unanimity and 0.756 under mean-over-pairs, because unanimity is an INTERSECTION over pairs and
pairwise is a UNION -- so the two extrapolate from the harness's 0.21 in opposite directions and
a two-donor harness cannot separate them (with two donors they are the same function).

h1-test CAN separate them: `orion_hct116` is a legal third donor there (it is not that harness).
So build both keys with three donors on h1-test and score them.  orion cannot be used for this --
its only legal third donor covers 12 of 132 targets.
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
real = f"{ROOT}/data/harness/h1-test.h5ad"
pred = f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
for k, c in t["cols"].items():
    t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
dh, d, w = t["dh"], t["d"], t["w"]
tg = np.array([str(x) for x in pp[t["keep"]]]); P, G = dh.shape
M, H = {}, {}
for nd in ("k562_gwps", "orion_hek293t", "orion_hct116"):
    dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cols = np.array([gi.get(str(g), -1) for g in rg]); ok = cols >= 0
    pi = {p: i for i, p in enumerate(dn.perts)}
    Mx = np.zeros((P, G)); h = np.zeros(P, bool)
    for r, xx in enumerate(tg):
        if xx in pi: Mx[r, ok] = s[pi[xx]][cols[ok]]; h[r] = True
    M[nd] = Mx / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
    H[nd] = h; del dn, s

def unanimity(dons):
    out = np.zeros((P, G))
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        if len(have) < 2: continue
        st = np.vstack([M[n][i] for n in have]); sg = np.sign(st)
        un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
        out[i] = np.where(un, np.exp(np.log(np.abs(st)+1e-300).mean(0)), 0.0)
    return out

def pairwise(dons):
    out = np.zeros((P, G))
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        if len(have) < 2: continue
        tot = np.zeros(G); cnt = 0
        for a in range(len(have)):
            for b in range(a+1, len(have)):
                x, y = M[have[a]][i], M[have[b]][i]
                g = np.sqrt(np.abs(x)*np.abs(y)) * ((np.sign(x)==np.sign(y)) & (x!=0) & (y!=0))
                tot += g; cnt += 1
        out[i] = tot/max(cnt, 1)
    return out

def alpha(agr):
    al = np.empty_like(agr)
    for i in range(P):
        al[i] = rankdata(agr[i], method="average")/G
        if np.unique(agr[i]).size < 2: al[i] = 1.0
    return al

def kap(DH):
    X = float((np.einsum("ij,ij->i", DH, d)*w).sum())
    B = float(_across_pert_budget(DH, t["cols"], w))
    return X/np.sqrt(B*t["sumD"])
k0 = kap(dh)
D2 = ("k562_gwps", "orion_hek293t"); D3 = D2 + ("orion_hct116",)
print(f"# h1-test, {P} targets.  donor coverage "
      + ", ".join(f"{n} {int(H[n].sum())}" for n in D3))
print(f"\n  {'key':<38}{'graded frac':>13}" + "".join(f"{'eta='+str(e):>9}" for e in (0.5,1.0,2.0)))
for lbl, agr in (("2-donor (= unanimity = pairwise)", unanimity(D2)),
                 ("3-donor UNANIMITY (production `agreement`)", unanimity(D3)),
                 ("3-donor PAIRWISE (`agreement-pairwise`)", pairwise(D3))):
    al = alpha(agr)
    frac = float(np.median([(agr[i] > 0).mean() for i in range(P)]))
    row = [kap((al**e)*dh)/k0 for e in (0.5, 1.0, 2.0)]
    print(f"  {lbl:<38}{frac:>13.3f}" + "".join(f"{v:>9.3f}" for v in row))
# and stacked with the per-perturbation weight, at the recommended setting
A_, B_ = M[D2[0]], M[D2[1]]
na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
dc = np.where(H[D2[0]] & H[D2[1]],
              np.einsum("ij,ij->i", A_, B_)/np.maximum(na*nb, 1e-12), 0.0)
z = np.maximum((rankdata(dc, method="average")-0.5)/P, 0.3) ** 2
print(f"\n  stacked with z=max(rank01(donor_cos),0.3)^2, eta=1:")
for lbl, agr in (("2-donor", unanimity(D2)), ("3-donor unanimity", unanimity(D3)),
                 ("3-donor pairwise", pairwise(D3))):
    print(f"    {lbl:<22}{kap(z[:, None]*alpha(agr)*dh)/k0:>8.3f}x")
