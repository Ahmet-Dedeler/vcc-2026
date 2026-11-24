"""Are 097's `gamma_gene` and the builder's `--shrink-gamma` the same knob? NO. Measured.

  097 (this document):   alpha_g = ( ascending_midrank_g / G ) ** gg          power law on a FRACTION
  builder:               alpha_g = exp( -gamma * descending_midrank_g / 2000 ) exponential on a COUNT

Different families, so there is no scalar conversion -- and the gene axis size enters one and not
the other. This measures both on both harnesses and reports the quantity that decides everything
(kappa) plus two interpretable summaries: the mean alpha, and the alpha assigned to the TIE BLOCK,
which is the ~79-89 % of genes whose donor agreement is exactly 0 and which therefore all receive
one shared value.
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
RANK_SCALE = 2000.0


def desc_midrank(x):
    n = x.size
    order = np.argsort(-x, kind="stable"); od = x[order]
    edges = np.flatnonzero(np.concatenate(([True], od[1:] != od[:-1], [True])))
    st, en = edges[:-1], edges[1:]
    gr = (st + en - 1) / 2.0
    out = np.empty(n); out[order] = np.repeat(gr, en - st)
    return out


for nm, real, pred, dons in [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
     ["k562_gwps", "orion_hek293t", "orion_hct116"]),
    ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
     ["k562_gwps", "orion_hek293t"])]:
    rp, rs, rc, rj, rg = group_stats(real); pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    dh0, w, cols, sD = t["dh"], t["w"], t["cols"], t["sumD"]
    tg = np.array([str(x) for x in pp[t["keep"]]]); P, G = dh0.shape
    M, H = {}, {}
    for nd in dons:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s_ = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cc = np.array([gi.get(str(g), -1) for g in rg]); ok = cc >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        Mx = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(tg):
            if x in pi: Mx[r, ok] = s_[pi[x]][cc[ok]]; h[r] = True
        M[nd] = Mx / (float(np.median(np.linalg.norm(s_[:, cc[ok]], axis=1))) or 1.0)
        H[nd] = h; del dn, s_
    ag = np.zeros((P, G))
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        if len(have) < 2: continue
        st = np.vstack([M[n][i] for n in have]); sg = np.sign(st)
        un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
        ag[i] = np.where(un, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)

    def kap(A):
        DH = A * dh0
        X = float((np.einsum("ij,ij->i", DH, t["d"]) * w).sum())
        B = float(_across_pert_budget(DH, cols, w))
        return X / np.sqrt(B * sD)
    k0 = kap(np.ones_like(ag))
    nz = float(np.mean([(ag[i] > 0).mean() for i in range(P)]))
    print(f"\n##### {nm}: G={G}, {len(dons)} donors, nonzero-agreement genes {100*nz:.1f} % "
          f"(so the TIE BLOCK is {100*(1-nz):.1f} % of every row)   base kappa {k0:.5f}")
    print(f"  {'form':<44}{'mean alpha':>12}{'alpha on tie block':>20}{'kappa':>10}{'M':>8}")
    for gg in (0.25, 1.0, 2.0):
        A = np.empty_like(ag)
        for i in range(P):
            A[i] = (rankdata(ag[i], method="average") / G) ** gg
            if np.unique(ag[i]).size < 2: A[i] = 1.0
        tb = float(np.median([A[i][ag[i] == 0].mean() for i in range(P) if (ag[i] == 0).any()]))
        print(f"  {'097  (rank/G)^gg      gg=' + f'{gg}':<44}{A.mean():>12.4f}{tb:>20.4f}"
              f"{kap(A):>10.5f}{kap(A)/k0:>8.3f}")
    for gam in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0):
        A = np.empty_like(ag)
        for i in range(P):
            A[i] = (np.ones(G) if np.unique(ag[i]).size < 2
                    else np.exp(-gam * desc_midrank(ag[i]) / RANK_SCALE))
        tb = float(np.median([A[i][ag[i] == 0].mean() for i in range(P) if (ag[i] == 0).any()]))
        print(f"  {'BLDR exp(-g*rank/2000)  gamma=' + f'{gam}':<44}{A.mean():>12.4f}{tb:>20.4f}"
              f"{kap(A):>10.5f}{kap(A)/k0:>8.3f}")
