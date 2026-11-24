"""kappa and pds for the BUILDER's exponential form -- i.e. a13 as actually built.

097's tables all use alpha = (ascending_midrank/G)**gg. The shipped builder uses
alpha = exp(-gamma * descending_midrank / 2000). a13 was built with `--shrink-gamma 0.25` on the
builder, so a13's alpha is the EXPONENTIAL one at 0.25 -- not 097's power law at 0.25, which is a
far milder object (mean alpha 0.834 against 0.353).
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import rankdata
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor"); sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats, CTRL, TS
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.discrimination import discrimination_score
from cell_eval2.metrics.delta import _across_pert_budget
from vcc2026.donor_transfer import Donor


def desc_midrank(x):
    n = x.size; o = np.argsort(-x, kind="stable"); od = x[o]
    e = np.flatnonzero(np.concatenate(([True], od[1:] != od[:-1], [True])))
    st, en = e[:-1], e[1:]
    out = np.empty(n); out[o] = np.repeat((st + en - 1) / 2.0, en - st)
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
    ci = int(np.flatnonzero(rp == CTRL)[0])
    bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
    rmeans = bulk_lognorm_means(rs.astype(np.float64), TS)
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
    ag = np.zeros((P, G)); dc = np.zeros(P); ndon = np.zeros(P, int)
    for i in range(P):
        have = [n for n in dons if H[n][i]]; ndon[i] = len(have)
        if len(have) < 2: continue
        st = np.vstack([M[n][i] for n in have]); sg = np.sign(st)
        un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
        ag[i] = np.where(un, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)
        tot, cnt = 0.0, 0
        for a in range(len(have)):
            for b in range(a + 1, len(have)):
                x, y = M[have[a]][i], M[have[b]][i]
                tot += float(x @ y) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12); cnt += 1
        dc[i] = tot / cnt
    z = (rankdata(dc, method="average") - 0.5) / P; z[ndon < 2] = 0.0
    z = np.maximum(z, 0.3)
    print(f"\n##### {nm}  (builder's exponential form, gamma_pert=1.0, floor 0.3)")
    print(f"  {'gamma_gene':<12}{'mean alpha':>12}{'kappa':>10}{'M':>8}{'pds':>10}")
    for gam in (0.0, 0.05, 0.25, 1.0):
        A = np.ones((P, G)) if gam == 0 else np.vstack(
            [np.ones(G) if np.unique(ag[i]).size < 2
             else np.exp(-gam * desc_midrank(ag[i]) / 2000.0) for i in range(P)])
        dh = (z if gam > 0 else np.ones(P))[:, None] * A * dh0
        X = float((np.einsum("ij,ij->i", dh, t["d"]) * w).sum())
        B = float(_across_pert_budget(dh, cols, w))
        kap = X / np.sqrt(B * sD)
        pm = np.tile(bctrl, (P + 1, 1)); pm[:P] += dh
        pv = discrimination_score(pred_bulk=(np.append(tg, CTRL), pm), real_bulk=(rp, rmeans),
                                  pert_col="target_gene", control=CTRL, distance="cosine",
                                  exclude_target_gene=True, exclusion_scope="panel",
                                  control_source="pred", genes=rg)
        pds = float(np.mean(list(pv.values()))) if isinstance(pv, dict) else float(pv)
        if gam == 0: k0 = kap
        print(f"  {gam:<12}{A.mean():>12.4f}{kap:>10.5f}{kap/k0:>8.3f}{pds:>10.4f}")
