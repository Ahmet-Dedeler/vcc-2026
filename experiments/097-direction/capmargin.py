"""The right `claim` diagnostic: per-perturbation CAP MARGIN, not the panel total.

miller's correction, and it is right. Total `claim` sitting exactly at the #247 cap means
`min(C_hat_p, K*C_p)` is SATURATED -- claim is being set by the reference's jackknife, not by
ours -- so a flat panel total cannot falsify the concern. The failure mode is per-row: the
per-perturbation weight takes some rows to z = 0.09 (gamma 2, floor 0.3), an 11x cut on those
rows, and shrinkage flattens the profile as well as scaling it. If any row's C_hat_p falls under
K*C_p, that row leaves the cap, `claim` falls, and `s* = sqrt(claim/B')` moves under every
multiplier quoted anywhere in this thread.

So this reports, for the base and for the recommended arm, at each arm's own cliff amplitude:
  * how many rows are AT the cap (C_hat_p >= K*C_p),
  * the distribution of the margin C_hat_p / (K*C_p), especially its minimum,
  * the margin restricted to the LOW-z rows, which is where it is being spent.
"""
from __future__ import annotations
import sys, time
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.moments import jackknife_correction
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.metrics.delta import _across_pert_budget, PRED_TRACE_CAP_K
from decompose import terms
from msekit import CTRL, TS, group_stats
from vcc2026.donor_transfer import emit_counts, Donor
which = sys.argv[1]
REAL, PRED = {
 "h1": (f"{ROOT}/data/harness/h1-test.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
 "orion": (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
}[which]
DON = ("k562_gwps", "orion_hek293t")
rp, rs, rc, rj, rg = group_stats(REAL); pp, ps, pc, pj, pg = group_stats(PRED)
G = rg.size; ci = int(np.flatnonzero(rp == CTRL)[0])
bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
dh0 = bulk_lognorm_means(ps.astype(np.float64), TS) - bctrl[None, :]
S = ps.sum(1).astype(np.float64); NC = int(round(float(np.median(pc))))
t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
keep, rows, w, cols = t["keep"], t["rows"], t["w"], t["cols"]
Cp = rj[rows]                                   # the reference's own jackknife per scored row
tg = np.array([str(x) for x in pp[keep]]); P = keep.size
Ms, hv = {}, {}
for nd in DON:
    dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, sD_ = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cc = np.array([gi.get(str(g), -1) for g in rg]); ok = cc >= 0
    pi = {p: i for i, p in enumerate(dn.perts)}
    M = np.zeros((P, G)); h = np.zeros(P, bool)
    for r, xx in enumerate(tg):
        if xx in pi: M[r, ok] = sD_[pi[xx]][cc[ok]]; h[r] = True
    Ms[nd] = M / (float(np.median(np.linalg.norm(sD_[:, cc[ok]], axis=1))) or 1.0)
    hv[nd] = h; del dn, sD_
A_, B_ = Ms[DON[0]], Ms[DON[1]]; both = hv[DON[0]] & hv[DON[1]]
na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
dcos = np.where(both, np.einsum("ij,ij->i", A_, B_)/np.maximum(na*nb, 1e-12), 0.0)
ag = np.sqrt(np.abs(A_)*np.abs(B_)) * ((np.sign(A_)==np.sign(B_)) & (A_!=0) & (B_!=0))
al = np.empty((P, G))
for i in range(P):
    al[i] = rankdata(ag[i], method="average")/G
    if np.unique(ag[i]).size < 2: al[i] = 1.0
zk = np.maximum((rankdata(dcos, method="average") - 0.5)/P, 0.3) ** 2
Wfull = {"base": np.ones_like(dh0)}
Wm = np.ones_like(dh0); Wm[keep] = zk[:, None]*al; Wfull["STACK p2f3g1"] = Wm
print(f"# {which}: P={P} scored rows, NCELL={NC}, cap K={PRED_TRACE_CAP_K}", flush=True)
for tag, Wa in Wfull.items():
    dh = Wa * dh0
    B1 = float(_across_pert_budget(dh[keep], cols, w))
    s = float(np.sqrt(t["claim"]/B1))          # provisional cliff; claim re-read below
    for it in range(2):
        sums = np.zeros_like(ps); jks = np.zeros_like(pj)
        t0 = time.time()
        for i in range(pp.size):
            prof = np.expm1(bctrl + s*dh[i]); np.clip(prof, 0.0, None, out=prof)
            tot = prof.sum()
            if tot <= 0: continue
            prof = prof/tot
            v, col, ip = emit_counts(prof, np.full(NC, S[i]/NC), 0.0,
                                     np.random.default_rng(7919+i), mode="quasi",
                                     spread=0.5, overdispersion=4.0)
            X = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                           shape=(NC, G))
            sums[i] = np.asarray(X.sum(0)).ravel().astype(np.int64)
            jks[i] = jackknife_correction(X, np.zeros(NC, np.intp), 1, TS, chunk=64)[0]
        Chat = jks[keep]
        claim = float((np.minimum(Chat, PRED_TRACE_CAP_K*Cp)*w).sum())
        s = float(np.sqrt(claim/B1))
    marg = Chat/(PRED_TRACE_CAP_K*Cp)
    lowz = zk <= np.quantile(zk, 0.25)
    print(f"  {tag:<14} s={s:.4f}  claim={claim:.6f}  rows AT the cap "
          f"{int((marg>=1).sum())}/{P}   margin C_hat/(K*C_p): min {marg.min():.3f} "
          f"p05 {np.percentile(marg,5):.3f} median {np.median(marg):.2f}   "
          f"LOW-z quartile: min {marg[lowz].min():.3f} median {np.median(marg[lowz]):.2f}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
