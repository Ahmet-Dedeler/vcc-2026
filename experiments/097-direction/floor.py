"""What does a FLOOR on the per-perturbation weight cost?

Two hazards both argue for a floor, and they are the same hazard seen twice:
  * `make_submission_pm_cd4.py:1285` raises on an empty row, and a row emitting the control
    profile calls nothing, which is an `n_pred` loss.
  * `pds_cosine` is SCALE-INVARIANT, so a graded weight is free for it in principle -- but the
    scorer's own docstring records that a ZERO-norm predicted effect ties its whole row at
    distance 1.0 and scores exactly 0.5 under the midrank policy.  A row whose effect falls
    below COUNT QUANTIZATION after emission is that row.  At pds raw 0.8333 and span 0.452,
    each such target costs 0.00041 overall, so 66 of them would cost -0.027 -- more than the
    whole mse gain.
Priced here on the mse side; the pds side is bounded above by the count above.
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
OUT = {}
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
            if x in pi:
                M[r, ok] = s[pi[x]][cols[ok]]; h[r] = True
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
    zr = (rankdata(dcos, method="average") - 0.5)/P; zr[~both] = 1.0

    def kap(DH):
        X = float((np.einsum("ij,ij->i", DH, d)*w).sum())
        B = float(_across_pert_budget(DH, t["cols"], w))
        return X/np.sqrt(B*t["sumD"])
    k0 = kap(dh)
    print(f"\n##### {nm}: floor on the FINAL per-perturbation weight z = max(rank01^g, f)")
    print(f"  {'':<12}" + "".join(f"{'f='+str(f):>10}" for f in (0.0, 0.1, 0.2, 0.3, 0.5, 1.0)))
    for g in (1.0, 2.0, 3.0):
        for lbl, mat in (("pert only", None), ("x gene^1", al)):
            row = []
            for f in (0.0, 0.1, 0.2, 0.3, 0.5, 1.0):
                z = np.maximum(zr ** g, f)
                DH = z[:, None]*dh if mat is None else z[:, None]*mat*dh
                row.append(kap(DH)/k0)
            nlow = int((np.maximum(zr**g, 0.0) < 0.05).sum())
            print(f"  g={g:<3.1f} {lbl:<7}" + "".join(f"{v:>10.3f}" for v in row)
                  + f"   ({nlow} rows below z=0.05 unfloored)")
