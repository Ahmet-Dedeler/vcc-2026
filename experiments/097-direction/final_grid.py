"""The final arm table, with the single-donor policy and the floor both settled by measurement.

singledonor.py established that h1-test's 4 rows with fewer than two donors carry 13.1 % of the
predicted amplitude A and +0.5 % of X (mean cos_p -0.0012 against +0.0369 for the rest), so they
must be DAMPED, not held at z = 1.  The policy used here:

    z_p = max( rank01(donor_cos_p) , floor ) ** gamma     with donor_cos_p = 0 when < 2 donors
          (a target with no cross-donor evidence therefore ranks low, but never below the floor)
    alpha_pg = midrank of the donor-agreement vector within row p, /G ; = 1 if the row has <2
          distinct values

Everything through the exact identity with B recomputed.
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
H = {}
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
    H[nm] = dict(dh=dh, d=d, w=w, cols=t["cols"], sumD=t["sumD"],
                 r=(rankdata(dcos, method="average") - 0.5)/P, al=al)

def shape(h, DH):
    X = float((np.einsum("ij,ij->i", DH, h["d"])*h["w"]).sum())
    A = float((np.einsum("ij,ij->i", DH, DH)*h["w"]).sum())
    B = float(_across_pert_budget(DH, h["cols"], h["w"]))
    return X/np.sqrt(B*h["sumD"]), (A-B)/B
for nm in H:
    H[nm]["k0"], H[nm]["e0"] = shape(H[nm], H[nm]["dh"])
    print(f"# {nm}: base kappa {H[nm]['k0']:.5f} eps {H[nm]['e0']:+.5f}")

print(f"\n{'arm':<40}{'M h1':>8}{'M ori':>8}{'min':>8}{'eps h1':>9}{'eps ori':>9}")
rows = []
for g in (0.5, 1.0, 2.0, 3.0):
    for f in (0.0, 0.2, 0.3, 0.5):
        for gene in (False, True):
            out = []
            for nm in ("h1-test", "orion"):
                h = H[nm]
                z = np.maximum(h["r"], f) ** g
                DH = z[:, None] * (h["al"] if gene else 1.0) * h["dh"]
                k, e = shape(h, DH)
                out.append((k/h["k0"], e))
            lbl = f"pert g={g} floor={f}" + (" x gene^1" if gene else "")
            rows.append((lbl, out[0][0], out[1][0], out[0][1], out[1][1]))
# gene-only reference row
out = []
for nm in ("h1-test", "orion"):
    h = H[nm]; k, e = shape(h, h["al"]*h["dh"]); out.append((k/h["k0"], e))
rows.append(("gene^1 ALONE (no pert weight)", out[0][0], out[1][0], out[0][1], out[1][1]))
rows.sort(key=lambda r: -min(r[1], r[2]))
for lbl, m1, m2, e1, e2 in rows[:16]:
    print(f"{lbl:<40}{m1:>8.3f}{m2:>8.3f}{min(m1,m2):>8.3f}{e1:>+9.5f}{e2:>+9.5f}")
print("  ... (best 16 of %d by the rule-1 minimum)" % len(rows))
