"""Score every arm through the EXACT identity, not the sqrt(W) shortcut.

miller (exp096 builder) is right and it changes the accounting.  Writing the metric at its own
optimum over the global amplitude s,

    MSE(s) - 1 = [ s^2 A - 2 s X + Kc - min(claim, s^2 B) ] / sumD ,   B = _across_pert_budget

the minimum sits on the rho cliff s_c = sqrt(claim/B) and evaluates to

    MSE_min = 1 + [ Kc + claim*(A-B)/B - 2*sqrt(claim/B)*X ] / sumD
            = 1 + [ Kc + claim*eps - 2*sqrt(claim*sumD)*kappa ] / sumD

with the two DIMENSIONLESS shape quantities

    kappa = X / sqrt(B*sumD)        alignment, budget-referenced   (096's cosine, with B for A)
    eps   = (A - B)/B               the common-mode share the sqrt(W)=X^2/A shortcut DROPS

`MSE_MECHANISM.md`'s master formula keeps only kappa.  That is fine for families that leave eps
alone (a global scalar: eps is invariant) and WRONG for any family that moves the common mode
relative to the differential -- which is exactly Miller's shrink-toward-the-panel-mean.  And the
penalty is worse on the board than on either harness, because it is multiplied by `claim` and
claim/B is 0.18 (h1) / 0.37 (orion) but 0.88 on the challenge at the #247 cap.

So each arm is reported as (M = kappa ratio, eps) and then mapped onto the challenge panel with
the board's own solved constants.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
from vcc2026.donor_transfer import Donor
from cell_eval2.metrics.delta import _across_pert_budget

# challenge constants, all from MSE-MECHANISM.md's board-anchored solve
CH = dict(claim=0.24653, Kc=0.005410, sumD=0.228572, kappa0=0.07402, B=0.280503)
B_ANCHOR, SPAN = 0.987863, 0.950163

HARNESS = [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
]
DONORS = ("k562_gwps", "orion_hek293t")


def prep(real, pred):
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    tg = np.array([str(x) for x in pp[t["keep"]]])
    P, G = t["dh"].shape
    Ms, hv = {}, {}
    for nmd in DONORS:
        dn = Donor.load(f"{ROOT}/data/donors/{nmd}")
        m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in rg]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        M = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(tg):
            if x in pi:
                M[r, ok] = s[pi[x]][cols[ok]]; h[r] = True
        Ms[nmd] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nmd] = h
        del dn, s
    A_, B_ = Ms[DONORS[0]], Ms[DONORS[1]]
    both = hv[DONORS[0]] & hv[DONORS[1]]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    t["dcos"] = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na * nb, 1e-12), 0.0)
    ag = np.sqrt(np.abs(A_) * np.abs(B_)) * ((np.sign(A_) == np.sign(B_)) & (A_ != 0) & (B_ != 0))
    al = np.empty_like(ag)
    for i in range(P):
        al[i] = (np.argsort(np.argsort(ag[i])).astype(float) + 0.5) / G
    t["al"] = al
    t["targets"] = tg
    return t


def shape(DH, t):
    """(kappa, eps, X, A, B) for a candidate effect matrix."""
    w, d = t["w"], t["d"]
    X = float((np.einsum("ij,ij->i", DH, d) * w).sum())
    A = float((np.einsum("ij,ij->i", DH, DH) * w).sum())
    B = float(_across_pert_budget(DH, t["cols"], w))
    return X / np.sqrt(B * t["sumD"]), (A - B) / B, X, A, B


def rank01(s):
    return (np.argsort(np.argsort(s)).astype(float) + 0.5) / len(s)


def arms(t):
    dh = t["dh"]
    z = rank01(t["dcos"]); al = t["al"]
    Sw = float(t["w"].sum())
    dhbar = (t["w"][:, None] * dh).sum(0) / Sw
    out = [("base (as shipped)", dh)]
    for g in (1.0, 2.0, 3.0):
        out.append((f"per-PERT donor_cos^{g:g}", (z ** g)[:, None] * dh))
    for e in (0.5, 1.0, 2.0):
        out.append((f"per-GENE agreement^{e:g} (to ZERO)", (al ** e) * dh))
    for e in (0.5, 1.0):
        out.append((f"per-GENE agreement^{e:g} (to MEAN)",
                    (al ** e) * dh + (1 - al ** e) * dhbar[None, :]))
    for g, e in ((2.0, 1.0), (2.0, 2.0), (3.0, 1.0)):
        out.append((f"STACK pert^{g:g} x gene^{e:g}", (z ** g)[:, None] * (al ** e) * dh))
    # oracles, for the ceiling column
    x = np.einsum("ij,ij->i", dh, t["d"]); a = np.einsum("ij,ij->i", dh, dh)
    out.append(("ORACLE per-pert z=x/a", (x / np.maximum(a, 1e-30))[:, None] * dh))
    U = (t["w"][:, None] * dh * t["d"]).sum(0); Q = (t["w"][:, None] * dh * dh).sum(0)
    sg = np.where(Q > 0, U / np.maximum(Q, 1e-30), 0.0)
    out.append(("ORACLE per-gene s_g", sg[None, :] * dh))
    return out


def challenge_raw(M, eps):
    c, sD = CH["claim"], CH["sumD"]
    return 1 + (CH["Kc"] + c * eps - 2 * np.sqrt(c * sD) * CH["kappa0"] * M) / sD


res = {}
for nm, r, p in HARNESS:
    t = prep(r, p)
    k0, e0, X0, A0, B0 = shape(t["dh"], t)
    print(f"\n########## {nm}  n={t['n']}  claim/B = {t['claim']/B0:.3f} "
          f"(challenge at the cap: {CH['claim']/CH['B']:.3f}) ##########")
    print(f"  base kappa={k0:.5f}  eps=(A-B)/B={e0:+.5f}   "
          f"[096 reported deconfounded cos X/sqrt(A*sumD) = {X0/np.sqrt(A0*t['sumD']):.5f}]")
    print(f"\n  {'arm':<38}{'kappa':>9}{'M':>8}{'eps':>9}{'sqrtW/A shortcut':>18}")
    res[nm] = {}
    for lbl, DH in arms(t):
        k, e, X, A, B = shape(DH, t)
        short = np.sqrt((X * X / A) / (X0 * X0 / A0))
        res[nm][lbl] = (k / k0, e)
        print(f"  {lbl:<38}{k:>9.5f}{k/k0:>8.3f}{e:>+9.5f}{short:>17.3f}x")

print("\n\n########## mapped onto the CHALLENGE panel at the #247 cap "
      "(a12's emission policy) ##########")
print("  transferable M = min(h1, orion);  eps = max(h1, orion) i.e. the pessimistic pair")
print(f"\n  {'arm':<38}{'M h1':>7}{'M ori':>7}{'M use':>7}{'eps use':>9}"
      f"{'raw mse':>9}{'score_mse':>10}{'overall':>9}{'d vs a12':>9}")
base_raw = challenge_raw(1.0, res["h1-test"]["base (as shipped)"][1])
for lbl in res["h1-test"]:
    m1, e1 = res["h1-test"][lbl]
    m2, e2 = res["orion"][lbl]
    M, E = min(m1, m2), max(e1, e2)
    raw = challenge_raw(M, E)
    sc = max(0.0, min(1.0, (B_ANCHOR - raw) / SPAN))
    b_sc = max(0.0, min(1.0, (B_ANCHOR - base_raw) / SPAN))
    print(f"  {lbl:<38}{m1:>7.3f}{m2:>7.3f}{M:>7.3f}{E:>+9.5f}"
          f"{raw:>9.4f}{sc:>10.4f}{sc/6:>9.4f}{(sc-b_sc)/6:>+9.4f}")
print(f"\n  base row reproduces MSE-MECHANISM's 0.8699 at the cap: {base_raw:.4f}")
