"""Global vs per-gene vs per-perturbation shrinkage, and what survives estimation noise.

At the optimum the metric collapses to one scalar. Writing the shrunk prediction as `s . dh_p`
(elementwise per gene, or a scalar per perturbation), with `U_g = sum_p w_p dh_pg^2`,
`V_g = sum_p w_p dh_pg d_pg`, and `B(s) ~= sum_g s_g^2 U_g` (the common mode is <1 % of A here):

    numerator(s) = sum_g s_g^2 U_g - 2 sum_g s_g V_g + Kc - min(B(s), claim)

Below the rho cliff the first and last terms cancel to the common mode, so the problem is
    maximise  sum_g s_g V_g   subject to   sum_g s_g^2 U_g = claim
whose solution is `s_g = lambda V_g / U_g` and whose value is `sqrt(claim * W)` with

    W = sum_g V_g^2 / U_g        per-gene family
    W = X^2 / A                  global-scalar family
    W = sum_p w_p <dh_p,d_p>^2 / ||dh_p||^2      per-perturbation family

    ==>  MSE_min  =  1 + ( Kc - 2 sqrt(claim * W) ) / sumD

so the whole shrinkage question reduces to how much `sqrt(W)` a family can buy -- and W is an
ORACLE quantity, fitted against the withheld truth. The honest number is the CROSS-VALIDATED
one: fit `s` on one half of the perturbations, score it on the other.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from decompose import terms  # noqa: E402
from msekit import group_stats  # noqa: E402

REAL = f"{ROOT}/data/harness/h1-test.h5ad"
PRED = f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"


def W_of(dh, d, w, s):
    """Explained sum of squares actually delivered by shrinkage pattern `s` [G] or [P,1]."""
    sdh = s * dh
    num = float((np.einsum("ij,ij->i", sdh, d) * w).sum())
    den = float((np.einsum("ij,ij->i", sdh, sdh) * w).sum())
    return num * num / den if den > 0 else 0.0


def main():
    rp, rs, rc, rj, rg = group_stats(REAL)
    pp, ps, pc, pj, pg = group_stats(PRED)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    dh, d, w = t["dh"], t["d"], t["w"]
    # zero each row's own target coordinate so every sum below matches the scored gene set
    for k, c in t["cols"].items():
        dh[k, c] = 0.0
        d[k, c] = 0.0
    A, X, sD, Kc = t["A"], t["X"], t["sumD"], t["Kc"]
    P, G = dh.shape
    U = np.einsum("ij,ij->j", dh * w[:, None], dh)
    V = np.einsum("ij,ij->j", dh * w[:, None], d)
    ok = U > 0

    W_glob = X * X / A
    W_gene = float((V[ok] ** 2 / U[ok]).sum())
    nr = np.einsum("ij,ij->i", dh, dh)
    xp = np.einsum("ij,ij->i", dh, d)
    W_pert = float((w[nr > 0] * xp[nr > 0] ** 2 / nr[nr > 0]).sum())

    print(f"h1-test: {P} perturbations, {G} genes,  sumD={sD:.6f}  Kc={Kc:.6f}  A={A:.6f}")
    print(f"\n  ORACLE explained sum of squares W (fitted on the truth -- an upper bound only)")
    for name, W in (("global scalar", W_glob), ("per-gene", W_gene),
                    ("per-perturbation", W_pert)):
        print(f"    {name:<20}{W:>12.6g}   sqrt(W/W_global) = {np.sqrt(W/W_glob):>7.3f}x")

    print(f"\n  CROSS-VALIDATED (fit on half the perturbations, evaluate on the other half)")
    rng = np.random.default_rng(0)
    res = {k: [] for k in ("global", "gene", "gene_JS", "pert_oracle")}
    for rep in range(20):
        perm = rng.permutation(P)
        fit, ev = perm[: P // 2], perm[P // 2:]
        Uf = np.einsum("ij,ij->j", dh[fit] * w[fit, None], dh[fit])
        Vf = np.einsum("ij,ij->j", dh[fit] * w[fit, None], d[fit])
        Af = float((np.einsum("ij,ij->i", dh[fit], dh[fit]) * w[fit]).sum())
        Xf = float((np.einsum("ij,ij->i", dh[fit], d[fit]) * w[fit]).sum())
        g_glob = Xf / Af
        okf = Uf > 0
        s_gene = np.zeros(G)
        s_gene[okf] = Vf[okf] / Uf[okf]
        # James-Stein / ridge shrink of the per-gene coefficients toward the global one.
        # tau^2 = between-gene variance of the TRUE coefficients, estimated by subtracting the
        # sampling variance of s_gene from its observed spread; lam_g = tau^2/(tau^2+v_g).
        resid = d[fit] - s_gene[None, :] * dh[fit]
        v_g = np.zeros(G)
        v_g[okf] = np.einsum("ij,ij->j", resid[:, okf] * w[fit, None], resid[:, okf]) / Uf[okf] \
            / max(fit.size - 1, 1)
        spread = float(np.average((s_gene[okf] - g_glob) ** 2, weights=Uf[okf]))
        tau2 = max(spread - float(np.average(v_g[okf], weights=Uf[okf])), 0.0)
        lam = np.zeros(G)
        lam[okf] = tau2 / (tau2 + v_g[okf]) if tau2 > 0 else 0.0
        s_js = g_glob + lam * (s_gene - g_glob)
        res["global"].append(W_of(dh[ev], d[ev], w[ev], np.full(G, g_glob)))
        res["gene"].append(W_of(dh[ev], d[ev], w[ev], s_gene))
        res["gene_JS"].append(W_of(dh[ev], d[ev], w[ev], s_js))
        nre = np.einsum("ij,ij->i", dh[ev], dh[ev])
        xpe = np.einsum("ij,ij->i", dh[ev], d[ev])
        res["pert_oracle"].append(
            float((w[ev][nre > 0] * xpe[nre > 0] ** 2 / nre[nre > 0]).sum()))
    base = float(np.mean(res["global"]))
    for k in ("global", "gene", "gene_JS", "pert_oracle"):
        v = np.asarray(res[k])
        print(f"    {k:<20}W = {v.mean():>11.6g} +- {v.std():>9.3g}   "
              f"sqrt(W/W_global) = {np.sqrt(v.mean()/base):>6.3f}x")
    print("    (`pert_oracle` cannot cross-validate -- a per-perturbation scalar has no\n"
          "     out-of-sample form when the held-out unit IS the perturbation; it is shown as\n"
          "     the oracle bound on the evaluation half only.)")

    print(f"\n  What each W is worth in raw mse on the CHALLENGE panel"
          f" (claim at the #247 cap):")
    bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
    sDc, Kcc, Kcap, Ac, Xc = (float(bs["sumD"]), float(bs["v11_Kc"]), float(bs["K"]),
                              float(bs["v11_A"]), float(bs["X"]))
    Wc = Xc * Xc / Ac
    for name, mult in (("global (measured)", 1.0),
                       ("per-gene, cross-validated", np.sqrt(np.mean(res["gene_JS"]) / base)),
                       ("per-gene, oracle", np.sqrt(W_gene / W_glob))):
        W = Wc * mult ** 2
        for cl, cn in ((float(bs["v11_claim"]), "phi=0"), (float(bs["a8_claim"]), "phi=4"),
                       (Kcap, "at the cap")):
            v = 1 + (Kcc - 2 * np.sqrt(cl * W)) / sDc
            print(f"    {name:<28}{cn:<12}raw mse {v:>7.4f}   "
                  f"score {min(max((0.987862879-v)/0.950162566,0),1):>6.4f}")


if __name__ == "__main__":
    main()
