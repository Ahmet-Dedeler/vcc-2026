"""Decompose `expr_mse_unbiased_capped_norm` into named terms and verify the identity.

    N_p   = ( ||b_hat_p - b_p||^2 - rho*min(C_hat_p, C_p) - C_p ) / G_p
    D_p   = ( ||b_p - b_ctrl||^2  - C_p - C_ctrl )               / G_p
    MSE   = sum_p N_p / sum_p D_p

Write d_p = b_p - b_ctrl (real effect) and dh_p = b_hat_p - b_ctrl (predicted effect), both on
row p's OWN scored gene set (all genes bar p's target). Then
||b_hat_p - b_p||^2 = ||dh_p||^2 - 2<dh_p, d_p> + ||d_p||^2 and, with w_p = 1/G_p,

    A     = sum_p w_p ||dh_p||^2          predicted amplitude   (needs only pred + real ctrl)
    X     = sum_p w_p <dh_p, d_p>         alignment             (needs the reference)
    R     = sum_p w_p ||d_p||^2           real amplitude        (reference only)
    K     = sum_p w_p C_p                 reference correction
    Kc    = sum_p w_p C_ctrl              control correction
    claim = sum_p w_p min(C_hat_p, C_p)   what the submission claims
    rho   = min(1, B / claim)

    sum N = A - 2X + R - rho*claim - K          sum D = R - K - Kc

    ==> MSE - 1 = ( A - 2X + Kc - rho*claim ) / sum D

which is the identity everything in this experiment is read through.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from msekit import CTRL, TS, group_stats, score  # noqa: E402

from cell_eval2.metrics.delta import (  # noqa: E402
    PRED_TRACE_CAP_K, _across_pert_budget, _row_weights,
)
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402


def terms(pred_perts, pred_sums, pred_counts, pred_jk,
          real_perts, real_sums, real_counts, real_jk, genes):
    """Every named term above, plus the scorer's own MSE for cross-checking."""
    pm = bulk_lognorm_means(pred_sums.astype(np.float64), TS)
    rm = bulk_lognorm_means(real_sums.astype(np.float64), TS)
    gpos = {g: i for i, g in enumerate(genes)}
    rindex = {str(p): i for i, p in enumerate(real_perts)}
    ci = rindex[CTRL]
    G = genes.size
    bctrl = rm[ci]
    Cctrl = float(real_jk[ci])

    keep, rows, cols = [], [], {}
    for i, p in enumerate(pred_perts):
        p = str(p)
        if p == CTRL or p not in rindex:
            continue
        k = len(keep)
        keep.append(i)
        rows.append(rindex[p])
        if p in gpos:
            cols[k] = gpos[p]
    keep = np.asarray(keep, np.intp)
    rows = np.asarray(rows, np.intp)
    w = _row_weights(keep.size, G, cols)          # 1/G_p, exactly the scorer's

    dh = pm[keep] - bctrl[None, :]
    d = rm[rows] - bctrl[None, :]
    a = np.einsum("ij,ij->i", dh, dh)
    x = np.einsum("ij,ij->i", dh, d)
    r = np.einsum("ij,ij->i", d, d)
    for k, c in cols.items():                      # drop each row's own target coordinate
        a[k] -= dh[k, c] ** 2
        x[k] -= dh[k, c] * d[k, c]
        r[k] -= d[k, c] ** 2

    Chat = pred_jk[keep]
    Cp = real_jk[rows]
    capped = np.minimum(Chat, PRED_TRACE_CAP_K * Cp)
    claim = float((capped * w).sum())
    budget = _across_pert_budget(pm[keep], cols, w)
    rho = min(1.0, budget / claim) if claim > 0 else 1.0

    A = float((a * w).sum())
    X = float((x * w).sum())
    R = float((r * w).sum())
    K = float((Cp * w).sum())
    Kc = float((np.full(keep.size, Cctrl) * w).sum())
    sN = A - 2 * X + R - rho * claim - K
    sD = R - K - Kc
    return dict(n=keep.size, G=G, A=A, X=X, R=R, K=K, Kc=Kc, claim=claim,
                budget=budget, rho=rho, sumN=sN, sumD=sD, mse=sN / sD,
                # per-row pieces, for the per-gene / per-pert questions
                w=w, a=a, x=x, r=r, Chat=Chat, Cp=Cp, capped=capped,
                dh=dh, d=d, cols=cols, keep=keep, rows=rows, bctrl=bctrl,
                Cctrl=Cctrl, pm=pm, rm=rm)


def report(name, t, mse_ref=None):
    sD = t["sumD"]
    print(f"\n===== {name}  ({t['n']} perturbations, {t['G']} genes) =====")
    print(f"  MSE (from terms)          {t['mse']:.10f}")
    if mse_ref is not None:
        print(f"  MSE (cell_eval2 scorer)   {mse_ref:.10f}   identity residual "
              f"{t['mse'] - mse_ref:+.3e}")
    print(f"  sum N = {t['sumN']:.8g}      sum D = {sD:.8g}")
    print("  --- terms, and each as a fraction of sum D (i.e. in MSE units) ---")
    for k in ("A", "X", "R", "K", "Kc", "claim", "budget"):
        print(f"    {k:<7}{t[k]:>14.8g}   / sumD = {t[k]/sD:>+10.5f}")
    print(f"    rho    {t['rho']:>14.6f}   (budget/claim = {t['budget']/t['claim']:.4f})")
    print("  --- the identity  MSE - 1 = (A - 2X + Kc - rho*claim)/sumD ---")
    for lbl, v in (("+A      (predicted amplitude)", t["A"]),
                   ("-2X     (alignment credit)", -2 * t["X"]),
                   ("+Kc     (control correction)", t["Kc"]),
                   ("-rho*claim (correction claimed)", -t["rho"] * t["claim"])):
        print(f"    {lbl:<34}{v/sD:>+10.5f}")
    print(f"    {'= MSE - 1':<34}{t['mse']-1:>+10.5f}")
    # alignment geometry
    cos = t["X"] / np.sqrt(t["A"] * t["R"])
    print(f"  --- geometry ---")
    print(f"    aggregate cosine  X/sqrt(A*R)      {cos:.5f}")
    print(f"    amplitude ratio   sqrt(A/R)        {np.sqrt(t['A']/t['R']):.5f}")
    print(f"    optimal global alpha*  X/A         {t['X']/t['A']:.5f}")
    gain = (t["A"] - t["X"]) ** 2 / t["A"] / sD
    print(f"    MSE at alpha* (corrections frozen) {t['mse']-gain:.5f}   (gain {gain:.5f})")


if __name__ == "__main__":
    real, pred = sys.argv[1], sys.argv[2]
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    assert np.array_equal(rg, pg)
    m, _, _ = score(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    report(os.path.basename(real), t, m)
