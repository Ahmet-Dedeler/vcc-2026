"""Solve the challenge panel's two hidden terms from two of our own board points, then falsify.

    MSE(s) - 1 = [ s^2 A  +  Kc  -  min(1, s^2 B / claim) * claim  -  2 s X ] / sumD

`A`, `B`, `claim`, `Kc` are MEASURED on the real challenge panel (see mech_challenge.py).
`X` (alignment with the withheld truth) and `sumD` (the reference's corrected effect size) are
the only unknowns. Two board raws at two different `s` on the SAME model determine both:

    v11  VCC_ALPHA 1.0796  ->  s = 1.00000  ->  raw 1.5801   (docs/STRATEGY.md addendum 5)
    v12  VCC_ALPHA 0.40    ->  s = 0.37050  ->  raw 0.9640   (leaderboard, "Alpha 0.40")

and a THIRD point, which the fit never sees, falsifies it:

    a8   VCC_ALPHA 1.0796, phi=4  ->  raw 1.3627   ("Wildflower", the current entry)

a8 differs from v11 only in its EMISSION (gamma-Poisson phi=4) plus the pruned whitener, so it
moves `claim` by +0.048 and `A` by -0.0007 and nothing else. If the frame is right it must land
on 1.3627 without being fitted to.
"""
from __future__ import annotations

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
CTX = ("A", "B", "C")
B_ANCHOR, SPAN = 0.987862879, 0.950162566     # 6-point fit on the board, max resid 2.1e-5
ALPHA_PROD = 1.0796


def measured(stats):
    """Context-mean A, B, claim, Kc for one submission -- all from mech_challenge.main."""
    import sys
    sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
    from mech_challenge import main
    cp = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy.npz").items())
    rows = main(stats, stats, cp=cp)
    return {k: float(np.mean([r[k] for r in rows])) for k in ("A", "B", "claim", "Kc")}, rows


def mse_of(s, m, X, sD):
    rho_claim = min(m["claim"], s * s * m["B"])          # min(1, s^2 B/claim) * claim
    return 1.0 + (s * s * m["A"] + m["Kc"] - rho_claim - 2 * s * X) / sD


def solve(m, pts):
    """Two (s, raw) points -> (X, sumD). Linear in both once the rho branch is fixed."""
    (s1, y1), (s2, y2) = pts

    def branch(s):
        """(const, coef_of_X) for the numerator, given which side of the rho cliff s is on."""
        rho_claim = min(m["claim"], s * s * m["B"])
        return s * s * m["A"] + m["Kc"] - rho_claim, -2 * s
    c1, k1 = branch(s1)
    c2, k2 = branch(s2)
    # (y-1)*sD = c + k*X  ->  solve for X and sD
    M = np.array([[k1, -(y1 - 1.0)], [k2, -(y2 - 1.0)]])
    v = np.array([-c1, -c2])
    X, sD = np.linalg.solve(M, v)
    return float(X), float(sD)


def curve(m, X, sD, label, grid=None, cp_cap=None):
    grid = grid if grid is not None else np.linspace(0.05, 1.4, 28)
    print(f"\n  {label}: raw mse vs amplitude")
    print(f"    {'s':>7}{'VCC_ALPHA':>11}{'rho':>7}{'raw mse':>10}{'score_mse':>11}{'overall':>10}")
    best = None
    for s in grid:
        v = mse_of(s, m, X, sD)
        sc = min(max((B_ANCHOR - v) / SPAN, 0.0), 1.0)
        rho = min(1.0, s * s * m["B"] / m["claim"])
        star = ""
        if best is None or v < best[1]:
            best = (s, v)
        print(f"    {s:>7.3f}{ALPHA_PROD*s:>11.4f}{rho:>7.3f}{v:>10.4f}{sc:>11.4f}"
              f"{sc/6:>10.5f}{star}")
    s_c = np.sqrt(m["claim"] / m["B"])
    v = mse_of(s_c, m, X, sD)
    sc = min(max((B_ANCHOR - v) / SPAN, 0.0), 1.0)
    print(f"    rho cliff at s = sqrt(claim/B) = {s_c:.4f}  (VCC_ALPHA {ALPHA_PROD*s_c:.4f})"
          f"  ->  raw {v:.4f}, score {sc:.4f}, overall {sc/6:+.5f}")
    return s_c, v
