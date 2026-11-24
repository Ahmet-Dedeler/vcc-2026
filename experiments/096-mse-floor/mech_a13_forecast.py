"""a13's raw `mse`, forecast from the SHIPPED a13 artifact and the corrected board solve.

Everything except the alignment `X` is measured on `submissions/a13-shrink.vcc` itself:
`A`, `B`, `Kc`, `C_hat_p`, the per-row cap margin and `claim`, all per context, with #348's
budget evaluated per context as the scorer does.

The one assumed quantity is the DIRECTION multiplier

    kappa = X_a13 / X_prod ,   X_prod = the alignment at VCC_ALPHA 1.0796, unshrunk

which no surface we hold measures on the challenge panel. It is reported as a curve, with the
three readings that bracket it named, rather than collapsed into one number.

Also settles the check A13-PREDICTION.md left open: is a13's `claim` still at the #247 cap
after the shrinkage, per row?

    uv run python mech_a13_forecast.py
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from mech_a12_audit import BOARD, CTX, cp_proxies, measure, panel, score_mse  # noqa: E402

W = f"{ROOT}/experiments/096-mse-floor/work"


def main():
    cps = cp_proxies()
    cp = {c: cps["ntc"][c][0] for c in CTX}
    arms = {}
    for arm in ("a8", "a12", "a13"):
        rows = measure(f"{W}/{arm}-stats.npz", cp)
        m = panel(rows)
        m["rows"] = rows
        arms[arm] = m

    a8, a12, a13 = arms["a8"], arms["a12"], arms["a13"]
    mu = float(np.sqrt(a13["A"] / a8["A"]))     # a13's amplitude relative to production
    print("===== a13 as SHIPPED, measured on the real challenge panel =====")
    print(f"{'ctx':>4}{'A_c':>10}{'B_c':>10}{'claim_c':>10}{'K_c':>10}{'Kc_c':>9}"
          f"{'rho_c':>8}{'min(cl,B)':>11}{'num_c':>10}{'C_hat mean':>12}{'C/C_p':>8}"
          f"{'cap binds':>11}")
    for r in a13["rows"]:
        rc = min(r["claim"], r["B"])
        print(f"{r['ctx']:>4}{r['A']:>10.6f}{r['B']:>10.6f}{r['claim']:>10.6f}{r['K']:>10.6f}"
              f"{r['Kc']:>9.6f}{min(1.0, r['B']/r['claim']):>8.4f}{rc:>11.6f}"
              f"{r['A']+r['Kc']-rc:>10.6f}{r['Chat'].mean():>12.3f}"
              f"{r['Chat'].mean()/r['Cp']:>8.3f}{r['cap_binds']:>7d}/{r['n']:<4d}")
    num13 = a13["A"] + a13["Kc"] - a13["rho_claim"]
    print(f"{'mean':>4}{a13['A']:>10.6f}{a13['B']:>10.6f}{a13['claim']:>10.6f}{a13['K']:>10.6f}"
          f"{a13['Kc']:>9.6f}{a13['rho_min']:>8.4f}{a13['rho_claim']:>11.6f}{num13:>10.6f}")
    print(f"\n  amplitude vs production (VCC_ALPHA 1.0796, unshrunk): "
          f"sqrt(A_a13/A_a8) = {mu:.4f}   (a12 sits at {np.sqrt(a12['A']/a8['A']):.4f})")
    print(f"  claim as a fraction of the #247 cap: {100*a13['claim']/a13['K']:.2f}% "
          f"(a12: {100*a12['claim']/a12['K']:.2f}%)")
    print("  A13-PREDICTION.md left this open ('the check is per-row and has not landed').")
    print("  It has now landed, on the shipped artifact.")

    print("\n===== raw mse vs the one assumed quantity, the direction multiplier kappa =====")
    Xp, sD = 0.019240, 0.229141        # LS on a8+a10+a11, the three arms whose `claim`
    #                                    is their OWN measured C_hat and so carries no C_p
    #                                    proxy risk. Residual <= 0.0001 raw on all three.
    print(f"  X_prod = {Xp:.6f}, sumD = {sD:.6f}   (LS on a8+a10+a11; residual <= 0.0001)")
    print("  a13's claim is at the #247 cap, so it inherits the C_p proxy band measured in")
    print("  mech_cp_bias.py: +-0.007 raw on top of everything below.")
    print(f"  {'kappa':>8}{'X_a13':>10}{'cos_a13':>9}{'raw mse':>10}{'score_mse':>11}"
          f"{'overall d':>11}   reading")
    # a12's own kappa is exactly its amplitude, because a12 only rescales.
    k12 = np.sqrt(a12["A"] / a8["A"])
    # the per-gene half at gamma 0.25 measured on the harnesses (work/a13_resolve.npz):
    z = np.load(f"{W}/a13_resolve.npz")
    i = int(np.flatnonzero(np.isclose(z["gam"], 0.25))[0])
    surv = float(np.sqrt(z["A"][i] / z["A"][0]))     # norm survival of the per-gene half
    ch1, cor = float(z["mh1"][i]) / surv, float(z["mor"][i]) / surv
    print(f"  (the per-gene half at gamma=0.25 raises the COSINE by {ch1:.3f}x on h1-test and "
          f"{cor:.3f}x on orion;\n   rule 1 takes the smaller. It says nothing about the "
          f"per-perturbation half or about the challenge panel.)")
    a12_raw_board = BOARD["a12"]
    # The 097 DE gate's non-mse net for a13 against a12, per surface, `reach` dropped as
    # unreadable (same-policy twins differ by 0.0207 / 0.0245 raw). Rule 1 takes the worse.
    DE = {"h1-test": +0.0044 - 0.0011 + 0.0026 + 0.0006,
          "orion": -0.0086 + 0.0057 - 0.0145 + 0.0000}
    de_worse = min(DE.values())
    sc12 = score_mse(a12_raw_board)
    cases = [
        (1.000, "no directional gain at all -- X scales with the norm"),
        (ch1, f"per-gene half only, h1-test ({ch1:.3f}x); the per-perturbation half assumed 1.0"),
        (ch1 * 1.035, "per-gene h1 x per-perturbation orion -- both rule-1 minima"),
        (ch1 * 1.144, "per-gene h1 x per-perturbation h1 -- the optimistic same-surface pair"),
        (0.8212 not in (0,) and 1.4902, "what A13-PREDICTION.md's published raw 0.8212 requires"),
        (cor, f"per-gene half only, orion ({cor:.3f}x)"),
    ]
    for c, note in sorted(cases):
        k = mu * c
        raw = 1.0 + (num13 - 2 * k * Xp) / sD
        cos = k * Xp / np.sqrt(a13["A"] * sD)
        d_mse = (score_mse(raw) - sc12) / 6.0
        print(f"  {k:>8.4f}{k*Xp:>10.6f}{cos:>9.5f}{raw:>10.4f}{score_mse(raw):>11.4f}"
              f"{d_mse:>+11.5f}   c={c:.3f}  {note}")
    print(f"\n  `overall d` above is the mse leg ONLY. The 097 DE gate's non-mse net for a13")
    print(f"  against a12 is {DE['h1-test']:+.4f} (h1-test) / {DE['orion']:+.4f} (orion), with")
    print(f"  `reach` dropped as unreadable; rule 1 takes {de_worse:+.4f} scaled = "
          f"{de_worse/6:+.5f} overall.")
    print(f"  {'cosine mult c':>14}{'raw mse':>10}{'overall':>10}{'vs a12':>10}")
    for c in (1.0, ch1, ch1 * 1.035, ch1 * 1.144, 1.4902, cor):
        k = mu * c
        raw = 1.0 + (num13 - 2 * k * Xp) / sD
        ov = BOARD["a12_overall"] + (score_mse(raw) - sc12) / 6.0 + de_worse / 6.0
        print(f"  {c:>14.3f}{raw:>10.4f}{ov:>10.5f}{ov-BOARD['a12_overall']:>+10.5f}")
    print("\n  For reference, the same table run on a12's own measured terms reproduces the")
    print(f"  board: kappa = {k12:.4f} -> raw "
          f"{1.0 + (a12['A']+a12['Kc']-a12['rho_claim'] - 2*k12*Xp)/sD:.4f} against "
          f"{a12_raw_board:.4f} measured.")

    print("\n===== what a13 needs on the OTHER members, at the board's own exchange rate =====")
    print("  Two board amplitude moves are now available and they agree in sign on `pds` and")
    print("  `nmae`, which the 096 model held invariant:")
    print(f"    {'member':>8}{'a11 (dalpha -0.2416)':>24}{'a12 (dalpha +0.1728)':>24}"
          f"{'slope /unit alpha':>20}")
    d11 = dict(pds=+0.003328, mse=+0.082041, nmae=+0.002644, fid=-0.003638,
               reach=-0.011375, jac=-0.001414)
    d12 = dict(pds=-0.002004, mse=+0.009122, nmae=-0.003008, fid=+0.000062,
               reach=+0.004569, jac=+0.002494)
    for k in ("pds", "nmae", "fid", "reach", "jac"):
        s1, s2 = d11[k] / -0.2416, d12[k] / 0.1728
        agree = "consistent" if s1 * s2 > 0 else "SIGN FLIP"
        print(f"    {k:>8}{d11[k]:>24.6f}{d12[k]:>24.6f}"
              f"{f'{s1:+.4f} / {s2:+.4f}':>20}  {agree}")


if __name__ == "__main__":
    main()
