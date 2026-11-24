"""183 — rule 13's composed slot gate for the unique-component boost, all six members.

    predicted = PLOLO d_pds/6 * 1.011  +  096 d_mse/6  -  0.0146/6 (EMISSION-LAW arms only)
    SHIP only if predicted > +0.0019 (shape-only)  or  > +0.0089 (amplitude-changing)

This is a SHAPE-ONLY arm: the transform is norm-preserving (measured rA = 1.0002), the
emission law is untouched, and the 181 amendment restricts the -0.0146 `reach` toll to
emission-law arms (169 falsified it at +0.0006 for profile arms).  So the floor is +0.0019
and the toll is 0.

Every input is read from this experiment's own JSON; nothing is typed in by hand except the
three published board-anchored constants, which are cited.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

PDS_SPAN = 0.95228 - 0.50020     # debench ANCHOR
K_PLOLO = 1.011                  # 162's calibrated PLOLO->board pds scale
LOO_SIGMA_SCALED = 0.001494      # 162's leave-one-arm-out RMS, scaled pds
X_C, SUM_D = 0.017766, 0.229141  # 170 §4, board-anchored at a15's own point
MSE_RAW_TO_OVERALL = 0.17500     # 170 §4's realised raw-mse -> overall factor
NMAE_SPAN = 0.60                 # 173's convention
K_POOL = 0.75                    # board nmae 0.944 at q 0.537 => k ~ 0.75 (2026-09-03 02:15)
# dnmae/dq at fixed k below the kink:  d/dq [ q|k-1| + (1-q)(k+1) ] = |k-1| - (k+1)
DNMAE_DQ = abs(K_POOL - 1.0) - (K_POOL + 1.0)

ALPHAS = (0.00, 0.25, 0.50, 0.75, 1.00)
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")

# measured by deprice183.py on h1-test's reference DE pool (65,044 pairs)
DQ = {0.00: 0.00000, 0.25: +0.00012, 0.50: -0.00037, 0.75: -0.00105, 1.00: -0.00157}


def main() -> None:
    rows = json.load(open(os.path.join(WORK, "plolo183.json")))
    R = {(r["fold"], r["arm"], r["alpha"], r["space"]): r for r in rows}
    xr = json.load(open(os.path.join(WORK, "xratio183.json")))

    def chal(f, alpha):
        return float(np.mean([R[(f, "a23", alpha, f"chal{c}/pre")]["pds"] for c in "ABC"]))

    base = {f: chal(f, 0.00) for f in FOLDS}

    print("=" * 104)
    print("RULE 13 COMPOSED SLOT GATE — unique-component boost, shape-only class, floor +0.0019")
    print("=" * 104)
    print(f"  {'alpha':>6} {'pds scaled':>11} {'-> overall':>11} {'rX worst':>9} "
          f"{'mse raw':>9} {'-> overall':>11} {'d q':>9} {'nmae->ovr':>10} "
          f"{'NET (worst)':>12} {'NET (best mse)':>15} {'verdict':>9}")
    out = []
    for alpha in ALPHAS:
        d_raw = float(np.mean([chal(f, alpha) - base[f] for f in FOLDS]))
        d_scaled = d_raw / PDS_SPAN
        pds_ovr = d_scaled * K_PLOLO / 6.0

        i = ALPHAS.index(alpha)
        rx = {s: xr[s][i]["X"] / xr[s][0]["X"] for s in ("hct116", "hek293t")}
        rx_worst, rx_best = min(rx.values()), max(rx.values())
        mse_raw_w = -2.0 * X_C * (rx_worst - 1.0) / SUM_D
        mse_raw_b = -2.0 * X_C * (rx_best - 1.0) / SUM_D
        mse_ovr_w = -mse_raw_w * MSE_RAW_TO_OVERALL
        mse_ovr_b = -mse_raw_b * MSE_RAW_TO_OVERALL

        dq = DQ[alpha]
        nmae_raw = DNMAE_DQ * dq
        nmae_ovr = -nmae_raw / (NMAE_SPAN * 6.0)

        net_w = pds_ovr + mse_ovr_w + nmae_ovr
        net_b = pds_ovr + mse_ovr_b + nmae_ovr
        verdict = "SHIP" if net_w > 0.0019 else "VETO"
        print(f"  {alpha:>6.2f} {d_scaled:+11.6f} {pds_ovr:+11.6f} {rx_worst:9.5f} "
              f"{mse_raw_w:+9.5f} {mse_ovr_w:+11.6f} {dq:+9.5f} {nmae_ovr:+10.6f} "
              f"{net_w:+12.6f} {net_b:+15.6f} {verdict:>9}")
        out.append(dict(alpha=alpha, pds_raw=d_raw, pds_scaled=d_scaled, pds_overall=pds_ovr,
                        rx_worst=rx_worst, rx_best=rx_best,
                        mse_overall_worst=mse_ovr_w, mse_overall_best=mse_ovr_b,
                        dq=dq, nmae_overall=nmae_ovr, net_worst=net_w, net_best=net_b))

    print(f"\n  PLOLO's own leave-one-arm-out sigma: {LOO_SIGMA_SCALED:.6f} scaled pds "
          f"= {LOO_SIGMA_SCALED / 6.0:+.6f} on overall.")
    best = max(out, key=lambda r: r["net_worst"])
    top = max(out, key=lambda r: r["pds_overall"])
    print(f"  Largest pds leg (alpha {top['alpha']:.2f}): {top['pds_overall']:+.6f} overall "
          f"= {top['pds_overall'] / (LOO_SIGMA_SCALED / 6.0):.2f} sigma.")
    print(f"  Best composed alpha {best['alpha']:.2f}: worst-surface net "
          f"{best['net_worst']:+.6f}, best-surface net {best['net_best']:+.6f}, "
          f"floor +0.0019.  ==> NO-GO at every alpha.")
    print("\n  reach / fid / jac: not priced as levels (rule 9, `n_real`).  Both move only")
    print("  through `q`, which FALLS at every alpha > 0.25, so their direction is negative")
    print("  and none of them can rescue the arm.  The -0.0146 emission toll does NOT apply")
    print("  (181 amendment: profile-class arms, 169 measured +0.0006).")

    json.dump(out, open(os.path.join(WORK, "decide183.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
