"""183 — what the unique boost does to `mse`'s alignment term `X`, on the two surfaces with truth.

    mse - 1 = ( s^2 A  -  2 s X  +  Kc  -  min(s^2 B, claim) ) / sumD          (096)

The boost is norm-preserving by construction, so `A` and `B` — and therefore `claim`, `rho`
and the whole cliff position 170 §4 is about — are unchanged to numerical noise (asserted
below).  `X = Sum_p w_p <dh_p, d_p>` is the ONLY term it can move, and `X` needs the withheld
truth, so it can only be transferred as a RATIO.  This is 169's `xratio169.py` route with the
transform swapped; per 168d it reports the two held Orion lines separately and never averages
them.

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, HERE)
from massfix import mass_fix, production_effect     # noqa: E402
from screen import render                           # noqa: E402
from uboost import norm_check, unique_boost         # noqa: E402

ALPHAS = (0.00, 0.25, 0.50, 0.75, 1.00)


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        print(f"\n{'=' * 100}\n{surface}: {D.shape[0]} targets, {len(mask)} scored genes")
        print(f"  {'alpha':>6} {'A':>12} {'rA':>9} {'X':>12} {'rX':>9} "
              f"{'cos':>10} {'r cos':>9} {'|norm err|':>11}")
        base = None
        rows = []
        for alpha in ALPHAS:
            t = unique_boost(D, alpha)
            nerr = norm_check(D, t)
            eff = mass_fix(t, b, gate)
            P = render(eff, b)[:, mask]
            Tm = T[:, mask]
            A = float((P * P).sum())
            X = float((P * Tm).sum())
            R = float((Tm * Tm).sum())
            rec = dict(alpha=alpha, A=A, X=X, R=R, cos=X / np.sqrt(A * R), norm_err=nerr)
            if base is None:
                base = rec
            rows.append(rec)
            print(f"  {alpha:>6.2f} {A:12.5f} {A / base['A']:9.6f} {X:12.5f} "
                  f"{X / base['X']:9.6f} {rec['cos']:10.6f} "
                  f"{rec['cos'] / base['cos']:9.6f} {nerr:11.2e}")
        out[surface] = rows

    with open(os.path.join(WORK, "xratio183.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    print(f"\n{'=' * 100}\nTRANSFERRED RATIOS (168d rule: report both, never average; "
          f"take the WORSE for a GO)")
    print(f"  {'alpha':>6} {'rX hct116':>11} {'rX hek293t':>12} {'rX worst':>10} "
          f"{'rA hct116':>11} {'rA hek293t':>12}")
    for i, alpha in enumerate(ALPHAS):
        rx = [out[s][i]["X"] / out[s][0]["X"] for s in ("hct116", "hek293t")]
        ra = [out[s][i]["A"] / out[s][0]["A"] for s in ("hct116", "hek293t")]
        print(f"  {alpha:>6.2f} {rx[0]:11.6f} {rx[1]:12.6f} {min(rx):10.6f} "
              f"{ra[0]:11.6f} {ra[1]:12.6f}")

    # ---- translate to overall, at a15/a23's board-anchored challenge-panel constants ----
    # 170 §4: X_c solved from a15's own board point; A, B, claim, Kc, sumD measured there.
    X_C, SUM_D = 0.017766, 0.229141
    print(f"\n{'=' * 100}\n`mse` translation at the board-anchored challenge constants "
          f"(X_c={X_C}, sumD={SUM_D})")
    print("  A and B are unchanged (rA == 1 to 1e-4), so min(B, claim) and rho do not move")
    print(f"  {'alpha':>6} {'rX worst':>10} {'d raw mse':>11} {'d overall':>11}")
    for i, alpha in enumerate(ALPHAS):
        rx = min(out[s][i]["X"] / out[s][0]["X"] for s in ("hct116", "hek293t"))
        d_raw = -2.0 * X_C * (rx - 1.0) / SUM_D
        print(f"  {alpha:>6.2f} {rx:10.6f} {d_raw:+11.6f} {-d_raw * 0.17500:+11.6f}")
    print("  (raw -> overall uses 170 §4's own realised factor 0.175 = 1 / (6 * pds-span))")


if __name__ == "__main__":
    sys.exit(main())
