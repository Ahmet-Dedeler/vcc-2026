"""209 Phase 3 — the amplitude term and the `VCC_ALPHA` re-solve.

`analyse209.py` prices `mse` through 197's identity holding the delivered amplitude `A`
fixed, which is 204 §6's convention.  A deflation does not hold `A` fixed exactly: the
norm restoration acts on the DONOR's own median row norm, while the delivered row norm is
measured after the operator, the CD4 residual and the render.  Measured, the delivered
median row norm moves by up to 1.5% on the norm-restored arms and 2.7% on the `_raw` one.

This script prices that second-order term explicitly, both ways:

  (i)  `VCC_ALPHA` RE-SOLVED, alpha -> 0.982 / s, restoring the delivered amplitude
       exactly.  Then the arm is pure shape and `analyse209.py`'s cosine-only price is the
       whole `mse` story.  187 measured the global-alpha optimum at a15's 0.982 to four
       decimals FOR THE SHIPPED SHAPE, and it is degree-(-1) homogeneous in the delivered
       amplitude, so the re-solve is a division and costs nothing to find.
  (ii) `VCC_ALPHA` LEFT AT 0.982.  Then A -> s^2 A, B -> s^2 B, X -> s X in
       MSE - 1 = (A + Kc - min(claim, B) - 2X) / SumD, and because a15 sits on the `B` leg
       in two of three contexts (rho 0.992 / 0.940 / 0.881), shrinking B shrinks the
       credited jackknife term.  That is a real extra cost and it is charged here.

SATURATION.  `claim` is a functional of the EMITTED counts' jackknife noise.  This arm class
touches neither the emission law, the spread policy, `keep_top`, nor the per-cell library
sizes, so `claim` is unchanged; what moves is `B`, the prediction's across-perturbation
spread.  We ship `claim` at 100.0% of `PRED_TRACE_CAP_K` (197 §1), so there is no headroom
above and the only direction saturation can move is down, through `B`.  Both branches are
reported.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

SUM_D = 0.22914103411323902
KC = 0.0054102686242422674
CLAIM = 0.24614814267726884
A_TERM = 0.2331583274832182
B_TERM = 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
OURS_COS = 0.079249
ALPHA_SHIPPED = 0.982
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")


def raw_general(s: float, cos: float) -> float:
    """MSE - 1 = (A + Kc - min(claim, B) - 2X)/SumD with A->s^2 A, B->s^2 B, X->s X."""
    A, B = s * s * A_TERM, s * s * B_TERM
    X = s * cos * np.sqrt(A_TERM * SUM_D)
    return 1.0 + (A + KC - min(CLAIM, B) - 2.0 * X) / SUM_D


def score_mse(raw: float) -> float:
    return min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0)


def main() -> None:
    a = json.load(open(WORK / "analyse209.json"))
    base_raw = raw_general(1.0, OURS_COS)
    base_score = score_mse(base_raw)
    print(f"baseline: raw_general(1, {OURS_COS}) = {base_raw:.6f}  "
          f"(197's per-context board value 0.873805; the {base_raw - 0.873805:+.6f} offset "
          f"is the aggregate-vs-per-context approximation and CANCELS in every delta)")
    print(f"a15 rho per context (178 §3): 0.992 / 0.940 / 0.881  ->  on the `B` leg in "
          f"2 of 3, so shrinking B is charged\n")

    hdr = (f"{'arm':18s} {'s (delivered)':>13s} {'alpha re-solve':>14s} "
           f"{'Δcos(board)':>12s} | {'Δmse ovr, alpha re-solved':>26s} "
           f"{'Δmse ovr, alpha FIXED':>22s} {'extra cost of not re-solving':>30s}")
    print(hdr)
    print("-" * len(hdr))
    out = {}
    keys = [k for k in a["arms"] if k != "equal"]
    keys.sort(key=lambda k: -a["arms"][k]["d_cos"])
    for k in ["equal"] + keys:
        r = a["arms"][k]
        s = float(np.mean([r["row_norm_ratio"][f] for f in FOLDS]))
        cos_arm = r["implied_board_cos_rel"]
        d_resolved = (score_mse(raw_general(1.0, cos_arm)) - base_score) / 6.0
        d_fixed = (score_mse(raw_general(s, cos_arm)) - base_score) / 6.0
        out[k] = dict(s_delivered=s, alpha_resolved=ALPHA_SHIPPED / s,
                      d_board_cos=cos_arm - OURS_COS,
                      d_mse_overall_alpha_resolved=d_resolved,
                      d_mse_overall_alpha_fixed=d_fixed,
                      extra_cost_of_not_resolving=d_fixed - d_resolved,
                      total_alpha_resolved=r["overall_pds"] + d_resolved + r["overall_nmae"],
                      total_alpha_fixed=r["overall_pds"] + d_fixed + r["overall_nmae"])
        print(f"{k:18s} {s:13.5f} {ALPHA_SHIPPED / s:14.5f} "
              f"{cos_arm - OURS_COS:+12.6f} | {d_resolved:+26.5f} {d_fixed:+22.5f} "
              f"{d_fixed - d_resolved:+30.5f}")

    print(f"\n{'arm':18s} {'ovr pds':>10s} {'ovr mse':>10s} {'ovr nmae':>10s} "
          f"{'TOTAL (alpha re-solved)':>24s} {'TOTAL (alpha fixed)':>21s} {'vs gate':>12s}")
    for k in ["equal"] + keys:
        r, o = a["arms"][k], out[k]
        gate = 0.0019 if abs(o["s_delivered"] - 1.0) <= 0.015 else 0.0089
        print(f"{k:18s} {r['overall_pds']:+10.5f} "
              f"{o['d_mse_overall_alpha_resolved']:+10.5f} {r['overall_nmae']:+10.5f} "
              f"{o['total_alpha_resolved']:+24.5f} {o['total_alpha_fixed']:+21.5f} "
              f"{o['total_alpha_resolved'] / gate:+8.2f}x{gate:.4f}")
    (WORK / "price209.json").write_text(json.dumps(out, indent=2))
    print("\nwrote work/price209.json")


if __name__ == "__main__":
    main()
