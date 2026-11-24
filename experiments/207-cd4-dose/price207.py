"""207 Phase 3 (cheap pass) — price the whole dose grid on the priceable members.

`pds`   PLOLO, 4 folds, challenge frames, a23 render (162's pre-committed predictor).
        overall = mean_fold d pds_raw / 0.452 / 6 * 1.011.

`mse`   through 096/197's identity, NOT through a local `mse` (162 measured PLOLO's `mse` at
        2.5x worse than predicting zero).  A CD4-dose arm holds the emission law, the spread
        policy, `keep_top` and `min_cpm` fixed, so `claim`, `Kc` and `SigmaD` are unchanged
        and -- because `add_post_whiten_residual` restores the base's median row norm by one
        scalar -- `A` is unchanged too.  That is CHECKED (`median_row_norm_wcols`), not
        assumed.  The only term that moves is `X`:

            MSE_c - 1 = (A_c + Kc_c - min(claim_c, B_c) - 2X) / SigmaD,  X = cos*sqrt(SigmaD*A)
            d raw_mse / d cos = -2.0826   ->   +0.3653 overall per +1.0 cosine

`nmae`  187 section 5's gene-wise identity with `alpha` fixed: d raw_nmae = -2 W_hi d q_w,
        W_hi = 0.0153 (187's measured value; the board's implied 0.012 is more favourable, so
        this is the conservative choice).  The exact 188 identity is run separately at the
        selected doses by `six207.py`.

TRANSFER OF THE COSINE.  PLOLO's cosine is carried onto the board RELATIVELY.  The transfer
factor is not assumed: the a14 -> a15 step is a BOARD-MEASURED cosine move of +6.8% (197
section 1.2) and PLOLO reads +5.42% for the identical configuration change, so PLOLO
UNDERSTATES the board by x1.255 on this exact arm family.  Both the uncalibrated (x1.0,
conservative) and the calibrated (x1.255) prices are reported.

    uv run python price207.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")

SPAN_PDS, CALIB = 0.452, 1.011
SPAN_NMAE, W_HI = 0.5963, 0.0153
# 197 section 1 / 204: the withheld reference's own terms at our operating point.
SUM_D, KC, CLAIM = 0.22914103411323902, 0.0054102686242422674, 0.24614814267726884
A_TERM, B_TERM = 0.2331583274832182, 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
OURS_COS = 0.079249                       # 197: our board transfer cosine at a15/rc_hybrid
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089     # rule 13


def mse_at(c: float) -> float:
    q, eps = CLAIM / SUM_D, A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw: float) -> float:
    return float(min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0))


def main() -> None:
    by = {(r["fold"], r["arm"]): r
          for r in json.loads((WORK / "sweep207.json").read_text())}
    doses = sorted({float(a[1:]) for (f, a) in by if a.startswith("f") and a[1:].replace(
        ".", "").isdigit()})
    base = "f0"

    def m(arm, key, folds=FOLDS):
        return float(np.mean([by[(f, arm)][key] for f in folds]))

    cos_base = m(base, "cos_agg")
    cosm_base = m(base, "cos")
    # the transfer calibration, measured not assumed
    plolo_rel_a14a15 = m("a15", "cos_agg") / m("a14", "cos_agg") - 1.0
    k_cos = 0.0683 / plolo_rel_a14a15
    print(f"PLOLO a14->a15 relative cos_agg  {plolo_rel_a14a15:+.5f}   "
          f"board (197 §1.2) +0.06830   =>  transfer factor k_cos = {k_cos:.3f}")
    print(f"PLOLO cos_agg at f0 (CD4 OFF) = {cos_base:.6f};  board cosine at a15 = "
          f"{OURS_COS:.6f}\n")

    rows = []
    for v in doses:
        arm = f"f{v:g}"
        if (FOLDS[0], arm) not in by:
            continue
        d_pds = m(arm, "pds") - m(base, "pds")
        d_cos = m(arm, "cos_agg") - m(base, "cos_agg")
        d_qw = m(arm, "a15_sgn_w") - m(base, "a15_sgn_w")
        rel = d_cos / cos_base
        # board cosine implied by the arm, relative to the SHIPPED point f=0.10
        rel_ship = (m(arm, "cos_agg") - m("f0.1", "cos_agg")) / m("f0.1", "cos_agg")
        out = {}
        for tag, k in (("k1", 1.0), ("kcal", k_cos)):
            bc = OURS_COS * (1.0 + k * rel_ship)
            o_mse = (score_mse(mse_at(bc)) - score_mse(mse_at(OURS_COS))) / 6.0
            out[tag] = dict(board_cos=bc, mse=o_mse)
        d_pds_ship = m(arm, "pds") - m("f0.1", "pds")
        d_qw_ship = m(arm, "a15_sgn_w") - m("f0.1", "a15_sgn_w")
        o_pds = d_pds_ship / SPAN_PDS / 6.0 * CALIB
        o_nmae = (2.0 * W_HI * d_qw_ship) / SPAN_NMAE / 6.0
        rows.append(dict(arm=arm, f=v, cos_agg=m(arm, "cos_agg"), cos=m(arm, "cos"),
                         pds=m(arm, "pds"), q_w=m(arm, "a15_sgn_w"),
                         d_cos_vs0=d_cos, rel_vs0=rel, rel_vs_ship=rel_ship,
                         d_pds_ship=d_pds_ship, d_qw_ship=d_qw_ship,
                         o_pds=o_pds, o_nmae=o_nmae,
                         o_mse_k1=out["k1"]["mse"], o_mse_kcal=out["kcal"]["mse"],
                         board_cos_k1=out["k1"]["board_cos"],
                         board_cos_kcal=out["kcal"]["board_cos"],
                         total_k1=o_pds + o_nmae + out["k1"]["mse"],
                         total_kcal=o_pds + o_nmae + out["kcal"]["mse"],
                         rownorm=m(arm, "median_row_norm_wcols")))

    print(f"{'f':>7}{'cos_agg':>10}{'rel vs ship':>12}{'board cos':>11}{'pds':>10}"
          f"{'q_w':>9} | {'ovr pds':>9}{'ovr mse':>9}{'ovr nmae':>10}{'TOTAL':>10}"
          f"{'TOTAL cal':>11}")
    print("-" * 118)
    for r in rows:
        mark = "  <-- SHIPPED" if abs(r["f"] - 0.10) < 1e-9 else ""
        print(f"{r['f']:>7.3g}{r['cos_agg']:>10.6f}{r['rel_vs_ship']:>+12.4f}"
              f"{r['board_cos_k1']:>11.6f}{r['pds']:>10.6f}{r['q_w']:>9.6f} | "
              f"{r['o_pds']:>+9.5f}{r['o_mse_k1']:>+9.5f}{r['o_nmae']:>+10.6f}"
              f"{r['total_k1']:>+10.5f}{r['total_kcal']:>+11.5f}{mark}")

    best_k1 = max(rows, key=lambda r: r["total_k1"])
    best_cal = max(rows, key=lambda r: r["total_kcal"])
    print(f"\nbest priced total (uncalibrated k=1.0):   f = {best_k1['f']:g}  "
          f"{best_k1['total_k1']:+.5f}   vs rule-13 gates {GATE_SHAPE:+.4f} / "
          f"{GATE_AMP:+.4f}")
    print(f"best priced total (calibrated  k={k_cos:.3f}):  f = {best_cal['f']:g}  "
          f"{best_cal['total_kcal']:+.5f}")
    (WORK / "price207.json").write_text(json.dumps(
        dict(k_cos=k_cos, plolo_rel_a14a15=plolo_rel_a14a15, cos_base=cos_base,
             cosm_base=cosm_base, rows=rows), indent=2))
    print(f"wrote {WORK/'price207.json'}")


if __name__ == "__main__":
    main()
