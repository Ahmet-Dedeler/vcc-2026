"""216 — price the extended dose grid on the priceable members, 207's identities verbatim.

Nothing here is a new pricing model.  `price207.py`'s route is imported by value so that the
extended grid is priced on exactly the instrument 207's shipped operating point was chosen
with, and so its published rows are known-answer rows of this one:

    `pds`   PLOLO 4 folds, a23 render:  d raw pds / 0.452 / 6 * 1.011
    `mse`   096/197's identity through `X` only (the emission law, spread, `keep_top` and
            `min_cpm` are untouched and `add_post_whiten_residual` restores the base's median
            row norm, so `claim`, `Kc`, `SigmaD` and `A` are fixed)
    `nmae`  187 §5's gene-wise identity at fixed alpha, `d raw_nmae = -2 W_hi d q_w`

and the PLOLO->board cosine transfer `k_cos` is re-derived here from this experiment's own
a14 / a15 rows rather than copied.

    uv run python price216.py
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
SUM_D, KC, CLAIM = 0.22914103411323902, 0.0054102686242422674, 0.24614814267726884
A_TERM, B_TERM = 0.2331583274832182, 0.23104586978915895
B_ANCHOR, SPAN_MSE = 0.987862879, 0.950162566
OURS_COS = 0.079249
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
SIGMA_PDS, GATE1 = 0.001494, 1.645
LEADER_COS = 0.13211


def mse_at(c: float) -> float:
    q, eps = CLAIM / SUM_D, A_TERM / B_TERM - 1.0
    return 1.0 + KC / SUM_D + q * eps - 2.0 * c * np.sqrt(q * (1.0 + eps))


def score_mse(raw: float) -> float:
    return float(min(max((B_ANCHOR - raw) / SPAN_MSE, 0.0), 1.0))


def main() -> None:
    sweep = json.loads((WORK / "sweep216.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in sweep}
    gate = json.loads((WORK / "gate216.json").read_text())["arms"]

    def m(arm, key):
        return float(np.mean([by[(f, arm)][key] for f in FOLDS]))

    def mg(arm, key):
        return float(np.mean([gate[f"{f}/{arm}"][key] for f in FOLDS]))

    # k_cos, re-derived from this experiment's own a14 (fraction 0.05) and a15 (f 0.10) rows
    plolo_rel = mg("f010", "cos_agg") / mg("a14", "cos_agg") - 1.0
    k_cos = 0.0683 / plolo_rel
    ship = m("f0.1", "cos_agg")
    ship_pds = m("f0.1", "pds")
    ship_qw = m("f0.1", "a15_sgn_w")

    doses = sorted({r["f"] for r in sweep if r["kind"] == "cd4" and r["f"] > 0})
    rows = []
    for v in doses:
        arm = "f%g" % v
        rel_ship = (m(arm, "cos_agg") - ship) / ship
        d_pds = m(arm, "pds") - ship_pds
        d_qw = m(arm, "a15_sgn_w") - ship_qw
        o_pds = d_pds / SPAN_PDS / 6.0 * CALIB
        o_nmae = (2.0 * W_HI * d_qw) / SPAN_NMAE / 6.0
        out = {}
        for tag, k in (("k1", 1.0), ("kcal", k_cos)):
            bc = OURS_COS * (1.0 + k * rel_ship)
            out[tag] = dict(board_cos=bc,
                            mse=(score_mse(mse_at(bc)) - score_mse(mse_at(OURS_COS))) / 6.0)
        rows.append(dict(
            f=v, cos_agg=m(arm, "cos_agg"), pds=m(arm, "pds"), q_w=m(arm, "a15_sgn_w"),
            rel_vs_ship=rel_ship, d_pds=d_pds, pds_sigma=d_pds / SIGMA_PDS,
            o_pds=o_pds, o_nmae=o_nmae, o_mse_k1=out["k1"]["mse"],
            o_mse_kcal=out["kcal"]["mse"], board_cos_k1=out["k1"]["board_cos"],
            board_cos_kcal=out["kcal"]["board_cos"],
            total_k1=o_pds + o_nmae + out["k1"]["mse"],
            total_kcal=o_pds + o_nmae + out["kcal"]["mse"],
            within_gate1=bool(abs(d_pds / SIGMA_PDS) < GATE1)))

    print(f"PLOLO a14->a15 relative cos_agg {plolo_rel:+.5f}  =>  k_cos = {k_cos:.4f} "
          f"(207 published 1.2598)")
    print(f"\n{'f':>7}{'cos_agg':>10}{'rel/ship':>10}{'board cos':>11}{'pds':>10}"
          f"{'d_pds/sig':>10} | {'ovr pds':>9}{'ovr mse':>9}{'ovr nmae':>10}"
          f"{'TOTAL':>10}{'TOTAL cal':>11}  gate1")
    print("-" * 122)
    for r in rows:
        mark = " <- SHIPPED" if abs(r["f"] - 0.10) < 1e-9 else (
            " <- 207" if abs(r["f"] - 0.34) < 1e-9 else "")
        print(f"{r['f']:>7.3g}{r['cos_agg']:>10.6f}{r['rel_vs_ship']:>+10.4f}"
              f"{r['board_cos_k1']:>11.6f}{r['pds']:>10.6f}{r['pds_sigma']:>10.2f} | "
              f"{r['o_pds']:>+9.5f}{r['o_mse_k1']:>+9.5f}{r['o_nmae']:>+10.6f}"
              f"{r['total_k1']:>+10.5f}{r['total_kcal']:>+11.5f}"
              f"{'  Y' if r['within_gate1'] else '  N'}{mark}")

    best = max(rows, key=lambda r: r["total_k1"])
    best_g = max((r for r in rows if r["within_gate1"]), key=lambda r: r["total_k1"])
    best_cal = max(rows, key=lambda r: r["total_kcal"])
    top_cos = max(rows, key=lambda r: r["cos_agg"])
    summary = dict(
        k_cos=k_cos, shipped_cos_agg=ship,
        best_unconstrained=best, best_within_gate1=best_g, best_calibrated=best_cal,
        cosine_argmax_row=top_cos,
        board_cosine_at_cosine_argmax=top_cos["board_cos_k1"],
        board_cosine_at_cosine_argmax_cal=top_cos["board_cos_kcal"],
        leader_cosine=LEADER_COS,
        gap_closed_at_argmax_k1=(top_cos["board_cos_k1"] - OURS_COS)
        / (LEADER_COS - OURS_COS),
        gap_closed_at_gate1_best_k1=(best_g["board_cos_k1"] - OURS_COS)
        / (LEADER_COS - OURS_COS),
        gates=dict(shape=GATE_SHAPE, amplitude=GATE_AMP))
    print(f"\nbest priced total, unconstrained : f = {best['f']:g}  {best['total_k1']:+.5f} "
          f"(cal {best['total_kcal']:+.5f})")
    print(f"best priced total, inside Gate-1 : f = {best_g['f']:g}  {best_g['total_k1']:+.5f} "
          f"(cal {best_g['total_kcal']:+.5f})")
    print(f"cosine argmax f = {top_cos['f']:g}: implied board cosine "
          f"{top_cos['board_cos_k1']:.5f} (cal {top_cos['board_cos_kcal']:.5f}) "
          f"against the leaders' {LEADER_COS}")
    print(f"fraction of the direction gap to the leaders closed: "
          f"{summary['gap_closed_at_argmax_k1']:.1%} at the cosine argmax, "
          f"{summary['gap_closed_at_gate1_best_k1']:.1%} inside Gate-1")
    (WORK / "price216.json").write_text(json.dumps(dict(rows=rows, **summary), indent=2))


if __name__ == "__main__":
    main()
