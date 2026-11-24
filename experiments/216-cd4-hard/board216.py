"""216 — INBOX-2: re-price the dose ladder on the REALISED BOARD RATES, not on PLOLO's.

`cd4f034` landed at **+0.0130557** (entry `QSe3JlNUYZ2BzvlcsL14`), 3.0x its own preregistered
+0.00437.  Both of the transfer factors this experiment's §4.3/§6.0 tables were built on are
now measurable directly against a board point, so the grid is re-priced against them.

ONE CORRECTION TO INBOX-2, and it matters
-----------------------------------------
INBOX-2 derives the PLOLO->board factor as **1.80** from "212 measured +13.83%, the board
delivered +24.9%".  Those are two DIFFERENT COSINE FUNCTIONALS.  212's +13.83% is the **macro
per-target** cosine on the a15 tag; 197's board floor inverts the **energy-weighted aggregate**
`X / (sqrt(SigmaD) * mean_c sqrt(A_c))`, which is `cos_agg`.  This experiment's Phase 0 hit the
same trap from the other side (four gates looked like failures purely because 204/212 publish
the a15 tag and 207 publishes a23).

On the MATCHED functional the factor is:

    PLOLO cos_agg   f=0.10 -> 0.34 :  0.020086 -> 0.023698   = +17.98%
    board floor     rc_callset -> cd4f034 : 0.072280 -> 0.090268 = +24.89%
    k_cos(board) = 24.89 / 17.98 = 1.384          (not 1.80)

The `pds` factor is measured the same way and INBOX-2's reading of it is confirmed:

    PLOLO d raw pds  f=0.10 -> 0.34 : -0.001837
    board d raw pds  rc_callset -> cd4f034 : -0.00365
    k_pds(board) = 1.987                          (207/162 assumed 1.011)

MSE IS PRICED WITHOUT A MODEL.  Both board points sit exactly on 197's line
`raw_mse = 1.023611 - 2.074495 * cos` (0.873667 vs 0.87367 and 0.836351 vs 0.83635), so the
`mse` leg is an affine function of the projected cosine and needs none of 096's constants.

    uv run python board216.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")

# --- the realised board point (experiments/210-cd4-build/BOARD-RESULT.md) ----------------
BOARD = {
    "rc_callset": dict(overall=0.1918638367, pds_s=0.745871, mse_s=0.120057,
                       nmae_s=0.094290, fid_s=-0.029268, reach_s=0.192898, jac_s=0.027335,
                       pds_raw=0.83712, mse_raw=0.87367, nmae_raw=0.94326),
    "cd4f034": dict(overall=0.2049195830, pds_s=0.737800, mse_s=0.159276,
                    nmae_s=0.104861, fid_s=-0.020127, reach_s=0.219300, jac_s=0.028408,
                    pds_raw=0.83347, mse_raw=0.83635, nmae_raw=0.93687),
}
FLOOR_A, FLOOR_B = 1.023611, 2.074495       # raw_mse = A - B * cos
SPAN_PDS, SPAN_MSE = 0.452, 0.950162566
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
SIGMA_PDS_PLOLO = 0.001494
LEADER_COS = 0.13211


def floor_cos(raw_mse):
    return (FLOOR_A - raw_mse) / FLOOR_B


def main() -> None:
    sweep = json.loads((WORK / "sweep216.json").read_text())
    by = {(r["fold"], r["arm"]): r for r in sweep}

    def m(arm, key):
        return float(np.mean([by[(f, arm)][key] for f in FOLDS]))

    inc_cos_p = m("f0.34", "cos_agg")
    ship_cos_p = m("f0.1", "cos_agg")
    inc_pds_p = m("f0.34", "pds")
    ship_pds_p = m("f0.1", "pds")

    c_inc = floor_cos(BOARD["cd4f034"]["mse_raw"])
    c_ship = floor_cos(BOARD["rc_callset"]["mse_raw"])
    rel_plolo = inc_cos_p / ship_cos_p - 1.0
    rel_board = c_inc / c_ship - 1.0
    k_cos = rel_board / rel_plolo
    d_pds_plolo = inc_pds_p - ship_pds_p
    d_pds_board = BOARD["cd4f034"]["pds_raw"] - BOARD["rc_callset"]["pds_raw"]
    k_pds = d_pds_board / d_pds_plolo

    # the two identity checks that make the mse leg model-free
    id_ship = FLOOR_A - FLOOR_B * c_ship
    id_inc = FLOOR_A - FLOOR_B * c_inc
    print(f"197 floor identity reproduces both board points: "
          f"{id_ship:.6f} vs {BOARD['rc_callset']['mse_raw']}, "
          f"{id_inc:.6f} vs {BOARD['cd4f034']['mse_raw']}")
    print(f"board cosine  {c_ship:.6f} -> {c_inc:.6f}  = {rel_board:+.4%}")
    print(f"PLOLO cos_agg {ship_cos_p:.6f} -> {inc_cos_p:.6f}  = {rel_plolo:+.4%}")
    print(f"k_cos (matched functional) = {k_cos:.4f}   [INBOX-2 quoted 1.80 on a MIXED pair]")
    print(f"k_pds = {d_pds_board:+.5f} / {d_pds_plolo:+.5f} = {k_pds:.4f}   "
          f"[207/162 used 1.011]\n")

    doses = sorted({r["f"] for r in sweep if r["kind"] == "cd4"})
    rows = []
    for v in doses:
        arm = "f%g" % v
        rel = (m(arm, "cos_agg") - inc_cos_p) / inc_cos_p          # vs the INCUMBENT
        cb = c_inc * (1.0 + k_cos * rel)
        d_mse_raw = -FLOOR_B * (cb - c_inc)
        o_mse = (-d_mse_raw / SPAN_MSE) / 6.0
        dp_plolo = m(arm, "pds") - inc_pds_p
        dp_board = k_pds * dp_plolo
        o_pds = dp_board / SPAN_PDS / 6.0
        rows.append(dict(
            f=v, cos_agg_plolo=m(arm, "cos_agg"), rel_vs_inc_plolo=rel,
            board_cos=cb, rel_vs_inc_board=cb / c_inc - 1.0,
            board_mse_raw=BOARD["cd4f034"]["mse_raw"] + d_mse_raw,
            d_pds_plolo=dp_plolo, d_pds_board=dp_board,
            pds_sigma_plolo=dp_plolo / SIGMA_PDS_PLOLO,
            o_pds=o_pds, o_mse=o_mse, total=o_pds + o_mse,
            board_cos_vs_leader=cb / LEADER_COS))

    print(f"{'f':>7}{'PLOLO cos':>11}{'rel/inc':>9}{'board cos':>11}{'rel/inc':>9}"
          f"{'raw mse':>10}{'d_pds raw':>11}{'ovr pds':>10}{'ovr mse':>10}"
          f"{'TOTAL vs inc':>14}{'x leader':>10}")
    print("-" * 122)
    for r in rows:
        mk = " <- INCUMBENT" if abs(r["f"] - 0.34) < 1e-9 else ""
        print(f"{r['f']:>7.3g}{r['cos_agg_plolo']:>11.6f}{r['rel_vs_inc_plolo']:>+9.3f}"
              f"{r['board_cos']:>11.6f}{r['rel_vs_inc_board']:>+9.3f}"
              f"{r['board_mse_raw']:>10.5f}{r['d_pds_board']:>+11.5f}{r['o_pds']:>+10.5f}"
              f"{r['o_mse']:>+10.5f}{r['total']:>+14.5f}{r['board_cos_vs_leader']:>10.3f}{mk}")

    best = max(rows, key=lambda r: r["total"])
    top = max(rows, key=lambda r: r["board_cos"])
    clears = [r for r in rows if r["total"] > GATE_SHAPE]
    print(f"\nargmax of the board-rate price: f = {best['f']:g}  {best['total']:+.5f} "
          f"against the shape gate {GATE_SHAPE:+.4f} -> "
          f"{'CLEARS' if best['total'] > GATE_SHAPE else 'VETO'}")
    print(f"doses clearing the shape gate: "
          f"{[r['f'] for r in clears] if clears else 'NONE'}")
    print(f"board-cosine ceiling: f = {top['f']:g} at {top['board_cos']:.6f} "
          f"= {top['board_cos']/c_ship:.4f}x the shipped f=0.10 point, "
          f"{top['board_cos']/LEADER_COS:.3f}x the leaders' 0.13211")
    print(f"the leaders' move was x2.056; the dose ceiling on BOARD rates is "
          f"x{top['board_cos']/c_ship:.4f}")

    # what the board says about the members this grid cannot price
    unp = (BOARD["cd4f034"]["fid_s"] - BOARD["rc_callset"]["fid_s"]
           + BOARD["cd4f034"]["reach_s"] - BOARD["rc_callset"]["reach_s"]
           + BOARD["cd4f034"]["jac_s"] - BOARD["rc_callset"]["jac_s"]) / 6.0
    nm = (BOARD["cd4f034"]["nmae_s"] - BOARD["rc_callset"]["nmae_s"]) / 6.0
    print(f"\nrealised at f=0.10 -> 0.34, for reference (NOT extrapolated):")
    print(f"  nmae            {nm:+.6f} overall  (PLOLO's 188 identity predicted -0.000152)")
    print(f"  fid+reach+jac   {unp:+.6f} overall  (unpriced; falsifier 5 fired favourably)")
    print(f"  total realised  {BOARD['cd4f034']['overall']-BOARD['rc_callset']['overall']:+.7f}")

    out = dict(k_cos=k_cos, k_pds=k_pds, c_ship=c_ship, c_inc=c_inc,
               rel_board=rel_board, rel_plolo=rel_plolo,
               identity_check=[id_ship, id_inc], rows=rows, best=best,
               cosine_ceiling=top, clears_shape_gate=[r["f"] for r in clears],
               realised_nmae_overall=nm, realised_unpriced_overall=unp,
               inbox2_factor_quoted=1.80,
               inbox2_correction="INBOX-2's 1.80 compares 212's MACRO per-target cosine "
                                 "with 197's ENERGY-WEIGHTED aggregate; on the matched "
                                 "functional the factor is 1.384.")
    (WORK / "board216.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
