"""231 sec6 — Arm A priced in units of 230's board `reach` NOISE, not only against the floor.

230 (`experiments/230-toll-magnitude/RESULTS.md`) landed after the measurements here and changed
two things:

  1. **Rule 13's `-0.0146` `reach` toll is MIS-KEYED** -- its own calibration points a21/a22 are
     digest-identical to the incumbent on 900/900 blocks, so by the table's own test they belong
     in the zero-toll row, and the row the toll is charged to contains zero board points.
     **Charge 0.** (Arm A was already toll 0 by rule 13's third row; this removes the "both
     ways" caveat entirely.)

  2. **The unexplained spread of board d`reach` among DIGEST-IDENTICAL arms is
     sd 0.027011 scaled = 0.004502 on overall.**  That is the noise floor on the channel 219
     says carries every point of Arm A, and it is 2.4x rule 13's shape-only floor.

The asymmetry is the point: `pds` and `mse` are DETERMINISTIC on this board (`rc_hybrid`
reproduced both to sixteen significant figures), so the direction half of any arm is measurable
at a precision the ordering half is not.

    uv run python noise231.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

SD_REACH_OVERALL = 0.004502      # 230, sd of board d`reach` among digest-identical arms
SD_REACH_SCALED = 0.027011       # the same, on the scaled member
GATE_AMPLITUDE = 0.0089
GATE_SHAPE = 0.0019
SL_REACH = 1.1292
FAILURES = []


def gate(name, got, exp, tol, note=""):
    ok = abs(got - exp) <= tol
    if not ok:
        FAILURES.append(f"{name}: got {got!r} expected {exp!r} tol {tol}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<48} got {got:>11.6f} exp {exp:>10.6f}  {note}")


def main():
    P0 = json.load(open(os.path.join(WORK, "phase0_231.json")))
    PR = json.load(open(os.path.join(WORK, "price231.json")))
    PD = json.load(open(os.path.join(WORK, "pds231.json")))

    print("=== Phase 6 known-answer gates ===")
    gate("KA1 230's scaled sd -> overall (/6)", SD_REACH_SCALED / 6.0, SD_REACH_OVERALL, 5e-6)
    gate("KA2 230's '2.4x the entire floor'", SD_REACH_OVERALL / GATE_SHAPE, 2.4, 0.05)
    gap_raw = P0["gap_to_field_best_raw"]
    slope = P0["raw_to_scaled_slope"]
    gate("KA3 Phase 0 gap reloads", gap_raw, 0.014290, 5e-6)
    if FAILURES:
        print("\nPHASE 6 GATES FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    prize = gap_raw * slope / 6.0
    print("\n" + "=" * 96)
    print("THE CLOSURE, IN UNITS OF THE BOARD'S OWN `reach` NOISE")
    print("=" * 96)
    print(f"  230's board `reach` noise, overall                       sd = {SD_REACH_OVERALL:.6f}")
    print(f"  the ENTIRE remaining `reach` prize (ours -> world best)     = {prize:+.6f}"
          f"  = **{prize / SD_REACH_OVERALL:.2f} sd**")
    print(f"  rule 13's amplitude floor (Arm A's class)                   = {GATE_AMPLITUDE:+.6f}"
          f"  = {GATE_AMPLITUDE / SD_REACH_OVERALL:.2f} sd")
    print(f"  rule 13's shape-only floor                                  = {GATE_SHAPE:+.6f}"
          f"  = {GATE_SHAPE / SD_REACH_OVERALL:.2f} sd")
    print("\n  -> Becoming the BEST `reach` team in the world, from 7th of 656, is worth")
    print(f"     {prize / SD_REACH_OVERALL:.2f} sd of the board's own noise on that member. "
          f"**The whole channel is\n     under one standard deviation.**")

    print("\n" + "=" * 96)
    print("214's PUBLISHED RANGE AND THIS EXPERIMENT'S CELLS, IN sd")
    print("=" * 96)
    print(f"{'claim':<52} {'overall':>10} {'sd':>7}  {'vs floor':>9}  status")
    print("-" * 96)
    rows = [("214's headline, native x1.290 (no deflator)", 0.0088, "REFUTED on arithmetic"),
            ("214's x1.533-deflated headline", 0.0054, "REFUTED on arithmetic"),
            ("214's low end, reach at 164's 0.4277", 0.0046, "not refuted")]
    for lbl, v, st in rows:
        print(f"{lbl:<52} {v:>+10.6f} {v / SD_REACH_OVERALL:>7.2f}  "
              f"{v / GATE_AMPLITUDE:>8.2f}x  {st}")
    best = max(r["total"] for r in PR["route_A"] if not r["refuted_on_arithmetic"])
    bestrow = [r for r in PR["route_A"] if r["total"] == best][0]
    print(f"{'231 route A, best UNREFUTED cell (' + bestrow['pool'] + ' x0.4277)':<52} "
          f"{best:>+10.6f} {best / SD_REACH_OVERALL:>7.2f}  {best / GATE_AMPLITUDE:>8.2f}x  "
          f"fails floor AND noise")
    tb = PR["route_B"]["total"]
    print(f"{'231 route B (shipped-stack measurements only)':<52} {tb:>+10.6f} "
          f"{tb / SD_REACH_OVERALL:>7.2f}  {tb / GATE_AMPLITUDE:>8.2f}x  NEGATIVE")

    # the reach member alone, per cell
    print("\n  the `reach` MEMBER alone, per cell, in sd of its own noise:")
    for r in PR["route_A"]:
        if not r["refuted_on_arithmetic"]:
            print(f"    {r['pool']:<7} f={r['factor']:.4f}  reach {r['reach']:+.6f} = "
                  f"{r['reach'] / SD_REACH_OVERALL:.2f} sd")

    # ---- the measurable half vs the unmeasurable half ------------------------------------
    arm = PD["arms"]["h1union|4-fold"]
    d_pds = arm["overall"]
    d_cos = arm["cos_armA"] - arm["cos_inc"]
    five = 0.008007 * (d_cos / 0.01)
    print("\n" + "=" * 96)
    print("THE SPLIT 230 ASKS FOR — the MEASURABLE half against the UNMEASURABLE half")
    print("=" * 96)
    print(f"  MEASURABLE (deterministic on this board -- `pds` and `mse` reproduced to 16 s.f.):")
    print(f"    `pds`, measured on the shipped stack                    {d_pds:+.6f}")
    print(f"    transfer cosine {arm['rel_cos_pct']:+.4f}% rel = {d_cos:+.6f} abs, priced")
    print(f"      through the five-member rate card                     {five:+.6f}")
    print(f"    ------------------------------------------------------------------")
    print(f"    measurable half TOTAL                                   {d_pds + five:+.6f}")
    print(f"  UNMEASURABLE (the `reach`/`nmae` half, sd {SD_REACH_OVERALL:.6f}):")
    print(f"    best unrefuted estimate                                 "
          f"{best - d_pds - bestrow['mse']:+.6f}  = "
          f"{(best - d_pds - bestrow['mse']) / SD_REACH_OVERALL:.2f} sd")
    print(f"\n  **The half of Arm A that this board can measure is NEGATIVE. The half that would")
    print(f"    be positive is under {(best - d_pds - bestrow['mse']) / SD_REACH_OVERALL:.1f} sd "
          f"of the noise on its own channel.**")

    out = dict(sd_reach_overall=SD_REACH_OVERALL, remaining_prize_overall=prize,
               remaining_prize_sd=prize / SD_REACH_OVERALL,
               floor_amplitude_sd=GATE_AMPLITUDE / SD_REACH_OVERALL,
               floor_shape_sd=GATE_SHAPE / SD_REACH_OVERALL,
               claims_214={f"{v}": v / SD_REACH_OVERALL for _, v, _ in rows},
               best_unrefuted=best, best_unrefuted_sd=best / SD_REACH_OVERALL,
               route_B=tb, route_B_sd=tb / SD_REACH_OVERALL,
               measurable_half=d_pds + five, measurable_pds=d_pds, measurable_cosine=five,
               unmeasurable_half=best - d_pds - bestrow["mse"],
               unmeasurable_half_sd=(best - d_pds - bestrow["mse"]) / SD_REACH_OVERALL,
               toll_charged=0.0, toll_source="230: rule 13's -0.0146 is mis-keyed; charge 0")
    json.dump(out, open(os.path.join(WORK, "noise231.json"), "w"), indent=2)
    print(f"\nwrote {os.path.join(WORK, 'noise231.json')}")


if __name__ == "__main__":
    main()
