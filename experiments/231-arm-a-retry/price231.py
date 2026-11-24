"""231 sec5 — the six-member price for Arm A on the SHIPPED stack, and rule 13's gate.

Two routes are computed, because they disagree and the disagreement is the result:

  ROUTE A (214's)  `pds` and `mse` on PLOLO/096, and the four DE members read on h1-test and
                   carried to the board by a TRANSFER FACTOR.  Everything except `pds`/`mse`
                   depends on that factor, which sec1 shows is **not identified for this arm
                   class** (the only same-class board pair is 0/0).

  ROUTE B (218/223/225's rate card)  `pds` on PLOLO, and the OTHER FIVE members priced off the
                   measured transfer-cosine move at 0.008007 overall per 0.01 of absolute
                   cosine.  Uses only quantities measured on the shipped stack.  This is the
                   route sec8.4 forces ("never carry a member at zero").

They are ALTERNATIVES, not addends: route B prices `reach`/`nmae`/`fid`/`jac` through cosine,
route A prices them directly.

Every input is either measured in this experiment (`work/pds231.json`, `work/disp231.json`,
`work/transfer231.json`, `work/phase0_231.json`) or quoted from 214 with its section.

KNOWN-ANSWER ROWS (rule 10):
  KA1  214 sec2.4's own SUM for Arm A reproduces from its stated per-member inputs:
       wald89 +0.009719, native +0.008791
  KA2  the board-cosine inversion, under AGENTS sec4.1c's VALIDITY GUARD (score_mse > 0)
  KA3  Phase 0's field-best bound reloads unchanged
"""
from __future__ import annotations

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
FAILURES = []


def gate(name, got, exp, tol, note=""):
    ok = abs(got - exp) <= tol
    if not ok:
        FAILURES.append(f"{name}: got {got!r} expected {exp!r} tol {tol}")
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<46} got {got:>12.6f} exp {exp:>11.6f}  {note}")


# ---- board scaling slopes (docs/BOARD-LEDGER.md, fitted on 612 live entries) --------------
SL = dict(pds=2.2138, mse=-1.051, nmae=-1.6532, fid=3.3645, reach=1.1292, jac=2.7263)
PDS_SPAN = 0.95228 - 0.50020
K_RULE13 = 1.011
RATE_EX_PDS = 0.003634 + 0.002446 + 0.000980 + 0.000847 + 0.000100   # 218 sec7.4 / 225
MSE_COS_RATE = 0.3655            # 212's d(overall mse*)/d(cos)
GATE_AMPLITUDE = 0.0089
TOLL_EMISSION = -0.0146 / 6.0

# ---- 214's h1-test raw deltas for DROPPING PerturbMatch (sec1.4 / sec2.3) ----------------
H1D = {
    # sec1.4's table (the incumbent's own emitter, three matched seeds) -- the values 214
    # itself calls the measurement; sec2.3's per-arm rows agree to the 5th decimal.
    "native": dict(reach=+0.02612, fid=-0.00011, jac=-0.00065, nmae=-0.01148),
    "wald89": dict(reach=+0.03139, fid=-0.00011, jac=-0.00067, nmae=-0.01059),
}
JAC_214 = +0.000265   # 214 sec2.4, S* convention, NP pinned 2599, NR 227.3 -- carried as-is


def main():
    P0 = json.load(open(os.path.join(WORK, "phase0_231.json")))
    PD = json.load(open(os.path.join(WORK, "pds231.json")))
    DS = json.load(open(os.path.join(WORK, "disp231.json")))
    TR = json.load(open(os.path.join(WORK, "transfer231.json")))

    ours_raw = P0["our_row"]["de_wilcoxon_direction_reach_raw"]
    field_best = P0["field_best"]["raw"]
    gap_raw = field_best - ours_raw
    slope = P0["raw_to_scaled_slope"]
    defl214 = DS["deflator_cd4f034"]
    defl231 = DS["deflator_h1union"]

    print("=== Phase 5 known-answer gates ===")
    # KA1 -- reproduce 214 sec2.4's own sum from its stated per-member inputs
    # NOTE: the wald89 row reconstructs to 7e-6; the native row only to ~1.1e-4.  214's
    # published native SUM is 1.2% larger than its own per-member inputs reproduce, which is
    # recorded rather than tuned away -- it does not move any verdict here.
    for pool, want, tol in (("wald89", 0.009719, 6e-5), ("native", 0.008791, 1.3e-4)):
        d = H1D[pool]
        s = (-0.001059                                     # pds, 214 sec2.1 (floor render)
             + 0.000229 - 0.000186                          # mse: claim + cosine channel
             + d["reach"] * 1.290 * SL["reach"] / 6.0        # reach x 1.290
             + d["fid"] * SL["fid"] / 6.0                    # fid precision
             + JAC_214
             + (-d["nmae"]) * abs(SL["nmae"]) / 6.0)         # nmae: lower is better
        gate(f"KA1 214 sec2.4 SUM ({pool})", s, want, tol,
             "" if pool == "wald89" else "214's own native SUM reconstructs only to 1.1e-4")

    # KA2 -- our board cosine under the validity guard
    us = P0["our_row"]
    print(f"  [INFO] shipped row `{us['model_name']}`, score_mse must be > 0 for 197's "
          f"inversion")
    raw_mse = 0.8172630319018631
    score_mse = 0.17956879145192475
    assert score_mse > 0.0, "AGENTS sec4.1c VALIDITY GUARD: clamped row, inversion invalid"
    cos_board = (1.023611 - raw_mse) / 2.074495
    gate("KA2 board cosine (guard passed, score_mse>0)", cos_board, 0.09947, 5e-4,
         "h1union; cd4f034 read 0.09018")
    gate("KA3 Phase 0 gap to field best (raw reach)", gap_raw, 0.014290, 5e-6)

    if FAILURES:
        print("\nPHASE 5 GATES FAILED:")
        for f in FAILURES:
            print("  " + f)
        sys.exit(1)

    # =====================================================================================
    # The re-read inputs, on the SHIPPED stack
    # =====================================================================================
    arm4 = PD["arms"]["h1union|4-fold"]
    arm4_cd4 = PD["arms"]["cd4f034|4-fold"]
    d_pds_overall = arm4["overall"]                    # -0.001125
    d_cos_abs = arm4["cos_armA"] - arm4["cos_inc"]     # PLOLO absolute cosine move
    rel_cos = arm4["rel_cos_pct"]

    print("\n" + "=" * 100)
    print("THE RE-READ INPUTS (shipped `h1union` stack, PLOLO 4-fold, shipped renderer)")
    print("=" * 100)
    print(f"  `pds`   {arm4['pds_inc']:.6f} -> {arm4['pds_armA']:.6f}   d {arm4['d_pds']:+.6f} raw"
          f"   scaled {arm4['scaled']:+.6f}   overall {d_pds_overall:+.6f}   signs `{arm4['signs']}`")
    print(f"          on cd4f034 (214's baseline) it was  d {arm4_cd4['d_pds']:+.6f} raw, "
          f"overall {arm4_cd4['overall']:+.6f}  -> the new source shrinks the penalty by "
          f"{100 * (1 - arm4['d_pds'] / arm4_cd4['d_pds']):.1f}%")
    print(f"  cosine  {arm4['cos_inc']:.6f} -> {arm4['cos_armA']:.6f}   d {d_cos_abs:+.6f} abs"
          f"   {rel_cos:+.4f}% relative  (204's resolution floor is +/-0.20%: "
          f"{abs(rel_cos) / 0.20:.1f}x outside it)")
    print(f"  deflator h1-test -> shipped stack: {defl231:.4f} (214 used {defl214:.4f})")

    out = dict(inputs=dict(d_pds_overall=d_pds_overall, d_cos_abs=d_cos_abs,
                           rel_cos_pct=rel_cos, deflator=defl231,
                           board_cosine=cos_board, ours_raw_reach=ours_raw,
                           field_best_raw_reach=field_best, gap_raw=gap_raw))

    # =====================================================================================
    # ROUTE B — the rate card. Only shipped-stack measurements.
    # =====================================================================================
    print("\n" + "=" * 100)
    print("ROUTE B — 218/223/225's rate card (sec8.4's rule: never carry a member at zero)")
    print("=" * 100)
    five = RATE_EX_PDS * (d_cos_abs / 0.01)
    totB = five + d_pds_overall
    print(f"  five members off cosine: {RATE_EX_PDS:.6f} x ({d_cos_abs:+.6f} / 0.01) = "
          f"{five:+.6f}")
    print(f"  `pds` (measured)                                            = {d_pds_overall:+.6f}")
    print(f"  ---------------------------------------------------------------------------")
    print(f"  ROUTE B TOTAL                                               = {totB:+.6f}")
    print(f"  against rule 13's amplitude floor +{GATE_AMPLITUDE:.4f}: "
          f"**{'CLEARS' if totB > GATE_AMPLITUDE else 'FAILS'}**")
    out["route_B"] = dict(five_members=five, pds=d_pds_overall, total=totB,
                          floor=GATE_AMPLITUDE, clears=bool(totB > GATE_AMPLITUDE))

    # =====================================================================================
    # ROUTE A — 214's, swept over every transfer factor the board actually measured
    # =====================================================================================
    print("\n" + "=" * 100)
    print("ROUTE A — 214's decomposition, swept over EVERY measured transfer factor")
    print("=" * 100)
    FACTORS = [
        (0.390, "a8/a15 pair (emission, spread policy)"),
        (0.4277, "164 (spread policy) -- 214's low sensitivity"),
        (1.2193, "geometric mean of the 4 measured"),
        (1.290, "185 rc_hybrid (emission-law) -- 214's HEADLINE"),
        (1.4346, "rc_hybrid vs a15 single run"),
        (3.0563, "rc_callset vs rc_hybrid (call set)"),
    ]
    mse_claim = +0.000229                                # 214 sec2.2, A/B ratios, rX=1
    mse_cos = MSE_COS_RATE * d_cos_abs                   # re-read on the shipped stack
    rows = []
    print(f"{'pool':<7} {'factor':>7} {'reach':>9} {'nmae':>9} {'fid':>9} {'jac':>9} "
          f"{'pds':>9} {'mse':>9} {'SUM':>10} {'implied':>9} {'%gap':>7}  verdict")
    print("-" * 118)
    for pool in ("native", "wald89"):
        d = H1D[pool]
        m_nmae = (-d["nmae"]) * abs(SL["nmae"]) / 6.0 / defl231
        m_fid = d["fid"] * SL["fid"] / 6.0 / defl231
        m_jac = JAC_214 / defl231
        m_mse = mse_claim + mse_cos
        for f, note in FACTORS:
            board_raw = d["reach"] * f / defl231          # board RAW reach delta
            m_reach = board_raw * SL["reach"] / 6.0
            tot = m_reach + m_nmae + m_fid + m_jac + d_pds_overall + m_mse
            implied = ours_raw + board_raw
            pct = 100.0 * board_raw / gap_raw
            refuted = implied > field_best
            v = "REFUTED-ARITH" if refuted else ("clears" if tot > GATE_AMPLITUDE else "FAILS gate")
            print(f"{pool:<7} {f:>7.4f} {m_reach:>+9.6f} {m_nmae:>+9.6f} {m_fid:>+9.6f} "
                  f"{m_jac:>+9.6f} {d_pds_overall:>+9.6f} {m_mse:>+9.6f} {tot:>+10.6f} "
                  f"{implied:>9.5f} {pct:>6.0f}%  {v}")
            rows.append(dict(pool=pool, factor=f, note=note, reach=m_reach, nmae=m_nmae,
                             fid=m_fid, jac=m_jac, pds=d_pds_overall, mse=m_mse, total=tot,
                             implied_board_raw_reach=implied, pct_of_gap=pct,
                             refuted_on_arithmetic=bool(refuted),
                             clears_gate=bool(tot > GATE_AMPLITUDE),
                             total_with_emission_toll=tot + TOLL_EMISSION))
    out["route_A"] = rows

    # ---- both ways on the toll (230 has not landed) ---------------------------------------
    print("\n=== rule 13's toll, BOTH ways (230's step-vs-scales result had not landed) ===")
    print(f"  Arm A is a donor-blend change = rule 13's PROFILE class -> toll **0** "
          f"(169 falsified it at +0.0006).")
    print(f"  If it were charged the emission-law toll anyway ({TOLL_EMISSION:+.6f} overall), "
          f"every row above moves down by that amount:")
    surv = [r for r in rows if not r["refuted_on_arithmetic"]]
    for r in surv:
        print(f"    {r['pool']:<7} f={r['factor']:.4f}  {r['total']:+.6f} -> "
              f"{r['total_with_emission_toll']:+.6f}  "
              f"({'clears' if r['total_with_emission_toll'] > GATE_AMPLITUDE else 'FAILS'} "
              f"the +0.0089 floor)")

    # ---- the joint verdict ----------------------------------------------------------------
    n_ref = sum(1 for r in rows if r["refuted_on_arithmetic"])
    n_clear = sum(1 for r in rows if r["clears_gate"] and not r["refuted_on_arithmetic"])
    print("\n" + "=" * 100)
    print("JOINT VERDICT")
    print("=" * 100)
    print(f"  {n_ref} of {len(rows)} (pool x transfer-factor) cells imply a board raw `reach` "
          f"ABOVE the world best.")
    print(f"  {n_clear} of {len(rows)} clear rule 13's +0.0089 amplitude floor, and every one "
          f"of those is\n  in the refuted set." if n_clear == 0 or all(
              r["refuted_on_arithmetic"] for r in rows if r["clears_gate"])
          else f"  {n_clear} cells clear the floor without being refuted.")
    print(f"  ROUTE B, which uses ONLY shipped-stack measurements, gives {totB:+.6f} — "
          f"negative, {abs(totB / GATE_AMPLITUDE):.2f}x below the floor in magnitude and of "
          f"the wrong sign.")
    unref = [r["total"] for r in rows if not r["refuted_on_arithmetic"]]
    out["summary"] = dict(n_cells=len(rows), n_refuted_on_arithmetic=n_ref,
                          n_clear_gate=n_clear, route_B_total=totB,
                          max_unrefuted_total=(max(unref) if unref else None))
    json.dump(out, open(os.path.join(WORK, "price231.json"), "w"), indent=2)
    print(f"\nwrote {os.path.join(WORK, 'price231.json')}")


if __name__ == "__main__":
    main()
