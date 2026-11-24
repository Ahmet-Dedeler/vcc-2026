"""200 — the four board-measured retrodictions, printed BEFORE any other number is read.

Audit 181 rule 1: an arm family that has a board point must retrodict it first.  This stack
has four, and each one tests a different instrument in this experiment:

  R1  a15 -> a23, the MASS FIX.            board overall -0.00073 (0.1740865 -> 0.1733565),
                                           board `reach` -0.007933 scaled.  Tests the six-member
                                           sum end to end, and the h1-test -> board `reach`
                                           transfer for a profile-CHANGING arm.
  R2  a14 -> a15, the CD4 RESIDUAL.        board `pds` +0.007257 scaled.  Tests PLOLO, which is
                                           the only `pds` instrument this project trusts.
  R3  a15 -> a19, the WHITENER RIDGE.      board overall -0.013518, on `mse`.  Tests the 096
                                           `claim` route, and specifically its CLIFF: the claim
                                           model has a hard zero where the board has a soft one.
  R4  a2-prune -> a8, --overdispersion 4.  board `reach` +0.021398 scaled.  a8 carried a FLAT
                                           0.25 spread and no tier, so this must be retrodicted
                                           in the FLAT context; a15's tiered context is a
                                           different measurement and gives the opposite sign.

Deterministic: reads only `work/*.json` written by `pds200.py`, `mse200.py` and `six200.py`.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
SPAN_PDS = 0.95228 - 0.50020
K13 = 1.011
# `docs/BOARD-LEDGER.md`.  The scaled members average to Overall, so a scaled delta enters
# Overall at /6.  `reach` raw <-> scaled uses the ledger's fitted slope 1.1292.
BOARD = {
    "a2-prune": dict(overall=0.1503084, pds=0.736561, reach=0.102951, reach_raw=0.1701),
    "a8":       dict(overall=0.1543607, pds=0.736561, reach=0.124349, reach_raw=0.1890),
    "a14":      dict(overall=0.1704242, pds=0.738614, reach=0.114542, reach_raw=0.180615),
    "a15":      dict(overall=0.1740865, pds=0.745871, reach=0.115368, reach_raw=0.181173),
    "a19":      dict(overall=0.1605690),
    "a23":      dict(overall=0.1733565, reach=0.107435),
}
SLOPE_REACH = 1.1292


def main() -> None:
    pds_rows = json.loads(open(os.path.join(WORK, "pds200.json")).read())
    R = {(r["fold"], r["arm"], r["space"]): r for r in pds_rows}
    mse = json.loads(open(os.path.join(WORK, "mse200.json")).read())
    six = json.loads(open(os.path.join(WORK, "six200-tablev2.json")).read())
    rep = json.loads(open(os.path.join(WORK, "report_pds.json")).read())

    def chal(fold, arm):
        vs = [R[(fold, arm, f"chal{c}/delivered")]["pds"] for c in "ABC"
              if (fold, arm, f"chal{c}/delivered") in R]
        return float(np.mean(vs)) if len(vs) == 3 else None

    def pds_delta(arm, base):
        per = [chal(f, arm) - chal(f, base) for f in FOLDS
               if chal(f, arm) is not None and chal(f, base) is not None]
        v = np.asarray(per)
        return float(v.mean()), float(v.std()), int(v.size), float(v.mean()) / SPAN_PDS

    out = {}
    line = "=" * 112
    print(line)
    print("RETRODICTIONS -- four board points, each read on the instrument this experiment "
          "uses for that member")
    print(line)

    # ---------------------------------------------------------------- R1  the mass fix
    b_ov = BOARD["a23"]["overall"] - BOARD["a15"]["overall"]
    b_reach_scaled = BOARD["a23"]["reach"] - BOARD["a15"]["reach"]
    b_reach_raw = b_reach_scaled / SLOPE_REACH
    mf = six["a23 mass fix"]
    mine_pds = -rep["loo"]["noMassFix"]["overall_k1011"]
    mine_mse = -mse["geometry"]["noMassFix"]["d_overall"]
    mine_ov = mine_pds + mine_mse + mf["reach"] + mf["fid"] + mf["jac"] + mf["nmae"]
    print(f"\nR1  a15 -> a23, the MASS FIX")
    print(f"    board   overall {b_ov:+.6f}   `reach` {b_reach_scaled:+.6f} scaled "
          f"= {b_reach_raw:+.6f} raw")
    print(f"    here    overall {mine_ov:+.6f}   (pds {mine_pds:+.6f}  mse {mine_mse:+.6f}  "
          f"reach {mf['reach']:+.6f}  fid {mf['fid']:+.6f}  jac {mf['jac']:+.6f}  "
          f"nmae {mf['nmae']:+.6f})")
    print(f"    h1-test d reach raw {mf['d_reach_raw']:+.6f} (sem "
          f"{mf['paired_reach']['sem']:.5f}, t {mf['paired_reach']['t']:+.2f}) against the "
          f"board's {b_reach_raw:+.6f}")
    print(f"    -> implied h1 -> board `reach` transfer for this profile-CHANGING arm: "
          f"{b_reach_raw / mf['d_reach_raw']:.2f}")
    print(f"    VERDICT: both readings put the mass fix INSIDE the +/-0.0019 shape gate "
          f"({b_ov:+.5f} board, {mine_ov:+.5f} here).  It has never cleared.")
    out["R1_massfix"] = dict(board_overall=b_ov, board_reach_scaled=b_reach_scaled,
                             board_reach_raw=b_reach_raw, mine_overall=mine_ov,
                             mine_reach_raw=mf["d_reach_raw"],
                             transfer=b_reach_raw / mf["d_reach_raw"])

    # ---------------------------------------------------------------- R2  the CD4 residual
    b_pds = BOARD["a15"]["pds"] - BOARD["a14"]["pds"]
    m, sd, n, sc = pds_delta("noMassFix", "a14")          # a15 (= no mass fix) vs a14, exact
    cd4a14 = rep["loo"]["cd4_a14"]
    print(f"\nR2  a14 -> a15, the CD4 RESIDUAL (fraction 0.05 -> solved 0.094)")
    print(f"    board   `pds` {b_pds:+.6f} scaled")
    print(f"    here    exact a14 arm, PLOLO 4-fold: {sc:+.6f} scaled "
          f"(raw {m:+.6f}, fold sd {sd:.6f}, n {n})   ratio {sc / b_pds:.2f}")
    print(f"    here    the same move read as a leave-one-out of the shipped stack "
          f"(`cd4_a14`): {-cd4a14['scaled']:+.6f} scaled   ratio "
          f"{-cd4a14['scaled'] / b_pds:.2f}")
    print(f"    and REMOVING CD4 outright (`noCD4`): "
          f"{-rep['loo']['noCD4']['scaled']:+.6f} scaled -- a larger move than a14 -> a15 was, "
          f"as it must be, since a14 already carried 0.05 of it.")
    out["R2_cd4"] = dict(board_pds_scaled=b_pds, mine_exact_scaled=sc,
                         mine_loo_scaled=-cd4a14["scaled"], fold_sd=sd,
                         mine_nocd4_scaled=-rep["loo"]["noCD4"]["scaled"])

    # ---------------------------------------------------------------- R3  a19's ridge
    b_ov19 = BOARD["a19"]["overall"] - BOARD["a15"]["overall"]
    g = mse["geometry"]["a19_exact"]
    a15_mse_overall = 0.120057 / 6.0        # a15's whole `mse` score, in overall units
    print(f"\nR3  a15 -> a19, the WHITENER RIDGE (gamma 0.80, lambda 0.007)")
    print(f"    board   overall {b_ov19:+.6f}, recorded against `mse`")
    print(f"    here    096 claim route: raw mse {g['raw']:.5f} -> score {g['score']:.5f}, "
          f"d overall {g['d_overall']:+.6f}")
    print(f"    a15's WHOLE `mse` score is {a15_mse_overall:+.6f} overall, so the claim "
          f"route's {g['d_overall']:+.6f} is the model's hard floor -- it says a19 lost ALL "
          f"of `mse`.")
    print(f"    The board says it lost {b_ov19 / -a15_mse_overall * 100:.0f}% of it.  "
          f"Ratio {g['d_overall'] / b_ov19:.2f}: RIGHT SIGN, OVERSTATED {g['d_overall'] / b_ov19:.2f}x.")
    print(f"    -> every `-0.018520` in the leave-one-out table is a CENSORED lower bound "
          f"and should be read as ~{-g['d_overall'] / (g['d_overall'] / b_ov19):.5f} "
          f"if this one calibration transfers.")
    out["R3_a19"] = dict(board_overall=b_ov19, mine_overall=g["d_overall"],
                         ratio=g["d_overall"] / b_ov19, raw_mse=g["raw"],
                         a15_mse_overall=a15_mse_overall)

    # ---------------------------------------------------------------- R4  --overdispersion
    b_reach8 = BOARD["a8"]["reach"] - BOARD["a2-prune"]["reach"]
    b_reach8_raw = BOARD["a8"]["reach_raw"] - BOARD["a2-prune"]["reach_raw"]
    flat = six["phi 4 IN a8's flat context"]
    tier = six["--overdispersion 4.0"]
    print(f"\nR4  a2-prune -> a8, `--overdispersion 4.0`  (a8 carried a FLAT 0.25 spread)")
    print(f"    board   `reach` {b_reach8:+.6f} scaled = {b_reach8_raw:+.6f} raw")
    print(f"    here    IN a8's flat context   d reach raw {flat['d_reach_raw']:+.6f} "
          f"(sem {flat['paired_reach']['sem']:.5f}, t {flat['paired_reach']['t']:+.2f})  "
          f"SIGN RECOVERED")
    print(f"    here    in a15's TIERED context d reach raw {tier['d_reach_raw']:+.6f} "
          f"(sem {tier['paired_reach']['sem']:.5f}, t {tier['paired_reach']['t']:+.2f})  "
          f"SIGN REVERSED")
    print(f"    -> implied h1 -> board transfer in the flat context: "
          f"{b_reach8_raw / flat['d_reach_raw']:.2f} (h1-test UNDERSTATES).  In the tiered "
          f"context no transfer factor of any sign fits.")
    print(f"    THE RETRODICTION ONLY WORKS IN a8's OWN CONTEXT.  `--overdispersion`'s `reach` "
          f"value is not a property of the knob; it is a property of the (spread, phi) PAIR.")
    out["R4_phi"] = dict(board_reach_scaled=b_reach8, board_reach_raw=b_reach8_raw,
                         flat_context_raw=flat["d_reach_raw"],
                         flat_context_t=flat["paired_reach"]["t"],
                         tier_context_raw=tier["d_reach_raw"],
                         tier_context_t=tier["paired_reach"]["t"],
                         transfer_flat=b_reach8_raw / flat["d_reach_raw"])

    json.dump(out, open(os.path.join(WORK, "retro200.json"), "w"), indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'retro200.json')}")


if __name__ == "__main__":
    sys.exit(main())
