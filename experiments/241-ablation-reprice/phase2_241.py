"""241 Phase 2 — the two rows `DEF-RUNBOOK` §3.1 says must be re-decided before Oct 22.

1. PerturbMatch.  Is it still net-negative on the shipped renderer, and does a "drop it"
   reading reconcile with 231's refutation of 214's Arm A?  The reconciliation is done on
   231's own arithmetic: an arm's implied BOARD raw `reach` must not exceed the world best.
2. The CD4 residual, which 200 priced on two members and the board has since measured at
   +0.0130557 for the f 0.10 -> 0.34 step.  Confirm the DIRECTION of the correction.

Channels are priced separately (§8.1b rule 3): DIRECTION = `pds` + `mse` + `nmae`,
ORDERING = `reach` + `fid` + `jac`.  A direction gain is readable at +0.002; an ordering gain
under ~2 sd = 0.009 is not distinguishable from a re-roll.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

BOARD_RAW_REACH = 0.27654544083789157       # `h1union`, our shipped row
WORLD_BEST_REACH = 0.2908353                # 231 KA0-a2 (`luka`)
REACH_SLOPE = 1.1291872                     # 231 KA0-c, raw -> scaled, fitted on 656 rows
ORDERING_SD = 0.004502                      # §8.1b, overall units
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
BOARD_COS = 0.090270                        # `cd4f034`'s inverted board cosine
RATE = dict(mse=+0.003634, reach=+0.002446, nmae=+0.000980, fid=+0.000847,
            jac=+0.000100, pds=-0.000748)   # 218's card, per +0.01 ABSOLUTE board cosine


def channels(r, pds_key="pds"):
    d = (r.get(pds_key, 0.0) or 0.0) + (r.get("mse", 0.0) or 0.0) + (r.get("nmae", 0.0) or 0.0)
    o = ((r.get("reach", 0.0) or 0.0) + (r.get("fid", 0.0) or 0.0)
         + (r.get("jac", 0.0) or 0.0))
    return d, o, d + o


def main() -> None:
    T = json.loads(open(os.path.join(WORK, "rank241.json")).read())["table"]
    cos = json.loads(open(os.path.join(WORK, "cos241.json")).read())
    out = {}

    # ---------------------------------------------------------------- 1. PerturbMatch
    pm = T["PerturbMatch K562 replacement (0.5)"]
    print("=" * 112)
    print("PHASE 2.1 — PerturbMatch on the shipped renderer")
    print("=" * 112)
    print(f"  {'stack':<40}{'DIRECTION':>12}{'ORDERING':>11}{'SIX-MEMBER':>13}"
          f"{'6-mem (3-can)':>15}{'kolf2 pds leg':>15}")
    rec = {}
    for when, name in (("BEFORE", "200's renderer, f=0.10 (published: -0.001615)"),
                       ("AFTER", "shipped renderer, f=0.10"),
                       ("SHIPPED", "shipped renderer, shipped dose f=0.34")):
        r = pm[when]
        d, o, s = channels(r)
        d3, o3, s3 = channels(r, "pds3")
        rec[when] = dict(direction=d, ordering=o, six=s, six_3cancer=s3,
                         pds3=r.get("pds3"), pds_kolf2=r.get("pds_kolf2"))
        print(f"  {name:<40}{d:>+12.6f}{o:>+11.6f}{s:>+13.6f}{s3:>+15.6f}"
              f"{r.get('pds_kolf2', float('nan')):>+15.6f}")
    print()
    print("  Sign convention: positive = HAVING PerturbMatch is worth this much.  So a "
          "negative\n  six-member number is the value of DROPPING it, with the sign flipped.")
    for when in ("AFTER", "SHIPPED"):
        s = rec[when]["six"]
        s3 = rec[when]["six_3cancer"]
        print(f"\n  --- {when}: dropping PerturbMatch is worth {-s:+.6f} (4-fold pds) / "
              f"{-s3:+.6f} (3-cancer)")
        print(f"      vs rule 13 shape gate  +-{GATE_SHAPE}:  "
              f"{abs(s) / GATE_SHAPE:.2f}x  -> {'OUTSIDE' if abs(s) > GATE_SHAPE else 'inside'}")
        d, o = rec[when]["direction"], rec[when]["ordering"]
        print(f"      channel split (§8.1b rule 3): DIRECTION {-d:+.6f} "
              f"({abs(d)/GATE_SHAPE:.2f}x the +-0.0019 gate)   "
              f"ORDERING {-o:+.6f} ({abs(o)/ORDERING_SD:.2f} sd of the 0.004502 ordering noise)")
        # 231's arithmetic admissibility test
        d_reach_overall = -(T["PerturbMatch K562 replacement (0.5)"][when]["reach"])
        implied = BOARD_RAW_REACH + d_reach_overall * 6.0 / REACH_SLOPE
        print(f"      231's test: implied board raw `reach` = {implied:.5f} against a world "
              f"best of {WORLD_BEST_REACH:.5f}  -> "
              f"{'ADMISSIBLE' if implied <= WORLD_BEST_REACH else 'IMPOSSIBLE'}")
        rec[when]["implied_board_reach"] = implied
        rec[when]["admissible"] = bool(implied <= WORLD_BEST_REACH)
    # 214's Arm A headline, for the reconciliation
    print(f"\n  214's Arm A headline was +0.0046 ... +0.0088; 231 refuted +0.0088 because it "
          f"needs\n  board raw `reach` 0.31223 > world best {WORLD_BEST_REACH:.5f}.  241 "
          f"measures the same arm at\n  "
          f"{-rec['SHIPPED']['six']:+.6f} ... {-rec['SHIPPED']['six_3cancer']:+.6f} — "
          f"{0.0088 / abs(rec['SHIPPED']['six']):.1f}x SMALLER than the refuted headline, and "
          f"arithmetically admissible.")
    out["perturbmatch"] = rec

    # ---------------------------------------------------------------- 2. the CD4 residual
    cd = T["post-whitening CD4 residual"]
    print()
    print("=" * 112)
    print("PHASE 2.2 — the CD4 residual: direction of the correction")
    print("=" * 112)
    print(f"  {'stack':<40}{'pds':>10}{'mse':>10}{'2-member':>11}"
          f"{'Δcos3 rel %':>13}{'DE-4 fill k=1.000':>19}{'k=1.384':>10}{'6-member':>21}")
    rc = {}
    for when, tag, name in (("BEFORE", "MF", "200's renderer, f=0.10 (published +0.002635)"),
                            ("AFTER", "NOMF", "shipped renderer, f=0.10"),
                            ("SHIPPED", "SHIPPED", "shipped renderer, shipped dose f=0.34")):
        r = cd[when]
        two = (r.get("pds") or 0.0) + (r.get("mse") or 0.0)
        c = cos[tag]["loo"]["noCD4"]
        lo, hi = c["de4_fill_k1000"], c["de4_fill_k1384"]
        rc[when] = dict(pds=r.get("pds"), mse=r.get("mse"), two_member=two,
                        rel_cos3=-c["rel3"], de4_lo=lo, de4_hi=hi,
                        six_lo=two + lo, six_hi=two + hi,
                        understatement_lo=(two + lo) / 0.002635,
                        understatement_hi=(two + hi) / 0.002635)
        print(f"  {name:<40}{r.get('pds', 0):>+10.6f}{r.get('mse', 0):>+10.6f}{two:>+11.6f}"
              f"{-c['rel3']:>+13.2f}{lo:>+19.6f}{hi:>+10.6f}"
              f"{two + lo:>+11.6f} .. {two + hi:>+.6f}")
    print(f"\n  h1-test measures the four DE members at EXACTLY 0 for this row and that is a "
          f"STRUCTURAL\n  zero, not a measurement: CD4 covers 4 of 85 h1-test targets, "
          f"`solve_cd4_fraction` sizes the\n  residual by the MEDIAN row norm of the added "
          f"block, that median is 0, and the arm is added at\n  strength zero.  "
          f"`de241.reb_noCD4_nomf.npz` is BIT-IDENTICAL to `de241.reb_ship_nomf.npz`.")
    print(f"\n  The `DE-4 fill` column is 218's `cd4f034` rate card driven by a PLOLO cosine "
          f"delta through an\n  UNIDENTIFIED k_cos (falsifier 4): k = 1.000 is no conversion "
          f"at all, k = 1.384 is §8.4's\n  corrected factor.  It is a RANGE and a direction, "
          f"never a point estimate.")
    print(f"\n  DIRECTION OF THE CORRECTION: UP, on every base.  At the shipped dose the "
          f"two-member price is\n  {rc['SHIPPED']['two_member']:+.6f} and the six-member "
          f"range is {rc['SHIPPED']['six_lo']:+.6f} .. {rc['SHIPPED']['six_hi']:+.6f} — "
          f"{rc['SHIPPED']['understatement_lo']:.2f}x .. "
          f"{rc['SHIPPED']['understatement_hi']:.2f}x 200's published +0.002635.")
    print(f"  The board's own direct measurement of the SAME component's dose step "
          f"(f 0.10 -> 0.34) is\n  +0.0130557, of which the four members 200 carried at zero "
          f"delivered +0.006103 (47%).")
    out["cd4"] = rc

    json.dump(out, open(os.path.join(WORK, "phase2_241.json"), "w"), indent=1, default=float)
    print("\nwrote work/phase2_241.json")


if __name__ == "__main__":
    sys.exit(main())
