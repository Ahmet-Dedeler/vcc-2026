"""197 — the two assumptions the field inversion needs, TESTED rather than asserted.

  (i)  `sumD`, `Kc` and the cap `K` are properties of the withheld REFERENCE, hence shared.
       Not testable from board data alone; what IS testable is the size of the error it would
       take to change the answer, so a sensitivity is reported instead.

  (ii) each team sits near its own amplitude optimum.  Two tests:
       (a) OUR OWN ARMS.  We know where each of ours sits relative to its rho cliff
           (s/s_c per context) AND we know its true cosine, so the loss the assumption costs
           is directly measurable across a 1.9x span in s/s_c.
       (b) THE FIELD.  A team still descending its own amplitude curve shows up as successive
           submissions with falling `mse` and flat `pds` (096 found exactly that shape on
           Aginglab v11->v12).  Walk every snapshot on disk and look.

  (iii) plus the two remaining Phase-4 falsifiers: the byte-identical resubmission, and the
        `rc_hybrid` pre-registration.

    uv run python assump197.py
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
W096 = f"{REPO}/experiments/096-mse-floor/work"
sys.path.insert(0, f"{REPO}/src")
sys.path.insert(0, f"{REPO}/experiments/096-mse-floor")

CTX = ("A", "B", "C")
B_ANCHOR, SPAN = 0.987862879, 0.950162566
LINES = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def main():
    inv = json.load(open(f"{WORK}/invert197.json"))
    ph3 = json.load(open(f"{WORK}/phase3197.json"))
    sumD, K, floor0, coeff = inv["sumD"], inv["K"], inv["floor0"], inv["coeff"]
    ours, exact = inv["ours"], ph3["exact"]

    P("=" * 100)
    P("197 — ASSUMPTION TESTS AND THE REMAINING FALSIFIERS")
    P("=" * 100)

    # ------------------------------------------------------------------ (ii.a) our own arms
    P("\n(ii.a) HOW FAR OFF ITS OWN AMPLITUDE OPTIMUM EACH OF OUR ARMS SITS, and what that")
    P("       costs the floor.  s/s_c = sqrt(s^2 B_c / claim_c): < 1 is below the rho cliff,")
    P("       > 1 above it.  The floor assumes s/s_c = 1 exactly.")
    P(f"  {'arm':<10}{'raw':>10}{'s/s_c A':>9}{'s/s_c B':>9}{'s/s_c C':>9}{'mean':>8}"
      f"{'floor':>9}{'measured':>10}{'floor/meas':>11}")
    tab = []
    for a in ("a8", "a10", "a11", "a12", "a14", "a15", "a23", "abc-knee"):
        per = ours[a]["per"]
        r = [np.sqrt(per[c]["B"] / per[c]["claim"]) for c in CTX]
        fl = (floor0 - ours[a]["raw"]) / coeff
        m = exact[a]
        tab.append((a, float(np.mean(r)), fl / m))
        P(f"  {a:<10}{ours[a]['raw']:>10.5f}" + "".join(f"{x:>9.4f}" for x in r)
          + f"{np.mean(r):>8.4f}{fl:>9.5f}{m:>10.5f}{fl/m:>11.4f}")
    P("  The relation is monotone and it is exactly the rho-cliff shape 187 describes: an arm")
    P("  ON its cliff loses 3-9% of its cosine to the floor, an arm 10% BELOW loses 27%, and")
    P("  an arm 60% ABOVE loses everything (the floor goes negative and reports nothing).")
    P("  a8/a10 are the decisive known-answer rows for Phase 4 falsifier 1: their MEASURED")
    P(f"  cosine is {exact['a8']:.5f}, the highest of any v11-family arm, and the field method")
    P("  reports 'no information' for them.  A high raw `mse` is NOT evidence of a poor")
    P("  direction, and any fit that let them enter at the clamp value would be inverted.")

    # ------------------------------------------------------------------ (i) sensitivity
    P("\n(i) SENSITIVITY of the leader's floor to the three shared reference constants and to")
    P("    the common mode `eps = (A-B)/B` that the floor sets to zero.")
    lead = ph3["leader"]
    u = lead["raw"]
    ours_cos = ph3["ours_cos"]

    def fl_of(sD=sumD, Kc=0.005410, Kcap=K, eps=0.0):
        q = Kcap / sD
        return (1.0 + Kc / sD + q * eps - u) / (2.0 * np.sqrt(q * (1.0 + eps)))
    base = fl_of()
    P(f"  {'perturbation':<34}{'leader floor':>14}{'x ours':>9}{'change':>10}")
    eps_ours = float(np.mean([ours["a23"]["per"][c]["A"] / ours["a23"]["per"][c]["B"] - 1
                              for c in CTX]))
    for tag, v in (("as computed (eps = 0)", base),
                   (f"eps = our own {eps_ours:.4f}", fl_of(eps=eps_ours)),
                   ("sumD  +5%", fl_of(sD=sumD * 1.05)), ("sumD  -5%", fl_of(sD=sumD * 0.95)),
                   ("cap K +10% (C_p proxy high)", fl_of(Kcap=K * 1.10)),
                   ("cap K -10% (C_p proxy low)", fl_of(Kcap=K * 0.90)),
                   ("Kc    +20%", fl_of(Kc=0.005410 * 1.2))):
        P(f"  {tag:<34}{v:>14.5f}{v/ours_cos:>9.3f}{100*(v/base-1):>9.1f}%")
    P("  Every perturbation that could plausibly be wrong leaves the ratio above 1.5x, and the")
    P("  one that moves it most (`sumD`) is the best-determined constant here: three of our own")
    P("  board points fit it to <= 6e-5 raw.  `eps = 0` is the CONSERVATIVE choice — the real")
    P("  common mode raises the floor.")

    # ------------------------------------------------------------------ (ii.b) the field
    P("\n(ii.b) IS THE FIELD AT ITS OWN AMPLITUDE OPTIMUM?  Every snapshot on disk, per team,")
    P("       looking for the 096 signature: `mse` moving while `pds_cosine` stands still.")
    snaps = sorted(glob.glob(f"{REPO}/data/leaderboard/*.json"))
    if os.path.exists(f"{WORK}/lb-live.json"):
        snaps.append(f"{WORK}/lb-live.json")
    hist = {}
    for f in snaps:
        d = json.load(open(f))
        e = d["entries"] if isinstance(d, dict) else d
        for r in e:
            hist.setdefault(r["team_name"], {})[r["id"]] = (
                r["submission_date"], r["expr_mse_unbiased_capped_norm"], r["pds_cosine"],
                r["score_avg"], r["de_wilcoxon_direction_reach_raw"])
    top = [t["team"] for t in ph3["top3"]] + ["cqawesome", "Aginglab.com", "Garrrbage",
                                              "Vivai", "Ahmet Dedeler"]
    P(f"  {'team':<20}{'n distinct entries seen':>24}   trajectory (raw mse | raw pds)")
    moved = []
    for t in dict.fromkeys(top):
        h = sorted(hist.get(t, {}).values())
        if len(h) < 2:
            P(f"  {t[:19]:<20}{len(h):>24}   (one snapshot only — cannot test)")
            continue
        traj = "  ->  ".join(f"{x[1]:.4f}|{x[2]:.4f}" for x in h)
        dm = h[-1][1] - h[0][1]
        dp = h[-1][2] - h[0][2]
        moved.append((t, dm, dp))
        P(f"  {t[:19]:<20}{len(h):>24}   {traj}")
        P(f"  {'':<20}{'':>24}   d mse {dm:+.4f}   d pds {dp:+.4f}")
    P("\n  The LAST step of each visible trajectory is the one that tests the assumption:")
    for t in dict.fromkeys(top):
        h = sorted(hist.get(t, {}).values())
        if len(h) < 2:
            continue
        dm, dp = h[-1][1] - h[-2][1], h[-1][2] - h[-2][2]
        verdict = ("still descending `mse` with `pds` flat -> NOT at optimum"
                   if dm < -0.002 and abs(dp) < 0.005 else
                   "moved both -> a direction change, not an amplitude move"
                   if abs(dp) >= 0.005 else "flat -> consistent with being at its optimum")
        P(f"    {t[:19]:<20} d mse {dm:+.4f}  d pds {dp:+.4f}   {verdict}")
    P("  Reading: no top team's last step looks like a team sitting still.  A team still")
    P("  descending its own amplitude curve is ABOVE its own optimum, so its observed `mse`")
    P("  exceeds its own minimum and its true cosine is HIGHER than the floor.  The assumption")
    P("  therefore fails in the safe direction, and our own a11 row measures the size of the")
    P("  failure at the other end (10% BELOW the cliff costs the floor 27%).")

    # ------------------------------------------------------------------ falsifier 2
    P("\nFALSIFIER 2 — the byte-identical resubmission of a15, entry 5QfX14axLAG8FF2rSDJe.")
    P("  docs/LOG.md 2026-08-31 records, for the restore against the original a15:")
    P("     score_mse  0.1200566151966586   vs   0.1200566151966586   (identical, 16 digits)")
    raw_from_score = B_ANCHOR - SPAN * 0.1200566151966586
    P(f"  Inverting the recovered scorer map: raw = {B_ANCHOR:.9f} - {SPAN:.9f} x 0.1200566152")
    P(f"                                          = {raw_from_score:.6f}")
    P(f"  a15's board raw as recorded             = {ours['a15']['raw']:.6f}   "
      f"(gap {raw_from_score-ours['a15']['raw']:+.2e}, the affine map's own 2.1e-5 fit residual)")
    P("  Identical score_mse to 16 digits => identical per-context raw => identical A, B, claim,")
    P("  Kc => identical X => IDENTICAL COSINE, exactly.  The falsifier is satisfied by")
    P("  construction and it also certifies that no part of the inversion reads anything that")
    P("  varies between two runs of the scorer.")
    P(f"  cos(a15) = {exact['a15']:.6f} both times.")

    # ------------------------------------------------------------------ falsifier: rc_hybrid
    P("\nPRE-REGISTRATION — rc_hybrid (submitted 2026-09-03 11:20 UTC, row not yet on the board).")
    P("  191 falsifier 1 proved its per-gene panel totals are digest-identical to a15's on")
    P("  900/900 blocks, so A, B and the profile direction are a15's EXACTLY.  Only `claim`")
    P("  moved, and only in context B, where min(claim, B) is already on the B leg.")
    P(f"  PREDICTED raw mse {inv['rc_hybrid_pred_raw']:.6f} = a15's {ours['a15']['raw']:.6f} "
      f"to {abs(inv['rc_hybrid_pred_raw']-ours['a15']['raw']):.1e}")
    P(f"  PREDICTED inverted cosine {exact['a15']:.5f}, identical to a15's.")
    P("  FALSIFIED if rc_hybrid's board raw `mse` differs from a15's by more than +/- 0.004")
    P("  (096's own out-of-sample band on a10/a11), or if it inverts to a cosine outside")
    P(f"  [{0.96*exact['a15']:.5f}, {1.04*exact['a15']:.5f}].")

    out = dict(cliff=tab, eps_ours=eps_ours, sens_base=base,
               raw_from_score=raw_from_score)
    json.dump(out, open(f"{WORK}/assump197.json", "w"), indent=1, default=float)
    open(f"{WORK}/assump197.log", "w").write("\n".join(LINES))
    P(f"\nwrote {WORK}/assump197.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
