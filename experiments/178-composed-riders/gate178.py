"""178 — apply AGENTS.md 11.1 rule 13's composed slot gate to every arm, and show the arithmetic.

    predicted = PLOLO d_pds/6 * 1.011
              + 096 d_mse/6            (if amplitude / spread / phi moves)
              - 0.0146/6               (if the EMITTED PROFILE or the emission law changes)
    SHIP only if predicted > +0.0019 (shape-only)  or  > +0.0089 (amplitude-changing)

Every arm here moves `WHITEN_GAMMA`, the donor mixture or the sign pattern, all of which
reshape the profile the renderer reads, so:
  * the `mse` term is IN (162 gate 2 calls gamma amplitude-class, and 097 measured gamma's
    `mse` move directly);
  * the emitted-profile `reach` toll is IN;
  * the binding threshold is the amplitude-class one, **+0.0089**.

PLOLO `d_pds` is the 4-fold mean of the challenge-frame (A/B/C) delta against `a23`, divided by
the board `pds` span 0.452 (162/retrodict.py). The `k562_gwps` fold cannot express the donor
weighting — there is no K562 in it to up-weight — so its contribution is exactly 0.0 by
construction, and the 3-fold mean is printed beside the 4-fold one rather than substituted for
it.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
SPAN_PDS = 0.452
K_CAL = 1.011
TOLL = -0.0146
GATE_SHAPE, GATE_AMP = 0.0019, 0.0089
FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
BASE = "a23"
ARMS = ("g045", "g030", "flip", "w211", "flip_w211", "g045_flip", "g045_w211",
        "all3_g045", "all3_g030")
PARTS = {"flip_w211": ("flip", "w211"),
         "g045_flip": ("g045", "flip"),
         "g045_w211": ("g045", "w211"),
         "all3_g045": ("g045", "flip", "w211"),
         "all3_g030": ("g030", "flip", "w211")}


def main() -> None:
    rows = json.loads((WORK / "plolo178.json").read_text())
    d = {(r["fold"], r["arm"], r["space"]): r for r in rows}

    def chal(fold, arm, field="pds"):
        return float(np.mean([d[(fold, arm, f"chal{c}/delivered")][field] for c in "ABC"]))

    # ---- known-answer rows against 162's own arms.json ----------------------------------
    ref = json.loads(
        (HERE.parents[0] / "162-retrodicting-validator/work/arms.json").read_text())
    rd = {(r["fold"], r["arm"], r["space"]): r for r in ref}
    print("=== KA: 162's own challenge-frame endpoints, reproduced ===")
    worst = 0.0
    for fold in FOLDS:
        for arm in ("a15", "a23"):
            for c in "ABC":
                a = d[(fold, arm, f"chal{c}/delivered")]["pds"]
                b = rd[(fold, arm, f"chal{c}/delivered")]["pds"]
                worst = max(worst, abs(a - b))
    print(f"  worst |d pds| over 4 folds x 2 arms x 3 contexts: {worst:.2e}  "
          f"{'PASS' if worst < 1e-9 else 'FAIL'}")

    # ---- per-fold and mean scaled d_pds --------------------------------------------------
    print(f"\n=== PLOLO d_pds vs {BASE}, challenge frame, scaled by span {SPAN_PDS} ===")
    print(f"{'arm':<12}" + "".join(f"{f[:12]:>14}" for f in FOLDS)
          + f"{'mean4':>10}{'mean3*':>10}{'sd':>9}")
    dpds = {}
    for arm in ARMS:
        per = {f: (chal(f, arm) - chal(f, BASE)) / SPAN_PDS for f in FOLDS}
        m4 = float(np.mean(list(per.values())))
        m3 = float(np.mean([per[f] for f in FOLDS if f != "k562_gwps"]))
        dpds[arm] = dict(per_fold=per, mean4=m4, mean3=m3,
                         sd=float(np.std(list(per.values()), ddof=1)))
        print(f"{arm:<12}" + "".join(f"{per[f]:>+14.6f}" for f in FOLDS)
              + f"{m4:>+10.6f}{m3:>+10.6f}{dpds[arm]['sd']:>9.6f}")
    print("  *mean3 drops the k562-held fold, the one in which the donor weighting is a "
          "no-op by construction")

    # ---- the gate ------------------------------------------------------------------------
    # `chal178.py` supersedes `mse178.py`: it evaluates `min(claim_c, B_c)` per context at
    # the SHIPPED amplitude, which is 170's correction and the convention the board reads.
    ch = json.loads((WORK / "chal178_solve.json").read_text())
    b0 = ch["a15"]
    mse = {a: {"d_overall": (r["scaled"] - b0["scaled"]) / 6.0,
               "d_overall_knee": (r["knee_scaled"] - b0["knee_scaled"]) / 6.0,
               "alpha_star": r["alpha_star"], "rho": r["rho"], "raw": r["raw"]}
           for a, r in ch.items()}
    print(f"\n=== rule 13 composed gate (amplitude class, threshold +{GATE_AMP}) ===")
    print(f"{'arm':<12}{'pds term':>11}{'mse term':>11}{'reach toll':>12}"
          f"{'PREDICTED':>12}{'gate':>10}{'verdict':>10}")
    pred = {}
    for arm in ARMS:
        p = K_CAL * dpds[arm]["mean4"] / 6.0
        m = (mse[arm]["d_overall"] if mse and arm in mse else 0.0)
        t = TOLL / 6.0
        total = p + m + t
        pred[arm] = dict(pds_term=p, mse_term=m, toll=t, predicted=total)
        print(f"{arm:<12}{p:>+11.6f}{m:>+11.6f}{t:>+12.6f}{total:>+12.6f}"
              f"{GATE_AMP:>10.4f}{'SHIP' if total > GATE_AMP else 'NO-GO':>10}")

    # ---- sensitivity: the mse term after re-tuning alpha onto the arm's own knee ------
    # Legitimate only if the arm actually ships at that alpha AND `nmae`/`fid`/`reach` are
    # re-read there too, which they are not -- a15's alpha is board-tuned and every arm here
    # would need `alpha` between 0.90 and 1.02. Reported as an upper bound, not the verdict.
    print("\n=== sensitivity: `mse` after re-tuning VCC_ALPHA onto each arm's own knee ===")
    print(f"{'arm':<12}{'rho A/B/C':>22}{'alpha*':>9}{'mse (shipped)':>15}"
          f"{'mse (knee)':>12}{'PREDICTED (knee)':>18}{'verdict':>9}")
    for arm in ARMS:
        r = mse[arm]
        total = pred[arm]["pds_term"] + r["d_overall_knee"] + TOLL / 6.0
        print(f"{arm:<12}{'/'.join(f'{v:.3f}' for v in r['rho']):>22}"
              f"{r['alpha_star']:>9.4f}{r['d_overall']:>+15.6f}"
              f"{r['d_overall_knee']:>+12.6f}{total:>+18.6f}"
              f"{'SHIP' if total > GATE_AMP else 'NO-GO':>9}")

    # ---- composition -------------------------------------------------------------------
    print(f"\n=== composition: is the composite the sum of its riders? ===")
    print(f"{'composite':<12}{'channel':<10}{'composite':>12}{'sum parts':>12}"
          f"{'ratio':>8}  reading")
    for comp, members in PARTS.items():
        for label, get in (("pds", lambda a: dpds[a]["mean4"]),
                           ("mse", lambda a: (mse[a]["d_overall"] * 6 if mse else 0.0)),
                           ("overall", lambda a: pred[a]["predicted"] - TOLL / 6.0)):
            c = get(comp)
            s = sum(get(m) for m in members)
            ratio = c / s if abs(s) > 1e-12 else float("nan")
            read = ("super-additive" if ratio > 1.05 else
                    "sub-additive" if ratio < 0.95 else "additive")
            print(f"{comp:<12}{label:<10}{c:>+12.6f}{s:>+12.6f}{ratio:>8.3f}  {read}")
    (WORK / "gate178.json").write_text(json.dumps(
        {"dpds": dpds, "predicted": pred}, indent=1))
    print(f"\nwrote {WORK/'gate178.json'}")


if __name__ == "__main__":
    main()
