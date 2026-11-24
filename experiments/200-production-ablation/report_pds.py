"""200 — read `work/pds200.json` into the leave-one-out / add-one-in table on PLOLO `pds`.

Predictor: 162's PRE-COMMITTED one -- the unweighted mean over all four donor folds, in the
CHALLENGE frame (mean of contexts A/B/C).  No fold is selected; the between-fold scatter is
reported as the instrument's own error bar, and any arm whose scatter exceeds its effect is
flagged rather than quoted.

Scaling: raw -> scaled by the board pds span 0.45208 (`debench.ANCHOR`), scaled -> overall by
rule 13's `k = 1.011` calibration and `/6`.  `k = 0.899` (162's all-arm fitted scale) is
reported as a sensitivity column, never as the headline.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
FOLDS = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
SPAN = 0.95228 - 0.50020
K13, K899 = 1.011, 0.899
SIGMA = 0.001494          # 162's leave-one-arm-out RMS on scaled `pds`

LOO = ["noMassFix", "noWhiten", "noPM", "noCD4", "cd4_a14", "alpha_v11", "weights_4_1_1",
       "weights_2_1_1", "noGate", "gate5", "gate20", "noPrune", "gamma_050", "gamma_080",
       "a19ridge"]
AOI = ["BARE+whiten", "BARE+pm", "BARE+cd4", "BARE+massfix", "BARE+whiten_unpruned"]


def main() -> None:
    rows = json.loads(open(os.path.join(WORK, "pds200.json")).read())
    R = {(r["fold"], r["arm"], r["space"]): r for r in rows}
    arms = sorted({r["arm"] for r in rows})

    def chal(fold, arm, field="pds"):
        vs = [R[(fold, arm, f"chal{c}/delivered")][field] for c in "ABC"
              if (fold, arm, f"chal{c}/delivered") in R]
        return float(np.mean(vs)) if len(vs) == 3 else None

    def delta(arm, base):
        per = {}
        for f in FOLDS:
            a, b = chal(f, arm), chal(f, base)
            if a is not None and b is not None:
                per[f] = a - b
        if len(per) == 0:
            return None
        v = np.asarray(list(per.values()))
        return dict(per=per, mean=float(v.mean()), sd=float(v.std()),
                    n=int(v.size), scaled=float(v.mean()) / SPAN)

    print("=" * 118)
    print("SHIP absolute `pds`, challenge frame, per fold  (162's pre-committed predictor)")
    print("=" * 118)
    print(f"  {'arm':<22}" + "".join(f"{f:>16}" for f in FOLDS) + f"{'4-fold mean':>14}")
    for a in ("SHIP", "SHIP-dup", "BARE"):
        v = [chal(f, a) for f in FOLDS]
        print(f"  {a:<22}" + "".join(f"{(x if x is not None else float('nan')):16.6f}"
                                     for x in v)
              + f"{np.mean([x for x in v if x is not None]):14.6f}")
    d = [abs(chal(f, "SHIP") - chal(f, "SHIP-dup")) for f in FOLDS]
    print(f"  control: SHIP vs SHIP-dup (field-identical configs) worst |d| = {max(d):.2e}")

    def block(title, names, base):
        print()
        print("=" * 118)
        print(title)
        print("=" * 118)
        print(f"  {'ablation':<22}{'raw d':>10}{'fold sd':>10}{'scaled':>11}"
              f"{'/sigma':>8}{'overall k=1.011':>17}{'k=0.899':>11}   per-fold signs")
        out = {}
        for a in names:
            if a not in arms:
                continue
            r = delta(a, base)
            if r is None:
                continue
            signs = "".join("+" if r["per"][f] > 0 else ("-" if r["per"][f] < 0 else "0")
                            for f in FOLDS if f in r["per"])
            ov13 = r["scaled"] * K13 / 6
            ov899 = r["scaled"] * K899 / 6
            out[a] = dict(**r, overall_k1011=ov13, overall_k0899=ov899, signs=signs)
            note = "" if r["n"] == 4 else f"  ({r['n']} fold only)"
            print(f"  {a:<22}{r['mean']:>+10.6f}{r['sd']:>10.6f}{r['scaled']:>+11.6f}"
                  f"{r['scaled'] / SIGMA:>8.1f}{ov13:>+17.6f}{ov899:>+11.6f}   {signs}{note}")
        return out

    loo = block("LEAVE ONE OUT of the shipped stack   (negative = the component EARNS its "
                "place on `pds`)", LOO, "SHIP")
    aoi = block("ADD ONE IN to a bare consensus   (positive = the component helps on its own)",
                AOI, "BARE")

    print()
    print("=" * 118)
    print("MARGINAL vs SOLO -- the interaction sequential addition could never see")
    print("=" * 118)
    print(f"  {'component':<22}{'LOO (-value in stack)':>24}{'AOI (solo value)':>20}"
          f"{'interaction':>14}")
    PAIRS = [("whitening", "noWhiten", "BARE+whiten"),
             ("PerturbMatch", "noPM", "BARE+pm"),
             ("CD4 residual", "noCD4", "BARE+cd4"),
             ("mass fix", "noMassFix", "BARE+massfix")]
    for label, l, a in PAIRS:
        if l in loo and a in aoi:
            marg = -loo[l]["scaled"]
            solo = aoi[a]["scaled"]
            print(f"  {label:<22}{marg:>+24.6f}{solo:>+20.6f}{marg - solo:>+14.6f}")

    json.dump(dict(loo=loo, aoi=aoi), open(os.path.join(WORK, "report_pds.json"), "w"),
              indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'report_pds.json')}")


if __name__ == "__main__":
    sys.exit(main())
