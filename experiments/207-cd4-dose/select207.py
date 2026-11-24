"""207 Phase 4 falsifier 3 — the dose must not be selected on the surface it is scored on.

Leave-one-FOLD-out.  For each held fold, choose the dose that maximises the priced total on
the OTHER THREE folds only, then report that dose's priced total on the held fold.  The
selection frame and the reporting frame never overlap.  The honest headline is the mean of
the four held-fold prices, not the in-sample maximum.

    uv run python select207.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from price207 import (CALIB, FOLDS, OURS_COS, SPAN_NMAE, SPAN_PDS, W_HI,  # noqa: E402
                      mse_at, score_mse)

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
BASE = "f0.1"


def main() -> None:
    by = {}
    for n in ("sweep207.json", "library207.json", "form207.json", "libnull207.json"):
        p = WORK / n
        if p.exists():
            by.update({(r["fold"], r["arm"]): r for r in json.loads(p.read_text())})
    k_cos = json.loads((WORK / "price207.json").read_text())["k_cos"]

    def price(arm, folds, k=1.0):
        """The rule-13 composite, restricted to a subset of folds."""
        d_pds = float(np.mean([by[(f, arm)]["pds"] - by[(f, BASE)]["pds"] for f in folds]))
        rel = float(np.mean([by[(f, arm)]["cos_agg"] for f in folds])) / float(
            np.mean([by[(f, BASE)]["cos_agg"] for f in folds])) - 1.0
        d_qw = float(np.mean([by[(f, arm)]["a15_sgn_w"] - by[(f, BASE)]["a15_sgn_w"]
                              for f in folds]))
        o_pds = d_pds / SPAN_PDS / 6.0 * CALIB
        o_mse = (score_mse(mse_at(OURS_COS * (1 + k * rel)))
                 - score_mse(mse_at(OURS_COS))) / 6.0
        o_nmae = (2.0 * W_HI * d_qw) / SPAN_NMAE / 6.0
        return o_pds + o_mse + o_nmae

    families = {
        "residual (f)": [a for (f, a) in by if f == FOLDS[0] and a.startswith("f")
                         and a[1:].replace(".", "").replace("-", "").isdigit()],
        "library (w) @a*": [a for (f, a) in by if f == FOLDS[0]
                            and a.startswith("lib") and a.endswith("@a*")
                            and "+res" not in a and "shuf" not in a and "kolf" not in a],
        "composed (w,f) @a*": [a for (f, a) in by if f == FOLDS[0]
                               and a.startswith("lib") and a.endswith("@a*")
                               and "+res" in a],
        "whitened residual": [a for (f, a) in by if f == FOLDS[0] and a.startswith("wres")],
    }
    out = {}
    for fam, arms in families.items():
        arms = sorted({a for a in arms if all((f, a) in by for f in FOLDS)})
        if not arms:
            continue
        print(f"\n{fam}  ({len(arms)} arms)")
        print(f"{'held fold':<18}{'chosen on the other 3':<26}{'in-sample':>11}"
              f"{'HELD-FOLD price':>17}")
        held_prices, chosen = [], []
        for held in FOLDS:
            other = [f for f in FOLDS if f != held]
            best = max(arms, key=lambda a: price(a, other))
            ins = price(best, other)
            hp = price(best, [held])
            held_prices.append(hp)
            chosen.append(best)
            print(f"{held:<18}{best:<26}{ins:>+11.5f}{hp:>+17.5f}")
        ins_all = max(arms, key=lambda a: price(a, FOLDS))
        print(f"{'MEAN of held-fold prices':<44}{'':>11}{np.mean(held_prices):>+17.5f}"
              f"   (in-sample argmax over all 4: {ins_all} at "
              f"{price(ins_all, FOLDS):+.5f})")
        out[fam] = dict(chosen=chosen, held_prices=held_prices,
                        mean_held=float(np.mean(held_prices)),
                        in_sample_argmax=ins_all, in_sample=price(ins_all, FOLDS),
                        mean_held_cal=float(np.mean(
                            [price(c, [h], k_cos) for c, h in zip(chosen, FOLDS)])))
    (WORK / "select207.json").write_text(json.dumps(out, indent=2))
    print(f"\nwrote {WORK/'select207.json'}")


if __name__ == "__main__":
    main()
