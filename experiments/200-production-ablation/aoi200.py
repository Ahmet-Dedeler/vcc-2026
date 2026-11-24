"""200 — ADD-ONE-IN from a bare consensus, and the interaction it exposes.

Sequential addition can only ever measure a component against the stack that existed BELOW it.
Leave-one-out measures it against the stack that exists ABOVE it too.  Running both on the same
instruments, in the same units, isolates the interaction term -- which for two of these four
components is larger than the component's own solo value.

Priced on `pds` (PLOLO, 4 folds, challenge frame) and `mse` (096 claim route) only.  There is no
add-one-in family on h1-test: a bare-consensus DE emission is three more seeds per rung and the
`reach` sem on this surface (0.003-0.007 raw) is larger than every solo value below.

The BARE consensus is `pds200.py`'s: gamma 0, PerturbMatch 0, no CD4, no mass fix, with the
amplitude (`alpha` 0.982) and the gate (`--min-cpm 10`) held at the shipped value so each rung
moves exactly one thing.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")

PAIRS = [
    # (component, LOO pds key, AOI pds key, LOO mse arm, AOI mse arm)
    ("union whitening (gamma 0.65)", "noWhiten",   "BARE+whiten",  "noWhiten",   "BARE+whiten"),
    ("PerturbMatch blend 0.5",       "noPM",       "BARE+pm",      "noPM",       "BARE+pm"),
    ("post-whitening CD4 (f = 0.10)", "noCD4",     "BARE+cd4",     "noCD4",      "BARE+cd4"),
    ("a23 mass fix (minnorm)",       "noMassFix",  "BARE+massfix", "noMassFix",  "BARE+massfix"),
]


def main() -> None:
    rep = json.loads(open(os.path.join(WORK, "report_pds.json")).read())
    mse = json.loads(open(os.path.join(WORK, "mse200.json")).read())
    G = mse["geometry"]
    bare_score = G["BARE"]["score"]

    print("=" * 118)
    print("ADD-ONE-IN vs LEAVE-ONE-OUT   (overall board units; positive = the component HELPS)")
    print("=" * 118)
    print(f"{'component':<32}{'LOO pds':>10}{'LOO mse':>10}{'LOO sum':>10}"
          f"{'AOI pds':>10}{'AOI mse':>10}{'AOI sum':>10}{'interaction':>13}")
    out = {}
    for label, lp, ap, lm, am in PAIRS:
        loo_pds = -rep["loo"][lp]["overall_k1011"]
        aoi_pds = rep["aoi"][ap]["overall_k1011"]
        # `mse` is a SCORE, so both directions are read off the score itself: the leave-one-out
        # against SHIP's score, the add-one-in against BARE's.
        loo_mse = (G["SHIP"]["score"] - G[lm]["score"]) / 6.0
        aoi_mse = (G[am]["score"] - bare_score) / 6.0
        ls, as_ = loo_pds + loo_mse, aoi_pds + aoi_mse
        out[label] = dict(loo_pds=loo_pds, loo_mse=loo_mse, loo_sum=ls,
                          aoi_pds=aoi_pds, aoi_mse=aoi_mse, aoi_sum=as_,
                          interaction=ls - as_)
        print(f"{label:<32}{loo_pds:>+10.6f}{loo_mse:>+10.6f}{ls:>+10.6f}"
              f"{aoi_pds:>+10.6f}{aoi_mse:>+10.6f}{as_:>+10.6f}{ls - as_:>+13.6f}")

    print()
    print("Reading:  LOO = what removing it from the FULL stack costs.  AOI = what adding it to")
    print("a BARE consensus buys.  interaction = LOO - AOI, i.e. how much of the component's")
    print("value only exists in the presence of the others.")
    print()
    print(f"BARE `mse` score {bare_score:.5f} (raw {G['BARE']['raw']:.5f}) -- the bare consensus")
    print("is already OFF the mse cliff, so every add-one-in rung except whitening reads 0 on")
    print("`mse` by censoring, not by being neutral.  Those three AOI sums are `pds` only.")

    json.dump(out, open(os.path.join(WORK, "aoi200.json"), "w"), indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'aoi200.json')}")


if __name__ == "__main__":
    sys.exit(main())
