"""518 — is the strength key a SUBSTITUTE for the LODO key or a complement? Exact PLOLO rows (score518), vs 446bj2's bulk:
  nokey          w = 1 (the unkeyed stack: what the LODO key buys, with the sign flipped)
  strength only  w = exp(theta z), theta 1 / 1.5 / 2 (no LODO key at all)
  LODO + strength (skey518's best, for reference)

    uv run python experiments/518-bulk-mse-gap/subst518.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402

import numpy as np  # noqa: E402

from score518 import Fold  # noqa: E402
from skey518 import zfeat  # noqa: E402

FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")


def main():
    out = {}
    for f in FOLDS:
        F = Fold(f)
        b, _ = F.score()
        rec = {}
        z = zfeat(F, "meta")
        arms = {"nokey": np.ones(F.n)} | {f"strength only th{t:g}": np.exp(t * z) for t in (1.0, 1.5, 2.0)} | \
               {"LODO + strength th1.25": F.w_base * np.exp(1.25 * z)}
        for nm, w in arms.items():
            s, _ = F.score(w=w)
            rec[nm] = dict(r=s["cos_agg"] / b["cos_agg"] - 1, dpds=s["pds"] - b["pds"])
        out[f] = rec
        print(f"{f:14s} " + " | ".join(f"{nm}: {100 * v['r']:+.2f} %" for nm, v in rec.items()), flush=True)
    names = list(out[FOLDS[0]])
    print("cancer mean: " + " | ".join(f"{nm}: {100 * np.mean([out[f][nm]['r'] for f in FOLDS[:3]]):+.2f} %" for nm in names))
    (HERE / "work/subst518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
