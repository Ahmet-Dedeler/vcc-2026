"""518 — the magnitude pot that REMAINS: the per-target oracle on 446bj2's fold bulk (cap 3), on cap 6 (495), and on
cap 6 + the strength key (exact rows via score518), per frame (context), and how concentrated the remaining loss is
(the smallest number of targets carrying half of it).

    uv run python experiments/518-bulk-mse-gap/remain518.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402

import numpy as np  # noqa: E402

from folds518 import frame_decomp  # noqa: E402
from score518 import Fold  # noqa: E402  (puts 420 on the path)
from lib420 import weights  # noqa: E402
from skey518 import zfeat  # noqa: E402

WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")


def main():
    out = {}
    for f in FOLDS:
        F = Fold(f)
        d = (F.truth * F.truth).sum(1)
        w6 = weights(F.kv1, F.carried, 0.5, 6.0)[0]
        arms = {"446bj2 (cap 3)": dict(), "cap 6 (495)": dict(w=w6), "cap 6 + strength th0.5": dict(w=w6 * np.exp(0.5 * zfeat(F, "meta")))}
        rec = {}
        for nm, kw in arms.items():
            _, arr = F.score(**kw)
            x, a = arr["x"], arr["a"]
            per = {}
            for j, c in enumerate("ABC"):
                fd = frame_decomp(x[j], a[j], d)
                loss = np.sort(fd["loss"])[::-1]
                half = int(np.searchsorted(np.cumsum(loss) / loss.sum(), 0.5) + 1)
                per[c] = dict(r_mag=fd["r_mag"], r_neg=fd["r_neg"], n_half=half, n=int(loss.size))
            rec[nm] = dict(frames=per, r_mag=float(np.mean([per[c]["r_mag"] for c in "ABC"])))
        out[f] = rec
        print(f"{f:14s} " + " | ".join(f"{nm}: oracle r {100 * v['r_mag']:+.1f} % (A {100 * v['frames']['A']['r_mag']:.1f} B {100 * v['frames']['B']['r_mag']:.1f} "
                                       f"C {100 * v['frames']['C']['r_mag']:.1f}; half the loss in {v['frames']['A']['n_half']}/{v['frames']['A']['n']} targets)"
                                       for nm, v in rec.items()))
    (WORK / "remain518.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
