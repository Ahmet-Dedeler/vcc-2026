#!/usr/bin/env python3
"""393 Task 3: price Delta-raw-mse for a separate shared scale c.

Numerator at (c, s=1), per context (187 identity + shared/centred split):
    num_c(c) - num_c(1) = (c^2-1)*S_c - 2*(c-1)*Xs_c
    Delta raw = mean_c Delta num_c / sumD   (sumD = 0.229141, common)
S_c = P||dbar||^2 from decomp393.json; Xs_c bracketed in {0, S/2, S}.
Also prints the score-space conversion (187's score()).
Writes only into this experiment's work/.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"

CTX = ("A", "B", "C")
SUMD = 0.229141
B_ANCHOR, SPAN = 0.987862879, 0.950162566


def score(raw_c):
    return float(np.mean([min(max((B_ANCHOR - r) / SPAN, 0.0), 1.0) for r in raw_c]))


def main() -> None:
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))
    dec = json.loads((WORK / "decomp393.json").read_text())
    S = {c: dec["per_context"][c]["P_dbar2"] for c in CTX}
    P(f"S_c = P||dbar||^2: { {c: round(S[c], 6) for c in CTX} }  "
      f"mean {np.mean(list(S.values())):.6f}; sumD {SUMD}")
    P("Delta raw(c) = mean_c [(c^2-1) S_c - 2(c-1) Xs_c]/sumD, "
      "Xs in {0, S/2, S} (c* in {0, .5, 1})")
    grid = [0, 0.25, 0.5, 0.75, 1, 1.25]
    P(f"\n{'Xs/S':>6}" + "".join(f"{'c=' + str(c):>12}" for c in grid))
    res = {}
    for frac in (0.0, 0.5, 1.0):
        row = []
        for c in grid:
            d = float(np.mean([(c * c - 1) * S[x] - 2 * (c - 1) * frac * S[x] for x in CTX]) / SUMD)
            row.append(d)
        res[str(frac)] = row
        P(f"{frac:>6}" + "".join(f"{v:>+12.6f}" for v in row))
    best = min(res["0.0"])
    P(f"\nbracket best case (Xs=0, c=0): Delta raw {best:+.6f}  "
      f"(|.| vs kill 0.01: {'KILL' if abs(best) < 0.01 else 'ALIVE'})")
    P(f"in score units that is {abs(best) / SPAN:+.6f} scaled ({abs(best) / SPAN / 6:+.6f} overall)")
    (WORK / "price393.json").write_text(json.dumps(
        dict(S=S, sumD=SUMD, grid=grid, delta_raw=res, best=best), indent=1))
    (WORK / "price393.log").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
