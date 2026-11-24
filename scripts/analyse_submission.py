"""Read a submission's board result and back out what it says about the harness.

Two things are worth extracting the moment a score lands, because both feed the next decision
and neither is obvious from the number itself.

**The LODO -> board slope.** Every configuration choice is made on a harness whose numbers are
not leaderboard numbers. The translation currently rests on two calibration points (the v1 and
v3 families) that are 0.010 apart in LODO and imply slopes anywhere from 0.33 to 0.86 -- and
the change just submitted is six times larger than the gap they span. A third point, far from
the first two, collapses that uncertainty.

**The implied transfer cosine.** `expr_mse ~ r^2 - 2r*cos + 1` where `r` is the ratio of
emitted to true effect magnitude. With `r` known from the build (alpha times the effect scale,
against the true magnitude implied by the previous submission), the board's raw `mse` inverts
to a cosine -- the one quantity that says whether `mse` will ever be reachable, and the one the
local harnesses cannot measure against the real contexts.

    uv run python scripts/analyse_submission.py --lodo-old 0.7378 --lodo-new 0.8078
"""
from __future__ import annotations

import argparse
import glob
import json

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--team", default="Ahmet Dedeler")
    ap.add_argument("--lodo-new", type=float, required=True,
                    help="macro PDS of the configuration just submitted")
    ap.add_argument("--lodo-old", type=float, default=0.7378,
                    help="macro PDS of the configuration the board currently shows")
    ap.add_argument("--r-assumed", type=float, default=0.92,
                    help="emitted / true effect magnitude, from the previous submission's mse")
    args = ap.parse_args()

    f = sorted(glob.glob("data/leaderboard/lb_*.json"))[-1]
    e = json.load(open(f))["entries"]
    e.sort(key=lambda x: -x["score_avg"])
    mine = next((r for r in e if args.team in r["team_name"]), None)
    if mine is None:
        raise SystemExit(f"{args.team} not in {f}")
    rank = e.index(mine) + 1
    print(f"{f}\n  rank {rank}/{len(e)}  overall {mine['score_avg']:.4f}  "
          f"model {mine.get('model_name')!r}")
    for k in ("pds", "mse", "nmae", "fid", "reach", "jac"):
        print(f"    {k:6} scaled {mine['score_' + k]:+.4f}")

    board_new = mine["pds_cosine"]
    print(f"\n  raw pds {board_new:.4f}")
    print(f"  LODO {args.lodo_old:.4f} -> {args.lodo_new:.4f} "
          f"(delta {args.lodo_new - args.lodo_old:+.4f})")
    print("  board delta depends on the previous entry's raw pds -- pass it in by hand and "
          "divide.")

    mse = mine["expr_mse_unbiased_capped_norm"]
    r = args.r_assumed
    cos = (r ** 2 + 1 - mse) / (2 * r)
    print(f"\n  raw mse {mse:.4f} at assumed r={r} implies cos = {cos:.4f}")
    print(f"  min mse at that cosine = {1 - cos ** 2:.4f}; official baseline b = 0.989")
    print(f"  scaled mse becomes positive only above cos = 0.105 -> "
          f"{'REACHABLE' if cos > 0.105 else 'still out of reach'}")


if __name__ == "__main__":
    main()
