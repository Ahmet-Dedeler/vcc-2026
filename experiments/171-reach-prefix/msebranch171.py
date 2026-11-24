"""171 — the `mse` price of each arm at a15's SHIPPED amplitude, on the correct branch.

`claim171.py`'s last table applies the linear form `d raw mse = -d claim / sumD`, which is
164's rule and is correct only BELOW the knee.  The scored numerator is

    raw mse = 1 + ( s^2 A  -  2 s X  +  Kc  -  min(s^2 B, claim) ) / sumD

so once `claim >= s^2 B` the credit saturates and raw `mse` stops responding to claim
ALTOGETHER.  a15 ships `VCC_ALPHA = 0.982` against v11's 1.0796, i.e. s = 0.90960, and
s^2 B = 0.23271 against a15's claim 0.24581 -- a15 is already 5.6% ABOVE the knee.

Consequence, and it is the one that matters for this experiment: an arm that RAISES claim
buys exactly nothing, and an arm that lowers it is free until it crosses 0.23271.  164 read
this as "a15 is on the cliff"; the sharper statement is that a15 has 0.01310 of claim
HEADROOM before `mse` moves at all.

This script prints the exact branch for every policy `claim171.py` measured.
"""
from __future__ import annotations

import re
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
bs = np.load(f"{ROOT}/experiments/096-mse-floor/work/board_solve.npz")
sD, A, B = float(bs["sumD"]), float(bs["v11_A"]), float(bs["v11_B"])
Kc, X, K = float(bs["v11_Kc"]), float(bs["X"]), float(bs["K"])

S15 = 0.982 / 1.0796           # a15's shipped amplitude in v11 units
CLAMP_HI, CLAMP_SPAN = 0.987862879, 0.950162566


def raw_mse(s, claim):
    return 1 + (s * s * A - 2 * s * X + Kc - min(s * s * B, claim)) / sD


def score(raw):
    return min(max((CLAMP_HI - raw) / CLAMP_SPAN, 0.0), 1.0)


def main():
    log = open(f"{ROOT}/experiments/171-reach-prefix/work/claim171.log").read()
    rows = []
    for line in log.splitlines():
        m = re.match(r"^(\S.{0,29}?)\s{2,}(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+"
                     r"(YES|no)\s", line)
        if m:
            rows.append((m.group(1).strip(), float(m.group(5)), m.group(6)))
    if not rows:
        sys.exit("no policy rows parsed from claim171.log")

    knee = S15 * S15 * B
    a15 = dict(rows and [(n, c) for n, c, _ in rows])["a15 AS SHIPPED"]
    print(f"a15 amplitude s = {S15:.5f}   knee s^2 B = {knee:.5f}   "
          f"a15 claim = {a15:.5f}   headroom = {a15 - knee:+.5f}")
    print(f"reference cap K = {K:.5f}   (claim can never exceed it)")
    print(f"a15 raw mse = {raw_mse(S15, a15):.5f}   score_mse = {score(raw_mse(S15, a15)):.5f}")
    print()
    print("HELD AT a15's SHIPPED VCC_ALPHA 0.982 — one-variable arms, correct branch")
    print(f"  {'policy':<30}{'claim':>9}{'branch':>12}{'raw mse':>10}"
          f"{'score_mse':>11}{'d score':>10}{'d OVERALL':>11}")
    base_raw = raw_mse(S15, a15)
    for name, claim, capped in rows:
        r = raw_mse(S15, claim)
        br = "SATURATED" if claim >= knee else "linear"
        print(f"  {name:<30}{claim:>9.5f}{br:>12}{r:>10.5f}{score(r):>11.5f}"
              f"{score(r)-score(base_raw):>+10.5f}{(score(r)-score(base_raw))/6:>+11.6f}")
    print()
    print("  Reading: every arm at or above claim 0.23271 has raw mse EXACTLY a15's, so its")
    print("  `mse` delta is 0.000000 -- not 'small', identically zero, because the metric")
    print("  reads min(s^2 B, claim) and both sides saturate.  claim171.py's linear column")
    print("  over-credits those arms by +0.000551 and is superseded here.")


if __name__ == "__main__":
    main()
