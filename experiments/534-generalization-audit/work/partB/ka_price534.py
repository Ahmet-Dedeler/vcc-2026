"""534 Part B -- KA-P on efp: my pricer (price534.Panel on my comp534/full534 outputs) == 438's own pricer (price438.Panel
on 438's stored efp jsons) for every step both can price, per fold, levelled and unlevelled, k = 1 and 1.407: |d| <= 1e-12.
Also: the H1 x2 half (M1 -> M2, levelled k = 1, cancer mean) == 438 section 3.1's efp value from summary438.json.
Light (reads jsons). Run after the efp known-answer runs; aborts on failure.

    uv run python ka_price534.py          -> work/partB/efp/ka_price534.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402

import numpy as np  # noqa: E402

import designs534 as D  # noqa: E402
import price534 as P5  # noqa: E402

sys.path.insert(0, str(P5.ROOT / "experiments/438-october-stack-repricing"))
import price438 as P8  # noqa: E402  (reads 438's work read-only; writes nothing on import)

TOL = 1e-12


def main():
    mine, theirs = P5.Panel("efp"), P8.Panel("efp")
    assert mine.folds == theirs.folds == list(P5.CANCER)
    pairs = [(D.M["M1"], D.M["M0"]), (D.M["M2"], D.M["M1"]), (D.M["M3"], D.M["M2"]), (D.M["C-a"], D.M["M1"]),
             (D.M["C-a"], D.M["M3"]), (D.M["M3"], D.M["M1"])]
    worst, rows = 0.0, {}
    for a, r in pairs:
        for level in (True, False):
            for k in (1.0, P5.K_KEY):
                for f in mine.folds:
                    x = mine.price(f, a, k, r, level)["ov"]
                    y = theirs.price(f, a, k, r, level)["ov"]
                    worst = max(worst, abs(x - y))
                    rows[f"{a} vs {r} | L={level} k={k:.4g} | {f}"] = [x, y, abs(x - y)]
    print(f"KA-P: {len(rows)} (pair, route, k, fold) prices; worst |d| {worst:.1e}")
    assert worst <= TOL, f"KA-P FAILED: {worst}"
    s438 = json.loads((P5.W438 / "summary438.json").read_text())
    key = [t for t in s438 if t.startswith("Stack and halves vs each panel's 417a analogue, LEVELLED")][0]
    want = s438[key]["H1 x2 half"]["efp"]["mean"]
    got = float(np.mean([mine.price(f, D.M["M2"], 1.0, D.M["M1"], True)["ov"] for f in mine.folds]))
    print(f"KA-P: H1 x2 half on efp, levelled k = 1: mine {got:+.12f} vs 438 summary {want:+.12f} (|d| {abs(got - want):.1e})")
    assert abs(got - want) <= TOL
    out = dict(worst=worst, n=len(rows), h1x2_half_efp=dict(mine=got, record=want), rows=rows)
    (P5.WORK / "efp" / "ka_price534.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
