#!/usr/bin/env python3
"""662 — which truth-power regime does the board match for CALL-SET changes? 537's H1 real-emitter runs re-scored on a fine
PowerSurface grid (f 0.35 … 1.0 in steps of 0.05, plus 521's cmp), with 537's exact member functions (528's rows_at on
lib521's state; realized order). Read-only: no emission, no cell-eval2, < 2 GB.

Why a finer grid: the K20 board read (10-04) has fid's f 1 SIGN but ~1/6 of f 1's size, and jac's f 0.5 sign. The board sits
between 521's f 1 and f 0.5 for call-set size, a range the four stored regimes do not resolve.

Known answer: at f 1 / 0.5 / 0.35 / cmp every comparison's d_mean equals 537's price537.json to 1e-12 (asserted).
Writes work/regime662.json: {comparison: {f: {d_mean, sd, value}}}.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/662-callset-regime/regime662.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import price528 as P5  # noqa: E402
import price537 as P7  # noqa: E402  (read-only: load(), the stats lookup)
from policy521 import SGN, means  # noqa: E402
from table521 import regime_pool  # noqa: E402

_F = sorted(set(np.round(np.arange(0.35, 1.0001, 0.05), 2)) | set(np.round(np.arange(0.55, 0.7001, 0.01), 2)))
FGRID = tuple(f"{x:.2f}".rstrip("0").rstrip(".") if x != 1 else "1.0" for x in _F)  # 0.05 steps, 0.01 steps over 0.55-0.70
REGS = FGRID + ("cmp",)
MEMBERS = ("fid", "jac", "reach", "nmae")
SEEDS = (0, 1, 2)
COMPARISONS = [  # (label in price537.json, arm, base)
    ("REF cp (528, kept) on 446bj2z", "lcap-kept", "base"),
    ("P2 K 2,000 alone (446bj2zk20's construction)", "cpm-n1600", "base"),
    ("P3 K 1,600 alone (446bj2zk16's construction)", "cpm-n1200", "base"),
    ("P2/3 K 2,000 vs K 1,600 (523's dose pair, K20 - k16)", "cpm-n1600", "cpm-n1200"),
    ("P1 CLEAN4 hook on 446bj2z (= 496a's knob on the zs base)", "cpm+c4h", "base"),
    ("P1 CLEAN4 hook on the cp base", "lcap-kept+c4h", "lcap-kept"),
    ("P1 cp + CLEAN4 hook, total vs 446bj2z", "lcap-kept+c4h", "base"),
    ("P8 CLEAN4 hook attenuated r 0.17 on 446bj2z", "cpm+c4h-r0.17", "base"),
    ("P8 CLEAN4 hook r 0.17 on the cp base", "lcap-kept+c4h-r0.17", "lcap-kept"),
    ("P8 cp + hook r 0.17, total vs 446bj2z", "lcap-kept+c4h-r0.17", "base"),
    ("P9 sign-shuffled CLEAN4 hook on the cp base", "lcap-kept+c4h-s", "lcap-kept"),
    ("P9 cp + sign-shuffled hook, total vs 446bj2z", "lcap-kept+c4h-s", "base"),
]


def main():
    P = json.load(open(ROOT / "experiments/537-arm-pipeline/work/price537.json"))["comparisons"]
    states = {s: L.load_state(s) for s in SEEDS}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in SEEDS for r in REGS}
    nreal = {r: [int(x) for x in pools[(0, r)].sum(1)] for r in REGS}
    cache = {}

    def rows(mode, s, r):
        if (mode, s, r) not in cache:
            A = P7.load(mode, s)
            assert A is not None, (mode, s)
            cache[(mode, s, r)] = P5.rows_at(states[s], A, pools[(s, r)])
        return cache[(mode, s, r)]

    out = {"FGRID": list(FGRID), "n_real_seed0": nreal, "comparisons": {}}
    ka = []
    for label, arm, b0 in COMPARISONS:
        assert P[label]["arm"] == arm and P[label]["base"] == b0, label
        rec = {"arm": arm, "base": b0}
        for r in REGS:
            dms = np.stack([means(rows(arm, s, r)) - means(rows(b0, s, r)) for s in SEEDS])
            dm = dms.mean(0)
            rec[r] = dict(d_mean=dict(zip(MEMBERS, map(float, dm))), sd=dict(zip(MEMBERS, map(float, dms.std(0, ddof=1)))),
                          value=float(np.mean(dms @ SGN)), per_seed=[dict(zip(MEMBERS, map(float, x))) for x in dms])
            if r in P[label]:
                ref = P[label][r]["d_mean"]
                err = max(abs(ref[m] - rec[r]["d_mean"][m]) for m in MEMBERS)
                ka.append((label, r, err))
                assert err < 1e-12, (label, r, err)
        out["comparisons"][label] = rec
        print(label[:60].ljust(60), " | ".join(f"f {r} fid {rec[r]['d_mean']['fid']:+.5f} jac {rec[r]['d_mean']['jac']:+.5f}"
                                               for r in ("1.0", "0.6", "0.5")), flush=True)
    out["KA"] = dict(n=len(ka), max_err=max(e for *_, e in ka))
    print(f"KA: {len(ka)} (comparison, regime) cells equal price537.json, max |err| {out['KA']['max_err']:.2e}")
    (HERE / "work/regime662.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
