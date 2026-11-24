#!/usr/bin/env python3
"""701 POST HOC (coordinator [2d6f46], 10-16): does the production emitter's REACH rise with the call-set size on H1?
This is condition (i) of the October K rule (711 PREREGISTRATION), fixed before these runs exist.

Runs: 537's real production-emitter H1 harness (run537.py cpm-n<CPM> base <seed> ... --zs dn_k250_zs), sizes CPM 2200
(the base, `cpm`, K 2,600), 2800, 3000, 3400, 3800 (K 3,200 / 3,400 / 3,800 / 4,200), seeds 0 / 1 / 2.
Members per seed and regime from the saved stats with price537's arithmetic (P5.rows_at on 521's regime pools; the KA:
the base's members recomputed from stats == its metrics json, as price537 asserts).
Statistic: per seed, the least-squares slope of reach on CPM (thousands) across the five sizes at the board's size regime
f 0.87 (FINDINGS 0b: the board sits near H1's f 1 regime for SIZE changes, dial f ≈ 0.87); f 1.0 is reported beside it.
Condition (i) HOLDS iff, at f 0.87, all 3 per-seed slopes are > 0 AND the mean slope > 2 × sd(slopes)/sqrt(3) (paired-seed
significance), AND the mean paired Δreach at CPM 3800 ≥ the one at CPM 3000 (still rising at 711's size).
usage: uv run python experiments/701-kup-callset/reach701.py
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
import price537 as P7  # noqa: E402
from policy521 import means  # noqa: E402
from table521 import regime_pool  # noqa: E402

SIZES = (("base", 2200), ("cpm-n2800", 2800), ("cpm-n3000", 3000), ("cpm-n3400", 3400), ("cpm-n3800", 3800))
REGS = ("0.87", "1.0")


def main():
    states = {s: L.load_state(s) for s in (0, 1, 2)}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in (0, 1, 2) for r in REGS}
    for s in (0, 1, 2):          # KA (price537's): the base's members from stats == its metrics json at f 1.0
        A = P7.load("base", s, False)
        R = P5.rows_at(states[s], A, pools[(s, "1.0")])
        got = [float(np.nanmean(R[m])) for m in range(4)]
        assert all(abs(g - A["raw"][k]) < 1e-12 for g, k in zip(got, ("fid", "jac", "reach", "nmae"))), got
    print("KA: base members from stats == metrics json on seeds 0/1/2", flush=True)
    out = {"sizes": [c for _, c in SIZES], "regimes": {}}
    for reg in REGS:
        reach = {}
        for mode, cpm in SIZES:
            for s in (0, 1, 2):
                X = P7.load(mode, s)
                if X is None:
                    raise SystemExit(f"missing H1 run: {mode} seed {s}")
                reach[(cpm, s)] = float(means(P5.rows_at(states[s], X, pools[(s, reg)]))[2])
        x = np.array([c for _, c in SIZES], float) / 1000.0
        slopes = []
        for s in (0, 1, 2):
            y = np.array([reach[(c, s)] for _, c in SIZES])
            slopes.append(float(np.polyfit(x, y, 1)[0]))
        sl = np.array(slopes)
        se = float(sl.std(ddof=1) / np.sqrt(3))
        d = {c: float(np.mean([reach[(c, s)] - reach[(2200, s)] for s in (0, 1, 2)])) for _, c in SIZES}
        holds = bool((sl > 0).all() and sl.mean() > 2 * se and d[3800] >= d[3000])
        out["regimes"][reg] = dict(slopes_per_seed=slopes, slope_mean=float(sl.mean()), slope_se=se,
                                   dreach_vs_2200_mean=d, condition_i=holds,
                                   reach=[[reach[(c, s)] for s in (0, 1, 2)] for _, c in SIZES])
        print(f"f {reg}: reach slope per 1,000 CPM per seed {[round(v, 5) for v in slopes]} mean {sl.mean():+.5f} "
              f"(se {se:.5f}); mean dreach vs CPM 2200: " + ", ".join(f"{c} {d[c]:+.5f}" for _, c in SIZES[1:])
              + f" -> condition (i) {'HOLDS' if holds else 'FAILS'}", flush=True)
    out["verdict_i_at_f087"] = out["regimes"]["0.87"]["condition_i"]
    (HERE / "work" / "reach701.json").write_text(json.dumps(out, indent=1))
    print("CONDITION (i) at f 0.87:", "HOLDS" if out["verdict_i_at_f087"] else "FAILS")


if __name__ == "__main__":
    main()
