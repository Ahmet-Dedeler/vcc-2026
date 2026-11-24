#!/usr/bin/env python3
"""701 EXPLORATORY (science [d80a10], 10-16): reach701.py's arithmetic on six sizes, adding CPM 4600 (K 5,000), so a
possible K 5,000 dose partner of 711 has its held-out price. NOT the October K rule's condition (i), which stays
reach701.py's five sizes and its verdict.
Reports, per regime (f 0.87 and 1.0): per-seed reach at each size, the paired step 3800 -> 4600 per seed with its mean and
seed standard error, and the six-point slope per 1,000 CPM.
usage: uv run python experiments/701-kup-callset/reach701d.py
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

SIZES = (("base", 2200), ("cpm-n2800", 2800), ("cpm-n3000", 3000), ("cpm-n3400", 3400), ("cpm-n3800", 3800),
         ("cpm-n4600", 4600))
REGS = ("0.87", "1.0")


def main():
    states = {s: L.load_state(s) for s in (0, 1, 2)}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in (0, 1, 2) for r in REGS}
    for s in (0, 1, 2):          # the same known answer as reach701: the base's members from stats == its metrics json
        A = P7.load("base", s, False)
        R = P5.rows_at(states[s], A, pools[(s, "1.0")])
        got = [float(np.nanmean(R[m])) for m in range(4)]
        assert all(abs(g - A["raw"][k]) < 1e-12 for g, k in zip(got, ("fid", "jac", "reach", "nmae"))), got
    print("KA: base members from stats == metrics json on seeds 0/1/2 (as reach701)", flush=True)
    out = {"sizes": [c for _, c in SIZES], "regimes": {}, "label": "EXPLORATORY (not condition (i))"}
    for reg in REGS:
        reach = {}
        for mode, cpm in SIZES:
            for s in (0, 1, 2):
                X = P7.load(mode, s)
                if X is None:
                    raise SystemExit(f"missing H1 run: {mode} seed {s}")
                reach[(cpm, s)] = float(means(P5.rows_at(states[s], X, pools[(s, reg)]))[2])
        step = np.array([reach[(4600, s)] - reach[(3800, s)] for s in (0, 1, 2)])
        x = np.array([c for _, c in SIZES], float) / 1000.0
        slopes = [float(np.polyfit(x, np.array([reach[(c, s)] for _, c in SIZES]), 1)[0]) for s in (0, 1, 2)]
        se = float(step.std(ddof=1) / np.sqrt(3))
        out["regimes"][reg] = dict(reach=[[reach[(c, s)] for s in (0, 1, 2)] for _, c in SIZES],
                                   step_3800_to_4600=step.tolist(), step_mean=float(step.mean()), step_se=se,
                                   slope6_per_seed=slopes)
        print(f"f {reg}: step CPM 3800 -> 4600 per seed {[round(v, 5) for v in step]} mean {step.mean():+.5f} (se {se:.5f}); "
              f"six-point slope per 1,000 CPM {[round(v, 5) for v in slopes]}", flush=True)
    (HERE / "work" / "reach701d.json").write_text(json.dumps(out, indent=1))
    print("wrote work/reach701d.json (EXPLORATORY)")


if __name__ == "__main__":
    main()
