#!/usr/bin/env python3
"""701 for 711 Amendment 2 (coordinator [2d6f46], 10-16 12:37Z; science [d80a10] writes it BEFORE ku38's read):
condition (i') at the H1 dial f 0.61, with reach701.py's exact arithmetic (same states, regime pools, rows_at, means and known
answer), regime "0.61" instead of "0.87". No new choices. Per K size: the paired-seed delta of reach against K 2,600 (CPM 2200),
its mean, se = sd(ddof=1) / sqrt(3), and pass = all 3 seeds > 0 AND mean > 2 * se.
K 3,400 = CPM 3000 (ku30), K 4,200 = CPM 3800 (ku38); K = 400 + CPM.
usage: uv run python experiments/701-kup-callset/reach701_f061.py   -> work/reach701_f061.json
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

REG = "0.61"
SIZES = {3400: "cpm-n3000", 4200: "cpm-n3800"}


def main():
    states = {s: L.load_state(s) for s in (0, 1, 2)}
    for s in (0, 1, 2):          # reach701's known answer, unchanged: the base's members from stats == its metrics json at f 1.0
        A = P7.load("base", s, False)
        R = P5.rows_at(states[s], A, regime_pool(states[s], "1.0")[0])
        got = [float(np.nanmean(R[m])) for m in range(4)]
        assert all(abs(g - A["raw"][k]) < 1e-12 for g, k in zip(got, ("fid", "jac", "reach", "nmae"))), got
    pools = {s: regime_pool(states[s], REG)[0] for s in (0, 1, 2)}
    reach = lambda mode, s: float(means(P5.rows_at(states[s], P7.load(mode, s), pools[s]))[2])
    base = {s: reach("base", s) for s in (0, 1, 2)}
    out = {"regime": 0.61, "statistic": "reach (raw) at H1 regime f 0.61, paired by seed, vs K 2,600 (CPM 2200)",
           "code": "experiments/701-kup-callset/reach701_f061.py (reach701.py's arithmetic, regime 0.61)",
           "known_answer": "base members from stats == metrics json on seeds 0/1/2 (reach701's KA) PASS", "per_size": {}}
    for k, mode in SIZES.items():
        d = np.array([reach(mode, s) - base[s] for s in (0, 1, 2)])
        se = float(d.std(ddof=1) / np.sqrt(3))
        ok = bool((d > 0).all() and d.mean() > 2 * se)
        out["per_size"][str(k)] = {"cpm": int(mode.split("-n")[1]), "dreach_per_seed": [float(v) for v in d],
                                   "mean": float(d.mean()), "se": se, "pass": ok}
        print(f"K {k} (CPM {mode.split('-n')[1]}): dreach per seed {[round(float(v), 5) for v in d]} mean {d.mean():+.5f} "
              f"se {se:.5f} -> {'PASS' if ok else 'FAIL'}", flush=True)
    (HERE / "work" / "reach701_f061.json").write_text(json.dumps(out, indent=1))
    print("wrote work/reach701_f061.json")


if __name__ == "__main__":
    main()
