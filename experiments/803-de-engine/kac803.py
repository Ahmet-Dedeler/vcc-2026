#!/usr/bin/env python3
"""803 KA (c): the engine's member deltas for 618's gate pairs vs 618's gate (gate618.json: cell-eval2's raw members at f 1,
mean of seeds 0-2) and vs 619's score619.json (528's rows_at at f 1 / 0.5 / 0.35 / cmp, per seed). On the Mac's own stats
these are equal by 619's KA; here the stats are the Linux re-runs, so the residual is the platform difference (§3).
    uv run python experiments/803-de-engine/kac803.py <out dir with h1tv/...>  -> <out dir>/kac803.json
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REG = {"1.0": "f1", "0.5": "f0.5", "0.35": "f0.35", "cmp": "cmp"}
MEM = ("fid", "jac", "reach", "nmae")
PAIRS = {"G1": ("zs", "g1pre", "zs (446bj2z - 446bj2)"), "G2": ("mde", "nomde", "coarse MDE (433a - 417e)"),
         "G3": ("cp", "zs", "cp call set (446bj2zcp - 446bj2z), WORN on 618")}


def main():
    out = Path(sys.argv[1])
    gate = json.loads((ROOT / "experiments/618-h1-trainval-surface/work/gate618.json").read_text())
    s619 = json.loads((ROOT / "experiments/619-board-multipliers-618/work/score619.json").read_text())
    rec, worst = {}, 0.0
    for g, (arm, base, name619) in PAIRS.items():
        per_seed = []
        for s in (0, 1, 2):
            fa, fb = out / "h1tv" / arm / f"s{s}/members.npz", out / "h1tv" / base / f"s{s}/members.npz"
            if not (fa.exists() and fb.exists()):
                continue
            a, b = np.load(fa), np.load(fb)
            ce = {m: float(json.loads(str(a["ce_raw"]))[m] - json.loads(str(b["ce_raw"]))[m]) for m in MEM}
            lib = {r: {m: float(np.nanmean(a[f"{m}_{k}"]) - np.nanmean(b[f"{m}_{k}"])) for m in MEM} for r, k in REG.items()}
            per_seed.append((s, ce, lib))
        if len(per_seed) < 3:
            rec[g] = dict(seeds=len(per_seed), note="incomplete")
            continue
        r = dict(seeds=3, vs_gate618={}, vs_score619={})
        for m in MEM:
            mine = float(np.mean([x[1][m] for x in per_seed]))
            r["vs_gate618"][m] = dict(engine=mine, gate618=gate[g][m], diff=mine - gate[g][m])
            worst = max(worst, abs(mine - gate[g][m]))
        for rg in REG:
            for m in MEM:
                ref = [x[m] for x in s619[name619]["regimes"][rg]["per_seed"]]
                mine = [x[2][rg][m] for x in per_seed]
                d = float(np.max(np.abs(np.array(mine) - np.array(ref))))
                r["vs_score619"][f"{rg}|{m}"] = dict(engine_mean=float(np.mean(mine)), score619_mean=float(np.mean(ref)), max_seed_diff=d)
                worst = max(worst, d)
        rec[g] = r
    rec["worst_abs_diff"] = worst
    rec["pass_2e-5"] = worst < 2e-5
    dest = Path(sys.argv[2]) if len(sys.argv) > 2 else out / "kac803.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
