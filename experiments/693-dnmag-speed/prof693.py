#!/usr/bin/env python3
"""693 — where does src/vcc2026/dn_mag.py spend its time? Light profile on real H1 blocks (630's fixtures, the KA surface of
631's test_dn_mag631.py), in the shipped order: zs pass -> 630's UP pass (b_k100_q20r3p10) -> dn_mag.apply("d_m200").

Per block: wall time of apply, the number of selected (called DOWN) genes, the inner-loop iterations, ranks_of calls,
de_stats calls and guard rounds; then a cProfile of the slowest block sorted by own time.
usage: uv run python experiments/693-dnmag-speed/prof693.py [seed=0] [n_blocks=8] [--cprofile]
"""
from __future__ import annotations

import cProfile
import glob
import io
import json
import pstats
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
sys.path.insert(0, str(ROOT / "src"))
import rnup630 as U  # noqa: E402
import src_rank_up as UP  # noqa: E402
from vcc2026 import dn_mag as DN  # noqa: E402
from vcc2026 import rank_mag as ZS  # noqa: E402

seed = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 0
n_b = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 8
do_prof = "--cprofile" in sys.argv

# instrumentation (wrappers count calls; they do not change results)
CNT = {"ranks_of": 0, "de_stats": 0, "dn_pass": 0}
_ro, _ds, _dp = ZS.ranks_of, DN.de_stats, DN.dn_pass
def ranks_of(*a, **k):
    CNT["ranks_of"] += 1
    return _ro(*a, **k)
def de_stats(*a, **k):
    CNT["de_stats"] += 1
    return _ds(*a, **k)
def dn_pass(*a, **k):
    CNT["dn_pass"] += 1
    return _dp(*a, **k)
ZS.ranks_of, DN.de_stats, DN.dn_pass = ranks_of, de_stats, dn_pass

r, c, _ = U.E.load_reference()
keep = np.asarray(c.keep, dtype=bool)
kpos = np.cumsum(keep) - 1
gpos = {str(g): i for i, g in enumerate(r["genes"])}
lib_med = float(np.median(r["libraries"]))
tie = DN.scorer_tie(c)
files = sorted(glob.glob(str(ROOT / f"experiments/630-rank-neutral-up/work/blocks/{U.TAG.format(seed=seed)}_*.npz")))[:n_b]
zspec = U.ZS.parse_spec("dn_k250_zs")
rows, prepared = [], []
for f in files:
    X0, n_car, target = U.load_block(f)
    _, _, _, _, l0 = U.E.block_stats(X0, c)
    Xz, _ = U.ZS.postpass(X0, n_car, c, l0, zspec)
    full = gpos.get(target, -1)
    own = int(kpos[full]) if full >= 0 and keep[full] else -1
    Xb, _ = UP.apply(Xz, n_car, c, "b_k100_q20r3p10", tie=tie, lib_med=lib_med, own_col=own)
    prepared.append((target, Xb, n_car, own))
    for k in CNT:
        CNT[k] = 0
    t0 = time.perf_counter()
    t_ds = time.perf_counter()
    _ds(Xb, c, tie)
    t_ds = time.perf_counter() - t_ds
    t0 = time.perf_counter()
    Xd, d = DN.apply(Xb, n_car, c, "d_m200", tie=tie, own_col=own)
    dt = time.perf_counter() - t0
    row = {"target": target, "seconds": round(dt, 2), "de_stats_s": round(t_ds, 3), "n_car": int(n_car),
           "n_bal": int(Xb.shape[0] - n_car), "sel": int(d["dn_sel"]), "moved_genes": int(d["dn_moved"]),
           "counts": int(d["dn_counts"]), "guard_rounds": int(d["dn_guard_rounds"]), **{f"n_{k}": v for k, v in CNT.items()}}
    rows.append(row)
    print(json.dumps(row), flush=True)
tot = sum(x["seconds"] for x in rows)
print(f"\n{len(rows)} blocks: total {tot:.1f} s, mean {tot / max(len(rows), 1):.2f} s/block; "
      f"ranks_of calls/block mean {np.mean([x['n_ranks_of'] for x in rows]):.0f}; "
      f"de_stats {np.mean([x['de_stats_s'] for x in rows]):.3f} s/call")
if do_prof and prepared:
    i = int(np.argmax([x["seconds"] for x in rows]))
    target, Xb, n_car, own = prepared[i]
    pr = cProfile.Profile()
    pr.enable()
    DN.apply(Xb, n_car, c, "d_m200", tie=tie, own_col=own)
    pr.disable()
    s = io.StringIO()
    pstats.Stats(pr, stream=s).sort_stats("tottime").print_stats(14)
    print(f"\ncProfile of the slowest block ({target}):\n" + s.getvalue())
json.dump(rows, open(HERE / f"work/prof693_s{seed}.json", "w"), indent=1)
