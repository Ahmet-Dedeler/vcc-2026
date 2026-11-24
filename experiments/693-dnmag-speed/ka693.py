#!/usr/bin/env python3
"""693 KA: dn_mag v2 (this dir's dn_mag2.py) == v1 (src/vcc2026/dn_mag.py, 0a8cc560), byte for byte.
On 630's real H1 fixture blocks (test_dn_mag631's surface), for a seed and a mode:
  dn   : zs pass -> DN                      (631's DN-only route)
  updn : zs pass -> 630's UP (b_k100_q20r3p10) -> DN   (the shipped order)
both modules get the identical input block; the output blocks must be np.array_equal and the diagnostics dicts equal.
Also times v1 and v2 per block.
usage: uv run python experiments/693-dnmag-speed/ka693.py <seed> <dn|updn> [n_blocks=85]
"""
from __future__ import annotations

import glob
import hashlib
import importlib.util
import json
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
from vcc2026 import dn_mag as V1  # noqa: E402

spec = importlib.util.spec_from_file_location("dn_mag2", HERE / "dn_mag2.py")
V2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V2)

seed, mode = int(sys.argv[1]), sys.argv[2]
n_b = int(sys.argv[3]) if len(sys.argv) > 3 else 85
r, c, _ = U.E.load_reference()
keep = np.asarray(c.keep, dtype=bool)
kpos = np.cumsum(keep) - 1
gpos = {str(g): i for i, g in enumerate(r["genes"])}
lib_med = float(np.median(r["libraries"]))
tie = V1.scorer_tie(c)
files = sorted(glob.glob(str(ROOT / f"experiments/630-rank-neutral-up/work/blocks/{U.TAG.format(seed=seed)}_*.npz")))[:n_b]
zspec = U.ZS.parse_spec("dn_k250_zs")
t1 = t2 = 0.0
n_ok, fs = 0, 0
for i, f in enumerate(files):
    X0, n_car, target = U.load_block(f)
    _, _, _, _, l0 = U.E.block_stats(X0, c)
    Xb, _ = U.ZS.postpass(X0, n_car, c, l0, zspec)
    full = gpos.get(target, -1)
    own = int(kpos[full]) if full >= 0 and keep[full] else -1
    if mode == "updn":
        Xb, _ = UP.apply(Xb, n_car, c, "b_k100_q20r3p10", tie=tie, lib_med=lib_med, own_col=own)
    a = time.perf_counter()
    X1, d1 = V1.apply(Xb, n_car, c, "d_m200", tie=tie, own_col=own)
    b = time.perf_counter()
    X2, d2 = V2.apply(Xb, n_car, c, "d_m200", tie=tie, own_col=own)
    e = time.perf_counter()
    t1 += b - a
    t2 += e - b
    if X1.dtype != X2.dtype or not np.array_equal(X1, X2):
        raise SystemExit(f"FAIL block {i} {target}: v2's block differs from v1's ({int((X1 != X2).sum())} entries)")
    if d1 != d2:
        raise SystemExit(f"FAIL block {i} {target}: diagnostics differ {d1} vs {d2}")
    fs += int(d1["dn_invariant_failed"] or d1["dn_guard_failed"])
    n_ok += 1
out = {"seed": seed, "mode": mode, "blocks": n_ok, "v1_s": round(t1, 1), "v2_s": round(t2, 1),
       "speedup": round(t1 / max(t2, 1e-9), 2), "failsafe_blocks": fs,
       "v1_sha": hashlib.sha256((ROOT / "src/vcc2026/dn_mag.py").read_bytes()).hexdigest(),
       "v2_sha": hashlib.sha256((HERE / "dn_mag2.py").read_bytes()).hexdigest()}
(HERE / "work").mkdir(exist_ok=True)
json.dump(out, open(HERE / f"work/ka693_s{seed}_{mode}.json", "w"), indent=1)
print(f"PASS: v2 == v1 byte for byte on {n_ok} blocks (seed {seed}, {mode}); v1 {t1:.1f} s, v2 {t2:.1f} s, "
      f"speedup x{t1 / max(t2, 1e-9):.2f}; fail-safe blocks {fs}; v2 sha {out['v2_sha'][:16]}")
