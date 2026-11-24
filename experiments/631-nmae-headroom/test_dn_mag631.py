#!/usr/bin/env python3
"""KA for the DN builder module (src_dn_mag631.py) against the H1 pricing run (rnup630.py with shipping semantics).
For every seed-<s> block: re-apply the production zs pass; for UPDN also 630's module (src_rank_up.apply, b_k100_q20r3p10);
then src_dn_mag631.apply(..., "d_m200") with the scorer's tie sums and the target's kept-axis column. The output's evaluator
stats (sums, jk, p, q, lfc) must equal rnup630's stats-rnup_s<s>_<variant>.npz bit for bit, with totals and libraries exact.
  uv run python experiments/631-nmae-headroom/test_dn_mag631.py <seed> <dn|updn> [<n_targets>]
"""
from __future__ import annotations

import glob
import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/630-rank-neutral-up"))
sys.path.insert(0, str(HERE))
import rnup630 as U  # noqa: E402
import src_rank_up as UP  # noqa: E402
import src_dn_mag631 as DN  # noqa: E402

seed, mode = int(sys.argv[1]), sys.argv[2]
n_t = int(sys.argv[3]) if len(sys.argv) > 3 else 85
variant = {"dn": "b_k000_q20r3p10_D200", "updn": "b_k100_q20r3p10_D200"}[mode]
r, c, _ = U.E.load_reference()
keep = np.asarray(c.keep, dtype=bool)
kpos = np.cumsum(keep) - 1
gpos = {str(g): i for i, g in enumerate(r["genes"])}
lib_med = float(np.median(r["libraries"]))
tie = DN.scorer_tie(c)
assert np.array_equal(tie, c.tie)
ref = np.load(U.WORK / f"stats-rnup_s{seed}_{variant}.npz", allow_pickle=False)
files = sorted(glob.glob(str(ROOT / f"experiments/630-rank-neutral-up/work/blocks/{U.TAG.format(seed=seed)}_*.npz")))[:n_t]
zspec = U.ZS.parse_spec("dn_k250_zs")
n_ok = 0
for i, f in enumerate(files):
    X0, n_car, target = U.load_block(f)
    _, _, _, _, l0 = U.E.block_stats(X0, c)
    Xz, _ = U.ZS.postpass(X0, n_car, c, l0, zspec)
    full = gpos.get(target, -1)
    own = int(kpos[full]) if full >= 0 and keep[full] else -1
    Xb = Xz
    if mode == "updn":
        Xb, _ = UP.apply(Xz, n_car, c, "b_k100_q20r3p10", tie=tie, lib_med=lib_med, own_col=own)
    Xd, d = DN.apply(Xb, n_car, c, "d_m200", tie=tie, own_col=own)
    if d["dn_invariant_failed"] or d["dn_guard_failed"]:
        raise SystemExit(f"FAIL target {i} {target}: a fail-safe fired {d}")
    su, ju, pu, qu, lu = U.E.block_stats(Xd, c)
    for nm, a, b in (("sums", su, ref["sums"][i]), ("p", pu, ref["p"][i]), ("q", qu, ref["q"][i]), ("lfc", lu, ref["lfc"][i])):
        if not np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True):
            raise SystemExit(f"FAIL target {i} {target}: module {nm} differs from rnup630's {variant} stats")
    if ju != float(ref["jk"][i]):
        raise SystemExit(f"FAIL target {i} jk")
    if not np.array_equal(Xd.sum(1), Xz.sum(1)) or not np.array_equal(Xd.sum(0), Xz.sum(0)):
        raise SystemExit(f"FAIL target {i}: library or per-gene total moved")
    n_ok += 1
print(f"PASS: src_dn_mag631.apply ({mode}) reproduces rnup630 {variant} seed {seed} on {n_ok} targets bit for bit "
      f"(module sha {hashlib.sha256((HERE / 'src_dn_mag631.py').read_bytes()).hexdigest()[:16]})")
