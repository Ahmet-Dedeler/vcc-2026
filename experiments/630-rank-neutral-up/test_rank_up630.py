"""KA for the builder module (src_rank_up.py) against the H1 pricing run (rnup630.py, v2).

For every dumped seed-<s> block: re-apply the production zs pass, then src_rank_up.apply(..., variant) with the scorer's
tie sums, H1's median control library and the target's kept-axis column. Its output's evaluator stats (sums, jk, p, q,
lfc) must equal rnup630's stats-rnup_s<s>_<variant>.npz row for row, bit for bit, and every block must show 0 call flips,
exact per-gene totals and libraries.

  uv run python experiments/630-rank-neutral-up/test_rank_up630.py <seed> <variant> [<n_targets>]
"""
from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rnup630 as U  # noqa: E402
import src_rank_up as M  # noqa: E402

seed, variant = int(sys.argv[1]), sys.argv[2]
n_t = int(sys.argv[3]) if len(sys.argv) > 3 else 85
r, c, _ = U.E.load_reference()
keep = np.asarray(c.keep, dtype=bool)
G = r["genes"][keep]
kpos = np.cumsum(keep) - 1
gpos = {str(g): i for i, g in enumerate(r["genes"])}
lib_med = float(np.median(r["libraries"]))
tie = M.scorer_tie(c)
assert np.array_equal(tie, c.tie), "module tie sums differ from the evaluator's control block"
ref = np.load(U.WORK / f"stats-rnup_s{seed}_{variant}.npz", allow_pickle=False)
files = sorted(glob.glob(str(HERE / f"work/blocks/{U.TAG.format(seed=seed)}_*.npz")))[:n_t]
spec = U.ZS.parse_spec("dn_k250_zs")
n_ok = 0
for i, f in enumerate(files):
    X0, n_car, target = U.load_block(f)
    _, _, _, _, l0 = U.E.block_stats(X0, c)
    Xz, _ = U.ZS.postpass(X0, n_car, c, l0, spec)
    full = gpos.get(target, -1)
    own = int(kpos[full]) if full >= 0 and keep[full] else -1
    Xu, d = M.apply(Xz, n_car, c, variant, tie=tie, lib_med=lib_med, own_col=own)
    su, ju, pu, qu, lu = U.E.block_stats(Xu, c)
    for nm, a, b in (("sums", su, ref["sums"][i]), ("p", pu, ref["p"][i]), ("q", qu, ref["q"][i]), ("lfc", lu, ref["lfc"][i])):
        if not np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True):
            raise SystemExit(f"FAIL target {i} {target}: module {nm} differs from rnup630's {variant} stats")
    if ju != float(ref["jk"][i]):
        raise SystemExit(f"FAIL target {i} jk")
    if not np.array_equal(Xu.sum(1), Xz.sum(1)) or not np.array_equal(Xu.sum(0), Xz.sum(0)):
        raise SystemExit(f"FAIL target {i}: library or per-gene total moved")
    n_ok += 1
print(f"PASS: src_rank_up.apply reproduces rnup630 {variant} seed {seed} on {n_ok} targets bit for bit "
      f"(module sha {__import__('hashlib').sha256((HERE / 'src_rank_up.py').read_bytes()).hexdigest()[:16]})")
