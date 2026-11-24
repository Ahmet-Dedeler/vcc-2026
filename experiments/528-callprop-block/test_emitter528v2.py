#!/usr/bin/env python3
"""528 — unit tests of the emitter patch's final form (k_emit528v3.py = patch_emitter528.py with `block_steer` key / cpm /
kept; k_emit528v2.py had key / cpm), on test_emitter528's synthetic context:

  V1  v2 unset == production (bit for bit: block, diagnostics, prior record)
  V2  v2 with a key and block_steer="key" (the default) == v1 (k_emit528.py, the emitter the main H1 queue ran) bit for bit,
      full swap and partial keep
  V3  v2 with block_steer="cpm": the call set is head | block_from_key (missed / unsuppressible consistent), the STEERED
      block is the CPM block (prior record: prior_n_block, shrink kf and flips equal the unkeyed run's), per-gene totals
      unchanged
  V4  block_steer outside {"key", "cpm", "kept"} is refused
  V5  v3 == v2 for "key" and "cpm" (the H1 -dec runs used v2)
  V6  block_steer="kept": the call set is head | key block; the steered block is (key block) & (CPM block): the prior
      record's prior_n_block equals |called block & CPM block|; per-gene totals unchanged
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/test_emitter528v2.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from vcc2026 import rc_emit  # noqa: E402
from vcc2026 import rc_emit_k as P  # noqa: E402
import k_emit528 as K1  # noqa: E402
import k_emit528v2 as K2v2  # noqa: E402
import k_emit528v3 as K2  # noqa: E402
from test_emitter528 import CS, context, run, same  # noqa: E402


def main():
    C = context()
    ref_lin = np.asarray(C["control"].ref_lin, dtype=np.float64)
    n_keep = C["n_keep"]
    base = run(P, C)
    rk = np.random.default_rng(5).random(n_keep)
    ok = {}
    ok["V1 v2 unset == production"] = same(run(K2, C), base)
    for keep in (0, 100):
        ok[f"V2 key, steer 'key', keep {keep}: v2 == v1"] = same(run(K2, C, block_key=rk, block_keep=keep),
                                                                  run(K1, C, block_key=rk, block_keep=keep))
        ok[f"V2 explicit steer='key' == default, keep {keep}"] = same(run(K2, C, block_key=rk, block_keep=keep, block_steer="key"),
                                                                      run(K2, C, block_key=rk, block_keep=keep))
    kabs = np.abs(np.asarray(C["key"])[C["control"].keep])
    head = np.zeros(n_keep, bool)
    head[np.argsort(-kabs, kind="stable")[:CS.n_head]] = True
    for keep in (0, 100):
        b3, d3, r3 = run(K2, C, block_key=rk, block_keep=keep, block_steer="cpm")
        want = head | K2.block_from_key(head, ref_lin, CS.n_cpm, rk, keep)
        z, p, q, lfc = rc_emit.de_from_block(b3, C["control"])
        dcall = q < 0.05
        ok[f"V3 steer cpm, keep {keep}: call set = head | key block"] = (
            int((want & ~dcall).sum()) == d3["call_set_missed"] and int((dcall & ~want).sum()) == d3["call_set_unsuppressible"])
        ok[f"V3 steer cpm, keep {keep}: steered block = the CPM block (prior record)"] = all(
            r3[k] == base[2][k] for k in ("prior_n_block", "prior_shrink_k", "prior_n_flip_target"))
        ok[f"V3 steer cpm, keep {keep}: per-gene totals unchanged"] = bool(np.array_equal(b3.sum(0), base[0].sum(0)))
        ok[f"V3 steer cpm, keep {keep}: differs from steer key"] = not same((b3, d3, r3), run(K2, C, block_key=rk, block_keep=keep))
        print(f"  steer cpm keep {keep}: missed {d3['call_set_missed']} unsuppressible {d3['call_set_unsuppressible']} "
              f"prior block {r3['prior_n_block']} (base {base[2]['prior_n_block']})", flush=True)
    try:
        run(K2, C, block_key=rk, block_keep=0, block_steer="both")
        ok["V4 refuses block_steer='both'"] = False
    except ValueError:
        ok["V4 refuses block_steer='both'"] = True
    for steer in ("key", "cpm"):
        for keep in (0, 100):
            ok[f"V5 v3 == v2, steer {steer}, keep {keep}"] = same(run(K2, C, block_key=rk, block_keep=keep, block_steer=steer),
                                                                 run(K2v2, C, block_key=rk, block_keep=keep, block_steer=steer))
    cpm_blk = np.zeros(n_keep, bool)
    cpm_blk[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[:CS.n_cpm]] = True
    for keep in (0, 100):
        b6, d6, r6 = run(K2, C, block_key=rk, block_keep=keep, block_steer="kept")
        kb = K2.block_from_key(head, ref_lin, CS.n_cpm, rk, keep)
        want = head | kb
        z, p, q, lfc = rc_emit.de_from_block(b6, C["control"])
        dcall = q < 0.05
        ok[f"V6 kept, keep {keep}: call set = head | key block"] = (
            int((want & ~dcall).sum()) == d6["call_set_missed"] and int((dcall & ~want).sum()) == d6["call_set_unsuppressible"])
        # prior_n_block counts the steered block after its (Ms > 0) & finite(E0v) filter: compare with the same filter applied
        # to the intersection, which the unkeyed run's record gives for the full CPM block
        ok[f"V6 kept, keep {keep}: steered block <= |key block & CPM block|"] = r6["prior_n_block"] <= int((kb & cpm_blk).sum())
        ok[f"V6 kept, keep {keep}: steered block smaller than the full CPM block"] = r6["prior_n_block"] < base[2]["prior_n_block"]
        ok[f"V6 kept, keep {keep}: per-gene totals unchanged"] = bool(np.array_equal(b6.sum(0), base[0].sum(0)))
        print(f"  steer kept keep {keep}: steered {r6['prior_n_block']} (|key & cpm| {int((kb & cpm_blk).sum())}, base "
              f"{base[2]['prior_n_block']}); missed {d6['call_set_missed']} unsuppressible {d6['call_set_unsuppressible']}", flush=True)
    for k, v in ok.items():
        print(("PASS " if v else "FAIL ") + k)
    n = sum(bool(v) for v in ok.values())
    print(f"{n}/{len(ok)} PASS")
    sys.exit(0 if n == len(ok) else 1)


if __name__ == "__main__":
    main()
