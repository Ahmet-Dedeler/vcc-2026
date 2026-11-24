#!/usr/bin/env python3
"""389 smoke: kernel equivalence + single-target K=1 bit-identity + K=3 run.

Fast checks before any full panel. Writes nothing (no WORK outputs).
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ.setdefault(
    "NUMBA_CACHE_DIR",
    str(Path(__file__).resolve().parent / "work" / "numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run323 as R  # noqa: E402
import evaluate323 as E  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import emit_counts as EMIT_COUNTS  # noqa: E402
import k_emit389 as K  # noqa: E402

CLASSES3 = ((1 / 3, 16078.75), (1 / 3, 20190.5), (1 / 3, 25330.5))


def target_inputs(i, r, emit_c, builder, effects):
    target = r["targets"][i]
    profile = builder.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
    spread = builder.emission_spread_vector(
        effects[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
        keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
        spread_lo=0., spread_hi=.60)
    gate = builder.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
    key = np.where(gate, effects[i], 0.)
    return target, profile, spread, key


def main():
    t0 = time.monotonic()
    r, c, emit_c = E.load_reference()
    with R.model_env("uniform"):
        builder = R.load_module("builder389smoke",
                                ROOT / "scripts/make_submission_pm_cd4.py")
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        effects = R.map_axis(z["cap15"], z["genes"], r["genes"])
    print(f"loaded ({time.monotonic() - t0:.1f}s)", flush=True)

    w = K.check_kernel_equivalence(emit_c.cs, emit_c.tie)
    print(f"kernel equivalence max|dz| = {w:.3e}", flush=True)

    pick_rng, _ = builder.make_rngs(0)
    picks0 = pick_rng.choice(len(r["libraries"]), 400, replace=False)
    target, profile, spread, key = target_inputs(0, r, emit_c, builder, effects)
    cs = rc_emit.CallSet(400, 2200, 1.9)
    kw = dict(spread=spread, overdispersion=4., control=emit_c, key=key,
              jitter_seed=1000, call_set_first=True, call_set=cs,
              emit_counts=EMIT_COUNTS)
    b1, d1 = rc_emit.rc_hybrid_block(
        profile, r["libraries"][picks0], .1, np.random.default_rng(7), **kw)
    bk, dk = K.rc_hybrid_block_k(
        profile, r["libraries"][picks0], .1, np.random.default_rng(7),
        classes=None, lfc_target=None, **kw)
    print("K=1 bit-identical:", bool(np.array_equal(b1, bk)),
          "diagKeys:", sorted(set(d1) ^ set(dk)) or "same", flush=True)
    assert np.array_equal(b1, bk)

    b3, d3 = K.rc_hybrid_block_k(
        profile, r["libraries"][picks0], .1, np.random.default_rng(7),
        classes=list(CLASSES3), lfc_target=None, **kw)
    print("K=3 totals identical:", bool(np.array_equal(b3.sum(0), b1.sum(0))),
          "cellTotalsSum:", int(b3.sum(1).sum()), "S:", int(b1.sum()),
          flush=True)
    assert np.array_equal(b3.sum(0), b1.sum(0))
    cc, Lk, Tc, nk = K.class_cells(400, list(CLASSES3), int(b1.sum()))
    print("classes:", Lk.tolist(), nk.tolist(), "perCell_T:", int(b1.sum() // 400),
          flush=True)
    assert np.array_equal(np.sort(b3.sum(1)), np.sort(Tc))
    z1 = rc_emit.wilcoxon_z(
        emit_c.cs, emit_c.tie,
        np.ascontiguousarray(
            np.log1p(b1 / (b1.sum(1, keepdims=True) / 1e6))[:, emit_c.keep].T))
    z3 = rc_emit.wilcoxon_z(
        emit_c.cs, emit_c.tie,
        np.ascontiguousarray(
            np.log1p(b3 / (b3.sum(1, keepdims=True) / 1e6))[:, emit_c.keep].T))
    dz = np.abs(np.abs(z3) - np.abs(z1))
    print(f"|z| drift K=3/plain target0: median {np.median(dz):.3f} "
          f"p90 {np.percentile(dz, 90):.3f} max {dz.max():.3f}", flush=True)
    print(f"smoke OK ({time.monotonic() - t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
