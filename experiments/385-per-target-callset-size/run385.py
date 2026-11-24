#!/usr/bin/env python3
"""385 harness: known-answer cap15 emit+score in this experiment's own work/.

Same frozen call sequence as 323/380 (profile_from_effect /
emission_spread_vector / emission_rank_axis / rc_hybrid_block with
CallSet(400, 2200, 1.9) / block_stats / score); reads 323's reference and
cap15 effects read-only; writes only into experiments/385-*/work/.
Step 2+ pricing arms reuse emit_arm() with per-target CallSet sizes.

Usage: uv run python run385.py cap15
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

# Redirect numba's file cache into this experiment: a concurrent worker's 319o edit
# changed rc_emit.py's source tag, and an aliased module load from equiv385's gate
# poisoned the shared src/ cache lookup (ModuleNotFoundError at first wilcoxon_z
# compile). Recompiling the textually-unchanged jitted functions from the same
# numba version is bit-identical; src/ itself is never touched. See RESULTS.md.
os.environ.setdefault(
    "NUMBA_CACHE_DIR",
    str(Path(__file__).resolve().parent / "work" / "numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run323 as R  # noqa: E402 -- read-only import, never edited
import evaluate323 as E  # noqa: E402 -- read-only import, never edited
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import emit_counts as EMIT_COUNTS  # noqa: E402
import equiv385 as EQ  # noqa: E402 -- this experiment's drift gate (read-only pins)


def verify_with_known_drift() -> dict:
    """equiv385.verify385: all frozen inputs exact, except the two pinned files whose
    shas must equal the gate-time pins (proven inert: builder fns AST-identical, default
    rc_hybrid_block path bit-identical on a real target). Any further edit -> ValueError."""
    rec = EQ.verify385()
    eval_frozen = json.loads((BASE / "evaluation-inputs.json").read_text())
    from importlib.metadata import version as pkg_version
    for name, expected in eval_frozen["inputs_sha256"].items():
        if R.sha(ROOT / name) != expected:
            raise ValueError(f"evaluation input drift: {name}")
    if pkg_version("cell-eval2") != eval_frozen["cell_eval2_version"]:
        raise ValueError("scorer version drift")
    print("drift gate pins hold (see work/drift385.json)", flush=True)
    return rec


def builder385():
    """Direct builder load (R.builder() enforces the frozen gate, which is what drifted).
    Safety comes from drift gate A, not from this loader."""
    with R.model_env("uniform"):
        b = R.load_module("builder385", ROOT / "scripts/make_submission_pm_cd4.py")
    assert b.H1_UNION_G is None and b.MASS_FIX is None
    return b


def native_cap15(builder_module):
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        comp = json.loads((BASE / "work/compose.json").read_text())
        if builder_module.array_sha256(z["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
            raise ValueError("cap15 final effect hash drift vs compose.json")
        return z["cap15"].copy(), z["genes"].copy(), z["targets"].copy()


def emit_arm(name, effects_native, genes, callset_fn=None) -> dict:
    """Emit+score one arm. callset_fn(i, target) -> rc_emit.CallSet (default constant)."""
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError(f"refusing to overwrite metrics-{name}.json")
    b = builder385()
    r, c, emit_c = E.load_reference()
    effects = R.map_axis(effects_native, genes, r["genes"])
    pick_rng, emit_rng = b.make_rngs(0)
    sums, jks, ps, qs, lfcs, digests = [], [], [], [], [], []
    counts_sha = hashlib.sha256()
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        profile = b.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
        cs = callset_fn(i, str(target)) if callset_fn else rc_emit.CallSet(400, 2200, 1.9)
        block, _ = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], .1, emit_rng, spread=spread,
            overdispersion=4., control=emit_c, key=np.where(gate, effects[i], 0.),
            jitter_seed=1000 + i, call_set_first=True,
            call_set=cs, emit_counts=EMIT_COUNTS)
        s, j, p, q, l = E.block_stats(block, c)
        sums.append(s)
        jks.append(j)
        ps.append(p)
        qs.append(q)
        lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        if (i + 1) % 20 == 0:
            print(f"{name} {i + 1}/{len(r['targets'])} ({time.monotonic() - t0:.1f}s)", flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{name}.npz", targets=r["targets"], sums=np.asarray(sums),
             jk=np.asarray(jks), p=np.asarray(ps), q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{name}.csv")
    out = {"arm": name, "raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(),
           "mse_legs": legs, "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(WORK / f"stats-{name}.npz")}
    R.write_json(WORK / f"metrics-{name}.json", out)
    print(name, json.dumps({"raw": raw}), flush=True)
    return out


def known_answer(receipt: dict) -> dict:
    expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
    diffs = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
    checks = {
        "emitted_count_digest_equal": receipt["emitted_count_digest_sha256"]
        == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal": receipt["panel_total_digest_sha256"]
        == expected["panel_total_digest_sha256"],
        "raw_metrics_within_tolerance": {m: diffs[m] <= E.TOLERANCE[m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and all(checks["raw_metrics_within_tolerance"].values())
                  and checks["mse_legs_exact"])
    result = {"pass": passed, "checks": checks, "absolute_differences": diffs,
              "tolerance": E.TOLERANCE,
              "expected_source": "experiments/323-current-family-validation/work/metrics-cap15.json"}
    R.write_json(WORK / "known-answer.json", result)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", choices=["cap15"])
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    verify_with_known_drift()
    b = builder385()
    eff, genes, _ = native_cap15(b)
    receipt = emit_arm(args.arm, eff, genes)
    ka = known_answer(receipt)
    print("known-answer:", json.dumps(ka, indent=2))
    if not ka["pass"]:
        raise AssertionError("cap15 did not reproduce 323; see work/known-answer.json")


if __name__ == "__main__":
    main()
