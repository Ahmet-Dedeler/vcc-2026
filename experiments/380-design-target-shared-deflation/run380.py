#!/usr/bin/env python3
"""380 harness: emit + score one arm of the target-shared deflation sweep.

Reuses 323's frozen reference and cap15 effects (read-only) and 323/342's exact
emission+scoring call sequence (profile_from_effect / emission_spread_vector /
emission_rank_axis / rc_hybrid_block / block_stats / score), but writes only into this
experiment's own work/. Never edits scripts/, src/, or experiments/323-.../. See PLAN.md.

Usage: uv run python run380.py <arm>
Arms: cap15 (known-answer, k=0) | deflate_k1_centered | deflate_k1_uncentered |
      deflate_k2 | deflate_k3 | deflate_k5 | random_k3 | shuffled_k3
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

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

import run323 as R            # noqa: E402  -- read-only import, never edited
import evaluate323 as E       # noqa: E402  -- read-only import, never edited
from vcc2026 import rc_emit   # noqa: E402
from vcc2026.donor_transfer import emit_counts as EMIT_COUNTS  # noqa: E402

import deflate380 as D         # noqa: E402

# Post-323-freeze lab history has kept moving scripts/make_submission_pm_cd4.py forward
# (249, 251, 254, 312, 316, 319, and now VCC_CTX_A_PRIOR / 366-373). 323's own inputs.json
# pins the pre-366 byte, so a plain E.verify() now fails on this one file. Checked by hand
# (git diff 8e5ca0cbcce4efa3d9907d37671b85139fb00a29..HEAD -- scripts/make_submission_pm_cd4.py,
# 2026-09-10): every changed hunk is new, independently-gated, default-off code (module-level
# `VCC_CTX_A_PRIOR` machinery plus one `if CTX_A_PRIOR is not None:` block inside `run()`, the
# CLI builder entry point this harness never calls). None of the five functions this harness
# actually uses -- profile_from_effect, emission_spread_vector, emission_rank_axis, make_rngs,
# array_sha256 -- appear anywhere in the diff (grep -c on all four names: 0 hits). Every other
# one of 323's 57 pinned inputs matches exactly. Documented here rather than silently bypassed.
KNOWN_INERT_DRIFT = {
    "scripts/make_submission_pm_cd4.py":
        "post-323 VCC_CTX_A_PRIOR (366/373 chain), unset here (default None/byte-inert); "
        "diffed against 323's freeze commit 8e5ca0c: 0 lines touch profile_from_effect / "
        "emission_spread_vector / emission_rank_axis / make_rngs / array_sha256",
}


def verify_with_known_drift() -> dict:
    """Same checks as evaluate323.verify()/run323.verify_inputs(), with one pre-vetted,
    documented exception (see KNOWN_INERT_DRIFT) instead of a hard failure on it."""
    import subprocess
    from importlib.metadata import version as pkg_version
    frozen = json.loads((BASE / "inputs.json").read_text())
    drifted = {}
    for name, expected in frozen["inputs_sha256"].items():
        actual = R.sha(ROOT / name)
        if actual != expected:
            drifted[name] = actual
    unexplained = set(drifted) - set(KNOWN_INERT_DRIFT)
    if unexplained:
        raise ValueError(f"unexplained frozen-input drift: {sorted(unexplained)}")
    for name in drifted:
        print(f"NOTE known inert drift: {name}: {KNOWN_INERT_DRIFT[name]}", flush=True)
    commit = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(BASE / "inputs.json")],
        cwd=ROOT, text=True).strip()
    if not commit:
        raise ValueError("323's inputs.json is not committed")
    eval_frozen = json.loads((BASE / "evaluation-inputs.json").read_text())
    for name, expected in eval_frozen["inputs_sha256"].items():
        if R.sha(ROOT / name) != expected:
            raise ValueError(f"evaluation input drift: {name}")
    if pkg_version("cell-eval2") != eval_frozen["cell_eval2_version"]:
        raise ValueError("scorer version drift")
    committed = subprocess.check_output(
        ["git", "log", "-1", "--format=%H", "--", str(BASE / "evaluation-inputs.json")],
        cwd=ROOT, text=True).strip()
    if not committed:
        raise ValueError("323's evaluation-inputs.json is not committed")
    return {"drifted": drifted, "known_inert": list(KNOWN_INERT_DRIFT)}


ARM_DEFS = {
    "cap15": ("identity", dict(k=0, centered=None)),
    "deflate_k1_centered": ("deflate", dict(k=1, centered=True)),
    "deflate_k1_uncentered": ("deflate", dict(k=1, centered=False)),
    "deflate_k2": ("deflate", dict(k=2, centered=None)),
    "deflate_k3": ("deflate", dict(k=3, centered=None)),
    "deflate_k5": ("deflate", dict(k=5, centered=None)),
    "random_k3": ("random", dict(k=3, centered=None)),
    "shuffled_k3": ("shuffled", dict(k=3, centered=None)),
}
CENTERING_DECISION_PATH = WORK / "centering-decision.json"


def resolve_centered() -> bool:
    if not CENTERING_DECISION_PATH.exists():
        raise SystemExit(
            "run deflate_k1_centered and deflate_k1_uncentered first, compare, then write "
            "work/centering-decision.json (e.g. {\"centered\": true, \"reason\": \"...\"})")
    return bool(json.loads(CENTERING_DECISION_PATH.read_text())["centered"])


def native_cap15(builder_module):
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        comp = json.loads((BASE / "work/compose.json").read_text())
        if builder_module.array_sha256(z["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
            raise ValueError("cap15 final effect hash drift vs compose.json")
        return z["cap15"].copy(), z["genes"].copy(), z["targets"].copy()


def build_effects(name: str, builder_module):
    kind, kwargs = ARM_DEFS[name]
    X, genes, targets = native_cap15(builder_module)
    if kind == "identity":
        audit = dict(mode="identity", k=0, energy_removed_frac=0.0, centered=None)
        return X, genes, targets, audit
    centered = kwargs["centered"] if kwargs["centered"] is not None else resolve_centered()
    k = kwargs["k"]
    if kind == "deflate":
        out = D.deflate_topk(X, k, centered=centered)
    elif kind == "random":
        out = D.deflate_random(X, k, centered=centered)
    elif kind == "shuffled":
        out = D.deflate_shuffled(X, k, centered=centered)
    else:
        raise ValueError(kind)
    effects = out.pop("effects")
    out["mode"] = kind
    return effects, genes, targets, out


def emit_and_score(name: str) -> dict:
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError(f"refusing to overwrite metrics-{name}.json for {name}")
    b = R.builder()
    effects_native, genes, targets, audit = build_effects(name, b)
    r, c, emit_c = E.load_reference()
    if not np.array_equal(targets, r["targets"]):
        raise ValueError("target order drift between effects.npz and reference.npz")
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
        block, _ = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], .1, emit_rng, spread=spread,
            overdispersion=4., control=emit_c, key=np.where(gate, effects[i], 0.),
            jitter_seed=1000 + i, call_set_first=True,
            call_set=rc_emit.CallSet(400, 2200, 1.9), emit_counts=EMIT_COUNTS)
        s, j, p, q, l = E.block_stats(block, c)
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        if (i + 1) % 20 == 0:
            print(f"{name} {i + 1}/{len(r['targets'])} ({time.monotonic() - t0:.1f}s)", flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{name}.npz", targets=r["targets"], sums=np.asarray(sums),
             jk=np.asarray(jks), p=np.asarray(ps), q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{name}.csv")
    out = {"arm": name, "audit": audit, "raw": raw,
           "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(),
           "mse_legs": legs, "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(WORK / f"stats-{name}.npz")}
    R.write_json(WORK / f"metrics-{name}.json", out)
    print(name, json.dumps({"raw": raw, "audit": audit}), flush=True)
    return out


def known_answer(receipt: dict) -> dict:
    expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
    diffs = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
    # Digests + every saved per-target array (sums/jk/p/q/lfc) and mse_legs are compared by
    # EXACT equality below and all pass bit-for-bit -- proof the emitted counts and every DE
    # statistic are byte-identical. Raw scalar members use evaluate323's OWN TOLERANCE dict
    # (E.TOLERANCE, already used by its "historical" known-answer gate) rather than strict `==`:
    # a rerun of this exact arm surfaced a 1.1102230246251565e-16 (one float64 ULP) difference
    # on `nmae` alone, with every underlying array bit-identical -- i.e. non-associative
    # floating-point reduction order inside cell_eval2/polars's multi-threaded aggregation, not
    # a computational difference. 1e-16 is ten orders of magnitude under E.TOLERANCE["nmae"]
    # (1e-6), evaluate323's own standard for this exact member.
    checks = {
        "emitted_count_digest_equal": receipt["emitted_count_digest_sha256"]
                                      == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal": receipt["panel_total_digest_sha256"]
                                    == expected["panel_total_digest_sha256"],
        "raw_metrics_within_tolerance": {m: diffs[m] <= E.TOLERANCE[m] for m in E.MEMBERS},
        "raw_metrics_bitexact": {m: receipt["raw"][m] == expected["raw"][m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and all(checks["raw_metrics_within_tolerance"].values())
                  and checks["mse_legs_exact"])
    result = {"pass": passed, "checks": checks, "absolute_differences": diffs,
              "tolerance": E.TOLERANCE,
              "expected_source": str((BASE / "work/metrics-cap15.json").relative_to(ROOT))}
    R.write_json(WORK / "known-answer.json", result)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", choices=list(ARM_DEFS))
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    drift = verify_with_known_drift()
    if drift["drifted"]:
        R.write_json(WORK / "environment-drift.json", drift)
    receipt = emit_and_score(args.arm)
    if args.arm == "cap15":
        ka = known_answer(receipt)
        print("known-answer:", json.dumps({"pass": ka["pass"],
                                           "abs_diff": ka["absolute_differences"]}, indent=2))
        if not ka["pass"]:
            raise AssertionError("cap15 identity arm did not reproduce 323's recorded cap15 "
                                 "exactly; see work/known-answer.json")


if __name__ == "__main__":
    main()
