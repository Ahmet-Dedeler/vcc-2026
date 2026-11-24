#!/usr/bin/env python3
"""342 — Hamilton integer-total rounding, remeasured on the current family.

Entry points: `freeze` | `verify` | `run` | `attest`. Nothing outside this directory is
written. `run` re-emits all 85 cap15 H1 blocks under three rounding rules through
production's `rc_hybrid_block` with ONE seamed statement (seam342.py), scores every
block with 323's exact evaluator, requires the floor arm to reproduce 323's recorded
cap15 outputs EXACTLY before any changed arm is interpreted, then writes the paired
report (report342.py). Freeze hash-binds code, registration, seam attestation, 323's
effects/reference/evaluator artifacts and the production sources; verify refuses a
dirty or uncommitted registered implementation.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

# One memory-guarded process with bounded threads, never parallel scorers.
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
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
import run323 as R            # noqa: E402
import evaluate323 as E       # noqa: E402
from vcc2026 import rc_emit   # noqa: E402
import rounding342 as RD      # noqa: E402
import seam342 as SM          # noqa: E402

ARMS = {"floor": "floor", "hamilton": "hamilton", "random_up": "random_up"}
ORDER = ("floor", "hamilton", "random_up")
OWN = ("PREREGISTRATION.md", "rounding342.py", "seam342.py", "run342.py",
       "report342.py", "audit342.py", "test342.py")
BUDGET_SECONDS = 2400
START_CUTOFF_UTC = datetime(2026, 9, 9, 19, 35, tzinfo=timezone.utc)
EMISSION = {
    "prediction_seed": 0, "emission_seed_offset": 1_000_003, "jitter_seed_base": 1000,
    "cells": 400, "min_cpm": 10.0, "tail_spread": 1.0, "keep_top": 8000,
    "keep_spread": 0.25, "spread_map": "block", "spread_lo": 0.0, "spread_hi": 0.60,
    "overdispersion": 4.0, "dispersion": 0.1, "call_set": [400, 2200, 1.9],
    "call_set_first": True, "h1_union": "OFF", "mass_fix": None,
    "effects": "323 cap15 final effects, identical for all three arms",
}
GATES = {
    "integrity": "known-answer exact, rounding invariants, identity checks, seam attestation, "
                 "no non-finite resamples",
    "point_gain": "hamilton local H1 overall point >= +0.001 versus floor",
    "lower_bound_vs_floor": "paired overall 95% lower bound > 0 versus floor",
    "lower_bound_vs_random_up": "paired overall 95% lower bound > 0 versus random_up",
    "no_resolved_mse_nmae_harm": "no raw MSE or NMAE contrast (hamilton versus floor or "
                                 "random_up) whose entire 95% interval is harmful",
    "no_resolved_fid_reach_jac_harm": "no raw fidelity/reach/Jaccard contrast (same two "
                                      "comparators) whose entire 95% interval is harmful",
    "pass_meaning": "LOCAL_GAIN_REQUIRES_SEPARATE_ABC_REVIEW only; never a cleared package or slot",
}
CELL_EVAL2_FILES = ("moments.py", "prep.py", "run.py", "de.py", "config.py", "catalog.py",
                    "configs/vcc2026.yaml", "metrics/delta.py", "metrics/discrimination.py",
                    "metrics/de.py", "metrics/direction.py")


def rel(path: Path) -> str:
    return str(Path(path).resolve().relative_to(ROOT))


def input_paths() -> list[Path]:
    files = [HERE / n for n in OWN]
    files += [BASE / n for n in ("inputs.json", "evaluation-inputs.json", "anchor-inputs.json",
                                 "run323.py", "evaluate323.py", "exact_de323.py",
                                 "report323.py")]
    files += [BASE / "work" / n for n in (
        "compose.json", "effects.npz", "reference.npz", "reference.json",
        "reference_de.parquet", "h1-anchors.json", "metrics-known.json",
        "metrics-cap15.json", "stats-cap15.npz", "per-target-cap15.csv",
        "whitener_h1_pruned.npz", "whitener_h1_pruned.json")]
    files += [ROOT / n for n in (
        "src/vcc2026/donor_transfer.py", "src/vcc2026/rc_emit.py",
        "src/vcc2026/whiten.py", "src/vcc2026/multi_donor.py",
        "scripts/make_submission_pm_cd4.py",
        "experiments/061-perturbmatch-hct-count/build_candidate.py",
        "experiments/062-cd4-post-whiten/run.py",
        "experiments/096-mse-floor/msekit.py")]
    return files


def cell_eval2_manifest() -> dict:
    import cell_eval2
    package = Path(cell_eval2.__file__).resolve().parent
    return {"version": version("cell-eval2"), "package_dir": str(package),
            "sha256": {f: R.sha(package / f) for f in CELL_EVAL2_FILES}}


def _plain(value):
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value)
    return value


def freeze() -> None:
    manifest_path = HERE / "inputs.json"
    if manifest_path.exists():
        raise FileExistsError("already frozen; do not overwrite inputs.json")
    if WORK.exists() and (list(WORK.glob("metrics-*.json")) or (WORK / "results.json").exists()):
        raise RuntimeError("generation already happened; freeze must precede it")
    seam = SM.Seam()
    hashes = {}
    for f in input_paths():
        hashes[rel(f)] = R.sha(f)
        print(f"pin {rel(f)}", flush=True)
    R.write_json(manifest_path, {
        "experiment": "342-hamilton-current-family",
        "classification": "reused H1 development remeasurement of an unsubmitted rounding "
                          "rider; not A/B/C confirmation, not a cleared package or slot",
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "inputs_sha256": hashes,
        "cell_eval2": cell_eval2_manifest(),
        "seam_attestation": seam.attestation,
        "arms": ARMS, "arm_order": list(ORDER),
        "baseline_effects": "experiments/323-current-family-validation/work/effects.npz::cap15",
        "control_seed_base": RD.CONTROL_SEED_BASE,
        "control_seed_rule": "np.random.default_rng(3420909 + target_index), target_index in "
                             "323's reference target order; the emission stream is never used",
        "emission": EMISSION,
        "bootstrap": {"draws": 10000, "seed": 20260909, "unit": "target", "paired": True,
                      "mse": "whole-panel ratio with the claim cap recomputed per resample"},
        "gates": GATES,
        "run_budget_seconds": BUDGET_SECONDS,
        "start_cutoff_utc": START_CUTOFF_UTC.isoformat(),
        "not_implemented_by_design": ["empirical-null call-set proposal", "block-sign repair"],
    })
    print(json.dumps({"seam": seam.attestation["line_differences"],
                      "pinned": len(hashes)}, indent=2), flush=True)


def verify(run_tests: bool = True) -> dict:
    manifest = json.loads((HERE / "inputs.json").read_text())
    expected = {rel(p) for p in input_paths()}
    if set(manifest["inputs_sha256"]) != expected:
        raise ValueError("input allowlist changed since freeze")
    for path, digest in manifest["inputs_sha256"].items():
        if R.sha(ROOT / path) != digest:
            raise ValueError(f"frozen input drift: {path}")
    own = [rel(HERE / n) for n in OWN] + [rel(HERE / "inputs.json")]
    for path in own:
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if tracked.returncode != 0:
            raise ValueError(f"registered file is not committed: {path}")
    dirty = subprocess.check_output(["git", "status", "--porcelain", "--", *own],
                                    cwd=ROOT, text=True)
    if dirty:
        raise ValueError(f"registered implementation is dirty; commit before measurement:\n{dirty}")
    if SM.Seam().attestation != manifest["seam_attestation"]:
        raise ValueError("seam attestation drift against the frozen manifest")
    if cell_eval2_manifest() != manifest["cell_eval2"]:
        raise ValueError("scorer package drift")
    E.verify()
    if run_tests:
        import test342
        test342.run_all()
    return manifest


def _rows_equal(a: Path, b: Path) -> bool:
    fa, fb = pl.read_csv(a), pl.read_csv(b)
    if fa.columns != fb.columns or fa.height != fb.height:
        return False
    def table(frame):
        return {(str(p), str(m)): float(v) for p, m, v in
                zip(frame["perturbation"], frame["metric"], frame["value"])}
    ka, kb = table(fa), table(fb)
    if ka.keys() != kb.keys():
        return False
    return all((np.isnan(ka[k]) and np.isnan(kb[k])) or ka[k] == kb[k] for k in ka)


def known_answer(receipt: dict) -> None:
    """The floor arm must reproduce 323's cap15 outputs EXACTLY; hashes are informational."""
    expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
    checks = {
        "emitted_count_digest_equal": receipt["emitted_count_digest_sha256"]
                                      == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal": receipt["panel_total_digest_sha256"]
                                    == expected["panel_total_digest_sha256"],
        "raw_metrics_exact": {m: receipt["raw"][m] == expected["raw"][m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    differences = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
    arrays = {}
    with np.load(BASE / "work/stats-cap15.npz", allow_pickle=False) as a, \
            np.load(WORK / "stats-floor.npz", allow_pickle=False) as b:
        for key in ("targets", "sums", "jk", "p", "q", "lfc"):
            x, y = a[key], b[key]
            same = x.shape == y.shape and x.dtype == y.dtype
            if same:
                same = bool(np.array_equal(x, y, equal_nan=(x.dtype.kind == "f")))
            arrays[key] = bool(same)
    checks["statistics_arrays_exact"] = arrays
    checks["per_gene_totals_exact"] = arrays["sums"]
    checks["per_target_rows_exact"] = _rows_equal(BASE / "work/per-target-cap15.csv",
                                                  WORK / "per-target-floor.csv")
    checks["statistics_file_sha256_equal_informational"] = (
        receipt["statistics_sha256"] == expected["statistics_sha256"])
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and all(checks["raw_metrics_exact"].values()) and checks["mse_legs_exact"]
                  and all(arrays.values()) and checks["per_target_rows_exact"])
    R.write_json(WORK / "known-answer.json", {
        "pass": passed, "checks": checks, "raw_absolute_differences": differences,
        "requirement": "exact equality with 323 cap15 on the emitted-count digest, the panel "
                       "digest, per-gene totals, every saved statistic, every per-target row "
                       "and all six raw metrics; the statistics file hash is informational",
        "expected_source": rel(BASE / "work/metrics-cap15.json")})
    if not passed:
        raise AssertionError("floor arm did not reproduce 323 cap15 exactly; "
                             "see work/known-answer.json (preserved, not tuned)")


def rounding_invariants(rule: str, records: list[dict]) -> dict:
    """Floor retains its truncation deficit; only changed arms restore rounded mass."""
    out = {
        "every_block_called_twice": all(t["calls"] == 2 for t in records),
        "every_block_mass_matches_rule": all(
            t["panel_total"] == t["floor_sum" if rule == "floor" else "round_sum_x"]
            for t in records),
        "every_block_zero_one_differences": all(t["zero_one_differences"] for t in records),
        "every_block_chosen_fractional": all(t["chosen_all_fractional"] for t in records),
        "every_block_moved_equals_expected": all(
            t["moved"] == (0 if rule == "floor" else t["residual_k"]) for t in records),
        "delivered_totals_equal_rule_totals": True,
        "scored_totals_equal_rule_totals": True,
    }
    out["all_pass"] = bool(records) and all(out.values())
    return out


def emit_arm(name: str, rule: str, r, c, emit_c, effects, builder, seam: SM.Seam,
             *, deadline: float) -> dict:
    pick_rng, emit_rng = builder.make_rngs(EMISSION["prediction_seed"])
    call_set = rc_emit.CallSet(*EMISSION["call_set"])
    sums, jks, ps, qs, lfcs, digests, rounding = [], [], [], [], [], [], []
    counts_sha = hashlib.sha256()
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        if time.monotonic() > deadline:
            raise TimeoutError(f"40-minute budget exceeded before {name}/{target}")
        profile = builder.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=EMISSION["min_cpm"])
        picks = pick_rng.choice(len(r["libraries"]), EMISSION["cells"], replace=False)
        spread = builder.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=EMISSION["min_cpm"],
            tail_spread=EMISSION["tail_spread"], keep_top=EMISSION["keep_top"],
            keep_top_frac=None, keep_spread=EMISSION["keep_spread"], rank_score=None,
            spread_map=EMISSION["spread_map"], spread_lo=EMISSION["spread_lo"],
            spread_hi=EMISSION["spread_hi"])
        gate = builder.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=EMISSION["min_cpm"])[0]
        rounder = RD.BlockRounder(rule, i)
        block, diagnostics = rc_emit.rc_hybrid_block(
            profile, r["libraries"][picks], EMISSION["dispersion"], emit_rng,
            spread=spread, overdispersion=EMISSION["overdispersion"], control=emit_c,
            key=np.where(gate, effects[i], 0.), jitter_seed=EMISSION["jitter_seed_base"] + i,
            call_set_first=EMISSION["call_set_first"], call_set=call_set,
            emit_counts=seam.emit_counts_with(rounder))
        if rounder.calls != 2:
            raise AssertionError(f"{name}/{target}: rounder called {rounder.calls} times, expected 2")
        if not np.array_equal(block.sum(0).astype(np.int64), rounder.totals):
            raise AssertionError(f"{name}/{target}: delivered totals differ from the rounding rule")
        s, j, p, q, l = E.block_stats(block, c)
        if not np.array_equal(s, rounder.totals):
            raise AssertionError(f"{name}/{target}: scored totals differ from the rounding rule")
        sums.append(s); jks.append(j); ps.append(p); qs.append(q); lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        summary = rounder.summary(gate=emit_c.keep)
        summary["target"] = str(target)
        summary["panel_total_digest_sha256"] = digests[-1]
        if rule == "random_up":
            ham = RD.hamilton_choice(rounder.x)
            union = np.union1d(ham, rounder.chosen).size
            summary["jaccard_with_hamilton_choice"] = (
                float(np.intersect1d(ham, rounder.chosen).size / union) if union else 1.0)
        summary["rc"] = {k: _plain(v) for k, v in diagnostics.items()}
        rounding.append(summary)
        del block
        if (i + 1) % 10 == 0:
            print(f"{name} {i + 1}/{len(r['targets'])} ({time.monotonic() - t0:.1f}s)", flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{name}.npz", targets=r["targets"], sums=np.asarray(sums),
             jk=np.asarray(jks), p=np.asarray(ps), q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{name}.csv")
    invariants = rounding_invariants(rule, rounding)
    keys = ("residual_k", "moved", "moved_in_scorer_gate", "moved_outside_scorer_gate",
            "n_fraction_positive", "panel_total", "max_abs_delta")
    summary = {k: {"min": int(min(t[k] for t in rounding)),
                   "median": float(np.median([t[k] for t in rounding])),
                   "max": int(max(t[k] for t in rounding))} for k in keys}
    delivered = [t["rc"]["call_set_delivered"] for t in rounding if "call_set_delivered" in t["rc"]]
    summary["call_set_delivered_median"] = float(np.median(delivered)) if delivered else None
    if rule == "random_up":
        summary["jaccard_with_hamilton_choice_median"] = float(np.median(
            [t["jaccard_with_hamilton_choice"] for t in rounding]))
    R.write_json(WORK / f"rounding-{name}.json", {
        "arm": name, "rule": rule, "invariants": invariants, "summary": summary,
        "targets": rounding})
    out = {"arm": name, "rule": rule, "raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(), "mse_legs": legs,
           "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(WORK / f"stats-{name}.npz"),
           "rounding_invariants": invariants, "rounding_summary": summary}
    R.write_json(WORK / f"metrics-{name}.json", out)
    print(name, json.dumps(raw), flush=True)
    return out


def run(decision: str | None = None) -> None:
    started = datetime.now(timezone.utc)
    if started > START_CUTOFF_UTC and not decision:
        raise SystemExit("starting after 19:35 UTC on 2026-09-09 requires a new coordinator "
                         "decision: pass --decision TEXT")
    manifest = verify()
    WORK.mkdir(parents=True, exist_ok=True)
    for name in ORDER:
        if (WORK / f"metrics-{name}.json").exists():
            raise FileExistsError(f"refusing to overwrite metrics-{name}.json")
    if (WORK / "results.json").exists():
        raise FileExistsError("results.json exists; preserved, not overwritten")
    timer = time.monotonic()
    record = {"mode": "run", "started_utc": started.isoformat(), "pid": os.getpid(),
              "coordinator_decision": decision, "status": "running",
              "input_manifest_sha256": R.sha(HERE / "inputs.json")}
    try:
        comp = json.loads((BASE / "work/compose.json").read_text())
        if not comp["gates_pass"]:
            raise ValueError("323 composition gates did not pass")
        if R.sha(BASE / "work/effects.npz") != comp["effects_sha256"]:
            raise ValueError("323 effects archive changed")
        builder = R.builder()
        if builder.H1_UNION_G is not None or builder.MASS_FIX is not None:
            raise AssertionError("H1 union or mass fix unexpectedly active")
        r, c, emit_c = E.load_reference()
        with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
            if not np.array_equal(z["targets"], r["targets"]):
                raise ValueError("target order drift between effects and reference")
            if builder.array_sha256(z["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
                raise ValueError("cap15 final effect hash drift")
            effects = R.map_axis(z["cap15"], z["genes"], r["genes"])
        seam = SM.Seam()
        print("seam:", json.dumps(seam.attestation["line_differences"]), flush=True)
        receipts = {}
        for name in ORDER:
            elapsed = time.monotonic() - timer
            if elapsed > BUDGET_SECONDS:
                raise TimeoutError(f"40-minute budget exceeded before {name} ({elapsed:.0f}s)")
            print(f"== arm {name} / rule {ARMS[name]} ({elapsed:.0f}s elapsed)", flush=True)
            receipts[name] = emit_arm(name, ARMS[name], r, c, emit_c, effects, builder, seam,
                                      deadline=timer+BUDGET_SECONDS)
            if name == "floor":
                known_answer(receipts[name])
                print("known-answer: floor reproduces 323 cap15 exactly", flush=True)
            gc.collect()
        del r, c, emit_c, effects
        gc.collect()
        if time.monotonic() - timer > BUDGET_SECONDS:
            raise TimeoutError("40-minute budget exceeded before the report")
        import report342
        result = report342.report()
        record.update(status="complete", result_status=result["status"],
                      gates=result["gates"])
    except BaseException as exc:
        record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        record.update(finished_utc=datetime.now(timezone.utc).isoformat(),
                      elapsed_seconds=time.monotonic() - timer,
                      process_max_rss_bytes_macos=resource.getrusage(
                          resource.RUSAGE_SELF).ru_maxrss,
                      frozen_utc=manifest["frozen_utc"])
        R.write_json(WORK / "execution.json", record)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("freeze", "verify", "run", "attest"))
    parser.add_argument("--decision", default=None,
                        help="coordinator decision text; required to start after the cutoff")
    args = parser.parse_args()
    if args.mode == "attest":
        print(json.dumps(SM.Seam().attestation, indent=2, sort_keys=True))
    elif args.mode == "run":
        run(args.decision)
    elif args.mode == "verify":
        verify()
        print("verify: frozen inputs, committed implementation, seam, scorer and tests OK")
    else:
        freeze()


if __name__ == "__main__":
    main()
