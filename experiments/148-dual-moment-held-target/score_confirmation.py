#!/usr/bin/env python3
"""Score the fixed Experiment 148 candidate and evaluate the held 45 targets."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUTHORIZATION_PATH = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
GUARD_PATH = ROOT / "scripts/run_with_memory_guard.py"
REAL_PATH = ROOT / "data/harness/h1-test.h5ad"
CANDIDATE_PATH = HERE / "work/candidate.h5ad"
INCUMBENT_RESULTS = ROOT / (
    "experiments/095-clean-harness/out-h1test/runs/"
    "g0.65_q1.0_a1.08_spread0.25/run/results.csv"
)
INCUMBENT_AGG = ROOT / (
    "experiments/095-clean-harness/out-h1test/runs/"
    "g0.65_q1.0_a1.08_spread0.25/run/agg_results.csv"
)
OUTDIR = HERE / "work/scoring-v1"
PERT_COL = "target_gene"
CONTROL = "non-targeting"
GIB = 1024**3
MAX_RSS_GIB = 24.0
MIN_FREE_PERCENT = 15

ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}
MEMBERS = tuple(ANCHORS)

GUARD_SPEC = importlib.util.spec_from_file_location("memory_guard_148", GUARD_PATH)
assert GUARD_SPEC and GUARD_SPEC.loader
GUARD = importlib.util.module_from_spec(GUARD_SPEC)
sys.modules[GUARD_SPEC.name] = GUARD
GUARD_SPEC.loader.exec_module(GUARD)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def verify_authorization() -> dict[str, Any]:
    authorization = json.loads(AUTHORIZATION_PATH.read_text())
    if authorization.get("status") != "authorized_once":
        raise ValueError("development access is not authorized")
    scorer = authorization["scorer"]
    if scorer["sha256"] != sha256(Path(__file__).resolve()):
        raise ValueError("scorer source does not match authorization")
    bindings = authorization["bindings_sha256"]
    observed = {
        str((ROOT / relative).relative_to(ROOT)): sha256(ROOT / relative)
        for relative in bindings
    }
    if observed != bindings:
        raise ValueError("frozen input bindings do not match authorization")
    if authorization["held_targets"] != authorization["protocol_held_targets"]:
        raise ValueError("authorization and protocol held cohorts differ")
    return authorization


def process_group_stats(pgid: int) -> tuple[int, int, list[int]]:
    """Return RSS bytes, process count, and PIDs for one process group."""
    output = subprocess.check_output(
        ["ps", "-axo", "pid=,pgid=,rss="],
        text=True,
    )
    pids: list[int] = []
    total_kib = 0
    for line in output.splitlines():
        try:
            pid_text, pgid_text, rss_text = line.split()
            pid = int(pid_text)
            observed_pgid = int(pgid_text)
            rss_kib = int(rss_text)
        except ValueError:
            continue
        if observed_pgid == pgid:
            pids.append(pid)
            total_kib += rss_kib
    return total_kib * 1024, len(pids), sorted(pids)


def terminate_group_and_wait(pgid: int) -> list[int]:
    """Terminate a scorer process group and return any surviving PIDs."""
    for sig, timeout in ((signal.SIGTERM, 15.0), (signal.SIGKILL, 5.0)):
        try:
            os.killpg(pgid, sig)
        except ProcessLookupError:
            return []
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            _, count, pids = process_group_stats(pgid)
            if count == 0:
                return []
            time.sleep(0.25)
    return process_group_stats(pgid)[2]


def run_guarded(command: list[str], outdir: Path) -> dict[str, Any]:
    started = time.monotonic()
    peak = 0
    samples = 0
    minimum_free: int | None = None
    stdout_path = outdir / "cell-eval2.stdout.log"
    stderr_path = outdir / "cell-eval2.stderr.log"
    report_path = outdir / "memory-guard.json"
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        try:
            while True:
                root_returncode = process.poll()
                rss, count, pids = process_group_stats(process.pid)
                peak = max(peak, rss)
                samples += 1
                free = GUARD.free_percent()
                if free is not None:
                    minimum_free = (
                        free
                        if minimum_free is None
                        else min(minimum_free, free)
                    )
                if samples == 1 or samples % 10 == 0:
                    print(
                        f"[score-guard] rss={rss / GIB:.2f} GiB "
                        f"peak={peak / GIB:.2f} GiB procs={count} "
                        f"free={free if free is not None else '?'}%",
                        flush=True,
                    )
                reason: str | None = None
                if rss > MAX_RSS_GIB * GIB:
                    reason = (
                        f"process group exceeded {MAX_RSS_GIB:.1f} GiB"
                    )
                elif free is not None and free < MIN_FREE_PERCENT:
                    reason = (
                        "system free memory fell below "
                        f"{MIN_FREE_PERCENT}%"
                    )
                if reason is not None:
                    survivors = terminate_group_and_wait(process.pid)
                    report = {
                        "status": "limit_exceeded",
                        "reason": reason,
                        "peak_rss_gib": peak / GIB,
                        "minimum_free_percent": minimum_free,
                        "samples": samples,
                        "elapsed_seconds": time.monotonic() - started,
                        "last_observed_pids": pids,
                        "surviving_pids": survivors,
                        "process_group_extinct": not survivors,
                    }
                    atomic_write_json(report_path, report)
                    raise RuntimeError(reason)
                if root_returncode is not None and count == 0:
                    break
                time.sleep(1.0)
        except BaseException:
            _, count, _ = process_group_stats(process.pid)
            if count:
                terminate_group_and_wait(process.pid)
            raise
        returncode = process.wait()
    _, final_count, final_pids = process_group_stats(process.pid)
    report = {
        "status": "complete" if returncode == 0 else "command_failed",
        "returncode": returncode,
        "peak_rss_gib": peak / GIB,
        "minimum_free_percent": minimum_free,
        "samples": samples,
        "elapsed_seconds": time.monotonic() - started,
        "command": command,
        "final_process_count": final_count,
        "final_process_ids": final_pids,
        "process_group_extinct": final_count == 0,
    }
    atomic_write_json(report_path, report)
    if final_count:
        raise RuntimeError(
            f"scorer process group still has descendants: {final_pids}"
        )
    if returncode != 0:
        raise RuntimeError(
            f"cell-eval2 failed with {returncode}; see {stderr_path}"
        )
    return report


def long_metric(frame: pd.DataFrame, metric: str, targets: set[str]) -> pd.Series:
    rows = frame[
        (frame["metric"] == metric) & frame["perturbation"].isin(targets)
    ]
    return rows.set_index("perturbation")["value"].astype(float)


def aggregate_subset(frame: pd.DataFrame, targets: set[str]) -> dict[str, float]:
    result = {
        metric: float(long_metric(frame, metric, targets).mean())
        for metric in (
            "pds_cosine",
            "de_wilcoxon_lfc_nmae",
            "de_wilcoxon_direction_fidelity_yield_raw",
            "de_wilcoxon_direction_reach_raw",
            "de_wilcoxon_sig_jaccard",
        )
    }
    numerator = long_metric(frame, "expr_mse_unbiased_capped", targets)
    denominator = long_metric(frame, "expr_distance_unbiased", targets)
    result["expr_mse_unbiased_capped_norm"] = float(
        numerator.sum() / denominator.sum()
    )
    return result


def aggregate_full(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path).set_index("statistic")
    return {metric: float(frame.loc["mean", metric]) for metric in MEMBERS}


def official(raw: dict[str, float]) -> dict[str, float]:
    scores: dict[str, float] = {}
    for metric, value in raw.items():
        baseline, replicate = ANCHORS[metric]
        score = (value - baseline) / (replicate - baseline)
        if metric == "expr_mse_unbiased_capped_norm":
            score = max(0.0, min(1.0, score))
        scores[metric] = float(score)
    scores["avg_score"] = float(np.mean(list(scores.values())))
    return scores


def deltas(
    incumbent: dict[str, float],
    candidate: dict[str, float],
) -> dict[str, float]:
    return {metric: candidate[metric] - incumbent[metric] for metric in incumbent}


def main() -> int:
    if OUTDIR.exists():
        raise FileExistsError(f"immutable scoring output exists: {OUTDIR}")
    authorization = verify_authorization()
    OUTDIR.mkdir(parents=True)
    run_dir = OUTDIR / "run"
    command = [
        "uv",
        "run",
        "cell-eval2",
        "run",
        "-ap",
        str(CANDIDATE_PATH),
        "-ar",
        str(REAL_PATH),
        "-o",
        str(run_dir),
        "--preset",
        "vcc2026",
        "--pert-col",
        PERT_COL,
        "--control",
        CONTROL,
    ]
    memory = run_guarded(command, OUTDIR)

    candidate_results_path = run_dir / "results.csv"
    candidate_agg_path = run_dir / "agg_results.csv"
    incumbent_frame = pd.read_csv(INCUMBENT_RESULTS)
    candidate_frame = pd.read_csv(candidate_results_path)
    held = set(authorization["held_targets"])

    incumbent_held = aggregate_subset(incumbent_frame, held)
    candidate_held = aggregate_subset(candidate_frame, held)
    held_raw_delta = deltas(incumbent_held, candidate_held)
    incumbent_held_scaled = official(incumbent_held)
    candidate_held_scaled = official(candidate_held)
    held_scaled_delta = deltas(incumbent_held_scaled, candidate_held_scaled)

    incumbent_full = aggregate_full(INCUMBENT_AGG)
    candidate_full = aggregate_full(candidate_agg_path)
    full_raw_delta = deltas(incumbent_full, candidate_full)
    incumbent_full_scaled = official(incumbent_full)
    candidate_full_scaled = official(candidate_full)
    full_scaled_delta = deltas(incumbent_full_scaled, candidate_full_scaled)

    incumbent_pds = long_metric(incumbent_frame, "pds_cosine", held).sort_index()
    candidate_pds = long_metric(candidate_frame, "pds_cosine", held).sort_index()
    pds_exact = bool(
        incumbent_pds.index.equals(candidate_pds.index)
        and np.array_equal(incumbent_pds.to_numpy(), candidate_pds.to_numpy())
    )

    checks = {
        "all_held_pds_exact": pds_exact,
        "held_projected_avg_gain_at_least_0.004": (
            held_scaled_delta["avg_score"] >= 0.004
        ),
        "held_raw_nmae_delta_at_most_minus_0.004": (
            held_raw_delta["de_wilcoxon_lfc_nmae"] <= -0.004
        ),
        "held_raw_reach_delta_at_least_0.010": (
            held_raw_delta["de_wilcoxon_direction_reach_raw"] >= 0.010
        ),
        "held_raw_fidelity_delta_at_least_minus_0.003": (
            held_raw_delta["de_wilcoxon_direction_fidelity_yield_raw"] >= -0.003
        ),
        "held_raw_jaccard_delta_at_least_minus_0.010": (
            held_raw_delta["de_wilcoxon_sig_jaccard"] >= -0.010
        ),
        "held_raw_mse_delta_at_most_0.010": (
            held_raw_delta["expr_mse_unbiased_capped_norm"] <= 0.010
        ),
        "full_projected_avg_gain_at_least_0.005": (
            full_scaled_delta["avg_score"] >= 0.005
        ),
    }
    passed = all(checks.values())
    result = {
        "experiment": "148-dual-moment-held-target",
        "stage": "held_target_confirmation",
        "provenance_utc": "2026-08-31T18:45:00+00:00",
        "authorization_sha256": sha256(AUTHORIZATION_PATH),
        "candidate_results_sha256": sha256(candidate_results_path),
        "candidate_aggregate_sha256": sha256(candidate_agg_path),
        "held_target_count": len(held),
        "incumbent_held_raw": incumbent_held,
        "candidate_held_raw": candidate_held,
        "held_raw_delta": held_raw_delta,
        "incumbent_held_official_projection": incumbent_held_scaled,
        "candidate_held_official_projection": candidate_held_scaled,
        "held_scaled_delta": held_scaled_delta,
        "incumbent_full_raw": incumbent_full,
        "candidate_full_raw": candidate_full,
        "full_raw_delta": full_raw_delta,
        "incumbent_full_official_projection": incumbent_full_scaled,
        "candidate_full_official_projection": candidate_full_scaled,
        "full_scaled_delta": full_scaled_delta,
        "checks": checks,
        "passes": passed,
        "decision": "ADVANCE_TO_CROSS_CONTEXT" if passed else "TERMINAL_NO_GO",
        "memory_guard": memory,
        "parallel_scorers_used": False,
        "challenge_build_created": False,
        "challenge_submission_created": False,
    }
    atomic_write_json(OUTDIR / "results.json", result)
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "held_raw_delta": held_raw_delta,
                "held_scaled_delta": held_scaled_delta,
                "full_scaled_delta": full_scaled_delta,
                "checks": checks,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
