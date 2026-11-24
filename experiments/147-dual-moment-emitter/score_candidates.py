#!/usr/bin/env python3
"""Score Experiment 147's four frozen H1-fast artifacts serially."""

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

import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCORING_PROTOCOL_PATH = HERE / "DEVELOPMENT_SCORING_PROTOCOL.v1.json"
GUARD_PATH = ROOT / "scripts/run_with_memory_guard.py"
SPEC = importlib.util.spec_from_file_location("memory_guard", GUARD_PATH)
assert SPEC and SPEC.loader
GUARD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = GUARD
SPEC.loader.exec_module(GUARD)

REAL_PATH = ROOT / "data/harness/h1-fast.h5ad"
BUNDLE_PATH = ROOT / "experiments/002-calibration-sweep/fast/bundle"
ARTIFACT_ROOT = ROOT / "experiments/147-dual-moment-emitter/work/h1-fast-v2"
SCORING_ROOT = ARTIFACT_ROOT / "scoring-v1"
ARMS = (
    ("incumbent", ARTIFACT_ROOT / "incumbent.h5ad"),
    ("lambda-0.25", ARTIFACT_ROOT / "lambda-0.25.h5ad"),
    ("lambda-0.50", ARTIFACT_ROOT / "lambda-0.50.h5ad"),
    ("lambda-1.00", ARTIFACT_ROOT / "lambda-1.00.h5ad"),
)
PERT_COL = "target_gene"
CONTROL = "non-targeting"
GIB = 1024**3
MAX_RSS_GIB = 24.0
MIN_FREE_PERCENT = 15


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
    authorization = json.loads(SCORING_PROTOCOL_PATH.read_text())
    if authorization.get("status") != "authorized_once":
        raise ValueError("development scoring is not authorized")
    scorer = authorization.get("scorer", {})
    if scorer.get("path") != str(Path(__file__).resolve().relative_to(ROOT)):
        raise ValueError("authorization names a different scorer")
    if scorer.get("sha256") != sha256(Path(__file__).resolve()):
        raise ValueError("scorer source does not match authorization")
    if authorization.get("memory_guard_sha256") != sha256(GUARD_PATH):
        raise ValueError("memory guard does not match authorization")
    for section_name in ("development_protocol",):
        section = authorization.get(section_name, {})
        path = ROOT / str(section.get("path", ""))
        if not path.is_file() or section.get("sha256") != sha256(path):
            raise ValueError(f"{section_name} does not match authorization")
    build = authorization.get("build_evidence", {})
    for path_key, hash_key in (
        ("report", "report_sha256"),
        ("memory_report", "memory_report_sha256"),
    ):
        path = ROOT / str(build.get(path_key, ""))
        if not path.is_file() or build.get(hash_key) != sha256(path):
            raise ValueError(f"build evidence {path_key} does not match authorization")
    if build.get("development_perturbation_truth_opened") is not False:
        raise ValueError("build evidence does not preserve the truth-access boundary")
    if build.get("all_pseudobulk_exact") is not True:
        raise ValueError("build evidence does not certify exact pseudobulk")
    expected = authorization.get("candidate_sha256", {})
    observed = {
        str(path.relative_to(ROOT)): sha256(path)
        for _, path in ARMS
    }
    if observed != expected:
        raise ValueError("candidate artifacts do not match authorization")
    return authorization


def run_guarded(
    command: list[str],
    *,
    stdout_path: Path,
    stderr_path: Path,
    report_path: Path,
) -> None:
    started = time.monotonic()
    peak = 0
    samples = 0
    minimum_free: int | None = None
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        try:
            while process.poll() is None:
                rss, count = GUARD.process_rss(process.pid)
                peak = max(peak, rss)
                samples += 1
                free = GUARD.free_percent()
                if free is not None:
                    minimum_free = free if minimum_free is None else min(minimum_free, free)
                if samples == 1 or samples % 10 == 0:
                    print(
                        f"[score-guard] rss={rss / GIB:.2f} GiB "
                        f"peak={peak / GIB:.2f} GiB procs={count} "
                        f"free={free if free is not None else '?'}%",
                        flush=True,
                    )
                reason: str | None = None
                if rss > MAX_RSS_GIB * GIB:
                    reason = f"process tree exceeded {MAX_RSS_GIB:.1f} GiB"
                elif free is not None and free < MIN_FREE_PERCENT:
                    reason = f"system free memory fell below {MIN_FREE_PERCENT}%"
                if reason is not None:
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    atomic_write_json(
                        report_path,
                        {
                            "status": "limit_exceeded",
                            "reason": reason,
                            "peak_rss_gib": peak / GIB,
                            "minimum_free_percent": minimum_free,
                            "samples": samples,
                            "elapsed_seconds": time.monotonic() - started,
                        },
                    )
                    raise RuntimeError(reason)
                time.sleep(1.0)
        finally:
            returncode = process.wait()
    atomic_write_json(
        report_path,
        {
            "status": "complete" if returncode == 0 else "command_failed",
            "returncode": returncode,
            "peak_rss_gib": peak / GIB,
            "minimum_free_percent": minimum_free,
            "samples": samples,
            "elapsed_seconds": time.monotonic() - started,
            "command": command,
        },
    )
    if returncode != 0:
        raise RuntimeError(
            f"guarded command failed with {returncode}; see {stderr_path}"
        )


def read_metrics(agg_path: Path, score_path: Path) -> dict[str, Any]:
    raw_frame = pd.read_csv(agg_path).set_index("statistic")
    raw = {
        column: float(raw_frame.loc["mean", column])
        for column in raw_frame.columns
        if column != "statistic"
    }
    scaled_frame = pd.read_csv(score_path).set_index("metric")
    scaled = {
        str(metric): float(value)
        for metric, value in scaled_frame["from_replicate"].items()
    }
    return {"raw": raw, "scaled": scaled}


def main() -> int:
    if SCORING_ROOT.exists():
        raise FileExistsError(f"immutable scoring directory exists: {SCORING_ROOT}")
    authorization = verify_authorization()
    SCORING_ROOT.mkdir(parents=True)

    results: dict[str, dict[str, Any]] = {}
    for index, (label, prediction_path) in enumerate(ARMS, 1):
        print(f"[{index}/{len(ARMS)}] scoring {label}", flush=True)
        arm_dir = SCORING_ROOT / label
        run_dir = arm_dir / "run"
        arm_dir.mkdir()
        command = [
            "uv",
            "run",
            "cell-eval2",
            "run",
            "-ap",
            str(prediction_path),
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
        run_guarded(
            command,
            stdout_path=arm_dir / "cell-eval2.stdout.log",
            stderr_path=arm_dir / "cell-eval2.stderr.log",
            report_path=arm_dir / "memory-guard.json",
        )
        aggregate = run_dir / "agg_results.csv"
        score_path = arm_dir / "score.csv"
        score_process = subprocess.run(
            [
                "uv",
                "run",
                "cell-eval2",
                "score",
                "--user-agg",
                str(aggregate),
                "--real-bundle",
                str(BUNDLE_PATH),
                "-o",
                str(score_path),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        (arm_dir / "score.stdout.log").write_text(score_process.stdout)
        (arm_dir / "score.stderr.log").write_text(score_process.stderr)
        if score_process.returncode != 0:
            raise RuntimeError(
                f"cell-eval2 score failed for {label}: {score_process.stderr[-2000:]}"
            )
        results[label] = read_metrics(aggregate, score_path)
        atomic_write_json(SCORING_ROOT / "partial-results.json", results)
        print(
            f"[{index}/{len(ARMS)}] {label} "
            f"avg={results[label]['scaled']['avg_score']:.9f}",
            flush=True,
        )

    baseline = results["incumbent"]
    pds_name = "pds_cosine"
    mse_name = "expr_mse_unbiased_capped_norm"
    nmae_name = "de_wilcoxon_lfc_nmae"
    fid_name = "de_wilcoxon_direction_fidelity_yield_raw"
    reach_name = "de_wilcoxon_direction_reach_raw"
    comparisons: dict[str, dict[str, Any]] = {}
    for label, _ in ARMS[1:]:
        candidate = results[label]
        raw_delta = {
            metric: candidate["raw"][metric] - baseline["raw"][metric]
            for metric in baseline["raw"]
        }
        scaled_delta = {
            metric: candidate["scaled"][metric] - baseline["scaled"][metric]
            for metric in baseline["scaled"]
        }
        checks = {
            "raw_pds_exact": abs(raw_delta[pds_name]) <= 1e-12,
            "scaled_avg_gain_at_least_0.002": scaled_delta["avg_score"] >= 0.002,
            "raw_nmae_improves": raw_delta[nmae_name] < 0.0,
            "raw_fidelity_regression_at_most_0.001": raw_delta[fid_name] >= -0.001,
            "raw_reach_regression_at_most_0.001": raw_delta[reach_name] >= -0.001,
            "raw_mse_regression_at_most_0.01": raw_delta[mse_name] <= 0.01,
        }
        comparisons[label] = {
            "raw_delta_candidate_minus_incumbent": raw_delta,
            "scaled_delta_candidate_minus_incumbent": scaled_delta,
            "checks": checks,
            "passes": all(checks.values()),
        }

    passing = [
        label
        for label, comparison in comparisons.items()
        if comparison["passes"]
    ]
    selected = (
        max(
            passing,
            key=lambda label: comparisons[label][
                "scaled_delta_candidate_minus_incumbent"
            ]["avg_score"],
        )
        if passing
        else None
    )
    final = {
        "experiment": "147-dual-moment-emitter",
        "stage": "h1_fast_development_scoring",
        "provenance_utc": "2026-08-31T18:45:00+00:00",
        "authorization_sha256": sha256(SCORING_PROTOCOL_PATH),
        "authorization": authorization,
        "results": results,
        "comparisons": comparisons,
        "selected": selected,
        "decision": "ADVANCE" if selected else "TERMINAL_NO_GO",
        "parallel_scorers_used": False,
        "challenge_build_created": False,
        "challenge_submission_created": False,
    }
    atomic_write_json(SCORING_ROOT / "results.json", final)
    print(
        json.dumps(
            {
                "decision": final["decision"],
                "selected": selected,
                "comparisons": comparisons,
            },
            indent=2,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
