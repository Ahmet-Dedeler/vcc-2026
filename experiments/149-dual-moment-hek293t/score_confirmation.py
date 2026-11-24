#!/usr/bin/env python3
"""Score the frozen held-HEK incumbent and dual-moment candidate serially."""

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
REAL_PATH = ROOT / "data/harness/orion-hek293t-thinned.h5ad"
BUILD_DIR = HERE / "work/source-build-v1"
EFFECTS_PATH = BUILD_DIR / "effects.npz"
PREDICTIONS = {
    "incumbent": BUILD_DIR / "incumbent.h5ad",
    "candidate": BUILD_DIR / "candidate.h5ad",
}
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

GUARD_SPEC = importlib.util.spec_from_file_location("memory_guard_149", GUARD_PATH)
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


def sha256_strings(values: np.ndarray) -> str:
    payload = "\n".join(map(str, values)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def verify_authorization() -> dict[str, Any]:
    authorization = json.loads(AUTHORIZATION_PATH.read_text())
    if authorization.get("status") != "authorized_once":
        raise ValueError("held-context scoring is not authorized")
    if authorization["scorer"]["sha256"] != sha256(Path(__file__).resolve()):
        raise ValueError("scorer source does not match authorization")
    bindings = authorization["bindings_sha256"]
    observed = {
        str((ROOT / relative).relative_to(ROOT)): sha256(ROOT / relative)
        for relative in bindings
    }
    if observed != bindings:
        raise ValueError("frozen scoring bindings do not match authorization")
    return authorization


def process_group_stats(pgid: int) -> tuple[int, int, list[int]]:
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


def run_guarded(command: list[str], arm_dir: Path) -> dict[str, Any]:
    started = time.monotonic()
    peak = 0
    samples = 0
    minimum_free: int | None = None
    stdout_path = arm_dir / "cell-eval2.stdout.log"
    stderr_path = arm_dir / "cell-eval2.stderr.log"
    report_path = arm_dir / "memory-guard.json"
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
                        f"[{arm_dir.name}] rss={rss / GIB:.2f} GiB "
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
                    atomic_write_json(
                        report_path,
                        {
                            "status": "limit_exceeded",
                            "reason": reason,
                            "peak_rss_gib": peak / GIB,
                            "minimum_free_percent": minimum_free,
                            "samples": samples,
                            "elapsed_seconds": time.monotonic() - started,
                            "last_observed_pids": pids,
                            "surviving_pids": survivors,
                            "process_group_extinct": not survivors,
                        },
                    )
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
    if returncode != 0:
        raise RuntimeError(
            f"cell-eval2 failed with {returncode}; see {stderr_path}"
        )
    if final_count:
        raise RuntimeError(
            f"scorer process group still has descendants: {final_pids}"
        )
    return report


def aggregate(path: Path) -> dict[str, float]:
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


def delta(left: dict[str, float], right: dict[str, float]) -> dict[str, float]:
    return {key: right[key] - left[key] for key in left}


def metric_arrays(
    frame: pd.DataFrame,
    targets: np.ndarray,
) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    target_pos = {str(target): index for index, target in enumerate(targets)}
    for metric in MEMBERS:
        values = np.full(len(targets), np.nan, dtype=np.float64)
        if metric == "expr_mse_unbiased_capped_norm":
            continue
        rows = frame[frame["metric"] == metric]
        for row in rows.itertuples(index=False):
            position = target_pos.get(str(row.perturbation))
            if position is not None:
                values[position] = float(row.value)
        result[metric] = values
    for metric in ("expr_mse_unbiased_capped", "expr_distance_unbiased"):
        values = np.full(len(targets), np.nan, dtype=np.float64)
        rows = frame[frame["metric"] == metric]
        for row in rows.itertuples(index=False):
            position = target_pos.get(str(row.perturbation))
            if position is not None:
                values[position] = float(row.value)
        result[metric] = values
    return result


def aggregate_indices(
    arrays: dict[str, np.ndarray],
    indices: np.ndarray,
) -> dict[str, float]:
    raw = {
        metric: float(np.nanmean(arrays[metric][indices]))
        for metric in MEMBERS
        if metric != "expr_mse_unbiased_capped_norm"
    }
    numerator = arrays["expr_mse_unbiased_capped"][indices]
    denominator = arrays["expr_distance_unbiased"][indices]
    keep = np.isfinite(numerator) & np.isfinite(denominator)
    raw["expr_mse_unbiased_capped_norm"] = float(
        numerator[keep].sum() / denominator[keep].sum()
    )
    return raw


def paired_bootstrap(
    incumbent_frame: pd.DataFrame,
    candidate_frame: pd.DataFrame,
    targets: np.ndarray,
    *,
    draws: int,
    seed: int,
) -> dict[str, float | int]:
    incumbent = metric_arrays(incumbent_frame, targets)
    candidate = metric_arrays(candidate_frame, targets)
    full_indices = np.arange(len(targets), dtype=np.int64)
    point = (
        official(aggregate_indices(candidate, full_indices))["avg_score"]
        - official(aggregate_indices(incumbent, full_indices))["avg_score"]
    )
    rng = np.random.default_rng(seed)
    samples = np.empty(draws, dtype=np.float64)
    for draw in range(draws):
        indices = rng.integers(0, len(targets), len(targets))
        candidate_score = official(
            aggregate_indices(candidate, indices)
        )["avg_score"]
        incumbent_score = official(
            aggregate_indices(incumbent, indices)
        )["avg_score"]
        samples[draw] = candidate_score - incumbent_score
    return {
        "draws": draws,
        "seed": seed,
        "estimate": float(point),
        "lower_2_5pct": float(np.percentile(samples, 2.5)),
        "upper_97_5pct": float(np.percentile(samples, 97.5)),
        "positive_draw_fraction": float(np.mean(samples > 0)),
    }


def main() -> int:
    if OUTDIR.exists():
        raise FileExistsError(f"immutable scoring output exists: {OUTDIR}")
    authorization = verify_authorization()
    OUTDIR.mkdir(parents=True)
    memory: dict[str, dict[str, Any]] = {}
    for label, prediction in PREDICTIONS.items():
        arm_dir = OUTDIR / label
        arm_dir.mkdir()
        print(f"[score] {label}", flush=True)
        memory[label] = run_guarded(
            [
                "uv",
                "run",
                "cell-eval2",
                "run",
                "-ap",
                str(prediction),
                "-ar",
                str(REAL_PATH),
                "-o",
                str(arm_dir / "run"),
                "--preset",
                "vcc2026",
                "--pert-col",
                PERT_COL,
                "--control",
                CONTROL,
            ],
            arm_dir,
        )

    incumbent_agg = aggregate(
        OUTDIR / "incumbent/run/agg_results.csv"
    )
    candidate_agg = aggregate(
        OUTDIR / "candidate/run/agg_results.csv"
    )
    raw_delta = delta(incumbent_agg, candidate_agg)
    incumbent_scaled = official(incumbent_agg)
    candidate_scaled = official(candidate_agg)
    scaled_delta = delta(incumbent_scaled, candidate_scaled)

    incumbent_frame = pd.read_csv(
        OUTDIR / "incumbent/run/results.csv"
    )
    candidate_frame = pd.read_csv(
        OUTDIR / "candidate/run/results.csv"
    )
    with np.load(EFFECTS_PATH, allow_pickle=False) as archive:
        targets = np.asarray(archive["targets"]).astype(str)
    if len(targets) != int(authorization["target_count"]):
        raise ValueError("authorized target count differs")
    if sha256_strings(targets) != authorization["target_order_sha256"]:
        raise ValueError("authorized target order differs")
    incumbent_pds = (
        incumbent_frame[incumbent_frame["metric"] == "pds_cosine"]
        .set_index("perturbation")
        .loc[targets, "value"]
        .to_numpy(dtype=np.float64)
    )
    candidate_pds = (
        candidate_frame[candidate_frame["metric"] == "pds_cosine"]
        .set_index("perturbation")
        .loc[targets, "value"]
        .to_numpy(dtype=np.float64)
    )
    pds_exact = bool(np.array_equal(incumbent_pds, candidate_pds))
    bootstrap = paired_bootstrap(
        incumbent_frame,
        candidate_frame,
        targets,
        draws=int(authorization["bootstrap"]["draws"]),
        seed=int(authorization["bootstrap"]["seed"]),
    )
    checks = {
        "all_per_target_pds_exact": pds_exact,
        "projected_avg_gain_at_least_0.004": (
            scaled_delta["avg_score"] >= 0.004
        ),
        "raw_nmae_delta_at_most_minus_0.004": (
            raw_delta["de_wilcoxon_lfc_nmae"] <= -0.004
        ),
        "raw_reach_delta_at_least_0.010": (
            raw_delta["de_wilcoxon_direction_reach_raw"] >= 0.010
        ),
        "raw_fidelity_delta_at_least_minus_0.003": (
            raw_delta["de_wilcoxon_direction_fidelity_yield_raw"] >= -0.003
        ),
        "raw_jaccard_delta_at_least_minus_0.010": (
            raw_delta["de_wilcoxon_sig_jaccard"] >= -0.010
        ),
        "raw_mse_delta_at_most_0.010": (
            raw_delta["expr_mse_unbiased_capped_norm"] <= 0.010
        ),
        "bootstrap_lower_bound_above_zero": (
            bootstrap["lower_2_5pct"] > 0.0
        ),
    }
    passed = all(checks.values())
    result = {
        "experiment": "149-dual-moment-hek293t",
        "stage": "held_context_confirmation",
        "provenance_utc": "2026-08-31T18:45:00+00:00",
        "authorization_sha256": sha256(AUTHORIZATION_PATH),
        "incumbent_results_sha256": sha256(
            OUTDIR / "incumbent/run/results.csv"
        ),
        "candidate_results_sha256": sha256(
            OUTDIR / "candidate/run/results.csv"
        ),
        "incumbent_aggregate_sha256": sha256(
            OUTDIR / "incumbent/run/agg_results.csv"
        ),
        "candidate_aggregate_sha256": sha256(
            OUTDIR / "candidate/run/agg_results.csv"
        ),
        "target_count": len(targets),
        "incumbent_raw": incumbent_agg,
        "candidate_raw": candidate_agg,
        "raw_delta": raw_delta,
        "incumbent_official_projection": incumbent_scaled,
        "candidate_official_projection": candidate_scaled,
        "scaled_delta": scaled_delta,
        "paired_target_bootstrap": bootstrap,
        "checks": checks,
        "passes": passed,
        "decision": (
            "ADVANCE_TO_PRODUCTION_BUILD" if passed else "TERMINAL_NO_GO"
        ),
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
                "raw_delta": raw_delta,
                "scaled_delta": scaled_delta,
                "bootstrap": bootstrap,
                "checks": checks,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
