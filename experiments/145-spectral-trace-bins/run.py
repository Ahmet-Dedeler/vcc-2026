#!/usr/bin/env python3
"""Fail-closed runner for Experiment 145.

``prepare`` creates the frozen spectral candidates and null plans without
loading development truth. ``develop`` performs one authorized 613-target
nested evaluation. ``audit`` never loads truth.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import Any, Iterator, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = Path(__file__).resolve()
TEST_PATH = HERE / "test_run.py"
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
PREPARE_AUTHORIZATION = HERE / "SOURCE_ONLY_PREPARE_AUTHORIZATION.v2.json"
DEVELOPMENT_AUTHORIZATION = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v2.json"
MEMORY_GUARD = ROOT / "scripts/run_with_memory_guard.py"

EXP140 = ROOT / "experiments/140-batch-matched-controls"
EXP140_IDENTITY = EXP140 / "work/identity-reproduction-v5"
EXP140_IDENTITY_MANIFEST = EXP140_IDENTITY / "artifact-manifest.json"

EXP142 = ROOT / "experiments/142-batch-reliability-weighting"
EXP142_DEVELOP = EXP142 / "develop.py"
EXP142_DEVELOPMENT_AUTHORIZATION = EXP142 / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
EXP142_PREOUTCOME = EXP142 / "work/development-preoutcome-v1"
EXP142_PREOUTCOME_MANIFEST = EXP142_PREOUTCOME / "artifact-manifest.json"
EXP142_TRUTH_LEDGER = EXP142_PREOUTCOME / "truth-source-hashes.json"
EXP142_PLAN = EXP142_PREOUTCOME / "development-plan.npz"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
EXP132_RUN = EXP132 / "run.py"
EXP132_AUTHORIZATION = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
EXP132_PLAN = EXP132 / "work/preregistered-freeze-v2"
EXP132_PLAN_MANIFEST = EXP132_PLAN / "plan-manifest.json"
DEVELOPMENT_TARGETS = EXP132_PLAN / "selection_targets.txt"
RESPONSE_GENES = EXP132_PLAN / "response_genes.txt"

GROUPS = ("k562", "non_k562")
SCORED_DATASETS = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat")
FAMILIES = ("k562", "rpe1", "jurkat")

PREOUTCOME_DIR = HERE / "work/preoutcome-v2"
RESULT_DIR = HERE / "work/development-result-v2"
PREPARE_LOCK = HERE / "work/preoutcome-v2.lock"
DEVELOP_LOCK = HERE / "work/development-result-v2.lock"
OUTCOME_ACCESS_LOCK = HERE / "work/development-outcome-access.v2.lock"

PLAN_ARTIFACT = "development-plan.npz"
NULL_PLAN_ARTIFACT = "null-plans.npz"
TRUTH_LEDGER_ARTIFACT = "truth-source-hashes.json"
RESULT_ARTIFACT = "result.json"
ARRAY_ARTIFACT = "per-target-bootstrap.npz"
MANIFEST_ARTIFACT = "artifact-manifest.json"

DEV_COUNT = 613
RESPONSE_COUNT = 4992
FOLD_COUNT = 5
EXPECTED_FOLD_SIZES = (123, 123, 123, 122, 122)
BOOTSTRAP_DRAWS = 20_000
SPECTRAL_NULL_DRAWS = 99
TARGET_NULL_DRAWS = 99
SPECTRAL_NULL_SEED = 145_202_601
TARGET_NULL_SEED = 145_202_602
BETAS = (0.0, 0.25, 0.50, 0.75, 1.00)
TRACE_BOUNDARIES = {
    "k562": (93, 806, 2006, 4992),
    "non_k562": (198, 994, 2252, 4992),
}
MIN_GAIN = 0.003
MAX_CONCENTRATION = 0.5
EPS = 1e-12
MAX_RSS_GIB = 8.0
MIN_FREE_PERCENT = 25


class ContractError(RuntimeError):
    """Raised when a frozen experimental boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    require(bool(values), f"empty axis: {path}")
    require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while piece := handle.read(block):
            digest.update(piece)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    return {"bytes": path.stat().st_size, "sha256": sha256_file(path)}


def directory_manifest(directory: Path) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "files": {
            path.name: file_record(path)
            for path in sorted(directory.iterdir())
            if path.is_file() and path.name != MANIFEST_ARTIFACT
        },
    }


def validate_directory_manifest(
    directory: Path, manifest_name: str = MANIFEST_ARTIFACT
) -> dict[str, Any]:
    manifest_path = directory / manifest_name
    require(manifest_path.is_file(), f"missing artifact manifest: {directory}")
    manifest = read_json(manifest_path)
    expected_names = set(manifest["files"])
    actual_names = {
        path.name
        for path in directory.iterdir()
        if path.is_file() and path.name != manifest_name
    }
    require(actual_names == expected_names, f"artifact set changed: {directory}")
    for name, expected in manifest["files"].items():
        path = directory / name
        require(path.stat().st_size == expected["bytes"], f"bytes changed: {path}")
        require(sha256_file(path) == expected["sha256"], f"hash changed: {path}")
    return manifest


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".npz"
    )
    os.close(descriptor)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_copy(source: Path, destination: Path) -> None:
    require(not destination.exists(), f"refuse overwrite: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    try:
        shutil.copyfile(source, temporary)
        with Path(temporary).open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


@contextmanager
def atomic_stage_directory(path: Path) -> Iterator[Path]:
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{path.name}.", dir=path.parent))
    try:
        yield staging
        os.replace(staging, path)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


@contextmanager
def exclusive_lock(path: Path) -> Iterator[None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor: int | None = None
    owner: tuple[int, int] | None = None
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        stat = os.fstat(descriptor)
        owner = (stat.st_dev, stat.st_ino)
        os.write(descriptor, f"{os.getpid()}\n".encode())
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None
        yield
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if owner is not None:
            try:
                current = path.stat()
            except FileNotFoundError:
                pass
            else:
                if (current.st_dev, current.st_ino) == owner:
                    path.unlink()


@contextmanager
def process_tree_memory_guard() -> Iterator[dict[str, Any]]:
    require(MEMORY_GUARD.is_file(), f"missing memory guard: {MEMORY_GUARD}")
    (HERE / "work").mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".development-memory-v2.", dir=HERE / "work"
    ) as temporary:
        temp = Path(temporary)
        ready = temp / "ready.json"
        stop = temp / "stop"
        report = temp / "report.json"
        command = [
            sys.executable,
            str(MEMORY_GUARD),
            "--max-rss-gib",
            str(MAX_RSS_GIB),
            "--min-free-percent",
            str(MIN_FREE_PERCENT),
            "--poll-seconds",
            "0.25",
            "--monitor-pid",
            str(os.getpid()),
            "--ready-file",
            str(ready),
            "--stop-file",
            str(stop),
            "--report-file",
            str(report),
            "--quiet",
        ]
        monitor = subprocess.Popen(command)
        deadline = time.monotonic() + 30.0
        while not ready.exists():
            require(monitor.poll() is None, "memory monitor exited before ready")
            require(time.monotonic() < deadline, "memory monitor readiness timed out")
            time.sleep(0.05)
        observed: dict[str, Any] = {}
        try:
            yield observed
        finally:
            if monitor.poll() is None:
                stop.touch(exist_ok=False)
            try:
                returncode = monitor.wait(timeout=30.0)
            except subprocess.TimeoutExpired:
                monitor.kill()
                monitor.wait(timeout=10.0)
                raise ContractError("memory monitor did not stop")
            require(returncode == 0, f"memory monitor failed: {returncode}")
            payload = read_json(report)
            require(payload.get("status") == "complete", "memory monitor incomplete")
            require(
                float(payload["peak_rss_gib"]) <= MAX_RSS_GIB,
                "process-tree RSS exceeded limit",
            )
            minimum_free = payload.get("minimum_free_percent")
            require(
                minimum_free is None or int(minimum_free) >= MIN_FREE_PERCENT,
                "system free memory fell below limit",
            )
            observed.update(payload)


def protocol() -> dict[str, Any]:
    value = read_json(PROTOCOL_PATH)
    require(value.get("schema_version") == 1, "protocol schema changed")
    require(value.get("experiment") == HERE.name, "protocol experiment changed")
    require(tuple(value["candidate"]["betas"]) == BETAS, "beta grid changed")
    for group in GROUPS:
        require(
            tuple(value["candidate"]["descending_trace_boundaries_one_based"][group])
            == TRACE_BOUNDARIES[group],
            f"{group}: trace boundaries changed",
        )
    return value


def resolve_repo_path(value: str) -> Path:
    path = (ROOT / value).resolve()
    require(ROOT == path or ROOT in path.parents, f"path escapes repository: {value}")
    return path


def authorization_bindings(authorization: Mapping[str, Any]) -> dict[Path, str]:
    raw = authorization.get("bindings")
    require(isinstance(raw, Mapping), "authorization lacks bindings")
    result = {}
    for name, value in raw.items():
        path = Path(str(name))
        resolved = (
            path.resolve() if path.is_absolute() else resolve_repo_path(str(name))
        )
        digest = value if isinstance(value, str) else value.get("sha256")
        require(isinstance(digest, str), f"invalid authorization binding: {name}")
        require(resolved not in result, f"duplicate authorization binding: {name}")
        result[resolved] = digest
    return result


def validate_bound_paths(
    authorization: Mapping[str, Any],
    required: Sequence[Path],
    *,
    hash_files: bool,
) -> dict[Path, str]:
    bindings = authorization_bindings(authorization)
    for path in required:
        resolved = path.resolve()
        require(resolved in bindings, f"authorization missing binding: {path}")
        require(resolved.is_file(), f"authorized file missing: {path}")
        if hash_files:
            require(
                sha256_file(resolved) == bindings[resolved],
                f"authorization hash mismatch: {path}",
            )
    return bindings


def require_safety_flags(authorization: Mapping[str, Any]) -> None:
    for key in (
        "untouched_410_test_authorized",
        "reserved_confirmation_expression_authorized",
        "count_generation_authorized",
        "cell_eval2_authorized",
        "challenge_build_authorized",
        "submission_authorized",
        "null_stage_authorized",
    ):
        require(authorization.get(key) is False, f"unsafe authorization flag: {key}")


def source_paths() -> tuple[Path, ...]:
    return (
        RUN_PATH,
        TEST_PATH,
        PROTOCOL_PATH,
        EXP140_IDENTITY_MANIFEST,
        EXP140_IDENTITY / "whitener-k562.npz",
        EXP140_IDENTITY / "whitener-non_k562.npz",
        EXP142_DEVELOP,
        EXP142_PREOUTCOME_MANIFEST,
        EXP142_TRUTH_LEDGER,
        EXP142_PLAN,
        EXP142_PREOUTCOME / "geometry-k562.npz",
        EXP142_PREOUTCOME / "geometry-non_k562.npz",
        EXP132_PLAN_MANIFEST,
        DEVELOPMENT_TARGETS,
        RESPONSE_GENES,
    )


def validate_prepare_authorization() -> dict[str, Any]:
    require(
        PREPARE_AUTHORIZATION.is_file(),
        "SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json is required",
    )
    authorization = read_json(PREPARE_AUTHORIZATION)
    require(
        authorization.get("source_only_prepare_authorized") is True,
        "source-only prepare is not authorized",
    )
    require(
        authorization.get("development_613_authorized") is False,
        "prepare authorization must not open development truth",
    )
    require_safety_flags(authorization)
    validate_bound_paths(authorization, source_paths(), hash_files=True)
    return authorization


def load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_exp142() -> ModuleType:
    return load_module(EXP142_DEVELOP, "experiment142_for_145")


def trace_boundaries(evals: np.ndarray) -> tuple[int, int, int, int]:
    values = np.asarray(evals, dtype=np.float64)
    require(values.shape == (RESPONSE_COUNT,), "eigenvalue shape changed")
    require(np.all(np.isfinite(values)) and np.all(values > 0.0), "bad eigenvalues")
    require(np.all(np.diff(values) >= 0.0), "eigenvalues are not nondecreasing")
    descending = values[::-1]
    cumulative = np.cumsum(descending) / descending.sum()
    boundaries = tuple(
        int(np.searchsorted(cumulative, fraction, side="left") + 1)
        for fraction in (0.25, 0.50, 0.75)
    )
    return (*boundaries, RESPONSE_COUNT)


def band_labels_ascending(boundaries: Sequence[int]) -> np.ndarray:
    b1, b2, b3, end = map(int, boundaries)
    require(0 < b1 < b2 < b3 < end == RESPONSE_COUNT, "invalid boundaries")
    descending = np.empty(RESPONSE_COUNT, dtype=np.uint8)
    descending[:b1] = 0
    descending[b1:b2] = 1
    descending[b2:b3] = 2
    descending[b3:end] = 3
    return descending[::-1].copy()


def multipliers(labels: np.ndarray, beta: float) -> np.ndarray:
    lab = np.asarray(labels, dtype=np.uint8)
    require(lab.shape == (RESPONSE_COUNT,), "band label shape changed")
    values = np.asarray(
        [1.0, 1.0, np.exp(-beta / 2.0), np.exp(-beta)],
        dtype=np.float64,
    )
    return values[lab]


def row_norm_preserving_spectral_transform(
    baseline: np.ndarray,
    evecs: np.ndarray,
    labels: np.ndarray,
    beta: float,
    *,
    chunk: int = 64,
) -> np.ndarray:
    if beta == 0.0:
        return np.asarray(baseline).copy()
    base = np.asarray(baseline, dtype=np.float64)
    vectors = np.asarray(evecs, dtype=np.float32)
    require(vectors.shape == (RESPONSE_COUNT, RESPONSE_COUNT), "eigenvector shape")
    weight = multipliers(labels, beta)
    output = np.empty(base.shape, dtype=np.float32)
    for lo in range(0, len(base), chunk):
        hi = min(lo + chunk, len(base))
        block = base[lo:hi]
        transformed = ((block @ vectors) * weight[None, :]) @ vectors.T
        old_norm = np.linalg.norm(block, axis=1, keepdims=True)
        new_norm = np.linalg.norm(transformed, axis=1, keepdims=True)
        transformed = old_norm * np.divide(
            transformed,
            np.maximum(new_norm, EPS),
            out=np.zeros_like(transformed),
        )
        output[lo:hi] = transformed.astype(np.float32)
    return output


def row_angles_degrees(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    left = np.asarray(a, dtype=np.float64)
    right = np.asarray(b, dtype=np.float64)
    numerator = np.sum(left * right, axis=1)
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    cosine = np.divide(
        numerator,
        np.maximum(denominator, EPS),
        out=np.ones_like(numerator),
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def stable_derangement(rng: np.random.Generator, values: np.ndarray) -> np.ndarray:
    source = np.asarray(values, dtype=np.int64)
    require(len(source) >= 2, "cannot derange fewer than two rows")
    for _ in range(10_000):
        candidate = source[rng.permutation(len(source))]
        if np.all(candidate != source):
            return candidate
    raise ContractError("could not construct derangement")


def make_null_plans(folds: np.ndarray) -> dict[str, np.ndarray]:
    spectral_rng = np.random.default_rng(SPECTRAL_NULL_SEED)
    target_rng = np.random.default_rng(TARGET_NULL_SEED)
    plans = {}
    for group in GROUPS:
        labels = band_labels_ascending(TRACE_BOUNDARIES[group])
        permutations = np.empty((SPECTRAL_NULL_DRAWS, RESPONSE_COUNT), dtype=np.uint8)
        for draw in range(SPECTRAL_NULL_DRAWS):
            permutations[draw] = spectral_rng.permutation(labels)
        require(
            all(
                np.array_equal(
                    np.bincount(permutations[draw], minlength=4),
                    np.bincount(labels, minlength=4),
                )
                for draw in range(SPECTRAL_NULL_DRAWS)
            ),
            f"{group}: spectral null changed bin sizes",
        )
        plans[f"spectral_band_labels__{group}"] = permutations

    target_permutations = np.empty((TARGET_NULL_DRAWS, DEV_COUNT), dtype=np.int16)
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    for draw in range(TARGET_NULL_DRAWS):
        permutation = np.empty(DEV_COUNT, dtype=np.int64)
        for fold in range(FOLD_COUNT):
            indices = all_indices[folds == fold]
            permutation[indices] = stable_derangement(target_rng, indices)
        require(np.all(permutation != all_indices), "target null has fixed points")
        require(
            np.array_equal(folds[permutation], folds),
            "target null crossed outer folds",
        )
        target_permutations[draw] = permutation.astype(np.int16)
    plans["target_alignment_permutations"] = target_permutations
    return plans


def validate_plan_arrays(
    folds: np.ndarray, blocks: np.ndarray, draws: np.ndarray
) -> None:
    fold_values = np.asarray(folds)
    block_values = np.asarray(blocks)
    draw_values = np.asarray(draws)
    require(fold_values.shape == (DEV_COUNT,), "fold vector shape changed")
    require(block_values.shape == (DEV_COUNT,), "block vector shape changed")
    require(
        np.array_equal(np.unique(fold_values), np.arange(FOLD_COUNT)),
        "fold labels are not exactly 0..4",
    )
    require(
        tuple(int(np.sum(fold_values == fold)) for fold in range(FOLD_COUNT))
        == EXPECTED_FOLD_SIZES,
        "fold sizes changed",
    )
    block_count = int(block_values.max()) + 1
    require(
        np.array_equal(np.unique(block_values), np.arange(block_count)),
        "block labels are not contiguous",
    )
    require(block_count == 531, "block count changed")
    for block in range(block_count):
        require(
            len(np.unique(fold_values[block_values == block])) == 1,
            f"block {block} crosses outer folds",
        )
    require(
        draw_values.shape == (BOOTSTRAP_DRAWS, block_count),
        "bootstrap plan changed",
    )
    require(
        np.all((draw_values >= 0) & (draw_values < block_count)),
        "bootstrap contains an out-of-range block",
    )


def prepare() -> None:
    validate_prepare_authorization()
    protocol()
    validate_directory_manifest(EXP140_IDENTITY)
    validate_directory_manifest(EXP142_PREOUTCOME)
    require(not PREOUTCOME_DIR.exists(), f"preoutcome exists: {PREOUTCOME_DIR}")
    require(not RESULT_DIR.exists(), "development result already exists")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")

    with np.load(EXP142_PLAN, allow_pickle=False) as data:
        folds = np.asarray(data["development_folds"], dtype=np.int8)
        blocks = np.asarray(data["development_blocks"], dtype=np.int16)
        draws = np.asarray(data["bootstrap_block_draws"], dtype=np.uint16)
    validate_plan_arrays(folds, blocks, draws)
    null_plans = make_null_plans(folds)

    with (
        exclusive_lock(PREPARE_LOCK),
        atomic_stage_directory(PREOUTCOME_DIR) as staging,
    ):
        atomic_npz(
            staging / PLAN_ARTIFACT,
            development_folds=folds,
            development_blocks=blocks,
            bootstrap_block_draws=draws,
        )
        atomic_npz(staging / NULL_PLAN_ARTIFACT, **null_plans)
        reports = {}
        for group in GROUPS:
            with np.load(
                EXP140_IDENTITY / f"whitener-{group}.npz", allow_pickle=False
            ) as operator:
                evals = np.asarray(operator["evals"], dtype=np.float64)
                evecs = np.asarray(operator["evecs"], dtype=np.float32)
            observed_boundaries = trace_boundaries(evals)
            require(
                observed_boundaries == TRACE_BOUNDARIES[group],
                f"{group}: trace boundary drift {observed_boundaries}",
            )
            gram = evecs.T @ evecs
            orthogonality_error = float(
                np.max(np.abs(gram - np.eye(RESPONSE_COUNT, dtype=np.float32)))
            )
            require(orthogonality_error <= 2e-5, f"{group}: eigenvectors drift")
            with np.load(
                EXP142_PREOUTCOME / f"geometry-{group}.npz", allow_pickle=False
            ) as geometry:
                baseline = np.asarray(geometry["baseline"], dtype=np.float32)
            require(
                baseline.shape == (DEV_COUNT, RESPONSE_COUNT),
                f"{group}: baseline shape",
            )
            labels = band_labels_ascending(observed_boundaries)
            predictions = np.empty(
                (len(BETAS), DEV_COUNT, RESPONSE_COUNT), dtype=np.float32
            )
            diagnostics = {}
            for index, beta in enumerate(BETAS):
                predictions[index] = row_norm_preserving_spectral_transform(
                    baseline, evecs, labels, beta
                )
                if beta == 0.0:
                    require(
                        np.array_equal(predictions[index], baseline),
                        f"{group}: beta-zero identity failed",
                    )
                old_norm = np.linalg.norm(baseline, axis=1)
                new_norm = np.linalg.norm(predictions[index], axis=1)
                relative = np.abs(new_norm - old_norm) / np.maximum(old_norm, EPS)
                angles = row_angles_degrees(baseline, predictions[index])
                diagnostics[f"{beta:.2f}"] = {
                    "median_angle_degrees": float(np.median(angles)),
                    "p95_angle_degrees": float(np.quantile(angles, 0.95)),
                    "max_angle_degrees": float(angles.max()),
                    "max_relative_row_norm_error": float(relative.max()),
                    "nonfinite_values": int(
                        np.size(predictions[index])
                        - np.isfinite(predictions[index]).sum()
                    ),
                    "zero_rows": int(np.sum(new_norm <= EPS)),
                }
                require(np.all(np.isfinite(predictions[index])), "nonfinite candidate")
                require(float(relative.max()) <= 2e-6, "row norm drift")
            atomic_npz(
                staging / f"candidates-{group}.npz",
                betas=np.asarray(BETAS, dtype=np.float64),
                band_labels_ascending=labels,
                baseline=baseline,
                predictions=predictions,
            )
            reports[group] = {
                "trace_boundaries_descending_one_based": list(observed_boundaries),
                "orthogonality_max_abs": orthogonality_error,
                "eigenvalue_min": float(evals.min()),
                "eigenvalue_median": float(np.median(evals)),
                "eigenvalue_max": float(evals.max()),
                "diagnostics": diagnostics,
                "baseline_sha256": array_sha256(baseline),
                "candidate_array_sha256": array_sha256(predictions),
            }

        atomic_copy(EXP142_TRUTH_LEDGER, staging / TRUTH_LEDGER_ARTIFACT)
        result = {
            "schema_version": 1,
            "stage": "source_only_spectral_trace_bin_prepare",
            "betas": list(BETAS),
            "candidate_count": len(BETAS),
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "spectral_assignment_null_draws": SPECTRAL_NULL_DRAWS,
            "target_alignment_null_draws": TARGET_NULL_DRAWS,
            "development_outcomes_opened": False,
            "sealed_410_effect_values_opened": False,
            "reserved_confirmation_expression_opened": False,
            "truth_arrays_opened": False,
            "null_outcomes_evaluated": False,
            "prepare_authorization_sha256": sha256_file(PREPARE_AUTHORIZATION),
            "run_sha256": sha256_file(RUN_PATH),
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "upstream_140_identity_manifest_sha256": sha256_file(
                EXP140_IDENTITY_MANIFEST
            ),
            "upstream_142_preoutcome_manifest_sha256": sha256_file(
                EXP142_PREOUTCOME_MANIFEST
            ),
            "upstream_142_truth_ledger_sha256": sha256_file(EXP142_TRUTH_LEDGER),
            "null_plan_sha256": {
                key: array_sha256(value) for key, value in null_plans.items()
            },
            "groups": reports,
        }
        atomic_json(staging / RESULT_ARTIFACT, result)
        atomic_json(staging / MANIFEST_ARTIFACT, directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_preoutcome() -> dict[str, Any]:
    validate_directory_manifest(PREOUTCOME_DIR)
    result = read_json(PREOUTCOME_DIR / RESULT_ARTIFACT)
    require(result["development_outcomes_opened"] is False, "prepare opened truth")
    require(
        result["sealed_410_effect_values_opened"] is False,
        "prepare opened sealed outcomes",
    )
    require(result["truth_arrays_opened"] is False, "prepare opened truth arrays")
    require(result["null_outcomes_evaluated"] is False, "prepare evaluated nulls")
    require(result["betas"] == list(BETAS), "beta grid changed")
    require(result["run_sha256"] == sha256_file(RUN_PATH), "run.py changed")
    require(
        result["protocol_sha256"] == sha256_file(PROTOCOL_PATH),
        "protocol changed",
    )
    require(
        result["prepare_authorization_sha256"] == sha256_file(PREPARE_AUTHORIZATION),
        "prepare authorization changed",
    )
    with np.load(PREOUTCOME_DIR / NULL_PLAN_ARTIFACT, allow_pickle=False) as data:
        for key in data.files:
            require(
                array_sha256(data[key]) == result["null_plan_sha256"][key],
                f"null plan changed: {key}",
            )
    return result


def truth_ledger_paths() -> tuple[Path, ...]:
    ledger = read_json(PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT)
    paths = []
    for records in ledger["datasets"].values():
        for record in records.values():
            paths.append(resolve_repo_path(record["path"]))
    return tuple(paths)


def development_control_paths() -> tuple[Path, ...]:
    return (
        RUN_PATH,
        TEST_PATH,
        PROTOCOL_PATH,
        PREPARE_AUTHORIZATION,
        MEMORY_GUARD,
        PREOUTCOME_DIR / MANIFEST_ARTIFACT,
        PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT,
        PREOUTCOME_DIR / NULL_PLAN_ARTIFACT,
        EXP142_DEVELOP,
        EXP142_DEVELOPMENT_AUTHORIZATION,
        EXP132_RUN,
        EXP132_AUTHORIZATION,
        EXP132_PLAN_MANIFEST,
        DEVELOPMENT_TARGETS,
        RESPONSE_GENES,
    )


def validate_development_authorization() -> dict[str, Any]:
    validate_preoutcome()
    require(
        DEVELOPMENT_AUTHORIZATION.is_file(),
        "DEVELOPMENT_ACCESS_AUTHORIZATION.v2.json is required",
    )
    authorization = read_json(DEVELOPMENT_AUTHORIZATION)
    require(
        authorization.get("development_613_authorized") is True,
        "613-target development access is not authorized",
    )
    require(
        authorization.get("source_only_prepare_authorized") is False,
        "development authorization must be stage-specific",
    )
    require_safety_flags(authorization)
    bindings = validate_bound_paths(
        authorization, development_control_paths(), hash_files=True
    )
    ledger = read_json(PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT)
    for records in ledger["datasets"].values():
        for record in records.values():
            path = resolve_repo_path(record["path"])
            require(path.resolve() in bindings, f"truth source not authorized: {path}")
            require(
                bindings[path.resolve()] == record["sha256"],
                f"authorization differs from truth ledger: {path}",
            )
            require(path.is_file(), f"truth source missing: {path}")
            require(path.stat().st_size == record["bytes"], f"bytes changed: {path}")
    return authorization


def create_outcome_access_lock(authorization: Mapping[str, Any]) -> None:
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    atomic_json(
        OUTCOME_ACCESS_LOCK,
        {
            "schema_version": 1,
            "development_targets": DEV_COUNT,
            "sealed_410_effect_values_accessed": False,
            "reserved_confirmation_expression_accessed": False,
            "truth_container_hashes_recomputed": False,
            "truth_access_mode": "selected_rows_and_response_columns_only",
            "null_stage_accessed": False,
            "authorization_sha256": sha256_file(DEVELOPMENT_AUTHORIZATION),
            "authorization_schema": authorization.get("schema_version"),
        },
    )


def load_plan() -> dict[str, np.ndarray]:
    with np.load(PREOUTCOME_DIR / PLAN_ARTIFACT, allow_pickle=False) as data:
        return {name: np.asarray(data[name]) for name in data.files}


def load_candidates(group: str) -> dict[str, np.ndarray]:
    with np.load(
        PREOUTCOME_DIR / f"candidates-{group}.npz", allow_pickle=False
    ) as data:
        result = {name: np.asarray(data[name]) for name in data.files}
    require(
        np.array_equal(result["betas"], np.asarray(BETAS, dtype=np.float64)),
        f"{group}: beta grid changed",
    )
    require(
        result["predictions"].shape == (len(BETAS), DEV_COUNT, RESPONSE_COUNT),
        f"{group}: candidate shape changed",
    )
    require(
        np.array_equal(result["predictions"][0], result["baseline"]),
        f"{group}: beta-zero identity changed",
    )
    return result


def load_development_truth_slice_only(
    parent: ModuleType, source: ModuleType
) -> dict[str, np.ndarray]:
    targets = read_lines(DEVELOPMENT_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    require(len(targets) == DEV_COUNT, "development target count changed")
    require(len(responses) == RESPONSE_COUNT, "response gene count changed")
    truth = {}
    for dataset in SCORED_DATASETS:
        bank = source.Bank.from_root(dataset, parent.TRUTH_ROOTS[dataset])
        truth[dataset] = np.asarray(
            source.load_effect_rows(bank, targets, responses), dtype=np.float64
        )
        require(
            truth[dataset].shape == (DEV_COUNT, RESPONSE_COUNT),
            f"{dataset}: truth shape changed",
        )
        parent.memory_check(f"loaded development truth {dataset}", MAX_RSS_GIB)
    return truth


def normalized_rows(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return np.divide(values, np.maximum(norms, EPS), out=np.zeros_like(values))


def pds_from_square_similarity(similarity: np.ndarray) -> np.ndarray:
    matrix = np.asarray(similarity, dtype=np.float64)
    require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1], "PDS square")
    own = np.diag(matrix)
    greater = np.sum(matrix > own[:, None], axis=1)
    equal = np.sum(matrix == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    denominator = max(len(matrix) - 1.0, 1.0)
    return 1.0 - (midrank - 1.0) / denominator


def build_score_cache(
    candidates: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    parent: ModuleType,
    folds: np.ndarray,
) -> dict[str, Any]:
    target_count = len(folds)
    truth_unit = {
        dataset: normalized_rows(truth[dataset]) for dataset in SCORED_DATASETS
    }
    identity_similarity = {}
    identity_credit = {}
    for dataset in SCORED_DATASETS:
        group = parent.DATASET_GROUP[dataset]
        baseline = candidates[group]["baseline"]
        identity_similarity[dataset] = normalized_rows(baseline) @ truth_unit[dataset].T
        identity_credit[dataset] = parent.signed_credit(baseline, truth[dataset])
    full = {}
    inner = {fold: {} for fold in range(FOLD_COUNT)}
    all_indices = np.arange(target_count, dtype=np.int64)
    for beta_index, beta in enumerate(BETAS):
        full[beta] = {}
        for dataset in SCORED_DATASETS:
            if beta == 0.0:
                similarity = identity_similarity[dataset]
                credit_delta = np.zeros(target_count, dtype=np.float64)
                full_pds_delta = np.zeros(target_count, dtype=np.float64)
            else:
                group = parent.DATASET_GROUP[dataset]
                prediction = candidates[group]["predictions"][beta_index]
                similarity = normalized_rows(prediction) @ truth_unit[dataset].T
                credit_delta = (
                    parent.signed_credit(prediction, truth[dataset])
                    - identity_credit[dataset]
                )
                full_pds_delta = pds_from_square_similarity(
                    similarity
                ) - pds_from_square_similarity(identity_similarity[dataset])
            full[beta][dataset] = {
                "pds": full_pds_delta,
                "signed_credit": credit_delta,
            }
            for fold in range(FOLD_COUNT):
                indices = all_indices[folds != fold]
                candidate_inner = similarity[np.ix_(indices, indices)]
                baseline_inner = identity_similarity[dataset][np.ix_(indices, indices)]
                record = inner[fold].setdefault(beta, {})
                record[dataset] = {
                    "pds": (
                        np.zeros(len(indices), dtype=np.float64)
                        if beta == 0.0
                        else pds_from_square_similarity(candidate_inner)
                        - pds_from_square_similarity(baseline_inner)
                    ),
                    "signed_credit": credit_delta[indices],
                }
    return {"full": full, "inner": inner}


def family_arrays(
    score: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, dict[str, np.ndarray]]:
    return {
        metric: {
            "k562": 0.5
            * (
                np.asarray(score["k562_gwps"][metric], dtype=np.float64)
                + np.asarray(score["k562_essential"][metric], dtype=np.float64)
            ),
            "rpe1": np.asarray(score["rpe1"][metric], dtype=np.float64),
            "jurkat": np.asarray(score["nadig_jurkat"][metric], dtype=np.float64),
        }
        for metric in ("pds", "signed_credit")
    }


def weighted_summary(
    score: Mapping[str, Mapping[str, np.ndarray]],
    weights: np.ndarray | None = None,
) -> dict[str, Any]:
    arrays = family_arrays(score)
    if weights is None:
        weight = np.ones(len(next(iter(arrays["pds"].values()))), dtype=np.float64)
    else:
        weight = np.asarray(weights, dtype=np.float64)
    require(weight.ndim == 1 and np.sum(weight) > 0.0, "invalid weights")
    families = {
        metric: {
            family: float(np.dot(values, weight) / np.sum(weight))
            for family, values in arrays[metric].items()
        }
        for metric in arrays
    }
    per_target = {
        metric: np.mean(np.stack(list(arrays[metric].values())), axis=0)
        for metric in arrays
    }
    return {
        "families": families,
        "macro": {
            metric: float(np.dot(values, weight) / np.sum(weight))
            for metric, values in per_target.items()
        },
        "per_target": per_target,
    }


def prefer_candidate(
    candidate_pds: np.ndarray | float,
    candidate_credit: np.ndarray | float,
    candidate_beta: np.ndarray | float,
    selected_pds: np.ndarray | float,
    selected_credit: np.ndarray | float,
    selected_beta: np.ndarray | float,
) -> np.ndarray:
    candidate_pds_array = np.asarray(candidate_pds)
    candidate_credit_array = np.asarray(candidate_credit)
    candidate_beta_array = np.asarray(candidate_beta)
    selected_pds_array = np.asarray(selected_pds)
    selected_credit_array = np.asarray(selected_credit)
    selected_beta_array = np.asarray(selected_beta)
    return (candidate_pds_array > selected_pds_array) | (
        (candidate_pds_array == selected_pds_array)
        & (
            (candidate_credit_array > selected_credit_array)
            | (
                (candidate_credit_array == selected_credit_array)
                & (candidate_beta_array < selected_beta_array)
            )
        )
    )


def select_beta(
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    weights: np.ndarray | None = None,
) -> tuple[float, dict[float, dict[str, Any]], list[float]]:
    summaries = {beta: weighted_summary(scores[beta], weights) for beta in BETAS}
    eligible = [
        beta
        for beta in BETAS
        if all(
            value >= -EPS
            for metric in ("pds", "signed_credit")
            for value in summaries[beta]["families"][metric].values()
        )
    ]
    require(0.0 in eligible, "beta-zero identity became ineligible")
    selected = 0.0
    for beta in eligible[1:]:
        if bool(
            prefer_candidate(
                summaries[beta]["macro"]["pds"],
                summaries[beta]["macro"]["signed_credit"],
                beta,
                summaries[selected]["macro"]["pds"],
                summaries[selected]["macro"]["signed_credit"],
                selected,
            )
        ):
            selected = beta
    return selected, summaries, eligible


def observed_nested(cache: Mapping[str, Any], folds: np.ndarray) -> dict[str, Any]:
    all_indices = np.arange(len(folds), dtype=np.int64)
    oof = {
        dataset: {
            metric: np.empty(len(folds), dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }
    records = []
    for fold in range(FOLD_COUNT):
        inner_indices = all_indices[folds != fold]
        selected, summaries, eligible = select_beta(cache["inner"][fold])
        held = all_indices[folds == fold]
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                oof[dataset][metric][held] = cache["full"][selected][dataset][metric][
                    held
                ]
        records.append(
            {
                "fold": fold,
                "selected_beta": selected,
                "eligible_betas": eligible,
                "inner_targets": len(inner_indices),
                "selected_inner_macro": summaries[selected]["macro"],
            }
        )
    final_beta, summaries, eligible = select_beta(cache["full"])
    fixed = cache["full"][final_beta]
    return {
        "folds": records,
        "oof": oof,
        "oof_summary": weighted_summary(oof),
        "final_beta": final_beta,
        "final_eligible_betas": eligible,
        "all_candidate_summaries": summaries,
        "fixed": fixed,
        "fixed_summary": weighted_summary(fixed),
    }


def block_count_matrix(draws: np.ndarray, block_count: int) -> np.ndarray:
    plan = np.asarray(draws, dtype=np.int64)
    require(plan.ndim == 2, "bootstrap plan shape")
    output = np.zeros((plan.shape[0], block_count), dtype=np.uint16)
    rows = np.repeat(np.arange(plan.shape[0], dtype=np.int64), plan.shape[1])
    np.add.at(output, (rows, plan.ravel()), 1)
    return output


def block_sums(
    values: np.ndarray,
    blocks: np.ndarray,
    mask: np.ndarray,
    block_count: int,
) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    selected = np.asarray(mask, dtype=bool)
    return np.bincount(
        blocks[selected],
        weights=data[selected],
        minlength=block_count,
    ).astype(np.float64)


def bootstrap_reselected(
    cache: Mapping[str, Any],
    folds: np.ndarray,
    blocks: np.ndarray,
    draws: np.ndarray,
) -> tuple[dict[str, dict[str, float]], dict[str, np.ndarray]]:
    block_count = int(blocks.max()) + 1
    draw_counts = block_count_matrix(draws, block_count).astype(np.float64)
    selected_by_fold = np.zeros((len(draw_counts), FOLD_COUNT), dtype=np.uint8)
    beta_to_index = {beta: index for index, beta in enumerate(BETAS)}
    held_numerators = {
        metric: {
            family: np.zeros(len(draw_counts), dtype=np.float64) for family in FAMILIES
        }
        for metric in ("pds", "signed_credit")
    }
    held_denominators = {
        family: np.zeros(len(draw_counts), dtype=np.float64) for family in FAMILIES
    }

    for fold in range(FOLD_COUNT):
        inner_mask = folds != fold
        held_mask = folds == fold
        inner_block_counts = np.bincount(
            blocks[inner_mask], minlength=block_count
        ).astype(np.float64)
        inner_denominator = draw_counts @ inner_block_counts
        require(np.all(inner_denominator > 0.0), "empty bootstrap inner fold")

        inner_family = {
            beta: family_arrays(cache["inner"][fold][beta]) for beta in BETAS
        }
        draw_family = {
            beta: {
                metric: {
                    family: (
                        draw_counts
                        @ block_sums(
                            inner_family[beta][metric][family],
                            blocks[inner_mask],
                            np.ones(int(inner_mask.sum()), dtype=bool),
                            block_count,
                        )
                    )
                    / inner_denominator
                    for family in FAMILIES
                }
                for metric in ("pds", "signed_credit")
            }
            for beta in BETAS
        }

        selected = np.zeros(len(draw_counts), dtype=np.uint8)
        selected_pds = np.zeros(len(draw_counts), dtype=np.float64)
        selected_credit = np.zeros(len(draw_counts), dtype=np.float64)
        selected_beta = np.zeros(len(draw_counts), dtype=np.float64)
        for beta in BETAS[1:]:
            index = beta_to_index[beta]
            eligible = np.ones(len(draw_counts), dtype=bool)
            for metric in ("pds", "signed_credit"):
                for family in FAMILIES:
                    eligible &= draw_family[beta][metric][family] >= -EPS
            macro_pds = np.mean(
                np.stack([draw_family[beta]["pds"][family] for family in FAMILIES]),
                axis=0,
            )
            macro_credit = np.mean(
                np.stack(
                    [draw_family[beta]["signed_credit"][family] for family in FAMILIES]
                ),
                axis=0,
            )
            better = eligible & prefer_candidate(
                macro_pds,
                macro_credit,
                beta,
                selected_pds,
                selected_credit,
                selected_beta,
            )
            selected[better] = index
            selected_pds[better] = macro_pds[better]
            selected_credit[better] = macro_credit[better]
            selected_beta[better] = beta
        selected_by_fold[:, fold] = selected

        held_block_counts = np.bincount(
            blocks[held_mask], minlength=block_count
        ).astype(np.float64)
        denominator = draw_counts @ held_block_counts
        require(np.all(denominator > 0.0), "empty bootstrap held fold")
        for metric in ("pds", "signed_credit"):
            full_family = {
                beta: family_arrays(cache["full"][beta])[metric] for beta in BETAS
            }
            for family in FAMILIES:
                candidate_numerators = np.stack(
                    [
                        draw_counts
                        @ block_sums(
                            full_family[beta][family],
                            blocks,
                            held_mask,
                            block_count,
                        )
                        for beta in BETAS
                    ],
                    axis=1,
                )
                chosen = np.take_along_axis(
                    candidate_numerators, selected[:, None], axis=1
                )[:, 0]
                held_numerators[metric][family] += chosen
        for family in FAMILIES:
            held_denominators[family] += denominator

    arrays = {}
    for metric in ("pds", "signed_credit"):
        family_values = []
        for family in FAMILIES:
            value = held_numerators[metric][family] / held_denominators[family]
            arrays[f"{metric}__{family}"] = value
            family_values.append(value)
        arrays[metric] = np.mean(np.stack(family_values), axis=0)
    arrays["selected_beta_index_by_fold"] = selected_by_fold
    summary = {
        metric: {
            "lower_95": float(np.quantile(arrays[metric], 0.025)),
            "upper_95": float(np.quantile(arrays[metric], 0.975)),
            "mean": float(np.mean(arrays[metric])),
        }
        for metric in ("pds", "signed_credit")
    }
    return summary, arrays


def fold_family_cells(
    oof: Mapping[str, Mapping[str, np.ndarray]],
    folds: np.ndarray,
) -> dict[str, dict[str, list[float]]]:
    arrays = family_arrays(oof)
    return {
        metric: {
            family: [float(values[folds == fold].mean()) for fold in range(FOLD_COUNT)]
            for family, values in arrays[metric].items()
        }
        for metric in arrays
    }


def concentration_summary(
    oof: Mapping[str, Mapping[str, np.ndarray]],
    folds: np.ndarray,
    blocks: np.ndarray,
) -> dict[str, dict[str, float]]:
    arrays = family_arrays(oof)
    result = {}
    for metric in ("pds", "signed_credit"):
        family_gain = np.asarray(
            [max(0.0, float(arrays[metric][family].mean())) for family in FAMILIES]
        )
        per_target = np.mean(np.stack(list(arrays[metric].values())), axis=0)
        block_gain = np.asarray(
            [
                max(0.0, float(per_target[blocks == block].sum()))
                for block in np.unique(blocks)
            ]
        )
        fold_gain = np.asarray(
            [
                max(0.0, float(per_target[folds == fold].sum()))
                for fold in range(FOLD_COUNT)
            ]
        )

        def fraction(values: np.ndarray) -> float:
            total = float(values.sum())
            return float(values.max() / total) if total > 0.0 else 1.0

        result[metric] = {
            "family": fraction(family_gain),
            "target_block": fraction(block_gain),
            "outer_fold": fraction(fold_gain),
        }
    return result


def promotion_checks(
    nested: Mapping[str, Any],
    cells: Mapping[str, Mapping[str, Sequence[float]]],
    bootstrap: Mapping[str, Mapping[str, float]],
    concentration: Mapping[str, Mapping[str, float]],
) -> dict[str, bool]:
    checks = {
        "final_beta_nonzero": float(nested["final_beta"]) > 0.0,
        "nested_macro_pds_at_least_0.003": (
            nested["oof_summary"]["macro"]["pds"] >= MIN_GAIN
        ),
        "nested_macro_signed_credit_at_least_0.003": (
            nested["oof_summary"]["macro"]["signed_credit"] >= MIN_GAIN
        ),
        "fixed_macro_pds_at_least_0.003": (
            nested["fixed_summary"]["macro"]["pds"] >= MIN_GAIN
        ),
        "fixed_macro_signed_credit_at_least_0.003": (
            nested["fixed_summary"]["macro"]["signed_credit"] >= MIN_GAIN
        ),
        "pds_bootstrap_lower_positive": bootstrap["pds"]["lower_95"] > 0.0,
        "signed_credit_bootstrap_lower_positive": (
            bootstrap["signed_credit"]["lower_95"] > 0.0
        ),
    }
    for metric in ("pds", "signed_credit"):
        checks[f"{metric}_every_family_nonnegative"] = all(
            value >= -EPS
            for value in nested["oof_summary"]["families"][metric].values()
        )
        checks[f"{metric}_every_family_fold_cell_nonnegative"] = all(
            value >= -EPS for family in FAMILIES for value in cells[metric][family]
        )
        for kind in ("family", "target_block", "outer_fold"):
            checks[f"{metric}_{kind}_concentration_at_most_0.5"] = (
                concentration[metric][kind] <= MAX_CONCENTRATION
            )
    return checks


def develop() -> None:
    authorization = validate_development_authorization()
    require(not RESULT_DIR.exists(), f"development result exists: {RESULT_DIR}")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    with exclusive_lock(DEVELOP_LOCK):
        with process_tree_memory_guard() as memory_report:
            create_outcome_access_lock(authorization)
            parent = load_exp142()
            source = load_module(EXP132_RUN, "experiment132_for_145")
            truth = load_development_truth_slice_only(parent, source)
            plan = load_plan()
            folds = np.asarray(plan["development_folds"], dtype=np.int8)
            blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
            draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.int64)
            validate_plan_arrays(folds, blocks, draws)
            candidates = {group: load_candidates(group) for group in GROUPS}
            cache = build_score_cache(candidates, truth, parent, folds)
            nested = observed_nested(cache, folds)
            cells = fold_family_cells(nested["oof"], folds)
            bootstrap, bootstrap_arrays = bootstrap_reselected(
                cache, folds, blocks, draws
            )
            concentration = concentration_summary(nested["oof"], folds, blocks)
            checks = promotion_checks(nested, cells, bootstrap, concentration)
            passes = all(checks.values())
        result = {
            "schema_version": 1,
            "stage": "authorized_613_target_spectral_trace_bins",
            "candidate_family": "post_whitening_spectral_trace_bin_attenuation",
            "folds": nested["folds"],
            "final_beta": nested["final_beta"],
            "final_eligible_betas": nested["final_eligible_betas"],
            "nested": {
                "macro": nested["oof_summary"]["macro"],
                "families": nested["oof_summary"]["families"],
            },
            "fixed": {
                "macro": nested["fixed_summary"]["macro"],
                "families": nested["fixed_summary"]["families"],
            },
            "all_candidate_summaries": {
                str(beta): {
                    "macro": summary["macro"],
                    "families": summary["families"],
                }
                for beta, summary in nested["all_candidate_summaries"].items()
            },
            "fold_family_cells": cells,
            "bootstrap_reselected": bootstrap,
            "positive_gain_concentration": concentration,
            "checks": checks,
            "passes": passes,
            "decision": ("PASS_FOR_FROZEN_NULL_STAGE" if passes else "TERMINAL_NO_GO"),
            "null_plans_frozen_before_truth": True,
            "null_stage_run": False,
            "development_613_outcomes_opened": True,
            "sealed_410_effect_values_opened": False,
            "reserved_confirmation_expression_opened": False,
            "truth_container_hashes_recomputed": False,
            "truth_access_mode": "selected_rows_and_response_columns_only",
            "memory_guard": memory_report,
            "count_generation": False,
            "cell_eval2": False,
            "challenge_build": False,
            "submission": False,
            "development_authorization_sha256": sha256_file(DEVELOPMENT_AUTHORIZATION),
            "preoutcome_manifest_sha256": sha256_file(
                PREOUTCOME_DIR / MANIFEST_ARTIFACT
            ),
        }
        with atomic_stage_directory(RESULT_DIR) as staging:
            arrays = {}
            for dataset in SCORED_DATASETS:
                for metric in ("pds", "signed_credit"):
                    arrays[f"oof__{dataset}__{metric}"] = nested["oof"][dataset][metric]
                    arrays[f"fixed__{dataset}__{metric}"] = nested["fixed"][dataset][
                        metric
                    ]
            for name, values in bootstrap_arrays.items():
                arrays[f"bootstrap__{name}"] = values
            atomic_npz(staging / ARRAY_ARTIFACT, **arrays)
            atomic_json(staging / RESULT_ARTIFACT, result)
            atomic_json(staging / MANIFEST_ARTIFACT, directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    protocol()
    payload: dict[str, Any] = {
        "schema_version": 1,
        "experiment": HERE.name,
        "candidate_count": len(BETAS),
        "prepare_authorization_exists": PREPARE_AUTHORIZATION.exists(),
        "preoutcome_exists": PREOUTCOME_DIR.exists(),
        "development_authorization_exists": DEVELOPMENT_AUTHORIZATION.exists(),
        "outcome_access_recorded": OUTCOME_ACCESS_LOCK.exists(),
        "development_result_exists": RESULT_DIR.exists(),
        "truth_arrays_opened_by_audit": False,
        "sealed_test_loader_exists": False,
        "null_stage_exists": False,
    }
    if PREPARE_AUTHORIZATION.exists():
        validate_prepare_authorization()
        payload["prepare_authorization_valid"] = True
    if PREOUTCOME_DIR.exists():
        validate_preoutcome()
        payload["preoutcome_valid_without_truth_read"] = True
    if DEVELOPMENT_AUTHORIZATION.exists() and PREOUTCOME_DIR.exists():
        validate_development_authorization()
        payload["development_authorization_valid_without_truth_read"] = True
    if RESULT_DIR.exists():
        validate_directory_manifest(RESULT_DIR)
        result = read_json(RESULT_DIR / RESULT_ARTIFACT)
        require(
            result["sealed_410_effect_values_opened"] is False,
            "result reports sealed outcome access",
        )
        payload["development_result"] = {
            "passes": result["passes"],
            "decision": result["decision"],
        }
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "develop", "audit"))
    args = parser.parse_args()
    if args.command == "prepare":
        prepare()
    elif args.command == "develop":
        develop()
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
