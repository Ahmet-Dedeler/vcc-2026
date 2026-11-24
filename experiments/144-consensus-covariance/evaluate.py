#!/usr/bin/env python3
"""Frozen development evaluation for Experiment 144.

``prepare`` constructs the single fixed consensus-covariance candidate and
freezes bootstrap/null plans without opening development truth. ``evaluate``
performs one separately authorized 613-target evaluation. ``audit`` never
opens biological outcome arrays.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = Path(__file__).resolve()
TEST_PATH = HERE / "test_evaluate.py"
PROTOCOL_PATH = HERE / "development_protocol.json"
CORE_PATH = HERE / "run.py"
MEMORY_GUARD = HERE / "memory_guard.py"

PREPARE_AUTHORIZATION = HERE / "DEVELOPMENT_PREFREEZE_AUTHORIZATION.v3.json"
DEVELOPMENT_AUTHORIZATION = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
PREOUTCOME_DIR = HERE / "work/development-preoutcome-v2"
RESULT_DIR = HERE / "work/development-result-v2"
PREPARE_LOCK = HERE / "work/development-preoutcome-v2.lock"
DEVELOP_LOCK = HERE / "work/development-result-v2.lock"
OUTCOME_ACCESS_LOCK = HERE / "work/development-outcome-access.v1.lock"

SOURCE_AUTHORIZATION = HERE / "SOURCE_ONLY_PREPARE_AUTHORIZATION.v2.json"
SOURCE_DIR = HERE / "work/source-only-v2"
SOURCE_MANIFEST = SOURCE_DIR / "artifact-manifest.json"
SOURCE_RESULT = SOURCE_DIR / "result.json"
SOURCE_MEMORY_REPORT = HERE / "work/source-only-v2-memory.json"

EXP140 = ROOT / "experiments/140-batch-matched-controls"
EXP140_IDENTITY = EXP140 / "work/identity-reproduction-v5"
EXP140_MANIFEST = EXP140_IDENTITY / "artifact-manifest.json"

EXP142 = ROOT / "experiments/142-batch-reliability-weighting"
EXP142_DEVELOP = EXP142 / "develop.py"
EXP142_DEVELOPMENT_AUTHORIZATION = EXP142 / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
EXP142_PREOUTCOME = EXP142 / "work/development-preoutcome-v1"
EXP142_MANIFEST = EXP142_PREOUTCOME / "artifact-manifest.json"
EXP142_PLAN = EXP142_PREOUTCOME / "development-plan.npz"
EXP142_TRUTH_LEDGER = EXP142_PREOUTCOME / "truth-source-hashes.json"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
EXP132_RUN = EXP132 / "run.py"
EXP132_AUTHORIZATION = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
EXP132_PLAN = EXP132 / "work/preregistered-freeze-v2"
EXP132_PLAN_MANIFEST = EXP132_PLAN / "plan-manifest.json"
DEVELOPMENT_TARGETS = EXP132_PLAN / "selection_targets.txt"
RESPONSE_GENES = EXP132_PLAN / "response_genes.txt"

COMMANDS = ("prepare", "evaluate", "audit")
GROUPS = ("k562", "non_k562")
SCORED_DATASETS = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat")
FAMILIES = ("k562", "rpe1", "jurkat")
DATASET_GROUP = {
    "k562_gwps": "k562",
    "k562_essential": "k562",
    "rpe1": "non_k562",
    "nadig_jurkat": "non_k562",
}

CANDIDATE_NAME = "consensus:g=0.65:eta=0.25"
OLD_GAMMA = 0.65
CONSENSUS_GAMMA = 0.65
BLEND_ETA = 0.25
DEV_COUNT = 613
RESPONSE_COUNT = 4_992
FOLD_COUNT = 5
EXPECTED_FOLD_SIZES = (123, 123, 123, 122, 122)
BLOCK_COUNT = 531
BOOTSTRAP_DRAWS = 20_000
NULL_DRAWS = 99
NULL_SEED = 144_202_603
MAX_RSS_GIB = 8.0
MIN_FREE_PERCENT = 25.0
MIN_PDS_GAIN = 0.003
MIN_SIGNED_CREDIT_GAIN = 0.0
MIN_FOLD_PDS = -0.0005
EPS = 1e-12

PLAN_NAME = "development-plan.npz"
TRUTH_LEDGER_NAME = "truth-source-hashes.json"
RESULT_NAME = "result.json"
ARRAYS_NAME = "per-target-bootstrap-null.npz"
MANIFEST_NAME = "artifact-manifest.json"

PREPARE_CONFIRMATION = "AUTHORIZE_DEVELOPMENT_PREFREEZE_V3"
DEVELOP_CONFIRMATION = "AUTHORIZE_613_TARGET_DEVELOPMENT_EVALUATION_V1"


def load_module(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


CORE = load_module(CORE_PATH, "experiment144_source_core")
ContractError = CORE.ContractError
require = CORE.require


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


def protocol() -> dict[str, Any]:
    value = read_json(PROTOCOL_PATH)
    require(value.get("schema_version") == 1, "development protocol schema changed")
    require(
        value.get("experiment") == HERE.name, "development protocol experiment changed"
    )
    require(
        value.get("status") == "UNAUTHORIZED", "development protocol status changed"
    )
    candidate = value["candidate"]
    require(candidate["name"] == CANDIDATE_NAME, "candidate name changed")
    require(candidate["new_grid_or_selection"] is False, "candidate selection widened")
    require(
        candidate["old_operator_inverse_gamma"] == -OLD_GAMMA, "old inverse changed"
    )
    require(
        candidate["consensus_operator_gamma"] == CONSENSUS_GAMMA,
        "consensus gamma changed",
    )
    require(candidate["blend_eta"] == BLEND_ETA, "blend eta changed")
    null = value["target_alignment_null"]
    require(null["draws"] == NULL_DRAWS, "null draw count changed")
    require(null["seed"] == NULL_SEED, "null seed changed")
    return value


def resolve_repo_path(value: str) -> Path:
    path = (ROOT / value).resolve()
    require(ROOT == path or ROOT in path.parents, f"path escapes repository: {value}")
    return path


def authorization_bindings(authorization: Mapping[str, Any]) -> dict[Path, str]:
    raw = authorization.get("bindings")
    require(isinstance(raw, Mapping), "authorization lacks bindings")
    bindings = {}
    for name, digest in raw.items():
        path = Path(str(name))
        resolved = (
            path.resolve() if path.is_absolute() else resolve_repo_path(str(name))
        )
        require(isinstance(digest, str), f"invalid authorization binding: {name}")
        require(resolved not in bindings, f"duplicate authorization binding: {name}")
        bindings[resolved] = digest
    return bindings


def validate_upstream_manifest(
    directory: Path,
    guard: Any,
    *,
    expected_schema_version: int | None,
) -> dict[str, Any]:
    manifest_path = directory / MANIFEST_NAME
    require(manifest_path.is_file(), f"manifest missing: {directory}")
    manifest = read_json(manifest_path)
    require(
        manifest.get("schema_version") == expected_schema_version,
        f"manifest schema changed: {directory}",
    )
    expected = set(manifest["files"])
    actual = {
        path.name
        for path in directory.iterdir()
        if path.is_file() and path.name != MANIFEST_NAME
    }
    require(actual == expected, f"artifact set changed: {directory}")
    for name, record in manifest["files"].items():
        path = directory / name
        require(path.stat().st_size == record["bytes"], f"bytes changed: {path}")
        require(
            CORE.sha256_file(path, guard) == record["sha256"],
            f"hash changed: {path}",
        )
    return manifest


def source_package_files() -> tuple[Path, ...]:
    return (
        SOURCE_DIR / "candidate-spec.json",
        SOURCE_DIR / "consensus-operator-k562.npz",
        SOURCE_DIR / "consensus-operator-k562.basis.json",
        SOURCE_DIR / "consensus-operator-non_k562.npz",
        SOURCE_DIR / "consensus-operator-non_k562.basis.json",
        SOURCE_RESULT,
        SOURCE_MANIFEST,
    )


def prepare_control_paths() -> tuple[Path, ...]:
    return (
        RUN_PATH,
        TEST_PATH,
        PROTOCOL_PATH,
        CORE_PATH,
        HERE / "test_run.py",
        HERE / "preregistered_protocol.json",
        MEMORY_GUARD,
        SOURCE_AUTHORIZATION,
        SOURCE_MEMORY_REPORT,
        *source_package_files(),
        EXP140_MANIFEST,
        EXP140_IDENTITY / "whitener-k562.npz",
        EXP140_IDENTITY / "whitener-non_k562.npz",
        EXP142_MANIFEST,
        EXP142_PLAN,
        EXP142_PREOUTCOME / "geometry-k562.npz",
        EXP142_PREOUTCOME / "geometry-non_k562.npz",
        EXP142_TRUTH_LEDGER,
        DEVELOPMENT_TARGETS,
        RESPONSE_GENES,
    )


def require_safety_flags(
    authorization: Mapping[str, Any],
    *,
    development_authorized: bool,
) -> None:
    require(
        authorization.get("development_613_authorized") is development_authorized,
        "development authorization flag changed",
    )
    required_false = (
        "untouched_410_test_authorized",
        "reserved_confirmation_expression_authorized",
        "count_generation_authorized",
        "cell_eval2_authorized",
        "challenge_build_authorized",
        "submission_authorized",
        "commit_authorized",
        "push_authorized",
    )
    for key in required_false:
        require(authorization.get(key) is False, f"unsafe authorization flag: {key}")


def validate_control_bindings(
    authorization: Mapping[str, Any],
    required: Sequence[Path],
    guard: Any,
    *,
    exact: bool = False,
) -> None:
    bindings = authorization_bindings(authorization)
    required_paths = {path.resolve() for path in required}
    if exact:
        require(set(bindings) == required_paths, "authorization binding set changed")
    else:
        require(
            required_paths.issubset(bindings), "authorization misses control bindings"
        )
    for path in required:
        resolved = path.resolve()
        require(resolved.is_file(), f"bound control file missing: {path}")
        require(
            CORE.sha256_file(resolved, guard) == bindings[resolved],
            f"control binding changed: {path}",
        )


def validate_prepare_authorization(guard: Any) -> dict[str, Any]:
    require(PREPARE_AUTHORIZATION.is_file(), "development prefreeze is not authorized")
    authorization = read_json(PREPARE_AUTHORIZATION)
    require(authorization.get("schema_version") == 1, "prepare authorization schema")
    require(
        authorization.get("experiment") == HERE.name, "prepare authorization experiment"
    )
    require(
        authorization.get("stage") == "development_prefreeze_v3",
        "prepare authorization stage",
    )
    require(
        authorization.get("confirmation") == PREPARE_CONFIRMATION,
        "prepare authorization confirmation",
    )
    require(
        authorization.get("development_prefreeze_authorized") is True,
        "development prefreeze is not authorized",
    )
    require_safety_flags(authorization, development_authorized=False)
    validate_control_bindings(
        authorization,
        prepare_control_paths(),
        guard,
        exact=True,
    )
    return authorization


def truth_ledger_records() -> tuple[dict[str, Any], ...]:
    ledger = read_json(PREOUTCOME_DIR / TRUTH_LEDGER_NAME)
    records = []
    for dataset in SCORED_DATASETS:
        for record in ledger["datasets"][dataset].values():
            records.append(record)
    records.extend(
        (
            ledger["development_targets"],
            ledger["response_genes"],
            ledger["source_authorization"],
            ledger["source_plan_manifest"],
        )
    )
    return tuple(records)


def development_control_paths() -> tuple[Path, ...]:
    return (
        RUN_PATH,
        TEST_PATH,
        PROTOCOL_PATH,
        CORE_PATH,
        MEMORY_GUARD,
        PREPARE_AUTHORIZATION,
        PREOUTCOME_DIR / MANIFEST_NAME,
        PREOUTCOME_DIR / RESULT_NAME,
        PREOUTCOME_DIR / PLAN_NAME,
        PREOUTCOME_DIR / TRUTH_LEDGER_NAME,
        PREOUTCOME_DIR / "candidates-k562.npz",
        PREOUTCOME_DIR / "candidates-non_k562.npz",
        EXP142_DEVELOP,
        EXP142_DEVELOPMENT_AUTHORIZATION,
        EXP132_RUN,
        EXP132_AUTHORIZATION,
        EXP132_PLAN_MANIFEST,
        DEVELOPMENT_TARGETS,
        RESPONSE_GENES,
    )


def validate_development_authorization(guard: Any) -> dict[str, Any]:
    require(
        DEVELOPMENT_AUTHORIZATION.is_file(), "development evaluation is not authorized"
    )
    authorization = read_json(DEVELOPMENT_AUTHORIZATION)
    require(
        authorization.get("schema_version") == 1, "development authorization schema"
    )
    require(
        authorization.get("experiment") == HERE.name,
        "development authorization experiment",
    )
    require(
        authorization.get("stage") == "development_613_evaluate_v1",
        "development authorization stage",
    )
    require(
        authorization.get("confirmation") == DEVELOP_CONFIRMATION,
        "development authorization confirmation",
    )
    require(
        authorization.get("development_prefreeze_authorized") is False,
        "development authorization must be stage-specific",
    )
    require_safety_flags(authorization, development_authorized=True)
    bindings = authorization_bindings(authorization)
    controls = {path.resolve() for path in development_control_paths()}
    truth_records = truth_ledger_records()
    truth_paths = {resolve_repo_path(record["path"]) for record in truth_records}
    require(set(bindings) == controls | truth_paths, "development binding set changed")
    validate_control_bindings(authorization, development_control_paths(), guard)
    for record in truth_records:
        path = resolve_repo_path(record["path"])
        require(bindings[path] == record["sha256"], f"truth ledger mismatch: {path}")
        require(path.is_file(), f"truth source missing: {path}")
        require(path.stat().st_size == record["bytes"], f"truth bytes changed: {path}")
    return authorization


def validate_plan(
    folds: np.ndarray,
    blocks: np.ndarray,
    draws: np.ndarray,
) -> None:
    require(folds.shape == (DEV_COUNT,), "fold vector shape")
    require(blocks.shape == (DEV_COUNT,), "block vector shape")
    require(
        tuple(int(np.sum(folds == fold)) for fold in range(FOLD_COUNT))
        == EXPECTED_FOLD_SIZES,
        "fold sizes changed",
    )
    require(
        np.array_equal(np.unique(blocks), np.arange(BLOCK_COUNT)), "block ids changed"
    )
    require(draws.shape == (BOOTSTRAP_DRAWS, BLOCK_COUNT), "bootstrap shape changed")
    require(np.all((draws >= 0) & (draws < BLOCK_COUNT)), "bootstrap range changed")
    for block in range(BLOCK_COUNT):
        require(
            len(np.unique(folds[blocks == block])) == 1,
            f"target block {block} crosses folds",
        )


def stable_derangement(
    rng: np.random.Generator,
    values: np.ndarray,
) -> np.ndarray:
    source = np.asarray(values, dtype=np.int64)
    require(len(source) >= 2, "cannot derange fewer than two targets")
    for _ in range(10_000):
        candidate = source[rng.permutation(len(source))]
        if np.all(candidate != source):
            return candidate
    raise ContractError("could not construct target derangement")


def make_null_permutations(folds: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(NULL_SEED)
    indices = np.arange(DEV_COUNT, dtype=np.int64)
    output = np.empty((NULL_DRAWS, DEV_COUNT), dtype=np.int16)
    for draw in range(NULL_DRAWS):
        permutation = np.empty(DEV_COUNT, dtype=np.int64)
        for fold in range(FOLD_COUNT):
            selected = indices[folds == fold]
            permutation[selected] = stable_derangement(rng, selected)
        require(np.all(permutation != indices), "null permutation has fixed points")
        require(
            np.array_equal(folds[permutation], folds), "null permutation crossed folds"
        )
        output[draw] = permutation.astype(np.int16)
    return output


def row_angles_degrees(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    a = normalized_rows(first)
    b = normalized_rows(second)
    cosine = np.sum(a * b, axis=1)
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def atomic_copy(source: Path, destination: Path, guard: Any) -> None:
    require(not destination.exists(), f"refuse overwrite: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    try:
        guard.check(f"before copy {source.name}")
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
        guard.check(f"after copy {source.name}")
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def build_candidate(
    group: str, guard: Any
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    with np.load(
        EXP142_PREOUTCOME / f"geometry-{group}.npz",
        allow_pickle=False,
    ) as geometry:
        baseline = np.asarray(geometry["baseline"], dtype=np.float32)
    require(baseline.shape == (DEV_COUNT, RESPONSE_COUNT), f"{group}: baseline shape")
    with np.load(
        EXP140_IDENTITY / f"whitener-{group}.npz",
        allow_pickle=False,
    ) as old_operator:
        old_evals = np.asarray(old_operator["evals"], dtype=np.float64)
        old_evecs = np.asarray(old_operator["evecs"], dtype=np.float32)
    with np.load(
        SOURCE_DIR / f"consensus-operator-{group}.npz",
        allow_pickle=False,
    ) as consensus_operator:
        consensus_evals = np.asarray(consensus_operator["evals"], dtype=np.float64)
        consensus_evecs = np.asarray(consensus_operator["evecs"], dtype=np.float32)
    for name, values in (
        ("old evals", old_evals),
        ("consensus evals", consensus_evals),
    ):
        require(values.shape == (RESPONSE_COUNT,), f"{group}: {name} shape")
        require(
            np.all(np.isfinite(values)) and np.all(values > 0.0), f"{group}: {name}"
        )
    for name, values in (
        ("old evecs", old_evecs),
        ("consensus evecs", consensus_evecs),
    ):
        require(
            values.shape == (RESPONSE_COUNT, RESPONSE_COUNT),
            f"{group}: {name} shape",
        )
    guard.check(f"{group} candidate inputs")
    pre_whitening = CORE.apply_operator(
        baseline,
        old_evals,
        old_evecs,
        -OLD_GAMMA,
    )
    reconstructed = CORE.apply_operator(
        pre_whitening,
        old_evals,
        old_evecs,
        OLD_GAMMA,
    )
    reconstruction_max_abs = float(np.max(np.abs(reconstructed - baseline)))
    reconstruction_angle = row_angles_degrees(baseline, reconstructed)
    require(reconstruction_max_abs <= 2e-6, f"{group}: inverse reconstruction error")
    require(
        float(reconstruction_angle.max()) <= 1e-4,
        f"{group}: inverse reconstruction angle",
    )
    consensus = CORE.apply_operator(
        pre_whitening,
        consensus_evals,
        consensus_evecs,
        CONSENSUS_GAMMA,
    )
    candidate = CORE.row_norm_preserving_blend(baseline, consensus, BLEND_ETA)
    baseline_norm = np.linalg.norm(baseline, axis=1)
    candidate_norm = np.linalg.norm(candidate, axis=1)
    norm_error = np.abs(candidate_norm - baseline_norm) / np.maximum(baseline_norm, EPS)
    require(np.all(np.isfinite(candidate)), f"{group}: nonfinite candidate")
    require(float(norm_error.max()) <= 2e-6, f"{group}: row norm drift")
    angles = row_angles_degrees(baseline, candidate)
    guard.check(f"{group} candidate complete")
    return (
        baseline,
        candidate,
        {
            "baseline_sha256": CORE.array_sha256(baseline),
            "candidate_sha256": CORE.array_sha256(candidate),
            "inverse_reconstruction_max_abs": reconstruction_max_abs,
            "inverse_reconstruction_max_angle_degrees": float(
                reconstruction_angle.max()
            ),
            "candidate_median_angle_degrees": float(np.median(angles)),
            "candidate_p95_angle_degrees": float(np.quantile(angles, 0.95)),
            "candidate_max_angle_degrees": float(angles.max()),
            "candidate_max_relative_row_norm_error": float(norm_error.max()),
        },
    )


def prepare(max_rss_gib: float, min_free_percent: float) -> None:
    protocol()
    require(not PREOUTCOME_DIR.exists(), f"preoutcome exists: {PREOUTCOME_DIR}")
    require(not RESULT_DIR.exists(), f"development result exists: {RESULT_DIR}")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    effective_max = min(float(max_rss_gib), MAX_RSS_GIB)
    effective_floor = max(float(min_free_percent), MIN_FREE_PERCENT)
    external = CORE.validate_external_guard(effective_max, effective_floor)
    guard = CORE.MemoryGuard(effective_max, effective_floor)
    validate_prepare_authorization(guard)
    CORE.validate_directory_manifest(SOURCE_DIR, guard)
    validate_upstream_manifest(
        EXP140_IDENTITY,
        guard,
        expected_schema_version=None,
    )
    validate_upstream_manifest(
        EXP142_PREOUTCOME,
        guard,
        expected_schema_version=1,
    )

    with np.load(EXP142_PLAN, allow_pickle=False) as plan:
        folds = np.asarray(plan["development_folds"], dtype=np.int8)
        blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
        draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.uint16)
    validate_plan(folds, blocks, draws)
    null_permutations = make_null_permutations(folds)

    with (
        CORE.exclusive_lock(PREPARE_LOCK),
        CORE.atomic_stage_directory(PREOUTCOME_DIR, guard) as staging,
    ):
        CORE.atomic_npz(
            staging / PLAN_NAME,
            guard,
            development_folds=folds,
            development_blocks=blocks,
            bootstrap_block_draws=draws,
            target_alignment_permutations=null_permutations,
        )
        group_reports = {}
        for group in GROUPS:
            baseline, candidate, report = build_candidate(group, guard)
            CORE.atomic_npz(
                staging / f"candidates-{group}.npz",
                guard,
                baseline=baseline,
                candidate=candidate,
            )
            group_reports[group] = report
        atomic_copy(EXP142_TRUTH_LEDGER, staging / TRUTH_LEDGER_NAME, guard)
        result = {
            "schema_version": 1,
            "stage": "development_prefreeze_v3",
            "status": "PREPARED_PENDING_INDEPENDENT_NO_TRUTH_AUDIT",
            "candidate": CANDIDATE_NAME,
            "candidate_count": 2,
            "new_grid_or_selection": False,
            "groups": group_reports,
            "fold_sizes": list(EXPECTED_FOLD_SIZES),
            "target_blocks": BLOCK_COUNT,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "target_alignment_null_draws": NULL_DRAWS,
            "target_alignment_null_sha256": CORE.array_sha256(null_permutations),
            "development_outcomes_opened": False,
            "sealed_410_effect_values_opened": False,
            "reserved_confirmation_expression_opened": False,
            "truth_arrays_opened": False,
            "count_generation": False,
            "cell_eval2": False,
            "challenge_build": False,
            "submission": False,
            "prepare_authorization_sha256": CORE.sha256_file(
                PREPARE_AUTHORIZATION,
                guard,
            ),
            "source_manifest_sha256": CORE.sha256_file(SOURCE_MANIFEST, guard),
            "protocol_sha256": CORE.sha256_file(PROTOCOL_PATH, guard),
            "run_sha256": CORE.sha256_file(RUN_PATH, guard),
            "external_memory_guard": external,
            "memory": guard.report(),
        }
        CORE.atomic_json(staging / RESULT_NAME, result, guard)
        CORE.atomic_json(
            staging / MANIFEST_NAME,
            CORE.directory_manifest(staging, guard),
            guard,
        )
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_preoutcome(guard: Any) -> dict[str, Any]:
    CORE.validate_directory_manifest(PREOUTCOME_DIR, guard)
    result = read_json(PREOUTCOME_DIR / RESULT_NAME)
    require(
        result["status"] == "PREPARED_PENDING_INDEPENDENT_NO_TRUTH_AUDIT",
        "prefreeze status changed",
    )
    require(result["development_outcomes_opened"] is False, "prefreeze opened truth")
    require(result["truth_arrays_opened"] is False, "prefreeze opened truth arrays")
    require(result["candidate"] == CANDIDATE_NAME, "candidate changed")
    require(result["new_grid_or_selection"] is False, "candidate family widened")
    with np.load(PREOUTCOME_DIR / PLAN_NAME, allow_pickle=False) as plan:
        folds = np.asarray(plan["development_folds"], dtype=np.int8)
        blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
        draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.uint16)
        permutations = np.asarray(
            plan["target_alignment_permutations"],
            dtype=np.int16,
        )
    validate_plan(folds, blocks, draws)
    require(
        permutations.shape == (NULL_DRAWS, DEV_COUNT),
        "null permutation shape changed",
    )
    require(
        CORE.array_sha256(permutations) == result["target_alignment_null_sha256"],
        "null permutations changed",
    )
    for group in GROUPS:
        with np.load(
            PREOUTCOME_DIR / f"candidates-{group}.npz",
            allow_pickle=False,
        ) as data:
            baseline = np.asarray(data["baseline"], dtype=np.float32)
            candidate = np.asarray(data["candidate"], dtype=np.float32)
        require(
            baseline.shape == (DEV_COUNT, RESPONSE_COUNT),
            f"{group}: baseline shape changed",
        )
        require(candidate.shape == baseline.shape, f"{group}: candidate shape changed")
        require(
            CORE.array_sha256(baseline) == result["groups"][group]["baseline_sha256"],
            f"{group}: baseline changed",
        )
        require(
            CORE.array_sha256(candidate) == result["groups"][group]["candidate_sha256"],
            f"{group}: candidate changed",
        )
    return result


def create_outcome_access_lock(
    authorization: Mapping[str, Any],
    guard: Any,
) -> None:
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    CORE.atomic_json(
        OUTCOME_ACCESS_LOCK,
        {
            "schema_version": 1,
            "experiment": HERE.name,
            "development_targets": DEV_COUNT,
            "sealed_410_effect_values_accessed": False,
            "reserved_confirmation_expression_accessed": False,
            "truth_container_hashes_recomputed": False,
            "truth_access_mode": "selected 613 rows by 4992 response columns only",
            "authorization_sha256": CORE.sha256_file(
                DEVELOPMENT_AUTHORIZATION,
                guard,
            ),
            "authorization_schema": authorization.get("schema_version"),
        },
        guard,
    )


def load_truth(guard: Any) -> tuple[dict[str, np.ndarray], ModuleType]:
    parent = load_module(EXP142_DEVELOP, "experiment142_for_experiment144")
    source = load_module(EXP132_RUN, "experiment132_for_experiment144")
    targets = read_lines(DEVELOPMENT_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    require(len(targets) == DEV_COUNT, "development target count changed")
    require(len(responses) == RESPONSE_COUNT, "response count changed")
    truth = {}
    for dataset in SCORED_DATASETS:
        bank = source.Bank.from_root(dataset, parent.TRUTH_ROOTS[dataset])
        truth[dataset] = np.asarray(
            source.load_effect_rows(bank, targets, responses),
            dtype=np.float64,
        )
        require(
            truth[dataset].shape == (DEV_COUNT, RESPONSE_COUNT),
            f"{dataset}: truth shape changed",
        )
        guard.check(f"loaded authorized truth {dataset}")
    return truth, parent


def normalized_rows(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return np.divide(values, np.maximum(norms, EPS), out=np.zeros_like(values))


def pds_from_square_similarity(similarity: np.ndarray) -> np.ndarray:
    matrix = np.asarray(similarity, dtype=np.float64)
    require(matrix.shape == (DEV_COUNT, DEV_COUNT), "PDS gallery shape")
    own = np.diag(matrix)
    greater = np.sum(matrix > own[:, None], axis=1)
    equal = np.sum(matrix == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (DEV_COUNT - 1.0)


def load_candidates() -> dict[str, dict[str, np.ndarray]]:
    output = {}
    for group in GROUPS:
        with np.load(
            PREOUTCOME_DIR / f"candidates-{group}.npz",
            allow_pickle=False,
        ) as data:
            output[group] = {
                "baseline": np.asarray(data["baseline"], dtype=np.float32),
                "candidate": np.asarray(data["candidate"], dtype=np.float32),
            }
    return output


def score_candidates(
    candidates: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    parent: ModuleType,
) -> dict[str, dict[str, np.ndarray]]:
    output = {}
    for dataset in SCORED_DATASETS:
        group = DATASET_GROUP[dataset]
        baseline = candidates[group]["baseline"]
        candidate = candidates[group]["candidate"]
        truth_values = truth[dataset]
        truth_unit = normalized_rows(truth_values)
        baseline_pds = pds_from_square_similarity(
            normalized_rows(baseline) @ truth_unit.T
        )
        candidate_pds = pds_from_square_similarity(
            normalized_rows(candidate) @ truth_unit.T
        )
        output[dataset] = {
            "pds": candidate_pds - baseline_pds,
            "signed_credit": (
                parent.signed_credit(candidate, truth_values)
                - parent.signed_credit(baseline, truth_values)
            ),
        }
    return output


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


def summarize(
    score: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, Any]:
    arrays = family_arrays(score)
    per_target = {
        metric: np.mean(np.stack(list(values.values())), axis=0)
        for metric, values in arrays.items()
    }
    return {
        "families": {
            metric: {
                family: float(values.mean()) for family, values in family_values.items()
            }
            for metric, family_values in arrays.items()
        },
        "macro": {
            metric: float(values.mean()) for metric, values in per_target.items()
        },
        "per_target": per_target,
    }


def block_count_matrix(draws: np.ndarray) -> np.ndarray:
    output = np.zeros((len(draws), BLOCK_COUNT), dtype=np.uint16)
    rows = np.repeat(np.arange(len(draws), dtype=np.int64), draws.shape[1])
    np.add.at(output, (rows, draws.ravel()), 1)
    return output


def bootstrap_summary(
    score: Mapping[str, Mapping[str, np.ndarray]],
    blocks: np.ndarray,
    draws: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    counts = block_count_matrix(np.asarray(draws, dtype=np.int64)).astype(np.float64)
    block_sizes = np.bincount(blocks, minlength=BLOCK_COUNT).astype(np.float64)
    denominator = counts @ block_sizes
    require(np.all(denominator > 0.0), "empty bootstrap draw")
    families = family_arrays(score)
    arrays: dict[str, np.ndarray] = {}
    for metric, family_values in families.items():
        values = []
        for family, per_target in family_values.items():
            sums = np.bincount(
                blocks,
                weights=per_target,
                minlength=BLOCK_COUNT,
            ).astype(np.float64)
            draws_for_family = (counts @ sums) / denominator
            arrays[f"{metric}__{family}"] = draws_for_family
            values.append(draws_for_family)
        arrays[metric] = np.mean(np.stack(values), axis=0)
    summary = {
        metric: {
            "lower_95": float(np.quantile(arrays[metric], 0.025)),
            "upper_95": float(np.quantile(arrays[metric], 0.975)),
            "mean": float(np.mean(arrays[metric])),
        }
        for metric in ("pds", "signed_credit")
    }
    return summary, arrays


def concentration_summary(
    score: Mapping[str, Mapping[str, np.ndarray]],
    folds: np.ndarray,
    blocks: np.ndarray,
) -> dict[str, dict[str, float]]:
    arrays = family_arrays(score)
    output = {}
    for metric, family_values in arrays.items():
        family_gain = np.asarray(
            [max(0.0, float(values.mean())) for values in family_values.values()]
        )
        per_target = np.mean(np.stack(list(family_values.values())), axis=0)
        block_gain = np.asarray(
            [
                max(0.0, float(per_target[blocks == block].sum()))
                for block in range(BLOCK_COUNT)
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

        output[metric] = {
            "family": fraction(family_gain),
            "target_block": fraction(block_gain),
            "outer_fold": fraction(fold_gain),
        }
    return output


def fold_summary(
    summary: Mapping[str, Any],
    folds: np.ndarray,
) -> dict[str, list[float]]:
    return {
        metric: [float(values[folds == fold].mean()) for fold in range(FOLD_COUNT)]
        for metric, values in summary["per_target"].items()
    }


def pre_null_checks(
    summary: Mapping[str, Any],
    bootstrap: Mapping[str, Any],
    folds: Mapping[str, Sequence[float]],
    concentration: Mapping[str, Mapping[str, float]],
) -> dict[str, bool]:
    pds_folds = np.asarray(folds["pds"], dtype=np.float64)
    return {
        "fixed_macro_pds_at_least_0.003": summary["macro"]["pds"] >= MIN_PDS_GAIN,
        "fixed_macro_signed_credit_nonnegative": summary["macro"]["signed_credit"]
        >= MIN_SIGNED_CREDIT_GAIN,
        "pds_bootstrap_lower_positive": bootstrap["pds"]["lower_95"] > 0.0,
        "signed_credit_bootstrap_lower_nonnegative": bootstrap["signed_credit"][
            "lower_95"
        ]
        >= 0.0,
        "every_family_pds_nonnegative": all(
            value >= 0.0 for value in summary["families"]["pds"].values()
        ),
        "every_family_signed_credit_nonnegative": all(
            value >= 0.0 for value in summary["families"]["signed_credit"].values()
        ),
        "at_least_four_outer_folds_positive_pds": int(np.sum(pds_folds > 0.0)) >= 4,
        "minimum_outer_fold_pds": float(pds_folds.min()) >= MIN_FOLD_PDS,
        "pds_family_concentration_at_most_0.6": concentration["pds"]["family"] <= 0.6,
        "pds_target_block_concentration_at_most_0.5": concentration["pds"][
            "target_block"
        ]
        <= 0.5,
        "pds_outer_fold_concentration_at_most_0.6": concentration["pds"]["outer_fold"]
        <= 0.6,
    }


def null_distribution(
    candidates: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    parent: ModuleType,
    permutations: np.ndarray,
    guard: Any,
) -> dict[str, np.ndarray]:
    identity_credit = {
        dataset: parent.signed_credit(
            candidates[DATASET_GROUP[dataset]]["baseline"],
            truth[dataset],
        )
        for dataset in SCORED_DATASETS
    }
    truth_unit = {
        dataset: normalized_rows(truth[dataset]) for dataset in SCORED_DATASETS
    }
    identity_pds = {
        dataset: pds_from_square_similarity(
            normalized_rows(candidates[DATASET_GROUP[dataset]]["baseline"])
            @ truth_unit[dataset].T
        )
        for dataset in SCORED_DATASETS
    }
    deltas = {
        group: candidates[group]["candidate"] - candidates[group]["baseline"]
        for group in GROUPS
    }
    output = {
        "pds": np.empty(len(permutations), dtype=np.float64),
        "signed_credit": np.empty(len(permutations), dtype=np.float64),
    }
    for draw, permutation in enumerate(permutations):
        null_candidates = {}
        for group in GROUPS:
            baseline = candidates[group]["baseline"]
            raw = baseline + deltas[group][permutation]
            null_candidates[group] = CORE.row_norm_preserving_blend(
                baseline,
                raw,
                1.0,
            )
        per_dataset = {}
        for dataset in SCORED_DATASETS:
            group = DATASET_GROUP[dataset]
            prediction = null_candidates[group]
            per_dataset[dataset] = {
                "pds": pds_from_square_similarity(
                    normalized_rows(prediction) @ truth_unit[dataset].T
                )
                - identity_pds[dataset],
                "signed_credit": (
                    parent.signed_credit(prediction, truth[dataset])
                    - identity_credit[dataset]
                ),
            }
        null_summary = summarize(per_dataset)
        for metric in output:
            output[metric][draw] = null_summary["macro"][metric]
        if draw % 5 == 0:
            guard.check(f"target-alignment null {draw}")
    return output


def evaluate(max_rss_gib: float, min_free_percent: float) -> None:
    protocol()
    require(PREOUTCOME_DIR.exists(), "development prefreeze is absent")
    require(not RESULT_DIR.exists(), f"development result exists: {RESULT_DIR}")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    effective_max = min(float(max_rss_gib), MAX_RSS_GIB)
    effective_floor = max(float(min_free_percent), MIN_FREE_PERCENT)
    external = CORE.validate_external_guard(effective_max, effective_floor)
    guard = CORE.MemoryGuard(effective_max, effective_floor)
    validate_preoutcome(guard)
    authorization = validate_development_authorization(guard)

    with CORE.exclusive_lock(DEVELOP_LOCK):
        create_outcome_access_lock(authorization, guard)
        truth, parent = load_truth(guard)
        candidates = load_candidates()
        with np.load(PREOUTCOME_DIR / PLAN_NAME, allow_pickle=False) as plan:
            folds = np.asarray(plan["development_folds"], dtype=np.int8)
            blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
            draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.int64)
            permutations = np.asarray(
                plan["target_alignment_permutations"],
                dtype=np.int64,
            )
        validate_plan(folds, blocks, draws)
        score = score_candidates(candidates, truth, parent)
        observed = summarize(score)
        fold_values = fold_summary(observed, folds)
        bootstrap, bootstrap_arrays = bootstrap_summary(score, blocks, draws)
        concentration = concentration_summary(score, folds, blocks)
        checks = pre_null_checks(observed, bootstrap, fold_values, concentration)
        pre_null_pass = all(checks.values())
        null_arrays: dict[str, np.ndarray] = {}
        null_summary: dict[str, Any] | None = None
        if pre_null_pass:
            null_arrays = null_distribution(
                candidates,
                truth,
                parent,
                permutations,
                guard,
            )
            null_summary = {
                metric: {
                    "p95": float(np.quantile(values, 0.95)),
                    "maximum": float(values.max()),
                    "mean": float(values.mean()),
                }
                for metric, values in null_arrays.items()
            }
            checks["observed_pds_above_null_p95"] = (
                observed["macro"]["pds"] > null_summary["pds"]["p95"]
            )
            checks["observed_signed_credit_above_null_p95"] = (
                observed["macro"]["signed_credit"]
                > null_summary["signed_credit"]["p95"]
            )
        else:
            checks["observed_pds_above_null_p95"] = False
            checks["observed_signed_credit_above_null_p95"] = False
        passes = all(checks.values())
        result = {
            "schema_version": 1,
            "stage": "authorized_613_target_consensus_covariance",
            "status": (
                "PASS_FOR_SEPARATELY_PREREGISTERED_CONFIRMATION"
                if passes
                else "TERMINAL_NO_GO"
            ),
            "candidate": CANDIDATE_NAME,
            "fixed": {
                "macro": observed["macro"],
                "families": observed["families"],
            },
            "outer_folds": fold_values,
            "bootstrap": bootstrap,
            "positive_gain_concentration": concentration,
            "pre_null_checks": {
                key: value
                for key, value in checks.items()
                if not key.startswith("observed_")
            },
            "pre_null_pass": pre_null_pass,
            "target_alignment_null_run": pre_null_pass,
            "target_alignment_null": null_summary,
            "checks": checks,
            "passes": passes,
            "development_613_outcomes_opened": True,
            "sealed_410_effect_values_opened": False,
            "reserved_confirmation_expression_opened": False,
            "truth_container_hashes_recomputed": False,
            "truth_access_mode": "selected 613 rows by 4992 response columns only",
            "count_generation": False,
            "cell_eval2": False,
            "challenge_build": False,
            "submission": False,
            "development_authorization_sha256": CORE.sha256_file(
                DEVELOPMENT_AUTHORIZATION,
                guard,
            ),
            "preoutcome_manifest_sha256": CORE.sha256_file(
                PREOUTCOME_DIR / MANIFEST_NAME,
                guard,
            ),
            "external_memory_guard": external,
            "memory": guard.report(),
        }
        arrays = {}
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                arrays[f"per_target__{dataset}__{metric}"] = score[dataset][metric]
        for name, values in bootstrap_arrays.items():
            arrays[f"bootstrap__{name}"] = values
        for name, values in null_arrays.items():
            arrays[f"null__{name}"] = values
        with CORE.atomic_stage_directory(RESULT_DIR, guard) as staging:
            CORE.atomic_npz(staging / ARRAYS_NAME, guard, **arrays)
            CORE.atomic_json(staging / RESULT_NAME, result, guard)
            CORE.atomic_json(
                staging / MANIFEST_NAME,
                CORE.directory_manifest(staging, guard),
                guard,
            )
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    protocol()
    guard = CORE.MemoryGuard(MAX_RSS_GIB, MIN_FREE_PERCENT)
    CORE.validate_directory_manifest(SOURCE_DIR, guard)
    payload: dict[str, Any] = {
        "schema_version": 1,
        "experiment": HERE.name,
        "prepare_authorization_exists": PREPARE_AUTHORIZATION.exists(),
        "development_prefreeze_exists": PREOUTCOME_DIR.exists(),
        "development_authorization_exists": DEVELOPMENT_AUTHORIZATION.exists(),
        "outcome_access_recorded": OUTCOME_ACCESS_LOCK.exists(),
        "development_result_exists": RESULT_DIR.exists(),
        "audit_opened_outcome_arrays": False,
    }
    if PREPARE_AUTHORIZATION.exists():
        validate_prepare_authorization(guard)
        payload["prepare_authorization_valid"] = True
    if PREOUTCOME_DIR.exists():
        validate_preoutcome(guard)
        payload["development_prefreeze_valid"] = True
    if DEVELOPMENT_AUTHORIZATION.exists() and PREOUTCOME_DIR.exists():
        validate_development_authorization(guard)
        payload["development_authorization_valid_without_truth_read"] = True
    if RESULT_DIR.exists():
        CORE.validate_directory_manifest(RESULT_DIR, guard)
        result = read_json(RESULT_DIR / RESULT_NAME)
        require(result["sealed_410_effect_values_opened"] is False, "410 opened")
        require(
            result["reserved_confirmation_expression_opened"] is False,
            "reserved opened",
        )
        require(result["count_generation"] is False, "counts generated")
        require(result["cell_eval2"] is False, "cell-eval2 run")
        require(result["challenge_build"] is False, "challenge build run")
        require(result["submission"] is False, "submission run")
        payload["development_result_valid"] = True
        payload["decision"] = result["status"]
    payload["memory"] = guard.report()
    print(json.dumps(payload, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=COMMANDS)
    parser.add_argument("--max-rss-gib", type=float, default=MAX_RSS_GIB)
    parser.add_argument("--min-free-percent", type=float, default=MIN_FREE_PERCENT)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "prepare":
        prepare(args.max_rss_gib, args.min_free_percent)
    elif args.command == "evaluate":
        evaluate(args.max_rss_gib, args.min_free_percent)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
