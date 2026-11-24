#!/usr/bin/env python3
"""Fail-closed runner for Experiment 143.

Commands
--------
prepare
    Source-only construction of the frozen normrel candidate and exact
    equal-weight baseline. This command cannot import or load development truth.
develop
    Authorized one-candidate evaluation on the 613-target development panel.
    It records the irreversible outcome-access lock before truth validation.
audit
    Read-only control-plane and artifact inspection. It never loads truth.

There is deliberately no confirmation, count, scorer, challenge-build, or
submission command.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from types import ModuleType
from typing import Any, Iterator, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = Path(__file__).resolve()
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
PREPARE_AUTHORIZATION = HERE / "SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json"
DEVELOPMENT_AUTHORIZATION = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"

EXP142 = ROOT / "experiments/142-batch-reliability-weighting"
EXP142_RUN = EXP142 / "develop.py"
EXP142_PREFREEZE_AUTHORIZATION = EXP142 / "DEVELOPMENT_PREFREEZE_AUTHORIZATION.v1.json"
EXP142_DEVELOPMENT_AUTHORIZATION = EXP142 / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
EXP142_PREOUTCOME = EXP142 / "work/development-preoutcome-v1"
EXP142_PREOUTCOME_MANIFEST = EXP142_PREOUTCOME / "artifact-manifest.json"
EXP142_TRUTH_LEDGER = EXP142_PREOUTCOME / "truth-source-hashes.json"
EXP142_PLAN = EXP142_PREOUTCOME / "development-plan.npz"
EXP142_SOURCE_BUILDS = {
    donor: EXP142 / "work/source-build-v4" / donor
    for donor in ("k562_gwps", "orion_hct116", "orion_hek293t")
}

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
EXP132_RUN = EXP132 / "run.py"
EXP132_AUTHORIZATION = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
EXP132_PLAN = EXP132 / "work/preregistered-freeze-v2"
DEVELOPMENT_TARGETS = EXP132_PLAN / "selection_targets.txt"
RESPONSE_GENES = EXP132_PLAN / "response_genes.txt"

PREOUTCOME_DIR = HERE / "work/preoutcome-v1"
RESULT_DIR = HERE / "work/development-result-v1"
PREPARE_LOCK = HERE / "work/preoutcome-v1.lock"
DEVELOP_LOCK = HERE / "work/development-result-v1.lock"
OUTCOME_ACCESS_LOCK = HERE / "work/development-outcome-access.v1.lock"

PLAN_ARTIFACT = "development-plan.npz"
EFFECT_ARTIFACT = {
    "k562": "effects-k562.npz",
    "non_k562": "effects-non_k562.npz",
}
TRUTH_LEDGER_ARTIFACT = "truth-source-hashes.json"
RESULT_ARTIFACT = "result.json"
ARRAY_ARTIFACT = "per-target-bootstrap.npz"
MANIFEST_ARTIFACT = "artifact-manifest.json"

GROUP_SLOTS = {
    "k562": ("orion_hct116", "orion_hek293t"),
    "non_k562": ("k562_composite", "orion_hct116", "orion_hek293t"),
}
SLOT_SOURCE = {
    "k562_composite": "k562_gwps",
    "orion_hct116": "orion_hct116",
    "orion_hek293t": "orion_hek293t",
}
SOURCE_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
SCORED_DATASETS = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat")
FAMILIES = ("k562", "rpe1", "jurkat")

DEV_COUNT = 613
RESPONSE_COUNT = 4992
FOLD_COUNT = 5
EXPECTED_FOLD_SIZES = (123, 123, 123, 122, 122)
BOOTSTRAP_DRAWS = 20_000
EXPONENT = -1.5
FLOOR_PERCENTILE = 10.0
EPS = 1e-12
MIN_GAIN = 0.003
MAX_CONCENTRATION = 0.5


class ContractError(RuntimeError):
    """Raised when a frozen experimental boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    require(bool(values), f"empty axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while piece := handle.read(block):
            digest.update(piece)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    return hashlib.sha256(value.tobytes()).hexdigest()


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


def protocol() -> dict[str, Any]:
    value = read_json(PROTOCOL_PATH)
    require(value.get("schema_version") == 1, "protocol schema changed")
    require(value.get("experiment") == HERE.name, "protocol experiment changed")
    candidate = value["candidate"]
    require(candidate["count"] == 1, "candidate count changed")
    formula = candidate["per_group_formula"]
    require("^-1.5" in formula["raw_weight"], "candidate exponent changed")
    require("10th percentile" in formula["floor"], "candidate floor percentile changed")
    require(
        formula["denominator"]
        == "median target norm over the same 613 rows and 4,992 response genes",
        "candidate denominator changed",
    )
    return value


def resolve_repo_path(value: str) -> Path:
    path = (ROOT / value).resolve()
    require(ROOT == path or ROOT in path.parents, f"path escapes repository: {value}")
    return path


def frozen_bindings() -> dict[Path, str]:
    bindings = protocol()["frozen_bindings"]
    return {resolve_repo_path(name): str(digest) for name, digest in bindings.items()}


def validate_hashes(bindings: Mapping[Path, str]) -> None:
    for path, expected in bindings.items():
        require(path.is_file(), f"bound source missing: {path}")
        require(sha256_file(path) == expected, f"bound source changed: {path}")


def authorization_bindings(authorization: Mapping[str, Any]) -> dict[Path, str]:
    raw = authorization.get("bindings")
    require(isinstance(raw, Mapping), "authorization lacks bindings")
    result: dict[Path, str] = {}
    for name, value in raw.items():
        path = Path(str(name))
        path = path.resolve() if path.is_absolute() else resolve_repo_path(str(name))
        digest = value if isinstance(value, str) else value.get("sha256")
        require(isinstance(digest, str), f"invalid authorization binding: {name}")
        require(path not in result, f"duplicate authorization binding: {name}")
        result[path] = digest
    return result


def validate_bound_paths(
    authorization: Mapping[str, Any],
    required_paths: Sequence[Path],
    *,
    hash_files: bool,
) -> dict[Path, str]:
    bindings = authorization_bindings(authorization)
    for path in required_paths:
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
    ):
        require(authorization.get(key) is False, f"unsafe authorization flag: {key}")


def validate_prepare_authorization() -> dict[str, Any]:
    require(
        PREPARE_AUTHORIZATION.is_file(),
        "SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json is required; do not self-authorize",
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
    required = (RUN_PATH, PROTOCOL_PATH, *frozen_bindings())
    bindings = validate_bound_paths(authorization, required, hash_files=True)
    for path, expected in frozen_bindings().items():
        require(
            bindings[path.resolve()] == expected,
            f"prepare authorization differs from protocol: {path}",
        )
    return authorization


def validate_142_lineage() -> None:
    parent = load_experiment142()
    parent.validate_development_authorization(verify_truth_files=False)
    preoutcome = read_json(EXP142_PREOUTCOME / RESULT_ARTIFACT)
    require(
        preoutcome["development_outcomes_opened"] is False, "142 prepare opened truth"
    )
    require(
        preoutcome["sealed_410_effect_values_opened"] is False,
        "142 prepare opened sealed outcomes",
    )


def validate_source_inputs() -> None:
    validate_hashes(frozen_bindings())
    validate_142_lineage()
    for directory in EXP142_SOURCE_BUILDS.values():
        validate_directory_manifest(directory)


def load_experiment142() -> ModuleType:
    expected = frozen_bindings()[EXP142_RUN.resolve()]
    require(sha256_file(EXP142_RUN) == expected, "Experiment 142 evaluator changed")
    spec = importlib.util.spec_from_file_location("experiment142_for_143", EXP142_RUN)
    require(
        spec is not None and spec.loader is not None, "cannot import Experiment 142"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def percentile_floor(
    keys: np.ndarray,
    available: np.ndarray,
    *,
    percentile: float = FLOOR_PERCENTILE,
) -> float:
    values = np.asarray(keys, dtype=np.float64)
    mask = np.asarray(available, dtype=bool)
    require(values.shape == mask.shape, "key/availability shape mismatch")
    require(np.all(np.isfinite(values[mask])), "non-finite available normrel key")
    require(np.all(values[mask] >= 0.0), "negative normrel key")
    positive = values[mask & (values > 0.0)]
    require(positive.size > 0, "group has no positive normrel keys")
    floor = float(np.percentile(positive, percentile, method="linear"))
    require(np.isfinite(floor) and floor > 0.0, "invalid normrel floor")
    return floor


def normrel_weights(
    keys: np.ndarray,
    available: np.ndarray,
    *,
    exponent: float = EXPONENT,
    percentile: float = FLOOR_PERCENTILE,
) -> tuple[np.ndarray, float]:
    values = np.asarray(keys, dtype=np.float64)
    mask = np.asarray(available, dtype=bool)
    require(values.ndim == 2, "normrel keys must be target by donor")
    require(values.shape == mask.shape, "key/availability shape mismatch")
    require(np.all(mask.sum(axis=1) >= 1), "target without an available donor")
    floor = percentile_floor(values, mask, percentile=percentile)
    raw = np.zeros_like(values)
    raw[mask] = np.maximum(values[mask], floor) ** exponent
    totals = raw.sum(axis=1, keepdims=True)
    require(np.all(totals > 0.0), "zero normrel weight total")
    weights = raw / totals
    require(
        np.allclose(weights.sum(axis=1), 1.0, rtol=0.0, atol=EPS),
        "normrel weights are not normalized",
    )
    require(np.all(weights[~mask] == 0.0), "unavailable donor received weight")
    return weights, floor


def equal_weights(available: np.ndarray) -> np.ndarray:
    mask = np.asarray(available, dtype=bool)
    require(mask.ndim == 2, "availability must be target by donor")
    counts = mask.sum(axis=1, keepdims=True)
    require(np.all(counts >= 1), "target without an available donor")
    return mask.astype(np.float64) / counts


def mix_slots_exact(
    slots: np.ndarray,
    weights: np.ndarray,
    available: np.ndarray,
    exact_equal_baseline: np.ndarray | None = None,
) -> np.ndarray:
    values = np.asarray(slots)
    active = np.asarray(available, dtype=bool)
    w = np.asarray(weights, dtype=np.float64)
    require(values.ndim == 3, "slots must be target by donor by gene")
    require(w.shape == values.shape[:2], "slot weight shape mismatch")
    require(active.shape == w.shape, "slot availability shape mismatch")
    equal = equal_weights(active)
    if exact_equal_baseline is not None and np.array_equal(w, equal):
        baseline = np.asarray(exact_equal_baseline)
        require(
            baseline.shape == (values.shape[0], values.shape[2]),
            "exact baseline shape mismatch",
        )
        return baseline.copy()
    return np.einsum("td,tdg->tg", w, values, optimize=True)


def group_normrel_keys(
    group: str,
    centered: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    slots = GROUP_SLOTS[group]
    keys = np.empty((DEV_COUNT, len(slots)), dtype=np.float64)
    available = np.ones_like(keys, dtype=bool)
    denominators = {}
    for column, slot in enumerate(slots):
        source = SLOT_SOURCE[slot]
        values = np.asarray(centered[source], dtype=np.float64)
        require(
            values.shape == (DEV_COUNT, RESPONSE_COUNT),
            f"{source}: centered response shape changed",
        )
        denominator = float(np.median(np.linalg.norm(values, axis=1)))
        require(denominator > 0.0, f"{source}: invalid development-panel median norm")
        denominators[slot] = denominator
        keys[:, column] = np.linalg.norm(values, axis=1) / denominator
    return keys, available, denominators


def load_source_centered() -> dict[str, np.ndarray]:
    output = {}
    for donor, directory in EXP142_SOURCE_BUILDS.items():
        values = np.load(directory / "pooled-centered.npy", mmap_mode="r")
        require(
            values.shape == (DEV_COUNT, RESPONSE_COUNT),
            f"{donor}: pooled-centered shape changed",
        )
        output[donor] = np.asarray(values, dtype=np.float64)
    return output


def prepare() -> None:
    validate_prepare_authorization()
    validate_source_inputs()
    require(not PREOUTCOME_DIR.exists(), f"preoutcome output exists: {PREOUTCOME_DIR}")
    require(not RESULT_DIR.exists(), "development result already exists")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    parent = load_experiment142()
    targets = read_lines(DEVELOPMENT_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    require(len(targets) == DEV_COUNT, "development target count changed")
    require(len(responses) == RESPONSE_COUNT, "response gene count changed")

    with (
        exclusive_lock(PREPARE_LOCK),
        atomic_stage_directory(PREOUTCOME_DIR) as staging,
    ):
        with np.load(EXP142_PLAN, allow_pickle=False) as plan:
            folds = np.asarray(plan["development_folds"], dtype=np.int8)
            blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
            draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.uint16)
        require(
            tuple(int(np.sum(folds == fold)) for fold in range(FOLD_COUNT))
            == EXPECTED_FOLD_SIZES,
            "frozen fold sizes changed",
        )
        require(draws.shape[0] == BOOTSTRAP_DRAWS, "bootstrap draw count changed")
        atomic_npz(
            staging / PLAN_ARTIFACT,
            development_folds=folds,
            development_blocks=blocks,
            bootstrap_block_draws=draws,
        )

        centered = load_source_centered()

        group_reports = {}
        for group in GROUP_SLOTS:
            geometry = parent.load_prepared_geometry(group)
            require(
                tuple(map(str, geometry["donors"])) == GROUP_SLOTS[group],
                f"{group}: frozen slot order changed",
            )
            keys, available, denominators = group_normrel_keys(group, centered)
            candidate_weights, floor = normrel_weights(keys, available)
            baseline_weights = equal_weights(available)
            baseline_final = mix_slots_exact(
                geometry["slots"],
                baseline_weights,
                available,
                exact_equal_baseline=geometry["baseline"],
            )
            require(
                np.array_equal(baseline_final, geometry["baseline"]),
                f"{group}: exact identity baseline failed",
            )
            candidate_final = mix_slots_exact(
                geometry["slots"], candidate_weights, available
            )
            atomic_npz(
                staging / EFFECT_ARTIFACT[group],
                donors=np.asarray(GROUP_SLOTS[group], dtype="U32"),
                available=available,
                normrel_keys=keys,
                baseline_weights=baseline_weights,
                candidate_weights=candidate_weights,
                baseline=baseline_final,
                candidate=candidate_final,
            )
            group_reports[group] = {
                "donors": list(GROUP_SLOTS[group]),
                "normrel_floor": floor,
                "normrel_key_quantiles": {
                    str(q): float(np.quantile(keys[available], q))
                    for q in (0.0, 0.1, 0.5, 0.9, 1.0)
                },
                "candidate_weight_quantiles": {
                    str(q): float(np.quantile(candidate_weights[available], q))
                    for q in (0.0, 0.1, 0.5, 0.9, 1.0)
                },
                "one_donor_targets": int(np.sum(available.sum(axis=1) == 1)),
                "k562_norm_source": (
                    "original_k562_gwps" if group == "non_k562" else None
                ),
                "same_axis_median_norms": denominators,
                "perturbmatch": (
                    "already embedded in frozen Experiment 142 k562_composite slot"
                    if group == "non_k562"
                    else None
                ),
                "positive_scalar_steps_omitted": (
                    "alpha, group-global scale, and whitening preserve do not alter "
                    "PDS or signed credit; CD4 has zero complete support on this panel"
                ),
                "exact_identity_baseline": True,
                "baseline_sha256": array_sha256(baseline_final),
                "candidate_sha256": array_sha256(candidate_final),
            }

        atomic_copy(EXP142_TRUTH_LEDGER, staging / TRUTH_LEDGER_ARTIFACT)
        result = {
            "schema_version": 1,
            "stage": "source_only_preoutcome_prepare",
            "development_outcomes_opened": False,
            "sealed_410_effect_values_opened": False,
            "reserved_confirmation_expression_opened": False,
            "truth_arrays_opened": False,
            "candidate_count": 1,
            "candidate": "normrel613^-1.5",
            "fold_sizes": [int(np.sum(folds == fold)) for fold in range(FOLD_COUNT)],
            "bootstrap_draws": int(draws.shape[0]),
            "prepare_authorization_sha256": sha256_file(PREPARE_AUTHORIZATION),
            "run_sha256": sha256_file(RUN_PATH),
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "upstream_142_preoutcome_manifest_sha256": sha256_file(
                EXP142_PREOUTCOME_MANIFEST
            ),
            "upstream_142_truth_ledger_sha256": sha256_file(EXP142_TRUTH_LEDGER),
            "cd4": {
                "active": False,
                "complete_targets": 0,
                "reason": "no complete mean-three CD4 support on the 613-target panel",
            },
            "groups": group_reports,
        }
        atomic_json(staging / RESULT_ARTIFACT, result)
        atomic_json(staging / MANIFEST_ARTIFACT, directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_preoutcome_control_plane() -> dict[str, Any]:
    validate_directory_manifest(PREOUTCOME_DIR)
    result = read_json(PREOUTCOME_DIR / RESULT_ARTIFACT)
    require(result["development_outcomes_opened"] is False, "prepare opened truth")
    require(
        result["sealed_410_effect_values_opened"] is False,
        "prepare opened sealed outcomes",
    )
    require(result["truth_arrays_opened"] is False, "prepare opened truth arrays")
    require(result["candidate_count"] == 1, "candidate count changed")
    require(result["candidate"] == "normrel613^-1.5", "candidate changed")
    require(
        result["run_sha256"] == sha256_file(RUN_PATH), "run.py changed after prepare"
    )
    require(
        result["protocol_sha256"] == sha256_file(PROTOCOL_PATH),
        "protocol changed after prepare",
    )
    require(
        result["prepare_authorization_sha256"] == sha256_file(PREPARE_AUTHORIZATION),
        "prepare authorization changed",
    )
    require(
        result["upstream_142_preoutcome_manifest_sha256"]
        == sha256_file(EXP142_PREOUTCOME_MANIFEST),
        "Experiment 142 preoutcome manifest changed",
    )
    require(
        result["upstream_142_truth_ledger_sha256"] == sha256_file(EXP142_TRUTH_LEDGER),
        "Experiment 142 truth ledger changed",
    )
    require(
        sha256_file(PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT)
        == sha256_file(EXP142_TRUTH_LEDGER),
        "copied truth ledger differs from Experiment 142",
    )
    return result


def truth_ledger_paths() -> tuple[Path, ...]:
    ledger = read_json(PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT)
    paths = []
    for records in ledger["datasets"].values():
        for record in records.values():
            paths.append(resolve_repo_path(record["path"]))
    return tuple(paths)


def validate_truth_ledger_files(bindings: Mapping[Path, str]) -> None:
    ledger = read_json(PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT)
    for records in ledger["datasets"].values():
        for record in records.values():
            path = resolve_repo_path(record["path"])
            require(path.resolve() in bindings, f"truth source not authorized: {path}")
            require(
                bindings[path.resolve()] == record["sha256"],
                f"authorization differs from truth ledger: {path}",
            )
            require(
                path.stat().st_size == record["bytes"], f"truth bytes changed: {path}"
            )
            require(
                sha256_file(path) == record["sha256"], f"truth hash changed: {path}"
            )


def development_control_paths() -> tuple[Path, ...]:
    return (
        RUN_PATH,
        PROTOCOL_PATH,
        PREPARE_AUTHORIZATION,
        PREOUTCOME_DIR / MANIFEST_ARTIFACT,
        PREOUTCOME_DIR / TRUTH_LEDGER_ARTIFACT,
        EXP142_RUN,
        EXP142_DEVELOPMENT_AUTHORIZATION,
        EXP132_RUN,
        EXP132_AUTHORIZATION,
        EXP132_PLAN / "plan-manifest.json",
        DEVELOPMENT_TARGETS,
        RESPONSE_GENES,
    )


def validate_development_authorization(*, verify_truth_files: bool) -> dict[str, Any]:
    validate_preoutcome_control_plane()
    require(
        DEVELOPMENT_AUTHORIZATION.is_file(),
        "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json is required; do not self-authorize",
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
    for path in truth_ledger_paths():
        require(path.resolve() in bindings, f"development authorization misses {path}")
    if verify_truth_files:
        validate_truth_ledger_files(bindings)
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
            "authorization_sha256": sha256_file(DEVELOPMENT_AUTHORIZATION),
            "authorization_schema": authorization.get("schema_version"),
        },
    )


def load_prepared_plan() -> dict[str, np.ndarray]:
    with np.load(PREOUTCOME_DIR / PLAN_ARTIFACT, allow_pickle=False) as data:
        return {name: np.asarray(data[name]) for name in data.files}


def load_prepared_effects(group: str) -> dict[str, np.ndarray]:
    with np.load(PREOUTCOME_DIR / EFFECT_ARTIFACT[group], allow_pickle=False) as data:
        return {name: np.asarray(data[name]) for name in data.files}


def score_candidate(
    effects: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    parent: ModuleType,
) -> dict[str, dict[str, np.ndarray]]:
    indices = np.arange(DEV_COUNT, dtype=np.int64)
    output = {}
    for dataset in SCORED_DATASETS:
        group = parent.DATASET_GROUP[dataset]
        baseline = effects[group]["baseline"]
        candidate = effects[group]["candidate"]
        output[dataset] = {
            "pds": (
                parent.pds_subset(candidate, truth[dataset], indices, indices)
                - parent.pds_subset(baseline, truth[dataset], indices, indices)
            ),
            "signed_credit": (
                parent.signed_credit(candidate, truth[dataset])
                - parent.signed_credit(baseline, truth[dataset])
            ),
        }
    return output


def promotion_checks(
    summary: Mapping[str, Any],
    fold_cells: Mapping[str, Mapping[str, Sequence[float]]],
    bootstrap: Mapping[str, Mapping[str, float]],
    concentration: Mapping[str, Mapping[str, float]],
) -> dict[str, bool]:
    checks = {}
    for metric in ("pds", "signed_credit"):
        checks[f"{metric}_macro_at_least_0.003"] = summary["macro"][metric] >= MIN_GAIN
        checks[f"{metric}_bootstrap_lower_positive"] = (
            bootstrap[metric]["lower_95"] > 0.0
        )
        checks[f"{metric}_every_family_nonnegative"] = all(
            value >= -EPS for value in summary["families"][metric].values()
        )
        checks[f"{metric}_every_family_fold_cell_nonnegative"] = all(
            value >= -EPS for family in FAMILIES for value in fold_cells[metric][family]
        )
        checks[f"{metric}_family_concentration_at_most_0.5"] = (
            concentration[metric]["family"] <= MAX_CONCENTRATION
        )
        checks[f"{metric}_target_block_concentration_at_most_0.5"] = (
            concentration[metric]["target_block"] <= MAX_CONCENTRATION
        )
    return checks


def develop() -> None:
    authorization = validate_development_authorization(verify_truth_files=False)
    require(not RESULT_DIR.exists(), f"development result exists: {RESULT_DIR}")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    with exclusive_lock(DEVELOP_LOCK):
        create_outcome_access_lock(authorization)
        validate_development_authorization(verify_truth_files=True)
        parent = load_experiment142()
        truth = parent.load_development_truth(0.0)
        plan = load_prepared_plan()
        folds = np.asarray(plan["development_folds"], dtype=np.int8)
        blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
        draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.int64)
        effects = {group: load_prepared_effects(group) for group in GROUP_SLOTS}
        score = score_candidate(effects, truth, parent)
        indices = np.arange(DEV_COUNT, dtype=np.int64)
        summary = parent.aggregate_scores(score, indices)
        fold_cells = parent.fold_family_cells(score, folds)
        bootstrap, bootstrap_arrays = parent.bootstrap_summary(score, blocks, draws)
        concentration = parent.concentration_summary(score, blocks)
        checks = promotion_checks(summary, fold_cells, bootstrap, concentration)
        passes = all(checks.values())
        result = {
            "schema_version": 1,
            "stage": "authorized_613_target_one_candidate_development",
            "candidate": "normrel613^-1.5",
            "gallery": "full_613_for_every_query",
            "macro": summary["macro"],
            "families": summary["families"],
            "fold_family_cells": fold_cells,
            "bootstrap": bootstrap,
            "positive_gain_concentration": concentration,
            "checks": checks,
            "passes": passes,
            "decision": "PASS_FOR_SEPARATE_REVIEW_ONLY" if passes else "TERMINAL_NO_GO",
            "null_stage_exists": False,
            "development_613_outcomes_opened": True,
            "sealed_410_effect_values_opened": False,
            "reserved_confirmation_expression_opened": False,
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
                    arrays[f"delta__{dataset}__{metric}"] = score[dataset][metric]
            for metric, values in bootstrap_arrays.items():
                arrays[f"bootstrap__{metric}"] = values
            atomic_npz(staging / ARRAY_ARTIFACT, **arrays)
            atomic_json(staging / RESULT_ARTIFACT, result)
            atomic_json(staging / MANIFEST_ARTIFACT, directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    value = protocol()
    validate_source_inputs()
    payload: dict[str, Any] = {
        "schema_version": 1,
        "experiment": value["experiment"],
        "candidate_count": value["candidate"]["count"],
        "prepare_authorization_exists": PREPARE_AUTHORIZATION.exists(),
        "preoutcome_exists": PREOUTCOME_DIR.exists(),
        "development_authorization_exists": DEVELOPMENT_AUTHORIZATION.exists(),
        "outcome_access_recorded": OUTCOME_ACCESS_LOCK.exists(),
        "development_result_exists": RESULT_DIR.exists(),
        "truth_arrays_opened_by_audit": False,
        "sealed_test_loader_exists": False,
    }
    if PREOUTCOME_DIR.exists():
        validate_preoutcome_control_plane()
        payload["preoutcome_valid_without_truth_read"] = True
    if DEVELOPMENT_AUTHORIZATION.exists() and PREOUTCOME_DIR.exists():
        validate_development_authorization(verify_truth_files=False)
        payload["development_authorization_valid_without_truth_read"] = True
    if RESULT_DIR.exists():
        validate_directory_manifest(RESULT_DIR)
        result = read_json(RESULT_DIR / RESULT_ARTIFACT)
        require(
            result["sealed_410_effect_values_opened"] is False,
            "result reports sealed-outcome access",
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
