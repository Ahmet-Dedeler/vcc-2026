#!/usr/bin/env python3
"""Development evaluator for Experiment 142.

Commands
--------
prepare
    Source-only, pre-outcome materialization. Requires the separately supplied
    DEVELOPMENT_PREFREEZE_AUTHORIZATION.v1.json and cannot load the four
    Experiment 132 development truth matrices.
develop
    Requires a separately supplied DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json,
    validates every frozen binding and truth-source hash, records an atomic
    outcome-access lock, then evaluates only the 613 development targets.
audit
    Read-only validation of existing control-plane and artifact state. It never
    imports Experiment 132's outcome loader or reads a truth matrix.

There is intentionally no sealed-test, count-generation, cell-eval, challenge
build, or submission command.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import resource
import shutil
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

DEVELOPMENT_INTERFACE_V1 = HERE / "DEVELOPMENT_INTERFACE.v1.json"
DEVELOPMENT_INTERFACE_V2 = HERE / "DEVELOPMENT_INTERFACE.v2.json"
DEVELOPMENT_INTERFACE_V3 = HERE / "DEVELOPMENT_INTERFACE.v3.json"
DEVELOPMENT_AMENDMENT = HERE / "DEVELOPMENT_PREFREEZE_AMENDMENT_001.json"
DEVELOPMENT_AMENDMENT_002 = HERE / "DEVELOPMENT_PREFREEZE_AMENDMENT_002.json"
DEVELOPMENT_PREFREEZE_AUTHORIZATION = (
    HERE / "DEVELOPMENT_PREFREEZE_AUTHORIZATION.v1.json"
)
DEVELOPMENT_AUTHORIZATION = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
SOURCE_AUTHORIZATION = HERE / "SOURCE_ONLY_FREEZE_AUTHORIZATION.v6.json"
SOURCE_PREOUTCOME_AUTHORIZATION = (
    HERE / "work/preoutcome-freeze-v4/PREOUTCOME_AUTHORIZATION.v1.json"
)
SOURCE_GATE_DIR = HERE / "work/source-gate-v4"
SOURCE_BUILD_DIR = HERE / "work/source-build-v4"
SOURCE_EVALUATION_DIR = HERE / "work/source-evaluation-v4"
SOURCE_PLAN_DIR = HERE / "work/preregistered-plan-v2"
SOURCE_PLAN_MANIFEST = SOURCE_PLAN_DIR / "plan-manifest.json"
SOURCE_PLAN_SCHEMA = SOURCE_PLAN_DIR / "source-schema.json"

EXP140 = ROOT / "experiments/140-batch-matched-controls"
EXP140_PLAN_ARRAYS = EXP140 / "work/preregistered-plan-v2/plan-arrays.npz"
EXP140_IDENTITY_DIR = EXP140 / "work/identity-reproduction-v5"
EXP140_SOURCE_MANIFEST = EXP140 / "SOURCE_HASH_MANIFEST.v5.json"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
EXP132_RUN = EXP132 / "run.py"
EXP132_AUTHORIZATION = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
EXP132_PLAN = EXP132 / "work/preregistered-freeze-v2"
DEVELOPMENT_TARGETS = EXP132_PLAN / "selection_targets.txt"
RESPONSE_GENES = EXP132_PLAN / "response_genes.txt"

EXP061_RUN = ROOT / "experiments/061-perturbmatch-hct-count/build_candidate.py"
PERTURBMATCH_ROOT = ROOT / "data/external/perturbmatch"
CANONICAL_GENES = ROOT / "data/challenge/gene_names.csv"

PREOUTCOME_DIR = HERE / "work/development-preoutcome-v1"
RESULT_DIR = HERE / "work/development-result-v1"
PREPARE_LOCK = HERE / "work/development-preoutcome-v1.lock"
DEVELOP_LOCK = HERE / "work/development-result-v1.lock"
OUTCOME_ACCESS_LOCK = HERE / "work/development-outcome-access.v1.lock"

PLAN_ARTIFACT = "development-plan.npz"
GEOMETRY_ARTIFACT = {
    "k562": "geometry-k562.npz",
    "non_k562": "geometry-non_k562.npz",
}
TARGET_NULL_ARTIFACT = {
    "k562": "target-pairing-null-weights-k562.npz",
    "non_k562": "target-pairing-null-weights-non_k562.npz",
}
DONOR_NULL_ARTIFACT = {
    "k562": "donor-assignment-null-weights-k562.npz",
    "non_k562": "donor-assignment-null-weights-non_k562.npz",
}
TRUTH_HASH_ARTIFACT = "truth-source-hashes.json"
RESULT_ARTIFACT = "result.json"
ARRAY_ARTIFACT = "per-target-nulls-bootstrap.npz"
MANIFEST_ARTIFACT = "artifact-manifest.json"

SOURCE_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
GROUP_SLOTS = {
    "k562": ("orion_hct116", "orion_hek293t"),
    "non_k562": ("k562_composite", "orion_hct116", "orion_hek293t"),
}
SLOT_SOURCE = {
    "k562_composite": "k562_gwps",
    "orion_hct116": "orion_hct116",
    "orion_hek293t": "orion_hek293t",
}
SCORED_DATASETS = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat")
DATASET_GROUP = {
    "k562_gwps": "k562",
    "k562_essential": "k562",
    "rpe1": "non_k562",
    "nadig_jurkat": "non_k562",
}
DATASET_FAMILY = {
    "k562_gwps": "k562",
    "k562_essential": "k562",
    "rpe1": "rpe1",
    "nadig_jurkat": "jurkat",
}
FAMILIES = ("k562", "rpe1", "jurkat")
SIGNALS = {
    "reliability": "reliability_z",
    "cell_count": "cell_count_z",
    "total_umi": "total_umi_z",
    "occupied_batch_count": "occupied_batch_count_z",
    "pooled_effect_norm": "pooled_norm_z",
}
ANCESTRY_CONTROLS = tuple(name for name in SIGNALS if name != "reliability")
TRUTH_ROOTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "k562_essential": ROOT / "data/donors/k562_essential",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
}
SOURCE_DONOR_ROOTS = {name: ROOT / "data/donors" / name for name in SOURCE_DONORS}

ETAS = (0.0, 0.25, 0.5)
EXPECTED_FOLD_SIZES = (123, 123, 123, 122, 122)
FOLDS = 5
DEV_COUNT = 613
RESPONSE_COUNT = 4992
NULL_DRAWS = 49
BOOTSTRAP_DRAWS = 20_000
WHITEN_GAMMA = 0.65
EPS = 1e-12
MIN_GAIN = 0.003
MIN_CONTROL_PDS_MARGIN = 0.001
MAX_CONCENTRATION = 0.5
MIN_DIRECTIONAL_VALUES = 12
PM_BLEND = 0.5
EXPECTED_PM_TARGET_SUPPORT = 605
EXPECTED_PM_RESPONSE_SUPPORT = 4750
GIB = 1024**3
USER_PROVENANCE_ANCHOR = "2026-08-31T18:45:00+00:00"


class ContractError(RuntimeError):
    """Raised when a frozen experimental boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    require(bool(values), f"empty axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def read_csv_first_column(path: Path) -> tuple[str, ...]:
    lines = read_lines(path, unique=False)
    require(len(lines) > 1, f"empty CSV: {path}")
    values = tuple(line.split(",", 1)[0].strip() for line in lines[1:])
    require(len(values) == len(set(values)), f"duplicate CSV axis: {path}")
    return values


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while piece := handle.read(block):
            digest.update(piece)
    return digest.hexdigest()


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def memory_check(label: str, limit: float) -> None:
    value = peak_rss_gib()
    print(f"[rss] {label}: {value:.3f} GiB", flush=True)
    if limit > 0 and value >= limit:
        raise MemoryError(f"{label}: peak RSS {value:.3f} GiB >= {limit:.3f} GiB")


def file_record(path: Path) -> dict[str, Any]:
    return {"bytes": path.stat().st_size, "sha256": sha256_file(path)}


def directory_manifest(directory: Path) -> dict[str, Any]:
    files = {}
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.name != MANIFEST_ARTIFACT:
            files[path.name] = file_record(path)
    return {"schema_version": 1, "files": files}


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


def binding_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        resolved = path.resolve()
    elif value.startswith(("experiments/", "data/", "scripts/", "src/")):
        resolved = (ROOT / path).resolve()
    else:
        resolved = (HERE / path).resolve()
    require(
        resolved == ROOT or ROOT in resolved.parents,
        f"binding escapes repository: {value}",
    )
    return resolved


def binding_digest(value: Any) -> str:
    if isinstance(value, str):
        return value
    require(
        isinstance(value, Mapping), "authorization binding must be a hash or object"
    )
    digest = value.get("sha256")
    require(isinstance(digest, str), "authorization binding lacks sha256")
    return digest


def normalized_authorization_bindings(
    authorization: Mapping[str, Any],
) -> dict[Path, str]:
    raw = authorization.get("bindings")
    require(isinstance(raw, Mapping), "authorization lacks bindings")
    result: dict[Path, str] = {}
    for name, value in raw.items():
        if str(name).endswith("_provenance_utc"):
            continue
        path = binding_path(str(name))
        require(path not in result, f"duplicate resolved binding: {name}")
        result[path] = binding_digest(value)
    return result


def authorized_provenance_utc(authorization: Mapping[str, Any], stage: str) -> str:
    key = f"{stage}_provenance_utc"
    bindings = authorization.get("bindings")
    require(isinstance(bindings, Mapping), "authorization lacks bindings")
    expected_text = bindings.get(key)
    require(isinstance(expected_text, str), f"authorization lacks {key}")
    actual_text = os.environ.get("VCC_PROVENANCE_UTC")
    require(
        actual_text is not None,
        "artifact stage requires explicit VCC_PROVENANCE_UTC",
    )
    try:
        expected = dt.datetime.fromisoformat(expected_text.replace("Z", "+00:00"))
        actual = dt.datetime.fromisoformat(actual_text.replace("Z", "+00:00"))
        anchor = dt.datetime.fromisoformat(
            USER_PROVENANCE_ANCHOR.replace("Z", "+00:00")
        )
    except ValueError as error:
        raise ContractError("invalid provenance timestamp") from error
    require(expected.tzinfo is not None, f"{key} must include a timezone")
    require(actual.tzinfo is not None, "VCC_PROVENANCE_UTC must include a timezone")
    expected_utc = expected.astimezone(dt.UTC)
    actual_utc = actual.astimezone(dt.UTC)
    anchor_utc = anchor.astimezone(dt.UTC)
    require(
        expected_utc == anchor_utc,
        f"{key} does not match the user-provided provenance anchor",
    )
    require(
        actual_utc == expected_utc,
        "VCC_PROVENANCE_UTC does not match the authorization binding",
    )
    return actual_utc.isoformat()


def validate_bound_paths(
    authorization: Mapping[str, Any], required: Sequence[Path]
) -> dict[Path, str]:
    bindings = normalized_authorization_bindings(authorization)
    for path in required:
        resolved = path.resolve()
        require(resolved in bindings, f"authorization missing binding: {path}")
        require(resolved.is_file(), f"bound file missing: {path}")
        require(
            sha256_file(resolved) == bindings[resolved],
            f"authorization hash mismatch: {path}",
        )
    return bindings


def required_prefreeze_paths() -> tuple[Path, ...]:
    interface_v2 = read_json(DEVELOPMENT_INTERFACE_V2)
    interface_v3 = read_json(DEVELOPMENT_INTERFACE_V3)
    names = (
        *interface_v2["required_bindings"],
        *interface_v3["additional_required_bindings"],
    )
    paths = tuple(binding_path(name) for name in names)
    require(len(paths) == len(set(paths)), "duplicate prefreeze binding path")
    return paths


def validate_prefreeze_authorization() -> dict[str, Any]:
    require(
        DEVELOPMENT_PREFREEZE_AUTHORIZATION.is_file(),
        "DEVELOPMENT_PREFREEZE_AUTHORIZATION.v1.json is required; do not self-authorize",
    )
    authorization = read_json(DEVELOPMENT_PREFREEZE_AUTHORIZATION)
    require(
        authorization.get("development_preoutcome_prepare_authorized") is True,
        "preoutcome prepare is not authorized",
    )
    require(
        authorization.get("development_613_authorized") is False,
        "prefreeze authorization must not authorize development outcomes",
    )
    for key in (
        "untouched_410_test_authorized",
        "count_generation_authorized",
        "cell_eval2_authorized",
        "challenge_build_authorized",
        "submission_authorized",
    ):
        require(authorization.get(key) is False, f"unsafe prefreeze flag: {key}")
    validate_bound_paths(authorization, required_prefreeze_paths())
    return authorization


def source_manifest_directories() -> tuple[Path, ...]:
    return (
        *(SOURCE_BUILD_DIR / donor for donor in SOURCE_DONORS),
        *(SOURCE_EVALUATION_DIR / donor for donor in SOURCE_DONORS),
        SOURCE_GATE_DIR,
        EXP140_IDENTITY_DIR,
    )


def validate_source_plan_axes() -> None:
    """Tie every prepared row to the exact frozen development axes."""
    manifest = validate_directory_manifest(SOURCE_PLAN_DIR, SOURCE_PLAN_MANIFEST.name)
    require(
        manifest["development_outcomes_opened"] is False,
        "source plan opened development outcomes",
    )
    require(
        manifest["sealed_410_effect_values_opened"] is False,
        "source plan opened sealed outcomes",
    )
    require(
        manifest["development_targets"] == DEV_COUNT,
        "source plan development count changed",
    )
    require(
        manifest["response_genes"] == RESPONSE_COUNT,
        "source plan response count changed",
    )
    schema = read_json(SOURCE_PLAN_SCHEMA)
    require(
        schema["development_outcomes_opened"] is False,
        "source schema opened development outcomes",
    )
    require(
        schema["sealed_410_effect_values_opened"] is False,
        "source schema opened sealed outcomes",
    )
    parent_hashes = schema["parent_hashes"]
    for path in (DEVELOPMENT_TARGETS, RESPONSE_GENES):
        key = str(path.relative_to(ROOT))
        require(key in parent_hashes, f"source schema lacks frozen axis: {key}")
        require(
            sha256_file(path) == parent_hashes[key],
            f"frozen development axis changed: {path}",
        )


def validate_source_chain() -> None:
    validate_source_plan_axes()
    for directory in source_manifest_directories():
        validate_directory_manifest(directory)
    gate = read_json(SOURCE_GATE_DIR / RESULT_ARTIFACT)
    require(gate["source_gate_passes"] is True, "Experiment 142 source gate failed")
    require(gate["development_outcomes_opened"] is False, "source gate opened outcomes")
    require(gate["sealed_410_effect_values_opened"] is False, "source gate opened test")


def validate_frozen_geometry_inputs() -> None:
    """Revalidate only the compact donor inputs used to build development slots."""
    manifest = read_json(EXP140_SOURCE_MANIFEST)
    records = manifest["files"]
    required = [CANONICAL_GENES]
    for root in SOURCE_DONOR_ROOTS.values():
        required.extend(root / name for name in ("delta.npy", "perts.txt", "genes.txt"))
    for path in required:
        key = str(path.relative_to(ROOT))
        require(key in records, f"Experiment 140 source manifest lacks {key}")
        expected = records[key]
        require(
            path.stat().st_size == expected["bytes"], f"source bytes changed: {path}"
        )
        require(sha256_file(path) == expected["sha256"], f"source hash changed: {path}")


def block_grouped_folds(
    blocks: np.ndarray,
    *,
    fold_count: int = FOLDS,
    expected_sizes: Sequence[int] | None = EXPECTED_FOLD_SIZES,
) -> np.ndarray:
    """Deterministically assign intact blocks to the currently smallest fold."""
    labels = np.asarray(blocks, dtype=np.int64)
    require(labels.ndim == 1 and len(labels) > 0, "invalid development blocks")
    block_ids = np.unique(labels)
    groups = [(int(block), np.flatnonzero(labels == block)) for block in block_ids]
    groups.sort(key=lambda item: (-len(item[1]), item[0]))
    loads = np.zeros(fold_count, dtype=np.int64)
    folds = np.full(len(labels), -1, dtype=np.int8)
    for _, indices in groups:
        fold = min(range(fold_count), key=lambda value: (loads[value], value))
        folds[indices] = fold
        loads[fold] += len(indices)
    require(np.all(folds >= 0), "unassigned target fold")
    for block in block_ids:
        require(
            len(np.unique(folds[labels == block])) == 1,
            f"development block {block} was split",
        )
    if expected_sizes is not None:
        require(
            tuple(int(value) for value in loads) == tuple(expected_sizes),
            f"unexpected balanced fold sizes: {tuple(loads)}",
        )
    return folds


def average_tie_ranks(values: np.ndarray) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    require(data.ndim == 1, "rank input must be one-dimensional")
    order = np.argsort(data, kind="mergesort")
    ranks = np.empty(len(data), dtype=np.float64)
    start = 0
    while start < len(data):
        end = start + 1
        while end < len(data) and data[order[end]] == data[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1) + 1.0
        start = end
    return ranks


def midrank_z(values: np.ndarray, valid: np.ndarray | None = None) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    mask = np.isfinite(data) if valid is None else np.asarray(valid, dtype=bool)
    mask &= np.isfinite(data)
    result = np.full(len(data), np.nan, dtype=np.float64)
    count = int(mask.sum())
    require(count >= 2, "fewer than two valid values for midrank z")
    ranks = average_tie_ranks(data[mask])
    result[mask] = 2.0 * (ranks - 1.0) / (count - 1.0) - 1.0
    return result


def gated_midrank_z(values: np.ndarray, gate: np.ndarray) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    mask = np.asarray(gate, dtype=bool) & np.isfinite(data)
    result = np.zeros(len(data), dtype=np.float64)
    count = int(mask.sum())
    if count == 1:
        result[mask] = 1.0
    elif count >= 2:
        result[mask] = midrank_z(data[mask])
    return result


def normalized_weights(z: np.ndarray, eta: float) -> np.ndarray:
    values = np.asarray(z, dtype=np.float64)
    require(values.ndim == 2, "z must be target by donor")
    complete = np.all(np.isfinite(values), axis=1)
    raw = np.ones_like(values)
    raw[complete] = 1.0 + eta * values[complete]
    require(np.all(raw >= 0.5 - EPS), "weight fell below frozen lower bound")
    result = raw / raw.sum(axis=1, keepdims=True)
    common = complete & np.all(values == values[:, :1], axis=1)
    result[~complete | common] = 1.0 / values.shape[1]
    return result


def weight_grid(z: np.ndarray) -> np.ndarray:
    return np.stack([normalized_weights(z, eta) for eta in ETAS], axis=0)


def donor_assignment_null_weight_grid(
    observed_z: np.ndarray, permutations: np.ndarray
) -> np.ndarray:
    z = np.asarray(observed_z, dtype=np.float64)
    perm = np.asarray(permutations, dtype=np.int64)
    require(perm.ndim == 3 and perm.shape[1:] == z.shape, "donor-null shape")
    expected = np.arange(z.shape[1])
    require(
        np.all(np.sort(perm, axis=2) == expected[None, None, :]),
        "invalid donor-null permutation",
    )
    output = np.empty((len(perm), len(ETAS), *z.shape), dtype=np.float64)
    for draw in range(len(perm)):
        permuted = np.take_along_axis(z, perm[draw], axis=1)
        output[draw] = weight_grid(permuted)
    return output


def target_pairing_null_weight_grid(null_z: np.ndarray) -> np.ndarray:
    values = np.asarray(null_z, dtype=np.float64)
    require(values.ndim == 3, "target-null z must be draw by target by donor")
    output = np.empty(
        (values.shape[0], len(ETAS), values.shape[1], values.shape[2]),
        dtype=np.float64,
    )
    for draw in range(values.shape[0]):
        output[draw] = weight_grid(values[draw])
    return output


def load_operator(group: str) -> tuple[np.ndarray, np.ndarray]:
    with np.load(
        EXP140_IDENTITY_DIR / f"whitener-{group}.npz", allow_pickle=False
    ) as data:
        evals = np.asarray(data["evals"], dtype=np.float64)
        evecs = np.asarray(data["evecs"], dtype=np.float32)
    require(evals.shape == (RESPONSE_COUNT,), f"{group}: eigenvalue shape")
    require(evecs.shape == (RESPONSE_COUNT, RESPONSE_COUNT), f"{group}: basis shape")
    return evals, evecs


def apply_whitener(
    matrix: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    *,
    chunk: int = 64,
) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    output = np.empty(values.shape, dtype=np.float32)
    spectral_weight = evals ** (-WHITEN_GAMMA / 2.0)
    for lo in range(0, len(values), chunk):
        hi = min(lo + chunk, len(values))
        spectral = (values[lo:hi] @ evecs) * spectral_weight[None, :]
        output[lo:hi] = (spectral @ evecs.T).astype(np.float32)
    return output


def spectral_coordinates(
    matrix: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    *,
    chunk: int = 64,
) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float32)
    output = np.empty(values.shape, dtype=np.float32)
    spectral_weight = evals ** (-WHITEN_GAMMA / 2.0)
    for lo in range(0, len(values), chunk):
        hi = min(lo + chunk, len(values))
        output[lo:hi] = (
            np.asarray(values[lo:hi], dtype=np.float64)
            @ evecs
            * spectral_weight[None, :]
        ).astype(np.float32)
    return output


def pds_matrix(query: np.ndarray, gallery: np.ndarray) -> np.ndarray:
    left = np.asarray(query, dtype=np.float64)
    right = np.asarray(gallery, dtype=np.float64)
    left /= np.maximum(np.linalg.norm(left, axis=1, keepdims=True), EPS)
    right /= np.maximum(np.linalg.norm(right, axis=1, keepdims=True), EPS)
    similarity = left @ right.T
    output = np.empty_like(similarity)
    denominator = max(similarity.shape[1] - 1, 1)
    for row in range(len(similarity)):
        values = similarity[row]
        order = np.argsort(-values, kind="mergesort")
        start = 0
        while start < len(values):
            end = start + 1
            while end < len(values) and values[order[end]] == values[order[start]]:
                end += 1
            midrank = 1.0 + start + 0.5 * (end - start - 1)
            output[row, order[start:end]] = 1.0 - (midrank - 1.0) / denominator
            start = end
    return output


def recompute_target_null_z(
    half_centered: np.ndarray,
    valid_pairs: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    permutations: np.ndarray,
    target_null_p95: np.ndarray,
) -> np.ndarray:
    """Recompute all 49 false-pair reliability vectors from frozen source rows."""
    centered = np.asarray(half_centered)
    valid = np.asarray(valid_pairs, dtype=bool)
    perm = np.asarray(permutations, dtype=np.int64)
    thresholds = np.asarray(target_null_p95, dtype=np.float64)
    splits, halves, targets, genes = centered.shape
    require(halves == 2, "half axis changed")
    require(valid.shape == (splits, targets), "valid-pair shape changed")
    require(perm.shape == (NULL_DRAWS, targets), "matched-null shape changed")
    require(thresholds.shape == (targets,), "target-null threshold shape changed")
    require(evals.shape == (genes,) and evecs.shape == (genes, genes), "operator shape")

    directional = np.full(
        (NULL_DRAWS, splits * halves, targets), np.nan, dtype=np.float64
    )
    all_indices = np.arange(targets, dtype=np.int64)
    for split in range(splits):
        indices = all_indices[valid[split]]
        require(len(indices) >= 2, f"split {split}: insufficient gallery")
        transformed = [
            spectral_coordinates(centered[split, half, indices], evals, evecs)
            for half in range(halves)
        ]
        local_position = np.full(targets, -1, dtype=np.int64)
        local_position[indices] = np.arange(len(indices), dtype=np.int64)
        for direction, (query, gallery) in enumerate(
            ((transformed[0], transformed[1]), (transformed[1], transformed[0]))
        ):
            matrix = pds_matrix(query, gallery)
            for draw in range(NULL_DRAWS):
                designated = local_position[perm[draw, indices]]
                usable = designated >= 0
                rows = np.flatnonzero(usable)
                directional[draw, 2 * split + direction, indices[rows]] = matrix[
                    rows, designated[usable]
                ]

    counts = np.sum(np.isfinite(directional), axis=1)
    reliability = np.nanquantile(directional, 0.25, axis=1, method="linear")
    reliability[counts < MIN_DIRECTIONAL_VALUES] = np.nan
    z = np.empty_like(reliability)
    for draw in range(NULL_DRAWS):
        gate = (
            np.isfinite(reliability[draw])
            & np.isfinite(thresholds)
            & (reliability[draw] > thresholds)
        )
        z[draw] = gated_midrank_z(reliability[draw], gate)
    return z


def last_columns(source: Sequence[str], wanted: Sequence[str]) -> np.ndarray:
    index = {gene: position for position, gene in enumerate(source)}
    require(all(gene in index for gene in wanted), "source misses a response gene")
    return np.asarray([index[gene] for gene in wanted], dtype=np.int64)


def centered_source_rows(
    donor: str, targets: Sequence[str], responses: Sequence[str]
) -> np.ndarray:
    root = SOURCE_DONOR_ROOTS[donor]
    perts = read_lines(root / "perts.txt")
    genes = read_lines(root / "genes.txt", unique=False)
    target_index = {target: row for row, target in enumerate(perts)}
    require(
        all(target in target_index for target in targets), f"{donor}: missing target"
    )
    rows = np.asarray([target_index[target] for target in targets], dtype=np.int64)
    columns = last_columns(genes, responses)
    delta = np.load(root / "delta.npy", mmap_mode="r")
    require(delta.shape == (len(perts), len(genes)), f"{donor}: delta shape")
    mean = np.asarray(delta.mean(axis=0), dtype=np.float32)[columns]
    selected = np.asarray(delta[np.ix_(rows, columns)], dtype=np.float32)
    return selected - mean


def prepare_perturbmatch(
    path: Path,
    target_genes: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, Any]:
    """Load only the frozen PerturbMatch arrays needed by the v3 evaluator."""
    genes = read_lines(path / "genes.txt")
    perts = read_lines(path / "perts.txt")
    metadata = read_json(path / "meta.json")
    require(
        metadata.get("source_md5") == "0017076a5cc7b9979febf3bd19bd2157",
        "PerturbMatch source MD5 changed",
    )
    target_position = {gene: index for index, gene in enumerate(target_genes)}
    source_for_target: dict[int, int] = {}
    for source_index, gene in enumerate(genes):
        target_index = target_position.get(gene)
        if target_index is not None:
            source_for_target[target_index] = source_index
    target_columns = np.asarray(sorted(source_for_target), dtype=np.int64)
    source_columns = np.asarray(
        [source_for_target[index] for index in target_columns], dtype=np.int64
    )
    require(len(target_columns) > 0, "PerturbMatch has no canonical genes")

    delta = np.load(path / "delta.npy", mmap_mode="r")
    require(delta.shape == (len(perts), len(genes)), "PerturbMatch shape changed")
    mean = np.zeros(len(source_columns), dtype=np.float64)
    for lo in range(0, len(perts), chunk_rows):
        hi = min(lo + chunk_rows, len(perts))
        mean += np.asarray(
            delta[np.ix_(np.arange(lo, hi), source_columns)], dtype=np.float64
        ).sum(axis=0)
        memory_check("PerturbMatch center", max_rss_gib)
    mean /= len(perts)

    norms = np.empty(len(perts), dtype=np.float64)
    for lo in range(0, len(perts), chunk_rows):
        hi = min(lo + chunk_rows, len(perts))
        block = np.asarray(
            delta[np.ix_(np.arange(lo, hi), source_columns)], dtype=np.float64
        )
        block -= mean
        norms[lo:hi] = np.linalg.norm(block, axis=1)
        memory_check("PerturbMatch scale", max_rss_gib)
    scale = float(np.median(norms)) or 1.0
    return {
        "genes": genes,
        "perts": perts,
        "pert_index": {pert: index for index, pert in enumerate(perts)},
        "target_columns": target_columns,
        "source_columns": source_columns,
        "mean": mean,
        "scale": scale,
        "delta": delta,
        "metadata": metadata,
    }


def composite_k562_standardized_rows(
    k562_centered: np.ndarray,
    targets: Sequence[str],
    responses: Sequence[str],
    canonical_genes: Sequence[str],
    k562_scale: float,
    pm: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, int]]:
    """Apply the production 50:50 PM replacement on response coordinates only."""
    base = np.asarray(k562_centered, dtype=np.float64) / float(k562_scale)
    output = base.copy()
    canonical_position = {gene: index for index, gene in enumerate(canonical_genes)}
    response_position = {gene: index for index, gene in enumerate(responses)}
    pm_target_columns = np.asarray(pm["target_columns"], dtype=np.int64)
    pm_source_columns = np.asarray(pm["source_columns"], dtype=np.int64)
    pm_by_canonical = {
        int(target): int(source)
        for target, source in zip(pm_target_columns, pm_source_columns, strict=True)
    }
    response_columns = []
    source_columns = []
    for gene, response_column in response_position.items():
        canonical_column = canonical_position.get(gene)
        if canonical_column is not None and canonical_column in pm_by_canonical:
            response_columns.append(response_column)
            source_columns.append(pm_by_canonical[canonical_column])
    response_columns_array = np.asarray(response_columns, dtype=np.int64)
    source_columns_array = np.asarray(source_columns, dtype=np.int64)
    pm_mean_by_source = {
        int(source): float(mean)
        for source, mean in zip(pm_source_columns, np.asarray(pm["mean"]), strict=True)
    }
    response_pm_mean = np.asarray(
        [pm_mean_by_source[int(source)] for source in source_columns_array],
        dtype=np.float64,
    )
    pm_scale = float(pm["scale"])
    pm_delta = pm["delta"]
    pm_pert_index = pm["pert_index"]
    modified = 0
    missing_target = 0
    for row, target in enumerate(targets):
        pm_row = pm_pert_index.get(str(target))
        if pm_row is None:
            missing_target += 1
            continue
        before = output[row].copy()
        pm_standardized = (
            np.asarray(pm_delta[pm_row, source_columns_array], dtype=np.float64)
            - response_pm_mean
        ) / pm_scale
        output[row, response_columns_array] = (1.0 - PM_BLEND) * base[
            row, response_columns_array
        ] + PM_BLEND * pm_standardized
        missing = np.ones(len(responses), dtype=bool)
        missing[response_columns_array] = False
        require(
            np.array_equal(output[row, missing], before[missing]),
            f"{target}: PM-missing response coordinate changed",
        )
        modified += 1
    return output, {
        "targets_modified": modified,
        "targets_preserved_exactly": missing_target,
        "response_genes_replaced": len(response_columns_array),
    }


def frozen_group_parameters() -> dict[str, dict[str, Any]]:
    result = read_json(EXP140_IDENTITY_DIR / RESULT_ARTIFACT)
    require(result["identity_reproduction_passes"] is True, "identity result failed")
    return result["groups"]


def build_geometry(
    group: str,
    centered: Mapping[str, np.ndarray],
    targets: Sequence[str],
    responses: Sequence[str],
    canonical: Sequence[str],
    pm: Mapping[str, Any],
    parameters: Mapping[str, Mapping[str, Any]],
    limit: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    donors = GROUP_SLOTS[group]
    source_parameters = parameters[group]
    renorm = source_parameters["donor_renormalization"]
    evals, evecs = load_operator(group)
    pre_slots = []
    pm_report: dict[str, Any] | None = None
    for slot in donors:
        source = SLOT_SOURCE[slot]
        if slot == "k562_composite":
            scale = float(source_parameters["donor_raw_median_norms"]["k562_gwps"])
            standardized, pm_report = composite_k562_standardized_rows(
                centered["k562_gwps"],
                targets,
                responses,
                canonical,
                scale,
                pm,
            )
            reference_norm = scale * float(renorm["k562_gwps"])
            values = reference_norm * standardized
        else:
            values = float(renorm[source]) * np.asarray(
                centered[source], dtype=np.float64
            )
        pre_slots.append(values)
    if group == "non_k562":
        require(pm_report is not None, "non-K562 geometry lacks PerturbMatch report")
        require(
            pm_report["targets_modified"] == EXPECTED_PM_TARGET_SUPPORT,
            "PerturbMatch target support changed",
        )
        require(
            pm_report["response_genes_replaced"] == EXPECTED_PM_RESPONSE_SUPPORT,
            "PerturbMatch response support changed",
        )
    pre = np.stack(pre_slots, axis=1)
    slots = np.empty(pre.shape, dtype=np.float32)
    for donor_index in range(len(donors)):
        slots[:, donor_index] = apply_whitener(pre[:, donor_index], evals, evecs)
        memory_check(f"{group} slot {donors[donor_index]}", limit)
    baseline = apply_whitener(np.mean(pre, axis=1), evals, evecs)
    return {
        "donors": np.asarray(donors, dtype="U32"),
        "slots": slots,
        "baseline": baseline,
    }, {"perturbmatch": pm_report}


def role_source(group: str, slot: str) -> tuple[str, str]:
    return SLOT_SOURCE[slot], group


def load_raw_signal(group: str, signal_key: str) -> np.ndarray:
    values = []
    for slot in GROUP_SLOTS[group]:
        donor, role = role_source(group, slot)
        with np.load(
            SOURCE_EVALUATION_DIR / donor / f"{role}.npz", allow_pickle=False
        ) as data:
            values.append(np.asarray(data[signal_key], dtype=np.float64))
    result = np.column_stack(values)
    require(result.shape == (DEV_COUNT, len(GROUP_SLOTS[group])), "signal shape")
    return result


def recompute_group_target_null_z(group: str, limit: float) -> np.ndarray:
    donor_values = []
    for slot in GROUP_SLOTS[group]:
        donor, role = role_source(group, slot)
        source_dir = SOURCE_BUILD_DIR / donor
        with np.load(
            SOURCE_EVALUATION_DIR / donor / f"{role}.npz", allow_pickle=False
        ) as evaluation:
            permutations = np.asarray(
                evaluation["matched_null_permutations"], dtype=np.int16
            )
            thresholds = np.asarray(evaluation["target_null_p95"], dtype=np.float64)
        centered = np.load(source_dir / "half-centered.npy", mmap_mode="r")
        valid = np.asarray(np.load(source_dir / "valid-pairs.npy"), dtype=bool)
        evals, evecs = load_operator(role)
        donor_values.append(
            recompute_target_null_z(
                centered, valid, evals, evecs, permutations, thresholds
            )
        )
        memory_check(f"{group} target-pairing null {slot}", limit)
    result = np.stack(donor_values, axis=2)
    require(
        result.shape == (NULL_DRAWS, DEV_COUNT, len(GROUP_SLOTS[group])),
        f"{group}: target-null z shape",
    )
    return result


def source_truth_hash_ledger() -> dict[str, Any]:
    """Copy frozen truth hashes without loading or hashing any truth matrix."""
    authorization = read_json(EXP132_AUTHORIZATION)
    inputs = authorization["inputs"]
    datasets = {}
    for dataset, root in TRUTH_ROOTS.items():
        records = {}
        for name in ("delta.npy", "perts.txt", "genes.txt"):
            path = root / name
            key = str(path.relative_to(ROOT))
            require(key in inputs, f"Experiment 132 authorization lacks {key}")
            expected = inputs[key]
            records[name] = {
                "path": key,
                "bytes": int(expected["bytes"]),
                "sha256": str(expected["sha256"]),
            }
        datasets[dataset] = records
    return {
        "schema_version": 1,
        "source_authorization": {
            "path": str(EXP132_AUTHORIZATION.relative_to(ROOT)),
            "bytes": EXP132_AUTHORIZATION.stat().st_size,
            "sha256": sha256_file(EXP132_AUTHORIZATION),
        },
        "source_plan_manifest": {
            "path": str((EXP132_PLAN / "plan-manifest.json").relative_to(ROOT)),
            "bytes": (EXP132_PLAN / "plan-manifest.json").stat().st_size,
            "sha256": sha256_file(EXP132_PLAN / "plan-manifest.json"),
        },
        "development_targets": {
            "path": str(DEVELOPMENT_TARGETS.relative_to(ROOT)),
            "bytes": DEVELOPMENT_TARGETS.stat().st_size,
            "sha256": sha256_file(DEVELOPMENT_TARGETS),
        },
        "response_genes": {
            "path": str(RESPONSE_GENES.relative_to(ROOT)),
            "bytes": RESPONSE_GENES.stat().st_size,
            "sha256": sha256_file(RESPONSE_GENES),
        },
        "datasets": datasets,
        "truth_arrays_opened": False,
    }


def prepare(limit: float) -> None:
    authorization = validate_prefreeze_authorization()
    created_utc = authorized_provenance_utc(authorization, "prepare")
    validate_source_chain()
    validate_frozen_geometry_inputs()
    require(not PREOUTCOME_DIR.exists(), f"preoutcome output exists: {PREOUTCOME_DIR}")
    require(not RESULT_DIR.exists(), "development result already exists")
    require(
        not OUTCOME_ACCESS_LOCK.exists(), "development outcomes were already accessed"
    )
    targets = read_lines(DEVELOPMENT_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    require(len(targets) == DEV_COUNT, "development target count changed")
    require(len(responses) == RESPONSE_COUNT, "response gene count changed")

    with (
        exclusive_lock(PREPARE_LOCK),
        atomic_stage_directory(PREOUTCOME_DIR) as staging,
    ):
        with np.load(EXP140_PLAN_ARRAYS, allow_pickle=False) as plan:
            blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
            bootstrap = np.asarray(plan["bootstrap_block_draws"], dtype=np.uint16)
        folds = block_grouped_folds(blocks)
        require(bootstrap.shape[0] == BOOTSTRAP_DRAWS, "bootstrap draw count changed")
        atomic_npz(
            staging / PLAN_ARTIFACT,
            development_folds=folds,
            development_blocks=blocks,
            bootstrap_block_draws=bootstrap,
            etas=np.asarray(ETAS, dtype=np.float64),
        )

        centered = {
            donor: centered_source_rows(donor, targets, responses)
            for donor in SOURCE_DONORS
        }
        canonical = read_csv_first_column(CANONICAL_GENES)
        pm = prepare_perturbmatch(
            PERTURBMATCH_ROOT,
            canonical,
            chunk_rows=128,
            max_rss_gib=limit,
        )
        parameters = frozen_group_parameters()
        geometry_reports = {}

        with np.load(
            SOURCE_GATE_DIR / "candidate-weights-k562.npz", allow_pickle=False
        ) as source_k562:
            frozen_donor_null_k562 = np.asarray(
                source_k562["donor_assignment_nulls"], dtype=np.int8
            )
        with np.load(
            SOURCE_GATE_DIR / "candidate-weights-non_k562.npz", allow_pickle=False
        ) as source_non_k562:
            frozen_donor_null_non_k562 = np.asarray(
                source_non_k562["donor_assignment_nulls"], dtype=np.int8
            )
        frozen_donor_nulls = {
            "k562": frozen_donor_null_k562,
            "non_k562": frozen_donor_null_non_k562,
        }

        for group in ("k562", "non_k562"):
            geometry, report = build_geometry(
                group,
                centered,
                targets,
                responses,
                canonical,
                pm,
                parameters,
                limit,
            )
            for signal, source_key in SIGNALS.items():
                raw_z = load_raw_signal(group, source_key)
                geometry[f"raw_z__{signal}"] = raw_z
                geometry[f"weights__{signal}"] = weight_grid(raw_z)
            require(
                np.array_equal(
                    mix_slots(
                        geometry["slots"],
                        geometry["weights__reliability"][0],
                        0.0,
                        geometry["baseline"],
                    ),
                    geometry["baseline"],
                ),
                f"{group}: eta-zero baseline is not exact",
            )
            atomic_npz(staging / GEOMETRY_ARTIFACT[group], **geometry)

            target_null_z = recompute_group_target_null_z(group, limit)
            target_null_weights = target_pairing_null_weight_grid(target_null_z)
            atomic_npz(
                staging / TARGET_NULL_ARTIFACT[group],
                donors=geometry["donors"],
                etas=np.asarray(ETAS, dtype=np.float64),
                null_z=target_null_z,
                weights=target_null_weights,
            )
            donor_null_weights = donor_assignment_null_weight_grid(
                geometry["raw_z__reliability"], frozen_donor_nulls[group]
            )
            atomic_npz(
                staging / DONOR_NULL_ARTIFACT[group],
                donors=geometry["donors"],
                etas=np.asarray(ETAS, dtype=np.float64),
                permutations=frozen_donor_nulls[group],
                weights=donor_null_weights,
            )
            geometry_reports[group] = {
                "donors": list(map(str, geometry["donors"])),
                "eta_zero_exact": True,
                **report,
            }

        atomic_json(staging / TRUTH_HASH_ARTIFACT, source_truth_hash_ledger())
        result = {
            "schema_version": 1,
            "stage": "development_preoutcome_prepare",
            "created_utc": created_utc,
            "development_outcomes_opened": False,
            "sealed_410_effect_values_opened": False,
            "truth_arrays_opened": False,
            "prefreeze_authorization_sha256": sha256_file(
                DEVELOPMENT_PREFREEZE_AUTHORIZATION
            ),
            "prefreeze_authorization_schema": authorization.get("schema_version"),
            "fold_sizes": [int(np.sum(folds == fold)) for fold in range(FOLDS)],
            "geometry": geometry_reports,
            "target_pairing_null_draws": NULL_DRAWS,
            "donor_assignment_null_draws": NULL_DRAWS,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "peak_rss_gib": peak_rss_gib(),
        }
        atomic_json(staging / RESULT_ARTIFACT, result)
        atomic_json(staging / MANIFEST_ARTIFACT, directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_preoutcome() -> dict[str, Any]:
    authorization = validate_prefreeze_authorization()
    validate_source_chain()
    validate_directory_manifest(PREOUTCOME_DIR)
    result = read_json(PREOUTCOME_DIR / RESULT_ARTIFACT)
    require(
        result["prefreeze_authorization_sha256"]
        == sha256_file(DEVELOPMENT_PREFREEZE_AUTHORIZATION),
        "preoutcome artifact was created under a different authorization",
    )
    require(
        result["prefreeze_authorization_schema"] == authorization.get("schema_version"),
        "preoutcome authorization schema changed",
    )
    require(result["development_outcomes_opened"] is False, "prepare opened outcomes")
    require(result["sealed_410_effect_values_opened"] is False, "prepare opened test")
    require(result["truth_arrays_opened"] is False, "prepare opened truth arrays")
    require(tuple(result["fold_sizes"]) == EXPECTED_FOLD_SIZES, "fold sizes changed")
    return result


def required_development_paths() -> tuple[Path, ...]:
    return (
        *required_prefreeze_paths(),
        DEVELOPMENT_PREFREEZE_AUTHORIZATION,
        PREOUTCOME_DIR / MANIFEST_ARTIFACT,
        PREOUTCOME_DIR / TRUTH_HASH_ARTIFACT,
        EXP132_RUN,
        EXP132_AUTHORIZATION,
    )


def validate_ledger_binding(
    record: Mapping[str, Any],
    path: Path,
    bindings: Mapping[Path, str],
    *,
    verify_file: bool,
) -> None:
    require(
        record["path"] == str(path.relative_to(ROOT)), f"ledger path changed: {path}"
    )
    require(path.resolve() in bindings, f"ledger source not bound: {path}")
    require(
        bindings[path.resolve()] == record["sha256"],
        f"authorization differs from frozen ledger: {path}",
    )
    if verify_file:
        require(path.stat().st_size == record["bytes"], f"bytes changed: {path}")
        require(sha256_file(path) == record["sha256"], f"hash changed: {path}")


def validate_development_authorization(*, verify_truth_files: bool) -> dict[str, Any]:
    validate_preoutcome()
    require(
        DEVELOPMENT_AUTHORIZATION.is_file(),
        "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json is required; do not self-authorize",
    )
    authorization = read_json(DEVELOPMENT_AUTHORIZATION)
    require(
        authorization.get("development_613_authorized") is True,
        "613-target development access is not authorized",
    )
    for key in (
        "untouched_410_test_authorized",
        "count_generation_authorized",
        "cell_eval2_authorized",
        "challenge_build_authorized",
        "submission_authorized",
    ):
        require(authorization.get(key) is False, f"unsafe development flag: {key}")
    bindings = validate_bound_paths(authorization, required_development_paths())
    ledger = read_json(PREOUTCOME_DIR / TRUTH_HASH_ARTIFACT)
    for key, path in (
        ("source_authorization", EXP132_AUTHORIZATION),
        ("source_plan_manifest", EXP132_PLAN / "plan-manifest.json"),
        ("development_targets", DEVELOPMENT_TARGETS),
        ("response_genes", RESPONSE_GENES),
    ):
        validate_ledger_binding(
            ledger[key],
            path,
            bindings,
            verify_file=verify_truth_files,
        )
    for records in ledger["datasets"].values():
        for record in records.values():
            path = ROOT / record["path"]
            validate_ledger_binding(
                record,
                path,
                bindings,
                verify_file=verify_truth_files,
            )
    return authorization


def load_prepared_plan() -> dict[str, np.ndarray]:
    with np.load(PREOUTCOME_DIR / PLAN_ARTIFACT, allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def load_prepared_geometry(group: str) -> dict[str, np.ndarray]:
    with np.load(PREOUTCOME_DIR / GEOMETRY_ARTIFACT[group], allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def load_null_weights(kind: str, group: str) -> np.ndarray:
    paths = {
        "target_pairing": TARGET_NULL_ARTIFACT,
        "donor_assignment": DONOR_NULL_ARTIFACT,
    }
    require(kind in paths, f"unknown null family: {kind}")
    with np.load(PREOUTCOME_DIR / paths[kind][group], allow_pickle=False) as data:
        return np.asarray(data["weights"], dtype=np.float64)


def mix_slots(
    slots: np.ndarray,
    weights: np.ndarray,
    eta: float,
    baseline: np.ndarray,
) -> np.ndarray:
    values = np.asarray(slots)
    w = np.asarray(weights, dtype=np.float64)
    frozen = np.asarray(baseline)
    require(values.ndim == 3, "slots must be target by donor by gene")
    require(w.shape == values.shape[:2], "weight/slot shape mismatch")
    require(frozen.shape == (values.shape[0], values.shape[2]), "baseline shape")
    if eta == 0.0:
        return frozen.copy()
    require(
        np.allclose(w.sum(axis=1), 1.0, rtol=0.0, atol=EPS), "weights not normalized"
    )
    return np.einsum("td,tdg->tg", w, values, optimize=True)


def pds_subset(
    prediction: np.ndarray,
    truth: np.ndarray,
    query: np.ndarray,
    gallery: np.ndarray,
) -> np.ndarray:
    pred = np.asarray(prediction[query], dtype=np.float64)
    observed = np.asarray(truth[gallery], dtype=np.float64)
    pred /= np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), EPS)
    observed /= np.maximum(np.linalg.norm(observed, axis=1, keepdims=True), EPS)
    similarity = pred @ observed.T
    positions = {int(value): index for index, value in enumerate(gallery)}
    own_columns = np.asarray([positions[int(value)] for value in query], dtype=np.int64)
    own = similarity[np.arange(len(query)), own_columns]
    greater = np.sum(similarity > own[:, None], axis=1)
    equal = np.sum(similarity == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (len(gallery) - 1.0)


def signed_credit(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    numerator = np.sum(np.sign(truth) * prediction, axis=1)
    denominator = np.sum(np.abs(prediction), axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=denominator > EPS,
    )


def predictions_from_weights(
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    weights_by_group: Mapping[str, np.ndarray],
) -> dict[float, dict[str, np.ndarray]]:
    output: dict[float, dict[str, np.ndarray]] = {eta: {} for eta in ETAS}
    for group, geometry in geometries.items():
        weights = np.asarray(weights_by_group[group], dtype=np.float64)
        require(
            weights.shape[:2] == (len(ETAS), DEV_COUNT),
            f"{group}: weight grid shape",
        )
        for eta_index, eta in enumerate(ETAS):
            output[eta][group] = mix_slots(
                geometry["slots"],
                weights[eta_index],
                eta,
                geometry["baseline"],
            )
    return output


def score_predictions(
    predictions: Mapping[float, Mapping[str, np.ndarray]],
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
) -> dict[float, dict[str, dict[str, np.ndarray]]]:
    indices = np.arange(DEV_COUNT, dtype=np.int64)
    identity_pds = {}
    identity_credit = {}
    for dataset in SCORED_DATASETS:
        group = DATASET_GROUP[dataset]
        identity = geometries[group]["baseline"]
        identity_pds[dataset] = pds_subset(identity, truth[dataset], indices, indices)
        identity_credit[dataset] = signed_credit(identity, truth[dataset])
    scores = {}
    for eta in ETAS:
        scores[eta] = {}
        for dataset in SCORED_DATASETS:
            group = DATASET_GROUP[dataset]
            prediction = predictions[eta][group]
            scores[eta][dataset] = {
                "pds": (
                    pds_subset(prediction, truth[dataset], indices, indices)
                    - identity_pds[dataset]
                ),
                "signed_credit": (
                    signed_credit(prediction, truth[dataset]) - identity_credit[dataset]
                ),
            }
    return scores


def aggregate_scores(
    score: Mapping[str, Mapping[str, np.ndarray]], indices: np.ndarray
) -> dict[str, Any]:
    idx = np.asarray(indices, dtype=np.int64)
    family_arrays: dict[str, dict[str, np.ndarray]] = {
        metric: {
            "k562": 0.5
            * (
                np.asarray(score["k562_gwps"][metric])[idx]
                + np.asarray(score["k562_essential"][metric])[idx]
            ),
            "rpe1": np.asarray(score["rpe1"][metric])[idx],
            "jurkat": np.asarray(score["nadig_jurkat"][metric])[idx],
        }
        for metric in ("pds", "signed_credit")
    }
    families = {
        metric: {
            family: float(values.mean())
            for family, values in family_arrays[metric].items()
        }
        for metric in family_arrays
    }
    per_target = {
        metric: np.mean(np.stack(list(family_arrays[metric].values())), axis=0)
        for metric in family_arrays
    }
    return {
        "families": families,
        "macro": {metric: float(per_target[metric].mean()) for metric in per_target},
        "per_target": per_target,
    }


def selection_summary(
    eta: float,
    indices: np.ndarray,
    predictions: Mapping[float, Mapping[str, np.ndarray]],
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    *,
    inner_gallery: bool,
) -> dict[str, Any]:
    idx = np.asarray(indices, dtype=np.int64)
    if not inner_gallery:
        return aggregate_scores(scores[eta], idx)
    local = {}
    for dataset in SCORED_DATASETS:
        group = DATASET_GROUP[dataset]
        local[dataset] = {
            "pds": (
                pds_subset(predictions[eta][group], truth[dataset], idx, idx)
                - pds_subset(geometries[group]["baseline"], truth[dataset], idx, idx)
            ),
            "signed_credit": scores[eta][dataset]["signed_credit"][idx],
        }
    return aggregate_scores(local, np.arange(len(idx), dtype=np.int64))


def select_eta(
    indices: np.ndarray,
    predictions: Mapping[float, Mapping[str, np.ndarray]],
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    *,
    inner_gallery: bool,
) -> tuple[float, dict[float, dict[str, Any]]]:
    summaries = {
        eta: selection_summary(
            eta,
            indices,
            predictions,
            scores,
            geometries,
            truth,
            inner_gallery=inner_gallery,
        )
        for eta in ETAS
    }
    eligible = [
        eta
        for eta in ETAS
        if all(
            value >= -EPS
            for metric in ("pds", "signed_credit")
            for value in summaries[eta]["families"][metric].values()
        )
    ]
    require(0.0 in eligible, "eta-zero identity became ineligible")
    selected = 0.0
    for eta in eligible:
        candidate = summaries[eta]["macro"]
        incumbent = summaries[selected]["macro"]
        if (
            candidate["pds"],
            candidate["signed_credit"],
            -eta,
        ) > (
            incumbent["pds"],
            incumbent["signed_credit"],
            -selected,
        ):
            selected = eta
    return selected, summaries


def nested_and_fixed_selection(
    predictions: Mapping[float, Mapping[str, np.ndarray]],
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    folds: np.ndarray,
) -> dict[str, Any]:
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    oof = {
        dataset: {
            metric: np.empty(DEV_COUNT, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }
    records = []
    for fold in range(FOLDS):
        held = all_indices[folds == fold]
        inner = all_indices[folds != fold]
        eta, summaries = select_eta(
            inner,
            predictions,
            scores,
            geometries,
            truth,
            inner_gallery=True,
        )
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                oof[dataset][metric][held] = scores[eta][dataset][metric][held]
        records.append(
            {
                "fold": fold,
                "selected_eta": eta,
                "inner_macro": summaries[eta]["macro"],
            }
        )
    final_eta, summaries = select_eta(
        all_indices,
        predictions,
        scores,
        geometries,
        truth,
        inner_gallery=False,
    )
    return {
        "folds": records,
        "oof": oof,
        "oof_summary": aggregate_scores(oof, all_indices),
        "final_eta": final_eta,
        "fixed_score": scores[final_eta],
        "fixed_summary": summaries[final_eta],
    }


def fold_family_cells(
    score: Mapping[str, Mapping[str, np.ndarray]], folds: np.ndarray
) -> dict[str, dict[str, list[float]]]:
    result = {
        metric: {family: [] for family in FAMILIES}
        for metric in ("pds", "signed_credit")
    }
    for fold in range(FOLDS):
        indices = np.flatnonzero(folds == fold)
        summary = aggregate_scores(score, indices)
        for metric in result:
            for family in FAMILIES:
                result[metric][family].append(summary["families"][metric][family])
    return result


def block_bootstrap(
    values: np.ndarray, blocks: np.ndarray, draws: np.ndarray
) -> np.ndarray:
    labels = np.asarray(blocks, dtype=np.int64)
    samples = np.asarray(draws, dtype=np.int64)
    block_count = int(labels.max()) + 1
    require(samples.ndim == 2 and samples.shape[1] == block_count, "bootstrap shape")
    sums = np.bincount(
        labels, weights=np.asarray(values, dtype=np.float64), minlength=block_count
    )
    sizes = np.bincount(labels, minlength=block_count).astype(np.float64)
    return sums[samples].sum(axis=1) / sizes[samples].sum(axis=1)


def bootstrap_summary(
    score: Mapping[str, Mapping[str, np.ndarray]],
    blocks: np.ndarray,
    draws: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    summary = aggregate_scores(score, np.arange(DEV_COUNT, dtype=np.int64))
    report = {}
    arrays = {}
    for metric in ("pds", "signed_credit"):
        samples = block_bootstrap(summary["per_target"][metric], blocks, draws)
        arrays[metric] = samples
        report[metric] = {
            "lower_95": float(np.quantile(samples, 0.025, method="linear")),
            "upper_95": float(np.quantile(samples, 0.975, method="linear")),
        }
    return report, arrays


def positive_share(values: Sequence[float] | np.ndarray) -> float:
    positive = np.maximum(np.asarray(values, dtype=np.float64), 0.0)
    return 1.0 if positive.sum() == 0.0 else float(positive.max() / positive.sum())


def concentration_summary(
    score: Mapping[str, Mapping[str, np.ndarray]], blocks: np.ndarray
) -> dict[str, Any]:
    summary = aggregate_scores(score, np.arange(DEV_COUNT, dtype=np.int64))
    output = {}
    for metric in ("pds", "signed_credit"):
        block_values = np.bincount(
            np.asarray(blocks, dtype=np.int64),
            weights=summary["per_target"][metric],
        )
        output[metric] = {
            "family": positive_share(list(summary["families"][metric].values())),
            "target_block": positive_share(block_values),
        }
    return output


def score_weight_family(
    weights_by_group: Mapping[str, np.ndarray],
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    folds: np.ndarray,
) -> dict[str, Any]:
    predictions = predictions_from_weights(geometries, weights_by_group)
    scores = score_predictions(predictions, geometries, truth)
    selected = nested_and_fixed_selection(predictions, scores, geometries, truth, folds)
    return {"predictions": predictions, "scores": scores, "selection": selected}


def evaluate_null_family(
    kind: str,
    geometries: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    folds: np.ndarray,
    limit: float,
) -> tuple[dict[str, dict[str, np.ndarray]], dict[str, Any]]:
    weights = {group: load_null_weights(kind, group) for group in GROUP_SLOTS}
    values = {
        scope: {
            metric: np.empty(NULL_DRAWS, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for scope in ("nested", "fixed")
    }
    for draw in range(NULL_DRAWS):
        result = score_weight_family(
            {group: weights[group][draw] for group in GROUP_SLOTS},
            geometries,
            truth,
            folds,
        )["selection"]
        for metric in ("pds", "signed_credit"):
            values["nested"][metric][draw] = result["oof_summary"]["macro"][metric]
            values["fixed"][metric][draw] = result["fixed_summary"]["macro"][metric]
        if draw % 7 == 6:
            memory_check(f"{kind} null {draw + 1}", limit)
    report = {
        scope: {
            metric: {
                "p95": float(np.quantile(values[scope][metric], 0.95, method="linear")),
                "draws": NULL_DRAWS,
            }
            for metric in ("pds", "signed_credit")
        }
        for scope in ("nested", "fixed")
    }
    return values, report


def paired_control_report(
    primary: Mapping[str, Any],
    control: Mapping[str, Any],
    blocks: np.ndarray,
    draws: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    output = {}
    arrays = {}
    for scope, primary_score, control_score in (
        ("nested", primary["oof"], control["oof"]),
        ("fixed", primary["fixed_score"], control["fixed_score"]),
    ):
        primary_summary = aggregate_scores(
            primary_score, np.arange(DEV_COUNT, dtype=np.int64)
        )
        control_summary = aggregate_scores(
            control_score, np.arange(DEV_COUNT, dtype=np.int64)
        )
        output[scope] = {}
        for metric in ("pds", "signed_credit"):
            difference = (
                primary_summary["per_target"][metric]
                - control_summary["per_target"][metric]
            )
            samples = block_bootstrap(difference, blocks, draws)
            arrays[f"{scope}__{metric}"] = samples
            output[scope][metric] = {
                "mean": float(difference.mean()),
                "lower_95": float(np.quantile(samples, 0.025, method="linear")),
                "upper_95": float(np.quantile(samples, 0.975, method="linear")),
            }
    return output, arrays


def promotion_checks(
    primary: Mapping[str, Any],
    fold_cells: Mapping[str, Mapping[str, Sequence[float]]],
    bootstrap: Mapping[str, Mapping[str, Mapping[str, float]]],
    nulls: Mapping[str, Mapping[str, Mapping[str, Mapping[str, float]]]],
    controls: Mapping[str, Mapping[str, Mapping[str, Mapping[str, float]]]],
    concentration: Mapping[str, Mapping[str, Mapping[str, float]]],
) -> dict[str, bool]:
    checks = {"final_eta_gt_zero": float(primary["final_eta"]) > 0.0}
    for metric in ("pds", "signed_credit"):
        checks[f"nested_fold_family_{metric}_nonnegative"] = all(
            value >= -EPS for family in FAMILIES for value in fold_cells[metric][family]
        )
    for scope, summary in (
        ("nested", primary["oof_summary"]),
        ("fixed", primary["fixed_summary"]),
    ):
        for metric in ("pds", "signed_credit"):
            checks[f"{scope}_{metric}_macro_at_least_0.003"] = (
                summary["macro"][metric] >= MIN_GAIN
            )
            checks[f"{scope}_{metric}_bootstrap_lower_positive"] = (
                bootstrap[scope][metric]["lower_95"] > 0.0
            )
            for kind in ("target_pairing", "donor_assignment"):
                checks[f"{scope}_{kind}_{metric}_above_p95"] = (
                    summary["macro"][metric] > nulls[kind][scope][metric]["p95"]
                )
            checks[f"{scope}_{metric}_family_concentration"] = (
                concentration[scope][metric]["family"] <= MAX_CONCENTRATION
            )
            checks[f"{scope}_{metric}_target_block_concentration"] = (
                concentration[scope][metric]["target_block"] <= MAX_CONCENTRATION
            )
            for control, report in controls.items():
                minimum = MIN_CONTROL_PDS_MARGIN if metric == "pds" else 0.0
                comparator = (
                    report[scope][metric]["mean"] >= minimum
                    if metric == "pds"
                    else report[scope][metric]["mean"] > minimum
                )
                checks[f"{scope}_{control}_{metric}_margin"] = comparator
                checks[f"{scope}_{control}_{metric}_bootstrap_lower_positive"] = (
                    report[scope][metric]["lower_95"] > 0.0
                )
    return checks


def import_experiment132() -> Any:
    """Import the truth loader only after development authorization and locking."""
    spec = importlib.util.spec_from_file_location(
        "experiment132_for_142_develop", EXP132_RUN
    )
    require(
        spec is not None and spec.loader is not None, "cannot import Experiment 132"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_development_truth(limit: float) -> dict[str, np.ndarray]:
    parent = import_experiment132()
    targets = read_lines(DEVELOPMENT_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    truth = {}
    for dataset in SCORED_DATASETS:
        bank = parent.Bank.from_root(dataset, TRUTH_ROOTS[dataset])
        truth[dataset] = np.asarray(
            parent.load_effect_rows(bank, targets, responses), dtype=np.float64
        )
        require(
            truth[dataset].shape == (DEV_COUNT, RESPONSE_COUNT),
            f"{dataset}: truth shape changed",
        )
        memory_check(f"loaded development truth {dataset}", limit)
    return truth


def create_outcome_access_lock(authorization: Mapping[str, Any]) -> None:
    require(
        not OUTCOME_ACCESS_LOCK.exists(), "development outcome access already recorded"
    )
    atomic_json(
        OUTCOME_ACCESS_LOCK,
        {
            "schema_version": 1,
            "created_utc": authorized_provenance_utc(authorization, "develop"),
            "development_targets": DEV_COUNT,
            "sealed_410_effect_values_accessed": False,
            "authorization_sha256": sha256_file(DEVELOPMENT_AUTHORIZATION),
            "authorization_schema": authorization.get("schema_version"),
        },
    )


def develop(limit: float) -> None:
    authorization = validate_development_authorization(verify_truth_files=False)
    created_utc = authorized_provenance_utc(authorization, "develop")
    require(not RESULT_DIR.exists(), f"development result exists: {RESULT_DIR}")
    require(not OUTCOME_ACCESS_LOCK.exists(), "development outcomes already accessed")
    with exclusive_lock(DEVELOP_LOCK):
        create_outcome_access_lock(authorization)
        validate_development_authorization(verify_truth_files=True)
        truth = load_development_truth(limit)
        plan = load_prepared_plan()
        folds = np.asarray(plan["development_folds"], dtype=np.int8)
        blocks = np.asarray(plan["development_blocks"], dtype=np.int16)
        draws = np.asarray(plan["bootstrap_block_draws"], dtype=np.int64)
        geometries = {group: load_prepared_geometry(group) for group in GROUP_SLOTS}

        primary_weights = {
            group: geometries[group]["weights__reliability"] for group in GROUP_SLOTS
        }
        primary_run = score_weight_family(primary_weights, geometries, truth, folds)
        primary = primary_run["selection"]

        controls = {}
        control_bootstrap_arrays = {}
        control_selection = {}
        for control in ANCESTRY_CONTROLS:
            run = score_weight_family(
                {
                    group: geometries[group][f"weights__{control}"]
                    for group in GROUP_SLOTS
                },
                geometries,
                truth,
                folds,
            )
            control_selection[control] = run["selection"]
            report, arrays = paired_control_report(
                primary, run["selection"], blocks, draws
            )
            controls[control] = report
            control_bootstrap_arrays[control] = arrays

        bootstrap = {}
        bootstrap_arrays = {}
        concentration = {}
        for scope, score in (
            ("nested", primary["oof"]),
            ("fixed", primary["fixed_score"]),
        ):
            bootstrap[scope], bootstrap_arrays[scope] = bootstrap_summary(
                score, blocks, draws
            )
            concentration[scope] = concentration_summary(score, blocks)

        null_values = {}
        null_reports = {}
        for kind in ("target_pairing", "donor_assignment"):
            null_values[kind], null_reports[kind] = evaluate_null_family(
                kind, geometries, truth, folds, limit
            )
            for scope in ("nested", "fixed"):
                observed = (
                    primary["oof_summary"]["macro"]
                    if scope == "nested"
                    else primary["fixed_summary"]["macro"]
                )
                for metric in ("pds", "signed_credit"):
                    values = null_values[kind][scope][metric]
                    null_reports[kind][scope][metric]["add_one_p"] = float(
                        (1 + np.sum(values >= observed[metric])) / (len(values) + 1)
                    )

        fold_cells = fold_family_cells(primary["oof"], folds)
        checks = promotion_checks(
            primary,
            fold_cells,
            bootstrap,
            null_reports,
            controls,
            concentration,
        )
        result = {
            "schema_version": 1,
            "stage": "authorized_613_target_development",
            "created_utc": created_utc,
            "selected_eta": primary["final_eta"],
            "nested_folds": primary["folds"],
            "nested_oof": {
                key: value
                for key, value in primary["oof_summary"].items()
                if key != "per_target"
            },
            "fixed": {
                key: value
                for key, value in primary["fixed_summary"].items()
                if key != "per_target"
            },
            "fold_family_cells": fold_cells,
            "bootstrap": bootstrap,
            "nulls": null_reports,
            "ancestry_controls": {
                name: {
                    "selected_eta": selection["final_eta"],
                    "nested_oof": {
                        key: value
                        for key, value in selection["oof_summary"].items()
                        if key != "per_target"
                    },
                    "fixed": {
                        key: value
                        for key, value in selection["fixed_summary"].items()
                        if key != "per_target"
                    },
                    "reliability_minus_control": controls[name],
                }
                for name, selection in control_selection.items()
            },
            "positive_gain_concentration": concentration,
            "checks": checks,
            "passes": all(checks.values()),
            "decision": "PASS_FOR_SEPARATE_REVIEW_ONLY"
            if all(checks.values())
            else "NO-GO",
            "development_613_outcomes_opened": True,
            "sealed_410_effect_values_opened": False,
            "count_generation": False,
            "cell_eval2": False,
            "challenge_build": False,
            "submission": False,
            "peak_rss_gib": peak_rss_gib(),
        }

        with atomic_stage_directory(RESULT_DIR) as staging:
            arrays_to_save = {}
            for dataset in SCORED_DATASETS:
                for metric in ("pds", "signed_credit"):
                    arrays_to_save[f"oof__{dataset}__{metric}"] = primary["oof"][
                        dataset
                    ][metric]
                    arrays_to_save[f"fixed__{dataset}__{metric}"] = primary[
                        "fixed_score"
                    ][dataset][metric]
            for kind, scopes in null_values.items():
                for scope, metrics in scopes.items():
                    for metric, values in metrics.items():
                        arrays_to_save[f"null__{kind}__{scope}__{metric}"] = values
            for scope, metrics in bootstrap_arrays.items():
                for metric, values in metrics.items():
                    arrays_to_save[f"bootstrap__{scope}__{metric}"] = values
            for control, arrays in control_bootstrap_arrays.items():
                for key, values in arrays.items():
                    arrays_to_save[f"control__{control}__{key}"] = values
            atomic_npz(staging / ARRAY_ARTIFACT, **arrays_to_save)
            atomic_json(staging / RESULT_ARTIFACT, result)
            atomic_json(staging / MANIFEST_ARTIFACT, directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    payload: dict[str, Any] = {
        "schema_version": 1,
        "prefreeze_authorization_exists": DEVELOPMENT_PREFREEZE_AUTHORIZATION.exists(),
        "preoutcome_exists": PREOUTCOME_DIR.exists(),
        "development_authorization_exists": DEVELOPMENT_AUTHORIZATION.exists(),
        "outcome_access_recorded": OUTCOME_ACCESS_LOCK.exists(),
        "development_result_exists": RESULT_DIR.exists(),
        "sealed_test_loader_exists": False,
        "truth_arrays_opened_by_audit": False,
    }
    if DEVELOPMENT_PREFREEZE_AUTHORIZATION.exists():
        validate_prefreeze_authorization()
        payload["prefreeze_authorization_valid"] = True
    if PREOUTCOME_DIR.exists():
        validate_preoutcome()
        payload["preoutcome_valid"] = True
    if DEVELOPMENT_AUTHORIZATION.exists():
        validate_development_authorization(verify_truth_files=False)
        payload["development_authorization_valid_without_truth_read"] = True
    if RESULT_DIR.exists():
        validate_directory_manifest(RESULT_DIR)
        result = read_json(RESULT_DIR / RESULT_ARTIFACT)
        require(
            result["sealed_410_effect_values_opened"] is False,
            "result reports sealed-test access",
        )
        payload["development_result"] = {
            "passes": result["passes"],
            "decision": result["decision"],
            "selected_eta": result["selected_eta"],
        }
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "develop", "audit"))
    parser.add_argument("--max-rss-gib", type=float, default=18.0)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.max_rss_gib)
    elif args.command == "develop":
        develop(args.max_rss_gib)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
