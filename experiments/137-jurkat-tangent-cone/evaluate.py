#!/usr/bin/env python3
"""Leakage-guarded source-only evaluation for Experiment 137."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
AUTHORIZATION = HERE / "PREVALUE_AUTHORIZATION.v1.json"
METADATA_NPZ = HERE / "work/metadata-freeze-v1.npz"
METADATA_JSON = HERE / "work/metadata-freeze-v1.json"
CANDIDATE_NPZ = HERE / "work/source-candidate-v1/candidate.npz"
CANDIDATE_JSON = HERE / "work/source-candidate-v1/candidate.json"
TANGENT_PATH = HERE / "tangent.py"
PRETRUTH_PLAN = HERE / "work/pretruth-plan-v1.json"
RUNTIME_AMENDMENT = HERE / "RUNTIME_AMENDMENT.v1.json"
DEVELOPMENT_LOCK = HERE / "work/development-outcome-access.lock"
DEVELOPMENT_DIR = HERE / "work/development-result-v1"
DEVELOPMENT_NPZ = DEVELOPMENT_DIR / "evaluation-arrays.npz"
DEVELOPMENT_JSON = DEVELOPMENT_DIR / "evaluation-summary.json"
CANDIDATE_FREEZE = HERE / "SOURCE_CANDIDATE_FREEZE.v1.json"
TEST_LOCK = HERE / "work/test-outcome-access.lock"
TEST_DIR = HERE / "work/test-result-v1"
TEST_NPZ = TEST_DIR / "evaluation-arrays.npz"
TEST_JSON = TEST_DIR / "evaluation-summary.json"
NADIG_DELTA = ROOT / "data/donors/nadig_jurkat/delta.npy"

EXPECTED_PROTOCOL_SHA256 = (
    "f1f085296213cfcbfab93212cc963982f34b683e5ce6e7945edf5bee65ea65f7"
)
EXPECTED_METADATA_NPZ_SHA256 = (
    "682afb9a32a2f2732486da04cf27359abf2d4b91c1ba92627b8d05f45dce9d1c"
)
EXPECTED_METADATA_JSON_SHA256 = (
    "45de70d1e9c687d8ab8967946f9c7f46a4982d41b6c107f9753905a0fe0ba1a0"
)
EXPECTED_NADIG_DELTA_SHA256 = (
    "1fa16024c2d7115fc4b7c06f0ab7d9afe8ba6ebd0ae6dcef62c4955636e29129"
)
EXPECTED_NADIG_DELTA_BYTES = 80_639_608
PERMUTATIONS = 99
BOOTSTRAPS = 20_000
MINIMUM_PDS_GAIN = 0.003
MAXIMUM_CLUSTER_FRACTION = 0.5


class EvaluationError(RuntimeError):
    """Raised when an evaluation boundary or frozen invariant changes."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise EvaluationError(message)


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    require(not array.dtype.hasobject, "object arrays cannot be hashed")
    digest = hashlib.sha256()
    digest.update(b"exp137-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(array.shape).encode("ascii"))
    digest.update(b"\0")
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == encoded, f"refusing overwrite: {path}")
        return
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".npz"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        np.savez_compressed(temporary, **arrays)
        if path.exists():
            require(path.read_bytes() == temporary.read_bytes(), f"refusing overwrite: {path}")
            temporary.unlink()
        else:
            os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def create_or_validate_lock(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if path.exists():
        require(path.read_bytes() == encoded, f"access lock changed: {path}")
        return
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o444)
    try:
        os.write(descriptor, encoded)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def cap_key(cap: float) -> str:
    return f"cap_{int(round(float(cap) * 100)):04d}"


def fix_unsupported_pds(
    scores: np.ndarray, transferable: np.ndarray
) -> np.ndarray:
    values = np.asarray(scores, dtype=np.float64).copy()
    supported = np.asarray(transferable)
    require(values.ndim == 1, "PDS scores must be one-dimensional")
    require(supported.dtype == np.bool_, "transferability mask must be boolean")
    require(supported.shape == values.shape, "transferability mask shape changed")
    require(np.all(np.isfinite(values)), "PDS scores are non-finite")
    values[~supported] = 0.5
    return values


def validate_authorization() -> dict[str, Any]:
    require(AUTHORIZATION.is_file(), "pre-value authorization is missing")
    payload = read_json(AUTHORIZATION)
    require(payload.get("status") == "approved_before_effect_access", "not authorized")
    bindings = payload.get("bindings")
    require(isinstance(bindings, dict), "authorization bindings are missing")
    for role, path in {
        "metadata_npz": METADATA_NPZ,
        "metadata_json": METADATA_JSON,
        "protocol": PROTOCOL,
        "builder": HERE / "build_source_candidate.py",
        "evaluator": Path(__file__).resolve(),
        "tangent": TANGENT_PATH,
        "builder_test": HERE / "test_build_source_candidate.py",
        "evaluator_test": HERE / "test_evaluate.py",
        "tangent_test": HERE / "test_tangent.py",
        "metadata_test": HERE / "test_freeze_metadata.py",
    }.items():
        require(role in bindings, f"authorization lacks {role}")
        require(
            bindings[role].get("path") == str(path.relative_to(ROOT)),
            f"authorization {role} path changed",
        )
        require(
            bindings[role].get("sha256") == sha256_file(path),
            f"authorization {role} hash changed",
        )
    return payload


def load_frozen_arrays() -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    require(
        sha256_file(PROTOCOL) == EXPECTED_PROTOCOL_SHA256, "protocol changed"
    )
    require(
        sha256_file(METADATA_NPZ) == EXPECTED_METADATA_NPZ_SHA256,
        "metadata NPZ changed",
    )
    require(
        sha256_file(METADATA_JSON) == EXPECTED_METADATA_JSON_SHA256,
        "metadata JSON changed",
    )
    candidate_report = read_json(CANDIDATE_JSON)
    require(
        candidate_report.get("status") == "source_candidate_built_without_nadig_truth",
        "candidate status changed",
    )
    require(
        candidate_report.get("artifact", {}).get("sha256") == sha256_file(CANDIDATE_NPZ),
        "candidate artifact hash changed",
    )
    require(candidate_report.get("nadig_truth_opened") is False, "candidate saw truth")
    with np.load(METADATA_NPZ, allow_pickle=False) as archive:
        metadata = {name: archive[name] for name in archive.files}
    with np.load(CANDIDATE_NPZ, allow_pickle=False) as archive:
        candidate = {name: archive[name] for name in archive.files}
    require(
        candidate_report.get("arrays")
        == {name: sha256_array(value) for name, value in candidate.items()},
        "candidate array bindings changed",
    )
    bindings = candidate_report.get("bindings")
    require(isinstance(bindings, dict), "candidate bindings are missing")
    for key, path in {
        "authorization_sha256": AUTHORIZATION,
        "metadata_npz_sha256": METADATA_NPZ,
        "metadata_json_sha256": METADATA_JSON,
        "protocol_sha256": PROTOCOL,
        "builder_sha256": HERE / "build_source_candidate.py",
        "tangent_sha256": TANGENT_PATH,
        "experiment_107_builder_sha256": (
            ROOT / "experiments/107-jurkat-high-moi/strict_artifact_builder.py"
        ),
        "experiment_107_strict_sha256": (
            ROOT / "experiments/107-jurkat-high-moi/strict_nadig_validation.py"
        ),
        "source_model_sha256": (
            ROOT
            / "experiments/107-jurkat-high-moi/work/frozen-source-model/source-model.npz"
        ),
        "source_manifest_sha256": (
            ROOT
            / "experiments/107-jurkat-high-moi/work/frozen-source-model/manifest.json"
        ),
        "operator_json_sha256": (
            HERE / "work/source-candidate-v1/operator.json"
        ),
    }.items():
        require(
            bindings.get(key) == sha256_file(path),
            f"candidate binding changed: {key}",
        )
    production_bindings = bindings.get("production_implementation_sha256")
    require(
        isinstance(production_bindings, dict) and production_bindings,
        "production implementation bindings are missing",
    )
    for relative, expected_hash in production_bindings.items():
        path = ROOT / relative
        require(
            path.is_file() and sha256_file(path) == expected_hash,
            f"production implementation changed: {relative}",
        )
    operator = candidate_report.get("operator")
    require(isinstance(operator, dict), "candidate operator binding is missing")
    operator_path = ROOT / str(operator.get("path", ""))
    require(
        operator_path.is_file()
        and operator.get("bytes") == operator_path.stat().st_size
        and operator.get("sha256") == sha256_file(operator_path),
        "candidate operator artifact changed",
    )
    transitive = candidate_report.get("transitive_inputs")
    require(
        isinstance(transitive, list) and transitive,
        "candidate transitive inputs are missing",
    )
    seen: set[str] = set()
    for record in transitive:
        require(isinstance(record, dict), "invalid transitive input record")
        relative = str(record.get("path", ""))
        require(relative and relative not in seen, "duplicate transitive input path")
        seen.add(relative)
        path = ROOT / relative
        require(
            path.is_file()
            and record.get("bytes") == path.stat().st_size
            and record.get("sha256") == sha256_file(path),
            f"candidate transitive input changed: {relative}",
        )
    require(
        np.array_equal(candidate["targets"], metadata["fresh_targets"]),
        "candidate target axis changed",
    )
    require(
        np.array_equal(candidate["score_genes"], metadata["score_genes"]),
        "candidate score axis changed",
    )
    return metadata, candidate


def pretruth_payload() -> dict[str, Any]:
    validate_authorization()
    metadata, candidate = load_frozen_arrays()
    require(NADIG_DELTA.stat().st_size == EXPECTED_NADIG_DELTA_BYTES, "truth size changed")
    require(
        sha256_file(NADIG_DELTA) == EXPECTED_NADIG_DELTA_SHA256,
        "truth file hash changed",
    )
    return {
        "schema_version": 1,
        "experiment": "137-jurkat-tangent-cone",
        "status": "frozen_before_experiment_137_fresh_target_outcome_access",
        "bindings": {
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "protocol_sha256": sha256_file(PROTOCOL),
            "metadata_npz_sha256": sha256_file(METADATA_NPZ),
            "metadata_json_sha256": sha256_file(METADATA_JSON),
            "candidate_npz_sha256": sha256_file(CANDIDATE_NPZ),
            "candidate_json_sha256": sha256_file(CANDIDATE_JSON),
            "evaluator_sha256": sha256_file(Path(__file__)),
            "tangent_sha256": sha256_file(TANGENT_PATH),
        },
        "truth": {
            "path": str(NADIG_DELTA.relative_to(ROOT)),
            "bytes": NADIG_DELTA.stat().st_size,
            "sha256": sha256_file(NADIG_DELTA),
            "array_values_opened": False,
        },
        "axes": {
            "fresh_targets_sha256_array": sha256_array(metadata["fresh_targets"]),
            "test_targets_sha256_array": sha256_array(metadata["test_targets"]),
            "development_targets_sha256_array": sha256_array(
                metadata["development_targets"]
            ),
            "development_fold_id_sha256_array": sha256_array(
                metadata["development_fold_id"]
            ),
            "score_genes_sha256_array": sha256_array(metadata["score_genes"]),
            "fresh_nadig_indices_sha256_array": sha256_array(
                metadata["fresh_nadig_indices"]
            ),
        },
        "candidate_arrays": {
            name: sha256_array(value)
            for name, value in candidate.items()
        },
    }


def require_pretruth_plan() -> dict[str, Any]:
    require(PRETRUTH_PLAN.is_file(), "pretruth plan is missing")
    observed = read_json(PRETRUTH_PLAN)
    expected = pretruth_payload()
    if observed != expected:
        authorization = validate_authorization()
        review = authorization.get("review")
        require(isinstance(review, dict), "authorization review is missing")
        binding = review.get("runtime_amendment")
        require(isinstance(binding, dict), "runtime amendment is not authorized")
        require(
            binding.get("path") == str(RUNTIME_AMENDMENT.relative_to(ROOT)),
            "runtime amendment path changed",
        )
        require(
            binding.get("sha256") == sha256_file(RUNTIME_AMENDMENT),
            "runtime amendment hash changed",
        )
        amendment = read_json(RUNTIME_AMENDMENT)
        require(
            amendment.get("status")
            == "approved_numerically_blind_runtime_repair_after_development_open",
            "runtime amendment status changed",
        )
        require(
            amendment.get("pretruth_plan_sha256") == sha256_file(PRETRUTH_PLAN),
            "runtime amendment pretruth binding changed",
        )
        require(
            amendment.get("development_lock_sha256")
            == sha256_file(DEVELOPMENT_LOCK),
            "runtime amendment development lock binding changed",
        )
        require(
            amendment.get("old_evaluator_sha256")
            == observed.get("bindings", {}).get("evaluator_sha256"),
            "runtime amendment old evaluator binding changed",
        )
        require(
            amendment.get("new_evaluator_sha256")
            == sha256_file(Path(__file__).resolve()),
            "runtime amendment new evaluator binding changed",
        )
        require(
            amendment.get("development_result_observed") is False
            and amendment.get("development_artifact_written") is False
            and amendment.get("test_rows_accessed") is False,
            "runtime amendment outcome-state claim changed",
        )
        require(
            {key: value for key, value in observed.items() if key != "bindings"}
            == {key: value for key, value in expected.items() if key != "bindings"},
            "pretruth scientific payload changed",
        )
        observed_bindings = observed.get("bindings")
        expected_bindings = expected.get("bindings")
        require(
            isinstance(observed_bindings, dict)
            and isinstance(expected_bindings, dict),
            "pretruth bindings are missing",
        )
        mutable_provenance = {
            "authorization_sha256",
            "candidate_json_sha256",
            "evaluator_sha256",
        }
        require(
            set(observed_bindings) == set(expected_bindings),
            "pretruth binding keys changed",
        )
        for key in set(expected_bindings) - mutable_provenance:
            require(
                observed_bindings[key] == expected_bindings[key],
                f"pretruth scientific binding changed: {key}",
            )
    return observed


def derangement(rng: np.random.Generator, size: int) -> np.ndarray:
    base = np.arange(size, dtype=np.int64)
    for _ in range(10_000):
        trial = rng.permutation(base)
        if np.all(trial != base):
            return trial
    raise EvaluationError("could not produce a derangement")


def reciprocal_knn_clusters(
    incumbent: np.ndarray,
    *,
    fold_one_based: int,
    k: int = 5,
    minimum_cosine: float = 0.80,
) -> np.ndarray:
    values = np.asarray(incumbent, dtype=np.float64)
    require(values.ndim == 2 and values.shape[0] >= 2, "cluster matrix changed")
    norms = np.linalg.norm(values, axis=1)
    unit = np.divide(
        values, norms[:, None], out=np.zeros_like(values), where=norms[:, None] > 0
    )
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    count = min(max(1, int(k)), len(values) - 1)
    neighbors = np.argsort(-similarity, axis=1, kind="stable")[:, :count]
    directed = np.zeros((len(values), len(values)), dtype=bool)
    directed[np.arange(len(values))[:, None], neighbors] = True
    reciprocal = directed & directed.T & (similarity >= minimum_cosine)
    parent = np.arange(len(values), dtype=np.int64)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    for left, right in zip(*np.where(np.triu(reciprocal, k=1)), strict=True):
        left_root, right_root = find(int(left)), find(int(right))
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)
    roots = [find(index) for index in range(len(values))]
    labels = {root: number for number, root in enumerate(sorted(set(roots)))}
    return np.asarray(
        [f"f{fold_one_based}:c{labels[root]}" for root in roots], dtype="U32"
    )


def development_scores(
    candidate: Mapping[str, np.ndarray],
    truth: np.ndarray,
    tangent: ModuleType,
) -> tuple[
    dict[float, tuple[np.ndarray, ...]],
    dict[float, tuple[np.ndarray, ...]],
    tuple[np.ndarray, ...],
]:
    fold_id = np.asarray(candidate["development_fold_id"], dtype=np.int8)
    require(fold_id.shape == (315,), "development fold axis changed")
    development_offset = 79
    pds: dict[float, list[np.ndarray]] = {cap: [] for cap in tangent.CAP_DEGREES}
    credit: dict[float, list[np.ndarray]] = {
        cap: [] for cap in tangent.CAP_DEGREES
    }
    clusters: list[np.ndarray] = []
    for fold in range(1, 6):
        local = np.flatnonzero(fold_id == fold)
        require(local.size == 63, f"development fold {fold} size changed")
        rows = development_offset + local
        fold_truth = truth[local]
        fold_transferable = np.asarray(candidate["transferable"], dtype=bool)[rows]
        clusters.append(
            reciprocal_knn_clusters(
                np.asarray(candidate["incumbent_score"])[rows],
                fold_one_based=fold,
            )
        )
        for cap in tangent.CAP_DEGREES:
            values = np.asarray(candidate[f"candidate_{cap_key(cap)}"])[rows]
            pds[cap].append(
                fix_unsupported_pds(
                    tangent.pds_per_target(values, fold_truth),
                    fold_transferable,
                )
            )
            credit[cap].append(tangent.signed_credit_per_target(values, fold_truth))
    return (
        {cap: tuple(values) for cap, values in pds.items()},
        {cap: tuple(values) for cap, values in credit.items()},
        tuple(clusters),
    )


def nested_statistics(
    pds: Mapping[float, Sequence[np.ndarray]],
    credit: Mapping[float, Sequence[np.ndarray]],
    tangent: ModuleType,
) -> dict[str, Any]:
    selected_caps, selected_pds = tangent.select_cap(dict(pds))
    baseline_pds = np.asarray(
        [np.asarray(values).mean() for values in pds[0.0]], dtype=np.float64
    )
    selected_credit = np.asarray(
        [
            np.asarray(credit[cap][fold]).mean()
            for fold, cap in enumerate(selected_caps)
        ],
        dtype=np.float64,
    )
    baseline_credit = np.asarray(
        [np.asarray(values).mean() for values in credit[0.0]], dtype=np.float64
    )
    final = tangent.final_cap(dict(pds))
    return {
        "selected_caps_degrees": list(selected_caps),
        "held_pds": selected_pds,
        "baseline_pds": baseline_pds,
        "held_pds_delta": selected_pds - baseline_pds,
        "macro_pds_delta": float(np.mean(selected_pds - baseline_pds)),
        "held_signed_credit": selected_credit,
        "baseline_signed_credit": baseline_credit,
        "held_signed_credit_delta": selected_credit - baseline_credit,
        "macro_signed_credit_delta": float(
            np.mean(selected_credit - baseline_credit)
        ),
        "final_cap_degrees": float(final),
    }


def nested_permutation_null(
    candidate: Mapping[str, np.ndarray],
    truth: np.ndarray,
    tangent: ModuleType,
) -> np.ndarray:
    fold_id = np.asarray(candidate["development_fold_id"], dtype=np.int8)
    source = np.asarray(candidate["source_score"], dtype=np.float64)
    incumbent = np.asarray(candidate["incumbent_score"], dtype=np.float64)
    transferable = np.asarray(candidate["transferable"], dtype=bool)
    coordinate_supported = np.asarray(candidate["coordinate_supported"], dtype=bool)
    draws = np.empty(PERMUTATIONS, dtype=np.float64)
    for replicate in range(PERMUTATIONS):
        pds: dict[float, list[np.ndarray]] = {
            cap: [] for cap in tangent.CAP_DEGREES
        }
        for fold in range(1, 6):
            local = np.flatnonzero(fold_id == fold)
            rows = 79 + local
            fold_source = source[rows].copy()
            fold_transferable = transferable[rows]
            movable = np.flatnonzero(fold_transferable)
            rng = np.random.Generator(
                np.random.PCG64(
                    np.random.SeedSequence([1_372_026, 137, replicate, fold])
                )
            )
            permutation = derangement(rng, len(movable))
            fold_source[movable] = fold_source[movable[permutation]]
            for cap in tangent.CAP_DEGREES:
                changed, _, _ = tangent.spherical_tangent_candidate(
                    incumbent[rows],
                    fold_source,
                    fold_transferable,
                    coordinate_supported,
                    cap,
                )
                pds[cap].append(
                    fix_unsupported_pds(
                        tangent.pds_per_target(changed, truth[local]),
                        fold_transferable,
                    )
                )
        selected_caps, held = tangent.select_cap(
            {cap: tuple(values) for cap, values in pds.items()}
        )
        del selected_caps
        baseline = np.asarray(
            [np.asarray(values).mean() for values in pds[0.0]], dtype=np.float64
        )
        draws[replicate] = float(np.mean(held - baseline))
    return draws


def cluster_bootstrap(
    pds: Mapping[float, Sequence[np.ndarray]],
    clusters: Sequence[np.ndarray],
    tangent: ModuleType,
) -> np.ndarray:
    draws = np.empty(BOOTSTRAPS, dtype=np.float64)
    for replicate in range(BOOTSTRAPS):
        means = np.empty((len(tangent.CAP_DEGREES), 5), dtype=np.float64)
        for fold in range(5):
            labels = np.asarray(clusters[fold]).astype(str)
            unique = np.unique(labels)
            rng = np.random.Generator(
                np.random.PCG64(
                    np.random.SeedSequence([1_372_026, 137, 20_000, replicate, fold + 1])
                )
            )
            sampled = rng.choice(unique, size=len(unique), replace=True)
            multiplicity = {
                label: int(np.count_nonzero(sampled == label)) for label in unique
            }
            weights = np.asarray(
                [multiplicity[label] for label in labels], dtype=np.int64
            )
            denominator = int(weights.sum())
            require(denominator > 0, "empty bootstrap draw")
            for cap_index, cap in enumerate(tangent.CAP_DEGREES):
                means[cap_index, fold] = float(
                    np.dot(np.asarray(pds[cap][fold]), weights) / denominator
                )
        held = np.empty(5, dtype=np.float64)
        for fold in range(5):
            objective = np.delete(means, fold, axis=1).mean(axis=1)
            chosen = int(np.flatnonzero(objective == objective.max())[0])
            held[fold] = means[chosen, fold]
        draws[replicate] = float(np.mean(held - means[0]))
    return draws


def gain_concentration(
    deltas: Sequence[np.ndarray], clusters: Sequence[np.ndarray], macro: float
) -> float:
    positive_total = float(
        math.fsum(max(float(value), 0.0) for array in deltas for value in array)
    )
    if macro <= 0 or positive_total <= 0:
        return math.inf
    sums: dict[str, float] = {}
    for fold, (values, labels) in enumerate(zip(deltas, clusters, strict=True), start=1):
        for value, label in zip(values, labels, strict=True):
            key = f"{fold}:{label}"
            sums[key] = sums.get(key, 0.0) + max(float(value), 0.0)
    return max(sums.values(), default=0.0) / positive_total


def development_stage() -> dict[str, Any]:
    require_pretruth_plan()
    metadata, candidate = load_frozen_arrays()
    create_or_validate_lock(
        DEVELOPMENT_LOCK,
        {
            "schema_version": 1,
            "experiment": "137-jurkat-tangent-cone",
            "stage": "development",
            "pretruth_plan_sha256": sha256_file(PRETRUTH_PLAN),
            "truth_sha256": sha256_file(NADIG_DELTA),
            "test_rows_accessed": False,
        },
    )
    truth_memmap = np.load(NADIG_DELTA, mmap_mode="r", allow_pickle=False)
    development_rows = np.asarray(metadata["fresh_nadig_indices"], dtype=np.int64)[79:]
    score_columns = np.asarray(metadata["score_indices_in_nadig"], dtype=np.int64)
    truth = np.asarray(
        truth_memmap[np.ix_(development_rows, score_columns)], dtype=np.float64
    )
    del truth_memmap

    tangent = load_module("exp137_tangent_development", TANGENT_PATH)
    pds, credit, clusters = development_scores(candidate, truth, tangent)
    nested = nested_statistics(pds, credit, tangent)
    null = nested_permutation_null(candidate, truth, tangent)
    bootstrap = cluster_bootstrap(pds, clusters, tangent)
    per_target_deltas = tuple(
        np.asarray(pds[cap][fold]) - np.asarray(pds[0.0][fold])
        for fold, cap in enumerate(nested["selected_caps_degrees"])
    )
    concentration = gain_concentration(
        per_target_deltas, clusters, nested["macro_pds_delta"]
    )
    maximum_angle = max(
        float(
            np.asarray(
                candidate[f"realized_degrees_{cap_key(cap)}"], dtype=np.float64
            ).max(initial=0.0)
        )
        for cap in tangent.CAP_DEGREES
    )
    incumbent = np.asarray(candidate["incumbent_score"])
    transferable = np.asarray(candidate["transferable"], dtype=bool)
    coordinate_supported = np.asarray(candidate["coordinate_supported"], dtype=bool)
    exact_fallback_by_cap: dict[str, bool] = {}
    for cap in tangent.CAP_DEGREES:
        key = cap_key(cap)
        changed = np.asarray(candidate[f"candidate_{key}"])
        effective = np.asarray(candidate[f"effective_{key}"], dtype=bool)
        tangent.assert_exact_fallback(
            incumbent, changed, effective, coordinate_supported
        )
        exact_fallback_by_cap[key] = bool(
            np.array_equal(
                np.ascontiguousarray(changed[~transferable]).view(np.uint8),
                np.ascontiguousarray(incumbent[~transferable]).view(np.uint8),
            )
            and np.array_equal(
                np.ascontiguousarray(changed[:, ~coordinate_supported]).view(
                    np.uint8
                ),
                np.ascontiguousarray(incumbent[:, ~coordinate_supported]).view(
                    np.uint8
                ),
            )
        )
    checks = {
        "selected_cap_nonzero": nested["final_cap_degrees"] > 0.0,
        "nested_pds_delta_each_fold_nonnegative": bool(
            np.all(nested["held_pds_delta"] >= 0.0)
        ),
        "nested_macro_pds_delta_at_least_0_003": (
            nested["macro_pds_delta"] >= MINIMUM_PDS_GAIN
        ),
        "nested_signed_credit_delta_each_fold_nonnegative": bool(
            np.all(nested["held_signed_credit_delta"] >= 0.0)
        ),
        "nested_macro_signed_credit_delta_nonnegative": (
            nested["macro_signed_credit_delta"] >= 0.0
        ),
        "bootstrap_lower_95_strictly_positive": (
            float(np.quantile(bootstrap, 0.025, method="linear")) > 0.0
        ),
        "observed_strictly_exceeds_all_permutations": (
            nested["macro_pds_delta"] > float(null.max())
        ),
        "cluster_concentration_at_most_0_5": (
            concentration <= MAXIMUM_CLUSTER_FRACTION
        ),
        "maximum_realized_angle_at_most_6_degrees": maximum_angle <= 6.000000001,
        "exact_identity": bool(
            np.array_equal(
                np.ascontiguousarray(candidate["candidate_cap_0000"]).view(np.uint8),
                np.ascontiguousarray(candidate["incumbent_score"]).view(np.uint8),
            )
        ),
        "exact_fallback_every_cap": all(exact_fallback_by_cap.values()),
    }
    passes = all(checks.values())
    arrays: dict[str, np.ndarray] = {
        "truth_development": truth,
        "selected_caps_degrees": np.asarray(
            nested["selected_caps_degrees"], dtype=np.float64
        ),
        "held_pds": np.asarray(nested["held_pds"], dtype=np.float64),
        "baseline_pds": np.asarray(nested["baseline_pds"], dtype=np.float64),
        "held_pds_delta": np.asarray(nested["held_pds_delta"], dtype=np.float64),
        "held_signed_credit": np.asarray(
            nested["held_signed_credit"], dtype=np.float64
        ),
        "baseline_signed_credit": np.asarray(
            nested["baseline_signed_credit"], dtype=np.float64
        ),
        "held_signed_credit_delta": np.asarray(
            nested["held_signed_credit_delta"], dtype=np.float64
        ),
        "permutation_null": null,
        "cluster_bootstrap": bootstrap,
    }
    for fold in range(5):
        arrays[f"fold_{fold + 1}_clusters"] = np.asarray(clusters[fold])
        arrays[f"fold_{fold + 1}_selected_pds_delta"] = per_target_deltas[fold]
    atomic_savez(DEVELOPMENT_NPZ, **arrays)
    summary = {
        "schema_version": 1,
        "experiment": "137-jurkat-tangent-cone",
        "status": "development_complete",
        "promotion_passes": passes,
        "checks": checks,
        "nested": {
            "selected_caps_degrees": nested["selected_caps_degrees"],
            "held_pds_delta": nested["held_pds_delta"].tolist(),
            "macro_pds_delta": nested["macro_pds_delta"],
            "held_signed_credit_delta": nested[
                "held_signed_credit_delta"
            ].tolist(),
            "macro_signed_credit_delta": nested["macro_signed_credit_delta"],
            "final_cap_degrees": nested["final_cap_degrees"],
        },
        "bootstrap": {
            "replicates": BOOTSTRAPS,
            "lower_95": float(np.quantile(bootstrap, 0.025, method="linear")),
            "upper_95": float(np.quantile(bootstrap, 0.975, method="linear")),
        },
        "permutation_null": {
            "replicates": PERMUTATIONS,
            "maximum": float(null.max()),
            "p95_higher": float(np.quantile(null, 0.95, method="higher")),
        },
        "gain_concentration_fraction": concentration,
        "maximum_realized_angle_degrees": maximum_angle,
        "exact_fallback_by_cap": exact_fallback_by_cap,
        "artifact": {
            "path": str(DEVELOPMENT_NPZ.relative_to(ROOT)),
            "bytes": DEVELOPMENT_NPZ.stat().st_size,
            "sha256": sha256_file(DEVELOPMENT_NPZ),
        },
        "test_outcome_opened": False,
    }
    atomic_write_json(DEVELOPMENT_JSON, summary)
    if passes:
        freeze = {
            "schema_version": 1,
            "experiment": "137-jurkat-tangent-cone",
            "status": "nonidentity_candidate_frozen_before_test_outcome_access",
            "selected_cap_degrees": nested["final_cap_degrees"],
            "bindings": {
                "development_summary_sha256": sha256_file(DEVELOPMENT_JSON),
                "development_arrays_sha256": sha256_file(DEVELOPMENT_NPZ),
                "candidate_npz_sha256": sha256_file(CANDIDATE_NPZ),
                "candidate_json_sha256": sha256_file(CANDIDATE_JSON),
                "pretruth_plan_sha256": sha256_file(PRETRUTH_PLAN),
                "protocol_sha256": sha256_file(PROTOCOL),
                "evaluator_sha256": sha256_file(Path(__file__)),
                "tangent_sha256": sha256_file(TANGENT_PATH),
            },
        }
        atomic_write_json(CANDIDATE_FREEZE, freeze)
    return summary


def validated_test_cap(
    freeze: Mapping[str, Any],
    development: Mapping[str, Any],
    candidate: Mapping[str, np.ndarray],
    tangent: ModuleType,
) -> float:
    selected = freeze.get("selected_cap_degrees")
    require(isinstance(selected, (int, float)), "frozen cap is not numeric")
    cap = float(selected)
    require(math.isfinite(cap), "frozen cap is non-finite")
    require(cap in tangent.CAP_DEGREES, "frozen cap is outside the fixed grid")
    require(cap > 0.0, "identity cannot enter the test")
    nested = development.get("nested")
    require(isinstance(nested, dict), "development nested summary is missing")
    development_selected = nested.get("final_cap_degrees")
    require(
        isinstance(development_selected, (int, float)),
        "development final cap is not numeric",
    )
    development_cap = float(development_selected)
    require(
        math.isfinite(development_cap)
        and development_cap in tangent.CAP_DEGREES,
        "development final cap is outside the fixed grid",
    )
    require(cap == development_cap, "frozen cap differs from development selection")
    key = f"candidate_{cap_key(cap)}"
    require(key in candidate, "frozen candidate array is missing")
    require(
        np.asarray(candidate[key]).shape
        == np.asarray(candidate["incumbent_score"]).shape,
        "frozen candidate array shape changed",
    )
    return cap


def test_stage() -> dict[str, Any]:
    require(CANDIDATE_FREEZE.is_file(), "passing candidate freeze is missing")
    freeze = read_json(CANDIDATE_FREEZE)
    require(
        freeze.get("status")
        == "nonidentity_candidate_frozen_before_test_outcome_access",
        "candidate freeze status changed",
    )
    development = read_json(DEVELOPMENT_JSON)
    require(development.get("promotion_passes") is True, "development did not pass")
    metadata, candidate = load_frozen_arrays()
    require(
        freeze["bindings"]["development_summary_sha256"]
        == sha256_file(DEVELOPMENT_JSON),
        "development summary changed",
    )
    require(
        freeze["bindings"]["development_arrays_sha256"]
        == sha256_file(DEVELOPMENT_NPZ),
        "development arrays changed",
    )
    for key, path in {
        "candidate_npz_sha256": CANDIDATE_NPZ,
        "candidate_json_sha256": CANDIDATE_JSON,
        "pretruth_plan_sha256": PRETRUTH_PLAN,
        "protocol_sha256": PROTOCOL,
        "evaluator_sha256": Path(__file__),
        "tangent_sha256": TANGENT_PATH,
    }.items():
        require(
            freeze["bindings"].get(key) == sha256_file(path),
            f"candidate freeze binding changed: {key}",
        )
    tangent = load_module("exp137_tangent_test", TANGENT_PATH)
    cap = validated_test_cap(freeze, development, candidate, tangent)
    create_or_validate_lock(
        TEST_LOCK,
        {
            "schema_version": 1,
            "experiment": "137-jurkat-tangent-cone",
            "stage": "untouched_test",
            "candidate_freeze_sha256": sha256_file(CANDIDATE_FREEZE),
            "truth_sha256": sha256_file(NADIG_DELTA),
        },
    )
    truth_memmap = np.load(NADIG_DELTA, mmap_mode="r", allow_pickle=False)
    test_rows = np.asarray(metadata["fresh_nadig_indices"], dtype=np.int64)[:79]
    score_columns = np.asarray(metadata["score_indices_in_nadig"], dtype=np.int64)
    truth = np.asarray(
        truth_memmap[np.ix_(test_rows, score_columns)], dtype=np.float64
    )
    del truth_memmap
    incumbent = np.asarray(candidate["incumbent_score"], dtype=np.float64)[:79]
    selected = np.asarray(candidate[f"candidate_{cap_key(cap)}"], dtype=np.float64)[
        :79
    ]
    baseline_pds = tangent.pds_per_target(incumbent, truth)
    selected_pds = tangent.pds_per_target(selected, truth)
    baseline_credit = tangent.signed_credit_per_target(incumbent, truth)
    selected_credit = tangent.signed_credit_per_target(selected, truth)
    pds_delta = selected_pds - baseline_pds
    credit_delta = selected_credit - baseline_credit
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([1_372_026, 137, 79])))
    bootstrap = np.empty(BOOTSTRAPS, dtype=np.float64)
    for replicate in range(BOOTSTRAPS):
        indices = rng.integers(0, 79, size=79)
        bootstrap[replicate] = float(pds_delta[indices].mean())
    null = np.empty(PERMUTATIONS, dtype=np.float64)
    source = np.asarray(candidate["source_score"], dtype=np.float64)[:79]
    transferable = np.asarray(candidate["transferable"], dtype=bool)[:79]
    coordinate_supported = np.asarray(candidate["coordinate_supported"], dtype=bool)
    require(np.all(transferable), "the frozen 79-target test lost transferability")
    for replicate in range(PERMUTATIONS):
        random = np.random.Generator(
            np.random.PCG64(np.random.SeedSequence([1_372_026, 137, 79, replicate]))
        )
        permutation = derangement(random, 79)
        permuted, _, _ = tangent.spherical_tangent_candidate(
            incumbent,
            source[permutation],
            transferable,
            coordinate_supported,
            cap,
        )
        null[replicate] = float(
            (
                tangent.pds_per_target(permuted, truth)
                - baseline_pds
            ).mean()
        )
    checks = {
        "pds_delta_at_least_0_003": float(pds_delta.mean()) >= MINIMUM_PDS_GAIN,
        "signed_credit_delta_nonnegative": float(credit_delta.mean()) >= 0.0,
        "bootstrap_lower_95_strictly_positive": float(
            np.quantile(bootstrap, 0.025, method="linear")
        )
        > 0.0,
        "observed_strictly_exceeds_all_permutations": float(pds_delta.mean())
        > float(null.max()),
        "exact_unsupported_coordinate_fallback": bool(
            np.array_equal(
                np.ascontiguousarray(selected[:, ~coordinate_supported]).view(
                    np.uint8
                ),
                np.ascontiguousarray(incumbent[:, ~coordinate_supported]).view(
                    np.uint8
                ),
            )
        ),
    }
    arrays = {
        "truth_test": truth,
        "baseline_pds": baseline_pds,
        "selected_pds": selected_pds,
        "pds_delta": pds_delta,
        "baseline_signed_credit": baseline_credit,
        "selected_signed_credit": selected_credit,
        "signed_credit_delta": credit_delta,
        "bootstrap": bootstrap,
        "permutation_null": null,
    }
    atomic_savez(TEST_NPZ, **arrays)
    summary = {
        "schema_version": 1,
        "experiment": "137-jurkat-tangent-cone",
        "status": "untouched_test_complete",
        "promotion_passes": all(checks.values()),
        "checks": checks,
        "selected_cap_degrees": cap,
        "pds_delta": float(pds_delta.mean()),
        "signed_credit_delta": float(credit_delta.mean()),
        "bootstrap": {
            "replicates": BOOTSTRAPS,
            "lower_95": float(np.quantile(bootstrap, 0.025, method="linear")),
            "upper_95": float(np.quantile(bootstrap, 0.975, method="linear")),
        },
        "permutation_null": {
            "replicates": PERMUTATIONS,
            "maximum": float(null.max()),
            "p95_higher": float(np.quantile(null, 0.95, method="higher")),
        },
        "artifact": {
            "path": str(TEST_NPZ.relative_to(ROOT)),
            "bytes": TEST_NPZ.stat().st_size,
            "sha256": sha256_file(TEST_NPZ),
        },
        "challenge_action_authorized": False,
    }
    atomic_write_json(TEST_JSON, summary)
    return summary


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    subparsers = result.add_subparsers(dest="command", required=True)
    subparsers.add_parser("pretruth")
    subparsers.add_parser("development")
    subparsers.add_parser("test")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "pretruth":
        payload = pretruth_payload()
        atomic_write_json(PRETRUTH_PLAN, payload)
        print(json.dumps({"status": payload["status"], "sha256": sha256_file(PRETRUTH_PLAN)}))
        return 0
    if args.command == "development":
        print(json.dumps(development_stage(), sort_keys=True))
        return 0
    if args.command == "test":
        print(json.dumps(test_stage(), sort_keys=True))
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
