#!/usr/bin/env python3
"""Assemble the compact response input and response-free design audit for Exp 107."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

import metadata
import modeling
import source_artifacts


HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CHANNEL_COUNT = 16
TARGET_SUM = 10_000.0
EXPECTED_ORIGINAL_RESPONSE_COUNT = 374
EXPECTED_VCC_MAPPED_RESPONSE_COUNT = 358
EXPECTED_EXPORTABLE_RESPONSE_COUNT = 349
DESIGN_AUDIT_SCHEMA_VERSION = 1
TARGET_SUPPORT_KEYS = frozenset(
    {
        "guide_cells_before_qc",
        "guide_cells_after_qc",
        "target_cells_before_qc",
        "target_cells_after_qc",
        "distinct_called_guides",
        "untreated_guides_with_at_least_three_cells",
        "untreated_support_channels",
        "eligible",
    }
)
DESIGN_AUDIT_KEYS = frozenset(
    {
        "schema_version",
        "held_channel_folds",
        "fold_support_training_channels",
        "fold_support_distinct_called_guides",
        "fold_support_guides_with_three_untreated_cells",
        "fold_support_untreated_carrier_cells",
        "fold_support_training_channel_count",
        "fold_support_supported",
        "fold_zero_variance",
        "fold_nuisance_aliased",
        "fold_fitted",
        "fold_target_scales",
        "fold_nuisance_target_slopes",
        "fold_nullspace_indptr",
        "fold_nullspace_basis",
        "fold_projection_diagonal",
        "fold_singular_indptr",
        "fold_singular_values",
        "fold_retained_singular_directions",
        "fold_effective_rank",
        "fold_condition_number",
        "final_support_training_channels",
        "final_support_distinct_called_guides",
        "final_support_guides_with_three_untreated_cells",
        "final_support_untreated_carrier_cells",
        "final_support_training_channel_count",
        "final_support_supported",
        "final_zero_variance",
        "final_nuisance_aliased",
        "final_fitted",
        "final_target_scales",
        "final_nuisance_target_slopes",
        "final_nullspace_basis",
        "final_projection_diagonal",
        "final_singular_values",
        "final_retained_singular_directions",
        "final_effective_rank",
        "final_condition_number",
    }
)


class ModelInputError(RuntimeError):
    """The frozen source-model input contract was violated."""


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ModelInputError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise ModelInputError(f"{path} must contain a JSON object")
    return value


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def sha256_arrays(**arrays: np.ndarray) -> str:
    """Hash named arrays with dtype, shape, and C-order bytes."""

    digest = hashlib.sha256()
    for name in sorted(arrays):
        value = np.ascontiguousarray(arrays[name])
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(value.dtype.str.encode("ascii"))
        digest.update(b"\0")
        digest.update(json.dumps(value.shape, separators=(",", ":")).encode("ascii"))
        digest.update(b"\0")
        digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def design_audit_paths(input_path: Path) -> tuple[Path, Path]:
    """Return the response-free design audit NPZ and receipt paths."""

    audit_path = input_path.with_name(f"{input_path.stem}-design-audit.npz")
    return audit_path, audit_path.with_suffix(".json")


def _pack_row_matrices(
    matrices: tuple[np.ndarray, ...],
    *,
    columns: int,
) -> tuple[np.ndarray, np.ndarray]:
    indptr = np.zeros(len(matrices) + 1, dtype=np.int64)
    rows: list[np.ndarray] = []
    for index, matrix in enumerate(matrices, start=1):
        value = np.asarray(matrix, dtype=np.float64)
        if value.ndim != 2 or value.shape[1] != columns:
            raise ModelInputError("design-audit row matrix axis changed")
        rows.append(value)
        indptr[index] = indptr[index - 1] + value.shape[0]
    packed = (
        np.concatenate(rows, axis=0)
        if indptr[-1]
        else np.zeros((0, columns), dtype=np.float64)
    )
    return indptr, packed


def _pack_vectors(
    vectors: tuple[np.ndarray, ...],
    *,
    dtype: np.dtype[Any],
) -> tuple[np.ndarray, np.ndarray]:
    indptr = np.zeros(len(vectors) + 1, dtype=np.int64)
    values: list[np.ndarray] = []
    for index, vector in enumerate(vectors, start=1):
        value = np.asarray(vector, dtype=dtype)
        if value.ndim != 1:
            raise ModelInputError("design-audit vector is not one-dimensional")
        values.append(value)
        indptr[index] = indptr[index - 1] + value.size
    packed = np.concatenate(values) if indptr[-1] else np.zeros(0, dtype=dtype)
    return indptr, packed


def _unpack_row_matrices(
    indptr: np.ndarray,
    values: np.ndarray,
    *,
    count: int,
    columns: int,
    name: str,
) -> tuple[np.ndarray, ...]:
    pointers = np.asarray(indptr, dtype=np.int64)
    matrix = np.asarray(values, dtype=np.float64)
    if (
        pointers.shape != (count + 1,)
        or pointers[0] != 0
        or np.any(np.diff(pointers) < 0)
        or matrix.ndim != 2
        or matrix.shape[1] != columns
        or pointers[-1] != matrix.shape[0]
    ):
        raise ModelInputError(f"{name} packed row matrix changed")
    return tuple(
        matrix[start:stop].copy()
        for start, stop in zip(pointers[:-1], pointers[1:], strict=True)
    )


def _unpack_vectors(
    indptr: np.ndarray,
    values: np.ndarray,
    *,
    count: int,
    dtype: np.dtype[Any],
    name: str,
) -> tuple[np.ndarray, ...]:
    pointers = np.asarray(indptr, dtype=np.int64)
    vector = np.asarray(values, dtype=dtype)
    if (
        pointers.shape != (count + 1,)
        or pointers[0] != 0
        or np.any(np.diff(pointers) < 0)
        or vector.ndim != 1
        or pointers[-1] != vector.size
    ):
        raise ModelInputError(f"{name} packed vector changed")
    return tuple(
        vector[start:stop].copy()
        for start, stop in zip(pointers[:-1], pointers[1:], strict=True)
    )


def csr_from_components(
    indptr: np.ndarray,
    indices: np.ndarray,
    *,
    rows: int,
    columns: int,
) -> sparse.csr_matrix:
    """Construct and validate a binary CSR matrix without dense expansion."""

    pointers = np.asarray(indptr, dtype=np.int64)
    column_indices = np.asarray(indices, dtype=np.int64)
    if (
        pointers.shape != (rows + 1,)
        or pointers[0] != 0
        or np.any(np.diff(pointers) < 0)
        or pointers[-1] != column_indices.size
    ):
        raise ModelInputError("CSR pointers do not partition the sparse indices")
    if np.any(column_indices < 0) or np.any(column_indices >= columns):
        raise ModelInputError("CSR column index is outside the matrix")
    for start, stop in zip(pointers[:-1], pointers[1:], strict=True):
        row = column_indices[start:stop]
        if row.size and np.any(np.diff(row) <= 0):
            raise ModelInputError("CSR rows must be canonical and duplicate-free")
    return sparse.csr_matrix(
        (
            np.ones(column_indices.size, dtype=np.float64),
            column_indices,
            pointers,
        ),
        shape=(rows, columns),
    )


def log_cpm_response(
    response_counts: np.ndarray,
    library_size: np.ndarray,
) -> np.ndarray:
    """Apply the frozen log1p(10,000 * count / full-library-size) transform."""

    counts = np.asarray(response_counts)
    totals = np.asarray(library_size, dtype=np.int64)
    if counts.ndim != 2 or totals.shape != (counts.shape[0],):
        raise ModelInputError("response counts and library sizes do not align")
    if np.any(counts < 0) or np.any(totals <= 0):
        raise ModelInputError("response counts and library sizes must be nonnegative")
    if np.any(counts.sum(axis=1, dtype=np.int64) > totals):
        raise ModelInputError("response counts exceed full library size")
    return np.log1p(
        TARGET_SUM * counts.astype(np.float64) / totals.astype(np.float64)[:, None]
    )


def _support_vector(
    values: np.ndarray,
    *,
    name: str,
    target_count: int,
    boolean: bool = False,
) -> np.ndarray:
    array = np.asarray(values)
    if array.shape != (target_count,):
        raise ModelInputError(f"{name} target axis changed")
    if boolean:
        if array.dtype != np.bool_:
            raise ModelInputError(f"{name} must be boolean")
        return array.astype(bool, copy=False)
    try:
        numeric = np.asarray(array, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ModelInputError(f"{name} must contain integer counts") from error
    if (
        not np.all(np.isfinite(numeric))
        or np.any(numeric < 0)
        or np.any(numeric != np.floor(numeric))
    ):
        raise ModelInputError(f"{name} must contain non-negative integer counts")
    return numeric.astype(np.int64)


def load_precomputed_target_support(
    path: Path,
    *,
    target_count: int,
) -> modeling.TargetSupportEvidence:
    """Load the assignment-stage untreated eligibility on the canonical axis."""

    try:
        with np.load(path, allow_pickle=False) as values:
            if set(values.files) != TARGET_SUPPORT_KEYS:
                raise ModelInputError("target-support artifact keys changed")
            arrays = {name: values[name] for name in values.files}
    except (OSError, ValueError) as error:
        raise ModelInputError(
            f"cannot load target-support artifact: {error}"
        ) from error

    target_cells = np.asarray(arrays["target_cells_after_qc"])
    if target_cells.shape != (2, target_count):
        raise ModelInputError("target_cells_after_qc axis changed")
    untreated_carriers = _support_vector(
        target_cells[0],
        name="untreated_carrier_cells",
        target_count=target_count,
    )
    distinct = _support_vector(
        arrays["distinct_called_guides"],
        name="distinct_called_guides",
        target_count=target_count,
    )
    guides_with_three = _support_vector(
        arrays["untreated_guides_with_at_least_three_cells"],
        name="guides_with_three_untreated_cells",
        target_count=target_count,
    )
    channel_count = _support_vector(
        arrays["untreated_support_channels"],
        name="training_channel_count",
        target_count=target_count,
    )
    supported = _support_vector(
        arrays["eligible"],
        name="eligible",
        target_count=target_count,
        boolean=True,
    )
    recomputed = (
        (distinct >= modeling.MIN_DISTINCT_CALLED_GUIDES)
        & (guides_with_three >= modeling.MIN_GUIDES_WITH_THREE_UNTREATED_CELLS)
        & (untreated_carriers >= modeling.MIN_UNTREATED_CARRIER_CELLS)
        & (channel_count >= modeling.MIN_TRAINING_CHANNELS)
    )
    if not np.array_equal(supported, recomputed):
        raise ModelInputError("precomputed untreated eligibility is inconsistent")
    return modeling.TargetSupportEvidence(
        training_channels=np.asarray(modeling.FROZEN_CHANNELS, dtype=np.int64),
        distinct_called_guides=distinct,
        guides_with_three_untreated_cells=guides_with_three,
        untreated_carrier_cells=untreated_carriers,
        training_channel_count=channel_count,
        supported=supported,
    )


def assert_target_support_equal(
    expected: modeling.TargetSupportEvidence,
    observed: modeling.TargetSupportEvidence,
    *,
    label: str,
) -> None:
    """Require exact equality for every support field, including eligibility."""

    for field in (
        "training_channels",
        "distinct_called_guides",
        "guides_with_three_untreated_cells",
        "untreated_carrier_cells",
        "training_channel_count",
        "supported",
    ):
        if not np.array_equal(getattr(expected, field), getattr(observed, field)):
            raise ModelInputError(f"{label} changed at {field}")


def _validate_stage_binding(
    channel: int,
    assignment_report: dict[str, Any],
    response_report: dict[str, Any],
    assignment_path: Path,
    response_path: Path,
    *,
    expected_assignment_report_sha256: str,
    source_verification: dict[str, Any],
) -> None:
    assignment_report_path = WORK / f"channel{channel:02d}" / "assignment-report.json"
    if expected_assignment_report_sha256 != sha256_file(assignment_report_path):
        raise ModelInputError(
            f"channel {channel} assignment changed after assay summarization"
        )
    if assignment_report.get("status") != "complete":
        raise ModelInputError(f"channel {channel} assignment is incomplete")
    if assignment_report.get("output_sha256") != sha256_file(assignment_path):
        raise ModelInputError(f"channel {channel} assignment output changed")
    if assignment_report.get("assignment_stage_sha256") != sha256_file(
        HERE / "assignment_stage.py"
    ):
        raise ModelInputError(f"channel {channel} assignment code changed")
    if response_report.get("status") != "complete":
        raise ModelInputError(f"channel {channel} response extraction is incomplete")
    if response_report.get("assignment_report_sha256") != sha256_file(
        assignment_report_path
    ):
        raise ModelInputError(f"channel {channel} response used stale assignments")
    if response_report.get("assignment_output_sha256") != sha256_file(assignment_path):
        raise ModelInputError(f"channel {channel} response assignment changed")
    if response_report.get("output_sha256") != sha256_file(response_path):
        raise ModelInputError(f"channel {channel} response output changed")
    if response_report.get("assignment_stage_sha256") != sha256_file(
        HERE / "assignment_stage.py"
    ):
        raise ModelInputError(f"channel {channel} response extraction code changed")
    if response_report.get("metadata_sha256") != sha256_file(HERE / "metadata.py"):
        raise ModelInputError(f"channel {channel} response metadata code changed")
    if response_report.get("source_artifacts_sha256") != sha256_file(
        HERE / "source_artifacts.py"
    ):
        raise ModelInputError(
            f"channel {channel} response source-verification code changed"
        )
    if response_report.get("source_raw_archive_sha256") != source_verification.get(
        "raw_archive_sha256"
    ):
        raise ModelInputError(
            f"channel {channel} response used a different source archive"
        )


def _validate_assignment_summary(
    summary: dict[str, Any],
    *,
    source_verification: dict[str, Any],
) -> dict[str, str]:
    required_gates = {
        "retained_cells",
        "untreated_fraction",
        "median_called_guides",
    }
    gates = summary.get("gates")
    if (
        summary.get("status") != "pass"
        or not isinstance(gates, dict)
        or set(gates) != required_gates
        or any(value is not True for value in gates.values())
    ):
        raise ModelInputError("assignment assay-compatibility gate did not pass")
    expected_lineage = {
        "assignment_stage_sha256": sha256_file(HERE / "assignment_stage.py"),
        "source_artifacts_sha256": sha256_file(HERE / "source_artifacts.py"),
        "source_raw_archive_sha256": source_verification.get("raw_archive_sha256"),
        "source_oracle_receipt_sha256": source_verification.get(
            "oracle_receipt_sha256"
        ),
    }
    for field, expected in expected_lineage.items():
        if summary.get(field) != expected:
            raise ModelInputError(f"assignment summary lineage changed at {field}")

    report_hashes = summary.get("channel_report_sha256")
    expected_keys = {str(channel) for channel in range(1, CHANNEL_COUNT + 1)}
    if not isinstance(report_hashes, dict) or set(report_hashes) != expected_keys:
        raise ModelInputError("assignment summary channel hash set changed")
    for channel in range(1, CHANNEL_COUNT + 1):
        path = WORK / f"channel{channel:02d}" / "assignment-report.json"
        if report_hashes[str(channel)] != sha256_file(path):
            raise ModelInputError(
                f"channel {channel} assignment changed after assay summarization"
            )

    target_support_path = WORK / "target-support.npz"
    if summary.get("target_support_sha256") != sha256_file(target_support_path):
        raise ModelInputError("assignment target-support artifact changed")
    metadata_report_path = WORK / "metadata-crosswalk-report.json"
    metadata_output_path = WORK / "metadata-crosswalk.npz"
    if summary.get("metadata_crosswalk_report_sha256") != sha256_file(
        metadata_report_path
    ):
        raise ModelInputError("assignment metadata-crosswalk report changed")
    metadata_report = load_json(metadata_report_path)
    if summary.get("metadata_crosswalk_output_sha256") != metadata_report.get(
        "output_sha256"
    ) or summary.get("metadata_crosswalk_output_sha256") != sha256_file(
        metadata_output_path
    ):
        raise ModelInputError("assignment metadata-crosswalk output changed")
    return {str(key): str(value) for key, value in report_hashes.items()}


def validate_model_input_lineage(
    report: dict[str, Any],
    *,
    source_verification: dict[str, Any],
) -> None:
    """Re-open every upstream report recorded by a source-model input."""

    summary_path = WORK / "assignment-summary.json"
    if report.get("assignment_summary_sha256") != sha256_file(summary_path):
        raise ModelInputError("source-model input assignment summary is stale")
    summary = load_json(summary_path)
    assignment_hashes = _validate_assignment_summary(
        summary,
        source_verification=source_verification,
    )
    target_support_path = WORK / "target-support.npz"
    if report.get("target_support_sha256") != sha256_file(target_support_path):
        raise ModelInputError("source-model input target-support artifact is stale")
    recorded_assignment_hashes = report.get("assignment_report_sha256")
    if recorded_assignment_hashes != assignment_hashes:
        raise ModelInputError("source-model input assignment-report lineage changed")

    response_hashes = report.get("response_report_sha256")
    expected_keys = {str(channel) for channel in range(1, CHANNEL_COUNT + 1)}
    if not isinstance(response_hashes, dict) or set(response_hashes) != expected_keys:
        raise ModelInputError("source-model input response-report hash set changed")
    for channel in range(1, CHANNEL_COUNT + 1):
        directory = WORK / f"channel{channel:02d}"
        assignment_path = directory / "assignments.npz"
        assignment_report_path = directory / "assignment-report.json"
        response_path = directory / "response-counts.npz"
        response_report_path = directory / "response-report.json"
        if response_hashes[str(channel)] != sha256_file(response_report_path):
            raise ModelInputError(
                f"channel {channel} response report changed after model input"
            )
        _validate_stage_binding(
            channel,
            load_json(assignment_report_path),
            load_json(response_report_path),
            assignment_path,
            response_path,
            expected_assignment_report_sha256=assignment_hashes[str(channel)],
            source_verification=source_verification,
        )

    expected_source = {
        "source_raw_archive_sha256": source_verification.get("raw_archive_sha256"),
        "source_oracle_receipt_sha256": source_verification.get(
            "oracle_receipt_sha256"
        ),
        "source_preflight_receipt_sha256": source_verification.get(
            "preflight_receipt_sha256"
        ),
        "numpy_version": source_verification.get("numpy_version"),
        "scipy_version": source_verification.get("scipy_version"),
    }
    for field, expected in expected_source.items():
        if report.get(field) != expected:
            raise ModelInputError(f"source-model input lineage changed at {field}")

    input_relative = report.get("output_path")
    audit_relative = report.get("design_audit_path")
    if not isinstance(input_relative, str) or not isinstance(audit_relative, str):
        raise ModelInputError("source-model input design-audit paths are missing")
    input_path = (HERE / input_relative).resolve()
    audit_path = (HERE / audit_relative).resolve()
    if HERE not in input_path.parents or HERE not in audit_path.parents:
        raise ModelInputError("source-model input design-audit path is unsafe")
    expected_audit_path, audit_report_path = design_audit_paths(input_path)
    if audit_path != expected_audit_path:
        raise ModelInputError("source-model input design-audit path changed")
    if report.get("design_audit_output_sha256") != sha256_file(audit_path):
        raise ModelInputError("source-model input design-audit payload is stale")
    if report.get("design_audit_report_sha256") != sha256_file(audit_report_path):
        raise ModelInputError("source-model input design-audit report is stale")
    audit_report = load_json(audit_report_path)
    if audit_report.get("output_sha256") != report.get(
        "design_audit_output_sha256"
    ) or audit_report.get("design_binding_sha256") != report.get(
        "design_audit_design_binding_sha256"
    ):
        raise ModelInputError("source-model input design-audit binding changed")


def _design_binding_sha256(
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    guide_evidence: modeling.GuideEvidence,
) -> str:
    return sha256_arrays(
        target_indptr=target_design.indptr.astype(np.int64, copy=False),
        target_indices=target_design.indices.astype(np.int32, copy=False),
        nuisance=np.asarray(nuisance, dtype=np.float64),
        channels=np.asarray(channels),
        guide_indptr=guide_evidence.guide_presence.indptr.astype(
            np.int64,
            copy=False,
        ),
        guide_indices=guide_evidence.guide_presence.indices.astype(
            np.int32,
            copy=False,
        ),
        guide_to_target=np.asarray(
            guide_evidence.guide_target_indices,
            dtype=np.int64,
        ),
    )


def write_production_design_audit(
    path: Path,
    audit: modeling.ProductionDesignAudit,
    *,
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    guide_evidence: modeling.GuideEvidence,
) -> dict[str, Any]:
    """Serialize a response-free identifiability audit and hash-bound receipt."""

    path = path.resolve()
    if HERE not in path.parents:
        raise ModelInputError(
            f"refusing to write design audit outside experiment: {path}"
        )
    target_count = target_design.shape[1]
    fold_count = len(audit.fold_geometry)
    if fold_count != len(modeling.FROZEN_HELD_CHANNEL_FOLDS):
        raise ModelInputError("design audit fold count changed")
    fold_nullspace_indptr, fold_nullspace_basis = _pack_row_matrices(
        tuple(item.nullspace_basis for item in audit.fold_geometry),
        columns=target_count,
    )
    fold_singular_indptr, fold_singular_values = _pack_vectors(
        tuple(item.singular_values for item in audit.fold_geometry),
        dtype=np.dtype(np.float64),
    )
    retained_indptr, fold_retained = _pack_vectors(
        tuple(item.retained_singular_directions for item in audit.fold_geometry),
        dtype=np.dtype(np.bool_),
    )
    if not np.array_equal(fold_singular_indptr, retained_indptr):
        raise ModelInputError("design audit singular-value masks changed")

    arrays = {
        "schema_version": np.asarray(DESIGN_AUDIT_SCHEMA_VERSION, dtype=np.int16),
        "held_channel_folds": np.asarray(audit.held_channel_folds, dtype=np.int16),
        "fold_support_training_channels": np.stack(
            [item.training_channels for item in audit.fold_support]
        ),
        "fold_support_distinct_called_guides": np.stack(
            [item.distinct_called_guides for item in audit.fold_support]
        ),
        "fold_support_guides_with_three_untreated_cells": np.stack(
            [item.guides_with_three_untreated_cells for item in audit.fold_support]
        ),
        "fold_support_untreated_carrier_cells": np.stack(
            [item.untreated_carrier_cells for item in audit.fold_support]
        ),
        "fold_support_training_channel_count": np.stack(
            [item.training_channel_count for item in audit.fold_support]
        ),
        "fold_support_supported": np.stack(
            [item.supported for item in audit.fold_support]
        ),
        "fold_zero_variance": np.stack(
            [item.zero_variance for item in audit.fold_geometry]
        ),
        "fold_nuisance_aliased": np.stack(
            [item.nuisance_aliased for item in audit.fold_geometry]
        ),
        "fold_fitted": np.stack([item.fitted for item in audit.fold_geometry]),
        "fold_target_scales": np.stack(
            [item.target_scales for item in audit.fold_geometry]
        ),
        "fold_nuisance_target_slopes": np.stack(
            [item.nuisance_target_slopes for item in audit.fold_geometry]
        ),
        "fold_nullspace_indptr": fold_nullspace_indptr,
        "fold_nullspace_basis": fold_nullspace_basis,
        "fold_projection_diagonal": np.stack(
            [item.projection_diagonal for item in audit.fold_geometry]
        ),
        "fold_singular_indptr": fold_singular_indptr,
        "fold_singular_values": fold_singular_values,
        "fold_retained_singular_directions": fold_retained,
        "fold_effective_rank": np.asarray(
            [item.effective_rank for item in audit.fold_geometry],
            dtype=np.int64,
        ),
        "fold_condition_number": np.asarray(
            [item.condition_number for item in audit.fold_geometry],
            dtype=np.float64,
        ),
        "final_support_training_channels": audit.final_support.training_channels,
        "final_support_distinct_called_guides": (
            audit.final_support.distinct_called_guides
        ),
        "final_support_guides_with_three_untreated_cells": (
            audit.final_support.guides_with_three_untreated_cells
        ),
        "final_support_untreated_carrier_cells": (
            audit.final_support.untreated_carrier_cells
        ),
        "final_support_training_channel_count": (
            audit.final_support.training_channel_count
        ),
        "final_support_supported": audit.final_support.supported,
        "final_zero_variance": audit.final_geometry.zero_variance,
        "final_nuisance_aliased": audit.final_geometry.nuisance_aliased,
        "final_fitted": audit.final_geometry.fitted,
        "final_target_scales": audit.final_geometry.target_scales,
        "final_nuisance_target_slopes": (audit.final_geometry.nuisance_target_slopes),
        "final_nullspace_basis": audit.final_geometry.nullspace_basis,
        "final_projection_diagonal": audit.final_geometry.projection_diagonal,
        "final_singular_values": audit.final_geometry.singular_values,
        "final_retained_singular_directions": (
            audit.final_geometry.retained_singular_directions
        ),
        "final_effective_rank": np.asarray(
            audit.final_geometry.effective_rank,
            dtype=np.int64,
        ),
        "final_condition_number": np.asarray(
            audit.final_geometry.condition_number,
            dtype=np.float64,
        ),
    }
    if set(arrays) != DESIGN_AUDIT_KEYS:
        raise ModelInputError("internal design-audit schema changed")
    atomic_savez(path, **arrays)
    report = {
        "schema_version": DESIGN_AUDIT_SCHEMA_VERSION,
        "status": "complete_before_response_fit",
        "response_values_used": False,
        "effect_semantics": (
            "identifiable_high_moi_conditional_target_contrast_"
            "minimum_norm_row_space_component"
        ),
        "target_count": target_count,
        "nuisance_columns": int(np.asarray(nuisance).shape[1]),
        "held_channel_folds": [list(fold) for fold in audit.held_channel_folds],
        "fold_effective_rank": [item.effective_rank for item in audit.fold_geometry],
        "fold_condition_number": [
            item.condition_number for item in audit.fold_geometry
        ],
        "final_effective_rank": audit.final_geometry.effective_rank,
        "final_condition_number": audit.final_geometry.condition_number,
        "design_binding_sha256": _design_binding_sha256(
            target_design,
            nuisance,
            channels,
            guide_evidence,
        ),
        "modeling_implementation_sha256": sha256_file(HERE / "modeling.py"),
        "model_input_implementation_sha256": sha256_file(Path(__file__)),
        "output_path": str(path.relative_to(HERE)),
        "output_sha256": sha256_file(path),
    }
    atomic_write_json(path.with_suffix(".json"), report)
    return report


def _support_from_arrays(
    arrays: dict[str, np.ndarray],
    *,
    prefix: str,
    index: int | None = None,
) -> modeling.TargetSupportEvidence:
    def value(name: str) -> np.ndarray:
        current = np.asarray(arrays[f"{prefix}{name}"])
        return current if index is None else current[index]

    return modeling.TargetSupportEvidence(
        training_channels=np.asarray(
            value("training_channels"),
            dtype=np.int64,
        ),
        distinct_called_guides=np.asarray(
            value("distinct_called_guides"),
            dtype=np.int64,
        ),
        guides_with_three_untreated_cells=np.asarray(
            value("guides_with_three_untreated_cells"),
            dtype=np.int64,
        ),
        untreated_carrier_cells=np.asarray(
            value("untreated_carrier_cells"),
            dtype=np.int64,
        ),
        training_channel_count=np.asarray(
            value("training_channel_count"),
            dtype=np.int64,
        ),
        supported=np.asarray(value("supported"), dtype=bool),
    )


def load_production_design_audit(
    path: Path,
    *,
    target_design: sparse.csr_matrix,
    nuisance: np.ndarray,
    channels: np.ndarray,
    guide_evidence: modeling.GuideEvidence,
) -> modeling.ProductionDesignAudit:
    """Load a design-only audit and bind it to the exact current model inputs."""

    path = path.resolve()
    report = load_json(path.with_suffix(".json"))
    expected_report = {
        "schema_version": DESIGN_AUDIT_SCHEMA_VERSION,
        "status": "complete_before_response_fit",
        "response_values_used": False,
        "modeling_implementation_sha256": sha256_file(HERE / "modeling.py"),
        "model_input_implementation_sha256": sha256_file(Path(__file__)),
        "output_path": str(path.relative_to(HERE)),
        "output_sha256": sha256_file(path),
        "design_binding_sha256": _design_binding_sha256(
            target_design,
            nuisance,
            channels,
            guide_evidence,
        ),
    }
    for field, expected in expected_report.items():
        if report.get(field) != expected:
            raise ModelInputError(f"design audit changed at {field}")
    try:
        with np.load(path, allow_pickle=False) as values:
            if set(values.files) != DESIGN_AUDIT_KEYS:
                raise ModelInputError("design-audit artifact keys changed")
            arrays = {name: values[name] for name in values.files}
    except (OSError, ValueError) as error:
        raise ModelInputError(f"cannot load design-audit artifact: {error}") from error

    if int(np.asarray(arrays["schema_version"]).reshape(())) != (
        DESIGN_AUDIT_SCHEMA_VERSION
    ):
        raise ModelInputError("design-audit schema version changed")
    folds = tuple(
        tuple(int(value) for value in row)
        for row in np.asarray(arrays["held_channel_folds"]).tolist()
    )
    fold_count = len(modeling.FROZEN_HELD_CHANNEL_FOLDS)
    target_count = target_design.shape[1]
    nullspaces = _unpack_row_matrices(
        arrays["fold_nullspace_indptr"],
        arrays["fold_nullspace_basis"],
        count=fold_count,
        columns=target_count,
        name="fold nullspace",
    )
    singular_values = _unpack_vectors(
        arrays["fold_singular_indptr"],
        arrays["fold_singular_values"],
        count=fold_count,
        dtype=np.dtype(np.float64),
        name="fold singular values",
    )
    retained = _unpack_vectors(
        arrays["fold_singular_indptr"],
        arrays["fold_retained_singular_directions"],
        count=fold_count,
        dtype=np.dtype(np.bool_),
        name="fold retained singular directions",
    )
    fold_support = tuple(
        _support_from_arrays(arrays, prefix="fold_support_", index=index)
        for index in range(fold_count)
    )
    fold_geometry = tuple(
        modeling.DesignGeometryAudit(
            training_channels=np.asarray(
                arrays["fold_support_training_channels"][index],
                dtype=np.int64,
            ),
            evidence_support=fold_support[index].supported,
            zero_variance=np.asarray(
                arrays["fold_zero_variance"][index],
                dtype=bool,
            ),
            nuisance_aliased=np.asarray(
                arrays["fold_nuisance_aliased"][index],
                dtype=bool,
            ),
            fitted=np.asarray(arrays["fold_fitted"][index], dtype=bool),
            target_scales=np.asarray(
                arrays["fold_target_scales"][index],
                dtype=np.float64,
            ),
            nuisance_target_slopes=np.asarray(
                arrays["fold_nuisance_target_slopes"][index],
                dtype=np.float64,
            ),
            nullspace_basis=nullspaces[index],
            projection_diagonal=np.asarray(
                arrays["fold_projection_diagonal"][index],
                dtype=np.float64,
            ),
            singular_values=singular_values[index],
            retained_singular_directions=retained[index],
            effective_rank=int(arrays["fold_effective_rank"][index]),
            condition_number=float(arrays["fold_condition_number"][index]),
        )
        for index in range(fold_count)
    )
    final_support = _support_from_arrays(arrays, prefix="final_support_")
    final_geometry = modeling.DesignGeometryAudit(
        training_channels=final_support.training_channels,
        evidence_support=final_support.supported,
        zero_variance=np.asarray(arrays["final_zero_variance"], dtype=bool),
        nuisance_aliased=np.asarray(arrays["final_nuisance_aliased"], dtype=bool),
        fitted=np.asarray(arrays["final_fitted"], dtype=bool),
        target_scales=np.asarray(arrays["final_target_scales"], dtype=np.float64),
        nuisance_target_slopes=np.asarray(
            arrays["final_nuisance_target_slopes"],
            dtype=np.float64,
        ),
        nullspace_basis=np.asarray(
            arrays["final_nullspace_basis"],
            dtype=np.float64,
        ),
        projection_diagonal=np.asarray(
            arrays["final_projection_diagonal"],
            dtype=np.float64,
        ),
        singular_values=np.asarray(
            arrays["final_singular_values"],
            dtype=np.float64,
        ),
        retained_singular_directions=np.asarray(
            arrays["final_retained_singular_directions"],
            dtype=bool,
        ),
        effective_rank=int(np.asarray(arrays["final_effective_rank"]).reshape(())),
        condition_number=float(
            np.asarray(arrays["final_condition_number"]).reshape(())
        ),
    )
    return modeling.ProductionDesignAudit(
        held_channel_folds=folds,
        fold_support=fold_support,
        fold_geometry=fold_geometry,
        final_support=final_support,
        final_geometry=final_geometry,
    )


def build_model_input(
    output_path: Path = WORK / "source-model-input.npz",
) -> dict[str, Any]:
    """Assemble all untreated cells without opening any held outcome."""

    output_path = output_path.resolve()
    if HERE not in output_path.parents:
        raise ModelInputError(
            f"refusing to write outside the experiment directory: {output_path}"
        )
    summary_path = WORK / "assignment-summary.json"
    summary = load_json(summary_path)
    source_verification = source_artifacts.verify_source_artifacts()
    assignment_report_hashes = _validate_assignment_summary(
        summary,
        source_verification=source_verification,
    )

    crosswalks = metadata.load_pinned_crosswalks()
    target_count = len(crosswalks.guide.target_symbols)
    guide_count = len(crosswalks.guide.guide_ids)
    original_response_symbols = np.asarray(
        [record.symbol for record in crosswalks.response.records],
        dtype="U",
    )
    original_response_vcc_indices = np.asarray(
        [
            -1 if record.vcc_gene_index is None else record.vcc_gene_index
            for record in crosswalks.response.records
        ],
        dtype=np.int32,
    )
    original_response_exportable = np.asarray(
        [
            record.vcc_gene_index is not None and not record.is_panel_target
            for record in crosswalks.response.records
        ],
        dtype=bool,
    )
    response_original_indices = np.flatnonzero(original_response_exportable).astype(
        np.int32
    )
    if original_response_symbols.shape != (EXPECTED_ORIGINAL_RESPONSE_COUNT,):
        raise ModelInputError(
            "original response axis is not exactly "
            f"{EXPECTED_ORIGINAL_RESPONSE_COUNT} genes"
        )
    if (
        np.count_nonzero(original_response_vcc_indices >= 0)
        != EXPECTED_VCC_MAPPED_RESPONSE_COUNT
    ):
        raise ModelInputError(
            "original response-to-VCC coverage is not exactly "
            f"{EXPECTED_VCC_MAPPED_RESPONSE_COUNT} genes"
        )
    if response_original_indices.size != EXPECTED_EXPORTABLE_RESPONSE_COUNT:
        raise ModelInputError(
            "modeled response axis is not exactly "
            f"{EXPECTED_EXPORTABLE_RESPONSE_COUNT} genes"
        )
    response_symbols = original_response_symbols[response_original_indices]
    response_vcc_indices = original_response_vcc_indices[response_original_indices]
    response_exportable = np.ones(
        EXPECTED_EXPORTABLE_RESPONSE_COUNT,
        dtype=bool,
    )
    target_support_path = WORK / "target-support.npz"
    precomputed_support = load_precomputed_target_support(
        target_support_path,
        target_count=target_count,
    )
    responses: list[np.ndarray] = []
    nuisances: list[np.ndarray] = []
    channels: list[np.ndarray] = []
    barcodes: list[np.ndarray] = []
    target_designs: list[sparse.csr_matrix] = []
    guide_designs: list[sparse.csr_matrix] = []
    input_hashes: dict[str, str] = {}

    for channel in range(1, CHANNEL_COUNT + 1):
        directory = WORK / f"channel{channel:02d}"
        assignment_path = directory / "assignments.npz"
        assignment_report_path = directory / "assignment-report.json"
        response_path = directory / "response-counts.npz"
        response_report_path = directory / "response-report.json"
        assignment_report = load_json(assignment_report_path)
        response_report = load_json(response_report_path)
        _validate_stage_binding(
            channel,
            assignment_report,
            response_report,
            assignment_path,
            response_path,
            expected_assignment_report_sha256=assignment_report_hashes[str(channel)],
            source_verification=source_verification,
        )

        with np.load(assignment_path, allow_pickle=False) as assignment:
            cell_indices = assignment["cell_indices"].astype(np.int64, copy=False)
            condition_codes = assignment["condition_codes"].astype(
                np.uint8,
                copy=False,
            )
            library_size = assignment["library_size"].astype(np.int64, copy=False)
            called_burden = assignment["called_guide_burden"].astype(
                np.float64,
                copy=False,
            )
            ntc_burden = assignment["non_targeting_guide_burden"].astype(
                np.float64,
                copy=False,
            )
            current_barcodes = assignment["barcodes"].astype(str)
            targets = csr_from_components(
                assignment["target_indptr"],
                assignment["target_indices"],
                rows=cell_indices.size,
                columns=target_count,
            )
            guides = csr_from_components(
                assignment["called_guide_indptr"],
                assignment["called_guide_indices"],
                rows=cell_indices.size,
                columns=guide_count,
            )
        with np.load(response_path, allow_pickle=False) as response:
            response_cells = response["cell_indices"].astype(np.int64, copy=False)
            response_counts = response["response_counts"]
            response_library = response["library_size"].astype(
                np.int64,
                copy=False,
            )
        if (
            not np.array_equal(response_cells, cell_indices)
            or not np.array_equal(response_library, library_size)
            or response_counts.shape
            != (cell_indices.size, EXPECTED_ORIGINAL_RESPONSE_COUNT)
        ):
            raise ModelInputError(
                f"channel {channel} assignment and response axes disagree"
            )

        untreated = condition_codes == 0
        if not np.any(untreated):
            raise ModelInputError(f"channel {channel} has no untreated cells")
        responses.append(
            log_cpm_response(
                response_counts[untreated][:, response_original_indices],
                library_size[untreated],
            ).astype(np.float32)
        )
        nuisances.append(
            np.column_stack(
                (
                    np.log1p(library_size[untreated].astype(np.float64)),
                    ntc_burden[untreated],
                    called_burden[untreated],
                )
            )
        )
        channels.append(np.full(np.count_nonzero(untreated), channel, dtype=np.int16))
        barcodes.append(current_barcodes[untreated])
        target_designs.append(targets[untreated])
        guide_designs.append(guides[untreated])
        input_hashes[str(channel)] = sha256_file(response_report_path)

    response = np.concatenate(responses, axis=0)
    nuisance = np.concatenate(nuisances, axis=0)
    channel_labels = np.concatenate(channels)
    barcode_axis = np.concatenate(barcodes)
    target_design = sparse.vstack(target_designs, format="csr")
    guide_design = sparse.vstack(guide_designs, format="csr")
    if (
        response.shape[0] != nuisance.shape[0]
        or response.shape[0] != channel_labels.size
        or response.shape[0] != barcode_axis.size
        or response.shape[0] != target_design.shape[0]
        or response.shape[0] != guide_design.shape[0]
    ):
        raise ModelInputError("assembled source-model row axes disagree")
    if not np.all(np.isfinite(response)) or not np.all(np.isfinite(nuisance)):
        raise ModelInputError("assembled source-model arrays are non-finite")
    observed_support = modeling.derive_target_support(
        modeling.GuideEvidence(
            guide_presence=guide_design,
            guide_target_indices=crosswalks.guide.guide_to_target,
        ),
        target_design,
        channel_labels,
        training_channels=modeling.FROZEN_CHANNELS,
    )
    assert_target_support_equal(
        precomputed_support,
        observed_support,
        label="assembled untreated target support",
    )
    guide_evidence = modeling.GuideEvidence(
        guide_presence=guide_design,
        guide_target_indices=crosswalks.guide.guide_to_target,
    )
    design_audit = modeling.audit_production_design(
        target_design,
        nuisance,
        channel_labels,
        guide_evidence,
    )
    assert_target_support_equal(
        precomputed_support,
        design_audit.final_support,
        label="design-audit all-channel untreated target support",
    )
    design_audit_path, design_audit_report_path = design_audit_paths(output_path)
    design_audit_report = write_production_design_audit(
        design_audit_path,
        design_audit,
        target_design=target_design,
        nuisance=nuisance,
        channels=channel_labels,
        guide_evidence=guide_evidence,
    )
    atomic_savez(
        output_path,
        response=response,
        nuisance=nuisance,
        channels=channel_labels,
        barcodes=barcode_axis,
        target_indptr=target_design.indptr.astype(np.int64),
        target_indices=target_design.indices.astype(np.int32),
        guide_indptr=guide_design.indptr.astype(np.int64),
        guide_indices=guide_design.indices.astype(np.int32),
        target_symbols=np.asarray(crosswalks.guide.target_symbols, dtype="U"),
        guide_ids=np.asarray(crosswalks.guide.guide_ids, dtype="U"),
        guide_to_target=crosswalks.guide.guide_to_target,
        response_symbols=response_symbols,
        response_vcc_indices=response_vcc_indices,
        response_exportable=response_exportable,
        original_response_symbols=original_response_symbols,
        original_response_vcc_indices=original_response_vcc_indices,
        original_response_exportable=original_response_exportable,
        response_original_indices=response_original_indices,
        precomputed_training_channels=precomputed_support.training_channels,
        precomputed_distinct_called_guides=(precomputed_support.distinct_called_guides),
        precomputed_guides_with_three_untreated_cells=(
            precomputed_support.guides_with_three_untreated_cells
        ),
        precomputed_untreated_carrier_cells=(
            precomputed_support.untreated_carrier_cells
        ),
        precomputed_training_channel_count=(precomputed_support.training_channel_count),
        precomputed_eligible=precomputed_support.supported,
    )
    report = {
        "status": "complete",
        "primary_condition": "untreated",
        "cells": int(response.shape[0]),
        "channels": CHANNEL_COUNT,
        "response_genes": int(response.shape[1]),
        "modeled_response_genes": int(response.shape[1]),
        "exportable_response_genes": int(response.shape[1]),
        "original_response_genes": int(original_response_symbols.size),
        "original_vcc_mapped_response_genes": int(
            np.count_nonzero(original_response_vcc_indices >= 0)
        ),
        "excluded_response_genes": int(
            original_response_symbols.size - response.shape[1]
        ),
        "original_response_axis_sha256": sha256_arrays(
            symbols=original_response_symbols,
            vcc_indices=original_response_vcc_indices,
            exportable=original_response_exportable,
        ),
        "modeled_response_axis_sha256": sha256_arrays(
            original_indices=response_original_indices,
            symbols=response_symbols,
            vcc_indices=response_vcc_indices,
        ),
        "targets": target_count,
        "guides": guide_count,
        "target_design_nnz": int(target_design.nnz),
        "guide_design_nnz": int(guide_design.nnz),
        "response_dtype": str(response.dtype),
        "nuisance_columns": [
            "log1p_full_transcriptome_library_size",
            "called_non_targeting_guide_count",
            "total_called_guide_count",
        ],
        "assignment_summary_sha256": sha256_file(summary_path),
        "assignment_report_sha256": assignment_report_hashes,
        "target_support_sha256": sha256_file(target_support_path),
        "source_raw_archive_sha256": source_verification["raw_archive_sha256"],
        "source_oracle_receipt_sha256": source_verification["oracle_receipt_sha256"],
        "source_preflight_receipt_sha256": source_verification[
            "preflight_receipt_sha256"
        ],
        "numpy_version": source_verification["numpy_version"],
        "scipy_version": source_verification["scipy_version"],
        "response_report_sha256": input_hashes,
        "design_audit_path": str(design_audit_path.relative_to(HERE)),
        "design_audit_report_sha256": sha256_file(design_audit_report_path),
        "design_audit_output_sha256": design_audit_report["output_sha256"],
        "design_audit_design_binding_sha256": design_audit_report[
            "design_binding_sha256"
        ],
        "implementation_sha256": sha256_file(Path(__file__)),
        "output_path": str(output_path.relative_to(HERE)),
        "output_sha256": sha256_file(output_path),
        "nadig_or_challenge_outcome_opened": False,
    }
    atomic_write_json(output_path.with_suffix(".json"), report)
    return report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument(
        "--output",
        type=Path,
        default=WORK / "source-model-input.npz",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        report = build_model_input(args.output.resolve())
    except (
        ModelInputError,
        source_artifacts.SourceArtifactError,
        OSError,
        ValueError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
