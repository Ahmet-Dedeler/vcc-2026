#!/usr/bin/env python3
"""Streaming assignment and response extraction for Experiment 107.

The source matrices have 6,794,880 columns per channel. This module only
materializes expression-called cells, then retained hash singlets, and finally
the frozen 374-gene response panel. Full transcriptome matrices and dense
83,401-guide matrices are never constructed.
"""

from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
import os
import sys
import tarfile
import tempfile
from collections import Counter
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np

import demultiplex
import guide_calling
import metadata
import preprocess
import source_artifacts


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "downloads" / "GSE249595_RAW.tar"
SCHEMA_REPORT = HERE / "audit" / "schema.json"
WORK = HERE / "work"
CHANNEL_COUNT = 16
EXPECTED_COLUMNS = 6_794_880
EXPECTED_TRANSCRIPTOME_ROWS = 20_606
EXPECTED_GUIDE_ROWS = 83_401
EXPECTED_LABEL_ROWS = 7
MINIMUM_DETECTED_GENES = 10
MAXIMUM_MITOCHONDRIAL_FRACTION = 0.30
MAXIMUM_ABSOLUTE_LOG_LIBRARY_MAD = 5.0
RETAINED_CELL_RANGE = (400_000, 700_000)
UNTREATED_FRACTION_RANGE = (0.05, 0.15)
MEDIAN_CALLED_GUIDE_RANGE = (8.0, 18.0)


class AssignmentStageError(RuntimeError):
    """A frozen assignment-stage invariant failed."""


@dataclass(frozen=True)
class ChannelAssignment:
    """Compact final-cell state for one channel."""

    cell_indices: np.ndarray
    barcodes: np.ndarray
    hash_label_indices: np.ndarray
    condition_codes: np.ndarray
    library_size: np.ndarray
    detected_genes: np.ndarray
    mitochondrial_counts: np.ndarray
    total_guide_umi: np.ndarray
    called_guide_burden: np.ndarray
    non_targeting_guide_burden: np.ndarray
    called_guide_indptr: np.ndarray
    called_guide_indices: np.ndarray
    target_indptr: np.ndarray
    target_indices: np.ndarray


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
        raise AssignmentStageError(f"cannot read {path}: {error}") from error
    if not isinstance(value, dict):
        raise AssignmentStageError(f"{path} must contain a JSON object")
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


@contextmanager
def exclusive_lock() -> Iterator[None]:
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / ".assignment-stage.lock"
    with path.open("a+b") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise AssignmentStageError(
                "another experiment-107 assignment stage is active"
            ) from error
        yield


@contextmanager
def gzip_member(
    archive: tarfile.TarFile,
    name: str,
) -> Iterator[BinaryIO]:
    extracted = archive.extractfile(name)
    if extracted is None:
        raise AssignmentStageError(f"archive member is missing: {name}")
    try:
        with gzip.GzipFile(fileobj=extracted, mode="rb") as stream:
            yield stream
    finally:
        extracted.close()


def member_name(
    accession: str,
    channel: int,
    modality: str,
    role: str,
) -> str:
    suffix = ".mtx.gz" if role == "matrix" else ".tsv.gz"
    return f"{accession}_channel{channel}_{modality}_{role}{suffix}"


def _validated_selected_columns(
    selected_columns: np.ndarray,
    column_count: int,
) -> np.ndarray:
    selected = np.asarray(selected_columns)
    if selected.ndim != 1 or not np.issubdtype(selected.dtype, np.integer):
        raise AssignmentStageError("selected columns must be an integer vector")
    selected = selected.astype(np.int64, copy=False)
    if np.any(selected < 0) or np.any(selected >= column_count):
        raise AssignmentStageError("selected column is outside the matrix")
    if selected.size and np.any(np.diff(selected) <= 0):
        raise AssignmentStageError(
            "selected columns must be unique and strictly increasing"
        )
    return selected


def iter_column_major_coordinates(
    stream: BinaryIO,
    *,
    expected_rows: int,
    expected_columns: int,
) -> Iterator[tuple[int, int, int]]:
    """Yield a validated strict column-major Matrix Market coordinate stream."""

    dimensions, coordinates = preprocess._iter_matrix_market_coordinates(stream)
    rows, columns, _ = dimensions
    if rows != expected_rows or columns != expected_columns:
        raise AssignmentStageError(
            f"matrix shape changed: {(rows, columns)} != "
            f"{(expected_rows, expected_columns)}"
        )
    previous = (-1, -1)
    for row, column, value in coordinates:
        key = (column, row)
        if key <= previous:
            raise AssignmentStageError(
                "matrix coordinates are not strict column-major order"
            )
        previous = key
        yield row, column, value


def extract_selected_dense(
    stream: BinaryIO,
    selected_columns: np.ndarray,
    *,
    expected_rows: int,
    expected_columns: int,
    dtype: np.dtype[Any] = np.int64,
) -> np.ndarray:
    """Extract selected columns into a cells-by-features dense matrix."""

    selected = _validated_selected_columns(selected_columns, expected_columns)
    output = np.zeros((selected.size, expected_rows), dtype=dtype)
    selected_position = 0
    for row, column, value in iter_column_major_coordinates(
        stream,
        expected_rows=expected_rows,
        expected_columns=expected_columns,
    ):
        while (
            selected_position < selected.size and selected[selected_position] < column
        ):
            selected_position += 1
        if selected_position < selected.size and selected[selected_position] == column:
            output[selected_position, row] = value
    return output


def _make_sparse_batch(
    cell_indices: list[int],
    indptr: list[int],
    feature_indices: list[int],
    values: list[int],
) -> guide_calling.SparseGuideBatch:
    return guide_calling.SparseGuideBatch(
        cell_indices=np.asarray(cell_indices, dtype=np.int64),
        indptr=np.asarray(indptr, dtype=np.int64),
        guide_indices=np.asarray(feature_indices, dtype=np.int64),
        observed_umi=np.asarray(values, dtype=np.int64),
    )


def iter_selected_sparse_batches(
    stream: BinaryIO,
    selected_columns: np.ndarray,
    *,
    expected_rows: int,
    expected_columns: int,
    batch_cells: int = 2_048,
) -> Iterator[guide_calling.SparseGuideBatch]:
    """Stream selected columns as CSR-like batches, including empty columns."""

    selected = _validated_selected_columns(selected_columns, expected_columns)
    if not isinstance(batch_cells, int) or isinstance(batch_cells, bool):
        raise AssignmentStageError("batch_cells must be a positive integer")
    if batch_cells <= 0:
        raise AssignmentStageError("batch_cells must be a positive integer")
    if selected.size == 0:
        return

    batch_indices: list[int] = []
    indptr = [0]
    feature_indices: list[int] = []
    values: list[int] = []
    selected_position = 0

    def append_selected() -> None:
        batch_indices.append(int(selected[selected_position]))
        indptr.append(len(feature_indices))

    for row, column, value in iter_column_major_coordinates(
        stream,
        expected_rows=expected_rows,
        expected_columns=expected_columns,
    ):
        while (
            selected_position < selected.size and selected[selected_position] < column
        ):
            append_selected()
            selected_position += 1
            if len(batch_indices) == batch_cells:
                yield _make_sparse_batch(
                    batch_indices,
                    indptr,
                    feature_indices,
                    values,
                )
                batch_indices = []
                indptr = [0]
                feature_indices = []
                values = []
        if selected_position == selected.size:
            continue
        if selected[selected_position] == column:
            feature_indices.append(row)
            values.append(value)

    while selected_position < selected.size:
        append_selected()
        selected_position += 1
        if len(batch_indices) == batch_cells:
            yield _make_sparse_batch(
                batch_indices,
                indptr,
                feature_indices,
                values,
            )
            batch_indices = []
            indptr = [0]
            feature_indices = []
            values = []
    if batch_indices:
        yield _make_sparse_batch(batch_indices, indptr, feature_indices, values)


def extract_selected_barcodes(
    stream: BinaryIO,
    selected_columns: np.ndarray,
    *,
    expected_columns: int,
) -> np.ndarray:
    """Read only selected barcode rows while validating the full row count."""

    selected = _validated_selected_columns(selected_columns, expected_columns)
    output: list[str] = []
    selected_position = 0
    observed_rows = 0
    for observed_rows, raw in enumerate(stream, start=1):
        column = observed_rows - 1
        if selected_position < selected.size and selected[selected_position] == column:
            fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
            if len(fields) != 1 or not fields[0]:
                raise AssignmentStageError("barcode row must contain one value")
            output.append(fields[0])
            selected_position += 1
    if observed_rows != expected_columns:
        raise AssignmentStageError(
            f"barcode row count changed: {observed_rows} != {expected_columns}"
        )
    if selected_position != selected.size:
        raise AssignmentStageError("one or more selected barcodes were not found")
    if len(set(output)) != len(output):
        raise AssignmentStageError("selected barcode identifiers are not unique")
    return np.asarray(output, dtype=f"<U{max(map(len, output), default=1)}")


def expression_qc_mask(
    library_size: np.ndarray,
    detected_genes: np.ndarray,
    mitochondrial_counts: np.ndarray,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Apply the frozen expression-QC thresholds within one hash-singlet channel."""

    totals = np.asarray(library_size, dtype=np.int64)
    detected = np.asarray(detected_genes, dtype=np.int64)
    mitochondrial = np.asarray(mitochondrial_counts, dtype=np.int64)
    if (
        totals.ndim != 1
        or detected.shape != totals.shape
        or mitochondrial.shape != totals.shape
    ):
        raise AssignmentStageError("expression-QC vectors must have matching shapes")
    if totals.size == 0:
        raise AssignmentStageError("expression QC has no hash-singlet cells")
    if np.any(totals <= 0) or np.any(detected < 0) or np.any(mitochondrial < 0):
        raise AssignmentStageError("expression-QC inputs are outside valid ranges")
    if np.any(mitochondrial > totals):
        raise AssignmentStageError("mitochondrial counts exceed library size")

    log_library = np.log1p(totals.astype(np.float64))
    center = float(np.median(log_library))
    raw_mad = float(np.median(np.abs(log_library - center)))
    scaled_mad = 1.4826 * raw_mad
    if not np.isfinite(scaled_mad) or scaled_mad <= 0:
        raise AssignmentStageError("within-channel log-library MAD is undefined")
    absolute_mad = np.abs(log_library - center) / scaled_mad
    mitochondrial_fraction = mitochondrial / totals
    keep = (
        (detected >= MINIMUM_DETECTED_GENES)
        & (mitochondrial_fraction <= MAXIMUM_MITOCHONDRIAL_FRACTION)
        & (absolute_mad <= MAXIMUM_ABSOLUTE_LOG_LIBRARY_MAD)
    )
    return keep, {
        "input_cells": int(totals.size),
        "retained_cells": int(np.count_nonzero(keep)),
        "minimum_detected_genes": MINIMUM_DETECTED_GENES,
        "maximum_mitochondrial_fraction": MAXIMUM_MITOCHONDRIAL_FRACTION,
        "maximum_absolute_log_library_mad": MAXIMUM_ABSOLUTE_LOG_LIBRARY_MAD,
        "log_library_median": center,
        "log_library_raw_mad": raw_mad,
        "log_library_scaled_mad": scaled_mad,
        "mad_scale_constant": 1.4826,
        "failed_detected_genes": int(
            np.count_nonzero(detected < MINIMUM_DETECTED_GENES)
        ),
        "failed_mitochondrial_fraction": int(
            np.count_nonzero(mitochondrial_fraction > MAXIMUM_MITOCHONDRIAL_FRACTION)
        ),
        "failed_log_library_mad": int(
            np.count_nonzero(absolute_mad > MAXIMUM_ABSOLUTE_LOG_LIBRARY_MAD)
        ),
    }


def extract_response_counts(
    stream: BinaryIO,
    selected_columns: np.ndarray,
    response_records: Sequence[metadata.ResponseGeneCrosswalk],
    *,
    expected_rows: int = EXPECTED_TRANSCRIPTOME_ROWS,
    expected_columns: int = EXPECTED_COLUMNS,
) -> np.ndarray:
    """Extract and symbol-collapse exactly the frozen 374 response genes."""

    selected = _validated_selected_columns(selected_columns, expected_columns)
    records = tuple(response_records)
    if len(records) != metadata.EXPECTED_RESPONSE_GENES:
        raise AssignmentStageError("response axis is not exactly 374 records")
    row_to_response: dict[int, int] = {}
    for record in records:
        for source_row in record.source_feature_indices:
            if source_row in row_to_response:
                raise AssignmentStageError(
                    "one source feature maps to multiple response records"
                )
            row_to_response[source_row] = record.response_index

    output = np.zeros((selected.size, len(records)), dtype=np.int32)
    selected_position = 0
    for row, column, value in iter_column_major_coordinates(
        stream,
        expected_rows=expected_rows,
        expected_columns=expected_columns,
    ):
        response_index = row_to_response.get(row)
        if response_index is None:
            continue
        while (
            selected_position < selected.size and selected[selected_position] < column
        ):
            selected_position += 1
        if selected_position < selected.size and selected[selected_position] == column:
            current = int(output[selected_position, response_index]) + value
            if current > np.iinfo(np.int32).max:
                raise AssignmentStageError("response count exceeds int32 range")
            output[selected_position, response_index] = current
    return output


def _csr_from_calls(
    calls: Sequence[guide_calling.CellGuideCall],
    keep: np.ndarray,
    attribute: str,
) -> tuple[np.ndarray, np.ndarray]:
    pointers = [0]
    values: list[np.ndarray] = []
    for call, retained in zip(calls, keep, strict=True):
        if not retained:
            continue
        current = np.asarray(getattr(call, attribute), dtype=np.int32)
        values.append(current)
        pointers.append(pointers[-1] + current.size)
    concatenated = (
        np.concatenate(values).astype(np.int32, copy=False)
        if values and pointers[-1]
        else np.empty(0, dtype=np.int32)
    )
    return np.asarray(pointers, dtype=np.int64), concatenated


def validate_assignment_arrays(
    arrays: dict[str, np.ndarray],
    *,
    guide_to_target: np.ndarray,
    non_targeting_guides: np.ndarray,
    hash_condition_codes: np.ndarray,
) -> None:
    """Validate every per-cell and CSR invariant before publishing a channel."""

    cells = np.asarray(arrays["cell_indices"], dtype=np.int64)
    if cells.ndim != 1 or (cells.size and np.any(np.diff(cells) <= 0)):
        raise AssignmentStageError("final cell indices must be strictly increasing")
    row_count = cells.size
    per_cell = (
        "barcodes",
        "hash_label_indices",
        "condition_codes",
        "library_size",
        "detected_genes",
        "mitochondrial_counts",
        "total_guide_umi",
        "called_guide_burden",
        "non_targeting_guide_burden",
    )
    for name in per_cell:
        if np.asarray(arrays[name]).shape != (row_count,):
            raise AssignmentStageError(f"{name} does not match the final cell axis")
    barcodes = np.asarray(arrays["barcodes"])
    if barcodes.dtype.kind != "U" or len(set(barcodes.tolist())) != row_count:
        raise AssignmentStageError("final barcodes must be unique strings")

    hash_labels = np.asarray(arrays["hash_label_indices"], dtype=np.int64)
    conditions = np.asarray(arrays["condition_codes"], dtype=np.int64)
    hash_conditions = np.asarray(hash_condition_codes, dtype=np.int64)
    if np.any(hash_labels < 0) or np.any(hash_labels >= hash_conditions.size):
        raise AssignmentStageError("hash-label index is outside the canonical axis")
    if not np.array_equal(conditions, hash_conditions[hash_labels]):
        raise AssignmentStageError("hash labels and condition codes disagree")

    library_size = np.asarray(arrays["library_size"], dtype=np.int64)
    detected = np.asarray(arrays["detected_genes"], dtype=np.int64)
    mitochondrial = np.asarray(arrays["mitochondrial_counts"], dtype=np.int64)
    if (
        np.any(library_size <= 0)
        or np.any(detected < MINIMUM_DETECTED_GENES)
        or np.any(mitochondrial < 0)
        or np.any(mitochondrial > library_size)
        or np.any(mitochondrial / library_size > MAXIMUM_MITOCHONDRIAL_FRACTION)
    ):
        raise AssignmentStageError("published cells violate expression QC")

    guide_mapping = np.asarray(guide_to_target, dtype=np.int64)
    ntc = np.asarray(non_targeting_guides)
    if guide_mapping.shape != (EXPECTED_GUIDE_ROWS,) or ntc.shape != (
        EXPECTED_GUIDE_ROWS,
    ):
        raise AssignmentStageError("canonical guide annotation axis changed")
    if ntc.dtype != np.bool_:
        raise AssignmentStageError("non-targeting guide annotation is not boolean")

    guide_indptr = np.asarray(arrays["called_guide_indptr"], dtype=np.int64)
    guide_indices = np.asarray(arrays["called_guide_indices"], dtype=np.int64)
    target_indptr = np.asarray(arrays["target_indptr"], dtype=np.int64)
    target_indices = np.asarray(arrays["target_indices"], dtype=np.int64)
    for name, pointers, values in (
        ("called guides", guide_indptr, guide_indices),
        ("targets", target_indptr, target_indices),
    ):
        if (
            pointers.shape != (row_count + 1,)
            or pointers[0] != 0
            or np.any(np.diff(pointers) < 0)
            or pointers[-1] != values.size
        ):
            raise AssignmentStageError(f"{name} CSR pointers are invalid")
    if np.any(guide_indices < 0) or np.any(guide_indices >= EXPECTED_GUIDE_ROWS):
        raise AssignmentStageError("called guide is outside the canonical axis")
    target_count = int(guide_mapping.max(initial=-1)) + 1
    if np.any(target_indices < 0) or np.any(target_indices >= target_count):
        raise AssignmentStageError("called target is outside the canonical axis")

    called_burden = np.asarray(arrays["called_guide_burden"], dtype=np.int64)
    ntc_burden = np.asarray(
        arrays["non_targeting_guide_burden"],
        dtype=np.int64,
    )
    if not np.array_equal(called_burden, np.diff(guide_indptr)):
        raise AssignmentStageError("called-guide burden disagrees with guide CSR")
    for row in range(row_count):
        guide_start, guide_stop = guide_indptr[row : row + 2]
        target_start, target_stop = target_indptr[row : row + 2]
        row_guides = guide_indices[guide_start:guide_stop]
        row_targets = target_indices[target_start:target_stop]
        if row_guides.size and np.any(np.diff(row_guides) <= 0):
            raise AssignmentStageError("called guides are not canonical and unique")
        if row_targets.size and np.any(np.diff(row_targets) <= 0):
            raise AssignmentStageError("called targets are not canonical and unique")
        if ntc_burden[row] != np.count_nonzero(ntc[row_guides]):
            raise AssignmentStageError("non-targeting burden disagrees with guide CSR")
        mapped = np.unique(guide_mapping[row_guides][guide_mapping[row_guides] >= 0])
        if not np.array_equal(row_targets, mapped):
            raise AssignmentStageError("target CSR disagrees with called guides")


def _histogram(values: np.ndarray) -> dict[str, int]:
    return {
        str(int(value)): int(count)
        for value, count in zip(*np.unique(values, return_counts=True), strict=True)
    }


def derive_target_support(
    guide_cells_after_qc: np.ndarray,
    target_cells_after_qc: np.ndarray,
    target_untreated_channels: np.ndarray,
    guide_to_target: np.ndarray,
) -> dict[str, np.ndarray]:
    """Derive the frozen all-channel target support diagnostics."""

    guide_counts = np.asarray(guide_cells_after_qc, dtype=np.int64)
    target_counts = np.asarray(target_cells_after_qc, dtype=np.int64)
    channel_counts = np.asarray(target_untreated_channels, dtype=np.int64)
    mapping = np.asarray(guide_to_target, dtype=np.int64)
    if guide_counts.ndim != 2 or guide_counts.shape[0] != 2:
        raise AssignmentStageError("guide support counts must have two conditions")
    if target_counts.ndim != 2 or target_counts.shape[0] != 2:
        raise AssignmentStageError("target support counts must have two conditions")
    if mapping.shape != (guide_counts.shape[1],):
        raise AssignmentStageError("guide-to-target axis disagrees with support counts")
    if channel_counts.shape != (target_counts.shape[1],):
        raise AssignmentStageError("target channel-support axis changed")
    if np.any(guide_counts < 0) or np.any(target_counts < 0):
        raise AssignmentStageError("support counts cannot be negative")
    if np.any(mapping < metadata.NO_UNIQUE_TARGET) or np.any(
        mapping >= target_counts.shape[1]
    ):
        raise AssignmentStageError("guide-to-target mapping is outside target axis")

    supported_guides = mapping >= 0
    target_count = target_counts.shape[1]
    distinct_called_guides = np.bincount(
        mapping[supported_guides],
        weights=(guide_counts[0, supported_guides] > 0).astype(np.int64),
        minlength=target_count,
    ).astype(np.int32)
    untreated_guides_with_three_cells = np.bincount(
        mapping[supported_guides],
        weights=(guide_counts[0, supported_guides] >= 3).astype(np.int64),
        minlength=target_count,
    ).astype(np.int32)
    eligible = (
        (distinct_called_guides >= 3)
        & (untreated_guides_with_three_cells >= 2)
        & (target_counts[0] >= 20)
        & (channel_counts >= 4)
    )
    return {
        "distinct_called_guides": distinct_called_guides,
        "untreated_guides_with_at_least_three_cells": (
            untreated_guides_with_three_cells
        ),
        "untreated_support_channels": channel_counts.astype(np.int16),
        "eligible": eligible,
    }


def write_metadata_crosswalk_report(
    crosswalks: metadata.MetadataCrosswalks,
) -> dict[str, Any]:
    """Publish a machine-checkable exact metadata crosswalk."""

    response_indptr = [0]
    response_source_indices: list[int] = []
    for record in crosswalks.response.records:
        response_source_indices.extend(record.source_feature_indices)
        response_indptr.append(len(response_source_indices))
    output_path = WORK / "metadata-crosswalk.npz"
    atomic_savez(
        output_path,
        guide_ids=np.asarray(crosswalks.guide.guide_ids, dtype="U"),
        guide_labels=np.asarray(crosswalks.guide.guide_labels, dtype="U"),
        guide_status=np.asarray(crosswalks.guide.guide_status, dtype="U"),
        target_symbols=np.asarray(crosswalks.guide.target_symbols, dtype="U"),
        guide_to_target=crosswalks.guide.guide_to_target,
        non_targeting_guides=crosswalks.guide.non_targeting_guides,
        hash_labels=np.asarray(
            [record.label for record in crosswalks.hashes.records],
            dtype="U",
        ),
        hash_condition_codes=np.asarray(
            [
                0 if record.condition == "untreated" else 1
                for record in crosswalks.hashes.records
            ],
            dtype=np.uint8,
        ),
        response_symbols=np.asarray(
            [record.symbol for record in crosswalks.response.records],
            dtype="U",
        ),
        response_source_indptr=np.asarray(response_indptr, dtype=np.int64),
        response_source_indices=np.asarray(
            response_source_indices,
            dtype=np.int32,
        ),
        response_vcc_indices=np.asarray(
            [
                -1 if record.vcc_gene_index is None else record.vcc_gene_index
                for record in crosswalks.response.records
            ],
            dtype=np.int32,
        ),
        response_is_panel_target=np.asarray(
            [record.is_panel_target for record in crosswalks.response.records],
            dtype=bool,
        ),
    )
    report_path = WORK / "metadata-crosswalk-report.json"
    try:
        displayed_output_path = str(output_path.relative_to(HERE))
    except ValueError:
        displayed_output_path = str(output_path)
    report = {
        "status": "complete",
        "guides": len(crosswalks.guide.guide_ids),
        "guide_status_counts": crosswalks.guide.status_counts,
        "source_targets": len(crosswalks.guide.target_symbols),
        "non_targeting_guides": int(
            np.count_nonzero(crosswalks.guide.non_targeting_guides)
        ),
        "panel_covered_targets": len(crosswalks.guide.panel_covered_targets),
        "panel_missing_targets": list(crosswalks.guide.panel_missing_targets),
        "hash_labels": [
            {
                "feature_index": record.feature_index,
                "label": record.label,
                "condition": record.condition,
            }
            for record in crosswalks.hashes.records
        ],
        "response_genes": len(crosswalks.response.records),
        "response_genes_on_vcc_axis": int(
            np.count_nonzero(
                [
                    record.vcc_gene_index is not None
                    for record in crosswalks.response.records
                ]
            )
        ),
        "pds_visible_response_genes": len(crosswalks.response.visible_vcc_indices),
        "missing_source_symbols": list(crosswalks.response.missing_source_symbols),
        "ambiguous_source_symbols": list(crosswalks.response.ambiguous_source_symbols),
        "missing_vcc_symbols": list(crosswalks.response.missing_vcc_symbols),
        "metadata_sha256": sha256_file(HERE / "metadata.py"),
        "output_path": displayed_output_path,
        "output_sha256": sha256_file(output_path),
    }
    atomic_write_json(report_path, report)
    return report


def _channel_record(schema: dict[str, Any], channel: int) -> dict[str, Any]:
    matches = [
        record
        for record in schema.get("channels", [])
        if int(record.get("channel", -1)) == channel
    ]
    if len(matches) != 1:
        raise AssignmentStageError(
            f"schema does not contain exactly one channel {channel}"
        )
    return matches[0]


def _load_production_cell_calls(channel: int) -> tuple[np.ndarray, dict[str, Any]]:
    directory = WORK / f"channel{channel:02d}"
    report_path = directory / "cell-calls-production.json"
    output_path = directory / "cell-calls-production.npz"
    report = load_json(report_path)
    if report.get("status") != "complete" or report.get("iterations") != 10_000:
        raise AssignmentStageError(f"channel {channel} production calls are incomplete")
    if report.get("output_sha256") != sha256_file(output_path):
        raise AssignmentStageError(f"channel {channel} cell-call output changed")
    with np.load(output_path, allow_pickle=False) as values:
        accepted = values["accepted_indices"].astype(np.int64, copy=False)
    _validated_selected_columns(accepted, EXPECTED_COLUMNS)
    if report.get("accepted_barcodes") != int(accepted.size):
        raise AssignmentStageError(
            f"channel {channel} accepted-cell count disagrees with report"
        )
    return accepted, report


def _load_prepared_qc(channel: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = WORK / f"channel{channel:02d}" / "cell-call-inputs.npz"
    with np.load(path, allow_pickle=False) as values:
        return (
            values["totals"].astype(np.int64, copy=False),
            values["detected_features"].astype(np.int64, copy=False),
            values["mitochondrial_counts"].astype(np.int64, copy=False),
        )


def _load_schema() -> dict[str, Any]:
    schema = load_json(SCHEMA_REPORT)
    if schema.get("status") != "pass":
        raise AssignmentStageError("schema audit did not pass")
    return schema


def _hash_round_report(result: demultiplex.DemultiplexResult) -> list[dict[str, Any]]:
    return [
        {
            "round": round_result.round_number,
            "input_cells": len(round_result.active_cell_ids),
            "selected_q": round_result.selected_q,
            "selected_q_index": round_result.selected_q_index,
            "selected_call_counts": dict(
                sorted(Counter(round_result.selected_calls).items())
            ),
        }
        for round_result in result.rounds
    ]


def _verified_source(
    verification: dict[str, Any] | None = None,
) -> dict[str, Any]:
    verified = (
        source_artifacts.verify_source_artifacts()
        if verification is None
        else verification
    )
    if verified.get("status") != "pass":
        raise AssignmentStageError("source-artifact verification did not pass")
    if len(verified.get("channels", ())) != CHANNEL_COUNT:
        raise AssignmentStageError(
            "source-artifact verification did not cover all 16 channels"
        )
    return verified


def _validate_assignment_lineage(
    channel: int,
    report: dict[str, Any],
    verified: dict[str, Any],
) -> None:
    """Bind one assignment artifact to current source outputs and code."""

    directory = WORK / f"channel{channel:02d}"
    cell_call_report_path = directory / "cell-calls-production.json"
    cell_call_output_path = directory / "cell-calls-production.npz"
    cell_call_report = load_json(cell_call_report_path)
    expected_code_hashes = {
        "assignment_stage_sha256": sha256_file(Path(__file__)),
        "metadata_sha256": sha256_file(HERE / "metadata.py"),
        "demultiplex_sha256": sha256_file(HERE / "demultiplex.py"),
        "guide_calling_sha256": sha256_file(HERE / "guide_calling.py"),
        "source_artifacts_sha256": sha256_file(HERE / "source_artifacts.py"),
    }
    for field, expected in expected_code_hashes.items():
        if report.get(field) != expected:
            raise AssignmentStageError(
                f"channel {channel} assignment lineage changed at {field}"
            )
    if report.get("cell_call_report_sha256") != sha256_file(cell_call_report_path):
        raise AssignmentStageError(
            f"channel {channel} assignment used a stale cell-call report"
        )
    if (
        report.get("cell_call_output_sha256")
        != cell_call_report.get("output_sha256")
        or report.get("cell_call_output_sha256") != sha256_file(cell_call_output_path)
    ):
        raise AssignmentStageError(
            f"channel {channel} assignment used stale cell-call output"
        )
    if report.get("source_raw_archive_sha256") != verified.get(
        "raw_archive_sha256"
    ):
        raise AssignmentStageError(
            f"channel {channel} assignment used a different source archive"
        )
    if report.get("source_oracle_receipt_sha256") != verified.get(
        "oracle_receipt_sha256"
    ):
        raise AssignmentStageError(
            f"channel {channel} assignment used a different R-oracle receipt"
        )


def build_channel_assignment(
    channel: int,
    *,
    source_verification: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run hash assignment, guide calling, and expression QC for one channel."""

    if not 1 <= channel <= CHANNEL_COUNT:
        raise AssignmentStageError("channel must lie in 1..16")
    verified = _verified_source(source_verification)
    schema = _load_schema()
    record = _channel_record(schema, channel)
    accession = str(record["accession"])
    accepted, cell_call_report = _load_production_cell_calls(channel)
    crosswalks = metadata.load_pinned_crosswalks(channel=channel)
    label_names = tuple(record.label for record in crosswalks.hashes.records)

    with tarfile.open(ARCHIVE, mode="r:") as archive:
        with gzip_member(
            archive,
            member_name(accession, channel, "labels", "matrix"),
        ) as stream:
            label_counts = extract_selected_dense(
                stream,
                accepted,
                expected_rows=EXPECTED_LABEL_ROWS,
                expected_columns=EXPECTED_COLUMNS,
                dtype=np.int32,
            )
    demultiplexed = demultiplex.demultiplex(
        label_counts,
        tuple(str(index) for index in accepted),
        label_names,
    )
    call_array = np.asarray(demultiplexed.calls)
    label_to_index = {label: index for index, label in enumerate(label_names)}
    singlet_positions = np.flatnonzero(
        np.fromiter(
            (call in label_to_index for call in demultiplexed.calls),
            dtype=bool,
            count=len(demultiplexed.calls),
        )
    )
    singlet_columns = accepted[singlet_positions]
    if singlet_columns.size == 0:
        raise AssignmentStageError(f"channel {channel} has no retained hash singlets")
    singlet_hash_labels = call_array[singlet_positions]
    singlet_condition_codes = np.asarray(
        [
            0 if crosswalks.hashes.label_to_condition[label] == "untreated" else 1
            for label in singlet_hash_labels
        ],
        dtype=np.uint8,
    )

    guide_matrix_name = member_name(accession, channel, "guides", "matrix")

    def guide_batches() -> Iterator[guide_calling.SparseGuideBatch]:
        with tarfile.open(ARCHIVE, mode="r:") as archive:
            with gzip_member(archive, guide_matrix_name) as stream:
                yield from iter_selected_sparse_batches(
                    stream,
                    singlet_columns,
                    expected_rows=EXPECTED_GUIDE_ROWS,
                    expected_columns=EXPECTED_COLUMNS,
                )

    channel_null = guide_calling.fit_channel_guide_null(guide_batches())
    calls = tuple(
        guide_calling.iter_channel_guide_calls(
            guide_batches(),
            channel_null=channel_null,
            guide_to_target=crosswalks.guide.guide_to_target,
            non_targeting_guides=crosswalks.guide.non_targeting_guides,
        )
    )
    if tuple(call.cell_index for call in calls) != tuple(singlet_columns):
        raise AssignmentStageError("guide-call cell axis differs from hash singlets")

    all_totals, all_detected, all_mitochondrial = _load_prepared_qc(channel)
    qc_keep, qc_report = expression_qc_mask(
        all_totals[singlet_columns],
        all_detected[singlet_columns],
        all_mitochondrial[singlet_columns],
    )
    final_columns = singlet_columns[qc_keep]
    final_positions = singlet_positions[qc_keep]
    final_hash_labels = call_array[final_positions]
    hash_label_indices = np.asarray(
        [label_to_index[label] for label in final_hash_labels],
        dtype=np.int8,
    )
    condition_codes = singlet_condition_codes[qc_keep]

    target_count = len(crosswalks.guide.target_symbols)
    guide_cells_before_qc = np.zeros(
        (2, EXPECTED_GUIDE_ROWS),
        dtype=np.int32,
    )
    guide_cells_after_qc = np.zeros_like(guide_cells_before_qc)
    target_cells_before_qc = np.zeros((2, target_count), dtype=np.int32)
    target_cells_after_qc = np.zeros_like(target_cells_before_qc)
    for call, condition, retained in zip(
        calls,
        singlet_condition_codes,
        qc_keep,
        strict=True,
    ):
        guide_cells_before_qc[condition, call.called_guide_indices] += 1
        target_cells_before_qc[condition, call.target_indices] += 1
        if retained:
            guide_cells_after_qc[condition, call.called_guide_indices] += 1
            target_cells_after_qc[condition, call.target_indices] += 1

    called_guide_indptr, called_guide_indices = _csr_from_calls(
        calls,
        qc_keep,
        "called_guide_indices",
    )
    target_indptr, target_indices = _csr_from_calls(
        calls,
        qc_keep,
        "target_indices",
    )
    total_guide_umi = np.asarray(
        [
            call.total_guide_umi
            for call, keep in zip(calls, qc_keep, strict=True)
            if keep
        ],
        dtype=np.int64,
    )
    called_guide_burden = np.asarray(
        [
            call.total_called_guide_burden
            for call, keep in zip(calls, qc_keep, strict=True)
            if keep
        ],
        dtype=np.uint16,
    )
    ntc_burden = np.asarray(
        [
            call.non_targeting_guide_burden
            for call, keep in zip(calls, qc_keep, strict=True)
            if keep
        ],
        dtype=np.uint16,
    )

    with tarfile.open(ARCHIVE, mode="r:") as archive:
        with gzip_member(
            archive,
            member_name(accession, channel, "transcriptome", "barcode"),
        ) as stream:
            barcodes = extract_selected_barcodes(
                stream,
                final_columns,
                expected_columns=EXPECTED_COLUMNS,
            )

    output_path = WORK / f"channel{channel:02d}" / "assignments.npz"
    support_path = WORK / f"channel{channel:02d}" / "support-counts.npz"
    report_path = WORK / f"channel{channel:02d}" / "assignment-report.json"
    assignment_arrays = {
        "cell_indices": final_columns.astype(np.int32),
        "barcodes": barcodes,
        "hash_label_indices": hash_label_indices,
        "condition_codes": condition_codes,
        "library_size": all_totals[final_columns],
        "detected_genes": all_detected[final_columns].astype(np.uint16),
        "mitochondrial_counts": all_mitochondrial[final_columns],
        "total_guide_umi": total_guide_umi,
        "called_guide_burden": called_guide_burden,
        "non_targeting_guide_burden": ntc_burden,
        "called_guide_indptr": called_guide_indptr,
        "called_guide_indices": called_guide_indices,
        "target_indptr": target_indptr,
        "target_indices": target_indices,
    }
    validate_assignment_arrays(
        assignment_arrays,
        guide_to_target=crosswalks.guide.guide_to_target,
        non_targeting_guides=crosswalks.guide.non_targeting_guides,
        hash_condition_codes=np.asarray(
            [
                0 if record.condition == "untreated" else 1
                for record in crosswalks.hashes.records
            ],
            dtype=np.uint8,
        ),
    )
    atomic_savez(output_path, **assignment_arrays)
    atomic_savez(
        support_path,
        guide_cells_before_qc=guide_cells_before_qc,
        guide_cells_after_qc=guide_cells_after_qc,
        target_cells_before_qc=target_cells_before_qc,
        target_cells_after_qc=target_cells_after_qc,
    )
    report = {
        "status": "complete",
        "channel": channel,
        "accession": accession,
        "cell_call_report_sha256": sha256_file(
            WORK / f"channel{channel:02d}" / "cell-calls-production.json"
        ),
        "cell_call_output_sha256": cell_call_report["output_sha256"],
        "assignment_stage_sha256": sha256_file(Path(__file__)),
        "metadata_sha256": sha256_file(HERE / "metadata.py"),
        "demultiplex_sha256": sha256_file(HERE / "demultiplex.py"),
        "guide_calling_sha256": sha256_file(HERE / "guide_calling.py"),
        "source_artifacts_sha256": sha256_file(HERE / "source_artifacts.py"),
        "source_raw_archive_sha256": verified["raw_archive_sha256"],
        "source_oracle_receipt_sha256": verified["oracle_receipt_sha256"],
        "accepted_expression_cells": int(accepted.size),
        "hash_call_counts": dict(sorted(Counter(demultiplexed.calls).items())),
        "hash_rounds": _hash_round_report(demultiplexed),
        "hash_singlet_cells": int(singlet_columns.size),
        "expression_qc": qc_report,
        "retained_cells": int(final_columns.size),
        "untreated_cells": int(np.count_nonzero(condition_codes == 0)),
        "activated_cells": int(np.count_nonzero(condition_codes == 1)),
        "guide_null_total_umi": channel_null.total_guide_umi,
        "guide_null_nonzero_features": int(np.count_nonzero(channel_null.guide_umi)),
        "called_guide_burden_histogram": _histogram(called_guide_burden),
        "median_called_guides": float(np.median(called_guide_burden)),
        "output_path": str(output_path.relative_to(HERE)),
        "output_sha256": sha256_file(output_path),
        "support_path": str(support_path.relative_to(HERE)),
        "support_sha256": sha256_file(support_path),
        "support_axes": {
            "condition_codes": {"0": "untreated", "1": "activated"},
            "guide_axis": "metadata-crosswalk canonical guide order",
            "target_axis": "metadata-crosswalk canonical source-target order",
            "before_qc": "all retained deMULTIplex hash singlets",
            "after_qc": "hash singlets passing frozen expression QC",
        },
    }
    atomic_write_json(report_path, report)
    return report


def build_channel_response(
    channel: int,
    *,
    source_verification: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Extract only the frozen response panel for one assigned channel."""

    if not 1 <= channel <= CHANNEL_COUNT:
        raise AssignmentStageError("channel must lie in 1..16")
    verified = _verified_source(source_verification)
    schema = _load_schema()
    record = _channel_record(schema, channel)
    accession = str(record["accession"])
    assignment_path = WORK / f"channel{channel:02d}" / "assignments.npz"
    assignment_report_path = WORK / f"channel{channel:02d}" / "assignment-report.json"
    assignment_report = load_json(assignment_report_path)
    if assignment_report.get("output_sha256") != sha256_file(assignment_path):
        raise AssignmentStageError(f"channel {channel} assignments changed")
    _validate_assignment_lineage(channel, assignment_report, verified)
    with np.load(assignment_path, allow_pickle=False) as assignment:
        selected = assignment["cell_indices"].astype(np.int64, copy=False)
        library_size = assignment["library_size"].astype(np.int64, copy=False)
    crosswalks = metadata.load_pinned_crosswalks(channel=channel)
    with tarfile.open(ARCHIVE, mode="r:") as archive:
        with gzip_member(
            archive,
            member_name(accession, channel, "transcriptome", "matrix"),
        ) as stream:
            response_counts = extract_response_counts(
                stream,
                selected,
                crosswalks.response.records,
            )
    if np.any(response_counts.sum(axis=1, dtype=np.int64) > library_size):
        raise AssignmentStageError("response-panel counts exceed full library size")

    output_path = WORK / f"channel{channel:02d}" / "response-counts.npz"
    report_path = WORK / f"channel{channel:02d}" / "response-report.json"
    atomic_savez(
        output_path,
        cell_indices=selected.astype(np.int32),
        response_counts=response_counts,
        library_size=library_size,
    )
    report = {
        "status": "complete",
        "channel": channel,
        "assignment_report_sha256": sha256_file(assignment_report_path),
        "assignment_output_sha256": sha256_file(assignment_path),
        "assignment_stage_sha256": sha256_file(Path(__file__)),
        "metadata_sha256": sha256_file(HERE / "metadata.py"),
        "source_artifacts_sha256": sha256_file(HERE / "source_artifacts.py"),
        "source_raw_archive_sha256": verified["raw_archive_sha256"],
        "retained_cells": int(selected.size),
        "response_genes": int(response_counts.shape[1]),
        "visible_vcc_response_genes": len(crosswalks.response.visible_vcc_indices),
        "response_panel_umi": int(response_counts.sum(dtype=np.int64)),
        "output_path": str(output_path.relative_to(HERE)),
        "output_sha256": sha256_file(output_path),
    }
    atomic_write_json(report_path, report)
    return report


def summarize_assignments() -> dict[str, Any]:
    """Aggregate all channels and enforce the frozen assay-compatibility gate."""

    verified = _verified_source()
    reports: list[dict[str, Any]] = []
    burdens: list[np.ndarray] = []
    guide_cells_before_qc = np.zeros((2, EXPECTED_GUIDE_ROWS), dtype=np.int64)
    guide_cells_after_qc = np.zeros_like(guide_cells_before_qc)
    target_cells_before_qc: np.ndarray | None = None
    target_cells_after_qc: np.ndarray | None = None
    target_untreated_channels: np.ndarray | None = None
    total_cells = 0
    untreated_cells = 0
    for channel in range(1, CHANNEL_COUNT + 1):
        directory = WORK / f"channel{channel:02d}"
        report_path = directory / "assignment-report.json"
        output_path = directory / "assignments.npz"
        support_path = directory / "support-counts.npz"
        report = load_json(report_path)
        if report.get("status") != "complete":
            raise AssignmentStageError(
                f"channel {channel} assignment report is incomplete"
            )
        if report.get("output_sha256") != sha256_file(output_path):
            raise AssignmentStageError(f"channel {channel} assignments changed")
        _validate_assignment_lineage(channel, report, verified)
        if report.get("support_sha256") != sha256_file(support_path):
            raise AssignmentStageError(f"channel {channel} support counts changed")
        with np.load(output_path, allow_pickle=False) as values:
            current_burden = values["called_guide_burden"].astype(
                np.int64,
                copy=False,
            )
            conditions = values["condition_codes"].astype(np.int64, copy=False)
        if current_burden.size != conditions.size:
            raise AssignmentStageError(f"channel {channel} assignment vectors disagree")
        with np.load(support_path, allow_pickle=False) as values:
            channel_guide_before = values["guide_cells_before_qc"].astype(
                np.int64,
                copy=False,
            )
            channel_guide_after = values["guide_cells_after_qc"].astype(
                np.int64,
                copy=False,
            )
            channel_target_before = values["target_cells_before_qc"].astype(
                np.int64,
                copy=False,
            )
            channel_target_after = values["target_cells_after_qc"].astype(
                np.int64,
                copy=False,
            )
        if (
            channel_guide_before.shape != guide_cells_before_qc.shape
            or channel_guide_after.shape != guide_cells_after_qc.shape
            or channel_target_before.ndim != 2
            or channel_target_before.shape[0] != 2
            or channel_target_after.shape != channel_target_before.shape
        ):
            raise AssignmentStageError(f"channel {channel} support-count axes changed")
        if target_cells_before_qc is None:
            target_cells_before_qc = np.zeros_like(
                channel_target_before,
                dtype=np.int64,
            )
            target_cells_after_qc = np.zeros_like(
                channel_target_after,
                dtype=np.int64,
            )
            target_untreated_channels = np.zeros(
                channel_target_after.shape[1],
                dtype=np.int16,
            )
        if channel_target_before.shape != target_cells_before_qc.shape:
            raise AssignmentStageError(f"channel {channel} target-support axis changed")
        guide_cells_before_qc += channel_guide_before
        guide_cells_after_qc += channel_guide_after
        target_cells_before_qc += channel_target_before
        target_cells_after_qc += channel_target_after
        target_untreated_channels += channel_target_after[0] > 0
        reports.append(report)
        burdens.append(current_burden)
        total_cells += int(current_burden.size)
        untreated_cells += int(np.count_nonzero(conditions == 0))

    if (
        target_cells_before_qc is None
        or target_cells_after_qc is None
        or target_untreated_channels is None
    ):
        raise AssignmentStageError("no channel support counts were loaded")
    crosswalks = metadata.load_pinned_crosswalks()
    metadata_report = write_metadata_crosswalk_report(crosswalks)
    support = derive_target_support(
        guide_cells_after_qc,
        target_cells_after_qc,
        target_untreated_channels,
        crosswalks.guide.guide_to_target,
    )
    target_support_path = WORK / "target-support.npz"
    atomic_savez(
        target_support_path,
        guide_cells_before_qc=guide_cells_before_qc,
        guide_cells_after_qc=guide_cells_after_qc,
        target_cells_before_qc=target_cells_before_qc,
        target_cells_after_qc=target_cells_after_qc,
        **support,
    )

    all_burdens = np.concatenate(burdens)
    untreated_fraction = untreated_cells / total_cells
    median_called_guides = float(np.median(all_burdens))
    gates = {
        "retained_cells": RETAINED_CELL_RANGE[0]
        <= total_cells
        <= RETAINED_CELL_RANGE[1],
        "untreated_fraction": UNTREATED_FRACTION_RANGE[0]
        <= untreated_fraction
        <= UNTREATED_FRACTION_RANGE[1],
        "median_called_guides": MEDIAN_CALLED_GUIDE_RANGE[0]
        <= median_called_guides
        <= MEDIAN_CALLED_GUIDE_RANGE[1],
    }
    result = {
        "status": "pass" if all(gates.values()) else "fail",
        "retained_cells": total_cells,
        "untreated_cells": untreated_cells,
        "untreated_fraction": untreated_fraction,
        "median_called_guides": median_called_guides,
        "called_guide_burden_histogram": _histogram(all_burdens),
        "eligible_source_targets": int(np.count_nonzero(support["eligible"])),
        "target_support_thresholds": {
            "minimum_distinct_called_guides": 3,
            "minimum_distinct_guides_with_at_least_3_untreated_cells": 2,
            "minimum_untreated_carrier_cells": 20,
            "minimum_untreated_channels": 4,
        },
        "target_support_path": str(target_support_path.relative_to(HERE)),
        "target_support_sha256": sha256_file(target_support_path),
        "metadata_crosswalk_report_sha256": sha256_file(
            WORK / "metadata-crosswalk-report.json"
        ),
        "metadata_crosswalk_output_sha256": metadata_report["output_sha256"],
        "assignment_stage_sha256": sha256_file(Path(__file__)),
        "source_artifacts_sha256": sha256_file(HERE / "source_artifacts.py"),
        "source_raw_archive_sha256": verified["raw_archive_sha256"],
        "source_oracle_receipt_sha256": verified["oracle_receipt_sha256"],
        "gates": gates,
        "frozen_ranges": {
            "retained_cells": list(RETAINED_CELL_RANGE),
            "untreated_fraction": list(UNTREATED_FRACTION_RANGE),
            "median_called_guides": list(MEDIAN_CALLED_GUIDE_RANGE),
        },
        "channel_report_sha256": {
            str(channel): sha256_file(
                WORK / f"channel{channel:02d}" / "assignment-report.json"
            )
            for channel in range(1, CHANNEL_COUNT + 1)
        },
    }
    atomic_write_json(WORK / "assignment-summary.json", result)
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument(
        "command",
        choices=("assign", "assign-all", "response", "response-all", "summarize"),
    )
    result.add_argument("--channel", type=int, choices=range(1, 17))
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command in {"assign", "response"} and args.channel is None:
        print("ERROR: --channel is required", file=sys.stderr)
        return 2
    try:
        with exclusive_lock():
            if args.command == "assign":
                report = build_channel_assignment(args.channel)
            elif args.command == "assign-all":
                verified = _verified_source()
                channels = []
                for channel in range(1, CHANNEL_COUNT + 1):
                    channel_report = build_channel_assignment(
                        channel,
                        source_verification=verified,
                    )
                    channels.append(channel_report)
                    print(
                        json.dumps(
                            {
                                "progress": "assignment-complete",
                                "channel": channel,
                                "retained_cells": channel_report["retained_cells"],
                                "median_called_guides": channel_report[
                                    "median_called_guides"
                                ],
                            },
                            sort_keys=True,
                        ),
                        flush=True,
                    )
                report = {
                    "status": "complete",
                    "stage": "assign-all",
                    "channels": channels,
                }
            elif args.command == "response":
                report = build_channel_response(args.channel)
            elif args.command == "response-all":
                verified = _verified_source()
                channels = []
                for channel in range(1, CHANNEL_COUNT + 1):
                    channel_report = build_channel_response(
                        channel,
                        source_verification=verified,
                    )
                    channels.append(channel_report)
                    print(
                        json.dumps(
                            {
                                "progress": "response-complete",
                                "channel": channel,
                                "retained_cells": channel_report["retained_cells"],
                                "response_panel_umi": channel_report[
                                    "response_panel_umi"
                                ],
                            },
                            sort_keys=True,
                        ),
                        flush=True,
                    )
                report = {
                    "status": "complete",
                    "stage": "response-all",
                    "channels": channels,
                }
            else:
                report = summarize_assignments()
    except (
        AssignmentStageError,
        demultiplex.DemultiplexError,
        guide_calling.GuideCallingError,
        metadata.MetadataError,
        preprocess.PreprocessError,
        source_artifacts.SourceArtifactError,
        OSError,
        tarfile.TarError,
    ) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report.get("status") != "fail" else 2


if __name__ == "__main__":
    raise SystemExit(main())
