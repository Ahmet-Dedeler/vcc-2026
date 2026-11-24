#!/usr/bin/env python3
"""Pure contracts for experiment 112's Orion control-label correction.

This module performs no file discovery or I/O. In particular, it cannot read
raw Parquet expression, held truth, donor artifacts, or whitening operators.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


EXPERIMENT = "112-orion-control-label-correction"
SCHEMA_VERSION = 1
TARGET_SUM = 50_000.0
EXACT_CONTROL_LABEL = "Non-Targeting"
MISCLASSIFIED_TARGET = "NEGR1"
METADATA_ONLY_COLUMNS = ("gene_target", "pass_guide_filter")
LEGACY_CONTROL_PATTERNS = (
    r"^non.?targeting",
    r"^ntc",
    r"^control",
    r"^safe.?harbou?r",
    r"^neg",
)
EXPECTED_LINES = {
    "HCT116": {
        "exact_control_count": 165_777,
        "negr1_count": 154,
        "legacy_control_count": 165_931,
        "row_groups": 109,
    },
    "HEK293T": {
        "exact_control_count": 218_838,
        "negr1_count": 233,
        "legacy_control_count": 219_071,
        "row_groups": 223,
    },
}
LEGACY_ARTIFACT_NAMES = (
    "delta.npy",
    "b_ctrl.npy",
    "genes.txt",
    "perts.txt",
    "meta.json",
)
SHA256_LENGTH = 64


class ContractError(ValueError):
    """Raised when a frozen experiment-112 invariant is violated."""


@dataclass(frozen=True)
class ControlCounts:
    """Eligible row counts from metadata-only Orion columns."""

    exact_control: int
    negr1: int
    legacy_control: int


def legacy_is_control(label: object) -> bool:
    """Reproduce the two shared builders' legacy control predicate."""

    normalized = str(label).lower()
    return any(re.match(pattern, normalized) for pattern in LEGACY_CONTROL_PATTERNS)


def corrected_is_control(label: object) -> bool:
    """Recognize only the exact raw Orion control label."""

    return isinstance(label, str) and label == EXACT_CONTROL_LABEL


def expression_filter(label: str) -> tuple[tuple[str, str, object], ...]:
    """Return the only allowed future Parquet expression predicate."""

    if label not in {EXACT_CONTROL_LABEL, MISCLASSIFIED_TARGET}:
        raise ContractError(f"expression projection is not authorized for {label!r}")
    return (
        ("pass_guide_filter", "=", 1),
        ("gene_target", "=", label),
    )


def count_metadata_rows(
    targets: Sequence[object],
    pass_guide_filter: Sequence[object],
) -> ControlCounts:
    """Count the bug from label/filter metadata without touching expression."""

    targets_array = np.asarray(targets, dtype=object)
    eligible_array = np.asarray(pass_guide_filter)
    if targets_array.ndim != 1 or eligible_array.ndim != 1:
        raise ContractError("metadata columns must be one-dimensional")
    if len(targets_array) != len(eligible_array):
        raise ContractError("metadata columns have different row counts")
    eligible = eligible_array == 1
    exact = np.fromiter(
        (corrected_is_control(value) for value in targets_array),
        dtype=bool,
        count=len(targets_array),
    )
    negr1 = targets_array == MISCLASSIFIED_TARGET
    legacy = np.fromiter(
        (legacy_is_control(value) for value in targets_array),
        dtype=bool,
        count=len(targets_array),
    )
    return ControlCounts(
        exact_control=int(np.sum(eligible & exact)),
        negr1=int(np.sum(eligible & negr1)),
        legacy_control=int(np.sum(eligible & legacy)),
    )


def validate_control_counts(line: str, counts: ControlCounts) -> None:
    """Require the frozen metadata-only counts and no extra legacy matches."""

    try:
        expected = EXPECTED_LINES[line]
    except KeyError as exc:
        raise ContractError(f"unknown Orion line {line!r}") from exc
    observed = {
        "exact_control_count": counts.exact_control,
        "negr1_count": counts.negr1,
        "legacy_control_count": counts.legacy_control,
    }
    for name, value in observed.items():
        if value != expected[name]:
            raise ContractError(
                f"{line} {name} changed: observed {value}, expected {expected[name]}"
            )
    if counts.legacy_control != counts.exact_control + counts.negr1:
        raise ContractError(
            f"{line} legacy controls contain labels beyond "
            f"{EXACT_CONTROL_LABEL!r} and {MISCLASSIFIED_TARGET!r}"
        )


def log1p_cpm_from_sum(
    count_sum: np.ndarray,
    *,
    target_sum: float = TARGET_SUM,
    dtype: np.dtype | type | None = None,
) -> np.ndarray:
    """Mirror the Orion donor builder's pseudobulk normalization."""

    values = np.asarray(count_sum, dtype=np.float64)
    if values.ndim != 1 or len(values) == 0:
        raise ContractError("count sum must be a nonempty one-dimensional array")
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ContractError("count sum must be finite and nonnegative")
    total = float(values.sum(dtype=np.float64))
    if not np.isfinite(total) or total <= 0:
        raise ContractError("count sum must have positive finite mass")
    if not np.isfinite(target_sum) or target_sum <= 0:
        raise ContractError("target_sum must be positive and finite")
    result = np.log1p(float(target_sum) * values / total)
    return result if dtype is None else result.astype(dtype)


def reconstruct_control_baselines(
    exact_control_sum: np.ndarray,
    negr1_sum: np.ndarray,
    *,
    target_sum: float = TARGET_SUM,
    dtype: np.dtype | type = np.float32,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the reconstructed legacy and corrected control baselines."""

    exact = np.asarray(exact_control_sum, dtype=np.float64)
    negr1 = np.asarray(negr1_sum, dtype=np.float64)
    if exact.shape != negr1.shape:
        raise ContractError("exact-control and NEGR1 sums have different shapes")
    corrected = log1p_cpm_from_sum(exact, target_sum=target_sum, dtype=dtype)
    legacy = log1p_cpm_from_sum(
        exact + negr1,
        target_sum=target_sum,
        dtype=dtype,
    )
    return legacy, corrected


def require_exact_legacy_baseline(
    existing_old_b_ctrl: np.ndarray,
    reconstructed_old_b_ctrl: np.ndarray,
) -> None:
    """Require a bit-identical reconstruction of the current donor baseline."""

    existing = np.asarray(existing_old_b_ctrl)
    reconstructed = np.asarray(reconstructed_old_b_ctrl)
    if existing.shape != reconstructed.shape or existing.dtype != reconstructed.dtype:
        raise ContractError("legacy baseline shape or dtype changed")
    if not np.array_equal(existing, reconstructed):
        maximum = float(
            np.max(
                np.abs(existing.astype(np.float64) - reconstructed.astype(np.float64))
            )
        )
        raise ContractError(
            "raw source sums do not reconstruct the existing legacy baseline "
            f"bit-for-bit; max_abs={maximum:.9g}"
        )


def translate_deltas(
    old_delta: np.ndarray,
    old_b_ctrl: np.ndarray,
    corrected_b_ctrl: np.ndarray,
    *,
    output_dtype: np.dtype | type | None = None,
) -> np.ndarray:
    """Apply the preregistered baseline translation to every target row."""

    delta = np.asarray(old_delta)
    old_baseline = np.asarray(old_b_ctrl)
    new_baseline = np.asarray(corrected_b_ctrl)
    if delta.ndim != 2:
        raise ContractError("old_delta must be a target-by-gene matrix")
    if old_baseline.ndim != 1 or new_baseline.ndim != 1:
        raise ContractError("control baselines must be one-dimensional")
    if old_baseline.shape != new_baseline.shape:
        raise ContractError("old and corrected baseline axes differ")
    if delta.shape[1] != len(old_baseline):
        raise ContractError("delta and baseline gene axes differ")
    for name, value in (
        ("old_delta", delta),
        ("old_b_ctrl", old_baseline),
        ("corrected_b_ctrl", new_baseline),
    ):
        if not np.issubdtype(value.dtype, np.number):
            raise ContractError(f"{name} is not numeric")
        if not np.isfinite(value).all():
            raise ContractError(f"{name} contains non-finite values")
    translated = (
        delta.astype(np.float64)
        + old_baseline.astype(np.float64)[None, :]
        - new_baseline.astype(np.float64)[None, :]
    )
    dtype = delta.dtype if output_dtype is None else np.dtype(output_dtype)
    return translated.astype(dtype)


def reconstruct_target_log_profiles(
    delta: np.ndarray,
    baseline: np.ndarray,
) -> np.ndarray:
    """Recover target log profiles from a delta and its control origin."""

    values = np.asarray(delta)
    origin = np.asarray(baseline)
    if values.ndim != 2 or origin.ndim != 1 or values.shape[1] != len(origin):
        raise ContractError("delta and baseline shapes are incompatible")
    return values.astype(np.float64) + origin.astype(np.float64)[None, :]


def validate_axis_identity(
    old_genes: Sequence[str],
    new_genes: Sequence[str],
    old_targets: Sequence[str],
    new_targets: Sequence[str],
) -> None:
    """Require byte-order-equivalent axes and unchanged target coverage."""

    old_gene_array = np.asarray(old_genes, dtype=str)
    new_gene_array = np.asarray(new_genes, dtype=str)
    old_target_array = np.asarray(old_targets, dtype=str)
    new_target_array = np.asarray(new_targets, dtype=str)
    if not np.array_equal(old_gene_array, new_gene_array):
        raise ContractError("gene axis or order changed")
    if not np.array_equal(old_target_array, new_target_array):
        raise ContractError("perturbation axis, order, or coverage changed")
    if (
        MISCLASSIFIED_TARGET in old_target_array
        or MISCLASSIFIED_TARGET in new_target_array
    ):
        raise ContractError("NEGR1 unexpectedly appears in an Orion donor target axis")


def _is_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == SHA256_LENGTH
        and all(character in "0123456789abcdef" for character in value)
    )


def _require_keys(
    value: Mapping[str, object],
    expected: set[str],
    *,
    label: str,
) -> None:
    observed = set(value)
    if observed != expected:
        raise ContractError(
            f"{label} keys changed; missing={sorted(expected - observed)}, "
            f"extra={sorted(observed - expected)}"
        )


def validate_source_metadata_audit(audit: Mapping[str, object]) -> None:
    """Validate the future metadata-only provenance/count report."""

    _require_keys(
        audit,
        {
            "schema_version",
            "experiment",
            "stage",
            "expression_column_read",
            "columns_read",
            "challenge_panel",
            "lines",
        },
        label="source audit",
    )
    if audit["schema_version"] != SCHEMA_VERSION or audit["experiment"] != EXPERIMENT:
        raise ContractError("source audit identity changed")
    if audit["stage"] != "metadata-only-source-audit":
        raise ContractError("source audit stage changed")
    if audit["expression_column_read"] is not False:
        raise ContractError("metadata audit opened gene_expression")
    if tuple(audit["columns_read"]) != METADATA_ONLY_COLUMNS:
        raise ContractError("metadata audit projected unauthorized columns")

    panel = audit["challenge_panel"]
    if not isinstance(panel, Mapping):
        raise ContractError("challenge_panel must be an object")
    _require_keys(panel, {"sha256", "contains_negr1"}, label="challenge panel")
    if not _is_sha256(panel["sha256"]) or panel["contains_negr1"] is not False:
        raise ContractError("challenge panel provenance or NEGR1 exclusion failed")

    lines = audit["lines"]
    if not isinstance(lines, Mapping) or set(lines) != set(EXPECTED_LINES):
        raise ContractError("source audit Orion line set changed")
    for line, expected in EXPECTED_LINES.items():
        entry = lines[line]
        if not isinstance(entry, Mapping):
            raise ContractError(f"{line} source audit must be an object")
        _require_keys(
            entry,
            {
                "raw_parquet_files",
                "row_groups",
                "accepted_label_counts",
                "negr1_false_control_count",
                "legacy_control_count",
                "donor_perturbations_contain_negr1",
                "legacy_artifacts",
            },
            label=f"{line} source audit",
        )
        if entry["row_groups"] != expected["row_groups"]:
            raise ContractError(f"{line} raw Parquet row-group count changed")
        if entry["accepted_label_counts"] != {
            EXACT_CONTROL_LABEL: expected["exact_control_count"]
        }:
            raise ContractError(f"{line} exact control-label counts changed")
        counts = ControlCounts(
            exact_control=expected["exact_control_count"],
            negr1=int(entry["negr1_false_control_count"]),
            legacy_control=int(entry["legacy_control_count"]),
        )
        validate_control_counts(line, counts)
        if entry["donor_perturbations_contain_negr1"] is not False:
            raise ContractError(f"{line} donor perturbation axis contains NEGR1")

        raw_files = entry["raw_parquet_files"]
        if not isinstance(raw_files, list) or not raw_files:
            raise ContractError(f"{line} raw Parquet provenance is empty")
        paths: list[str] = []
        row_groups = 0
        for item in raw_files:
            if not isinstance(item, Mapping):
                raise ContractError(f"{line} raw Parquet record must be an object")
            _require_keys(
                item,
                {"path", "size_bytes", "row_groups", "sha256"},
                label=f"{line} raw Parquet record",
            )
            path = item["path"]
            if not isinstance(path, str) or not path.endswith(".parquet"):
                raise ContractError(f"{line} raw Parquet path is invalid")
            if not isinstance(item["size_bytes"], int) or item["size_bytes"] <= 0:
                raise ContractError(f"{line} raw Parquet size is invalid")
            if not isinstance(item["row_groups"], int) or item["row_groups"] <= 0:
                raise ContractError(f"{line} raw Parquet row_groups is invalid")
            if not _is_sha256(item["sha256"]):
                raise ContractError(f"{line} raw Parquet SHA-256 is invalid")
            paths.append(path)
            row_groups += int(item["row_groups"])
        if paths != sorted(paths) or len(paths) != len(set(paths)):
            raise ContractError(f"{line} raw Parquet paths are not sorted and unique")
        if row_groups != expected["row_groups"]:
            raise ContractError(f"{line} per-file row groups do not match the total")

        artifacts = entry["legacy_artifacts"]
        if not isinstance(artifacts, Mapping):
            raise ContractError(f"{line} legacy_artifacts must be an object")
        _require_keys(
            artifacts,
            set(LEGACY_ARTIFACT_NAMES),
            label=f"{line} legacy artifacts",
        )
        if not all(_is_sha256(value) for value in artifacts.values()):
            raise ContractError(f"{line} legacy artifact SHA-256 is invalid")
