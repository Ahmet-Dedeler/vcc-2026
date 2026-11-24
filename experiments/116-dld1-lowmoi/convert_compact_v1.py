#!/usr/bin/env python3
"""Build the frozen source-only compact canonical artifact for Experiment 116."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import tempfile
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np

import acquire_lfc_v1 as acquisition


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "CONVERSION_PROTOCOL.v1.json"
FREEZE_PATH = HERE / "PREVALUE_FREEZE.v1.json"
ACCESS_LOCK_PATH = HERE / "MATRIX_VALUE_ACCESS_LOCK.v1.json"
ACQUISITION_RECEIPT_PATH = HERE / "ACQUISITION_RECEIPT.v1.json"
WORK_ROOT = HERE / "work"
ARTIFACT_ROOT = WORK_ROOT / "compact-source-v1"
CONFIRMATION = "I_CONFIRM_EXPERIMENT_116_SOURCE_ONLY_CONVERSION_V1"
MISSING_TOKENS = frozenset({"", "NA", "NaN", "nan"})
BLOCK_ROWS = 128


class ConversionError(RuntimeError):
    """The frozen conversion contract or source data failed closed."""


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
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


def load_protocol(path: Path = PROTOCOL_PATH) -> dict[str, Any]:
    protocol = json.loads(path.read_text(encoding="utf-8"))
    if protocol.get("schema_version") != 1:
        raise ConversionError("unsupported conversion protocol schema")
    if protocol.get("experiment") != "116-dld1-lowmoi":
        raise ConversionError("wrong experiment in conversion protocol")
    if protocol.get("status") != "frozen_before_matrix_value_access":
        raise ConversionError("conversion protocol status changed")
    if protocol.get("scope") != "source_only_compact_canonical_artifact":
        raise ConversionError("conversion scope changed")
    forbidden = protocol.get("forbidden", {})
    if not forbidden or not all(value is True for value in forbidden.values()):
        raise ConversionError("all forbidden actions must remain forbidden")
    policy = protocol["transfer_arrays"]["panel_response_coordinate_policy"]
    if not policy.startswith("All 300 VCC panel response coordinates"):
        raise ConversionError("panel response exclusion policy changed")
    if protocol["duplicate_promoter_policy"]["selection_by_observed_effect"]:
        raise ConversionError("effect-based promoter selection is forbidden")
    return protocol


def verify_prevalue_freeze() -> dict[str, Any]:
    try:
        return acquisition.verify_prevalue_freeze(FREEZE_PATH)
    except acquisition.AcquisitionError as error:
        raise ConversionError(str(error)) from error


def read_axis(spec: dict[str, Any]) -> list[str]:
    path = ROOT / spec["path"]
    if path.is_symlink() or not path.is_file():
        raise ConversionError(f"missing or symlinked axis: {path}")
    if path.stat().st_size <= 0 or sha256_file(path) != spec["sha256"]:
        raise ConversionError(f"axis hash drift: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = [row[spec["column"]].strip() for row in rows]
    if len(values) != spec["row_count"] or len(values) != len(set(values)):
        raise ConversionError(f"axis row-count or uniqueness drift: {path}")
    if any(not value for value in values):
        raise ConversionError(f"axis contains an empty symbol: {path}")
    return values


def read_target_labels(protocol: dict[str, Any]) -> list[str]:
    spec = protocol["inputs"]
    path = ROOT / spec["target_labels_path"]
    if path.is_symlink() or not path.is_file():
        raise ConversionError("frozen target-label file is missing")
    if sha256_file(path) != spec["target_labels_sha256"]:
        raise ConversionError("target-label hash drift")
    labels = path.read_text(encoding="utf-8").splitlines()
    if len(labels) != spec["target_label_count"] or len(labels) != len(set(labels)):
        raise ConversionError("target-label count or uniqueness drift")
    return labels


def map_promoter_labels(
    labels: list[str],
    canonical_genes: list[str],
) -> tuple[list[str | None], dict[str, list[int]], dict[str, int]]:
    canonical = set(canonical_genes)
    mapped: list[str | None] = []
    groups: dict[str, list[int]] = defaultdict(list)
    counts = {"ntc": 0, "mapped": 0, "unsupported": 0}
    for index, label in enumerate(labels):
        if label.startswith("NTC_"):
            mapped.append(None)
            counts["ntc"] += 1
            continue
        symbol = label.split("_", 1)[0] if "_" in label else label
        if symbol in canonical:
            mapped.append(symbol)
            groups[symbol].append(index)
            counts["mapped"] += 1
        else:
            mapped.append(None)
            counts["unsupported"] += 1
    return mapped, dict(groups), counts


def validate_metadata_mapping(
    protocol: dict[str, Any],
    labels: list[str],
    canonical_genes: list[str],
    groups: dict[str, list[int]],
    counts: dict[str, int],
    panel: list[str],
) -> None:
    expected = protocol["predeclared_mapping"]["expected_metadata_counts"]
    actual = {
        "all_promoter_labels": len(labels),
        "ntc_labels": counts["ntc"],
        "non_ntc_labels": len(labels) - counts["ntc"],
        "canonical_mapped_promoter_labels": counts["mapped"],
        "unsupported_promoter_labels": counts["unsupported"],
        "canonical_mapped_target_symbols": len(groups),
        "canonical_targets_with_multiple_promoters": sum(
            len(indices) > 1 for indices in groups.values()
        ),
        "vcc_panel_target_overlap": len(set(groups) & set(panel)),
        "vcc_panel_promoter_overlap": sum(
            len(groups[symbol]) for symbol in set(groups) & set(panel)
        ),
    }
    if actual != expected:
        raise ConversionError(f"pre-value metadata mapping drift: {actual!r}")
    if any("_" in gene for gene in canonical_genes):
        raise ConversionError("canonical gene symbols unexpectedly contain underscores")


def parse_value(token: str) -> tuple[float, bool]:
    if token in MISSING_TOKENS:
        return math.nan, False
    try:
        value = float(token)
    except ValueError as error:
        raise ConversionError(f"invalid numeric token: {token!r}") from error
    if not math.isfinite(value):
        raise ConversionError(f"non-finite numeric token: {token!r}")
    return value, True


def parse_values(tokens: Iterable[str], expected_count: int) -> np.ndarray:
    values = np.fromiter(
        (parse_value(token)[0] for token in tokens),
        dtype=np.float64,
        count=expected_count,
    )
    if values.size != expected_count:
        raise ConversionError("source row has an unexpected value count")
    return values


def aggregation_plan(
    target_symbols: list[str],
    groups: dict[str, list[int]],
) -> tuple[np.ndarray, np.ndarray, list[tuple[int, np.ndarray]]]:
    singleton_targets: list[int] = []
    singleton_columns: list[int] = []
    duplicates: list[tuple[int, np.ndarray]] = []
    for target_index, symbol in enumerate(target_symbols):
        columns = groups[symbol]
        if len(columns) == 1:
            singleton_targets.append(target_index)
            singleton_columns.append(columns[0])
        else:
            duplicates.append(
                (target_index, np.asarray(columns, dtype=np.int64))
            )
    return (
        np.asarray(singleton_targets, dtype=np.int64),
        np.asarray(singleton_columns, dtype=np.int64),
        duplicates,
    )


def aggregate_block(
    values: np.ndarray,
    singleton_targets: np.ndarray,
    singleton_columns: np.ndarray,
    duplicate_groups: list[tuple[int, np.ndarray]],
    target_count: int,
) -> np.ndarray:
    aggregated = np.full((values.shape[0], target_count), np.nan, dtype=np.float64)
    aggregated[:, singleton_targets] = values[:, singleton_columns]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for target_index, columns in duplicate_groups:
            aggregated[:, target_index] = np.nanmedian(values[:, columns], axis=1)
    return aggregated


def apply_transfer_exclusion(
    canonical_indices: np.ndarray,
    aggregated: np.ndarray,
    panel_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    values = np.nan_to_num(aggregated, nan=0.0).astype(np.float32)
    measured = np.isfinite(aggregated)
    excluded_rows = panel_mask[canonical_indices]
    if np.any(excluded_rows):
        values[excluded_rows, :] = 0.0
        measured[excluded_rows, :] = False
    return values, measured


def _write_lines(path: Path, values: Iterable[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in values), encoding="utf-8")


def stream_split(
    *,
    split: str,
    source_path: Path,
    expected_source_sha256: str,
    expected_labels: list[str],
    canonical_genes: list[str],
    panel_mask: np.ndarray,
    target_symbols: list[str],
    singleton_targets: np.ndarray,
    singleton_columns: np.ndarray,
    duplicate_groups: list[tuple[int, np.ndarray]],
    author_h5: h5py.File,
    artifact_root: Path,
) -> dict[str, Any]:
    if source_path.is_symlink() or not source_path.is_file():
        raise ConversionError(f"{split}: accepted source is missing or symlinked")
    if sha256_file(source_path) != expected_source_sha256:
        raise ConversionError(f"{split}: accepted source SHA-256 drift")

    n_genes = len(canonical_genes)
    n_promoters = len(expected_labels)
    n_targets = len(target_symbols)
    group = author_h5.create_group(split)
    author_lfc = group.create_dataset(
        "lfc",
        shape=(n_genes, n_promoters),
        dtype=np.float64,
        chunks=(128, 64),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
        fillvalue=np.nan,
    )
    author_measured = group.create_dataset(
        "measured",
        shape=(n_genes, n_promoters),
        dtype=np.bool_,
        chunks=(128, 64),
        compression="gzip",
        compression_opts=4,
        shuffle=True,
        fillvalue=False,
    )
    transfer_lfc_path = artifact_root / f"transfer_lfc_{split}.npy"
    transfer_mask_path = artifact_root / f"transfer_measured_{split}.npy"
    transfer_lfc = np.lib.format.open_memmap(
        transfer_lfc_path,
        mode="w+",
        dtype=np.float32,
        shape=(n_genes, n_targets),
    )
    transfer_lfc[:] = 0.0
    transfer_measured = np.lib.format.open_memmap(
        transfer_mask_path,
        mode="w+",
        dtype=np.bool_,
        shape=(n_genes, n_targets),
    )
    transfer_measured[:] = False

    canonical_index = {symbol: index for index, symbol in enumerate(canonical_genes)}
    seen_source_labels: set[str] = set()
    source_label_digest = hashlib.sha256()
    unsupported_label_digest = hashlib.sha256()
    unsupported_count = 0
    source_rows = 0
    canonical_rows = 0
    buffered_indices: list[int] = []
    buffered_values: list[np.ndarray] = []

    def flush() -> None:
        nonlocal canonical_rows
        if not buffered_indices:
            return
        indices = np.asarray(buffered_indices, dtype=np.int64)
        block = np.vstack(buffered_values)
        order = np.argsort(indices)
        indices = indices[order]
        block = block[order, :]
        if len(np.unique(indices)) != len(indices):
            raise ConversionError(f"{split}: duplicate canonical response row")
        author_lfc[indices, :] = block
        author_measured[indices, :] = np.isfinite(block)
        aggregated = aggregate_block(
            block,
            singleton_targets,
            singleton_columns,
            duplicate_groups,
            n_targets,
        )
        values, measured = apply_transfer_exclusion(
            indices,
            aggregated,
            panel_mask,
        )
        transfer_lfc[indices, :] = values
        transfer_measured[indices, :] = measured
        canonical_rows += len(indices)
        buffered_indices.clear()
        buffered_values.clear()

    csv.field_size_limit(16 << 20)
    with gzip.open(source_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ConversionError(f"{split}: empty source matrix") from error
        if not header or header[0] != "gene" or header[1:] != expected_labels:
            raise ConversionError(f"{split}: source header drift")
        for row in reader:
            source_rows += 1
            if len(row) != n_promoters + 1:
                raise ConversionError(
                    f"{split}: row {source_rows} width {len(row)} != {n_promoters + 1}"
                )
            label = row[0]
            if not label or label in seen_source_labels:
                raise ConversionError(f"{split}: empty or duplicate source row {label!r}")
            seen_source_labels.add(label)
            source_label_digest.update(label.encode("utf-8"))
            source_label_digest.update(b"\n")
            values = parse_values(row[1:], n_promoters)
            index = canonical_index.get(label)
            if index is None:
                unsupported_count += 1
                unsupported_label_digest.update(label.encode("utf-8"))
                unsupported_label_digest.update(b"\n")
                continue
            buffered_indices.append(index)
            buffered_values.append(values)
            if len(buffered_indices) >= BLOCK_ROWS:
                flush()
        flush()

    transfer_lfc.flush()
    transfer_measured.flush()
    transfer_lfc._mmap.close()
    transfer_measured._mmap.close()
    return {
        "split": split,
        "source_path": str(source_path.relative_to(ROOT)),
        "source_sha256": expected_source_sha256,
        "source_rows": source_rows,
        "source_row_labels_sha256": source_label_digest.hexdigest(),
        "canonical_response_rows": canonical_rows,
        "unsupported_response_rows": unsupported_count,
        "unsupported_response_labels_sha256": unsupported_label_digest.hexdigest(),
    }


def vector_metrics(
    low1: np.ndarray,
    low2: np.ndarray,
    measured1: np.ndarray,
    measured2: np.ndarray,
    *,
    top_k: int = 100,
) -> dict[str, float | int]:
    common = measured1 & measured2
    count = int(common.sum())
    if count == 0:
        return {
            "common_measured": 0,
            "cosine": math.nan,
            "pearson": math.nan,
            "top100_sign_agreement": math.nan,
            "normalized_rmse": math.nan,
        }
    a = np.asarray(low1[common], dtype=np.float64)
    b = np.asarray(low2[common], dtype=np.float64)
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    cosine = (
        float(np.dot(a, b) / (norm_a * norm_b))
        if norm_a > 0.0 and norm_b > 0.0
        else math.nan
    )
    centered_a = a - a.mean()
    centered_b = b - b.mean()
    centered_norm = float(np.linalg.norm(centered_a) * np.linalg.norm(centered_b))
    pearson = (
        float(np.dot(centered_a, centered_b) / centered_norm)
        if centered_norm > 0.0
        else math.nan
    )
    k = min(top_k, count)
    importance = (np.abs(a) + np.abs(b)) / 2.0
    selected = np.argpartition(importance, -k)[-k:] if k < count else np.arange(count)
    sign_agreement = float(np.mean(np.sign(a[selected]) == np.sign(b[selected])))
    denominator = math.sqrt(
        max((float(np.mean(a * a)) + float(np.mean(b * b))) / 2.0, 1e-30)
    )
    normalized_rmse = float(math.sqrt(float(np.mean((a - b) ** 2))) / denominator)
    return {
        "common_measured": count,
        "cosine": cosine,
        "pearson": pearson,
        "top100_sign_agreement": sign_agreement,
        "normalized_rmse": normalized_rmse,
    }


def pairwise_masked_cosine(
    query_values: np.ndarray,
    gallery_values: np.ndarray,
    query_measured: np.ndarray,
    gallery_measured: np.ndarray,
) -> np.ndarray:
    n_query = query_values.shape[0]
    n_gallery = gallery_values.shape[0]
    result = np.full((n_query, n_gallery), np.nan, dtype=np.float64)
    for i in range(n_query):
        for j in range(n_gallery):
            common = query_measured[i] & gallery_measured[j]
            if not np.any(common):
                continue
            a = query_values[i, common].astype(np.float64, copy=False)
            b = gallery_values[j, common].astype(np.float64, copy=False)
            denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
            if denominator > 0.0:
                result[i, j] = float(np.dot(a, b) / denominator)
    return result


def pds_scores(similarity: np.ndarray, mapping: np.ndarray | None = None) -> np.ndarray:
    rows, columns = similarity.shape
    if rows != columns or rows < 2:
        raise ConversionError("PDS requires a square gallery with at least two targets")
    if mapping is None:
        mapping = np.arange(rows)
    scores = np.empty(rows, dtype=np.float64)
    for row in range(rows):
        own_column = int(mapping[row])
        own = similarity[row, own_column]
        distractors = np.delete(similarity[row, :], own_column)
        if not np.isfinite(own) or not np.all(np.isfinite(distractors)):
            raise ConversionError("non-finite pairwise cosine in PDS gallery")
        scores[row] = (
            float(np.sum(own > distractors))
            + 0.5 * float(np.sum(own == distractors))
        ) / len(distractors)
    return scores


def _finite_summary(values: np.ndarray) -> dict[str, float | int]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"count": 0, "median": math.nan, "mean": math.nan}
    return {
        "count": int(finite.size),
        "median": float(np.median(finite)),
        "mean": float(np.mean(finite)),
        "p10": float(np.quantile(finite, 0.10)),
        "p90": float(np.quantile(finite, 0.90)),
    }


def compute_diagnostics(
    *,
    artifact_root: Path,
    protocol: dict[str, Any],
    canonical_genes: list[str],
    panel: list[str],
    labels: list[str],
    groups: dict[str, list[int]],
    target_symbols: list[str],
) -> dict[str, Any]:
    panel_set = set(panel)
    panel_mask = np.asarray([gene in panel_set for gene in canonical_genes])
    nonpanel_count = int((~panel_mask).sum())
    low1 = np.load(artifact_root / "transfer_lfc_low1.npy", mmap_mode="r")
    low2 = np.load(artifact_root / "transfer_lfc_low2.npy", mmap_mode="r")
    measured1 = np.load(
        artifact_root / "transfer_measured_low1.npy",
        mmap_mode="r",
    )
    measured2 = np.load(
        artifact_root / "transfer_measured_low2.npy",
        mmap_mode="r",
    )
    expected_shape = (len(canonical_genes), len(target_symbols))
    for array in (low1, low2, measured1, measured2):
        if array.shape != expected_shape:
            raise ConversionError("transfer array shape drift")
    for array in (low1, low2):
        if not np.all(array[panel_mask, :] == 0.0):
            raise ConversionError("panel response values are not exactly zero")
    for array in (measured1, measured2):
        if np.any(array[panel_mask, :]):
            raise ConversionError("panel response mask is not exactly false")

    thresholds = protocol["qc_policy"]
    target_rows: list[dict[str, Any]] = []
    cosine_values = np.empty(len(target_symbols), dtype=np.float64)
    reliable = np.zeros(len(target_symbols), dtype=np.bool_)
    for target_index, symbol in enumerate(target_symbols):
        metrics = vector_metrics(
            low1[:, target_index],
            low2[:, target_index],
            measured1[:, target_index],
            measured2[:, target_index],
        )
        fraction = metrics["common_measured"] / nonpanel_count
        is_reliable = (
            fraction
            >= thresholds["minimum_common_nonpanel_measured_fraction"]
            and math.isfinite(float(metrics["cosine"]))
            and metrics["cosine"] >= thresholds["minimum_low1_low2_cosine"]
            and metrics["top100_sign_agreement"]
            >= thresholds["minimum_top100_sign_agreement"]
        )
        cosine_values[target_index] = float(metrics["cosine"])
        reliable[target_index] = is_reliable
        target_rows.append(
            {
                "target_symbol": symbol,
                "promoter_count": len(groups[symbol]),
                "vcc_panel_target": symbol in panel_set,
                "common_measured": metrics["common_measured"],
                "common_measured_fraction": fraction,
                "cosine": metrics["cosine"],
                "pearson": metrics["pearson"],
                "top100_sign_agreement": metrics["top100_sign_agreement"],
                "normalized_rmse": metrics["normalized_rmse"],
                "reliable": bool(is_reliable),
            }
        )

    with (artifact_root / "target_reliability.csv").open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(target_rows[0]))
        writer.writeheader()
        writer.writerows(target_rows)
    np.save(artifact_root / "reliable_target_mask.npy", reliable)

    duplicate_pairs: list[tuple[str, int, int]] = []
    for symbol in target_symbols:
        columns = groups[symbol]
        for left_offset, left in enumerate(columns):
            for right in columns[left_offset + 1 :]:
                duplicate_pairs.append((symbol, left, right))
    promoter_rows: list[dict[str, Any]] = []
    if duplicate_pairs:
        with h5py.File(artifact_root / "author_lfc.h5", "r") as handle:
            for symbol, left, right in duplicate_pairs:
                row: dict[str, Any] = {
                    "target_symbol": symbol,
                    "left_promoter": labels[left],
                    "right_promoter": labels[right],
                }
                for split in ("low1", "low2"):
                    values_left = handle[f"{split}/lfc"][:, left]
                    values_right = handle[f"{split}/lfc"][:, right]
                    mask_left = handle[f"{split}/measured"][:, left]
                    mask_right = handle[f"{split}/measured"][:, right]
                    metrics = vector_metrics(
                        values_left,
                        values_right,
                        mask_left & ~panel_mask,
                        mask_right & ~panel_mask,
                    )
                    row[f"{split}_common_measured"] = metrics["common_measured"]
                    row[f"{split}_cosine"] = metrics["cosine"]
                promoter_rows.append(row)
        with (artifact_root / "duplicate_promoter_reliability.csv").open(
            "w",
            encoding="utf-8",
            newline="",
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(promoter_rows[0]))
            writer.writeheader()
            writer.writerows(promoter_rows)

    overlap_symbols = [symbol for symbol in target_symbols if symbol in panel_set]
    overlap_indices = np.asarray(
        [target_symbols.index(symbol) for symbol in overlap_symbols],
        dtype=np.int64,
    )
    query1 = np.asarray(low1[:, overlap_indices].T)
    query2 = np.asarray(low2[:, overlap_indices].T)
    mask1 = np.asarray(measured1[:, overlap_indices].T)
    mask2 = np.asarray(measured2[:, overlap_indices].T)
    similarity12 = pairwise_masked_cosine(query1, query2, mask1, mask2)
    similarity21 = pairwise_masked_cosine(query2, query1, mask2, mask1)
    scores12 = pds_scores(similarity12)
    scores21 = pds_scores(similarity21)
    per_target_macro = (scores12 + scores21) / 2.0
    directional12 = float(scores12.mean())
    directional21 = float(scores21.mean())
    macro_pds = float(per_target_macro.mean())

    diagnostic_spec = protocol["source_only_diagnostics"]["vcc_overlap_retrieval"]
    rng = np.random.default_rng(diagnostic_spec["seed"])
    null_values = np.empty(diagnostic_spec["permutations"], dtype=np.float64)
    for index in range(null_values.size):
        permutation = rng.permutation(len(overlap_symbols))
        null_values[index] = (
            pds_scores(similarity12, permutation).mean()
            + pds_scores(similarity21, np.argsort(permutation)).mean()
        ) / 2.0
    p_value = float((1 + np.sum(null_values >= macro_pds)) / (1 + len(null_values)))
    null_p95 = float(np.quantile(null_values, 0.95))
    bootstrap = np.empty(
        diagnostic_spec["bootstrap_replicates"],
        dtype=np.float64,
    )
    for index in range(bootstrap.size):
        sample = rng.integers(0, len(per_target_macro), len(per_target_macro))
        bootstrap[index] = float(np.mean(per_target_macro[sample]))
    bootstrap_low, bootstrap_high = (
        float(value) for value in np.quantile(bootstrap, [0.025, 0.975])
    )

    gate = protocol["source_only_diagnostics"][
        "preregistered_held_experiment_recommendation_gate"
    ]
    reliable_fraction = float(np.mean(reliable))
    finite_cosine = cosine_values[np.isfinite(cosine_values)]
    median_target_cosine = (
        float(np.median(finite_cosine)) if finite_cosine.size else math.nan
    )
    gate_checks = {
        "vcc_overlap_targets": len(overlap_symbols)
        >= gate["minimum_vcc_overlap_targets"],
        "macro_replicate_pds": macro_pds >= gate["minimum_macro_replicate_pds"],
        "low1_to_low2_pds": directional12
        >= gate["minimum_each_direction_pds"],
        "low2_to_low1_pds": directional21
        >= gate["minimum_each_direction_pds"],
        "bootstrap_lower_95": bootstrap_low
        >= gate["minimum_bootstrap_lower_95"],
        "permutation_p_value": p_value <= gate["maximum_permutation_p_value"],
        "macro_above_null_p95": macro_pds > null_p95,
        "median_target_cosine": median_target_cosine
        >= gate["minimum_median_target_cosine"],
        "reliable_target_fraction": reliable_fraction
        >= gate["minimum_reliable_target_fraction"],
    }

    return {
        "target_level": {
            "target_count": len(target_symbols),
            "nonpanel_response_coordinate_count": nonpanel_count,
            "reliable_target_count": int(reliable.sum()),
            "reliable_target_fraction": reliable_fraction,
            "split_cosine": _finite_summary(cosine_values),
            "median_target_cosine": median_target_cosine,
        },
        "duplicate_promoters": {
            "target_count": sum(len(groups[symbol]) > 1 for symbol in target_symbols),
            "pair_count": len(duplicate_pairs),
            "low1_pairwise_cosine": _finite_summary(
                np.asarray(
                    [row["low1_cosine"] for row in promoter_rows],
                    dtype=np.float64,
                )
            ),
            "low2_pairwise_cosine": _finite_summary(
                np.asarray(
                    [row["low2_cosine"] for row in promoter_rows],
                    dtype=np.float64,
                )
            ),
        },
        "vcc_overlap_retrieval": {
            "target_count": len(overlap_symbols),
            "target_symbols": overlap_symbols,
            "low1_to_low2_pds": directional12,
            "low2_to_low1_pds": directional21,
            "macro_pds": macro_pds,
            "bootstrap_95": [bootstrap_low, bootstrap_high],
            "null_p95": null_p95,
            "permutation_p_value": p_value,
            "permutations": diagnostic_spec["permutations"],
            "bootstrap_replicates": diagnostic_spec["bootstrap_replicates"],
            "seed": diagnostic_spec["seed"],
        },
        "recommend_separately_preregistered_held_context_experiment": all(
            gate_checks.values()
        ),
        "recommendation_gate_checks": gate_checks,
    }


def build_artifact_manifest(
    artifact_root: Path,
    *,
    protocol_sha256: str,
    freeze_sha256: str,
    acquisition_receipt_sha256: str,
    build_summary: dict[str, Any],
    diagnostics: dict[str, Any],
) -> dict[str, Any]:
    files: dict[str, Any] = {}
    for path in sorted(artifact_root.iterdir()):
        if path.name == "artifact-manifest.json" or not path.is_file():
            continue
        files[path.name] = {
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return {
        "schema_version": 1,
        "experiment": "116-dld1-lowmoi",
        "dataset": "GSE337988",
        "status": "source_only_compact_artifact_complete",
        "protocol_sha256": protocol_sha256,
        "prevalue_freeze_sha256": freeze_sha256,
        "acquisition_receipt_sha256": acquisition_receipt_sha256,
        "build_summary": build_summary,
        "diagnostics": diagnostics,
        "files": files,
        "held_truth_opened": False,
        "existing_outcomes_opened": False,
        "cell_eval2_run": False,
        "counts_generated": False,
        "submission_built_or_sent": False,
    }


def verify_artifact(path: Path = ARTIFACT_ROOT) -> dict[str, Any]:
    manifest_path = path / "artifact-manifest.json"
    if not manifest_path.is_file():
        raise ConversionError("artifact manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "source_only_compact_artifact_complete":
        raise ConversionError("artifact status changed")
    for name, expected in manifest["files"].items():
        candidate = path / name
        if (
            candidate.is_symlink()
            or not candidate.is_file()
            or candidate.stat().st_size != expected["size_bytes"]
            or sha256_file(candidate) != expected["sha256"]
        ):
            raise ConversionError(f"artifact file drift: {name}")
    return manifest


def ensure_access_lock(receipt: dict[str, Any]) -> dict[str, Any]:
    payload = {
        "schema_version": 1,
        "experiment": "116-dld1-lowmoi",
        "status": "matrix_value_access_authorized_from_frozen_protocol",
        "created_at": "2026-08-30",
        "prevalue_freeze_sha256": sha256_file(FREEZE_PATH),
        "conversion_protocol_sha256": sha256_file(PROTOCOL_PATH),
        "acquisition_receipt_sha256": sha256_file(ACQUISITION_RECEIPT_PATH),
        "source_assets": {
            name: {
                "path": receipt["assets"][name]["path"],
                "size_bytes": receipt["assets"][name]["size_bytes"],
                "sha256": receipt["assets"][name]["sha256"],
            }
            for name in acquisition.ASSET_NAMES
        },
    }
    if ACCESS_LOCK_PATH.exists():
        existing = json.loads(ACCESS_LOCK_PATH.read_text(encoding="utf-8"))
        if existing != payload:
            raise ConversionError("matrix-value access lock drift")
    else:
        atomic_json(ACCESS_LOCK_PATH, payload)
    return payload


def build() -> dict[str, Any]:
    protocol = load_protocol()
    verify_prevalue_freeze()
    amendment = acquisition.load_amendment()
    if not ACQUISITION_RECEIPT_PATH.is_file():
        raise ConversionError("accepted acquisition receipt is missing")
    receipt = json.loads(ACQUISITION_RECEIPT_PATH.read_text(encoding="utf-8"))
    if not acquisition._existing_receipt_is_valid(receipt, amendment):
        raise ConversionError("accepted acquisition receipt failed verification")
    ensure_access_lock(receipt)

    if ARTIFACT_ROOT.exists():
        return verify_artifact()
    canonical_genes = read_axis(protocol["inputs"]["canonical_gene_axis"])
    panel = read_axis(protocol["inputs"]["vcc_panel_axis"])
    labels = read_target_labels(protocol)
    mapped, groups, mapping_counts = map_promoter_labels(labels, canonical_genes)
    validate_metadata_mapping(
        protocol,
        labels,
        canonical_genes,
        groups,
        mapping_counts,
        panel,
    )
    target_symbols = [gene for gene in canonical_genes if gene in groups]
    singleton_targets, singleton_columns, duplicate_groups = aggregation_plan(
        target_symbols,
        groups,
    )
    panel_set = set(panel)
    panel_mask = np.asarray([gene in panel_set for gene in canonical_genes])

    WORK_ROOT.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(
        tempfile.mkdtemp(prefix=".compact-source-v1.", dir=WORK_ROOT)
    )
    try:
        _write_lines(temporary_root / "canonical_genes.txt", canonical_genes)
        _write_lines(temporary_root / "promoter_labels.txt", labels)
        _write_lines(
            temporary_root / "promoter_target_symbols.txt",
            [symbol or "" for symbol in mapped],
        )
        _write_lines(temporary_root / "target_symbols.txt", target_symbols)
        np.save(temporary_root / "panel_response_mask.npy", panel_mask)
        (temporary_root / "target_promoters.json").write_text(
            json.dumps(
                {
                    symbol: [labels[index] for index in groups[symbol]]
                    for symbol in target_symbols
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

        split_summaries: dict[str, Any] = {}
        with h5py.File(temporary_root / "author_lfc.h5", "w") as author_h5:
            author_h5.attrs["experiment"] = "116-dld1-lowmoi"
            author_h5.attrs["orientation"] = (
                "canonical_response_gene_by_original_promoter_label"
            )
            author_h5.attrs["protocol_sha256"] = sha256_file(PROTOCOL_PATH)
            for name, split in (
                ("low1_lfc_matrix", "low1"),
                ("low2_lfc_matrix", "low2"),
            ):
                item = receipt["assets"][name]
                split_summaries[split] = stream_split(
                    split=split,
                    source_path=ROOT / item["path"],
                    expected_source_sha256=item["sha256"],
                    expected_labels=labels,
                    canonical_genes=canonical_genes,
                    panel_mask=panel_mask,
                    target_symbols=target_symbols,
                    singleton_targets=singleton_targets,
                    singleton_columns=singleton_columns,
                    duplicate_groups=duplicate_groups,
                    author_h5=author_h5,
                    artifact_root=temporary_root,
                )
        if (
            split_summaries["low1"]["source_rows"]
            != split_summaries["low2"]["source_rows"]
            or split_summaries["low1"]["source_row_labels_sha256"]
            != split_summaries["low2"]["source_row_labels_sha256"]
        ):
            raise ConversionError("Low1/Low2 response-row sequence differs")

        diagnostics = compute_diagnostics(
            artifact_root=temporary_root,
            protocol=protocol,
            canonical_genes=canonical_genes,
            panel=panel,
            labels=labels,
            groups=groups,
            target_symbols=target_symbols,
        )
        build_summary = {
            "canonical_gene_count": len(canonical_genes),
            "panel_response_coordinate_count": int(panel_mask.sum()),
            "original_promoter_label_count": len(labels),
            "mapped_promoter_label_count": mapping_counts["mapped"],
            "unsupported_promoter_label_count": mapping_counts["unsupported"],
            "ntc_label_count": mapping_counts["ntc"],
            "canonical_target_count": len(target_symbols),
            "duplicate_promoter_target_count": len(duplicate_groups),
            "split_summaries": split_summaries,
        }
        (temporary_root / "result.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "experiment": "116-dld1-lowmoi",
                    "status": "source_only_diagnostics_complete",
                    "build_summary": build_summary,
                    "diagnostics": diagnostics,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        manifest = build_artifact_manifest(
            temporary_root,
            protocol_sha256=sha256_file(PROTOCOL_PATH),
            freeze_sha256=sha256_file(FREEZE_PATH),
            acquisition_receipt_sha256=sha256_file(ACQUISITION_RECEIPT_PATH),
            build_summary=build_summary,
            diagnostics=diagnostics,
        )
        (temporary_root / "artifact-manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary_root, ARTIFACT_ROOT)
    except BaseException:
        shutil.rmtree(temporary_root, ignore_errors=True)
        raise
    return verify_artifact()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("authenticate")
    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("--confirm", required=True)
    subparsers.add_parser("verify")
    args = parser.parse_args()

    if args.command == "authenticate":
        protocol = load_protocol()
        freeze = verify_prevalue_freeze()
        print(
            json.dumps(
                {
                    "status": "ready_for_source_only_value_access",
                    "protocol_sha256": sha256_file(PROTOCOL_PATH),
                    "prevalue_freeze_sha256": sha256_file(FREEZE_PATH),
                    "bound_file_count": len(freeze["bound_files"]),
                    "forbidden": protocol["forbidden"],
                },
                indent=2,
                sort_keys=True,
            )
        )
        return
    if args.command == "verify":
        print(json.dumps(verify_artifact(), indent=2, sort_keys=True))
        return
    if args.confirm != CONFIRMATION:
        raise ConversionError("explicit source-only conversion token is required")
    print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
