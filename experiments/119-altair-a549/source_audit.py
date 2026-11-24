#!/usr/bin/env python3
"""Authenticate and audit the official Altair A549 perturbation-level H5AD.

This module is source-only. It reads the public Altair file, the two public VCC
axis CSVs, and small frozen OPS API receipts. It has no code path for held
truth, model fitting, evaluation, count generation, or submission.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
REPORT_PATH = HERE / "SOURCE_AUDIT.json"
ACQUISITION_RECEIPT_PATH = HERE / "SOURCE_ACQUISITION_RECEIPT.json"

EXPECTED_ROOT_KEYS = {
    "X",
    "layers",
    "obs",
    "obsm",
    "obsp",
    "uns",
    "var",
    "varm",
    "varp",
}
EXPECTED_LAYERS = {
    "adj_p_value",
    "baseMean",
    "lfcSE",
    "neg_log10_fdr",
    "p_value",
    "stat",
}
EXPECTED_OBS_COLUMNS = [
    "role",
    "gene_id",
    "gene_symbol",
    "n_guides",
    "n_cells",
    "control_type",
    "cluster_group_low_res_clusters",
    "cluster_group_high_res_clusters",
]


class SourceAuditError(RuntimeError):
    """A source provenance, schema, or semantic invariant failed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise SourceAuditError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    payload = "".join(f"{value}\n" for value in values).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SourceAuditError(f"cannot read valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise SourceAuditError(f"JSON root must be an object: {path}")
    return value


def load_manifest(path: Path = MANIFEST_PATH) -> dict[str, Any]:
    manifest = load_json(path)
    _require(manifest.get("schema_version") == 1, "manifest schema changed")
    _require(
        manifest.get("experiment") == "119-altair-a549",
        "manifest experiment changed",
    )
    _require(
        manifest.get("scope") == "source_acquisition_and_compact_source_only",
        "manifest scope changed",
    )
    authorization = manifest.get("authorization")
    _require(isinstance(authorization, dict), "authorization is missing")
    for forbidden in (
        "additional_large_downloads",
        "held_truth_access",
        "existing_outcome_access",
        "model_fitting",
        "evaluation",
        "cell_eval2",
        "count_generation",
        "submission_build",
        "submission",
    ):
        _require(
            authorization.get(forbidden) is False,
            f"forbidden authorization enabled: {forbidden}",
        )
    _require(
        authorization.get("source_file_already_acquired") is True,
        "existing source acquisition is not declared",
    )
    return manifest


def safe_repo_file(relative: str, label: str) -> Path:
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()), f"{label} escapes repository")
    _require(path.is_file() and not path.is_symlink(), f"{label} missing or unsafe")
    return path


def safe_metadata_file(relative: str, label: str) -> Path:
    path = (HERE / relative).resolve()
    metadata_root = (HERE / "metadata").resolve()
    _require(path.is_relative_to(metadata_root), f"{label} escapes metadata/")
    _require(path.is_file() and not path.is_symlink(), f"{label} missing or unsafe")
    return path


def authenticate_source(manifest: dict[str, Any]) -> dict[str, Any]:
    source = manifest["source"]
    path = safe_repo_file(source["local_path"], "Altair H5AD")
    size = path.stat().st_size
    _require(size == source["size_bytes"], f"source size drift: {size}")
    digest = sha256_file(path)
    _require(digest == source["sha256"], f"source SHA-256 drift: {digest}")
    return {
        "path": str(path.relative_to(ROOT)),
        "size_bytes": size,
        "sha256": digest,
    }


def authenticate_metadata(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    verified: dict[str, dict[str, Any]] = {}
    for name, asset in manifest["metadata_assets"].items():
        path = safe_metadata_file(asset["path"], name)
        digest = sha256_file(path)
        _require(digest == asset["sha256"], f"{name} metadata SHA-256 drift")
        verified[name] = {
            "path": str(path.relative_to(HERE)),
            "size_bytes": path.stat().st_size,
            "sha256": digest,
        }
    return verified


def read_csv_axis(
    path: Path,
    column: str,
    expected_count: int,
) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        _require(reader.fieldnames == [column], f"{path.name} header changed")
        values = tuple(str(row[column]) for row in reader)
    validate_axis(values, path.name, expected_count)
    return values


def validate_axis(
    values: Sequence[str],
    label: str,
    expected_count: int | None = None,
) -> None:
    if expected_count is not None:
        _require(len(values) == expected_count, f"{label} count changed")
    _require(
        all(value and value == value.strip() for value in values),
        f"{label} contains empty or padded labels",
    )
    _require(len(values) == len(set(values)), f"{label} contains duplicates")


def exact_symbol_projection(
    source_symbols: Sequence[str],
    canonical_symbols: Sequence[str],
) -> np.ndarray:
    """Map each source symbol by exact equality; never resolve aliases."""
    validate_axis(source_symbols, "source response symbols")
    validate_axis(canonical_symbols, "canonical response symbols")
    destination = {symbol: index for index, symbol in enumerate(canonical_symbols)}
    return np.asarray(
        [destination.get(symbol, -1) for symbol in source_symbols],
        dtype=np.int32,
    )


def _attr_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_string_node(node: h5py.Dataset | h5py.Group, label: str) -> np.ndarray:
    if isinstance(node, h5py.Dataset):
        _require(node.ndim == 1, f"{label} must be one-dimensional")
        _require(node.dtype.kind in {"O", "S", "U"}, f"{label} is not text")
        return np.asarray(node.asstr()[...], dtype=object)

    encoding = _attr_text(node.attrs.get("encoding-type", ""))
    if encoding == "categorical":
        _require(set(node.keys()) == {"categories", "codes"}, f"{label} changed")
        categories = read_string_node(node["categories"], f"{label}/categories")
        codes = np.asarray(node["codes"][...])
        _require(codes.ndim == 1, f"{label} codes must be one-dimensional")
        _require(np.all(codes >= 0), f"{label} contains missing values")
        _require(np.all(codes < len(categories)), f"{label} codes are invalid")
        return categories[codes.astype(np.int64, copy=False)]
    if encoding == "nullable-string-array":
        _require(set(node.keys()) == {"mask", "values"}, f"{label} changed")
        mask = np.asarray(node["mask"][...], dtype=bool)
        values = read_string_node(node["values"], f"{label}/values")
        _require(mask.shape == values.shape, f"{label} axes disagree")
        _require(not np.any(mask), f"{label} contains missing values")
        return values
    raise SourceAuditError(f"{label} uses unsupported encoding {encoding!r}")


def read_nullable_integer(node: h5py.Group, label: str) -> np.ndarray:
    _require(
        _attr_text(node.attrs.get("encoding-type", "")) == "nullable-integer",
        f"{label} is not nullable integer",
    )
    _require(set(node.keys()) == {"mask", "values"}, f"{label} members changed")
    mask = np.asarray(node["mask"][...], dtype=bool)
    values = np.asarray(node["values"][...], dtype=np.int64)
    _require(mask.shape == values.shape, f"{label} axes disagree")
    _require(not np.any(mask), f"{label} unexpectedly contains missing values")
    return values


def read_dataframe_index(handle: h5py.File, frame_name: str) -> np.ndarray:
    frame = handle.get(frame_name)
    _require(isinstance(frame, h5py.Group), f"/{frame_name} is missing")
    _require(
        _attr_text(frame.attrs.get("encoding-type", "")) == "dataframe",
        f"/{frame_name} is not an AnnData dataframe",
    )
    index_name = _attr_text(frame.attrs.get("_index", ""))
    _require(index_name in frame, f"/{frame_name} index binding is missing")
    return read_string_node(frame[index_name], f"/{frame_name}/{index_name}")


def require_dense_float32(
    handle: h5py.File,
    path: str,
    shape: tuple[int, int],
) -> h5py.Dataset:
    node = handle.get(path)
    _require(isinstance(node, h5py.Dataset), f"/{path} must be dense")
    _require(node.shape == shape, f"/{path} shape changed: {node.shape}")
    _require(node.dtype == np.dtype(np.float32), f"/{path} dtype changed")
    _require(
        _attr_text(node.attrs.get("encoding-type", "")) == "array",
        f"/{path} AnnData encoding changed",
    )
    return node


def _json_float_vector(values: Sequence[object], label: str) -> np.ndarray:
    result = np.empty(len(values), dtype=np.float64)
    for index, value in enumerate(values):
        if value is None:
            result[index] = np.nan
        else:
            try:
                result[index] = float(value)
            except (TypeError, ValueError) as error:
                raise SourceAuditError(f"{label} contains a non-number") from error
    return result


def compare_api_vector(
    api_values: Sequence[object],
    source_values: np.ndarray,
    tolerance: float,
    label: str,
) -> float:
    api = _json_float_vector(api_values, label)
    source = np.asarray(source_values, dtype=np.float64)
    _require(api.shape == source.shape, f"{label} shape changed")
    _require(
        np.array_equal(np.isnan(api), np.isnan(source)),
        f"{label} missing-value pattern disagrees with H5AD",
    )
    finite = np.isfinite(api)
    maximum = (
        float(np.max(np.abs(api[finite] - source[finite]))) if finite.any() else 0.0
    )
    _require(maximum <= tolerance, f"{label} API mismatch: {maximum}")
    return maximum


def audit_dense_values(
    handle: h5py.File,
    row_names: np.ndarray,
    chunk_rows: int = 32,
) -> dict[str, Any]:
    shape = (len(row_names), len(read_dataframe_index(handle, "var")))
    x = require_dense_float32(handle, "X", shape)
    layers = {
        name: require_dense_float32(handle, f"layers/{name}", shape)
        for name in sorted(EXPECTED_LAYERS)
    }

    x_finite = 0
    x_nan = 0
    fdr_finite = 0
    fdr_nan = 0
    x_min = math.inf
    x_max = -math.inf
    fdr_min = math.inf
    fdr_max = -math.inf
    row_finite = np.zeros(shape[0], dtype=np.int64)
    stat_error = 0.0
    fdr_error = 0.0
    core_nan_consistent = True
    fdr_nan_matches_adjusted = True
    fdr_subset_of_x = True

    for start in range(0, shape[0], chunk_rows):
        stop = min(shape[0], start + chunk_rows)
        xv = np.asarray(x[start:stop], dtype=np.float64)
        base = np.asarray(layers["baseMean"][start:stop], dtype=np.float64)
        se = np.asarray(layers["lfcSE"][start:stop], dtype=np.float64)
        stat = np.asarray(layers["stat"][start:stop], dtype=np.float64)
        adjusted = np.asarray(
            layers["adj_p_value"][start:stop],
            dtype=np.float64,
        )
        neg_fdr = np.asarray(
            layers["neg_log10_fdr"][start:stop],
            dtype=np.float64,
        )

        finite_x = np.isfinite(xv)
        finite_fdr = np.isfinite(neg_fdr)
        row_finite[start:stop] = finite_x.sum(axis=1)
        x_finite += int(finite_x.sum())
        x_nan += int(np.isnan(xv).sum())
        fdr_finite += int(finite_fdr.sum())
        fdr_nan += int(np.isnan(neg_fdr).sum())
        if finite_x.any():
            x_min = min(x_min, float(xv[finite_x].min()))
            x_max = max(x_max, float(xv[finite_x].max()))
        if finite_fdr.any():
            fdr_min = min(fdr_min, float(neg_fdr[finite_fdr].min()))
            fdr_max = max(fdr_max, float(neg_fdr[finite_fdr].max()))

        x_nan_mask = np.isnan(xv)
        core_nan_consistent &= (
            np.array_equal(x_nan_mask, np.isnan(base))
            and np.array_equal(x_nan_mask, np.isnan(se))
            and np.array_equal(x_nan_mask, np.isnan(stat))
        )
        fdr_nan_matches_adjusted &= np.array_equal(
            np.isnan(neg_fdr),
            np.isnan(adjusted),
        )
        fdr_subset_of_x &= not np.any(finite_fdr & ~finite_x)

        stat_mask = finite_x & np.isfinite(se) & np.isfinite(stat) & (se > 0)
        if stat_mask.any():
            stat_error = max(
                stat_error,
                float(np.max(np.abs(xv[stat_mask] / se[stat_mask] - stat[stat_mask]))),
            )
        fdr_relation_mask = (
            np.isfinite(adjusted) & (adjusted > 0) & np.isfinite(neg_fdr)
        )
        if fdr_relation_mask.any():
            fdr_error = max(
                fdr_error,
                float(
                    np.max(
                        np.abs(
                            neg_fdr[fdr_relation_mask]
                            + np.log10(adjusted[fdr_relation_mask])
                        )
                    )
                ),
            )

    return {
        "source_x_finite": x_finite,
        "source_x_nan": x_nan,
        "source_x_min": x_min,
        "source_x_max": x_max,
        "source_fdr_finite": fdr_finite,
        "source_fdr_nan": fdr_nan,
        "source_fdr_min": fdr_min,
        "source_fdr_max": fdr_max,
        "all_nan_rows": row_names[row_finite == 0].tolist(),
        "row_finite_counts": row_finite,
        "core_nan_pattern_consistent": bool(core_nan_consistent),
        "fdr_nan_matches_adjusted_p_value": bool(fdr_nan_matches_adjusted),
        "fdr_support_is_subset_of_x_support": bool(fdr_subset_of_x),
        "max_abs_x_over_lfcse_minus_stat": stat_error,
        "max_abs_neglog10fdr_plus_log10padj": fdr_error,
    }


def require_expected(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    for key, expected_value in expected.items():
        _require(
            actual.get(key) == expected_value,
            f"audit drift for {key}: {actual.get(key)!r} != {expected_value!r}",
        )


def build_report(manifest_path: Path = MANIFEST_PATH) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    source_binding = authenticate_source(manifest)
    metadata_bindings = authenticate_metadata(manifest)

    gene_spec = manifest["canonical_axes"]["genes"]
    panel_spec = manifest["canonical_axes"]["panel"]
    gene_path = safe_repo_file(gene_spec["path"], "canonical genes")
    panel_path = safe_repo_file(panel_spec["path"], "VCC panel")
    _require(
        sha256_file(gene_path) == gene_spec["sha256"],
        "canonical gene axis SHA-256 drift",
    )
    _require(
        sha256_file(panel_path) == panel_spec["sha256"],
        "VCC panel SHA-256 drift",
    )
    canonical = read_csv_axis(gene_path, gene_spec["column"], gene_spec["count"])
    panel = read_csv_axis(panel_path, panel_spec["column"], panel_spec["count"])

    collection = load_json(HERE / manifest["metadata_assets"]["collection"]["path"])
    dataset = load_json(HERE / manifest["metadata_assets"]["dataset"]["path"])
    aggregates = load_json(HERE / manifest["metadata_assets"]["aggregates"]["path"])
    api_alg12 = load_json(
        HERE / manifest["metadata_assets"]["alg12_diff_expression"]["path"]
    )
    object_head = load_json(
        HERE / manifest["metadata_assets"]["public_object_head"]["path"]
    )

    official = manifest["official_sources"]
    _require(
        collection.get("collection_id") == official["collection_id"],
        "collection id changed",
    )
    _require(
        collection.get("perturbation_count") == 1000,
        "collection perturbation count changed",
    )
    _require(
        dataset.get("dataset_id") == official["dataset_id"],
        "dataset id changed",
    )
    _require(dataset.get("modality") == "cropseq", "dataset modality changed")
    _require(
        dataset.get("observation_unit") == ["perturbation_id"],
        "dataset observation unit changed",
    )
    _require(dataset.get("n_cells") == 606075, "dataset cell count changed")

    expected_head = manifest["source"]["public_object_head"]
    normalized_head = {
        "etag": str(object_head["ETag"]).strip('"'),
        "last_modified": object_head["LastModified"],
        "version_id": object_head["VersionId"],
        "content_type": object_head["ContentType"],
        "server_side_encryption": object_head["ServerSideEncryption"],
    }
    _require(
        int(object_head["ContentLength"]) == manifest["source"]["size_bytes"],
        "public object size changed",
    )
    _require(normalized_head == expected_head, "public object HEAD metadata changed")

    source_path = ROOT / source_binding["path"]
    with h5py.File(source_path, "r") as handle:
        _require(set(handle.keys()) == EXPECTED_ROOT_KEYS, "H5AD root schema changed")
        _require(
            _attr_text(handle.attrs.get("encoding-type", "")) == "anndata",
            "root is not AnnData",
        )
        _require(set(handle["layers"].keys()) == EXPECTED_LAYERS, "layers changed")
        _require(
            list(handle["obs"].attrs["column-order"]) == EXPECTED_OBS_COLUMNS,
            "obs column order changed",
        )
        _require(
            list(handle["var"].attrs["column-order"]) == ["gene_symbol"],
            "var columns changed",
        )

        obs_names = read_dataframe_index(handle, "obs")
        var_ids = read_dataframe_index(handle, "var")
        obs_symbols = read_string_node(handle["obs/gene_symbol"], "/obs/gene_symbol")
        var_symbols = read_string_node(handle["var/gene_symbol"], "/var/gene_symbol")
        roles = read_string_node(handle["obs/role"], "/obs/role")
        gene_ids = read_string_node(handle["obs/gene_id"], "/obs/gene_id")
        control_types = read_string_node(
            handle["obs/control_type"]["categories"],
            "/obs/control_type/categories",
        )
        n_cells = np.asarray(handle["obs/n_cells"][...], dtype=np.int64)
        n_guides = np.asarray(handle["obs/n_guides"][...], dtype=np.int64)
        low_clusters = read_nullable_integer(
            handle["obs/cluster_group_low_res_clusters"],
            "/obs/cluster_group_low_res_clusters",
        )
        high_clusters = read_nullable_integer(
            handle["obs/cluster_group_high_res_clusters"],
            "/obs/cluster_group_high_res_clusters",
        )

        validate_axis(obs_names.tolist(), "source perturbation axis", 1001)
        validate_axis(var_ids.tolist(), "source Ensembl axis", 9999)
        validate_axis(var_symbols.tolist(), "source symbol axis", 9999)
        _require(
            np.array_equal(obs_names, obs_symbols),
            "obs index and gene_symbol disagree",
        )
        _require(control_types.tolist() == ["non-targeting"], "control type changed")
        _require(len(gene_ids) == len(obs_names), "obs gene id length changed")
        _require(low_clusters.shape == (1001,), "low cluster axis changed")
        _require(high_clusters.shape == (1001,), "high cluster axis changed")

        target_mask = roles == "targeting"
        control_mask = roles == "control"
        _require(np.sum(target_mask) == 1000, "target row count changed")
        _require(np.sum(control_mask) == 1, "control row count changed")
        _require(
            obs_names[control_mask].tolist() == ["Nontargeting"],
            "control aggregate changed",
        )
        _require(
            np.unique(n_guides[target_mask]).tolist() == [4],
            "target guide count changed",
        )
        _require(
            n_guides[control_mask].tolist() == [211],
            "control guide count changed",
        )

        dense = audit_dense_values(handle, obs_names)
        _require(
            dense["core_nan_pattern_consistent"],
            "X/baseMean/lfcSE/stat missingness disagrees",
        )
        _require(
            dense["fdr_nan_matches_adjusted_p_value"],
            "FDR and adjusted-p missingness disagrees",
        )
        _require(
            dense["fdr_support_is_subset_of_x_support"],
            "FDR is finite where LFC is unsupported",
        )
        _require(
            dense["max_abs_x_over_lfcse_minus_stat"]
            <= manifest["source_semantics"]["deseq_stat_tolerance"],
            "X/lfcSE no longer reproduces stat",
        )
        _require(
            dense["max_abs_neglog10fdr_plus_log10padj"] <= 1e-4,
            "neg_log10_fdr no longer matches adjusted p-values",
        )

        aggregate_ids = aggregates["aggregate_id"]
        _require(aggregate_ids == obs_names.tolist(), "API aggregate axis changed")
        _require(
            aggregates["perturbation_id"] == obs_names.tolist(),
            "API perturbation axis changed",
        )
        _require(
            aggregates["n_cells"] == n_cells.tolist(),
            "API aggregate cell counts changed",
        )

        alg12_index = obs_names.tolist().index("ALG12")
        _require(
            api_alg12["ensembl_id"] == var_ids.tolist(),
            "API differential-expression gene axis changed",
        )
        api_tolerance = manifest["source_semantics"]["api_rounding_tolerance"]
        lfc_api_error = compare_api_vector(
            api_alg12["log2_fold_change"],
            np.asarray(handle["X"][alg12_index], dtype=np.float64),
            api_tolerance,
            "ALG12 log2_fold_change",
        )
        fdr_api_error = compare_api_vector(
            api_alg12["neg_log10_fdr"],
            np.asarray(
                handle["layers/neg_log10_fdr"][alg12_index],
                dtype=np.float64,
            ),
            api_tolerance,
            "ALG12 neg_log10_fdr",
        )

    projection = exact_symbol_projection(var_symbols.tolist(), canonical)
    covered = projection >= 0
    supported_rows = dense["row_finite_counts"] > 0
    source_targets = set(obs_names[target_mask].tolist())
    supported_targets = set(obs_names[target_mask & supported_rows].tolist())
    panel_set = set(panel)
    panel_targets = sorted(source_targets & panel_set)
    panel_supported = sorted(supported_targets & panel_set)

    summary = {
        "source_shape": [len(obs_names), len(var_ids)],
        "target_rows": int(np.sum(target_mask)),
        "control_rows": int(np.sum(control_mask)),
        "supported_target_rows": int(np.sum(target_mask & supported_rows)),
        "all_nan_rows": dense["all_nan_rows"],
        "cells_total": int(n_cells.sum()),
        "target_cells_total": int(n_cells[target_mask].sum()),
        "control_cells_total": int(n_cells[control_mask].sum()),
        "target_cells_min": int(n_cells[target_mask].min()),
        "target_cells_median": float(np.median(n_cells[target_mask])),
        "target_cells_max": int(n_cells[target_mask].max()),
        "source_response_genes": len(var_symbols),
        "canonical_exact_overlap": int(covered.sum()),
        "source_symbols_unmapped": int((~covered).sum()),
        "canonical_symbols_uncovered": len(canonical) - int(covered.sum()),
        "source_x_finite": dense["source_x_finite"],
        "source_x_nan": dense["source_x_nan"],
        "source_fdr_finite": dense["source_fdr_finite"],
        "source_fdr_nan": dense["source_fdr_nan"],
        "vcc_panel_target_overlap": len(panel_targets),
        "vcc_panel_supported_overlap": len(panel_supported),
        "vcc_panel_targets": panel_targets,
        "vcc_panel_supported_targets": panel_supported,
    }
    require_expected(summary, manifest["expected_audit"])

    return {
        "schema_version": 1,
        "experiment": manifest["experiment"],
        "status": "official_source_authenticated_and_semantics_verified",
        "source": source_binding,
        "metadata": metadata_bindings,
        "official_ids": {
            "collection_id": official["collection_id"],
            "dataset_id": official["dataset_id"],
        },
        "schema": {
            "source_shape": summary["source_shape"],
            "matrix_X": {
                "dtype": "float32",
                "meaning": "log2_fold_change",
                "finite": dense["source_x_finite"],
                "nan": dense["source_x_nan"],
                "minimum": dense["source_x_min"],
                "maximum": dense["source_x_max"],
            },
            "layers": sorted(EXPECTED_LAYERS),
            "neg_log10_fdr": {
                "finite": dense["source_fdr_finite"],
                "nan": dense["source_fdr_nan"],
                "minimum": dense["source_fdr_min"],
                "maximum": dense["source_fdr_max"],
            },
        },
        "semantics_evidence": {
            "official_api_field": "log2_fold_change",
            "api_probe_aggregate": "ALG12",
            "api_lfc_max_abs_error_after_five_decimal_serialization": lfc_api_error,
            "api_fdr_max_abs_error_after_five_decimal_serialization": fdr_api_error,
            "max_abs_x_over_lfcse_minus_stat": dense["max_abs_x_over_lfcse_minus_stat"],
            "max_abs_neglog10fdr_plus_log10padj": dense[
                "max_abs_neglog10fdr_plus_log10padj"
            ],
            "core_nan_pattern_consistent": dense["core_nan_pattern_consistent"],
            "fdr_support_is_subset_of_x_support": dense[
                "fdr_support_is_subset_of_x_support"
            ],
            "interpretation": manifest["source_semantics"],
        },
        "summary": summary,
        "axis_bindings": {
            "source_perturbations_sha256": sha256_lines(obs_names.tolist()),
            "source_ensembl_ids_sha256": sha256_lines(var_ids.tolist()),
            "source_gene_symbols_sha256": sha256_lines(var_symbols.tolist()),
            "canonical_genes_sha256": gene_spec["sha256"],
            "vcc_panel_sha256": panel_spec["sha256"],
        },
        "safety": {
            "held_truth_opened": False,
            "existing_outcome_opened": False,
            "model_fitted": False,
            "evaluation_run": False,
            "counts_generated": False,
            "submission_built": False,
        },
    }


def build_acquisition_receipt(
    report: dict[str, Any],
    manifest_path: Path = MANIFEST_PATH,
) -> dict[str, Any]:
    manifest = load_manifest(manifest_path)
    return {
        "schema_version": 1,
        "experiment": manifest["experiment"],
        "status": "verified_existing_official_public_source",
        "verification_date_utc": manifest["date"],
        "source_manifest_sha256": sha256_file(manifest_path),
        "source": report["source"],
        "public_object_head": manifest["source"]["public_object_head"],
        "official_ids": report["official_ids"],
        "semantics_status": report["status"],
        "safety": report["safety"],
    }


def write_json_atomic(path: Path, value: dict[str, Any]) -> None:
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def verify_frozen_reports() -> tuple[dict[str, Any], dict[str, Any]]:
    report = build_report()
    receipt = build_acquisition_receipt(report)
    _require(REPORT_PATH.is_file(), "SOURCE_AUDIT.json is missing")
    _require(
        load_json(REPORT_PATH) == report,
        "SOURCE_AUDIT.json does not match the authenticated source",
    )
    _require(
        ACQUISITION_RECEIPT_PATH.is_file(),
        "SOURCE_ACQUISITION_RECEIPT.json is missing",
    )
    _require(
        load_json(ACQUISITION_RECEIPT_PATH) == receipt,
        "SOURCE_ACQUISITION_RECEIPT.json does not match the authenticated source",
    )
    return report, receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("audit", "verify"))
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()

    if args.command == "verify":
        report, _ = verify_frozen_reports()
    else:
        report = build_report()
        if args.write_report:
            write_json_atomic(REPORT_PATH, report)
            write_json_atomic(
                ACQUISITION_RECEIPT_PATH,
                build_acquisition_receipt(report),
            )

    print(
        json.dumps(
            {
                "status": report["status"],
                "source_shape": report["summary"]["source_shape"],
                "target_rows": report["summary"]["target_rows"],
                "supported_target_rows": report["summary"]["supported_target_rows"],
                "canonical_exact_overlap": report["summary"]["canonical_exact_overlap"],
                "vcc_panel_target_overlap": report["summary"][
                    "vcc_panel_target_overlap"
                ],
                "vcc_panel_supported_overlap": report["summary"][
                    "vcc_panel_supported_overlap"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
