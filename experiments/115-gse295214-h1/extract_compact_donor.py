#!/usr/bin/env python3
"""Build a frozen compact donor from the pinned GSE295214 pseudobulk H5AD.

This program is deliberately source-only. It authenticates the completed
download and its acquisition receipt, reads only the public author-produced
pseudobulk, and writes only beneath this experiment directory. It never opens
held perturbation truth, fits a model, evaluates a prediction, or emits counts.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import numpy as np

import acquire


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
AMENDMENT_PATH = HERE / "PSEUDOBULK_ACQUISITION_AMENDMENT.json"
DOWNLOAD_RECEIPT_PATH = HERE / "work" / "pseudobulk-download-receipt.json"
DEFAULT_OUTPUT_DIR = HERE / "work" / "compact-donor"

EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "cd57ac5692ce3f3e472067e37226b8f2a786a8f8b0cf271ab08d0915f7169050"
)
EXPECTED_AMENDMENT_SHA256 = (
    "4df32fb28300350f3c668c1e6395e1ca69ec7417b2581cfb6ebc86924428019c"
)
EXPECTED_CANONICAL_GENES_SHA256 = (
    "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"
)
EXPECTED_PANEL_SHA256 = (
    "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
)
EXPECTED_TARGET_AXIS_SHA256 = (
    "7060076272d411d4ee00b898cc26eedcb1eaaac5f5e94e2c60f3b45978320747"
)

SOURCE_TARGET_SUM = 10_000.0
PROJECT_TARGET_SUM = 50_000.0
TARGET_SUM_ATOL = 1.0
DEFAULT_CHUNK_ROWS = 32

ARTIFACT_FILENAMES = (
    "b_ctrl.npy",
    "delta.npy",
    "genes.txt",
    "measured_mask.npy",
    "meta.json",
    "perts.txt",
    "rows.tsv",
    "source_covered_gene_mask.npy",
    "source_gene_map.tsv",
    "vcc_target_mask.npy",
    "vcc_targets.txt",
)
ARTIFACT_MANIFEST_NAME = "artifact-manifest.json"
COMPACT_RECEIPT_NAME = "compact-donor-receipt.json"
COMPACT_RECEIPT_HASH_NAME = "compact-donor-receipt.sha256"


class ExtractionError(RuntimeError):
    """A frozen source, schema, conversion, or output invariant failed."""


@dataclass(frozen=True)
class ExtractionContract:
    source_h5ad: Path
    download_receipt: Path
    source_manifest: Path
    amendment: Path
    target_axis: Path
    canonical_genes_csv: Path
    panel_csv: Path
    expected_source_size: int
    expected_source_sha256: str
    expected_source_manifest_sha256: str
    expected_amendment_sha256: str
    expected_target_axis_sha256: str
    expected_canonical_genes_sha256: str
    expected_panel_sha256: str
    expected_rows: int
    expected_vcc_targets: int
    expected_canonical_genes: int
    expected_panel_targets: int
    dataset_label: str = "Pan2025_hESC"


@dataclass(frozen=True)
class AuthenticatedSource:
    path: Path
    size_bytes: int
    sha256: str
    receipt_sha256: str
    source_manifest_sha256: str
    amendment_sha256: str


@dataclass(frozen=True)
class GeneProjection:
    source_ids: tuple[str, ...]
    source_symbols: tuple[str, ...]
    source_to_canonical: np.ndarray
    single_source_indices: np.ndarray
    single_destinations: np.ndarray
    duplicate_groups: tuple[tuple[int, tuple[int, ...]], ...]
    covered_gene_mask: np.ndarray


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def _write_bytes(path: Path, payload: bytes) -> None:
    path.write_bytes(payload)


def _write_lines(path: Path, values: Iterable[str]) -> None:
    path.write_text("".join(f"{value}\n" for value in values), encoding="utf-8")


def _attr_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _decode_scalar(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8")
    return str(value)


def _read_string_vector(node: h5py.Dataset | h5py.Group, label: str) -> np.ndarray:
    """Read AnnData string, nullable-string, or categorical encodings."""
    if isinstance(node, h5py.Dataset):
        if node.ndim != 1:
            raise ExtractionError(f"{label} must be one-dimensional")
        if node.dtype.kind not in {"O", "S", "U"}:
            raise ExtractionError(f"{label} is not a string array")
        values = node.asstr()[...]
        return np.asarray(values, dtype=object)

    encoding = _attr_text(node.attrs.get("encoding-type", ""))
    if encoding == "categorical":
        if set(node.keys()) != {"categories", "codes"}:
            raise ExtractionError(f"{label} categorical members changed")
        categories = _read_string_vector(node["categories"], f"{label}/categories")
        codes = np.asarray(node["codes"][...])
        if codes.ndim != 1 or codes.dtype.kind not in {"i", "u"}:
            raise ExtractionError(f"{label} categorical codes are invalid")
        if np.any(codes < 0) or np.any(codes >= len(categories)):
            raise ExtractionError(f"{label} contains missing or invalid categories")
        return categories[codes.astype(np.int64, copy=False)]

    if encoding == "nullable-string-array":
        if set(node.keys()) != {"mask", "values"}:
            raise ExtractionError(f"{label} nullable-string members changed")
        mask = np.asarray(node["mask"][...], dtype=bool)
        values = _read_string_vector(node["values"], f"{label}/values")
        if mask.shape != values.shape:
            raise ExtractionError(f"{label} nullable-string axes disagree")
        if np.any(mask):
            raise ExtractionError(f"{label} contains missing strings")
        return values

    raise ExtractionError(f"{label} has unsupported encoding {encoding!r}")


def _read_dataframe_column(
    handle: h5py.File,
    frame_name: str,
    column_name: str,
) -> np.ndarray:
    if frame_name not in handle or not isinstance(handle[frame_name], h5py.Group):
        raise ExtractionError(f"/{frame_name} dataframe is missing")
    frame = handle[frame_name]
    if _attr_text(frame.attrs.get("encoding-type", "")) != "dataframe":
        raise ExtractionError(f"/{frame_name} is not an AnnData dataframe")
    if column_name not in frame:
        raise ExtractionError(f"/{frame_name}/{column_name} is missing")
    return _read_string_vector(frame[column_name], f"/{frame_name}/{column_name}")


def _read_dataframe_index(handle: h5py.File, frame_name: str) -> np.ndarray:
    frame = handle.get(frame_name)
    if not isinstance(frame, h5py.Group):
        raise ExtractionError(f"/{frame_name} dataframe is missing")
    index_name = _attr_text(frame.attrs.get("_index", ""))
    if not index_name or index_name not in frame:
        raise ExtractionError(f"/{frame_name} index binding is missing")
    return _read_string_vector(frame[index_name], f"/{frame_name}/{index_name}")


def _require_dense_float32(
    handle: h5py.File,
    path: str,
    expected_shape: tuple[int, int],
) -> h5py.Dataset:
    node = handle.get(path)
    if not isinstance(node, h5py.Dataset):
        raise ExtractionError(f"/{path} must be a dense HDF5 dataset")
    if _attr_text(node.attrs.get("encoding-type", "")) != "array":
        raise ExtractionError(f"/{path} must use AnnData dense-array encoding")
    if _attr_text(node.attrs.get("encoding-version", "")) != "0.2.0":
        raise ExtractionError(f"/{path} dense-array encoding version changed")
    if node.shape != expected_shape:
        raise ExtractionError(
            f"/{path} shape changed: {node.shape!r} != {expected_shape!r}"
        )
    if node.dtype != np.dtype(np.float32):
        raise ExtractionError(f"/{path} dtype changed: {node.dtype} != float32")
    return node


def _read_csv_axis(path: Path, column: str, expected_count: int) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != [column]:
            raise ExtractionError(
                f"{path.name} header changed: expected only {column!r}"
            )
        values = tuple(str(row[column]) for row in reader)
    _validate_axis(values, path.name, expected_count)
    return values


def _read_text_axis(path: Path, expected_count: int) -> tuple[str, ...]:
    if path.read_bytes().endswith(b"\n") is False:
        raise ExtractionError(f"{path.name} must end with a newline")
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    _validate_axis(values, path.name, expected_count)
    return values


def _validate_axis(values: Sequence[str], label: str, expected_count: int) -> None:
    if len(values) != expected_count:
        raise ExtractionError(
            f"{label} row count changed: {len(values)} != {expected_count}"
        )
    if any(not value or value != value.strip() for value in values):
        raise ExtractionError(f"{label} contains empty or padded symbols")
    if len(values) != len(set(values)):
        raise ExtractionError(f"{label} contains duplicate symbols")


def _require_hash(path: Path, expected: str, label: str) -> str:
    if path.is_symlink() or not path.is_file():
        raise ExtractionError(f"{label} is missing or symlinked: {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise ExtractionError(f"{label} SHA-256 changed: {actual}")
    return actual


def production_contract() -> ExtractionContract:
    """Load and cross-check the immutable production extraction contract."""
    _require_hash(
        SOURCE_MANIFEST_PATH,
        EXPECTED_SOURCE_MANIFEST_SHA256,
        "SOURCE_MANIFEST.json",
    )
    _require_hash(
        AMENDMENT_PATH,
        EXPECTED_AMENDMENT_SHA256,
        "PSEUDOBULK_ACQUISITION_AMENDMENT.json",
    )
    manifest = acquire.load_manifest(SOURCE_MANIFEST_PATH)
    amendment = acquire.load_pseudobulk_amendment(manifest, AMENDMENT_PATH)
    asset = manifest["assets"]["public_pseudobulk_database"]
    target_asset = manifest["assets"]["post_qc_target_axis"]
    panel = manifest["local_axes"]["vcc_2026_panel"]
    if amendment["source_semantics"] != {
        "analysis_repository": "xingjiepan/SCMG_scripts",
        "analysis_commit": "aa848b9d517071d27b141e058e9bdae1966a869f",
        "pseudobulk_notebook_path": (
            "manuscript2025/hESC_genome_scale_pert_seq_analysis/"
            "basic_perturb_seq_analysis/generate_pseudo_bulk.ipynb"
        ),
        "pseudobulk_notebook_sha256": (
            "1b57480326b7605bf2f7f1abd540405f7d127fa1a557cc0eec0c3bae7e801f74"
        ),
        "effect_representation": ("log1p(CPM10000_perturbed)-log1p(CPM10000_control)"),
        "required_layers": ["control", "measure_mask"],
        "author_produced_derived_artifact": True,
        "byte_reproduction_from_public_notebooks_claimed": False,
    }:
        raise ExtractionError("pseudobulk source semantics changed")
    return ExtractionContract(
        source_h5ad=(ROOT / asset["local_path"]).resolve(),
        download_receipt=DOWNLOAD_RECEIPT_PATH.resolve(),
        source_manifest=SOURCE_MANIFEST_PATH.resolve(),
        amendment=AMENDMENT_PATH.resolve(),
        target_axis=(ROOT / target_asset["local_path"]).resolve(),
        canonical_genes_csv=(ROOT / "data/challenge/gene_names.csv").resolve(),
        panel_csv=(ROOT / panel["path"]).resolve(),
        expected_source_size=asset["size_bytes"],
        expected_source_sha256=asset["sha256"],
        expected_source_manifest_sha256=EXPECTED_SOURCE_MANIFEST_SHA256,
        expected_amendment_sha256=EXPECTED_AMENDMENT_SHA256,
        expected_target_axis_sha256=EXPECTED_TARGET_AXIS_SHA256,
        expected_canonical_genes_sha256=EXPECTED_CANONICAL_GENES_SHA256,
        expected_panel_sha256=EXPECTED_PANEL_SHA256,
        expected_rows=target_asset["row_count"],
        expected_vcc_targets=manifest["expected_audit"]["vcc_panel_post_qc_overlap"],
        expected_canonical_genes=18_533,
        expected_panel_targets=panel["row_count"],
    )


def authenticate_download(contract: ExtractionContract) -> AuthenticatedSource:
    """Authenticate receipt bindings and then independently hash the H5AD."""
    manifest_hash = _require_hash(
        contract.source_manifest,
        contract.expected_source_manifest_sha256,
        "source manifest",
    )
    amendment_hash = _require_hash(
        contract.amendment,
        contract.expected_amendment_sha256,
        "acquisition amendment",
    )
    if (
        contract.download_receipt.is_symlink()
        or not contract.download_receipt.is_file()
    ):
        raise ExtractionError("completed pseudobulk download receipt is missing")
    receipt_bytes = contract.download_receipt.read_bytes()
    receipt_hash = hashlib.sha256(receipt_bytes).hexdigest()
    try:
        receipt = json.loads(receipt_bytes)
    except json.JSONDecodeError as error:
        raise ExtractionError("download receipt is not valid JSON") from error

    expected_asset_path = str(contract.source_h5ad.relative_to(ROOT))
    expected_receipt = {
        "schema_version": 1,
        "experiment": "115-gse295214-h1",
        "dataset": "GSE295214",
        "status": "verified_public_pseudobulk_download",
        "source_manifest_sha256": manifest_hash,
        "acquisition_amendment_sha256": amendment_hash,
    }
    for key, expected in expected_receipt.items():
        if receipt.get(key) != expected:
            raise ExtractionError(f"download receipt field changed: {key}")
    completed_at = receipt.get("completed_at_utc")
    if not isinstance(completed_at, str) or not completed_at.endswith("Z"):
        raise ExtractionError("download receipt completion timestamp is invalid")

    asset = receipt.get("asset")
    if not isinstance(asset, dict):
        raise ExtractionError("download receipt asset binding is missing")
    expected_asset = {
        "path": expected_asset_path,
        "size_bytes": contract.expected_source_size,
        "sha256": contract.expected_source_sha256,
        "repository_revision": "1d3f28d966772dcf75870beac052923efdc078c5",
        "remote_lfs_sha256": contract.expected_source_sha256,
    }
    if asset != expected_asset:
        raise ExtractionError("download receipt asset binding changed")

    amendment = json.loads(contract.amendment.read_text(encoding="utf-8"))
    if receipt.get("authorization") != amendment.get("authorization"):
        raise ExtractionError("download receipt authorization binding changed")
    if contract.source_h5ad.is_symlink() or not contract.source_h5ad.is_file():
        raise ExtractionError("completed pseudobulk H5AD is missing or symlinked")
    size = contract.source_h5ad.stat().st_size
    if size != contract.expected_source_size:
        raise ExtractionError(
            f"pseudobulk H5AD byte count changed: {size} "
            f"!= {contract.expected_source_size}"
        )
    source_hash = sha256_file(contract.source_h5ad)
    if source_hash != contract.expected_source_sha256:
        raise ExtractionError(f"pseudobulk H5AD SHA-256 changed: {source_hash}")
    return AuthenticatedSource(
        path=contract.source_h5ad,
        size_bytes=size,
        sha256=source_hash,
        receipt_sha256=receipt_hash,
        source_manifest_sha256=manifest_hash,
        amendment_sha256=amendment_hash,
    )


def _build_gene_projection(
    source_ids: Sequence[str],
    source_symbols: Sequence[str],
    canonical_genes: Sequence[str],
) -> GeneProjection:
    if len(source_ids) != len(source_symbols):
        raise ExtractionError("source var index and gene-name axes disagree")
    if len(source_ids) != len(set(source_ids)):
        raise ExtractionError("source var index contains duplicates")
    if any(not value or value != value.strip() for value in source_ids):
        raise ExtractionError("source var index contains empty or padded values")
    if any(not value or value != value.strip() for value in source_symbols):
        raise ExtractionError("source gene names contain empty or padded values")

    canonical_position = {gene: index for index, gene in enumerate(canonical_genes)}
    source_to_canonical = np.asarray(
        [canonical_position.get(symbol, -1) for symbol in source_symbols],
        dtype=np.int32,
    )
    grouped: dict[int, list[int]] = {}
    for source_index, destination in enumerate(source_to_canonical):
        if destination >= 0:
            grouped.setdefault(int(destination), []).append(source_index)
    if not grouped:
        raise ExtractionError("source has no genes on the canonical axis")

    singles: list[tuple[int, int]] = []
    duplicates: list[tuple[int, tuple[int, ...]]] = []
    for destination, source_indices in sorted(grouped.items()):
        if len(source_indices) == 1:
            singles.append((source_indices[0], destination))
        else:
            duplicates.append((destination, tuple(source_indices)))
    covered = np.zeros(len(canonical_genes), dtype=np.bool_)
    covered[np.fromiter(grouped, dtype=np.int64)] = True
    return GeneProjection(
        source_ids=tuple(source_ids),
        source_symbols=tuple(source_symbols),
        source_to_canonical=source_to_canonical,
        single_source_indices=np.asarray([item[0] for item in singles], dtype=np.int64),
        single_destinations=np.asarray([item[1] for item in singles], dtype=np.int64),
        duplicate_groups=tuple(duplicates),
        covered_gene_mask=covered,
    )


def _collapse_linear(
    values: np.ndarray,
    measured_source: np.ndarray,
    projection: GeneProjection,
    canonical_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Collapse source columns; duplicate groups require complete measurement."""
    if values.shape != measured_source.shape or values.ndim != 2:
        raise ExtractionError("linear values and measurement mask disagree")
    rows = values.shape[0]
    collapsed = np.zeros((rows, canonical_count), dtype=np.float64)
    measured = np.zeros((rows, canonical_count), dtype=np.bool_)

    src = projection.single_source_indices
    dst = projection.single_destinations
    if len(src):
        collapsed[:, dst] = values[:, src]
        measured[:, dst] = measured_source[:, src]

    for destination, source_indices_tuple in projection.duplicate_groups:
        source_indices = np.asarray(source_indices_tuple, dtype=np.int64)
        collapsed[:, destination] = np.sum(
            values[:, source_indices],
            axis=1,
            dtype=np.float64,
        )
        measured[:, destination] = np.all(
            measured_source[:, source_indices],
            axis=1,
        )
    return collapsed, measured


def _collapse_control(
    control_cpm10k: np.ndarray,
    projection: GeneProjection,
    canonical_count: int,
) -> np.ndarray:
    values = np.asarray(control_cpm10k, dtype=np.float64)
    if values.ndim != 1 or values.shape[0] != len(projection.source_symbols):
        raise ExtractionError("control vector and source gene axis disagree")
    collapsed = np.zeros(canonical_count, dtype=np.float64)
    src = projection.single_source_indices
    dst = projection.single_destinations
    if len(src):
        collapsed[dst] = values[src]
    for destination, source_indices in projection.duplicate_groups:
        collapsed[destination] = np.sum(
            values[np.asarray(source_indices, dtype=np.int64)],
            dtype=np.float64,
        )
    return collapsed


def _validate_source_rows(
    effect: np.ndarray,
    control_log: np.ndarray,
    mask_values: np.ndarray,
    row_labels: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if effect.shape != control_log.shape or effect.shape != mask_values.shape:
        raise ExtractionError("X, control, and measure_mask chunk shapes disagree")
    if not np.all(np.isfinite(effect)):
        raise ExtractionError(f"X contains nonfinite values near {row_labels[0]}")
    if not np.all(np.isfinite(control_log)):
        raise ExtractionError(f"control contains nonfinite values near {row_labels[0]}")
    if not np.all(np.isfinite(mask_values)):
        raise ExtractionError(
            f"measure_mask contains nonfinite values near {row_labels[0]}"
        )
    if np.any((mask_values != 0.0) & (mask_values != 1.0)):
        raise ExtractionError(
            f"measure_mask is not exactly binary near {row_labels[0]}"
        )
    measured = mask_values.astype(np.bool_, copy=False)
    perturbed_log = control_log + effect
    minimum_control = float(np.min(control_log))
    minimum_perturbed = float(np.min(perturbed_log))
    if minimum_control < -1e-6 or minimum_perturbed < -1e-6:
        raise ExtractionError(
            "log1p source values imply negative CPM near "
            f"{row_labels[0]}: control={minimum_control}, "
            f"perturbed={minimum_perturbed}"
        )
    control_cpm = np.expm1(np.maximum(control_log, 0.0), dtype=np.float64)
    perturbed_cpm = np.expm1(np.maximum(perturbed_log, 0.0), dtype=np.float64)
    control_sums = np.sum(control_cpm, axis=1, dtype=np.float64)
    perturbed_sums = np.sum(perturbed_cpm, axis=1, dtype=np.float64)
    bad_control = np.flatnonzero(
        np.abs(control_sums - SOURCE_TARGET_SUM) > TARGET_SUM_ATOL
    )
    bad_perturbed = np.flatnonzero(
        np.abs(perturbed_sums - SOURCE_TARGET_SUM) > TARGET_SUM_ATOL
    )
    if len(bad_control):
        index = int(bad_control[0])
        raise ExtractionError(
            f"control CPM10k sum drift at {row_labels[index]}: "
            f"{control_sums[index]:.9f}"
        )
    if len(bad_perturbed):
        index = int(bad_perturbed[0])
        raise ExtractionError(
            f"perturbed CPM10k sum drift at {row_labels[index]}: "
            f"{perturbed_sums[index]:.9f}"
        )
    return perturbed_cpm, control_cpm, measured


def _safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved == HERE or not resolved.is_relative_to(HERE):
        raise ExtractionError(f"output must be a child of {HERE}: {resolved}")
    return resolved


def _artifact_binding(path: Path) -> dict[str, Any]:
    binding: dict[str, Any] = {
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }
    if path.suffix == ".npy":
        array = np.load(path, mmap_mode="r", allow_pickle=False)
        binding["shape"] = list(array.shape)
        binding["dtype"] = array.dtype.str
    return binding


def _write_source_gene_map(
    path: Path,
    projection: GeneProjection,
    canonical_genes: Sequence[str],
) -> None:
    multiplicity: dict[int, int] = {}
    for destination in projection.source_to_canonical:
        if destination >= 0:
            multiplicity[int(destination)] = multiplicity.get(int(destination), 0) + 1
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "source_index",
                "source_id",
                "gene_symbol",
                "canonical_index",
                "canonical_gene",
                "mapped",
                "symbol_multiplicity",
            ]
        )
        for source_index, (source_id, symbol, destination) in enumerate(
            zip(
                projection.source_ids,
                projection.source_symbols,
                projection.source_to_canonical,
                strict=True,
            )
        ):
            mapped = int(destination >= 0)
            writer.writerow(
                [
                    source_index,
                    source_id,
                    symbol,
                    int(destination),
                    canonical_genes[int(destination)] if mapped else "",
                    mapped,
                    multiplicity.get(int(destination), 0) if mapped else 0,
                ]
            )


def extract_compact_donor(
    authenticated: AuthenticatedSource,
    contract: ExtractionContract,
    output_dir: Path,
    *,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
) -> dict[str, Any]:
    """Extract all frozen H1 rows into canonical TS=50,000 donor artifacts."""
    output_dir = _safe_output(output_dir)
    if chunk_rows <= 0:
        raise ExtractionError("chunk_rows must be positive")
    if output_dir.exists():
        raise ExtractionError(f"output already exists: {output_dir}")
    implementation_sha256 = sha256_file(Path(__file__).resolve())

    _require_hash(
        contract.target_axis,
        contract.expected_target_axis_sha256,
        "frozen Pan2025_hESC target axis",
    )
    _require_hash(
        contract.canonical_genes_csv,
        contract.expected_canonical_genes_sha256,
        "canonical gene axis",
    )
    _require_hash(
        contract.panel_csv,
        contract.expected_panel_sha256,
        "VCC target panel",
    )
    frozen_targets = _read_text_axis(contract.target_axis, contract.expected_rows)
    canonical_genes = _read_csv_axis(
        contract.canonical_genes_csv,
        "gene_name",
        contract.expected_canonical_genes,
    )
    panel = _read_csv_axis(
        contract.panel_csv,
        "target_gene",
        contract.expected_panel_targets,
    )
    panel_set = set(panel)
    vcc_targets = tuple(target for target in frozen_targets if target in panel_set)
    if len(vcc_targets) != contract.expected_vcc_targets:
        raise ExtractionError(
            f"VCC target overlap changed: {len(vcc_targets)} "
            f"!= {contract.expected_vcc_targets}"
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(dir=output_dir.parent, prefix=f".{output_dir.name}.")
    )
    _safe_output(temporary)
    try:
        with h5py.File(authenticated.path, "r") as handle:
            if _attr_text(handle.attrs.get("encoding-type", "")) != "anndata":
                raise ExtractionError("source root is not an AnnData H5AD")
            datasets = _read_dataframe_column(handle, "obs", "dataset")
            target_names = _read_dataframe_column(
                handle,
                "obs",
                "perturbed_gene_name",
            )
            obs_ids = _read_dataframe_index(handle, "obs")
            source_ids = tuple(
                _decode_scalar(value) for value in _read_dataframe_index(handle, "var")
            )
            source_symbols = tuple(
                _decode_scalar(value)
                for value in _read_dataframe_column(handle, "var", "gene_name")
            )

            n_obs = len(datasets)
            n_vars = len(source_symbols)
            if len(target_names) != n_obs or len(obs_ids) != n_obs:
                raise ExtractionError("source obs columns disagree")
            if len(source_ids) != n_vars:
                raise ExtractionError("source var columns disagree")
            shape = (n_obs, n_vars)
            x = _require_dense_float32(handle, "X", shape)
            control = _require_dense_float32(handle, "layers/control", shape)
            measure_mask = _require_dense_float32(
                handle,
                "layers/measure_mask",
                shape,
            )

            dataset_values = np.asarray(
                [_decode_scalar(value) for value in datasets],
                dtype=object,
            )
            target_values = np.asarray(
                [_decode_scalar(value) for value in target_names],
                dtype=object,
            )
            obs_id_values = np.asarray(
                [_decode_scalar(value) for value in obs_ids],
                dtype=object,
            )
            selected = np.flatnonzero(dataset_values == contract.dataset_label)
            if len(selected) != contract.expected_rows:
                raise ExtractionError(
                    f"{contract.dataset_label} row count changed: "
                    f"{len(selected)} != {contract.expected_rows}"
                )
            selected_targets = [str(target_values[index]) for index in selected]
            _validate_axis(
                selected_targets,
                f"{contract.dataset_label} source target axis",
                contract.expected_rows,
            )
            if tuple(sorted(selected_targets)) != frozen_targets:
                raise ExtractionError(
                    "Pan2025_hESC source targets differ from the frozen target axis"
                )
            target_to_source = dict(zip(selected_targets, selected, strict=True))
            ordered_source_rows = np.asarray(
                [target_to_source[target] for target in frozen_targets],
                dtype=np.int64,
            )
            source_order = np.argsort(ordered_source_rows, kind="stable")
            sorted_source_rows = ordered_source_rows[source_order]
            sorted_destinations = source_order.astype(np.int64, copy=False)

            projection = _build_gene_projection(
                source_ids,
                source_symbols,
                canonical_genes,
            )
            np.save(
                temporary / "source_covered_gene_mask.npy",
                projection.covered_gene_mask,
                allow_pickle=False,
            )
            _write_source_gene_map(
                temporary / "source_gene_map.tsv",
                projection,
                canonical_genes,
            )
            _write_lines(temporary / "genes.txt", canonical_genes)
            _write_lines(temporary / "perts.txt", frozen_targets)
            _write_lines(temporary / "vcc_targets.txt", vcc_targets)
            vcc_target_mask = np.asarray(
                [target in panel_set for target in frozen_targets],
                dtype=np.bool_,
            )
            np.save(
                temporary / "vcc_target_mask.npy",
                vcc_target_mask,
                allow_pickle=False,
            )

            with (temporary / "rows.tsv").open(
                "w",
                newline="",
                encoding="utf-8",
            ) as handle_rows:
                writer = csv.writer(
                    handle_rows,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writerow(
                    [
                        "output_index",
                        "source_row_index",
                        "source_obs_id",
                        "target_gene",
                        "is_vcc_2026_target",
                    ]
                )
                for output_index, (source_row, target, is_vcc) in enumerate(
                    zip(
                        ordered_source_rows,
                        frozen_targets,
                        vcc_target_mask,
                        strict=True,
                    )
                ):
                    writer.writerow(
                        [
                            output_index,
                            int(source_row),
                            str(obs_id_values[source_row]),
                            target,
                            int(is_vcc),
                        ]
                    )

            delta = np.lib.format.open_memmap(
                temporary / "delta.npy",
                mode="w+",
                dtype=np.float32,
                shape=(contract.expected_rows, contract.expected_canonical_genes),
            )
            measured_output = np.lib.format.open_memmap(
                temporary / "measured_mask.npy",
                mode="w+",
                dtype=np.bool_,
                shape=(contract.expected_rows, contract.expected_canonical_genes),
            )
            delta[:] = 0.0
            measured_output[:] = False

            first_source_row = int(sorted_source_rows[0])
            reference_control_log = np.asarray(
                control[first_source_row, :],
                dtype=np.float64,
            )
            if not np.all(np.isfinite(reference_control_log)):
                raise ExtractionError("reference control row contains nonfinite values")
            if float(np.min(reference_control_log)) < -1e-6:
                raise ExtractionError("reference control row implies negative CPM")
            reference_control_cpm = np.expm1(
                np.maximum(reference_control_log, 0.0),
                dtype=np.float64,
            )
            reference_sum = float(np.sum(reference_control_cpm, dtype=np.float64))
            if abs(reference_sum - SOURCE_TARGET_SUM) > TARGET_SUM_ATOL:
                raise ExtractionError(
                    f"reference control CPM10k sum drift: {reference_sum:.9f}"
                )
            collapsed_control = _collapse_control(
                reference_control_cpm,
                projection,
                contract.expected_canonical_genes,
            )
            b_ctrl = np.log1p(
                (PROJECT_TARGET_SUM / SOURCE_TARGET_SUM) * collapsed_control
            ).astype(np.float32)
            np.save(temporary / "b_ctrl.npy", b_ctrl, allow_pickle=False)

            measured_counts = np.empty(contract.expected_rows, dtype=np.int32)
            for start in range(0, contract.expected_rows, chunk_rows):
                end = min(start + chunk_rows, contract.expected_rows)
                source_rows = sorted_source_rows[start:end]
                destinations = sorted_destinations[start:end]
                if len(source_rows) > 1 and np.any(np.diff(source_rows) <= 0):
                    raise ExtractionError("selected source rows are not unique")
                effect_chunk = np.asarray(x[source_rows, :], dtype=np.float64)
                control_chunk = np.asarray(
                    control[source_rows, :],
                    dtype=np.float64,
                )
                mask_chunk = np.asarray(
                    measure_mask[source_rows, :],
                    dtype=np.float64,
                )
                row_labels = [
                    f"{contract.dataset_label}/{target_values[index]}"
                    for index in source_rows
                ]
                perturbed_cpm, control_cpm, measured_source = _validate_source_rows(
                    effect_chunk,
                    control_chunk,
                    mask_chunk,
                    row_labels,
                )
                if not np.array_equal(
                    control_chunk,
                    np.broadcast_to(reference_control_log, control_chunk.shape),
                ):
                    raise ExtractionError(
                        f"{contract.dataset_label} control rows are not identical"
                    )
                collapsed_perturbed, measured = _collapse_linear(
                    perturbed_cpm,
                    measured_source,
                    projection,
                    contract.expected_canonical_genes,
                )
                collapsed_chunk_control, measured_again = _collapse_linear(
                    control_cpm,
                    measured_source,
                    projection,
                    contract.expected_canonical_genes,
                )
                if not np.array_equal(measured, measured_again):
                    raise ExtractionError("collapsed measurement masks disagree")
                converted = np.log1p(
                    (PROJECT_TARGET_SUM / SOURCE_TARGET_SUM) * collapsed_perturbed
                ) - np.log1p(
                    (PROJECT_TARGET_SUM / SOURCE_TARGET_SUM) * collapsed_chunk_control
                )
                converted[~measured] = 0.0
                if not np.all(np.isfinite(converted)):
                    raise ExtractionError(
                        f"converted effects are nonfinite near {row_labels[0]}"
                    )
                delta[destinations] = converted.astype(np.float32)
                measured_output[destinations] = measured
                measured_counts[destinations] = np.sum(
                    measured,
                    axis=1,
                    dtype=np.int64,
                ).astype(np.int32)

            delta.flush()
            measured_output.flush()
            del delta
            del measured_output

        duplicate_rows = len(source_symbols) - len(set(source_symbols))
        mapped_rows = int(np.count_nonzero(projection.source_to_canonical >= 0))
        metadata = {
            "format_version": 1,
            "created_by": ("experiments/115-gse295214-h1/extract_compact_donor.py"),
            "implementation_sha256": implementation_sha256,
            "source": {
                "dataset": "GSE295214",
                "dataset_label": contract.dataset_label,
                "h5ad_size_bytes": authenticated.size_bytes,
                "h5ad_sha256": authenticated.sha256,
                "download_receipt_sha256": authenticated.receipt_sha256,
                "source_manifest_sha256": (authenticated.source_manifest_sha256),
                "acquisition_amendment_sha256": authenticated.amendment_sha256,
                "effect_representation": (
                    "log1p(CPM10000_perturbed)-log1p(CPM10000_control)"
                ),
            },
            "conversion": {
                "source_target_sum": SOURCE_TARGET_SUM,
                "project_target_sum": PROJECT_TARGET_SUM,
                "method": (
                    "reconstruct linear perturbed/control CPM10k with expm1; "
                    "sum duplicate source symbols in linear space; multiply "
                    "linear CPM by 5; compute log1p(TS50000 perturbed) minus "
                    "log1p(TS50000 control)"
                ),
                "unmapped_source_symbols_are_not_renormalized": True,
                "unsupported_delta_is_exact_zero": True,
                "duplicate_group_support_requires_all_source_coordinates": True,
                "control_rows_identical": True,
            },
            "axes": {
                "rows": contract.expected_rows,
                "canonical_genes": contract.expected_canonical_genes,
                "vcc_targets": contract.expected_vcc_targets,
                "source_features": len(source_symbols),
                "source_unique_gene_symbols": len(set(source_symbols)),
                "source_duplicate_feature_rows": duplicate_rows,
                "mapped_source_feature_rows": mapped_rows,
                "covered_canonical_genes": int(
                    np.count_nonzero(projection.covered_gene_mask)
                ),
                "measured_coordinates_per_row": {
                    "min": int(np.min(measured_counts)),
                    "median": float(np.median(measured_counts)),
                    "max": int(np.max(measured_counts)),
                },
            },
            "array_contract": {
                "delta.npy": {
                    "shape": [
                        contract.expected_rows,
                        contract.expected_canonical_genes,
                    ],
                    "dtype": "<f4",
                    "semantics": "TS50000 log1p perturbation effect",
                },
                "b_ctrl.npy": {
                    "shape": [contract.expected_canonical_genes],
                    "dtype": "<f4",
                    "semantics": "TS50000 log1p shared H1 control",
                },
                "measured_mask.npy": {
                    "shape": [
                        contract.expected_rows,
                        contract.expected_canonical_genes,
                    ],
                    "dtype": "|b1",
                    "semantics": (
                        "row-coordinate support after exact-symbol mapping; "
                        "all duplicate source coordinates must be measured"
                    ),
                },
                "source_covered_gene_mask.npy": {
                    "shape": [contract.expected_canonical_genes],
                    "dtype": "|b1",
                    "semantics": "canonical gene has at least one source column",
                },
                "vcc_target_mask.npy": {
                    "shape": [contract.expected_rows],
                    "dtype": "|b1",
                    "semantics": "row target belongs to the VCC-2026 panel",
                },
            },
            "safety": {
                "held_truth_opened": False,
                "model_fit": False,
                "evaluation_run": False,
                "counts_generated": False,
                "submission_built": False,
            },
        }
        _write_bytes(temporary / "meta.json", _json_bytes(metadata))

        artifact_manifest = {
            "schema_version": 1,
            "experiment": "115-gse295214-h1",
            "artifact_set": "Pan2025_hESC_compact_source_donor",
            "files": {
                name: _artifact_binding(temporary / name) for name in ARTIFACT_FILENAMES
            },
        }
        _write_bytes(
            temporary / ARTIFACT_MANIFEST_NAME,
            _json_bytes(artifact_manifest),
        )
        receipt = {
            "schema_version": 1,
            "experiment": "115-gse295214-h1",
            "status": "frozen_compact_source_donor",
            "implementation": {
                "path": "experiments/115-gse295214-h1/extract_compact_donor.py",
                "sha256": implementation_sha256,
            },
            "source": {
                "h5ad_size_bytes": authenticated.size_bytes,
                "h5ad_sha256": authenticated.sha256,
                "download_receipt_sha256": authenticated.receipt_sha256,
                "source_manifest_sha256": (authenticated.source_manifest_sha256),
                "acquisition_amendment_sha256": authenticated.amendment_sha256,
            },
            "axes": {
                "rows": contract.expected_rows,
                "canonical_genes": contract.expected_canonical_genes,
                "vcc_targets": contract.expected_vcc_targets,
            },
            "artifact_manifest": {
                "path": ARTIFACT_MANIFEST_NAME,
                "size_bytes": (temporary / ARTIFACT_MANIFEST_NAME).stat().st_size,
                "sha256": sha256_file(temporary / ARTIFACT_MANIFEST_NAME),
            },
        }
        receipt_bytes = _json_bytes(receipt)
        _write_bytes(temporary / COMPACT_RECEIPT_NAME, receipt_bytes)
        receipt_hash = hashlib.sha256(receipt_bytes).hexdigest()
        _write_bytes(
            temporary / COMPACT_RECEIPT_HASH_NAME,
            f"{receipt_hash}  {COMPACT_RECEIPT_NAME}\n".encode(),
        )
        os.replace(temporary, output_dir)
        return receipt
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def verify_compact_donor(
    output_dir: Path,
    contract: ExtractionContract | None = None,
) -> dict[str, Any]:
    """Verify the compact donor's complete deterministic hash chain."""
    output_dir = _safe_output(output_dir)
    if output_dir.is_symlink() or not output_dir.is_dir():
        raise ExtractionError(f"compact donor is missing or symlinked: {output_dir}")
    hash_path = output_dir / COMPACT_RECEIPT_HASH_NAME
    receipt_path = output_dir / COMPACT_RECEIPT_NAME
    manifest_path = output_dir / ARTIFACT_MANIFEST_NAME
    if any(
        path.is_symlink() or not path.is_file()
        for path in (hash_path, receipt_path, manifest_path)
    ):
        raise ExtractionError("compact donor receipt chain is incomplete")

    fields = hash_path.read_text(encoding="utf-8").rstrip("\n").split("  ")
    if len(fields) != 2 or fields[1] != COMPACT_RECEIPT_NAME:
        raise ExtractionError("compact donor receipt checksum file is malformed")
    if sha256_file(receipt_path) != fields[0]:
        raise ExtractionError("compact donor receipt SHA-256 mismatch")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") != "frozen_compact_source_donor":
        raise ExtractionError("compact donor receipt status changed")
    if contract is not None:
        expected_implementation = {
            "path": "experiments/115-gse295214-h1/extract_compact_donor.py",
            "sha256": sha256_file(Path(__file__).resolve()),
        }
        if receipt.get("implementation") != expected_implementation:
            raise ExtractionError("compact donor implementation binding changed")
        expected_source = {
            "h5ad_size_bytes": contract.expected_source_size,
            "h5ad_sha256": contract.expected_source_sha256,
            "download_receipt_sha256": sha256_file(contract.download_receipt),
            "source_manifest_sha256": contract.expected_source_manifest_sha256,
            "acquisition_amendment_sha256": contract.expected_amendment_sha256,
        }
        if receipt.get("source") != expected_source:
            raise ExtractionError("compact donor source binding changed")
        if receipt.get("axes") != {
            "rows": contract.expected_rows,
            "canonical_genes": contract.expected_canonical_genes,
            "vcc_targets": contract.expected_vcc_targets,
        }:
            raise ExtractionError("compact donor axis binding changed")
    manifest_binding = receipt.get("artifact_manifest")
    if not isinstance(manifest_binding, dict):
        raise ExtractionError("compact donor artifact-manifest binding is missing")
    if manifest_binding != {
        "path": ARTIFACT_MANIFEST_NAME,
        "size_bytes": manifest_path.stat().st_size,
        "sha256": sha256_file(manifest_path),
    }:
        raise ExtractionError("compact donor artifact-manifest binding changed")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ExtractionError("unsupported compact artifact manifest schema")
    files = manifest.get("files")
    if not isinstance(files, dict) or set(files) != set(ARTIFACT_FILENAMES):
        raise ExtractionError("compact donor artifact inventory changed")
    for name in ARTIFACT_FILENAMES:
        path = output_dir / name
        if path.is_symlink() or not path.is_file():
            raise ExtractionError(f"compact artifact is missing or symlinked: {name}")
        actual = _artifact_binding(path)
        if actual != files[name]:
            raise ExtractionError(f"compact artifact binding changed: {name}")
    return receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="must be a new directory beneath this experiment",
    )
    parser.add_argument("--chunk-rows", type=int, default=DEFAULT_CHUNK_ROWS)
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="verify an existing compact donor without opening the source H5AD",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.verify_only:
        contract = production_contract()
        print(
            json.dumps(
                verify_compact_donor(args.output_dir, contract),
                indent=2,
                sort_keys=True,
            )
        )
        return
    contract = production_contract()
    authenticated = authenticate_download(contract)
    receipt = extract_compact_donor(
        authenticated,
        contract,
        args.output_dir,
        chunk_rows=args.chunk_rows,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
