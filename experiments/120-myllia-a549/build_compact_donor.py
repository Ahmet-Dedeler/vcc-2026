#!/usr/bin/env python3
"""Audit GSE337804 A549 inputs and build a frozen source-only compact donor.

This program never reads VCC perturbation truth or existing experiment outcome
directories. It authenticates only the files in SOURCE_MANIFEST.json, computes
guide effects independently in two source channels, and writes only below this
experiment directory.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import numpy as np
from scipy import sparse


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
PROVENANCE_PATH = HERE / "SOURCE_PROVENANCE.json"
DEFAULT_AUDIT_PATH = HERE / "work" / "source-audit.json"
DEFAULT_OUTPUT_DIR = HERE / "work" / "compact-donor"

EXPECTED_MANIFEST_SHA256 = (
    "eb6ad356ea88fee18335919e755b4720defb29dbf1b159bf5ad30927b5172bf5"
)
EXPECTED_PROVENANCE_SHA256 = (
    "5598343ca28e4d0e0d99cdcb7af90e6d6d189b704734144a922f33be4a7e45fd"
)

CHANNEL_ORDER = ("ch_1", "ch_2")
TARGET_SUM = 50_000.0
DEFAULT_CHUNK_ROWS = 512
CONTROL_SYMBOL = "non-targeting"
DESERT_SYMBOL = "gene_desert"

ARRAY_FILES = (
    "b_ctrl.npy",
    "channel_b_ctrl.npy",
    "channel_source_covered_gene_mask.npy",
    "delta.npy",
    "guide_channel_delta.npy",
    "guide_channel_measured_mask.npy",
    "guide_delta.npy",
    "guide_measured_mask.npy",
    "measured_mask.npy",
    "panel_coordinate_mask.npy",
    "source_covered_gene_mask.npy",
    "target_qc_mask.npy",
    "transferable_gene_mask.npy",
)
TEXT_FILES = (
    "genes.txt",
    "guide_channel_rows.tsv",
    "guide_qc.tsv",
    "guide_rows.tsv",
    "meta.json",
    "perts.txt",
    "source_targets.txt",
    "target_qc.tsv",
    "vcc_design_overlap.txt",
    "vcc_targets.txt",
)
ARTIFACT_MANIFEST_NAME = "artifact-manifest.json"
RECEIPT_NAME = "compact-donor-receipt.json"
RECEIPT_HASH_NAME = "compact-donor-receipt.sha256"


class SourceError(RuntimeError):
    """A source, schema, safety, QC, or output invariant failed."""


@dataclass(frozen=True)
class Contract:
    manifest: dict[str, Any]
    provenance_sha256: str
    canonical_genes: tuple[str, ...]
    panel: tuple[str, ...]
    channel_paths: dict[str, Path]
    library_path: Path
    pdf_path: Path
    min_cells: int
    min_paired_guides: int


@dataclass(frozen=True)
class ChannelMetadata:
    channel: str
    path: Path
    obs_ids: np.ndarray
    guide_ids: np.ndarray
    symbols: np.ndarray
    genes: tuple[str, ...]
    ncount: np.ndarray
    nfeature: np.ndarray
    cells: int
    nnz: int


@dataclass(frozen=True)
class ChannelAggregate:
    channel: str
    guide_ids: tuple[str, ...]
    guide_symbols: tuple[str, ...]
    guide_cells: np.ndarray
    delta_full: np.ndarray
    control_log: np.ndarray
    control_projected_counts: np.ndarray
    control_total_counts: int
    covered_mask: np.ndarray


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def axis_sha256(values: Sequence[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode()
    ).hexdigest()


def json_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def write_lines(path: Path, values: Iterable[str]) -> None:
    path.write_text(
        "".join(f"{value}\n" for value in values),
        encoding="utf-8",
    )


def attr_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def decode_scalar(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8")
    return str(value)


def read_string_vector(node: h5py.Dataset | h5py.Group, label: str) -> np.ndarray:
    if isinstance(node, h5py.Dataset):
        if node.ndim != 1 or node.dtype.kind not in {"O", "S", "U"}:
            raise SourceError(f"{label} is not a one-dimensional string array")
        return np.asarray(node.asstr()[...], dtype=object)
    if attr_text(node.attrs.get("encoding-type", "")) != "categorical":
        raise SourceError(f"{label} has unsupported string encoding")
    if set(node.keys()) != {"categories", "codes"}:
        raise SourceError(f"{label} categorical members changed")
    categories = read_string_vector(node["categories"], f"{label}/categories")
    codes = np.asarray(node["codes"][...])
    if codes.ndim != 1 or codes.dtype.kind not in {"i", "u"}:
        raise SourceError(f"{label} categorical codes are invalid")
    if np.any(codes < 0) or np.any(codes >= len(categories)):
        raise SourceError(f"{label} categorical codes are out of range")
    return categories[codes.astype(np.int64, copy=False)]


def read_frame_index(frame: h5py.Group, label: str) -> np.ndarray:
    index_name = attr_text(frame.attrs.get("_index", ""))
    if not index_name or index_name not in frame:
        raise SourceError(f"{label} index binding is missing")
    return read_string_vector(frame[index_name], f"{label}/{index_name}")


def read_csv_axis(path: Path, column: str, expected: int) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != [column]:
            raise SourceError(f"{path.name} header changed")
        values = tuple(str(row[column]) for row in reader)
    validate_axis(values, path.name, expected)
    return values


def validate_axis(values: Sequence[str], label: str, expected: int) -> None:
    if len(values) != expected:
        raise SourceError(f"{label} length changed: {len(values)} != {expected}")
    if any(not value or value != value.strip() for value in values):
        raise SourceError(f"{label} contains empty or padded values")
    if len(values) != len(set(values)):
        raise SourceError(f"{label} contains duplicate values")


def require_file(path: Path, size: int, digest: str, label: str) -> str:
    if path.is_symlink() or not path.is_file():
        raise SourceError(f"{label} is missing or symlinked: {path}")
    actual_size = path.stat().st_size
    if actual_size != size:
        raise SourceError(f"{label} size changed: {actual_size} != {size}")
    actual_digest = sha256_file(path)
    if actual_digest != digest:
        raise SourceError(f"{label} SHA-256 changed: {actual_digest}")
    return actual_digest


def safe_child(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if resolved == HERE or not resolved.is_relative_to(HERE):
        raise SourceError(f"output must be a child of {HERE}: {resolved}")
    return resolved


def load_contract() -> Contract:
    if sha256_file(MANIFEST_PATH) != EXPECTED_MANIFEST_SHA256:
        raise SourceError("SOURCE_MANIFEST.json SHA-256 changed")
    provenance_hash = sha256_file(PROVENANCE_PATH)
    if provenance_hash != EXPECTED_PROVENANCE_SHA256:
        raise SourceError("SOURCE_PROVENANCE.json SHA-256 changed")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    authorization = manifest.get("authorization")
    expected_authorization = {
        "build_compact_source_donor": True,
        "cell_eval2": False,
        "challenge_build": False,
        "commit": False,
        "count_generation": False,
        "evaluation": False,
        "open_existing_outcome_directories": False,
        "open_vcc_held_truth": False,
        "push": False,
        "submission": False,
    }
    if authorization != expected_authorization:
        raise SourceError("source-only authorization changed")
    qc = manifest["qc_contract"]
    if qc != {
        "channel_pooling": "equal_weight_mean",
        "control_symbol": CONTROL_SYMBOL,
        "exact_symbols_only": True,
        "guide_channel_min_cells": 10,
        "guide_requires_both_channels": True,
        "minimum_paired_guides_per_target": 2,
        "panel_coordinates_transferable": False,
        "pseudobulk_target_sum": 50000,
        "target_pooling": "coordinate_wise_median_of_channel_pooled_guides",
    }:
        raise SourceError("QC contract changed")

    inputs = manifest["inputs"]
    channel_paths = {
        "ch_1": (ROOT / inputs["channel_1_h5ad"]["local_path"]).resolve(),
        "ch_2": (ROOT / inputs["channel_2_h5ad"]["local_path"]).resolve(),
    }
    axes = manifest["local_axes"]
    canonical_path = (ROOT / axes["canonical_genes"]["local_path"]).resolve()
    panel_path = (ROOT / axes["vcc_panel"]["local_path"]).resolve()
    require_file(
        canonical_path,
        canonical_path.stat().st_size,
        axes["canonical_genes"]["sha256"],
        "canonical gene axis",
    )
    require_file(
        panel_path,
        panel_path.stat().st_size,
        axes["vcc_panel"]["sha256"],
        "VCC panel axis",
    )
    canonical = read_csv_axis(
        canonical_path,
        "gene_name",
        axes["canonical_genes"]["row_count"],
    )
    panel = read_csv_axis(
        panel_path,
        "target_gene",
        axes["vcc_panel"]["row_count"],
    )
    return Contract(
        manifest=manifest,
        provenance_sha256=provenance_hash,
        canonical_genes=canonical,
        panel=panel,
        channel_paths=channel_paths,
        library_path=(ROOT / inputs["sgrna_library"]["local_path"]).resolve(),
        pdf_path=(ROOT / inputs["description_pdf"]["local_path"]).resolve(),
        min_cells=int(qc["guide_channel_min_cells"]),
        min_paired_guides=int(qc["minimum_paired_guides_per_target"]),
    )


def authenticate_inputs(contract: Contract) -> dict[str, dict[str, Any]]:
    inputs = contract.manifest["inputs"]
    bindings: dict[str, dict[str, Any]] = {}
    for key, item in inputs.items():
        path = (ROOT / item["local_path"]).resolve()
        digest = require_file(path, item["size_bytes"], item["sha256"], key)
        bindings[key] = {
            "path": item["local_path"],
            "sha256": digest,
            "size_bytes": path.stat().st_size,
        }
    return bindings


def load_library(path: Path) -> tuple[list[dict[str, str]], dict[str, dict[str, str]]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected_fields = ["id", "symbol", "rank", "class", "entrez", "guide", "count"]
        if reader.fieldnames != expected_fields:
            raise SourceError(f"sgRNA library columns changed: {reader.fieldnames}")
        rows = [dict(row) for row in reader]
    if len(rows) != 694:
        raise SourceError(f"sgRNA library row count changed: {len(rows)}")
    ids = [row["id"] for row in rows]
    guides = [row["guide"] for row in rows]
    if len(set(ids)) != len(ids) or len(set(guides)) != len(guides):
        raise SourceError("sgRNA library IDs or sequences are not unique")
    if any(not row["id"] or not row["symbol"] for row in rows):
        raise SourceError("sgRNA library has empty IDs or symbols")
    class_counts = {
        name: sum(row["class"] == name for row in rows)
        for name in ("gene", "ctrl_T", "ctrl_NT")
    }
    if class_counts != {"gene": 654, "ctrl_T": 20, "ctrl_NT": 20}:
        raise SourceError(f"sgRNA library class counts changed: {class_counts}")
    target_symbols = {row["symbol"] for row in rows if row["class"] == "gene"}
    if len(target_symbols) != 218:
        raise SourceError("sgRNA target symbol count changed")
    if {row["symbol"] for row in rows if row["class"] == "ctrl_NT"} != {
        CONTROL_SYMBOL
    }:
        raise SourceError("non-targeting control symbol changed")
    if {row["symbol"] for row in rows if row["class"] == "ctrl_T"} != {
        DESERT_SYMBOL
    }:
        raise SourceError("gene-desert control symbol changed")
    return rows, {row["id"]: row for row in rows}


def inspect_channel(
    channel: str,
    path: Path,
    expected: dict[str, Any],
    library_by_id: dict[str, dict[str, str]],
) -> ChannelMetadata:
    with h5py.File(path, "r") as handle:
        required_root = {"X", "obs", "obsm", "raw", "var", "varm"}
        if set(handle.keys()) != required_root:
            raise SourceError(f"{channel} H5AD root members changed")
        obs = handle["obs"]
        var = handle["var"]
        raw = handle["raw"]
        if not isinstance(obs, h5py.Group) or not isinstance(var, h5py.Group):
            raise SourceError(f"{channel} obs/var schema changed")
        if not isinstance(raw, h5py.Group) or set(raw.keys()) != {"X", "var"}:
            raise SourceError(f"{channel} raw schema changed")

        obs_ids = read_frame_index(obs, f"{channel}/obs")
        guide_ids = read_string_vector(obs["sgrna_id"], f"{channel}/obs/sgrna_id")
        symbols = read_string_vector(
            obs["sgrna_symbol"],
            f"{channel}/obs/sgrna_symbol",
        )
        channels = read_string_vector(obs["channel"], f"{channel}/obs/channel")
        cell_lines = read_string_vector(obs["cl"], f"{channel}/obs/cl")
        genes = tuple(
            decode_scalar(value)
            for value in read_string_vector(var["features"], f"{channel}/var/features")
        )
        var_index = tuple(
            decode_scalar(value) for value in read_frame_index(var, f"{channel}/var")
        )
        raw_genes = tuple(
            decode_scalar(value)
            for value in read_frame_index(raw["var"], f"{channel}/raw/var")
        )
        ncount = np.asarray(obs["nCount_RNA"][...])
        nfeature = np.asarray(obs["nFeature_RNA"][...])

        cells = len(obs_ids)
        if cells != expected["cells"] or len(genes) != expected["genes"]:
            raise SourceError(f"{channel} H5AD dimensions changed")
        if not (
            len(guide_ids)
            == len(symbols)
            == len(channels)
            == len(cell_lines)
            == len(ncount)
            == len(nfeature)
            == cells
        ):
            raise SourceError(f"{channel} obs columns disagree")
        validate_axis(tuple(map(str, obs_ids)), f"{channel} obs axis", cells)
        validate_axis(genes, f"{channel} gene axis", expected["genes"])
        if genes != var_index or genes != raw_genes:
            raise SourceError(f"{channel} var and raw.var axes differ")
        if axis_sha256(tuple(map(str, obs_ids))) != expected["obs_axis_sha256"]:
            raise SourceError(f"{channel} obs axis SHA-256 changed")
        if axis_sha256(genes) != expected["gene_axis_sha256"]:
            raise SourceError(f"{channel} gene axis SHA-256 changed")
        if set(map(str, channels)) != {channel} or set(map(str, cell_lines)) != {"A549"}:
            raise SourceError(f"{channel} identity metadata changed")

        guide_ids = np.asarray(list(map(str, guide_ids)), dtype=object)
        symbols = np.asarray(list(map(str, symbols)), dtype=object)
        unknown = sorted(set(guide_ids) - set(library_by_id))
        if unknown:
            raise SourceError(f"{channel} contains unknown guide IDs: {unknown[:3]}")
        for guide_id, symbol in zip(guide_ids, symbols, strict=True):
            if library_by_id[str(guide_id)]["symbol"] != str(symbol):
                raise SourceError(f"{channel} guide-symbol binding changed: {guide_id}")
        if len(set(guide_ids)) != expected["observed_guides"]:
            raise SourceError(f"{channel} observed-guide count changed")
        if int(np.sum(symbols == CONTROL_SYMBOL)) != expected["ntc_cells"]:
            raise SourceError(f"{channel} NTC cell count changed")
        if int(np.sum(symbols == DESERT_SYMBOL)) != expected["gene_desert_cells"]:
            raise SourceError(f"{channel} gene-desert cell count changed")

        raw_x = raw["X"]
        x = handle["X"]
        for label, node in (("X", x), ("raw/X", raw_x)):
            if not isinstance(node, h5py.Group):
                raise SourceError(f"{channel} {label} is not CSR")
            if attr_text(node.attrs.get("encoding-type", "")) != "csr_matrix":
                raise SourceError(f"{channel} {label} encoding changed")
            if tuple(map(int, node.attrs["shape"])) != (cells, len(genes)):
                raise SourceError(f"{channel} {label} shape changed")
            if set(node.keys()) != {"data", "indices", "indptr"}:
                raise SourceError(f"{channel} {label} CSR members changed")
            if node["data"].dtype != np.dtype(np.float64):
                raise SourceError(f"{channel} {label} data dtype changed")
            if node["indices"].dtype != np.dtype(np.int32):
                raise SourceError(f"{channel} {label} index dtype changed")
            if node["indptr"].dtype != np.dtype(np.int32):
                raise SourceError(f"{channel} {label} indptr dtype changed")
        nnz = len(raw_x["data"])
        if nnz != expected["nnz"]:
            raise SourceError(f"{channel} nnz changed")
        if int(np.sum(ncount, dtype=np.float64)) != expected["ncount_sum"]:
            raise SourceError(f"{channel} nCount_RNA sum changed")
        if int(np.sum(nfeature, dtype=np.int64)) != nnz:
            raise SourceError(f"{channel} nFeature_RNA sum does not equal nnz")

    return ChannelMetadata(
        channel=channel,
        path=path,
        obs_ids=np.asarray(obs_ids, dtype=object),
        guide_ids=guide_ids,
        symbols=symbols,
        genes=genes,
        ncount=np.asarray(ncount, dtype=np.float64),
        nfeature=np.asarray(nfeature, dtype=np.int64),
        cells=cells,
        nnz=nnz,
    )


def audit_source(contract: Contract) -> tuple[dict[str, Any], dict[str, ChannelMetadata]]:
    bindings = authenticate_inputs(contract)
    library, library_by_id = load_library(contract.library_path)
    expected = contract.manifest["expected_audit"]
    input_manifest = contract.manifest["inputs"]
    metadata: dict[str, ChannelMetadata] = {}
    for channel, input_key in (
        ("ch_1", "channel_1_h5ad"),
        ("ch_2", "channel_2_h5ad"),
    ):
        channel_expected = {
            **expected["channels"][channel],
            "obs_axis_sha256": input_manifest[input_key]["obs_axis_sha256"],
            "gene_axis_sha256": input_manifest[input_key]["gene_axis_sha256"],
        }
        metadata[channel] = inspect_channel(
            channel,
            contract.channel_paths[channel],
            channel_expected,
            library_by_id,
        )

    observed_guides = set().union(
        *(set(item.guide_ids) for item in metadata.values())
    )
    observed_symbols = set().union(
        *(set(item.symbols) for item in metadata.values())
    )
    target_symbols = {
        row["symbol"] for row in library if row["class"] == "gene"
    }
    panel_overlap = tuple(sorted(target_symbols & set(contract.panel)))
    common_genes = set(metadata["ch_1"].genes) & set(metadata["ch_2"].genes)
    common_canonical = common_genes & set(contract.canonical_genes)
    transferable = common_canonical - set(contract.panel)

    guide_channel_rows: list[tuple[str, str, str, int]] = []
    for channel in CHANNEL_ORDER:
        item = metadata[channel]
        for guide_id in sorted(set(item.guide_ids)):
            indices = np.flatnonzero(item.guide_ids == guide_id)
            symbol = str(item.symbols[indices[0]])
            if library_by_id[guide_id]["class"] == "gene":
                guide_channel_rows.append(
                    (channel, guide_id, symbol, int(len(indices)))
                )
    qualifying = [
        row for row in guide_channel_rows if row[3] >= contract.min_cells
    ]
    qualifying_by_guide: dict[str, list[tuple[str, str, str, int]]] = {}
    for row in qualifying:
        qualifying_by_guide.setdefault(row[1], []).append(row)
    paired_guides = {
        guide_id: rows[0][2]
        for guide_id, rows in qualifying_by_guide.items()
        if {row[0] for row in rows} == set(CHANNEL_ORDER)
        and len(rows) == len(CHANNEL_ORDER)
    }
    paired_by_target: dict[str, list[str]] = {}
    for guide_id, symbol in paired_guides.items():
        paired_by_target.setdefault(symbol, []).append(guide_id)
    passing_targets = tuple(
        sorted(
            symbol
            for symbol, guides in paired_by_target.items()
            if len(guides) >= contract.min_paired_guides
        )
    )
    passing_vcc = tuple(sorted(set(passing_targets) & set(contract.panel)))

    actual = {
        "combined_cells": sum(item.cells for item in metadata.values()),
        "combined_ntc_cells": sum(
            int(np.sum(item.symbols == CONTROL_SYMBOL))
            for item in metadata.values()
        ),
        "combined_gene_desert_cells": sum(
            int(np.sum(item.symbols == DESERT_SYMBOL))
            for item in metadata.values()
        ),
        "observed_guides": len(observed_guides),
        "observed_symbols": len(observed_symbols),
        "observed_guide_channel_target_rows": len(guide_channel_rows),
        "guide_channel_rows_min_10_cells": len(qualifying),
        "paired_guides_min_10_each_channel": len(paired_guides),
        "qc_passing_targets": len(passing_targets),
        "canonical_common_covered_genes": len(common_canonical),
        "canonical_transferable_nonpanel_genes": len(transferable),
        "vcc_exact_overlap": len(panel_overlap),
        "vcc_exact_overlap_symbols": list(panel_overlap),
        "qc_passing_vcc_targets": len(passing_vcc),
        "vcc_qc_failed_targets": sorted(set(panel_overlap) - set(passing_vcc)),
    }
    for key, value in actual.items():
        if expected[key] != value:
            raise SourceError(
                f"source audit changed for {key}: {value!r} != {expected[key]!r}"
            )

    report = {
        "audit": actual,
        "channels": {
            channel: {
                "cells": item.cells,
                "gene_axis_sha256": axis_sha256(item.genes),
                "genes": len(item.genes),
                "gene_desert_cells": int(np.sum(item.symbols == DESERT_SYMBOL)),
                "nnz": item.nnz,
                "ntc_cells": int(np.sum(item.symbols == CONTROL_SYMBOL)),
                "observed_guides": len(set(item.guide_ids)),
            }
            for channel, item in metadata.items()
        },
        "experiment": "120-myllia-a549",
        "implementation_sha256": sha256_file(Path(__file__).resolve()),
        "input_bindings": bindings,
        "manifest_sha256": EXPECTED_MANIFEST_SHA256,
        "provenance_sha256": contract.provenance_sha256,
        "safety": {
            "cell_eval2": False,
            "counts_generated": False,
            "evaluation_run": False,
            "existing_outcome_directories_opened": False,
            "held_truth_opened": False,
            "submission_built": False,
        },
        "schema_version": 1,
    }
    return report, metadata


def atomic_write(path: Path, payload: bytes) -> None:
    safe_child(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def canonical_projection(
    genes: Sequence[str],
    canonical_position: dict[str, int],
) -> tuple[np.ndarray, np.ndarray]:
    destinations = np.asarray(
        [canonical_position.get(gene, -1) for gene in genes],
        dtype=np.int32,
    )
    source_indices = np.flatnonzero(destinations >= 0).astype(np.int64)
    canonical_indices = destinations[source_indices].astype(np.int64)
    if len(np.unique(canonical_indices)) != len(canonical_indices):
        raise SourceError("exact source symbols map multiply to canonical coordinates")
    return source_indices, canonical_indices


def aggregate_channel(
    item: ChannelMetadata,
    library_by_id: dict[str, dict[str, str]],
    canonical_genes: Sequence[str],
    panel_mask: np.ndarray,
    *,
    chunk_rows: int,
) -> ChannelAggregate:
    target_guides = tuple(
        sorted(
            guide_id
            for guide_id in set(item.guide_ids)
            if library_by_id[guide_id]["class"] == "gene"
        )
    )
    guide_position = {guide_id: index + 1 for index, guide_id in enumerate(target_guides)}
    group_for_cell = np.full(item.cells, -1, dtype=np.int32)
    group_for_cell[item.symbols == CONTROL_SYMBOL] = 0
    for guide_id, group in guide_position.items():
        group_for_cell[item.guide_ids == guide_id] = group
    group_sums = np.zeros(
        (len(target_guides) + 1, len(item.genes)),
        dtype=np.uint64,
    )

    with h5py.File(item.path, "r") as handle:
        raw = handle["raw/X"]
        x = handle["X"]
        raw_indptr = np.asarray(raw["indptr"][...], dtype=np.int64)
        x_indptr = np.asarray(x["indptr"][...], dtype=np.int64)
        if not np.array_equal(raw_indptr, x_indptr):
            raise SourceError(f"{item.channel} X/raw.X indptr differ")
        if raw_indptr[0] != 0 or raw_indptr[-1] != item.nnz:
            raise SourceError(f"{item.channel} CSR indptr is invalid")
        if np.any(np.diff(raw_indptr) < 0):
            raise SourceError(f"{item.channel} CSR indptr is not monotonic")

        for start in range(0, item.cells, chunk_rows):
            stop = min(item.cells, start + chunk_rows)
            left = int(raw_indptr[start])
            right = int(raw_indptr[stop])
            raw_data = np.asarray(raw["data"][left:right])
            raw_indices = np.asarray(raw["indices"][left:right])
            x_data = np.asarray(x["data"][left:right])
            x_indices = np.asarray(x["indices"][left:right])
            if not np.array_equal(raw_data, x_data) or not np.array_equal(
                raw_indices,
                x_indices,
            ):
                raise SourceError(f"{item.channel} X and raw.X differ")
            if (
                not np.all(np.isfinite(raw_data))
                or np.any(raw_data < 0)
                or np.any(raw_data != np.floor(raw_data))
            ):
                raise SourceError(f"{item.channel} raw.X is not integer count data")
            local_indptr = raw_indptr[start : stop + 1] - left
            block = sparse.csr_matrix(
                (raw_data, raw_indices, local_indptr),
                shape=(stop - start, len(item.genes)),
            )
            row_sums = np.asarray(block.sum(axis=1)).ravel()
            row_nnz = np.diff(local_indptr)
            if not np.array_equal(row_sums, item.ncount[start:stop]):
                raise SourceError(f"{item.channel} nCount_RNA disagrees with raw.X")
            if not np.array_equal(row_nnz, item.nfeature[start:stop]):
                raise SourceError(f"{item.channel} nFeature_RNA disagrees with raw.X")

            groups = group_for_cell[start:stop]
            selected = np.flatnonzero(groups >= 0)
            if not len(selected):
                continue
            selector = sparse.csr_matrix(
                (
                    np.ones(len(selected), dtype=np.float64),
                    (groups[selected], selected),
                ),
                shape=(len(target_guides) + 1, stop - start),
            )
            partial = (selector @ block).tocoo()
            if np.any(partial.data != np.floor(partial.data)):
                raise SourceError(f"{item.channel} aggregation lost integer counts")
            np.add.at(
                group_sums,
                (partial.row, partial.col),
                partial.data.astype(np.uint64),
            )

    canonical_position = {
        gene: index for index, gene in enumerate(canonical_genes)
    }
    source_indices, destinations = canonical_projection(
        item.genes,
        canonical_position,
    )
    projected = np.zeros(
        (len(target_guides) + 1, len(canonical_genes)),
        dtype=np.float64,
    )
    projected[:, destinations] = group_sums[:, source_indices]
    totals = np.sum(group_sums, axis=1, dtype=np.uint64)
    if np.any(totals == 0):
        raise SourceError(f"{item.channel} contains an empty guide/control pseudobulk")
    normalized = projected * (TARGET_SUM / totals.astype(np.float64))[:, None]
    log_profiles = np.log1p(normalized)
    delta_full = log_profiles[1:] - log_profiles[0]
    if not np.all(np.isfinite(delta_full)):
        raise SourceError(f"{item.channel} produced nonfinite guide effects")

    guide_cells = np.asarray(
        [np.sum(item.guide_ids == guide_id) for guide_id in target_guides],
        dtype=np.int32,
    )
    guide_symbols = tuple(
        library_by_id[guide_id]["symbol"] for guide_id in target_guides
    )
    covered = np.zeros(len(canonical_genes), dtype=np.bool_)
    covered[destinations] = True
    delta_full[:, panel_mask] = 0.0
    return ChannelAggregate(
        channel=item.channel,
        guide_ids=target_guides,
        guide_symbols=guide_symbols,
        guide_cells=guide_cells,
        delta_full=delta_full.astype(np.float32),
        control_log=log_profiles[0].astype(np.float32),
        control_projected_counts=projected[0].astype(np.uint64),
        control_total_counts=int(totals[0]),
        covered_mask=covered,
    )


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denominator == 0.0:
        return 0.0
    return float(np.dot(left, right) / denominator)


def median_pairwise_cosine(values: np.ndarray) -> float:
    if len(values) < 2:
        return 0.0
    scores = [
        cosine(values[left], values[right])
        for left in range(len(values))
        for right in range(left + 1, len(values))
    ]
    return float(np.median(scores))


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def artifact_binding(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if path.suffix == ".npy":
        array = np.load(path, mmap_mode="r", allow_pickle=False)
        result["dtype"] = array.dtype.str
        result["shape"] = list(array.shape)
    return result


def build_compact_donor(
    contract: Contract,
    metadata: dict[str, ChannelMetadata],
    output_dir: Path,
    *,
    chunk_rows: int = DEFAULT_CHUNK_ROWS,
) -> dict[str, Any]:
    output_dir = safe_child(output_dir)
    if output_dir.exists() or output_dir.is_symlink():
        raise SourceError(f"output already exists: {output_dir}")
    if chunk_rows <= 0:
        raise SourceError("chunk_rows must be positive")
    library, library_by_id = load_library(contract.library_path)
    canonical = contract.canonical_genes
    canonical_position = {gene: index for index, gene in enumerate(canonical)}
    panel_mask = np.asarray(
        [gene in set(contract.panel) for gene in canonical],
        dtype=np.bool_,
    )
    if int(np.sum(panel_mask)) != 300:
        raise SourceError("panel does not map exactly to canonical axis")

    aggregates = {
        channel: aggregate_channel(
            metadata[channel],
            library_by_id,
            canonical,
            panel_mask,
            chunk_rows=chunk_rows,
        )
        for channel in CHANNEL_ORDER
    }
    common_covered = np.logical_and.reduce(
        [aggregates[channel].covered_mask for channel in CHANNEL_ORDER]
    )
    transferable = common_covered & ~panel_mask
    expected = contract.manifest["expected_audit"]
    if int(np.sum(common_covered)) != expected["canonical_common_covered_genes"]:
        raise SourceError("common canonical coverage changed during build")
    if int(np.sum(transferable)) != expected[
        "canonical_transferable_nonpanel_genes"
    ]:
        raise SourceError("transferable canonical coverage changed during build")

    by_channel_guide = {
        channel: {
            guide_id: index
            for index, guide_id in enumerate(aggregates[channel].guide_ids)
        }
        for channel in CHANNEL_ORDER
    }
    target_guides = {
        row["symbol"]: []
        for row in library
        if row["class"] == "gene"
    }
    for row in library:
        if row["class"] == "gene":
            target_guides[row["symbol"]].append(row["id"])

    guide_qc_rows: list[list[Any]] = []
    qualifying_channel_rows: list[tuple[str, str, str, int, np.ndarray]] = []
    guide_status: dict[str, dict[str, Any]] = {}
    for target in sorted(target_guides):
        for guide_id in sorted(target_guides[target]):
            counts: dict[str, int] = {}
            observed: dict[str, bool] = {}
            passes: dict[str, bool] = {}
            for channel in CHANNEL_ORDER:
                position = by_channel_guide[channel].get(guide_id)
                is_observed = position is not None
                cells = (
                    int(aggregates[channel].guide_cells[position])
                    if is_observed
                    else 0
                )
                passed = is_observed and cells >= contract.min_cells
                counts[channel] = cells
                observed[channel] = is_observed
                passes[channel] = passed
                if passed:
                    qualifying_channel_rows.append(
                        (
                            channel,
                            guide_id,
                            target,
                            cells,
                            aggregates[channel].delta_full[position],
                        )
                    )
            paired = all(passes.values())
            channel_cosine = ""
            pooled: np.ndarray | None = None
            if paired:
                left = aggregates["ch_1"].delta_full[
                    by_channel_guide["ch_1"][guide_id]
                ][transferable]
                right = aggregates["ch_2"].delta_full[
                    by_channel_guide["ch_2"][guide_id]
                ][transferable]
                channel_cosine = f"{cosine(left, right):.9g}"
                pooled = (
                    aggregates["ch_1"].delta_full[
                        by_channel_guide["ch_1"][guide_id]
                    ].astype(np.float64)
                    + aggregates["ch_2"].delta_full[
                        by_channel_guide["ch_2"][guide_id]
                    ].astype(np.float64)
                ) / 2.0
                pooled[~transferable] = 0.0
            guide_status[guide_id] = {
                "target": target,
                "paired": paired,
                "pooled": pooled,
                "counts": counts,
                "channel_cosine": channel_cosine,
            }
            guide_qc_rows.append(
                [
                    guide_id,
                    target,
                    counts["ch_1"],
                    counts["ch_2"],
                    int(observed["ch_1"]),
                    int(observed["ch_2"]),
                    int(passes["ch_1"]),
                    int(passes["ch_2"]),
                    int(paired),
                    channel_cosine,
                ]
            )

    qualifying_channel_rows.sort(key=lambda row: (row[1], row[0]))
    paired_guide_ids = tuple(
        sorted(guide_id for guide_id, item in guide_status.items() if item["paired"])
    )
    if len(qualifying_channel_rows) != expected["guide_channel_rows_min_10_cells"]:
        raise SourceError("qualifying guide-channel row count changed")
    if len(paired_guide_ids) != expected["paired_guides_min_10_each_channel"]:
        raise SourceError("paired guide count changed")

    guide_delta = np.asarray(
        [guide_status[guide_id]["pooled"] for guide_id in paired_guide_ids],
        dtype=np.float32,
    )
    guide_measured = np.broadcast_to(
        transferable,
        guide_delta.shape,
    ).copy()
    guide_channel_delta = np.asarray(
        [row[4] for row in qualifying_channel_rows],
        dtype=np.float32,
    )
    guide_channel_measured = np.asarray(
        [
            aggregates[row[0]].covered_mask & ~panel_mask
            for row in qualifying_channel_rows
        ],
        dtype=np.bool_,
    )
    guide_channel_delta[~guide_channel_measured] = 0.0

    eligible_by_target: dict[str, list[str]] = {}
    for guide_id in paired_guide_ids:
        eligible_by_target.setdefault(
            str(guide_status[guide_id]["target"]),
            [],
        ).append(guide_id)
    passing_targets = tuple(
        sorted(
            target
            for target, guide_ids in eligible_by_target.items()
            if len(guide_ids) >= contract.min_paired_guides
        )
    )
    if len(passing_targets) != expected["qc_passing_targets"]:
        raise SourceError("target QC count changed")
    target_delta = np.zeros(
        (len(passing_targets), len(canonical)),
        dtype=np.float32,
    )
    target_qc_rows: list[list[Any]] = []
    panel_set = set(contract.panel)
    source_targets = tuple(sorted(target_guides))
    for target in source_targets:
        eligible = sorted(eligible_by_target.get(target, []))
        passed = len(eligible) >= contract.min_paired_guides
        pair_cos = ""
        self_effect = ""
        if eligible:
            values = np.asarray(
                [guide_status[guide_id]["pooled"] for guide_id in eligible],
                dtype=np.float64,
            )
            pair_cos = f"{median_pairwise_cosine(values[:, transferable]):.9g}"
            if target in canonical_position and not panel_mask[canonical_position[target]]:
                self_effect = f"{float(np.median(values[:, canonical_position[target]])):.9g}"
        if passed:
            position = passing_targets.index(target)
            target_delta[position] = np.median(values, axis=0).astype(np.float32)
            target_delta[position, ~transferable] = 0.0
        observed_guides = sum(
            any(guide_id in by_channel_guide[channel] for channel in CHANNEL_ORDER)
            for guide_id in target_guides[target]
        )
        target_qc_rows.append(
            [
                target,
                len(target_guides[target]),
                observed_guides,
                len(eligible),
                int(passed),
                int(target in panel_set),
                pair_cos,
                self_effect,
            ]
        )
    measured_mask = np.broadcast_to(
        transferable,
        target_delta.shape,
    ).copy()
    if np.any(target_delta[:, panel_mask] != 0.0):
        raise SourceError("target delta retained panel coordinates")
    if np.any(guide_delta[:, panel_mask] != 0.0):
        raise SourceError("guide delta retained panel coordinates")
    if np.any(guide_channel_delta[:, panel_mask] != 0.0):
        raise SourceError("guide-channel delta retained panel coordinates")

    channel_b_ctrl = np.asarray(
        [aggregates[channel].control_log for channel in CHANNEL_ORDER],
        dtype=np.float32,
    )
    pooled_control_counts = np.sum(
        np.asarray(
            [
                aggregates[channel].control_projected_counts
                for channel in CHANNEL_ORDER
            ],
            dtype=np.uint64,
        ),
        axis=0,
        dtype=np.uint64,
    )
    pooled_control_total = sum(
        aggregates[channel].control_total_counts for channel in CHANNEL_ORDER
    )
    b_ctrl = np.log1p(
        pooled_control_counts.astype(np.float64)
        * (TARGET_SUM / float(pooled_control_total))
    ).astype(np.float32)
    target_qc_mask = np.asarray(
        [target in set(passing_targets) for target in source_targets],
        dtype=np.bool_,
    )
    vcc_design_overlap = tuple(sorted(set(source_targets) & panel_set))
    vcc_targets = tuple(sorted(set(passing_targets) & panel_set))
    if len(vcc_targets) != expected["qc_passing_vcc_targets"]:
        raise SourceError("QC-passing VCC target count changed")

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(dir=output_dir.parent, prefix=f".{output_dir.name}.")
    )
    safe_child(temporary)
    try:
        np.save(temporary / "b_ctrl.npy", b_ctrl, allow_pickle=False)
        np.save(
            temporary / "channel_b_ctrl.npy",
            channel_b_ctrl,
            allow_pickle=False,
        )
        np.save(
            temporary / "channel_source_covered_gene_mask.npy",
            np.asarray(
                [aggregates[channel].covered_mask for channel in CHANNEL_ORDER],
                dtype=np.bool_,
            ),
            allow_pickle=False,
        )
        np.save(temporary / "delta.npy", target_delta, allow_pickle=False)
        np.save(
            temporary / "guide_channel_delta.npy",
            guide_channel_delta,
            allow_pickle=False,
        )
        np.save(
            temporary / "guide_channel_measured_mask.npy",
            guide_channel_measured,
            allow_pickle=False,
        )
        np.save(
            temporary / "guide_delta.npy",
            guide_delta,
            allow_pickle=False,
        )
        np.save(
            temporary / "guide_measured_mask.npy",
            guide_measured,
            allow_pickle=False,
        )
        np.save(
            temporary / "measured_mask.npy",
            measured_mask,
            allow_pickle=False,
        )
        np.save(
            temporary / "panel_coordinate_mask.npy",
            panel_mask,
            allow_pickle=False,
        )
        np.save(
            temporary / "source_covered_gene_mask.npy",
            common_covered,
            allow_pickle=False,
        )
        np.save(
            temporary / "target_qc_mask.npy",
            target_qc_mask,
            allow_pickle=False,
        )
        np.save(
            temporary / "transferable_gene_mask.npy",
            transferable,
            allow_pickle=False,
        )
        write_lines(temporary / "genes.txt", canonical)
        write_lines(temporary / "perts.txt", passing_targets)
        write_lines(temporary / "source_targets.txt", source_targets)
        write_lines(
            temporary / "vcc_design_overlap.txt",
            vcc_design_overlap,
        )
        write_lines(temporary / "vcc_targets.txt", vcc_targets)
        write_tsv(
            temporary / "guide_channel_rows.tsv",
            ["row", "channel", "guide_id", "target", "cells"],
            [
                [index, channel, guide_id, target, cells]
                for index, (channel, guide_id, target, cells, _) in enumerate(
                    qualifying_channel_rows
                )
            ],
        )
        write_tsv(
            temporary / "guide_rows.tsv",
            [
                "row",
                "guide_id",
                "target",
                "ch_1_cells",
                "ch_2_cells",
                "channel_cosine",
            ],
            [
                [
                    index,
                    guide_id,
                    guide_status[guide_id]["target"],
                    guide_status[guide_id]["counts"]["ch_1"],
                    guide_status[guide_id]["counts"]["ch_2"],
                    guide_status[guide_id]["channel_cosine"],
                ]
                for index, guide_id in enumerate(paired_guide_ids)
            ],
        )
        write_tsv(
            temporary / "guide_qc.tsv",
            [
                "guide_id",
                "target",
                "ch_1_cells",
                "ch_2_cells",
                "ch_1_observed",
                "ch_2_observed",
                "ch_1_min_cells_pass",
                "ch_2_min_cells_pass",
                "paired_pass",
                "channel_cosine",
            ],
            guide_qc_rows,
        )
        write_tsv(
            temporary / "target_qc.tsv",
            [
                "target",
                "library_guides",
                "observed_guides",
                "paired_guides",
                "target_qc_pass",
                "vcc_panel_overlap",
                "median_pairwise_guide_cosine",
                "median_self_effect_nonpanel_only",
            ],
            target_qc_rows,
        )
        implementation_hash = sha256_file(Path(__file__).resolve())
        meta = {
            "array_contract": {
                "b_ctrl.npy": {
                    "dtype": b_ctrl.dtype.str,
                    "shape": list(b_ctrl.shape),
                    "semantics": "pooled A549 NTC TS50000 log1p profile",
                },
                "delta.npy": {
                    "dtype": target_delta.dtype.str,
                    "shape": list(target_delta.shape),
                    "semantics": "median replicated-guide target effect; panel coordinates zero",
                },
                "guide_channel_delta.npy": {
                    "dtype": guide_channel_delta.dtype.str,
                    "shape": list(guide_channel_delta.shape),
                    "semantics": "qualifying source guide-channel effect; panel coordinates zero",
                },
                "guide_delta.npy": {
                    "dtype": guide_delta.dtype.str,
                    "shape": list(guide_delta.shape),
                    "semantics": "equal-channel pooled replicated-guide effect; panel coordinates zero",
                },
                "measured_mask.npy": {
                    "dtype": measured_mask.dtype.str,
                    "shape": list(measured_mask.shape),
                    "semantics": "common exact source coverage excluding all VCC panel coordinates",
                },
            },
            "axes": {
                "canonical_genes": len(canonical),
                "channel_order": list(CHANNEL_ORDER),
                "common_covered_genes": int(np.sum(common_covered)),
                "guide_channel_rows": len(qualifying_channel_rows),
                "guide_rows": len(paired_guide_ids),
                "panel_coordinates_excluded": int(np.sum(panel_mask)),
                "source_target_genes": len(target_guides),
                "target_rows": len(passing_targets),
                "target_qc_mask_axis": "source_targets.txt",
                "transferable_genes": int(np.sum(transferable)),
                "vcc_design_overlap": expected["vcc_exact_overlap"],
                "vcc_qc_passing_targets": len(vcc_targets),
            },
            "conversion": {
                "channel_pooling": "equal_weight_mean",
                "control_pooling": "raw-count pooled across both channel NTC sets",
                "exact_symbols_only": True,
                "guide_channel_min_cells": contract.min_cells,
                "guide_requires_both_channels": True,
                "minimum_paired_guides_per_target": contract.min_paired_guides,
                "panel_coordinates_transferable": False,
                "pseudobulk": (
                    "sum raw counts; normalize over full source feature axis "
                    "to TS50000; exact-symbol project without renormalizing "
                    "unmapped mass; log1p; subtract channel NTC"
                ),
                "target_pooling": "coordinate-wise median of channel-pooled guides",
            },
            "created_by": "experiments/120-myllia-a549/build_compact_donor.py",
            "format_version": 1,
            "implementation_sha256": implementation_hash,
            "safety": {
                "cell_eval2": False,
                "counts_generated": False,
                "evaluation_run": False,
                "existing_outcome_directories_opened": False,
                "held_truth_opened": False,
                "submission_built": False,
            },
            "source": {
                "dataset": "GSE337804",
                "manifest_sha256": EXPECTED_MANIFEST_SHA256,
                "provenance_sha256": contract.provenance_sha256,
                "source_h5ad_sha256": {
                    channel: contract.manifest["inputs"][
                        f"channel_{1 if channel == 'ch_1' else 2}_h5ad"
                    ]["sha256"]
                    for channel in CHANNEL_ORDER
                },
            },
        }
        (temporary / "meta.json").write_bytes(json_bytes(meta))

        artifact_names = sorted((*ARRAY_FILES, *TEXT_FILES))
        manifest = {
            "artifact_set": "GSE337804_A549_compact_source_donor",
            "experiment": "120-myllia-a549",
            "files": {
                name: artifact_binding(temporary / name)
                for name in artifact_names
            },
            "schema_version": 1,
        }
        (temporary / ARTIFACT_MANIFEST_NAME).write_bytes(json_bytes(manifest))
        receipt = {
            "artifact_manifest_sha256": sha256_file(
                temporary / ARTIFACT_MANIFEST_NAME
            ),
            "experiment": "120-myllia-a549",
            "implementation_sha256": implementation_hash,
            "manifest_sha256": EXPECTED_MANIFEST_SHA256,
            "provenance_sha256": contract.provenance_sha256,
            "schema_version": 1,
            "source_sha256": {
                key: value["sha256"]
                for key, value in contract.manifest["inputs"].items()
            },
            "status": "verified_source_only_compact_donor",
        }
        (temporary / RECEIPT_NAME).write_bytes(json_bytes(receipt))
        receipt_hash = sha256_file(temporary / RECEIPT_NAME)
        (temporary / RECEIPT_HASH_NAME).write_text(
            f"{receipt_hash}  {RECEIPT_NAME}\n",
            encoding="utf-8",
        )
        temporary.replace(output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    verify_compact_donor(output_dir, contract=contract)
    return meta


def verify_compact_donor(
    output_dir: Path,
    *,
    contract: Contract | None = None,
) -> dict[str, Any]:
    output_dir = safe_child(output_dir)
    if output_dir.is_symlink() or not output_dir.is_dir():
        raise SourceError(f"compact donor is missing or symlinked: {output_dir}")
    required = {
        *ARRAY_FILES,
        *TEXT_FILES,
        ARTIFACT_MANIFEST_NAME,
        RECEIPT_NAME,
        RECEIPT_HASH_NAME,
    }
    if {path.name for path in output_dir.iterdir()} != required:
        raise SourceError("compact donor file set changed")
    manifest = json.loads(
        (output_dir / ARTIFACT_MANIFEST_NAME).read_text(encoding="utf-8")
    )
    if manifest.get("artifact_set") != "GSE337804_A549_compact_source_donor":
        raise SourceError("artifact manifest identity changed")
    for name, expected in manifest["files"].items():
        actual = artifact_binding(output_dir / name)
        if actual != expected:
            raise SourceError(f"artifact binding changed: {name}")
    receipt_path = output_dir / RECEIPT_NAME
    receipt_hash_line = (output_dir / RECEIPT_HASH_NAME).read_text(
        encoding="utf-8"
    )
    expected_line = f"{sha256_file(receipt_path)}  {RECEIPT_NAME}\n"
    if receipt_hash_line != expected_line:
        raise SourceError("compact receipt hash sidecar changed")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("artifact_manifest_sha256") != sha256_file(
        output_dir / ARTIFACT_MANIFEST_NAME
    ):
        raise SourceError("receipt artifact-manifest binding changed")
    if receipt.get("manifest_sha256") != EXPECTED_MANIFEST_SHA256:
        raise SourceError("receipt source-manifest binding changed")
    if receipt.get("provenance_sha256") != EXPECTED_PROVENANCE_SHA256:
        raise SourceError("receipt provenance binding changed")
    implementation_hash = sha256_file(Path(__file__).resolve())
    if receipt.get("implementation_sha256") != implementation_hash:
        raise SourceError("receipt implementation binding changed")
    active_contract = contract if contract is not None else load_contract()
    expected_sources = {
        key: value["sha256"]
        for key, value in active_contract.manifest["inputs"].items()
    }
    if receipt.get("source_sha256") != expected_sources:
        raise SourceError("receipt source bindings changed")

    meta = json.loads((output_dir / "meta.json").read_text(encoding="utf-8"))
    if meta.get("implementation_sha256") != implementation_hash:
        raise SourceError("meta implementation binding changed")
    panel_mask = np.load(
        output_dir / "panel_coordinate_mask.npy",
        allow_pickle=False,
    )
    transferable = np.load(
        output_dir / "transferable_gene_mask.npy",
        allow_pickle=False,
    )
    if int(np.sum(panel_mask)) != 300 or np.any(transferable & panel_mask):
        raise SourceError("panel exclusion mask changed")
    for name in (
        "delta.npy",
        "guide_delta.npy",
        "guide_channel_delta.npy",
    ):
        values = np.load(output_dir / name, mmap_mode="r", allow_pickle=False)
        if not np.all(np.isfinite(values)):
            raise SourceError(f"{name} contains nonfinite values")
        if np.any(values[:, panel_mask] != 0.0):
            raise SourceError(f"{name} contains nonzero panel coordinates")
    for delta_name, mask_name in (
        ("delta.npy", "measured_mask.npy"),
        ("guide_delta.npy", "guide_measured_mask.npy"),
        ("guide_channel_delta.npy", "guide_channel_measured_mask.npy"),
    ):
        values = np.load(output_dir / delta_name, mmap_mode="r", allow_pickle=False)
        measured = np.load(output_dir / mask_name, mmap_mode="r", allow_pickle=False)
        if values.shape != measured.shape or measured.dtype != np.bool_:
            raise SourceError(f"{mask_name} contract changed")
        if np.any(values[~measured] != 0.0):
            raise SourceError(f"{delta_name} has nonzero unsupported coordinates")
    return meta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser("audit")
    audit.add_argument("--write-report", action="store_true")
    build = subparsers.add_parser("build")
    build.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    build.add_argument("--chunk-rows", type=int, default=DEFAULT_CHUNK_ROWS)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = load_contract()
    if args.command == "verify":
        meta = verify_compact_donor(args.output_dir, contract=contract)
        print(json.dumps(meta["axes"], indent=2, sort_keys=True))
        return
    report, metadata = audit_source(contract)
    if args.command == "audit":
        if args.write_report:
            atomic_write(DEFAULT_AUDIT_PATH, json_bytes(report))
        print(json.dumps(report["audit"], indent=2, sort_keys=True))
        return
    meta = build_compact_donor(
        contract,
        metadata,
        args.output_dir,
        chunk_rows=args.chunk_rows,
    )
    print(json.dumps(meta["axes"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
