#!/usr/bin/env python3
"""Build a canonical-axis compact source from the official Altair A549 H5AD.

The output preserves source log2-fold-change and FDR values. It does not
convert them into project-space effects, fit a model, inspect held truth,
evaluate predictions, generate counts, or build a submission.
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

import source_audit


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUTPUT_DIR = HERE / "work" / "compact-source"
SOURCE_MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
AUDIT_PATH = HERE / "SOURCE_AUDIT.json"
ACQUISITION_RECEIPT_PATH = HERE / "SOURCE_ACQUISITION_RECEIPT.json"

EXPECTED_SOURCE_MANIFEST_SHA256 = (
    "61c1b6338ad197d34f28162766c7e5cf8e4b1be0324fe2e974015152285fff88"
)
EXPECTED_CANONICAL_GENES_SHA256 = (
    "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"
)
EXPECTED_PANEL_SHA256 = (
    "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
)

ARTIFACT_MANIFEST_NAME = "artifact-manifest.json"
COMPACT_RECEIPT_NAME = "compact-source-receipt.json"
COMPACT_RECEIPT_HASH_NAME = "compact-source-receipt.sha256"
ARTIFACT_FILES = (
    "fdr_mask.npy",
    "genes.txt",
    "log2_fold_change.npy",
    "measured_mask.npy",
    "meta.json",
    "neg_log10_fdr.npy",
    "perts.txt",
    "rows.tsv",
    "source_covered_gene_mask.npy",
    "source_gene_map.tsv",
    "source_supported_row_mask.npy",
    "vcc_supported_target_mask.npy",
    "vcc_supported_targets.txt",
    "vcc_target_mask.npy",
    "vcc_targets.txt",
)


class ExtractionError(RuntimeError):
    """A compact-source provenance, schema, or output invariant failed."""


@dataclass(frozen=True)
class ExtractionContract:
    source_h5ad: Path
    source_manifest: Path
    source_audit: Path
    acquisition_receipt: Path
    canonical_genes_csv: Path
    panel_csv: Path
    expected_source_size: int
    expected_source_sha256: str
    expected_source_manifest_sha256: str
    expected_canonical_genes_sha256: str
    expected_panel_sha256: str
    expected_source_rows: int
    expected_target_rows: int
    expected_source_genes: int
    expected_canonical_genes: int
    expected_panel_targets: int
    expected_vcc_overlap: int
    expected_vcc_supported_overlap: int


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ExtractionError(message)


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
        raise ExtractionError(f"cannot read valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ExtractionError(f"JSON root must be an object: {path}")
    return value


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_lines(path: Path, values: Iterable[str]) -> None:
    path.write_text(
        "".join(f"{value}\n" for value in values),
        encoding="utf-8",
    )


def validate_axis(
    values: Sequence[str],
    label: str,
    expected_count: int,
) -> None:
    _require(len(values) == expected_count, f"{label} count changed")
    _require(
        all(value and value == value.strip() for value in values),
        f"{label} contains empty or padded labels",
    )
    _require(len(values) == len(set(values)), f"{label} contains duplicates")


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


def read_text_axis(path: Path, expected_count: int) -> tuple[str, ...]:
    _require(path.read_bytes().endswith(b"\n"), f"{path.name} lacks final newline")
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    validate_axis(values, path.name, expected_count)
    return values


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    _require(
        resolved != HERE and resolved.is_relative_to(HERE),
        f"output must be a child of {HERE}",
    )
    return resolved


def production_contract() -> ExtractionContract:
    manifest = source_audit.load_manifest(SOURCE_MANIFEST_PATH)
    manifest_hash = sha256_file(SOURCE_MANIFEST_PATH)
    if EXPECTED_SOURCE_MANIFEST_SHA256 != "TO_BE_PINNED":
        _require(
            manifest_hash == EXPECTED_SOURCE_MANIFEST_SHA256,
            "SOURCE_MANIFEST.json changed",
        )

    genes = manifest["canonical_axes"]["genes"]
    panel = manifest["canonical_axes"]["panel"]
    _require(
        genes["sha256"] == EXPECTED_CANONICAL_GENES_SHA256,
        "canonical gene hash declaration changed",
    )
    _require(
        panel["sha256"] == EXPECTED_PANEL_SHA256,
        "panel hash declaration changed",
    )
    expected = manifest["expected_audit"]
    return ExtractionContract(
        source_h5ad=(ROOT / manifest["source"]["local_path"]).resolve(),
        source_manifest=SOURCE_MANIFEST_PATH.resolve(),
        source_audit=AUDIT_PATH.resolve(),
        acquisition_receipt=ACQUISITION_RECEIPT_PATH.resolve(),
        canonical_genes_csv=(ROOT / genes["path"]).resolve(),
        panel_csv=(ROOT / panel["path"]).resolve(),
        expected_source_size=manifest["source"]["size_bytes"],
        expected_source_sha256=manifest["source"]["sha256"],
        expected_source_manifest_sha256=manifest_hash,
        expected_canonical_genes_sha256=genes["sha256"],
        expected_panel_sha256=panel["sha256"],
        expected_source_rows=expected["source_shape"][0],
        expected_target_rows=expected["target_rows"],
        expected_source_genes=expected["source_shape"][1],
        expected_canonical_genes=genes["count"],
        expected_panel_targets=panel["count"],
        expected_vcc_overlap=expected["vcc_panel_target_overlap"],
        expected_vcc_supported_overlap=expected["vcc_panel_supported_overlap"],
    )


def authenticate_inputs(
    contract: ExtractionContract,
) -> tuple[dict[str, Any], dict[str, Any]]:
    report, receipt = source_audit.verify_frozen_reports()
    _require(
        report["source"]["size_bytes"] == contract.expected_source_size,
        "audit source size binding changed",
    )
    _require(
        report["source"]["sha256"] == contract.expected_source_sha256,
        "audit source hash binding changed",
    )
    _require(
        receipt["source_manifest_sha256"] == contract.expected_source_manifest_sha256,
        "acquisition receipt manifest binding changed",
    )
    _require(
        sha256_file(contract.canonical_genes_csv)
        == contract.expected_canonical_genes_sha256,
        "canonical gene axis hash drift",
    )
    _require(
        sha256_file(contract.panel_csv) == contract.expected_panel_sha256,
        "VCC panel axis hash drift",
    )
    return report, receipt


def array_binding(array: np.ndarray) -> dict[str, Any]:
    view = np.asarray(array)
    return {
        "dtype": view.dtype.str,
        "shape": list(view.shape),
        "c_order": bool(view.flags.c_contiguous),
    }


def write_rows(
    path: Path,
    source_rows: np.ndarray,
    source_gene_ids: np.ndarray,
    n_guides: np.ndarray,
    n_cells: np.ndarray,
    measured_counts: np.ndarray,
    canonical_measured_counts: np.ndarray,
    vcc_mask: np.ndarray,
    supported_mask: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "compact_row",
                "source_row",
                "perturbation",
                "source_gene_id",
                "n_guides",
                "n_cells",
                "source_measured_genes",
                "canonical_measured_genes",
                "source_supported",
                "vcc_target",
                "vcc_supported_target",
            ]
        )
        for compact_row, source_row in enumerate(source_rows):
            writer.writerow(
                [
                    compact_row,
                    int(source_row),
                    source_gene_ids[compact_row, 0],
                    source_gene_ids[compact_row, 1],
                    int(n_guides[compact_row]),
                    int(n_cells[compact_row]),
                    int(measured_counts[compact_row]),
                    int(canonical_measured_counts[compact_row]),
                    int(supported_mask[compact_row]),
                    int(vcc_mask[compact_row]),
                    int(vcc_mask[compact_row] and supported_mask[compact_row]),
                ]
            )


def write_gene_map(
    path: Path,
    source_ids: Sequence[str],
    source_symbols: Sequence[str],
    source_to_canonical: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "source_column",
                "source_ensembl_id",
                "source_gene_symbol",
                "canonical_column",
                "mapping",
            ]
        )
        for source_column, (ensembl_id, symbol, destination) in enumerate(
            zip(source_ids, source_symbols, source_to_canonical, strict=True)
        ):
            writer.writerow(
                [
                    source_column,
                    ensembl_id,
                    symbol,
                    "" if destination < 0 else int(destination),
                    "unmapped" if destination < 0 else "exact_symbol",
                ]
            )


def build_artifact_manifest(
    directory: Path,
    contract: ExtractionContract,
    report: dict[str, Any],
) -> dict[str, Any]:
    files = {}
    for name in ARTIFACT_FILES:
        path = directory / name
        _require(path.is_file() and not path.is_symlink(), f"missing artifact {name}")
        files[name] = {
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return {
        "schema_version": 1,
        "experiment": "119-altair-a549",
        "artifact": "canonical_compact_source",
        "source": {
            "path": str(contract.source_h5ad.relative_to(ROOT)),
            "size_bytes": contract.expected_source_size,
            "sha256": contract.expected_source_sha256,
        },
        "bindings": {
            "source_manifest_sha256": contract.expected_source_manifest_sha256,
            "source_audit_sha256": sha256_file(contract.source_audit),
            "source_acquisition_receipt_sha256": sha256_file(
                contract.acquisition_receipt
            ),
            "source_audit_implementation_sha256": sha256_file(HERE / "source_audit.py"),
            "extractor_sha256": sha256_file(Path(__file__).resolve()),
            "canonical_genes_file_sha256": contract.expected_canonical_genes_sha256,
            "vcc_panel_file_sha256": contract.expected_panel_sha256,
            "compact_perturbations_sha256": sha256_file(directory / "perts.txt"),
            "compact_canonical_genes_sha256": sha256_file(directory / "genes.txt"),
            "source_perturbations_sha256": report["axis_bindings"][
                "source_perturbations_sha256"
            ],
            "source_ensembl_ids_sha256": report["axis_bindings"][
                "source_ensembl_ids_sha256"
            ],
            "source_gene_symbols_sha256": report["axis_bindings"][
                "source_gene_symbols_sha256"
            ],
        },
        "semantics": {
            "log2_fold_change": "unaltered source X values",
            "neg_log10_fdr": "unaltered source layers/neg_log10_fdr values",
            "measured_mask": "finite source X after exact-symbol projection",
            "fdr_mask": (
                "finite source layers/neg_log10_fdr after exact-symbol projection"
            ),
            "unsupported": "NaN with false mask",
            "mapping": "exact_symbol_only_no_alias_merge",
            "row_scope": "1000 targeting aggregates; pooled control excluded",
        },
        "files": files,
    }


def build_receipt(
    directory: Path,
    artifact_manifest: dict[str, Any],
    contract: ExtractionContract,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "experiment": "119-altair-a549",
        "status": "verified_compact_source",
        "source_manifest_sha256": contract.expected_source_manifest_sha256,
        "source_audit_sha256": sha256_file(contract.source_audit),
        "source_acquisition_receipt_sha256": sha256_file(contract.acquisition_receipt),
        "source_sha256": contract.expected_source_sha256,
        "artifact_manifest_sha256": sha256_file(directory / ARTIFACT_MANIFEST_NAME),
        "artifact_file_count": len(artifact_manifest["files"]),
        "authorization": {
            "source_only": True,
            "held_truth_accessed": False,
            "existing_outcomes_accessed": False,
            "model_fitted": False,
            "evaluation_run": False,
            "counts_generated": False,
            "submission_built": False,
        },
    }


def extract_compact_source(
    contract: ExtractionContract,
    output_dir: Path,
    chunk_rows: int = 32,
) -> Path:
    output = safe_output(output_dir)
    _require(chunk_rows > 0, "chunk_rows must be positive")
    _require(not output.exists(), f"refusing to overwrite existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    report, _ = authenticate_inputs(contract)
    canonical = read_csv_axis(
        contract.canonical_genes_csv,
        "gene_name",
        contract.expected_canonical_genes,
    )
    panel = read_csv_axis(
        contract.panel_csv,
        "target_gene",
        contract.expected_panel_targets,
    )
    panel_set = set(panel)

    staging = Path(
        tempfile.mkdtemp(prefix=".compact-source-", dir=output.parent)
    ).resolve()
    try:
        with h5py.File(contract.source_h5ad, "r") as handle:
            obs_names = source_audit.read_dataframe_index(handle, "obs")
            roles = source_audit.read_string_node(handle["obs/role"], "/obs/role")
            obs_gene_ids = source_audit.read_string_node(
                handle["obs/gene_id"],
                "/obs/gene_id",
            )
            source_ids = source_audit.read_dataframe_index(handle, "var")
            source_symbols = source_audit.read_string_node(
                handle["var/gene_symbol"],
                "/var/gene_symbol",
            )
            n_guides_all = np.asarray(handle["obs/n_guides"][...], dtype=np.int64)
            n_cells_all = np.asarray(handle["obs/n_cells"][...], dtype=np.int64)

            validate_axis(
                obs_names.tolist(),
                "source perturbations",
                contract.expected_source_rows,
            )
            validate_axis(
                source_ids.tolist(),
                "source Ensembl ids",
                contract.expected_source_genes,
            )
            validate_axis(
                source_symbols.tolist(),
                "source symbols",
                contract.expected_source_genes,
            )

            target_rows = np.flatnonzero(roles == "targeting")
            _require(
                len(target_rows) == contract.expected_target_rows,
                "target row count changed",
            )
            perturbations = tuple(str(obs_names[index]) for index in target_rows)
            validate_axis(
                perturbations,
                "target perturbations",
                contract.expected_target_rows,
            )
            target_gene_ids = np.asarray(
                [
                    [str(obs_names[index]), str(obs_gene_ids[index])]
                    for index in target_rows
                ],
                dtype=object,
            )
            n_guides = n_guides_all[target_rows]
            n_cells = n_cells_all[target_rows]

            source_to_canonical = source_audit.exact_symbol_projection(
                source_symbols.tolist(),
                canonical,
            )
            mapped_source = np.flatnonzero(source_to_canonical >= 0)
            mapped_destination = source_to_canonical[mapped_source].astype(
                np.int64,
                copy=False,
            )
            _require(
                len(mapped_destination) == len(set(mapped_destination.tolist())),
                "exact mapping unexpectedly merges source symbols",
            )

            shape = (contract.expected_target_rows, contract.expected_canonical_genes)
            lfc = np.lib.format.open_memmap(
                staging / "log2_fold_change.npy",
                mode="w+",
                dtype=np.float32,
                shape=shape,
            )
            fdr = np.lib.format.open_memmap(
                staging / "neg_log10_fdr.npy",
                mode="w+",
                dtype=np.float32,
                shape=shape,
            )
            measured = np.lib.format.open_memmap(
                staging / "measured_mask.npy",
                mode="w+",
                dtype=np.bool_,
                shape=shape,
            )
            fdr_mask = np.lib.format.open_memmap(
                staging / "fdr_mask.npy",
                mode="w+",
                dtype=np.bool_,
                shape=shape,
            )
            lfc[:] = np.nan
            fdr[:] = np.nan
            measured[:] = False
            fdr_mask[:] = False

            source_measured_counts = np.zeros(len(target_rows), dtype=np.int64)
            canonical_measured_counts = np.zeros(len(target_rows), dtype=np.int64)
            source_x = source_audit.require_dense_float32(
                handle,
                "X",
                (contract.expected_source_rows, contract.expected_source_genes),
            )
            source_fdr = source_audit.require_dense_float32(
                handle,
                "layers/neg_log10_fdr",
                (contract.expected_source_rows, contract.expected_source_genes),
            )

            for start in range(0, len(target_rows), chunk_rows):
                stop = min(len(target_rows), start + chunk_rows)
                selected_rows = target_rows[start:stop]
                source_lfc_chunk = np.asarray(
                    source_x[selected_rows, :],
                    dtype=np.float32,
                )
                source_fdr_chunk = np.asarray(
                    source_fdr[selected_rows, :],
                    dtype=np.float32,
                )
                source_measured_counts[start:stop] = np.isfinite(source_lfc_chunk).sum(
                    axis=1
                )
                projected_lfc = source_lfc_chunk[:, mapped_source]
                projected_fdr = source_fdr_chunk[:, mapped_source]
                lfc[start:stop, mapped_destination] = projected_lfc
                fdr[start:stop, mapped_destination] = projected_fdr
                measured_chunk = np.isfinite(projected_lfc)
                fdr_chunk_mask = np.isfinite(projected_fdr)
                measured[start:stop, mapped_destination] = measured_chunk
                fdr_mask[start:stop, mapped_destination] = fdr_chunk_mask
                canonical_measured_counts[start:stop] = measured_chunk.sum(axis=1)

            lfc.flush()
            fdr.flush()
            measured.flush()
            fdr_mask.flush()
            del lfc, fdr, measured, fdr_mask

        source_supported = source_measured_counts > 0
        source_symbol_set = set(source_symbols.tolist())
        source_covered_gene_mask = np.asarray(
            [symbol in source_symbol_set for symbol in canonical],
            dtype=bool,
        )
        vcc_target_mask = np.asarray(
            [perturbation in panel_set for perturbation in perturbations],
            dtype=bool,
        )
        vcc_supported_mask = vcc_target_mask & source_supported
        _require(
            int(source_covered_gene_mask.sum())
            == report["summary"]["canonical_exact_overlap"],
            "canonical coverage count changed",
        )
        _require(
            int(vcc_target_mask.sum()) == contract.expected_vcc_overlap,
            "VCC target overlap changed",
        )
        _require(
            int(vcc_supported_mask.sum()) == contract.expected_vcc_supported_overlap,
            "VCC supported target overlap changed",
        )

        np.save(
            staging / "source_covered_gene_mask.npy",
            source_covered_gene_mask,
            allow_pickle=False,
        )
        np.save(
            staging / "source_supported_row_mask.npy",
            source_supported,
            allow_pickle=False,
        )
        np.save(
            staging / "vcc_target_mask.npy",
            vcc_target_mask,
            allow_pickle=False,
        )
        np.save(
            staging / "vcc_supported_target_mask.npy",
            vcc_supported_mask,
            allow_pickle=False,
        )
        write_lines(staging / "genes.txt", canonical)
        write_lines(staging / "perts.txt", perturbations)
        write_lines(
            staging / "vcc_targets.txt",
            [p for p, keep in zip(perturbations, vcc_target_mask, strict=True) if keep],
        )
        write_lines(
            staging / "vcc_supported_targets.txt",
            [
                p
                for p, keep in zip(
                    perturbations,
                    vcc_supported_mask,
                    strict=True,
                )
                if keep
            ],
        )
        write_rows(
            staging / "rows.tsv",
            target_rows,
            target_gene_ids,
            n_guides,
            n_cells,
            source_measured_counts,
            canonical_measured_counts,
            vcc_target_mask,
            source_supported,
        )
        write_gene_map(
            staging / "source_gene_map.tsv",
            source_ids.tolist(),
            source_symbols.tolist(),
            source_to_canonical,
        )

        lfc_view = np.load(
            staging / "log2_fold_change.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        fdr_view = np.load(
            staging / "neg_log10_fdr.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        measured_view = np.load(
            staging / "measured_mask.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        fdr_mask_view = np.load(
            staging / "fdr_mask.npy",
            mmap_mode="r",
            allow_pickle=False,
        )
        _require(np.all(np.isnan(lfc_view[~measured_view])), "LFC NaNs not preserved")
        _require(np.all(np.isfinite(lfc_view[measured_view])), "LFC mask is invalid")
        _require(np.all(np.isnan(fdr_view[~fdr_mask_view])), "FDR NaNs not preserved")
        _require(np.all(np.isfinite(fdr_view[fdr_mask_view])), "FDR mask is invalid")
        _require(
            not np.any(fdr_mask_view & ~measured_view),
            "FDR support exceeds LFC support",
        )

        meta = {
            "schema_version": 1,
            "experiment": "119-altair-a549",
            "artifact": "canonical_compact_source",
            "source": {
                "dataset": "Altair - A549 1k scRNA",
                "cell_line": "A549",
                "perturbation_type": "CRISPR knockout",
                "matrix_semantics": "log2_fold_change",
                "source_sha256": contract.expected_source_sha256,
            },
            "axes": {
                "perturbations": len(perturbations),
                "canonical_genes": len(canonical),
                "source_response_genes": len(source_symbols),
                "source_exact_symbol_matches": int(source_covered_gene_mask.sum()),
                "source_unmapped_symbols": len(source_symbols)
                - int(source_covered_gene_mask.sum()),
                "vcc_target_rows": int(vcc_target_mask.sum()),
                "vcc_supported_target_rows": int(vcc_supported_mask.sum()),
                "source_supported_rows": int(source_supported.sum()),
            },
            "arrays": {
                "log2_fold_change.npy": array_binding(lfc_view),
                "neg_log10_fdr.npy": array_binding(fdr_view),
                "measured_mask.npy": array_binding(measured_view),
                "fdr_mask.npy": array_binding(fdr_mask_view),
            },
            "mapping": {
                "method": "exact_symbol_only",
                "alias_resolution": False,
                "alias_merge": False,
                "unmapped_source_values_discarded": True,
                "uncovered_canonical_values": "NaN",
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
        write_json(staging / "meta.json", meta)

        artifact_manifest = build_artifact_manifest(staging, contract, report)
        write_json(staging / ARTIFACT_MANIFEST_NAME, artifact_manifest)
        receipt = build_receipt(staging, artifact_manifest, contract)
        write_json(staging / COMPACT_RECEIPT_NAME, receipt)
        receipt_hash = sha256_file(staging / COMPACT_RECEIPT_NAME)
        write_lines(staging / COMPACT_RECEIPT_HASH_NAME, [receipt_hash])
        verify_compact_source(staging, contract)

        os.replace(staging, output)
        return output
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def verify_compact_source(
    directory: Path,
    contract: ExtractionContract | None = None,
) -> dict[str, Any]:
    directory = directory.resolve()
    _require(directory.is_dir(), f"compact source is missing: {directory}")
    manifest_path = directory / ARTIFACT_MANIFEST_NAME
    receipt_path = directory / COMPACT_RECEIPT_NAME
    receipt_hash_path = directory / COMPACT_RECEIPT_HASH_NAME
    for path in (manifest_path, receipt_path, receipt_hash_path):
        _require(path.is_file() and not path.is_symlink(), f"missing {path.name}")

    artifact = load_json(manifest_path)
    receipt = load_json(receipt_path)
    _require(
        artifact.get("experiment") == "119-altair-a549",
        "artifact experiment changed",
    )
    _require(
        set(artifact.get("files", {})) == set(ARTIFACT_FILES),
        "artifact file inventory changed",
    )
    for name, binding in artifact["files"].items():
        path = directory / name
        _require(path.is_file() and not path.is_symlink(), f"missing artifact {name}")
        _require(path.stat().st_size == binding["size_bytes"], f"{name} size drift")
        _require(sha256_file(path) == binding["sha256"], f"{name} hash drift")

    receipt_hash = sha256_file(receipt_path)
    _require(
        receipt_hash_path.read_text(encoding="utf-8") == f"{receipt_hash}\n",
        "compact receipt sidecar changed",
    )
    _require(
        receipt["artifact_manifest_sha256"] == sha256_file(manifest_path),
        "compact receipt artifact binding changed",
    )
    if contract is not None:
        bindings = artifact["bindings"]
        _require(
            artifact["source"]
            == {
                "path": str(contract.source_h5ad.relative_to(ROOT)),
                "size_bytes": contract.expected_source_size,
                "sha256": contract.expected_source_sha256,
            },
            "artifact source binding changed",
        )
        _require(
            bindings["source_manifest_sha256"]
            == contract.expected_source_manifest_sha256,
            "artifact source manifest binding changed",
        )
        _require(
            bindings["source_audit_sha256"] == sha256_file(contract.source_audit),
            "artifact source audit binding changed",
        )
        _require(
            bindings["source_acquisition_receipt_sha256"]
            == sha256_file(contract.acquisition_receipt),
            "artifact acquisition receipt binding changed",
        )
        _require(
            bindings["source_audit_implementation_sha256"]
            == sha256_file(HERE / "source_audit.py"),
            "artifact source audit implementation binding changed",
        )
        _require(
            bindings["extractor_sha256"] == sha256_file(Path(__file__).resolve()),
            "artifact extractor binding changed",
        )
        _require(
            bindings["canonical_genes_file_sha256"]
            == contract.expected_canonical_genes_sha256,
            "artifact canonical gene binding changed",
        )
        _require(
            bindings["vcc_panel_file_sha256"] == contract.expected_panel_sha256,
            "artifact panel binding changed",
        )
        _require(
            receipt["source_manifest_sha256"]
            == contract.expected_source_manifest_sha256,
            "receipt source manifest binding changed",
        )
        _require(
            receipt["source_audit_sha256"] == sha256_file(contract.source_audit),
            "receipt source audit binding changed",
        )
        _require(
            receipt["source_acquisition_receipt_sha256"]
            == sha256_file(contract.acquisition_receipt),
            "receipt acquisition binding changed",
        )
        _require(
            receipt["source_sha256"] == contract.expected_source_sha256,
            "receipt source binding changed",
        )

    meta = load_json(directory / "meta.json")
    rows = meta["axes"]["perturbations"]
    genes = meta["axes"]["canonical_genes"]
    if contract is not None:
        _require(
            rows == contract.expected_target_rows
            and genes == contract.expected_canonical_genes,
            "compact axes changed",
        )
    _require(
        meta["mapping"]["method"] == "exact_symbol_only"
        and meta["mapping"]["alias_merge"] is False,
        "mapping contract changed",
    )
    _require(
        all(value is False for value in meta["safety"].values()),
        "source-only safety declaration changed",
    )

    lfc = np.load(
        directory / "log2_fold_change.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    fdr = np.load(
        directory / "neg_log10_fdr.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    measured = np.load(
        directory / "measured_mask.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    fdr_mask = np.load(
        directory / "fdr_mask.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    expected_shape = (rows, genes)
    _require(lfc.shape == expected_shape and lfc.dtype == np.float32, "LFC changed")
    _require(fdr.shape == expected_shape and fdr.dtype == np.float32, "FDR changed")
    _require(
        measured.shape == expected_shape and measured.dtype == np.bool_,
        "measured mask changed",
    )
    _require(
        fdr_mask.shape == expected_shape and fdr_mask.dtype == np.bool_,
        "FDR mask changed",
    )
    _require(np.all(np.isnan(lfc[~measured])), "unsupported LFC is not NaN")
    _require(np.all(np.isfinite(lfc[measured])), "measured LFC is not finite")
    _require(np.all(np.isnan(fdr[~fdr_mask])), "unsupported FDR is not NaN")
    _require(np.all(np.isfinite(fdr[fdr_mask])), "available FDR is not finite")
    _require(not np.any(fdr_mask & ~measured), "FDR support exceeds LFC support")

    read_text_axis(directory / "perts.txt", rows)
    read_text_axis(directory / "genes.txt", genes)
    _require(
        sha256_file(directory / "perts.txt")
        == artifact["bindings"]["compact_perturbations_sha256"],
        "compact perturbation axis changed",
    )
    _require(
        sha256_file(directory / "genes.txt")
        == artifact["bindings"]["compact_canonical_genes_sha256"],
        "compact canonical gene axis changed",
    )
    return meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("build", "verify"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument("--chunk-rows", type=int, default=32)
    args = parser.parse_args()

    contract = production_contract()
    if args.command == "build":
        output = extract_compact_source(
            contract,
            args.output_dir,
            chunk_rows=args.chunk_rows,
        )
        meta = verify_compact_source(output, contract)
    else:
        output = safe_output(args.output_dir)
        meta = verify_compact_source(output, contract)
    print(
        json.dumps(
            {
                "status": "verified_compact_source",
                "output": str(output),
                "axes": meta["axes"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
