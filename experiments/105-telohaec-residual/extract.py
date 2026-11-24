#!/usr/bin/env python3
"""Memory-bounded extraction of guide-robust TeloHAEC perturbation effects."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENT_ROOT = (REPO_ROOT / "experiments/105-telohaec-residual").resolve()
DATA_ROOT = (REPO_ROOT / "data/external/igvf-telohaec").resolve()
ALLOWED_WRITE_ROOTS = (EXPERIMENT_ROOT, DATA_ROOT)
CONTROL_TYPES = frozenset({"safe-targeting", "non-targeting"})
SCALE = 50_000.0


@dataclass(frozen=True)
class Guide:
    guide_id: str
    spacer: str
    guide_type: str
    description: str
    clean_target: str | None
    exclusion_reason: str


def _decoded(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def ensure_allowed_output(path: Path) -> Path:
    """Resolve an output path and reject writes outside the two authorized roots."""
    resolved = path.expanduser().resolve()
    if not any(resolved == root or resolved.is_relative_to(root) for root in ALLOWED_WRITE_ROOTS):
        raise ValueError(f"output path is outside the authorized write roots: {resolved}")
    return resolved


def clean_target(description: str, guide_id: str = "") -> tuple[str | None, str]:
    """Map a targeting-guide description to a conservative gene symbol."""
    description = description.strip()
    low_description = description.lower()
    low_guide_id = guide_id.lower()
    if (
        low_guide_id.startswith("enhancer_")
        or low_guide_id.startswith("9p21_")
        or low_description.startswith("enhancer_")
        or low_description.startswith("enhancer-")
        or low_description.startswith("9p21_")
        or low_description.startswith("9p21-")
    ):
        return None, "regulatory_element"
    if "_and_" in low_description or "-and-" in low_description:
        return None, "compound_target"
    symbol = re.sub(r"[_-]TSS\d+$", "", description, flags=re.IGNORECASE).upper()
    if not symbol or not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]*", symbol):
        return None, "non_gene_label"
    return symbol, ""


def load_guides(path: Path) -> list[Guide]:
    guides: list[Guide] = []
    seen_spacers: set[str] = set()
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"guide_id", "spacer", "type", "description"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"guide manifest lacks required columns: {sorted(required)}")
        for row in reader:
            spacer = row["spacer"].strip().upper()
            if not spacer or spacer in seen_spacers:
                raise ValueError(f"empty or duplicate spacer in guide manifest: {spacer!r}")
            seen_spacers.add(spacer)
            guide_type = row["type"].strip().lower()
            description = row["description"].strip()
            if guide_type == "targeting":
                target, reason = clean_target(description, row["guide_id"])
            elif guide_type in CONTROL_TYPES:
                target, reason = None, ""
            else:
                target, reason = None, "unsupported_guide_type"
            guides.append(
                Guide(
                    guide_id=row["guide_id"].strip(),
                    spacer=spacer,
                    guide_type=guide_type,
                    description=description,
                    clean_target=target,
                    exclusion_reason=reason,
                )
            )
    return guides


def load_panel(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "target_gene" not in reader.fieldnames:
            raise ValueError("panel CSV must contain a target_gene column")
        panel = [row["target_gene"].strip().upper() for row in reader]
    if len(panel) != 300 or len(set(panel)) != 300:
        raise ValueError(f"expected 300 unique panel targets, found {len(panel)} rows")
    return panel


def feature_axis(raw_names: Iterable[str], panel: set[str]) -> tuple[list[str], np.ndarray, list[str]]:
    """Collapse duplicate symbols and remove panel coordinates from the working axis."""
    unique_genes: list[str] = []
    unique_index: dict[str, int] = {}
    raw_to_unique: list[int] = []
    for raw_name in raw_names:
        if ":" not in raw_name:
            raise ValueError(f"feature name lacks SYMBOL:ENSEMBL structure: {raw_name!r}")
        symbol, ensembl = raw_name.split(":", 1)
        symbol = symbol.strip().upper()
        if not symbol or not ensembl.startswith("ENSG"):
            raise ValueError(f"unexpected feature name: {raw_name!r}")
        if symbol not in unique_index:
            unique_index[symbol] = len(unique_genes)
            unique_genes.append(symbol)
        raw_to_unique.append(unique_index[symbol])

    model_genes = [gene for gene in unique_genes if gene not in panel]
    model_index = {gene: i for i, gene in enumerate(model_genes)}
    raw_to_model = np.asarray(
        [model_index.get(unique_genes[i], -1) for i in raw_to_unique],
        dtype=np.int32,
    )
    return unique_genes, raw_to_model, model_genes


def log_pseudobulk(counts: np.ndarray) -> np.ndarray:
    total = float(np.sum(counts, dtype=np.float64))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("guide pseudobulk has nonpositive or nonfinite non-panel total")
    return np.log1p(SCALE * np.asarray(counts, dtype=np.float64) / total)


def _atomic_text(path: Path, text: str) -> None:
    path = ensure_allowed_output(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    ensure_allowed_output(temporary).write_text(text)
    os.replace(temporary, path)


def _sha256(path: Path, block_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def _scan_cells(
    barcode_dataset: h5py.Dataset,
    guides: list[Guide],
    chunk_rows: int,
) -> tuple[np.ndarray, np.ndarray]:
    spacer_to_index = {guide.spacer: i for i, guide in enumerate(guides)}
    cell_guides = np.empty(barcode_dataset.shape[0], dtype=np.int32)
    for start in range(0, barcode_dataset.shape[0], chunk_rows):
        end = min(start + chunk_rows, barcode_dataset.shape[0])
        for offset, raw_barcode in enumerate(barcode_dataset[start:end]):
            barcode = _decoded(raw_barcode)
            parts = barcode.split(":", 2)
            if len(parts) != 3:
                raise ValueError(f"unexpected barcode structure at row {start + offset}: {barcode!r}")
            spacer = parts[1].upper()
            try:
                cell_guides[start + offset] = spacer_to_index[spacer]
            except KeyError as exc:
                raise ValueError(f"cell spacer absent from guide manifest: {spacer}") from exc
    counts = np.bincount(cell_guides, minlength=len(guides)).astype(np.int64, copy=False)
    return cell_guides, counts


def _eligibility(
    guides: list[Guide],
    guide_cell_counts: np.ndarray,
    min_guide_cells: int,
    min_target_cells: int,
) -> tuple[np.ndarray, list[int], dict[str, int], dict[str, int]]:
    target_total_cells: dict[str, int] = {}
    target_eligible_cells: dict[str, int] = {}
    for guide, count in zip(guides, guide_cell_counts, strict=True):
        if guide.clean_target is not None:
            target_total_cells[guide.clean_target] = target_total_cells.get(guide.clean_target, 0) + int(count)
            if count >= min_guide_cells:
                target_eligible_cells[guide.clean_target] = (
                    target_eligible_cells.get(guide.clean_target, 0) + int(count)
                )

    eligible_manifest_indices: list[int] = []
    for index, (guide, count) in enumerate(zip(guides, guide_cell_counts, strict=True)):
        enough_guide_cells = count >= min_guide_cells
        if guide.guide_type in CONTROL_TYPES:
            eligible = enough_guide_cells
        else:
            eligible = (
                guide.clean_target is not None
                and enough_guide_cells
                and target_eligible_cells.get(guide.clean_target, 0) >= min_target_cells
            )
        if eligible:
            eligible_manifest_indices.append(index)

    guide_to_aggregate = np.full(len(guides), -1, dtype=np.int32)
    guide_to_aggregate[eligible_manifest_indices] = np.arange(
        len(eligible_manifest_indices), dtype=np.int32
    )
    return (
        guide_to_aggregate,
        eligible_manifest_indices,
        target_total_cells,
        target_eligible_cells,
    )


def aggregate_csr_to_guides(
    x_group: h5py.Group,
    cell_aggregate_ids: np.ndarray,
    raw_to_model: np.ndarray,
    output_path: Path,
    n_guides: int,
    n_model_genes: int,
    row_block: int,
) -> tuple[np.memmap, dict[str, int | float]]:
    """Scan a CSR matrix once and aggregate eligible rows into a disk memmap."""
    output_path = ensure_allowed_output(output_path)
    shape = tuple(int(x) for x in x_group.attrs["shape"])
    if x_group.attrs.get("encoding-type") != "csr_matrix":
        raise ValueError("X must be encoded as CSR")
    if shape[0] != len(cell_aggregate_ids) or shape[1] != len(raw_to_model):
        raise ValueError("CSR shape is inconsistent with obs or var")

    indptr = np.asarray(x_group["indptr"][:], dtype=np.int64)
    if len(indptr) != shape[0] + 1 or indptr[0] != 0 or np.any(np.diff(indptr) < 0):
        raise ValueError("invalid CSR indptr")
    if indptr[-1] != len(x_group["data"]) or indptr[-1] != len(x_group["indices"]):
        raise ValueError("CSR arrays have inconsistent lengths")

    sums = np.memmap(
        output_path,
        mode="w+",
        dtype=np.float64,
        shape=(n_guides, n_model_genes),
    )
    sums[:] = 0.0
    flat_sums = sums.reshape(-1)
    explicit_zeros = 0
    max_value = 0.0
    checked_nnz = 0

    for start_row in range(0, shape[0], row_block):
        end_row = min(start_row + row_block, shape[0])
        start_nnz = int(indptr[start_row])
        end_nnz = int(indptr[end_row])
        data = np.asarray(x_group["data"][start_nnz:end_nnz], dtype=np.float64)
        indices = np.asarray(x_group["indices"][start_nnz:end_nnz], dtype=np.int64)
        if not np.all(np.isfinite(data)) or np.any(data < 0):
            raise ValueError(f"X contains nonfinite or negative values near row {start_row}")
        if np.any(data != np.rint(data)):
            raise ValueError(f"X is not count-like near row {start_row}")
        if np.any(indices < 0) or np.any(indices >= shape[1]):
            raise ValueError(f"X contains out-of-range column indices near row {start_row}")
        explicit_zeros += int(np.count_nonzero(data == 0))
        if data.size:
            max_value = max(max_value, float(data.max()))
        checked_nnz += int(data.size)

        row_lengths = np.diff(indptr[start_row : end_row + 1])
        aggregate_rows = np.repeat(cell_aggregate_ids[start_row:end_row], row_lengths)
        model_columns = raw_to_model[indices]
        keep = (aggregate_rows >= 0) & (model_columns >= 0)
        if np.any(keep):
            flat_indices = (
                aggregate_rows[keep].astype(np.int64, copy=False) * n_model_genes
                + model_columns[keep].astype(np.int64, copy=False)
            )
            np.add.at(flat_sums, flat_indices, data[keep])

    sums.flush()
    return sums, {
        "nnz_checked": checked_nnz,
        "explicit_zero_entries": explicit_zeros,
        "maximum_count_entry": max_value,
    }


def _write_guide_qc(
    path: Path,
    guides: list[Guide],
    guide_cell_counts: np.ndarray,
    guide_to_aggregate: np.ndarray,
    target_total_cells: dict[str, int],
    target_eligible_cells: dict[str, int],
    min_guide_cells: int,
    min_target_cells: int,
) -> None:
    path = ensure_allowed_output(path)
    temporary = ensure_allowed_output(path.with_name(f".{path.name}.tmp"))
    with temporary.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "guide_id",
                "spacer",
                "guide_type",
                "description",
                "clean_target",
                "exclusion_reason",
                "observed_cells",
                "target_total_cells",
                "target_eligible_cells",
                "eligible",
                "eligibility_note",
            ]
        )
        for index, (guide, count) in enumerate(zip(guides, guide_cell_counts, strict=True)):
            target_total = target_total_cells.get(guide.clean_target or "", 0)
            target_eligible = target_eligible_cells.get(guide.clean_target or "", 0)
            eligible = guide_to_aggregate[index] >= 0
            if count == 0:
                note = "unobserved"
            elif count < min_guide_cells:
                note = f"fewer_than_{min_guide_cells}_guide_cells"
            elif guide.guide_type in CONTROL_TYPES:
                note = "eligible_control"
            elif guide.clean_target is None:
                note = guide.exclusion_reason
            elif target_eligible < min_target_cells:
                note = f"fewer_than_{min_target_cells}_target_cells"
            else:
                note = "eligible_targeting"
            writer.writerow(
                [
                    guide.guide_id,
                    guide.spacer,
                    guide.guide_type,
                    guide.description,
                    guide.clean_target or "",
                    guide.exclusion_reason,
                    int(count),
                    target_total,
                    target_eligible,
                    int(eligible),
                    note,
                ]
            )
    os.replace(temporary, path)


def extract(
    h5ad_path: Path,
    guide_path: Path,
    panel_path: Path,
    output_dir: Path,
    *,
    min_guide_cells: int = 5,
    min_target_cells: int = 20,
    row_block: int = 512,
    hash_outputs: bool = True,
) -> dict[str, object]:
    output_dir = ensure_allowed_output(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    guides = load_guides(guide_path)
    panel_list = load_panel(panel_path)
    panel = set(panel_list)

    temporary_sums = ensure_allowed_output(output_dir / ".guide_sums.float64.memmap")
    delta_temp = ensure_allowed_output(output_dir / ".delta.npy.tmp")
    delta_path = ensure_allowed_output(output_dir / "delta.npy")
    try:
        with h5py.File(h5ad_path, "r") as handle:
            if not {"X", "obs", "var"}.issubset(handle):
                raise ValueError("H5AD lacks X, obs, or var")
            raw_features = [_decoded(x) for x in handle["var/_index"][:]]
            unique_genes, raw_to_model, model_genes = feature_axis(raw_features, panel)
            if handle["obs/barcodes"].shape[0] != int(handle["X"].attrs["shape"][0]):
                raise ValueError("barcode row count does not match X")
            cell_guides, guide_cell_counts = _scan_cells(
                handle["obs/barcodes"], guides, chunk_rows=max(row_block, 1024)
            )
            (
                guide_to_aggregate,
                eligible_manifest_indices,
                target_total_cells,
                target_eligible_cells,
            ) = _eligibility(
                guides,
                guide_cell_counts,
                min_guide_cells=min_guide_cells,
                min_target_cells=min_target_cells,
            )
            cell_aggregate_ids = guide_to_aggregate[cell_guides]
            guide_sums, matrix_qc = aggregate_csr_to_guides(
                handle["X"],
                cell_aggregate_ids,
                raw_to_model,
                temporary_sums,
                n_guides=len(eligible_manifest_indices),
                n_model_genes=len(model_genes),
                row_block=row_block,
            )

        eligible_guides = [guides[i] for i in eligible_manifest_indices]
        control_rows = [
            i for i, guide in enumerate(eligible_guides) if guide.guide_type in CONTROL_TYPES
        ]
        if not control_rows:
            raise ValueError("no eligible control guides")
        control_profiles = np.stack([log_pseudobulk(guide_sums[i]) for i in control_rows])
        control_model = np.median(control_profiles, axis=0)

        target_rows: dict[str, list[int]] = {}
        for row, guide in enumerate(eligible_guides):
            if guide.clean_target is not None:
                target_rows.setdefault(guide.clean_target, []).append(row)
        perts = sorted(target_rows)
        unique_index = {gene: i for i, gene in enumerate(unique_genes)}
        model_unique_indices = np.asarray(
            [unique_index[gene] for gene in model_genes], dtype=np.int64
        )

        delta = np.lib.format.open_memmap(
            delta_temp,
            mode="w+",
            dtype=np.float32,
            shape=(len(perts), len(unique_genes)),
        )
        delta[:] = 0.0
        b_ctrl = np.zeros(len(unique_genes), dtype=np.float32)
        b_ctrl[model_unique_indices] = control_model.astype(np.float32)
        for perturbation_index, target in enumerate(perts):
            profiles = np.stack(
                [log_pseudobulk(guide_sums[row]) for row in target_rows[target]]
            )
            effect = np.median(profiles, axis=0) - control_model
            delta[perturbation_index, model_unique_indices] = effect.astype(np.float32)
        delta.flush()
        os.replace(delta_temp, delta_path)
        np.save(ensure_allowed_output(output_dir / "b_ctrl.npy"), b_ctrl)
        _atomic_text(output_dir / "genes.txt", "".join(f"{gene}\n" for gene in unique_genes))
        _atomic_text(output_dir / "perts.txt", "".join(f"{target}\n" for target in perts))
        _write_guide_qc(
            output_dir / "guide_qc.csv",
            guides,
            guide_cell_counts,
            guide_to_aggregate,
            target_total_cells,
            target_eligible_cells,
            min_guide_cells,
            min_target_cells,
        )

        panel_axis_overlap = sorted(panel.intersection(unique_genes))
        retained_panel_targets = sorted(panel.intersection(perts))
        panel_indices = np.asarray(
            [unique_index[gene] for gene in panel_axis_overlap], dtype=np.int64
        )
        saved_delta = np.load(delta_path, mmap_mode="r")
        if panel_indices.size and np.any(saved_delta[:, panel_indices] != 0):
            raise AssertionError("panel coordinates are not exactly zero")
        if panel_indices.size and np.any(b_ctrl[panel_indices] != 0):
            raise AssertionError("control panel coordinates are not exactly zero")

        output_names = [
            "delta.npy",
            "b_ctrl.npy",
            "genes.txt",
            "perts.txt",
            "guide_qc.csv",
        ]
        output_hashes = (
            {name: _sha256(output_dir / name) for name in output_names}
            if hash_outputs
            else {}
        )
        report: dict[str, object] = {
            "source": {
                "h5ad": h5ad_path.name,
                "guide_manifest": guide_path.name,
            },
            "parameters": {
                "minimum_cells_per_guide": min_guide_cells,
                "minimum_cells_per_target": min_target_cells,
                "target_cell_count_definition": (
                    "sum of cells from guides meeting minimum_cells_per_guide"
                ),
                "row_block": row_block,
                "pseudobulk_scale": SCALE,
                "control_aggregation": "coordinate-wise median of equal-weight guide log-pseudobulks",
                "target_aggregation": "coordinate-wise median of guide effects relative to control",
            },
            "schema": {
                "cells": int(len(cell_guides)),
                "raw_features": len(raw_features),
                "unique_feature_symbols": len(unique_genes),
                "duplicate_feature_rows_collapsed": len(raw_features) - len(unique_genes),
                "nonpanel_working_features": len(model_genes),
            },
            "guides": {
                "manifest_total": len(guides),
                "observed_total": int(np.count_nonzero(guide_cell_counts)),
                "eligible_total": len(eligible_guides),
                "eligible_control": len(control_rows),
                "eligible_targeting": len(eligible_guides) - len(control_rows),
            },
            "targets": {
                "clean_manifest_targets": len(
                    {guide.clean_target for guide in guides if guide.clean_target is not None}
                ),
                "retained_targets": len(perts),
                "retained_target_eligible_cells_min": min(
                    target_eligible_cells[target] for target in perts
                ),
                "targets_excluded_after_guide_filter": sum(
                    0 < cells < min_target_cells
                    for cells in target_eligible_cells.values()
                ),
                "panel_targets_total": len(panel),
                "panel_targets_in_clean_manifest": len(
                    panel.intersection(
                        {guide.clean_target for guide in guides if guide.clean_target is not None}
                    )
                ),
                "panel_targets_retained": len(retained_panel_targets),
                "panel_targets_retained_symbols": retained_panel_targets,
            },
            "panel_coordinate_exclusion": {
                "panel_targets_total": len(panel),
                "panel_symbols_present_on_source_axis": len(panel_axis_overlap),
                "panel_symbols_absent_from_source_axis": sorted(panel.difference(unique_genes)),
                "excluded_before_library_size_normalization": True,
                "delta_exact_zero_verified": True,
                "control_exact_zero_verified": True,
            },
            "matrix_validation": matrix_qc,
            "output": {
                "delta_shape": [len(perts), len(unique_genes)],
                "delta_dtype": "float32",
                "control_shape": [len(unique_genes)],
                "control_dtype": "float32",
                "sha256": output_hashes,
            },
        }
        _atomic_text(
            output_dir / "extraction_report.json",
            json.dumps(report, indent=2, sort_keys=True) + "\n",
        )
        return report
    finally:
        if temporary_sums.exists():
            temporary_sums.unlink()
        if delta_temp.exists():
            delta_temp.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=DATA_ROOT / "IGVFFI8645YMGX.h5ad",
    )
    parser.add_argument(
        "--guides",
        type=Path,
        default=DATA_ROOT / "IGVFFI9754AGFB.tsv.gz",
    )
    parser.add_argument(
        "--panel",
        type=Path,
        default=REPO_ROOT / "data/challenge/pert_counts.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=DATA_ROOT)
    parser.add_argument("--min-guide-cells", type=int, default=5)
    parser.add_argument("--min-target-cells", type=int, default=20)
    parser.add_argument("--row-block", type=int, default=512)
    args = parser.parse_args()
    report = extract(
        args.h5ad,
        args.guides,
        args.panel,
        args.output_dir,
        min_guide_cells=args.min_guide_cells,
        min_target_cells=args.min_target_cells,
        row_block=args.row_block,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
