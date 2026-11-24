#!/usr/bin/env python3
"""Stream K562 GWPS into guide-pair and target pseudobulk effect tables.

Guide effects retain only guide categories meeting the minimum-cell threshold.
Compact-target reconstruction is a separate integrity estimand: it pools every
cell from targets whose total target-level cell count meets that threshold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import shutil
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import h5py
import numpy as np


TARGET_SUM = 50_000.0
CONTROL_LABEL = "control"


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_memory(limit_gib: float, stage: str) -> None:
    peak = peak_rss_gib()
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def decode(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode() if isinstance(value, bytes) else str(value) for value in values
    )


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_panel(path: Path) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "target_gene" not in reader.fieldnames:
            raise ValueError("panel CSV lacks target_gene")
        return tuple(str(row["target_gene"]).strip() for row in reader)


def assign_guide_targets(
    gene_codes_ds: h5py.Dataset,
    guide_codes_ds: h5py.Dataset,
    n_guides: int,
    *,
    block_rows: int = 250_000,
) -> tuple[np.ndarray, np.ndarray]:
    guide_target = np.full(n_guides, -1, dtype=np.int32)
    guide_cells = np.zeros(n_guides, dtype=np.int64)
    for lo in range(0, len(gene_codes_ds), block_rows):
        hi = min(lo + block_rows, len(gene_codes_ds))
        genes = np.asarray(gene_codes_ds[lo:hi], dtype=np.int32)
        guides = np.asarray(guide_codes_ds[lo:hi], dtype=np.int32)
        valid = (genes >= 0) & (guides >= 0)
        guide_cells += np.bincount(guides[valid], minlength=n_guides)
        pairs = np.unique(np.stack((guides[valid], genes[valid]), axis=1), axis=0)
        for guide, target in pairs:
            old = int(guide_target[guide])
            if old >= 0 and old != int(target):
                raise ValueError(f"guide category {guide} maps to multiple targets")
            guide_target[guide] = int(target)
    return guide_target, guide_cells


def category_counts(
    codes_ds: h5py.Dataset,
    n_categories: int,
    *,
    block_rows: int = 250_000,
) -> np.ndarray:
    """Count categorical rows without opening the expression matrix."""
    counts = np.zeros(n_categories, dtype=np.int64)
    for lo in range(0, len(codes_ds), block_rows):
        hi = min(lo + block_rows, len(codes_ds))
        codes = np.asarray(codes_ds[lo:hi], dtype=np.int32)
        if np.any(codes < 0) or np.any(codes >= n_categories):
            raise ValueError("target category codes are missing or out of range")
        counts += np.bincount(codes, minlength=n_categories)
    return counts


def compact_target_layout(
    targets: tuple[str, ...],
    target_cells: np.ndarray,
    compact_perts: tuple[str, ...],
    *,
    control_code: int,
    min_cells: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map every compact target to all of its source cells.

    Guide eligibility is deliberately absent from this function. The compact
    donor was built by filtering total cells per perturbation and then pooling
    every cell assigned to each retained perturbation.
    """
    counts = np.asarray(target_cells, dtype=np.int64)
    if counts.shape != (len(targets),):
        raise ValueError("target-cell counts do not match target categories")
    if len(targets) != len(set(targets)):
        raise ValueError("H5AD target categories are not unique")
    expected = tuple(
        sorted(
            target
            for index, target in enumerate(targets)
            if index != control_code and counts[index] >= min_cells
        )
    )
    if compact_perts != expected:
        raise ValueError(
            "compact K562 target axis does not match total-cell eligibility"
        )
    source_code = {target: index for index, target in enumerate(targets)}
    compact_codes = np.asarray(
        [source_code[target] for target in compact_perts],
        dtype=np.int32,
    )
    output_row = np.full(len(targets), -1, dtype=np.int32)
    output_row[compact_codes] = np.arange(len(compact_codes), dtype=np.int32)
    return compact_codes, output_row, counts[compact_codes]


def grouped_sum_block(
    block: np.ndarray,
    codes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    valid = codes >= 0
    if not np.any(valid):
        return np.empty(0, dtype=np.int64), np.empty((0, block.shape[1]))
    rows = np.flatnonzero(valid)
    order = np.argsort(codes[rows], kind="stable")
    sorted_rows = rows[order]
    sorted_codes = codes[sorted_rows]
    starts = np.r_[0, 1 + np.flatnonzero(sorted_codes[1:] != sorted_codes[:-1])]
    groups = sorted_codes[starts]
    sums = np.add.reduceat(
        np.asarray(block[sorted_rows], dtype=np.float64), starts, axis=0
    )
    return groups.astype(np.int64, copy=False), sums


def extract(args: argparse.Namespace) -> dict[str, object]:
    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    scratch = args.output / ".scratch"
    scratch.mkdir()
    panel = read_panel(args.panel_csv)
    compact_perts = read_lines(args.compact_k562 / "perts.txt")
    compact_genes = read_lines(args.compact_k562 / "genes.txt")
    compact_delta = np.load(
        args.compact_k562 / "delta.npy", mmap_mode="r", allow_pickle=False
    )
    if compact_delta.shape != (len(compact_perts), len(compact_genes)):
        raise ValueError("compact K562 shape mismatch")

    try:
        with h5py.File(args.h5ad, "r") as handle:
            x = handle["X"]
            if not isinstance(x, h5py.Dataset):
                raise ValueError("expected a dense HDF5 X dataset")
            genes = decode(handle["var/gene_name"][:])
            targets = decode(handle["obs/perturbation/categories"][:])
            guides = decode(handle["obs/guide_id/categories"][:])
            if tuple(genes) != compact_genes:
                raise ValueError("H5AD and compact K562 gene axes differ")
            target_codes_ds = handle["obs/perturbation/codes"]
            guide_codes_ds = handle["obs/guide_id/codes"]
            guide_target, guide_cells = assign_guide_targets(
                target_codes_ds, guide_codes_ds, len(guides)
            )
            control_code = targets.index(CONTROL_LABEL)
            all_target_cells = category_counts(target_codes_ds, len(targets))
            target_codes, target_output_rows, target_cells = compact_target_layout(
                targets,
                all_target_cells,
                compact_perts,
                control_code=control_code,
                min_cells=args.min_cells,
            )
            guide_sums_path = scratch / "guide_sums.npy"
            guide_sums = np.lib.format.open_memmap(
                guide_sums_path,
                mode="w+",
                dtype=np.float64,
                shape=(len(guides), len(genes)),
            )
            guide_sums[:] = 0
            target_sums_path = scratch / "target_sums.npy"
            target_sums = np.lib.format.open_memmap(
                target_sums_path,
                mode="w+",
                dtype=np.float64,
                shape=(len(target_codes), len(genes)),
            )
            target_sums[:] = 0
            control_sum = np.zeros(len(genes), dtype=np.float64)
            control_cells = 0
            for lo in range(0, x.shape[0], args.block_rows):
                hi = min(lo + args.block_rows, x.shape[0])
                block = np.asarray(x[lo:hi], dtype=np.float32)
                target_codes = np.asarray(target_codes_ds[lo:hi], dtype=np.int32)
                qcodes = np.asarray(guide_codes_ds[lo:hi], dtype=np.int32)
                controls = target_codes == control_code
                if np.any(controls):
                    control_sum += np.asarray(block[controls], dtype=np.float64).sum(
                        axis=0
                    )
                    control_cells += int(np.sum(controls))
                perturbed = ~controls
                groups, sums = grouped_sum_block(block[perturbed], qcodes[perturbed])
                if len(groups):
                    guide_sums[groups] += sums
                groups, sums = grouped_sum_block(
                    block,
                    target_output_rows[target_codes],
                )
                if len(groups):
                    target_sums[groups] += sums
                del block, groups, sums
                if lo % (args.block_rows * 20) == 0:
                    print(f"streamed {hi:,}/{x.shape[0]:,} cells", flush=True)
                enforce_memory(args.max_rss_gib, f"extract rows {lo}:{hi}")
            guide_sums.flush()
            target_sums.flush()

        valid_guides = np.flatnonzero(
            (guide_target >= 0)
            & (guide_target != control_code)
            & (guide_cells >= args.min_cells)
        )
        output_guides = tuple(guides[index] for index in valid_guides)
        output_targets = tuple(
            targets[int(guide_target[index])] for index in valid_guides
        )
        output_cells = guide_cells[valid_guides]
        b_ctrl = np.log1p(TARGET_SUM * control_sum / control_sum.sum())
        guide_delta_path = args.output / "guide_delta.npy"
        guide_delta = np.lib.format.open_memmap(
            guide_delta_path,
            mode="w+",
            dtype=np.float32,
            shape=(len(valid_guides), len(genes)),
        )
        for lo in range(0, len(valid_guides), args.effect_block_rows):
            hi = min(lo + args.effect_block_rows, len(valid_guides))
            sums = np.asarray(guide_sums[valid_guides[lo:hi]], dtype=np.float64)
            totals = sums.sum(axis=1, keepdims=True)
            guide_delta[lo:hi] = (
                np.log1p(TARGET_SUM * sums / np.maximum(totals, 1.0)) - b_ctrl
            ).astype(np.float32)
            enforce_memory(args.max_rss_gib, f"guide effects {lo}:{hi}")
        guide_delta.flush()

        grouped: dict[str, list[int]] = defaultdict(list)
        for row, target in enumerate(output_targets):
            grouped[target].append(row)
        target_names = compact_perts
        target_delta_path = args.output / "target_pooled_delta.npy"
        target_delta = np.lib.format.open_memmap(
            target_delta_path,
            mode="w+",
            dtype=np.float32,
            shape=(len(target_names), len(genes)),
        )
        for lo in range(0, len(target_names), args.effect_block_rows):
            hi = min(lo + args.effect_block_rows, len(target_names))
            sums = np.asarray(target_sums[lo:hi], dtype=np.float64)
            totals = sums.sum(axis=1, keepdims=True)
            if np.any(totals <= 0):
                raise ValueError("compact-eligible target has zero expression total")
            target_delta[lo:hi] = (
                np.log1p(TARGET_SUM * sums / totals) - b_ctrl
            ).astype(np.float32)
            enforce_memory(args.max_rss_gib, f"target effects {lo}:{hi}")
        target_delta.flush()

        errors = []
        for row in range(len(target_names)):
            errors.append(
                float(
                    np.max(
                        np.abs(
                            np.asarray(target_delta[row], dtype=np.float64)
                            - np.asarray(compact_delta[row], dtype=np.float64)
                        )
                    )
                )
            )
        max_error = max(errors, default=float("nan"))
        if not np.isfinite(max_error) or max_error > args.max_reconstruction_error:
            raise ValueError(
                f"target pseudobulk does not reproduce compact K562: {max_error}"
            )

        (args.output / "genes.txt").write_text("\n".join(genes) + "\n")
        (args.output / "guide_ids.txt").write_text("\n".join(output_guides) + "\n")
        (args.output / "guide_targets.txt").write_text("\n".join(output_targets) + "\n")
        (args.output / "target_names.txt").write_text("\n".join(target_names) + "\n")
        np.save(args.output / "guide_cells.npy", output_cells, allow_pickle=False)
        np.save(args.output / "target_cells.npy", target_cells, allow_pickle=False)
        multi = {target: rows for target, rows in grouped.items() if len(rows) >= 2}
        metadata: dict[str, object] = {
            "experiment": "083-guide-latent-decontamination",
            "source_h5ad": str(args.h5ad.resolve()),
            "source_h5ad_bytes": args.h5ad.stat().st_size,
            "source_h5ad_sha256": sha256_file(args.h5ad),
            "compact_k562": str(args.compact_k562.resolve()),
            "target_sum": TARGET_SUM,
            "minimum_cells_per_guide": args.min_cells,
            "minimum_cells_per_target": args.min_cells,
            "guide_effect_pooling": "only guides meeting minimum_cells_per_guide",
            "target_reconstruction_pooling": (
                "all cells for targets meeting minimum_cells_per_target"
            ),
            "cells": int(x.shape[0]),
            "genes": len(genes),
            "control_cells": control_cells,
            "guide_categories": len(guides),
            "retained_guide_pairs": len(valid_guides),
            "retained_targets": len(target_names),
            "multi_guide_targets": len(multi),
            "panel_covered_targets": sum(target in grouped for target in panel),
            "panel_multi_guide_targets": sum(target in multi for target in panel),
            "compact_reconstruction": {
                "compared_targets": len(errors),
                "maximum_absolute_error": max_error,
                "tolerance": args.max_reconstruction_error,
                "passes": max_error <= args.max_reconstruction_error,
            },
            "outputs": {
                name: sha256_file(args.output / name)
                for name in (
                    "guide_delta.npy",
                    "target_pooled_delta.npy",
                    "guide_cells.npy",
                    "target_cells.npy",
                    "genes.txt",
                    "guide_ids.txt",
                    "guide_targets.txt",
                    "target_names.txt",
                )
            },
            "peak_rss_gib": peak_rss_gib(),
        }
        (args.output / "metadata.json").write_text(
            json.dumps(metadata, indent=2) + "\n"
        )
        return metadata
    finally:
        try:
            del guide_sums
        except UnboundLocalError:
            pass
        try:
            del target_sums
        except UnboundLocalError:
            pass
        shutil.rmtree(scratch, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--compact-k562", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-rows", type=int, default=3_886)
    parser.add_argument("--effect-block-rows", type=int, default=64)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--max-rss-gib", type=float, default=24.0)
    parser.add_argument("--max-reconstruction-error", type=float, default=2e-5)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data extraction without --guarded")
    if min(args.block_rows, args.effect_block_rows, args.min_cells) < 1:
        raise SystemExit("block sizes and min cells must be positive")
    metadata = extract(args)
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
