#!/usr/bin/env python3
"""Stage leakage-safe frozen-arm-D queries for experiment 082.

Development and final-panel paths are intentionally separate:

* development queries zero all 300 VCC targets and all 600 non-panel
  development targets before the dense whitening operator, then score only
  genes outside both sets;
* final-panel queries zero exactly the 300 VCC targets before whitening and
  retain the 600 development-target genes on the final score axis.

This prevents historical self-target coordinates from mixing into safe
development genes without deleting valid final-panel coordinates.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import resource
import sys
from pathlib import Path
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEVELOPMENT_FOLDS = ("orion_hct116", "orion_hek293t")
KOLF_FOLD = "kolf2"
FINAL_FOLDS = (*DEVELOPMENT_FOLDS, KOLF_FOLD)
OUTER_DONORS = DEVELOPMENT_FOLDS
MAIN_DONORS = ("k562_gwps", *DEVELOPMENT_FOLDS)
TEACHER_DONORS = ("rpe1", "nadig_jurkat")
CD4_DONORS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CACHE_SCHEMA = 3
EXPECTED_FINAL_SCORE_GENES = 17_806


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP062_PATH = ROOT / "experiments/062-cd4-post-whiten/run.py"
EXP064_PATH = ROOT / "experiments/064-pm-cd4-composition/run.py"
EXP062 = load_module("exp082_exp062", EXP062_PATH)
EXP064 = load_module("exp082_exp064", EXP064_PATH)


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
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


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    return hashlib.sha256(array.view(np.uint8)).hexdigest()


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: axis is empty or contains duplicates")
    return values


def last_occurrence_lookup(
    labels: Sequence[str],
) -> tuple[dict[str, int], dict[str, object]]:
    """Build a deterministic unique projection without altering matrix columns."""
    lookup: dict[str, int] = {}
    duplicates: set[str] = set()
    for index, label in enumerate(labels):
        if label in lookup:
            duplicates.add(label)
        lookup[label] = index
    duplicate_symbols = tuple(sorted(duplicates))
    return lookup, {
        "axis_rows": len(labels),
        "unique_symbols": len(lookup),
        "duplicate_symbols": len(duplicate_symbols),
        "duplicate_extra_rows": len(labels) - len(lookup),
        "axis_rows_sha256": sha256_lines(labels),
        "duplicate_symbols_sha256": sha256_lines(duplicate_symbols),
        "projection_policy": (
            "preserve the full ordered matrix axis; project each unique symbol "
            "from its deterministic last occurrence"
        ),
    }


def verify_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 {actual} != expected {expected}")


def verify_inputs(args: argparse.Namespace, manifest: Mapping[str, object]) -> None:
    for filename, expected in manifest["challenge"].items():
        path = args.genes_csv if filename == "gene_names.csv" else args.panel_csv
        verify_hash(path, str(expected), filename)
    for donor, files in manifest["donors"].items():
        for filename, expected in files.items():
            verify_hash(
                args.donor_root / donor / filename,
                str(expected),
                f"{donor}/{filename}",
            )
    for filename, expected in manifest["perturbmatch"].items():
        verify_hash(args.perturbmatch / filename, str(expected), f"pm/{filename}")
    for donor, files in manifest["cd4"].items():
        for filename, expected in files.items():
            verify_hash(
                args.cd4_root / donor / filename,
                str(expected),
                f"{donor}/{filename}",
            )
    for filename, expected in manifest["truth_cache"]["files"].items():
        if filename.endswith("_raw_outer.npy"):
            continue
        verify_hash(
            args.truth_cache / filename,
            str(expected),
            f"truth-cache/{filename}",
        )
    for held, entry in manifest["operators"].items():
        verify_hash(
            args.operator_dir / f"held_{held}.npz",
            str(entry["npz_sha256"]),
            f"operator/{held}.npz",
        )
        verify_hash(
            args.operator_dir / f"held_{held}.json",
            str(entry["json_sha256"]),
            f"operator/{held}.json",
        )
    dependencies = manifest["dependencies"]
    verify_hash(
        EXP062_PATH,
        str(dependencies["exp062_run_sha256"]),
        "experiment 062 implementation",
    )
    verify_hash(
        EXP064_PATH,
        str(dependencies["exp064_run_sha256"]),
        "experiment 064 implementation",
    )
    verify_hash(
        args.exp064_results,
        str(manifest["exp064_reference"]["results_sha256"]),
        "experiment-064/results.json",
    )


def load_exp064_reference(
    path: Path,
    manifest: Mapping[str, object],
) -> tuple[str, ...]:
    result = json.loads(path.read_text())
    reference = manifest["exp064_reference"]
    targets = tuple(map(str, result["axis"]["target_order"]))
    if len(targets) != int(reference["gallery_targets"]):
        raise ValueError("experiment-064 target count changed")
    if sha256_lines(targets) != str(reference["target_order_sha256"]):
        raise ValueError("experiment-064 target order changed")
    return targets


def load_operator(
    path: Path,
    source_banks: Sequence[EXP062.Bank],
    official: Sequence[str],
    panel: set[str],
    entry: Mapping[str, object],
) -> tuple[EXP062.Operator, np.ndarray]:
    """Load a held-fold or production operator and verify its full contract."""
    expected_columns = EXP064.operator_columns(source_banks, official, panel)
    verify_hash(path, str(entry["npz_sha256"]), f"operator/{path.name}")
    metadata_path = path.with_suffix(".json")
    verify_hash(
        metadata_path,
        str(entry["json_sha256"]),
        f"operator/{metadata_path.name}",
    )
    metadata = json.loads(metadata_path.read_text())
    if tuple(map(str, metadata.get("donors", ()))) != tuple(
        bank.name for bank in source_banks
    ):
        raise ValueError(f"{path}: operator donor order mismatch")
    if int(metadata.get("n_cols", -1)) != len(expected_columns):
        raise ValueError(f"{path}: operator-column count mismatch")
    expected_axis_hash = sha256_array(expected_columns)
    metadata_axis_hash = metadata.get("wcols_sha256")
    if metadata_axis_hash is not None and metadata_axis_hash != expected_axis_hash:
        raise ValueError(f"{path}: operator-column hash mismatch")
    pinned_axis_hash = entry.get("columns_sha256")
    if pinned_axis_hash is not None and pinned_axis_hash != expected_axis_hash:
        raise ValueError(f"{path}: pinned operator-column hash mismatch")
    if "gamma" in metadata and float(metadata["gamma"]) != EXP064.WHITEN_GAMMA:
        raise ValueError(f"{path}: whitening gamma changed")
    with np.load(path) as archive:
        evals = np.asarray(archive["evals"], dtype=np.float64)
        evecs = np.asarray(archive["evecs"], dtype=np.float32)
    if evals.shape != (len(expected_columns),):
        raise ValueError(f"{path}: eigenvalue shape mismatch")
    if evecs.shape != (len(expected_columns), len(expected_columns)):
        raise ValueError(f"{path}: eigenvector shape mismatch")
    return (
        EXP062.Operator(evals, evecs, metadata, path, metadata_path),
        expected_columns,
    )


def apply_operator_restriction(
    matrix: np.ndarray,
    operator: EXP062.Operator,
    operator_columns: np.ndarray,
    active_columns: np.ndarray,
    *,
    preserve_rows: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    """Apply the principal restriction of the dense whitening operator.

    For zero forbidden inputs, the retained dense output is exactly
    ``x_R W_RR``. Evaluating that principal block directly prevents forbidden
    historical-target outputs from ever being materialized or fed downstream.
    """
    full = np.asarray(operator_columns, dtype=np.int64)
    active = np.asarray(active_columns, dtype=np.int64)
    keep = np.isin(full, active, assume_unique=True)
    if not np.array_equal(full[keep], active):
        raise ValueError("active operator columns are not an ordered subset")
    if np.any(np.asarray(matrix[:, full[~keep]], dtype=np.float64) != 0.0):
        raise AssertionError("forbidden coordinate is nonzero before whitening")
    selected_evecs = np.asarray(operator.evecs[keep], dtype=np.float64)
    weights = np.asarray(
        operator.evals ** (-EXP064.WHITEN_GAMMA / 2.0), dtype=np.float64
    )
    source = np.asarray(matrix[:, active], dtype=np.float64)
    transformed = ((source @ selected_evecs) * weights) @ selected_evecs.T
    preserve = np.asarray(preserve_rows, dtype=np.int64)
    before = float(np.median(np.linalg.norm(source[preserve], axis=1))) or 1.0
    after = float(np.median(np.linalg.norm(transformed[preserve], axis=1))) or 1.0
    scale = before / after
    output = np.asarray(matrix, dtype=np.float64).copy()
    output[:, full] = 0.0
    output[:, active] = transformed * scale
    return output, {
        "method": "principal restriction W_RR = E_R diag(lambda^-gamma/2) E_R.T",
        "full_operator_columns": int(len(full)),
        "active_operator_columns": int(len(active)),
        "forbidden_input_coordinates": int(len(full) - len(active)),
        "active_columns_sha256": sha256_array(active),
        "dense_zero_equivalence_proved": True,
        "proof": (
            "For x_F=0, (x W)_R = x_R W_RR; the implementation evaluates "
            "that principal block directly and never materializes (x W)_F."
        ),
        "before_operator_median_norm": before,
        "after_operator_median_norm": after,
        "applied_scale": scale,
        "preserve_rows": int(len(preserve)),
    }


def add_post_whiten_residual(
    baseline: np.ndarray,
    cd4: EXP062.CD4Prepared,
    active_columns: np.ndarray,
    *,
    preserve_rows: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    candidate = np.asarray(baseline, dtype=np.float64).copy()
    active_mask = np.zeros(candidate.shape[1], dtype=bool)
    active_mask[np.asarray(active_columns, dtype=np.int64)] = True
    available = (cd4.gene_votes > 0) & active_mask[None, :]
    candidate[available] += EXP064.CD4_FRACTION * cd4.mean_effect[available]
    preserve = np.asarray(preserve_rows, dtype=np.int64)
    active = np.asarray(active_columns, dtype=np.int64)
    baseline_norm = (
        float(np.median(np.linalg.norm(baseline[np.ix_(preserve, active)], axis=1)))
        or 1.0
    )
    candidate_norm = (
        float(np.median(np.linalg.norm(candidate[np.ix_(preserve, active)], axis=1)))
        or 1.0
    )
    scale = baseline_norm / candidate_norm
    candidate *= scale
    return candidate, {
        "fraction": EXP064.CD4_FRACTION,
        "baseline_output_median_norm": baseline_norm,
        "candidate_pre_preserve_median_norm": candidate_norm,
        "post_add_global_preserve_scale": scale,
    }


def validate_truth_contract(
    truth_cache: Path,
    official: Sequence[str],
    panel: Sequence[str],
    *,
    expected_training_targets: int,
    minimum_outer_targets: int,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], np.ndarray]:
    master_genes = read_lines(truth_cache / "master_genes.txt")
    train_targets = read_lines(truth_cache / "train_perts.txt")
    outer_targets = read_lines(truth_cache / "outer_perts.txt")
    teacher_mask = np.asarray(
        np.load(truth_cache / "train_mask.npy", allow_pickle=False),
        dtype=bool,
    )
    if len(train_targets) != expected_training_targets:
        raise ValueError(
            f"expected {expected_training_targets} non-panel targets, "
            f"got {len(train_targets)}"
        )
    if len(outer_targets) < minimum_outer_targets:
        raise ValueError(f"only {len(outer_targets)} shared VCC targets")
    if set(train_targets) & set(panel):
        raise ValueError("non-panel targets overlap the VCC panel")
    if not set(outer_targets) <= set(panel):
        raise ValueError("outer targets are not a VCC-panel subset")
    if not set(master_genes) <= set(official):
        raise ValueError("truth-cache genes are outside the official axis")
    if set(master_genes) & set(panel):
        raise ValueError("a VCC target coordinate entered the panel score axis")
    if teacher_mask.shape != (len(master_genes),):
        raise ValueError("teacher-common score mask width mismatch")
    if np.any(
        teacher_mask
        & np.asarray(
            [gene in set(train_targets) for gene in master_genes],
            dtype=bool,
        )
    ):
        raise ValueError("training target coordinate entered teacher-common score mask")
    return master_genes, train_targets, outer_targets, teacher_mask


def established_final_score_axis(
    official: Sequence[str],
    panel: Sequence[str],
    donor_root: Path,
    manifest: Mapping[str, object],
) -> tuple[tuple[str, ...], np.ndarray, dict[str, object]]:
    """Recover and authenticate the complete held-donor final score axis."""
    panel_set = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    if len(official_pos) != len(official):
        raise ValueError("official gene axis contains duplicate symbols")
    missing_panel = sorted(panel_set - set(official_pos))
    if missing_panel:
        raise ValueError(
            f"{len(missing_panel)} VCC targets are absent from the official axis"
        )

    axes: dict[str, tuple[str, ...]] = {}
    donor_axis_audits: dict[str, object] = {}
    for donor in DEVELOPMENT_FOLDS:
        gene_rows = EXP062.read_lines(
            donor_root / donor / "genes.txt",
            require_unique=False,
        )
        gene_lookup, gene_axis_audit = last_occurrence_lookup(gene_rows)
        donor_axis_audits[donor] = gene_axis_audit
        axes[donor] = tuple(
            gene for gene in official if gene in gene_lookup and gene not in panel_set
        )
    reference = axes[DEVELOPMENT_FOLDS[0]]
    if any(axis != reference for axis in axes.values()):
        raise ValueError("development held-donor score axes are not identical")

    frozen = manifest.get("final_score_axis", {})
    expected_count = int(frozen.get("gene_count", -1))
    expected_hash = str(frozen.get("gene_order_sha256", ""))
    expected_positions_hash = str(frozen.get("official_positions_sha256", ""))
    positions = np.asarray([official_pos[gene] for gene in reference], dtype=np.int64)
    if (
        len(official) == 18_533
        and len(panel) == 300
        and expected_count != EXPECTED_FINAL_SCORE_GENES
    ):
        raise ValueError(
            "manifest final score-axis count is not the established 17,806"
        )
    if len(reference) != expected_count:
        raise ValueError(
            f"held-donor score axis {len(reference)} != frozen {expected_count}"
        )
    if sha256_lines(reference) != expected_hash:
        raise ValueError("held-donor score-gene order/hash changed")
    if sha256_array(positions) != expected_positions_hash:
        raise ValueError("held-donor official-position order/hash changed")
    return (
        reference,
        positions,
        {
            "gene_count": len(reference),
            "gene_order_sha256": sha256_lines(reference),
            "official_positions_sha256": sha256_array(positions),
            "development_donor_axes_equal": True,
            "excludes_exactly_vcc_coordinates_within_held_support": True,
            "donor_gene_axes": donor_axis_audits,
        },
    )


def generalized_pm_slot(
    incumbent: EXP062.PreparedSource,
    perturbmatch: EXP062.PreparedSource,
    shared_columns: np.ndarray,
) -> tuple[EXP062.PreparedSource, dict[str, int]]:
    """Use PM where both target rows exist; preserve incumbent otherwise."""
    effects = np.asarray(incumbent.effects, dtype=np.float64).copy()
    rows = incumbent.target_available & perturbmatch.target_available
    row_indices = np.flatnonzero(rows)
    effects[np.ix_(row_indices, shared_columns)] = (
        0.5 * incumbent.effects[np.ix_(row_indices, shared_columns)]
        + 0.5 * perturbmatch.effects[np.ix_(row_indices, shared_columns)]
    )
    slot = EXP062.PreparedSource(
        name="k562_gwps_pm50",
        effects=effects,
        target_available=np.asarray(incumbent.target_available, dtype=bool).copy(),
        mean_on_axis=np.asarray(incumbent.mean_on_axis, dtype=np.float64).copy(),
        raw_median_norm=float(incumbent.raw_median_norm),
        mapped_gene_count=int(incumbent.mapped_gene_count),
        duplicate_symbols=tuple(incumbent.duplicate_symbols),
    )
    return slot, {
        "targets_total": int(len(rows)),
        "targets_pm_modified": int(np.sum(rows)),
        "targets_incumbent_preserved": int(np.sum(~rows & incumbent.target_available)),
    }


def prepare_main_source_safe(
    bank: EXP062.Bank,
    targets: Sequence[str],
    official_genes: Sequence[str],
    blocked: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> EXP062.PreparedSource:
    """Prepare a source with blocked coordinates absent from every statistic."""
    official_pos = {gene: index for index, gene in enumerate(official_genes)}
    all_source_columns: list[int] = []
    last_source: dict[int, int] = {}
    duplicates: set[str] = set()
    seen: set[int] = set()
    for source_column, gene in enumerate(bank.genes):
        destination = official_pos.get(gene)
        if destination is None or gene in blocked:
            continue
        all_source_columns.append(source_column)
        if destination in seen:
            duplicates.add(gene)
        seen.add(destination)
        last_source[destination] = source_column
    if not all_source_columns:
        raise ValueError(f"{bank.name}: no leakage-safe genes map to official axis")
    all_columns = np.asarray(all_source_columns, dtype=np.int64)
    unique_positions = np.asarray(sorted(last_source), dtype=np.int64)
    unique_columns = np.asarray(
        [last_source[position] for position in unique_positions], dtype=np.int64
    )
    mean_all = EXP062.stream_column_mean(
        bank,
        all_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    raw_median, _, _ = EXP062.centered_statistics(
        bank,
        all_columns,
        mean_all,
        need_gene_scale=False,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    all_position = {
        int(source): index for index, source in enumerate(all_columns.tolist())
    }
    unique_mean = np.asarray(
        [mean_all[all_position[int(source)]] for source in unique_columns],
        dtype=np.float64,
    )
    mean_on_axis = np.zeros(len(official_genes), dtype=np.float64)
    mean_on_axis[unique_positions] = unique_mean
    effects = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    available = np.zeros(len(targets), dtype=bool)
    delta = np.load(bank.delta_path, mmap_mode="r")
    for index, target in enumerate(targets):
        row = bank.pert_pos.get(target)
        if row is None:
            continue
        values = np.asarray(delta[row, unique_columns], dtype=np.float64)
        effects[index, unique_positions] = (values - unique_mean) * (
            EXP062.EFFECT_SCALE / raw_median
        )
        available[index] = True
    return EXP062.PreparedSource(
        name=bank.name,
        effects=effects,
        target_available=available,
        mean_on_axis=mean_on_axis,
        raw_median_norm=raw_median,
        mapped_gene_count=len(unique_positions),
        duplicate_symbols=tuple(sorted(duplicates)),
    )


def slice_source(
    source: EXP062.PreparedSource, rows: slice | np.ndarray
) -> EXP062.PreparedSource:
    return EXP062.PreparedSource(
        name=source.name,
        effects=np.asarray(source.effects[rows], dtype=np.float64),
        target_available=np.asarray(source.target_available[rows], dtype=bool),
        mean_on_axis=source.mean_on_axis,
        raw_median_norm=source.raw_median_norm,
        mapped_gene_count=source.mapped_gene_count,
        duplicate_symbols=source.duplicate_symbols,
    )


def slice_cd4(
    source: EXP062.CD4Prepared, rows: slice | np.ndarray
) -> EXP062.CD4Prepared:
    return EXP062.CD4Prepared(
        mode=source.mode,
        mean_effect=np.asarray(source.mean_effect[rows], dtype=np.float64),
        gene_votes=np.asarray(source.gene_votes[rows]),
        target_condition_count=np.asarray(source.target_condition_count[rows]),
        condition_details=source.condition_details,
    )


def prepare_arm_inputs(
    banks: Mapping[str, EXP062.Bank],
    perturbmatch_bank: EXP062.Bank,
    cd4_banks: Sequence[EXP062.Bank],
    source_names: Sequence[str],
    targets: Sequence[str],
    official: Sequence[str],
    blocked: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[
    dict[str, EXP062.PreparedSource],
    EXP062.PreparedSource,
    EXP062.PreparedSource,
    EXP062.CD4Prepared,
    np.ndarray,
    dict[str, object],
]:
    prepared = {
        name: prepare_main_source_safe(
            banks[name],
            targets,
            official,
            blocked,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        for name in source_names
    }
    prepared_pm = prepare_main_source_safe(
        perturbmatch_bank,
        targets,
        official,
        blocked,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    k562_genes = set(banks["k562_gwps"].genes)
    pm_genes = set(perturbmatch_bank.genes)
    shared_columns = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene not in blocked and gene in k562_genes and gene in pm_genes
        ],
        dtype=np.int64,
    )
    pm_slot, pm_support = generalized_pm_slot(
        prepared["k562_gwps"], prepared_pm, shared_columns
    )
    cd4 = EXP062.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        blocked,
        EXP064.CD4_MODE,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    missing_pm = prepared["k562_gwps"].target_available & ~prepared_pm.target_available
    if np.any(pm_slot.effects[missing_pm] != prepared["k562_gwps"].effects[missing_pm]):
        raise AssertionError("PerturbMatch fallback changed incumbent-only rows")
    return (
        prepared,
        prepared_pm,
        pm_slot,
        cd4,
        shared_columns,
        {
            "blocked_target_labels": len(blocked),
            "blocked_target_labels_sha256": sha256_lines(sorted(blocked)),
            "blocked_coordinates": int(sum(gene in set(official) for gene in blocked)),
            "main_source_statistics_exclude_blocked_coordinates": True,
            "cd4_statistics_exclude_blocked_coordinates": True,
            "pm_shared_unblocked_genes": int(len(shared_columns)),
            "pm_shared_columns_sha256": sha256_array(shared_columns),
            "pm_support": pm_support,
            "pm_fallback_rows": int(np.sum(missing_pm)),
            "pm_fallback_exact_incumbent": True,
            "cd4_target_condition_histogram": {
                str(value): int(np.sum(cd4.target_condition_count == value))
                for value in sorted(set(cd4.target_condition_count.tolist()))
            },
            "source_coverage": {
                name: int(np.sum(prepared[name].target_available))
                for name in source_names
            },
        },
    )


def exp064_reproduction_audit(
    *,
    held: str,
    reference_targets: Sequence[str],
    evaluated_targets: Sequence[str],
    held_genes: set[str],
    prepared: Mapping[str, EXP062.PreparedSource],
    prepared_pm: EXP062.PreparedSource,
    pm_slot: EXP062.PreparedSource,
    cd4: EXP062.CD4Prepared,
    operator: EXP062.Operator,
    operator_columns: np.ndarray,
    staged_query: np.ndarray,
    staged_score_positions: np.ndarray,
    staged_score_genes: Sequence[str],
    official: Sequence[str],
    shared_columns: np.ndarray,
    atol: float,
) -> dict[str, object]:
    row_pos = {target: index for index, target in enumerate(evaluated_targets)}
    missing = [target for target in reference_targets if target not in row_pos]
    if missing:
        raise ValueError(f"{held}: exp064 targets absent from panel: {missing[:5]}")
    rows = np.asarray([row_pos[target] for target in reference_targets], dtype=np.int64)
    companion = next(name for name in DEVELOPMENT_FOLDS if name != held)
    incumbent_ref = slice_source(prepared["k562_gwps"], rows)
    pm_ref = slice_source(prepared_pm, rows)
    companion_ref = slice_source(prepared[companion], rows)
    cd4_ref = slice_cd4(cd4, rows)
    exact_pm_slot = EXP064.build_pm_slot(incumbent_ref, pm_ref, shared_columns)
    exact_pre, exact_votes = EXP062.build_baseline([exact_pm_slot, companion_ref])
    exact_white, _ = EXP062.apply_operator(exact_pre, operator, operator_columns)
    exact_arm_d, _ = EXP062.add_post_whiten_residual(
        exact_white,
        cd4_ref,
        EXP064.CD4_FRACTION,
        operator_columns,
    )
    retrieval = np.asarray(
        [index for index in operator_columns if official[int(index)] in held_genes],
        dtype=np.int64,
    )
    if not np.array_equal(
        np.asarray(staged_score_positions, dtype=np.int64),
        retrieval,
    ):
        raise AssertionError(f"{held}: staged score axis differs from arm-D retrieval")
    retrieval_genes = tuple(official[int(index)] for index in retrieval)
    if tuple(staged_score_genes) != retrieval_genes:
        raise AssertionError(f"{held}: staged score-gene order differs from arm D")
    observed = np.asarray(staged_query[rows], dtype=np.float64)
    expected = np.asarray(exact_arm_d[:, retrieval], dtype=np.float64)
    max_abs_error = float(np.max(np.abs(observed - expected)))
    if not np.allclose(observed, expected, rtol=0.0, atol=atol):
        raise AssertionError(
            f"{held}: arm-D algebra differs from direct reference by {max_abs_error}"
        )
    if not np.array_equal(pm_slot.effects[rows], exact_pm_slot.effects):
        raise AssertionError(f"{held}: generalized PM slot differs on exp064 rows")
    extension = tuple(
        target for target in evaluated_targets if target not in set(reference_targets)
    )
    extension_rows = np.asarray(
        [row_pos[target] for target in extension], dtype=np.int64
    )
    return {
        "passed": True,
        "comparison": (
            "tight numeric equality of saved panel_query rows to the direct "
            "arm-D formula on the complete retrieval axis"
        ),
        "absolute_tolerance": float(atol),
        "max_abs_error": max_abs_error,
        "validated_targets": len(reference_targets),
        "validated_target_order_sha256": sha256_lines(reference_targets),
        "retrieval_genes": int(len(retrieval)),
        "retrieval_columns_sha256": sha256_array(retrieval),
        "retrieval_gene_order_sha256": sha256_lines(retrieval_genes),
        "saved_query_score_axis_exact": True,
        "leakage_safe_source_statistics": True,
        "source_vote_histogram": {
            str(value): int(np.sum(exact_votes == value))
            for value in sorted(set(exact_votes.tolist()))
        },
        "extension_targets": len(extension),
        "extension_target_order_sha256": sha256_lines(extension),
        "extension_audit": {
            "pm_missing_rows": int(
                np.sum(
                    prepared["k562_gwps"].target_available[extension_rows]
                    & ~prepared_pm.target_available[extension_rows]
                )
            ),
            "cd4_condition_histogram": {
                str(value): int(
                    np.sum(cd4.target_condition_count[extension_rows] == value)
                )
                for value in sorted(
                    set(cd4.target_condition_count[extension_rows].tolist())
                )
            },
            "fallback_policy": (
                "preserve incumbent K562 row when PerturbMatch lacks a target; "
                "average only available main donors; add CD4 only where measured"
            ),
        },
    }


def build_arm_d_panel(
    pm_slot: EXP062.PreparedSource,
    other_sources: Sequence[EXP062.PreparedSource],
    cd4: EXP062.CD4Prepared,
    operator: EXP062.Operator,
    operator_columns: np.ndarray,
    active_operator_columns: np.ndarray,
    preserve_rows: np.ndarray,
    score_positions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Apply arm D with a leakage-safe principal operator restriction."""
    pre, votes = EXP062.build_baseline([pm_slot, *other_sources])
    whitened, whitening = apply_operator_restriction(
        pre,
        operator,
        operator_columns,
        active_operator_columns,
        preserve_rows=preserve_rows,
    )
    arm_d, cd4_preserve = add_post_whiten_residual(
        whitened,
        cd4,
        active_operator_columns,
        preserve_rows=preserve_rows,
    )
    output = np.asarray(arm_d[:, score_positions], dtype=np.float32)
    if not np.isfinite(output).all():
        raise ValueError("arm-D query contains non-finite values")
    return (
        output,
        arm_d,
        {
            "rows": int(output.shape[0]),
            "columns": int(output.shape[1]),
            "source_vote_histogram": {
                str(value): int(np.sum(votes == value))
                for value in sorted(set(votes.tolist()))
            },
            "whitening": whitening,
            "operator_restriction": whitening,
            "cd4_post_whiten": cd4_preserve,
            "cd4_targets_with_any_condition": int(
                np.sum(cd4.target_condition_count > 0)
            ),
            "cd4_targets_with_all_conditions": int(
                np.sum(cd4.target_condition_count == len(CD4_DONORS))
            ),
        },
    )


def stage(
    args: argparse.Namespace, manifest: Mapping[str, object]
) -> dict[str, object]:
    official = EXP062.read_csv_column(args.genes_csv, "gene_name")
    panel = EXP062.read_csv_column(args.panel_csv, "target_gene")
    master_genes, train_targets, outer_targets, teacher_mask = validate_truth_contract(
        args.truth_cache,
        official,
        panel,
        expected_training_targets=args.expected_training_targets,
        minimum_outer_targets=args.minimum_outer_targets,
    )
    exp064_targets = load_exp064_reference(args.exp064_results, manifest)
    required_atol = float(manifest["exp064_reference"]["required_max_abs_error"])
    if args.exp064_atol > required_atol:
        raise ValueError(
            f"--exp064-atol {args.exp064_atol} exceeds frozen limit {required_atol}"
        )
    if not set(exp064_targets) <= set(outer_targets):
        raise ValueError("experiment-064 target gallery is not in outer target cache")
    official_pos = {gene: index for index, gene in enumerate(official)}
    train_set = set(train_targets)
    panel_set = set(panel)
    train_blocked = panel_set | train_set
    panel_blocked = panel_set
    missing_panel_coordinates = sorted(panel_blocked - set(official_pos))
    if missing_panel_coordinates:
        raise ValueError(
            f"{len(missing_panel_coordinates)} VCC targets are absent from the "
            f"official axis; first={missing_panel_coordinates[:5]}"
        )
    missing_nonpanel_coordinates = tuple(sorted(train_set - set(official_pos)))
    missing_nonpanel_master_coordinates = tuple(sorted(train_set - set(master_genes)))
    present_train_blocked = train_blocked & set(official_pos)
    teacher_gene_sets = {
        donor: set(read_lines(args.donor_root / donor / "genes.txt"))
        for donor in TEACHER_DONORS
    }
    expected_teacher_genes = (
        set(master_genes)
        .intersection(*(teacher_gene_sets[donor] for donor in TEACHER_DONORS))
        .difference(train_set)
    )
    dev_score_genes = tuple(
        gene for gene, keep in zip(master_genes, teacher_mask, strict=True) if keep
    )
    if set(dev_score_genes) != expected_teacher_genes:
        raise ValueError(
            "teacher-common score mask does not equal the measured-gene "
            "intersection after target-coordinate exclusion"
        )
    panel_score_genes, panel_score_positions, final_axis_audit = (
        established_final_score_axis(
            official,
            panel,
            args.donor_root,
            manifest,
        )
    )
    dev_score_positions = np.asarray(
        [official_pos[gene] for gene in dev_score_genes], dtype=np.int64
    )
    panel_position = {gene: index for index, gene in enumerate(panel_score_genes)}
    missing_dev_score = [gene for gene in dev_score_genes if gene not in panel_position]
    if missing_dev_score:
        raise ValueError(
            f"{len(missing_dev_score)} development score genes are absent from "
            "the frozen final score axis"
        )
    dev_to_panel_positions = np.asarray(
        [panel_position[gene] for gene in dev_score_genes], dtype=np.int64
    )
    if set(panel_score_genes) & panel_set:
        raise AssertionError("VCC coordinate entered final panel score axis")
    if set(dev_score_genes) & (panel_set | train_set):
        raise AssertionError("held target coordinate entered development score axis")
    historical_panel_coordinates = set(panel_score_genes) & train_set
    expected_historical_coordinates = train_set & set(official_pos)
    if historical_panel_coordinates != expected_historical_coordinates:
        raise AssertionError(
            "final panel score axis did not retain every available historical "
            "target coordinate"
        )

    if args.output.exists() and any(args.output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty cache {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    np.save(
        args.output / "dev_to_panel_positions.npy",
        dev_to_panel_positions,
        allow_pickle=False,
    )

    fold_metadata: dict[str, object] = {}
    fold_specs = {
        held: {
            "sources": (
                "k562_gwps",
                next(name for name in DEVELOPMENT_FOLDS if name != held),
            ),
        }
        for held in DEVELOPMENT_FOLDS
    }
    fold_specs[KOLF_FOLD] = {"sources": MAIN_DONORS}
    for held in FINAL_FOLDS:
        source_names = tuple(fold_specs[held]["sources"])
        banks = {
            name: EXP062.Bank.load(name, args.donor_root / name)
            for name in source_names
        }
        if held in DEVELOPMENT_FOLDS:
            held_target_set = set(read_lines(args.donor_root / held / "perts.txt"))
            evaluation_targets = tuple(
                target for target in panel if target in held_target_set
            )
            if len(evaluation_targets) < args.minimum_outer_targets:
                raise ValueError(
                    f"{held}: only {len(evaluation_targets)} VCC-panel targets"
                )
            if not set(exp064_targets) <= set(evaluation_targets):
                raise ValueError(
                    f"{held}: experiment-064 gallery is not a subset of the "
                    "complete held-donor panel"
                )
        else:
            evaluation_targets = panel
        perturbmatch_bank = EXP064.load_perturbmatch(args.perturbmatch)
        cd4_banks = tuple(
            EXP062.Bank.load(
                name,
                args.cd4_root / name,
                require_cd4_scale=True,
            )
            for name in CD4_DONORS
        )
        operator_entry = manifest["operators"][held]
        operator, operator_columns = load_operator(
            args.operator_dir / f"held_{held}.npz",
            tuple(banks[name] for name in source_names),
            official,
            panel_set,
            operator_entry,
        )
        train_active_columns = np.asarray(
            [
                position
                for position in operator_columns
                if official[int(position)] not in train_set
            ],
            dtype=np.int64,
        )
        if not set(map(int, dev_score_positions)) <= set(
            map(int, train_active_columns)
        ):
            raise ValueError(f"{held}: development score axis is outside operator")
        if not set(map(int, panel_score_positions)) <= set(map(int, operator_columns)):
            raise ValueError(f"{held}: final panel score axis is outside operator")

        fold_dir = args.output / held
        fold_dir.mkdir(parents=True, exist_ok=True)
        panel_targets_path = fold_dir / "panel_targets.txt"
        panel_targets_path.write_text("\n".join(evaluation_targets) + "\n")
        (
            train_prepared,
            _train_prepared_pm,
            train_pm_slot,
            train_cd4,
            _train_shared_columns,
            train_support,
        ) = prepare_arm_inputs(
            banks,
            perturbmatch_bank,
            cd4_banks,
            source_names,
            train_targets,
            official,
            train_blocked,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        train_query, _train_full, train_audit = build_arm_d_panel(
            train_pm_slot,
            [train_prepared[name] for name in source_names if name != "k562_gwps"],
            train_cd4,
            operator,
            operator_columns,
            train_active_columns,
            np.arange(len(train_targets), dtype=np.int64),
            dev_score_positions,
        )
        train_source_count = np.sum(
            np.stack(
                [
                    train_pm_slot.target_available,
                    *[
                        train_prepared[name].target_available
                        for name in source_names
                        if name != "k562_gwps"
                    ],
                ]
            ),
            axis=0,
            dtype=np.uint8,
        )
        train_active_genes = {
            official[int(position)] for position in train_active_columns
        }
        train_audit["operator_restriction"][
            "nonpanel_target_overlap_with_active_operator"
        ] = len(train_active_genes & train_set)
        train_audit["operator_restriction"][
            "all_nonpanel_target_coordinates_excluded"
        ] = not bool(train_active_genes & train_set)
        train_audit["operator_restriction"]["excluded_nonpanel_targets_sha256"] = (
            sha256_lines(sorted(train_set))
        )
        train_audit["operator_restriction"][
            "nonpanel_targets_without_official_coordinate"
        ] = len(missing_nonpanel_coordinates)
        train_audit["operator_restriction"][
            "nonpanel_targets_without_official_coordinate_sha256"
        ] = sha256_lines(missing_nonpanel_coordinates)
        train_audit["operator_restriction"][
            "nonpanel_targets_without_master_coordinate"
        ] = len(missing_nonpanel_master_coordinates)
        train_audit["operator_restriction"][
            "nonpanel_targets_without_master_coordinate_sha256"
        ] = sha256_lines(missing_nonpanel_master_coordinates)
        train_audit["support"] = train_support

        (
            panel_prepared,
            panel_prepared_pm,
            panel_pm_slot,
            panel_cd4,
            panel_shared_columns,
            panel_support,
        ) = prepare_arm_inputs(
            banks,
            perturbmatch_bank,
            cd4_banks,
            source_names,
            evaluation_targets,
            official,
            panel_blocked,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        if held in DEVELOPMENT_FOLDS:
            evaluation_pos = {
                target: index for index, target in enumerate(evaluation_targets)
            }
            preserve_rows = np.asarray(
                [evaluation_pos[target] for target in exp064_targets],
                dtype=np.int64,
            )
        else:
            preserve_rows = np.arange(len(evaluation_targets), dtype=np.int64)
        panel_query, panel_full, panel_audit = build_arm_d_panel(
            panel_pm_slot,
            [panel_prepared[name] for name in source_names if name != "k562_gwps"],
            panel_cd4,
            operator,
            operator_columns,
            operator_columns,
            preserve_rows,
            panel_score_positions,
        )
        panel_source_count = np.sum(
            np.stack(
                [
                    panel_pm_slot.target_available,
                    *[
                        panel_prepared[name].target_available
                        for name in source_names
                        if name != "k562_gwps"
                    ],
                ]
            ),
            axis=0,
            dtype=np.uint8,
        )
        panel_audit["support"] = panel_support

        train_path = fold_dir / "train_query.npy"
        train_count_path = fold_dir / "train_source_count.npy"
        panel_path = fold_dir / "panel_query.npy"
        panel_count_path = fold_dir / "panel_source_count.npy"
        np.save(train_path, train_query, allow_pickle=False)
        np.save(train_count_path, train_source_count, allow_pickle=False)
        np.save(panel_path, panel_query, allow_pickle=False)
        np.save(panel_count_path, panel_source_count, allow_pickle=False)
        train_audit.update(
            {
                "query_sha256": sha256_file(train_path),
                "source_count_sha256": sha256_file(train_count_path),
                "source_count_histogram": {
                    str(value): int(np.sum(train_source_count == value))
                    for value in range(len(source_names) + 1)
                },
                "held_truth_opened": False,
            }
        )
        panel_audit.update(
            {
                "query_sha256": sha256_file(panel_path),
                "source_count_sha256": sha256_file(panel_count_path),
                "source_count_histogram": {
                    str(value): int(np.sum(panel_source_count == value))
                    for value in range(len(source_names) + 1)
                },
                "target_order_sha256": sha256_lines(evaluation_targets),
                "score_axis_excludes_exactly_vcc_coordinates": True,
                "score_gene_count": len(panel_score_genes),
                "score_gene_order_sha256": sha256_lines(panel_score_genes),
                "score_official_positions_sha256": sha256_array(panel_score_positions),
                "saved_query_is_complete_scored_axis": True,
                "held_truth_opened": False,
            }
        )

        reproduction = None
        if held in DEVELOPMENT_FOLDS:
            reproduction = exp064_reproduction_audit(
                held=held,
                reference_targets=exp064_targets,
                evaluated_targets=evaluation_targets,
                held_genes={
                    line.strip()
                    for line in (args.donor_root / held / "genes.txt")
                    .read_text()
                    .splitlines()
                    if line.strip()
                },
                prepared=panel_prepared,
                prepared_pm=panel_prepared_pm,
                pm_slot=panel_pm_slot,
                cd4=panel_cd4,
                operator=operator,
                operator_columns=operator_columns,
                staged_query=panel_query,
                staged_score_positions=panel_score_positions,
                staged_score_genes=panel_score_genes,
                official=official,
                shared_columns=panel_shared_columns,
                atol=args.exp064_atol,
            )

        fold_metadata[held] = {
            "held_donor": held,
            "source_donors": list(source_names),
            "held_donor_used_in_prediction": False,
            "held_donor_used_in_operator": False,
            "held_donor_delta_loaded_in_own_fold": False,
            "operator_npz_sha256": operator_entry["npz_sha256"],
            "operator_json_sha256": operator_entry["json_sha256"],
            "operator_columns": int(len(operator_columns)),
            "operator_columns_sha256": sha256_array(operator_columns),
            "panel_targets": len(evaluation_targets),
            "panel_target_order_sha256": sha256_lines(evaluation_targets),
            "panel_targets_file_sha256": sha256_file(panel_targets_path),
            "complete_held_panel": held in DEVELOPMENT_FOLDS,
            "panels": {
                "train": train_audit,
                "panel": panel_audit,
            },
        }
        if reproduction is not None:
            fold_metadata[held]["exp064_reproduction"] = reproduction
        del (
            banks,
            perturbmatch_bank,
            cd4_banks,
            train_prepared,
            _train_prepared_pm,
            train_pm_slot,
            train_cd4,
            panel_prepared,
            panel_prepared_pm,
            panel_pm_slot,
            panel_cd4,
            operator,
            train_query,
            panel_query,
            panel_full,
            train_source_count,
            panel_source_count,
        )
        gc.collect()
        enforce_memory_limit(args.max_rss_gib, f"{held} completion")

    metadata: dict[str, object] = {
        "schema_version": CACHE_SCHEMA,
        "experiment": "082-biological-target-kernel",
        "status": "arm_d_query_cache_complete",
        "arm": "D",
        "arm_definition": (
            "50:50 incumbent/PerturbMatch K562 slot, source-only gamma=0.65 "
            "union whitening, then CD4 mean3 gene_norm residual at fraction 0.05"
        ),
        "implementation": {
            "exp062_run_sha256": sha256_file(EXP062_PATH),
            "exp064_run_sha256": sha256_file(EXP064_PATH),
            "exact_arm_d_formula": True,
            "leakage_safe_preparation": (
                "blocked coordinates are excluded from source means, scales, "
                "effects, CD4 statistics, and whitening inputs; development "
                "uses the principal restriction of the frozen dense operator"
            ),
        },
        "development_folds": list(DEVELOPMENT_FOLDS),
        "final_folds": list(FINAL_FOLDS),
        "valid_outer_donors": list(OUTER_DONORS),
        "invalid_outer_donor": (
            "k562_gwps omitted because PerturbMatch is a reanalysis of the same assay"
        ),
        "training_targets": len(train_targets),
        "outer_targets": len(outer_targets),
        "outer_targets_role": (
            "legacy shared HCT116/HEK293T reference used only to validate the "
            "experiment-064 reproduction subset"
        ),
        "training_outer_overlap": len(set(train_targets) & set(outer_targets)),
        "vcc_targets": len(panel),
        "training_blocked_targets": len(train_blocked),
        "training_blocked_target_coordinates": len(present_train_blocked),
        "training_blocked_target_coordinates_sha256": sha256_lines(
            sorted(present_train_blocked)
        ),
        "training_targets_without_official_coordinate": len(
            missing_nonpanel_coordinates
        ),
        "training_targets_without_official_coordinate_sha256": sha256_lines(
            missing_nonpanel_coordinates
        ),
        "training_targets_without_master_coordinate": len(
            missing_nonpanel_master_coordinates
        ),
        "training_targets_without_master_coordinate_sha256": sha256_lines(
            missing_nonpanel_master_coordinates
        ),
        "panel_blocked_target_coordinates": len(panel_blocked),
        "panel_blocked_target_coordinates_sha256": sha256_lines(sorted(panel_blocked)),
        "historical_target_coordinates_restored_for_panel": len(
            historical_panel_coordinates
        ),
        "master_genes": len(master_genes),
        "dev_score_genes": len(dev_score_genes),
        "panel_score_genes": len(panel_score_genes),
        "dev_score_gene_order_sha256": sha256_lines(dev_score_genes),
        "panel_score_gene_order_sha256": sha256_lines(panel_score_genes),
        "panel_score_official_positions_sha256": sha256_array(panel_score_positions),
        "dev_to_panel_positions_sha256": sha256_array(dev_to_panel_positions),
        "train_order_sha256": sha256_lines(train_targets),
        "outer_order_sha256": sha256_lines(outer_targets),
        "panel_order_sha256": sha256_lines(panel),
        "exp064_target_order_sha256": sha256_lines(exp064_targets),
        "coordinate_contract": {
            "development_blocked": {
                "vcc_targets": len(panel),
                "nonpanel_development_targets": len(train_targets),
                "nonpanel_development_target_coordinates": len(
                    train_set & set(official_pos)
                ),
                "before_dense_whitening": True,
                "before_source_normalization": True,
                "score_axis_excludes_both_sets": True,
            },
            "final_panel_blocked": {
                "vcc_targets": len(panel),
                "nonpanel_development_targets": 0,
                "before_dense_whitening": True,
                "before_source_normalization": True,
                "score_axis_excludes_exactly_vcc_targets": True,
            },
        },
        "folds": fold_metadata,
        "panel_query_target_counts": {
            held: int(fold_metadata[held]["panel_targets"]) for held in FINAL_FOLDS
        },
        "panel_query_target_order_sha256": {
            held: str(fold_metadata[held]["panel_target_order_sha256"])
            for held in FINAL_FOLDS
        },
        "final_score_axis": final_axis_audit,
        "leakage_audit": {
            "all_vcc_target_coordinates_removed": True,
            "all_nonpanel_target_coordinates_removed_from_training": True,
            "training_dense_whitener_principally_restricted": True,
            "development_target_coordinates_excluded_from_development_score_axis": True,
            "teacher_common_measured_gene_axis_used_for_all_residual_fits": True,
            "development_target_coordinates_retained_on_final_panel_score_axis": True,
            "final_panel_excludes_exactly_the_vcc_target_set": True,
            "coordinate_less_nonpanel_targets_treated_as_no_coordinate_to_remove": True,
            "source_normalization_excludes_blocked_coordinates": True,
            "saved_panel_query_uses_complete_17806_gene_score_axis": True,
            "held_truth_opened": False,
            "held_donor_used_in_query": False,
            "held_donor_used_in_operator": False,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-root", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--truth-cache", type=Path, required=True)
    parser.add_argument("--operator-dir", type=Path, required=True)
    parser.add_argument("--exp064-results", type=Path, required=True)
    parser.add_argument("--artifact-manifest", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=48.0)
    parser.add_argument("--expected-training-targets", type=int, default=600)
    parser.add_argument("--minimum-outer-targets", type=int, default=240)
    parser.add_argument("--exp064-atol", type=float, default=1e-6)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data cache staging without --guarded")
    if (
        args.chunk_rows < 1
        or args.max_rss_gib <= 0
        or args.expected_training_targets < 1
        or args.minimum_outer_targets < 2
        or args.exp064_atol <= 0
    ):
        raise SystemExit("chunk rows and memory limit must be positive")
    manifest = json.loads(args.artifact_manifest.read_text())
    if int(manifest.get("schema_version", -1)) != 1:
        raise SystemExit("unsupported artifact manifest")
    verify_inputs(args, manifest)
    metadata = stage(args, manifest)
    print(json.dumps(metadata["leakage_audit"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
