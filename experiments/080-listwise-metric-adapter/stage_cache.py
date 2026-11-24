#!/usr/bin/env python3
"""Stage leakage-safe frozen-arm-D query matrices for experiment 080.

This script reconstructs experiment-064 arm D for the two development
holdouts and the source-only three-donor arm-D query needed for independent
KOLF confirmation.  Training and final-panel queries deliberately use
different score axes:

* non-panel fitting/evaluation removes the 300 VCC and 600 historical target
  coordinates before a principal restriction of the frozen dense whitener;
* final VCC-panel evaluation reproduces experiment 064's fixed baseline:
  source centering and norm calibration use all mapped coordinates, then the
  300 VCC coordinates are zeroed. Historical target coordinates remain
  ordinary scored genes.

Held-donor effects, including KOLF, are never opened here.
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
CD4_DONORS = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CACHE_SCHEMA = 5
EXPECTED_FINAL_SCORE_GENES = 17_806


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP062 = load_module("exp080_exp062", ROOT / "experiments/062-cd4-post-whiten/run.py")
EXP064 = load_module(
    "exp080_exp064", ROOT / "experiments/064-pm-cd4-composition/run.py"
)


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
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def verify_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 {actual} != expected {expected}")


def verify_inputs(args: argparse.Namespace, manifest: Mapping[str, object]) -> None:
    reference = manifest["exp064_reference"]
    expected_reference_contract = {
        "normalization_contract": "all_mapped_then_block_vcc",
        "independent_raw_bank_reconstruction_required": True,
        "learned_statistics_use_separate_retained_coordinate_features": True,
    }
    for key, value in expected_reference_contract.items():
        if reference.get(key) != value:
            raise ValueError(f"exp064 reference contract {key!r} changed")
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
            # The cache builder needs labels and masks, never held panel truth.
            # Outer truth is authenticated only after run.py's pre-gate passes.
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
    """Load either an exp062 held-fold or the production three-donor operator."""
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

    If ``W = E diag(w) E.T`` and forbidden input coordinates are zero, then
    the retained output is exactly ``x_R W_RR`` where
    ``W_RR = E_R diag(w) E_R.T``.  Computing with ``E_R`` directly prevents
    forbidden historical target coordinates from being materialized as
    post-whitening outputs.
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
    baseline_norm = (
        float(
            np.median(
                np.linalg.norm(
                    baseline[np.ix_(preserve, np.asarray(active_columns))], axis=1
                )
            )
        )
        or 1.0
    )
    candidate_norm = (
        float(
            np.median(
                np.linalg.norm(
                    candidate[np.ix_(preserve, np.asarray(active_columns))], axis=1
                )
            )
        )
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
    train_mask = np.asarray(
        np.load(truth_cache / "train_mask.npy", allow_pickle=False), dtype=bool
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
    if len(master_genes) != len(train_mask):
        raise ValueError("truth-cache mask length mismatch")
    if not set(master_genes) <= set(official):
        raise ValueError("truth-cache genes are outside the official axis")
    if set(master_genes) & set(panel):
        raise ValueError("a VCC target coordinate entered the score axis")
    if np.any(
        train_mask
        & np.asarray([gene in set(train_targets) for gene in master_genes], dtype=bool)
    ):
        raise ValueError("a fitting-target coordinate entered the correction axis")
    if not np.any(train_mask):
        raise ValueError("empty correction axis")
    return master_genes, train_targets, outer_targets, train_mask


def established_final_score_axis(
    official: Sequence[str],
    panel: Sequence[str],
    donor_root: Path,
    manifest: Mapping[str, object],
) -> tuple[tuple[str, ...], np.ndarray, dict[str, object]]:
    """Recover and authenticate the established held-truth score axis."""
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
    for donor in DEVELOPMENT_FOLDS:
        genes = set(
            EXP062.read_lines(
                donor_root / donor / "genes.txt",
                require_unique=False,
            )
        )
        axes[donor] = tuple(
            gene for gene in official if gene in genes and gene not in panel_set
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
    official: Sequence[str],
    blocked_positions: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> EXP062.PreparedSource:
    """Center and scale only retained coordinates, then map to the official axis."""
    mapping = EXP062.map_bank_to_axis(bank, official)
    blocked_set = set(map(int, np.asarray(blocked_positions, dtype=np.int64)))
    keep = np.asarray(
        [
            int(position) not in blocked_set
            for position in mapping.unique_target_positions
        ],
        dtype=bool,
    )
    source_columns = np.asarray(mapping.unique_source_columns[keep], dtype=np.int64)
    target_positions = np.asarray(mapping.unique_target_positions[keep], dtype=np.int64)
    if not len(source_columns):
        raise ValueError(f"{bank.name}: no unblocked genes map to the official axis")
    mean = EXP062.stream_column_mean(
        bank,
        source_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    raw_median, _, _ = EXP062.centered_statistics(
        bank,
        source_columns,
        mean,
        need_gene_scale=False,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    mean_on_axis = np.zeros(len(official), dtype=np.float64)
    mean_on_axis[target_positions] = mean
    effects = np.zeros((len(targets), len(official)), dtype=np.float64)
    available = np.zeros(len(targets), dtype=bool)
    delta = np.load(bank.delta_path, mmap_mode="r", allow_pickle=False)
    for index, target in enumerate(targets):
        row = bank.pert_pos.get(target)
        if row is None:
            continue
        values = np.asarray(delta[row, source_columns], dtype=np.float64)
        effects[index, target_positions] = (values - mean) * (
            EXP062.EFFECT_SCALE / raw_median
        )
        available[index] = True
    return EXP062.PreparedSource(
        name=bank.name,
        effects=effects,
        target_available=available,
        mean_on_axis=mean_on_axis,
        raw_median_norm=raw_median,
        mapped_gene_count=len(target_positions),
        duplicate_symbols=mapping.duplicate_symbols,
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
    official_pos = {gene: index for index, gene in enumerate(official)}
    blocked_positions = np.asarray(
        sorted(official_pos[gene] for gene in blocked if gene in official_pos),
        dtype=np.int64,
    )
    prepared = {
        name: prepare_main_source_safe(
            banks[name],
            targets,
            official,
            blocked_positions,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        for name in source_names
    }
    prepared_pm = prepare_main_source_safe(
        perturbmatch_bank,
        targets,
        official,
        blocked_positions,
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
            "normalization_contract": (
                "leakage-safe: center and norm only retained coordinates; "
                "blocked target coordinates never enter a learned statistic"
            ),
            "used_for_adapter_fit_or_selection": True,
            "blocked_coordinates": len(blocked),
            "blocked_coordinates_sha256": sha256_lines(sorted(blocked)),
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


def prepare_exp064_panel_inputs(
    banks: Mapping[str, EXP062.Bank],
    perturbmatch_bank: EXP062.Bank,
    cd4_banks: Sequence[EXP062.Bank],
    source_names: Sequence[str],
    targets: Sequence[str],
    official: Sequence[str],
    panel: Sequence[str],
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
    """Build the frozen exp064 baseline, distinct from leakage-safe fit features.

    Experiment 064 estimated each main donor's centering vector and norm from
    all mapped source coordinates, then zeroed the VCC target coordinates.
    Reproducing that fixed upstream baseline is an inference-only operation:
    adapter bases, weights, hyperparameters, nulls, and clusters are all fitted
    from ``prepare_arm_inputs`` outputs, whose statistics never see either
    target panel.
    """
    panel_set = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[gene] for gene in panel],
        dtype=np.int64,
    )
    prepared = {
        name: EXP062.prepare_main_source(
            banks[name],
            targets,
            official,
            panel_positions,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        for name in source_names
    }
    prepared_pm = EXP062.prepare_main_source(
        perturbmatch_bank,
        targets,
        official,
        panel_positions,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    shared_columns = EXP064.shared_k562_pm_columns(
        banks["k562_gwps"],
        perturbmatch_bank,
        official,
        panel_set,
    )
    pm_slot, pm_support = generalized_pm_slot(
        prepared["k562_gwps"],
        prepared_pm,
        shared_columns,
    )
    cd4 = EXP062.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        panel_set,
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
            "normalization_contract": (
                "exp064 original: center and norm on all mapped source "
                "coordinate, then zero the 300 VCC target coordinates"
            ),
            "inference_only_fixed_baseline": True,
            "used_for_adapter_fit_or_selection": False,
            "blocked_coordinates": len(panel_set),
            "blocked_coordinates_sha256": sha256_lines(sorted(panel_set)),
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
    banks: Mapping[str, EXP062.Bank],
    perturbmatch_bank: EXP062.Bank,
    cd4_banks: Sequence[EXP062.Bank],
    operator: EXP062.Operator,
    operator_columns: np.ndarray,
    saved_query: np.ndarray,
    staged_score_positions: np.ndarray,
    staged_score_genes: Sequence[str],
    official: Sequence[str],
    panel: Sequence[str],
    chunk_rows: int,
    max_rss_gib: float,
    atol: float,
) -> dict[str, object]:
    row_pos = {target: index for index, target in enumerate(evaluated_targets)}
    missing = [target for target in reference_targets if target not in row_pos]
    if missing:
        raise ValueError(
            f"{held}: exp064 targets absent from staged panel: {missing[:5]}"
        )
    rows = np.asarray([row_pos[target] for target in reference_targets], dtype=np.int64)
    companion = next(name for name in DEVELOPMENT_FOLDS if name != held)
    panel_set = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[gene] for gene in panel],
        dtype=np.int64,
    )
    independent_prepared = {
        name: EXP062.prepare_main_source(
            banks[name],
            reference_targets,
            official,
            panel_positions,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        for name in ("k562_gwps", companion)
    }
    independent_pm = EXP062.prepare_main_source(
        perturbmatch_bank,
        reference_targets,
        official,
        panel_positions,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    independent_cd4 = EXP062.prepare_cd4_mean(
        cd4_banks,
        reference_targets,
        official,
        panel_set,
        EXP064.CD4_MODE,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    independent_shared_columns = EXP064.shared_k562_pm_columns(
        banks["k562_gwps"],
        perturbmatch_bank,
        official,
        panel_set,
    )
    exact_pm_slot = EXP064.build_pm_slot(
        independent_prepared["k562_gwps"],
        independent_pm,
        independent_shared_columns,
    )
    exact_pre, exact_votes = EXP062.build_baseline(
        [exact_pm_slot, independent_prepared[companion]]
    )
    exact_white, _ = EXP062.apply_operator(exact_pre, operator, operator_columns)
    exact_arm_d, _ = EXP062.add_post_whiten_residual(
        exact_white,
        independent_cd4,
        EXP064.CD4_FRACTION,
        operator_columns,
    )
    held_genes = set(banks[held].genes)
    retrieval = np.asarray(
        [index for index in operator_columns if official[int(index)] in held_genes],
        dtype=np.int64,
    )
    if not np.array_equal(
        np.asarray(staged_score_positions, dtype=np.int64),
        retrieval,
    ):
        raise AssertionError(f"{held}: staged score axis differs from exp064 retrieval")
    retrieval_genes = tuple(official[int(index)] for index in retrieval)
    if tuple(staged_score_genes) != retrieval_genes:
        raise AssertionError(f"{held}: staged score-gene order differs from exp064")
    observed = np.asarray(saved_query[rows], dtype=np.float64)
    expected = np.asarray(exact_arm_d[:, retrieval], dtype=np.float64)
    max_abs_error = float(np.max(np.abs(observed - expected)))
    if not np.allclose(observed, expected, rtol=0.0, atol=atol):
        raise AssertionError(
            f"{held}: arm-D differs from exp064 reference by {max_abs_error}"
        )
    extension = tuple(
        target for target in evaluated_targets if target not in set(reference_targets)
    )
    extension_condition_counts = np.asarray(
        [sum(target in bank.pert_pos for bank in cd4_banks) for target in extension],
        dtype=np.int64,
    )
    return {
        "passed": True,
        "comparison": (
            "saved panel_query compared to an independent raw-bank "
            "reconstruction through exp064's original all-mapped-coordinate "
            "normalization, target blocking, whitening, and CD4 composition"
        ),
        "independent_from_staged_panel_preparation": True,
        "normalization_contract": (
            "center and norm on all mapped source coordinates, then block the "
            "300 VCC target coordinates exactly as exp064"
        ),
        "learned_statistics_contract": (
            "this reconstruction is baseline verification only; adapter "
            "training, selection, and permutation nulls use the separate "
            "retained-coordinate-safe train_query, while clustering uses the "
            "separate retained-coordinate-safe panel-row cluster_query"
        ),
        "absolute_tolerance": float(atol),
        "max_abs_error": max_abs_error,
        "validated_targets": len(reference_targets),
        "validated_target_order_sha256": sha256_lines(reference_targets),
        "retrieval_genes": int(len(retrieval)),
        "retrieval_columns_sha256": sha256_array(retrieval),
        "retrieval_gene_order_sha256": sha256_lines(retrieval_genes),
        "saved_query_score_axis_exact": True,
        "source_vote_histogram": {
            str(value): int(np.sum(exact_votes == value))
            for value in sorted(set(exact_votes.tolist()))
        },
        "extension_targets": len(extension),
        "extension_target_order_sha256": sha256_lines(extension),
        "extension_audit": {
            "pm_missing_rows": int(
                sum(
                    target in banks["k562_gwps"].pert_pos
                    and target not in perturbmatch_bank.pert_pos
                    for target in extension
                )
            ),
            "cd4_condition_histogram": {
                str(value): int(np.sum(extension_condition_counts == value))
                for value in sorted(set(extension_condition_counts.tolist()))
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
    master_genes, train_targets, outer_targets, train_mask = validate_truth_contract(
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
    train_blocked = set(panel) | set(train_targets)
    panel_blocked = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    missing_panel_coordinates = sorted(panel_blocked - set(official_pos))
    if missing_panel_coordinates:
        raise ValueError(
            f"{len(missing_panel_coordinates)} VCC targets are absent from the "
            f"official axis; first={missing_panel_coordinates[:5]}"
        )
    missing_nonpanel_coordinates = tuple(sorted(set(train_targets) - set(official_pos)))
    present_train_blocked = train_blocked & set(official_pos)
    train_set = set(train_targets)
    train_score_mask = np.asarray(
        [gene not in train_set for gene in master_genes], dtype=bool
    )
    train_score_genes = tuple(
        gene for gene, keep in zip(master_genes, train_score_mask, strict=True) if keep
    )
    panel_score_genes, panel_score_positions, final_axis_audit = (
        established_final_score_axis(
            official,
            panel,
            args.donor_root,
            manifest,
        )
    )
    correction_genes = tuple(
        gene for gene, keep in zip(master_genes, train_mask, strict=True) if keep
    )
    train_score_positions = np.asarray(
        [official_pos[gene] for gene in train_score_genes], dtype=np.int64
    )
    train_score_local = {gene: index for index, gene in enumerate(train_score_genes)}
    panel_score_local = {gene: index for index, gene in enumerate(panel_score_genes)}
    correction_in_train_score = np.asarray(
        [train_score_local[gene] for gene in correction_genes], dtype=np.int64
    )
    correction_in_panel_score = np.asarray(
        [panel_score_local[gene] for gene in correction_genes], dtype=np.int64
    )
    if set(train_score_genes) & train_blocked:
        raise AssertionError("blocked target coordinate entered training score axis")
    if set(panel_score_genes) & panel_blocked:
        raise AssertionError("VCC coordinate entered final panel score axis")
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
    fold_metadata: dict[str, object] = {}
    fold_specs = {
        held: {
            "sources": (
                "k562_gwps",
                next(name for name in DEVELOPMENT_FOLDS if name != held),
            ),
            "panel_targets": outer_targets,
        }
        for held in DEVELOPMENT_FOLDS
    }
    fold_specs[KOLF_FOLD] = {
        "sources": MAIN_DONORS,
        "panel_targets": panel,
    }
    for held in FINAL_FOLDS:
        source_names = tuple(fold_specs[held]["sources"])
        evaluation_targets = tuple(fold_specs[held]["panel_targets"])
        banks = {
            name: EXP062.Bank.load(name, args.donor_root / name) for name in MAIN_DONORS
        }
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
            set(panel),
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
        if not set(map(int, train_score_positions)) <= set(
            map(int, train_active_columns)
        ):
            raise ValueError(
                f"{held}: training score axis is outside restricted operator"
            )
        if not set(map(int, panel_score_positions)) <= set(map(int, operator_columns)):
            raise ValueError(f"{held}: panel score axis is outside operator")

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
            train_score_positions,
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
        train_audit["support"] = train_support

        (
            cluster_prepared,
            _cluster_prepared_pm,
            cluster_pm_slot,
            cluster_cd4,
            _cluster_shared_columns,
            cluster_support,
        ) = prepare_arm_inputs(
            banks,
            perturbmatch_bank,
            cd4_banks,
            source_names,
            evaluation_targets,
            official,
            train_blocked,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        cluster_query, _cluster_full, cluster_audit = build_arm_d_panel(
            cluster_pm_slot,
            [cluster_prepared[name] for name in source_names if name != "k562_gwps"],
            cluster_cd4,
            operator,
            operator_columns,
            train_active_columns,
            np.arange(len(evaluation_targets), dtype=np.int64),
            train_score_positions,
        )
        cluster_audit["operator_restriction"][
            "all_vcc_and_historical_target_coordinates_excluded"
        ] = True
        cluster_audit["support"] = {
            **cluster_support,
            "used_for_adapter_fit_or_selection": False,
            "used_for_source_clustering": True,
            "used_for_promotion_inference": True,
        }

        (
            panel_prepared,
            panel_prepared_pm,
            panel_pm_slot,
            panel_cd4,
            panel_shared_columns,
            panel_support,
        ) = prepare_exp064_panel_inputs(
            banks,
            perturbmatch_bank,
            cd4_banks,
            source_names,
            evaluation_targets,
            official,
            panel,
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
        panel_audit["support"] = panel_support

        fold_dir = args.output / held
        fold_dir.mkdir(parents=True, exist_ok=True)
        train_path = fold_dir / "train_query.npy"
        cluster_path = fold_dir / "cluster_query.npy"
        panel_path = fold_dir / "panel_query.npy"
        np.save(train_path, train_query, allow_pickle=False)
        np.save(cluster_path, cluster_query, allow_pickle=False)
        np.save(panel_path, panel_query, allow_pickle=False)
        train_audit["query_sha256"] = sha256_file(train_path)
        cluster_audit["query_sha256"] = sha256_file(cluster_path)
        cluster_audit["target_order_sha256"] = sha256_lines(evaluation_targets)
        cluster_audit["score_gene_count"] = len(train_score_genes)
        cluster_audit["score_gene_order_sha256"] = sha256_lines(train_score_genes)
        cluster_audit["score_official_positions_sha256"] = sha256_array(
            train_score_positions
        )
        cluster_audit["saved_query_is_complete_cluster_axis"] = True
        panel_audit["query_sha256"] = sha256_file(panel_path)
        panel_audit["target_order_sha256"] = sha256_lines(evaluation_targets)
        panel_audit["score_axis_excludes_exactly_vcc_coordinates"] = True
        panel_audit["score_gene_count"] = len(panel_score_genes)
        panel_audit["score_gene_order_sha256"] = sha256_lines(panel_score_genes)
        panel_audit["score_official_positions_sha256"] = sha256_array(
            panel_score_positions
        )
        panel_audit["saved_query_is_complete_scored_axis"] = True

        reproduction = None
        if held in DEVELOPMENT_FOLDS:
            saved_panel_query = np.load(panel_path, mmap_mode="r", allow_pickle=False)
            reproduction = exp064_reproduction_audit(
                held=held,
                reference_targets=exp064_targets,
                evaluated_targets=evaluation_targets,
                banks=banks,
                perturbmatch_bank=perturbmatch_bank,
                cd4_banks=cd4_banks,
                operator=operator,
                operator_columns=operator_columns,
                saved_query=saved_panel_query,
                staged_score_positions=panel_score_positions,
                staged_score_genes=panel_score_genes,
                official=official,
                panel=panel,
                chunk_rows=args.chunk_rows,
                max_rss_gib=args.max_rss_gib,
                atol=args.exp064_atol,
            )
            del saved_panel_query

        fold_metadata[held] = {
            "held_donor": held,
            "source_donors": list(source_names),
            "held_donor_used_in_prediction": False,
            "held_donor_used_in_operator": False,
            "operator_npz_sha256": operator_entry["npz_sha256"],
            "operator_json_sha256": operator_entry["json_sha256"],
            "operator_columns": int(len(operator_columns)),
            "operator_columns_sha256": sha256_array(operator_columns),
            "panel_targets": len(evaluation_targets),
            "panels": {
                "train": train_audit,
                "cluster": cluster_audit,
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
            cluster_prepared,
            _cluster_prepared_pm,
            cluster_pm_slot,
            cluster_cd4,
            panel_prepared,
            panel_prepared_pm,
            panel_pm_slot,
            panel_cd4,
            operator,
            train_query,
            cluster_query,
            _cluster_full,
            panel_query,
            panel_full,
        )
        gc.collect()
        enforce_memory_limit(args.max_rss_gib, f"{held} completion")

    metadata: dict[str, object] = {
        "schema_version": CACHE_SCHEMA,
        "experiment": "080-listwise-metric-adapter",
        "status": "arm_d_query_cache_complete",
        "arm": "D",
        "arm_definition": (
            "50:50 incumbent/PerturbMatch K562 slot, source-only gamma=0.65 "
            "union whitening, then CD4 mean3 gene_norm residual at fraction 0.05"
        ),
        "baseline_reproduction_contract": (
            "panel_query uses exp064's original all-mapped-coordinate source "
            "normalization followed by VCC-coordinate blocking; each development "
            "fold is independently reconstructed from raw banks at <=1e-6"
        ),
        "learned_feature_contract": (
            "train_query and cluster_query centering, norm calibration, and "
            "whitening restriction exclude all VCC and historical target "
            "coordinates; only train_query enters adapter fitting/selection and "
            "only cluster_query enters source clustering"
        ),
        "development_folds": list(DEVELOPMENT_FOLDS),
        "final_folds": list(FINAL_FOLDS),
        "invalid_outer_donor": (
            "k562_gwps omitted because PerturbMatch is a reanalysis of the same assay"
        ),
        "training_targets": len(train_targets),
        "outer_targets": len(outer_targets),
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
        "panel_blocked_target_coordinates": len(panel_blocked),
        "panel_blocked_target_coordinates_sha256": sha256_lines(sorted(panel_blocked)),
        "historical_target_coordinates_restored_for_panel": len(
            historical_panel_coordinates
        ),
        "train_score_genes": len(train_score_genes),
        "cluster_score_genes": len(train_score_genes),
        "panel_score_genes": len(panel_score_genes),
        "correction_genes": len(correction_genes),
        "train_score_gene_order_sha256": sha256_lines(train_score_genes),
        "cluster_score_gene_order_sha256": sha256_lines(train_score_genes),
        "cluster_score_official_positions_sha256": sha256_array(train_score_positions),
        "panel_score_gene_order_sha256": sha256_lines(panel_score_genes),
        "panel_score_official_positions_sha256": sha256_array(panel_score_positions),
        "correction_gene_order_sha256": sha256_lines(correction_genes),
        "correction_in_train_score_sha256": sha256_array(correction_in_train_score),
        "correction_in_panel_score_sha256": sha256_array(correction_in_panel_score),
        "train_order_sha256": sha256_lines(train_targets),
        "outer_order_sha256": sha256_lines(outer_targets),
        "panel_order_sha256": sha256_lines(panel),
        "exp064_target_order_sha256": sha256_lines(exp064_targets),
        "folds": fold_metadata,
        "final_score_axis": final_axis_audit,
        "leakage_audit": {
            "all_vcc_target_coordinates_removed": True,
            "all_nonpanel_target_coordinates_removed_from_training": True,
            "nonpanel_target_coordinates_present_in_final_panel_scoring": True,
            "training_dense_whitener_principally_restricted": True,
            "final_panel_excludes_exactly_the_vcc_target_set": True,
            "coordinate_less_nonpanel_targets_treated_as_no_coordinate_to_remove": True,
            "saved_panel_query_uses_complete_17806_gene_score_axis": True,
            "cluster_queries_exclude_all_vcc_and_historical_target_coordinates": True,
            "exp064_baseline_reproduction_is_inference_only": True,
            "learned_statistics_use_only_retained_coordinates": True,
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
    if int(manifest.get("schema_version", -1)) != 2:
        raise SystemExit("unsupported artifact manifest")
    verify_inputs(args, manifest)
    metadata = stage(args, manifest)
    print(json.dumps(metadata["leakage_audit"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
