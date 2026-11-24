#!/usr/bin/env python3
"""069 — fit whitening covariance to the frozen PerturbMatch K562 slot.

The comparator is experiment 064 arm D:

* 50:50 incumbent/PerturbMatch K562 slot, with incumbent fallback;
* the existing candidate-independent gamma-0.65 union whitener;
* the frozen post-whitening CD4 mean3 gene-normalized residual at fraction 0.05.

The candidate changes exactly one scientific object: its whitening covariance
uses the same frozen K562 processing ensemble as the prediction, plus the
non-held Orion donor.  The two covariance sources have equal total mass.  All
300 challenge perturbation rows and target-gene coordinates are excluded from
the candidate covariance fit.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import importlib.util
import json
import math
import os
import resource
import shutil
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXP064_PATH = ROOT / "experiments/064-pm-cd4-composition/run.py"
SPEC = importlib.util.spec_from_file_location("exp069_exp064", EXP064_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(f"cannot import frozen experiment-064 code: {EXP064_PATH}")
EXP064 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = EXP064
SPEC.loader.exec_module(EXP064)
BASE = EXP064.BASE

FOLDS = ("orion_hct116", "orion_hek293t")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
PM_FRACTION = 0.50
CD4_MODE = "gene_norm"
CD4_FRACTION = 0.05
WHITEN_GAMMA = 0.65
BOOTSTRAP_REPS = 20_000
BOOTSTRAP_SEED = 6_902_026
MACRO_GATE = 0.003
EXPECTED_BASELINE_WHITENER_SHA256 = dict(EXP064.EXPECTED_WHITENER_SHA256)


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB >= {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def array_sha256(values: np.ndarray) -> str:
    canonical = np.ascontiguousarray(values)
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def panel_background_rows(bank: BASE.Bank, panel: set[str]) -> np.ndarray:
    rows = np.asarray(
        [
            index
            for index, perturbation in enumerate(bank.perts)
            if perturbation not in panel
        ],
        dtype=np.int64,
    )
    if not len(rows):
        raise ValueError(f"{bank.name}: no non-panel covariance rows")
    if any(bank.perts[int(index)] in panel for index in rows):
        raise AssertionError(f"{bank.name}: panel row entered covariance basis")
    return rows


def build_pm_slot_with_fallback(
    incumbent: BASE.PreparedSource,
    perturbmatch: BASE.PreparedSource,
    shared_columns: np.ndarray,
) -> tuple[BASE.PreparedSource, dict[str, object]]:
    """Apply the frozen PM blend only where both K562 estimates are available."""
    if incumbent.effects.shape != perturbmatch.effects.shape:
        raise ValueError("incumbent and PerturbMatch target matrices differ in shape")
    columns = np.asarray(shared_columns, dtype=np.int64)
    rows = incumbent.target_available & perturbmatch.target_available
    effects = np.asarray(incumbent.effects, dtype=np.float64).copy()
    if rows.any() and len(columns):
        effects[np.ix_(rows, columns)] = (1.0 - PM_FRACTION) * incumbent.effects[
            np.ix_(rows, columns)
        ] + PM_FRACTION * perturbmatch.effects[np.ix_(rows, columns)]
    return (
        BASE.PreparedSource(
            name="k562_gwps_pm50_fallback",
            effects=effects,
            target_available=np.asarray(incumbent.target_available, dtype=bool).copy(),
            mean_on_axis=np.asarray(incumbent.mean_on_axis, dtype=np.float64).copy(),
            raw_median_norm=incumbent.raw_median_norm,
            mapped_gene_count=incumbent.mapped_gene_count,
            duplicate_symbols=incumbent.duplicate_symbols,
        ),
        {
            "targets_total": int(len(rows)),
            "targets_blended": int(rows.sum()),
            "targets_incumbent_fallback": int((~rows).sum()),
            "shared_non_panel_genes": int(len(columns)),
        },
    )


def freeze_complete_cd4_mean3(cd4: BASE.CD4Prepared) -> BASE.CD4Prepared:
    """Retain the residual only for targets represented in all three conditions."""
    complete = np.asarray(cd4.target_condition_count) == len(CD4_NAMES)
    effects = np.asarray(cd4.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(cd4.gene_votes, dtype=np.uint8).copy()
    effects[~complete] = 0.0
    votes[~complete] = 0
    return BASE.CD4Prepared(
        mode=cd4.mode,
        mean_effect=effects,
        gene_votes=votes,
        target_condition_count=np.asarray(
            cd4.target_condition_count, dtype=np.uint8
        ).copy(),
        condition_details=list(cd4.condition_details),
    )


@dataclass(frozen=True)
class BackgroundStats:
    rows: np.ndarray
    source_columns: np.ndarray
    operator_columns: np.ndarray
    mean: np.ndarray
    median_norm: float


@dataclass(frozen=True)
class EnsembleStats:
    mean: np.ndarray
    median_norm: float
    matched_rows: int


def background_stats(
    bank: BASE.Bank,
    official: Sequence[str],
    wcols: np.ndarray,
    panel: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> BackgroundStats:
    """Center and scale one source using non-panel rows only."""
    rows = panel_background_rows(bank, panel)
    source_columns, operator_columns = BASE.last_source_columns(bank, official, wcols)
    if not len(source_columns):
        raise ValueError(f"{bank.name}: no genes map to the operator")
    delta = np.load(bank.delta_path, mmap_mode="r")
    mean = np.zeros(len(source_columns), dtype=np.float64)
    for lo in range(0, len(rows), chunk_rows):
        rr = rows[lo : lo + chunk_rows]
        mean += np.asarray(delta[np.ix_(rr, source_columns)], dtype=np.float64).sum(
            axis=0, dtype=np.float64
        )
        enforce_limit(max_rss_gib, f"{bank.name} non-panel center")
    mean /= len(rows)
    norms = np.empty(len(rows), dtype=np.float64)
    for lo in range(0, len(rows), chunk_rows):
        rr = rows[lo : lo + chunk_rows]
        block = np.asarray(delta[np.ix_(rr, source_columns)], dtype=np.float64)
        block -= mean
        norms[lo : lo + len(rr)] = np.linalg.norm(block, axis=1)
        enforce_limit(max_rss_gib, f"{bank.name} non-panel scale")
    median_norm = float(np.median(norms)) or 1.0
    return BackgroundStats(
        rows=rows,
        source_columns=source_columns,
        operator_columns=operator_columns,
        mean=mean,
        median_norm=median_norm,
    )


def source_mass_row_weight(n_rows: int) -> float:
    if n_rows <= 0:
        raise ValueError("source must have at least one covariance row")
    return math.sqrt(0.5 / n_rows)


def mapped_source_columns(
    bank: BASE.Bank,
    official: Sequence[str],
    wcols: np.ndarray,
) -> dict[int, tuple[int, int]]:
    """Map official gene position to (source column, operator column)."""
    source_columns, operator_columns = BASE.last_source_columns(bank, official, wcols)
    return {
        int(wcols[int(operator_column)]): (int(source_column), int(operator_column))
        for source_column, operator_column in zip(
            source_columns, operator_columns, strict=True
        )
    }


def normalized_source_block(
    bank: BASE.Bank,
    stats: BackgroundStats,
    rows: np.ndarray,
    n_operator_genes: int,
) -> np.ndarray:
    delta = np.load(bank.delta_path, mmap_mode="r")
    values = np.asarray(delta[np.ix_(rows, stats.source_columns)], dtype=np.float64)
    values -= stats.mean
    values *= BASE.EFFECT_SCALE / stats.median_norm
    block = np.zeros((len(rows), n_operator_genes), dtype=np.float64)
    block[:, stats.operator_columns] = values
    return block


def blend_background_block(
    incumbent: BASE.Bank,
    perturbmatch: BASE.Bank,
    incumbent_stats: BackgroundStats,
    perturbmatch_stats: BackgroundStats,
    rows: np.ndarray,
    official: Sequence[str],
    wcols: np.ndarray,
    panel: set[str],
) -> tuple[np.ndarray, int]:
    """Build K562 ensemble rows with target and gene fallback to incumbent."""
    block = normalized_source_block(incumbent, incumbent_stats, rows, len(wcols))
    incumbent_map = mapped_source_columns(incumbent, official, wcols)
    pm_map = mapped_source_columns(perturbmatch, official, wcols)
    common_official = sorted(set(incumbent_map) & set(pm_map))
    if not common_official:
        raise ValueError("incumbent and PerturbMatch have no common operator genes")
    common_operator = np.asarray(
        [incumbent_map[position][1] for position in common_official],
        dtype=np.int64,
    )
    pm_source = np.asarray(
        [pm_map[position][0] for position in common_official], dtype=np.int64
    )
    local_positions: list[int] = []
    pm_rows: list[int] = []
    for local, incumbent_row in enumerate(rows):
        target = incumbent.perts[int(incumbent_row)]
        if target in panel:
            raise AssertionError("panel perturbation entered K562 ensemble basis")
        pm_row = perturbmatch.pert_pos.get(target)
        if pm_row is not None and perturbmatch.perts[pm_row] not in panel:
            local_positions.append(local)
            pm_rows.append(pm_row)
    if local_positions:
        pm_delta = np.load(perturbmatch.delta_path, mmap_mode="r")
        pm_values = np.asarray(
            pm_delta[
                np.ix_(
                    np.asarray(pm_rows, dtype=np.int64),
                    pm_source,
                )
            ],
            dtype=np.float64,
        )
        pm_mean_by_source = {
            int(source): float(mean)
            for source, mean in zip(
                perturbmatch_stats.source_columns,
                perturbmatch_stats.mean,
                strict=True,
            )
        }
        pm_mean = np.asarray(
            [pm_mean_by_source[int(source)] for source in pm_source],
            dtype=np.float64,
        )
        pm_values -= pm_mean
        pm_values *= BASE.EFFECT_SCALE / perturbmatch_stats.median_norm
        local_array = np.asarray(local_positions, dtype=np.int64)
        block[np.ix_(local_array, common_operator)] = (1.0 - PM_FRACTION) * block[
            np.ix_(local_array, common_operator)
        ] + PM_FRACTION * pm_values
    return block, len(local_positions)


def ensemble_background_stats(
    incumbent: BASE.Bank,
    perturbmatch: BASE.Bank,
    incumbent_stats: BackgroundStats,
    perturbmatch_stats: BackgroundStats,
    official: Sequence[str],
    wcols: np.ndarray,
    panel: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> EnsembleStats:
    """Center and scale the blended K562 source on its frozen row universe."""
    rows = incumbent_stats.rows
    mean = np.zeros(len(wcols), dtype=np.float64)
    matched_rows = 0
    for lo in range(0, len(rows), chunk_rows):
        rr = rows[lo : lo + chunk_rows]
        block, matched = blend_background_block(
            incumbent,
            perturbmatch,
            incumbent_stats,
            perturbmatch_stats,
            rr,
            official,
            wcols,
            panel,
        )
        mean += block.sum(axis=0, dtype=np.float64)
        matched_rows += matched
        enforce_limit(max_rss_gib, "K562 ensemble center")
    mean /= len(rows)

    norms = np.empty(len(rows), dtype=np.float64)
    for lo in range(0, len(rows), chunk_rows):
        rr = rows[lo : lo + chunk_rows]
        block, _ = blend_background_block(
            incumbent,
            perturbmatch,
            incumbent_stats,
            perturbmatch_stats,
            rr,
            official,
            wcols,
            panel,
        )
        block -= mean
        norms[lo : lo + len(rr)] = np.linalg.norm(block, axis=1)
        enforce_limit(max_rss_gib, "K562 ensemble scale")
    return EnsembleStats(
        mean=mean,
        median_norm=float(np.median(norms)) or 1.0,
        matched_rows=matched_rows,
    )


def fit_aligned_whitener(
    incumbent: BASE.Bank,
    perturbmatch: BASE.Bank,
    companion: BASE.Bank,
    official: Sequence[str],
    panel: Sequence[str],
    output_path: Path,
    scratch: Path,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[BASE.Operator, np.ndarray]:
    """Fit equal-source-mass covariance from ensemble K562 and companion."""
    panel_set = set(panel)
    wcols = EXP064.operator_columns((incumbent, companion), official, panel_set)
    incumbent_stats = background_stats(
        incumbent,
        official,
        wcols,
        panel_set,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    perturbmatch_stats = background_stats(
        perturbmatch,
        official,
        wcols,
        panel_set,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    companion_stats = background_stats(
        companion,
        official,
        wcols,
        panel_set,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    ensemble_stats = ensemble_background_stats(
        incumbent,
        perturbmatch,
        incumbent_stats,
        perturbmatch_stats,
        official,
        wcols,
        panel_set,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    incumbent_only = len(incumbent_stats.rows) - ensemble_stats.matched_rows
    metadata: dict[str, object] = {
        "axis": "union",
        "candidate": "perturbmatch_aligned_whitener",
        "donors": ["k562_gwps_pm50_fallback", companion.name],
        "n_cols": int(len(wcols)),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(panel_set),
        "excluded_panel_rows_by_source": {
            incumbent.name: int(incumbent.shape[0] - len(incumbent_stats.rows)),
            perturbmatch.name: int(
                perturbmatch.shape[0] - len(perturbmatch_stats.rows)
            ),
            companion.name: int(companion.shape[0] - len(companion_stats.rows)),
        },
        "basis_rows_by_source": {
            "k562_gwps_pm50_fallback": int(len(incumbent_stats.rows)),
            companion.name: int(len(companion_stats.rows)),
        },
        "k562_background_rows_blended": int(ensemble_stats.matched_rows),
        "k562_background_rows_incumbent_fallback": int(incumbent_only),
        "k562_ensemble_median_norm_before_rescale": ensemble_stats.median_norm,
        "source_mass": {
            "k562_gwps_pm50_fallback": 0.5,
            companion.name: 0.5,
        },
        "basis_normalization": (
            "K562 views independently centered and median-norm scaled on "
            "non-panel rows, blended 50:50 with fallback, then the ensemble "
            "is centered and median-norm rescaled; Orion is independently "
            "centered and median-norm scaled"
        ),
        "perturbmatch_fraction": PM_FRACTION,
        "perturbmatch_row_universe": (
            "incumbent K562 non-panel rows; incumbent fallback for missing "
            "PerturbMatch targets and genes"
        ),
        "panel_rows_excluded": True,
        "panel_coordinates_excluded": True,
        "gamma": WHITEN_GAMMA,
    }
    metadata_path = output_path.with_suffix(".json")
    if output_path.exists() or metadata_path.exists():
        if not output_path.exists() or not metadata_path.exists():
            raise ValueError(f"{output_path}: incomplete cached operator")
        actual = json.loads(metadata_path.read_text())
        if actual != metadata:
            raise ValueError(f"{output_path}: stale candidate whitener metadata")
        with np.load(output_path) as archive:
            evals = np.asarray(archive["evals"], dtype=np.float64)
            evecs = np.asarray(archive["evecs"], dtype=np.float32)
        if evecs.shape != (len(wcols), len(wcols)):
            raise ValueError(f"{output_path}: cached operator shape mismatch")
        return BASE.Operator(evals, evecs, actual, output_path, metadata_path), wcols

    output_path.parent.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    total_rows = len(incumbent_stats.rows) + len(companion_stats.rows)
    basis_path = scratch / "equal_source_mass_basis.npy"
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float64,
        shape=(total_rows, len(wcols)),
    )
    k_weight = source_mass_row_weight(len(incumbent_stats.rows))
    matched_rows = 0
    cursor = 0
    for lo in range(0, len(incumbent_stats.rows), chunk_rows):
        rr = incumbent_stats.rows[lo : lo + chunk_rows]
        block, matched = blend_background_block(
            incumbent,
            perturbmatch,
            incumbent_stats,
            perturbmatch_stats,
            rr,
            official,
            wcols,
            panel_set,
        )
        block -= ensemble_stats.mean
        block *= BASE.EFFECT_SCALE / ensemble_stats.median_norm
        basis[cursor : cursor + len(rr)] = block * k_weight
        cursor += len(rr)
        matched_rows += matched
        enforce_limit(max_rss_gib, "K562 ensemble covariance basis")
    if matched_rows != ensemble_stats.matched_rows:
        raise AssertionError(
            f"matched K562 rows {matched_rows} != {ensemble_stats.matched_rows}"
        )

    companion_weight = source_mass_row_weight(len(companion_stats.rows))
    for lo in range(0, len(companion_stats.rows), chunk_rows):
        rr = companion_stats.rows[lo : lo + chunk_rows]
        block = normalized_source_block(companion, companion_stats, rr, len(wcols))
        basis[cursor : cursor + len(rr)] = block * companion_weight
        cursor += len(rr)
        enforce_limit(max_rss_gib, f"{companion.name} covariance basis")
    if cursor != total_rows:
        raise AssertionError("candidate basis row count mismatch")
    basis.flush()

    print(
        f"fit aligned whitener companion={companion.name} "
        f"basis={basis.shape} matched_k562_rows={matched_rows} "
        f"peak_rss={peak_rss_gib():.3f} GiB",
        flush=True,
    )
    covariance = basis.T @ basis
    del basis
    gc.collect()
    evals, evecs = eigh(covariance, overwrite_a=True, check_finite=False, driver="evr")
    del covariance
    gc.collect()
    evals = np.maximum(
        np.asarray(evals, dtype=np.float64),
        float(evals[-1]) * 1e-12,
    )
    evecs32 = np.asarray(evecs, dtype=np.float32)
    del evecs
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.savez(handle, evals=evals, evecs=evecs32)
    temporary.replace(output_path)
    metadata_temporary = metadata_path.with_suffix(metadata_path.suffix + ".tmp")
    metadata_temporary.write_text(json.dumps(metadata, indent=2) + "\n")
    metadata_temporary.replace(metadata_path)
    try:
        basis_path.unlink()
    except FileNotFoundError:
        pass
    enforce_limit(max_rss_gib, f"{companion.name} aligned whitener fit")
    return (
        BASE.Operator(evals, evecs32, metadata, output_path, metadata_path),
        wcols,
    )


def paired_target_cluster_bootstrap(
    per_fold_delta: dict[str, np.ndarray],
    *,
    reps: int,
    seed: int,
) -> dict[str, object]:
    arrays = [np.asarray(per_fold_delta[held], dtype=np.float64) for held in FOLDS]
    if len({len(values) for values in arrays}) != 1:
        raise ValueError("held folds do not share the frozen full gallery")
    cluster_delta = np.mean(np.stack(arrays, axis=1), axis=1)
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    for lo in range(0, reps, 1_000):
        hi = min(lo + 1_000, reps)
        indices = rng.integers(
            0, len(cluster_delta), size=(hi - lo, len(cluster_delta))
        )
        draws[lo:hi] = cluster_delta[indices].mean(axis=1)
    return {
        "method": (
            "paired target-cluster bootstrap over the common full gallery; "
            "one target resample preserves baseline/candidate pairing and both folds"
        ),
        "reps": int(reps),
        "seed": int(seed),
        "estimate": float(cluster_delta.mean()),
        "ci95_percentile": [
            float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)),
        ],
        "positive_probability": float(np.mean(draws > 0)),
        "targets": int(len(cluster_delta)),
    }


def stopping_gate(
    fold_deltas: dict[str, float],
    macro_delta: float,
    bootstrap: dict[str, object],
) -> dict[str, object]:
    lower = float(bootstrap["ci95_percentile"][0])
    checks = {
        "both_folds_nonnegative": all(
            float(fold_deltas[held]) >= 0.0 for held in FOLDS
        ),
        "macro_gain_at_least_0_003": float(macro_delta) >= MACRO_GATE,
        "bootstrap_lower_bound_positive": lower > 0.0,
    }
    return {
        "thresholds": {
            "minimum_fold_delta": 0.0,
            "minimum_macro_delta": MACRO_GATE,
            "minimum_bootstrap_lower_bound": 0.0,
        },
        "checks": checks,
        "passes": all(checks.values()),
        "action": (
            "eligible for one count-level confirmation"
            if all(checks.values())
            else "stop; do not run count-level confirmation"
        ),
    }


def evaluate_fold(
    held: str,
    banks: dict[str, BASE.Bank],
    perturbmatch: BASE.Bank,
    prepared: dict[str, BASE.PreparedSource],
    pm_slot: BASE.PreparedSource,
    cd4: BASE.CD4Prepared,
    official: Sequence[str],
    panel: Sequence[str],
    targets: Sequence[str],
    baseline_whitener_path: Path,
    candidate_whitener_path: Path,
    scratch: Path,
    *,
    expected_baseline_sha256: str,
    chunk_rows: int,
    max_rss_gib: float,
    min_operator_genes: int,
    min_retrieval_genes: int,
) -> tuple[dict[str, object], list[dict[str, object]], np.ndarray]:
    companion_name = next(name for name in FOLDS if name != held)
    source_banks = (banks["k562_gwps"], banks[companion_name])
    baseline_operator, baseline_wcols = EXP064.load_frozen_operator(
        baseline_whitener_path,
        source_banks,
        official,
        set(panel),
        expected_sha256=expected_baseline_sha256,
    )
    candidate_operator, candidate_wcols = fit_aligned_whitener(
        banks["k562_gwps"],
        perturbmatch,
        banks[companion_name],
        official,
        panel,
        candidate_whitener_path,
        scratch,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    if not np.array_equal(baseline_wcols, candidate_wcols):
        raise AssertionError("baseline and candidate operator axes differ")
    wcols = candidate_wcols
    if len(wcols) < min_operator_genes:
        raise ValueError(f"{held}: only {len(wcols)} operator genes")

    arm_d_pre, source_votes = BASE.build_baseline([pm_slot, prepared[companion_name]])
    baseline_post, baseline_preserve = BASE.apply_operator(
        arm_d_pre, baseline_operator, wcols
    )
    candidate_post, candidate_preserve = BASE.apply_operator(
        arm_d_pre, candidate_operator, wcols
    )
    baseline_d, baseline_cd4_preserve = BASE.add_post_whiten_residual(
        baseline_post, cd4, CD4_FRACTION, wcols
    )
    candidate_d, candidate_cd4_preserve = BASE.add_post_whiten_residual(
        candidate_post, cd4, CD4_FRACTION, wcols
    )

    held_genes = set(banks[held].genes)
    retrieval = np.asarray(
        [index for index in wcols if official[index] in held_genes],
        dtype=np.int64,
    )
    if len(retrieval) < min_retrieval_genes:
        raise ValueError(f"{held}: only {len(retrieval)} retrieval genes")
    if set(panel) & {official[index] for index in retrieval}:
        raise AssertionError("panel target coordinate entered retrieval")
    truth = BASE.raw_truth(prepared[held])
    baseline_per_target = BASE.pds_per_target(baseline_d, truth, retrieval)
    candidate_per_target = BASE.pds_per_target(candidate_d, truth, retrieval)
    paired_delta = candidate_per_target - baseline_per_target
    baseline_pds = float(np.mean(baseline_per_target))
    candidate_pds = float(np.mean(candidate_per_target))
    delta = candidate_pds - baseline_pds
    per_target_rows = [
        {
            "held": held,
            "target_gene": target,
            "baseline_arm_d_pds": float(base),
            "aligned_whitener_arm_d_pds": float(candidate),
            "delta": float(change),
        }
        for target, base, candidate, change in zip(
            targets,
            baseline_per_target,
            candidate_per_target,
            paired_delta,
            strict=True,
        )
    ]
    fold = {
        "held": held,
        "companion_source": companion_name,
        "gallery_targets": len(targets),
        "baseline_arm_d_pds": baseline_pds,
        "aligned_whitener_arm_d_pds": candidate_pds,
        "delta": delta,
        "improved_targets": int(np.sum(paired_delta > 0)),
        "worsened_targets": int(np.sum(paired_delta < 0)),
        "tied_targets": int(np.sum(paired_delta == 0)),
        "operator_genes": len(wcols),
        "retrieval_genes": len(retrieval),
        "operator_columns_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "retrieval_columns_sha256": hashlib.sha256(retrieval.tobytes()).hexdigest(),
        "source_vote_counts": {
            str(value): int(np.sum(source_votes == value))
            for value in sorted(set(source_votes.tolist()))
        },
        "baseline_whitener": {
            "path": str(baseline_whitener_path.resolve()),
            "sha256": sha256_file(baseline_whitener_path),
            "metadata_sha256": sha256_file(baseline_whitener_path.with_suffix(".json")),
            "preserve": baseline_preserve,
            "cd4_preserve": baseline_cd4_preserve,
        },
        "candidate_whitener": {
            "path": str(candidate_whitener_path.resolve()),
            "sha256": sha256_file(candidate_whitener_path),
            "metadata_sha256": sha256_file(
                candidate_whitener_path.with_suffix(".json")
            ),
            "metadata": candidate_operator.metadata,
            "preserve": candidate_preserve,
            "cd4_preserve": candidate_cd4_preserve,
        },
    }
    del (
        baseline_operator,
        candidate_operator,
        arm_d_pre,
        baseline_post,
        candidate_post,
        baseline_d,
        candidate_d,
        truth,
        baseline_per_target,
        candidate_per_target,
    )
    gc.collect()
    shutil.rmtree(scratch, ignore_errors=True)
    enforce_limit(max_rss_gib, f"{held} fold completion")
    return fold, per_target_rows, paired_delta


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    flattened = [
        {
            key: json.dumps(value, sort_keys=True)
            if isinstance(value, (dict, list))
            else value
            for key, value in row.items()
        }
        for row in rows
    ]
    fields = sorted({key for row in flattened for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flattened)


def bank_hashes(bank: BASE.Bank, *, include_b_ctrl: bool = True) -> dict[str, str]:
    names = ["delta.npy", "genes.txt", "perts.txt", "meta.json"]
    if include_b_ctrl:
        names.insert(1, "b_ctrl.npy")
    return {
        name: sha256_file(bank.path / name)
        for name in names
        if (bank.path / name).exists()
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = BASE.read_csv_column(args.panel_csv, "target_gene")
    official = BASE.read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != args.expected_panel_targets:
        raise ValueError(
            f"panel has {len(panel)} targets, expected {args.expected_panel_targets}"
        )
    official_pos = {gene: index for index, gene in enumerate(official)}
    if any(target not in official_pos for target in panel):
        raise ValueError("a panel target is absent from the official gene axis")
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    panel_set = set(panel)

    banks = {
        "k562_gwps": BASE.Bank.load("k562_gwps", args.k562),
        "orion_hct116": BASE.Bank.load("orion_hct116", args.hct116),
        "orion_hek293t": BASE.Bank.load("orion_hek293t", args.hek293t),
    }
    perturbmatch = EXP064.load_perturbmatch(args.perturbmatch)
    cd4_banks = tuple(
        BASE.Bank.load(name, args.cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    targets = tuple(
        target
        for target in panel
        if target in banks["orion_hct116"].pert_pos
        and target in banks["orion_hek293t"].pert_pos
    )
    if len(targets) != args.expected_full_gallery_targets:
        raise ValueError(
            f"full common held gallery {len(targets)} != "
            f"{args.expected_full_gallery_targets}"
        )

    prepared = {
        name: BASE.prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name, bank in banks.items()
    }
    prepared_pm = BASE.prepare_main_source(
        perturbmatch,
        targets,
        official,
        panel_positions,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    shared_columns = EXP064.shared_k562_pm_columns(
        banks["k562_gwps"], perturbmatch, official, panel_set
    )
    if len(shared_columns) < args.min_pm_shared_genes:
        raise ValueError(f"only {len(shared_columns)} shared K562/PM genes")
    pm_slot, pm_slot_audit = build_pm_slot_with_fallback(
        prepared["k562_gwps"], prepared_pm, shared_columns
    )
    cd4 = freeze_complete_cd4_mean3(
        BASE.prepare_cd4_mean(
            cd4_banks,
            targets,
            official,
            panel_set,
            CD4_MODE,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
    )
    if np.any(cd4.mean_effect[:, panel_positions] != 0):
        raise AssertionError("CD4 residual contains a panel target coordinate")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_whitener_dir.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    folds: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    per_fold_delta: dict[str, np.ndarray] = {}
    for held in FOLDS:
        expected_hash = EXPECTED_BASELINE_WHITENER_SHA256[held]
        if held == "orion_hct116" and args.expected_hct_whitener_sha256:
            expected_hash = args.expected_hct_whitener_sha256
        if held == "orion_hek293t" and args.expected_hek_whitener_sha256:
            expected_hash = args.expected_hek_whitener_sha256
        fold, rows, paired = evaluate_fold(
            held,
            banks,
            perturbmatch,
            prepared,
            pm_slot,
            cd4,
            official,
            panel,
            targets,
            args.baseline_whitener_dir / f"held_{held}.npz",
            args.candidate_whitener_dir / f"held_{held}.npz",
            args.scratch / f"held_{held}",
            expected_baseline_sha256=expected_hash,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
            min_operator_genes=args.min_operator_genes,
            min_retrieval_genes=args.min_retrieval_genes,
        )
        folds.append(fold)
        per_target_rows.extend(rows)
        per_fold_delta[held] = paired
        print(
            f"held={held} baseline_D={fold['baseline_arm_d_pds']:.9f} "
            f"aligned_D={fold['aligned_whitener_arm_d_pds']:.9f} "
            f"delta={fold['delta']:+.9f} "
            f"peak_rss={peak_rss_gib():.3f} GiB",
            flush=True,
        )

    baseline_macro = float(np.mean([fold["baseline_arm_d_pds"] for fold in folds]))
    candidate_macro = float(
        np.mean([fold["aligned_whitener_arm_d_pds"] for fold in folds])
    )
    macro_delta = candidate_macro - baseline_macro
    bootstrap = paired_target_cluster_bootstrap(
        per_fold_delta, reps=args.bootstrap_reps, seed=args.bootstrap_seed
    )
    if not np.isclose(bootstrap["estimate"], macro_delta, atol=1e-14):
        raise AssertionError("bootstrap estimate does not match macro fold delta")
    fold_deltas = {str(fold["held"]): float(fold["delta"]) for fold in folds}
    gate = stopping_gate(fold_deltas, macro_delta, bootstrap)

    result: dict[str, object] = {
        "experiment": "069-perturbmatch-aligned-whitener",
        "status": "completed",
        "protocol": {
            "metric": "effect-level rank PDS",
            "design": "one frozen candidate versus existing experiment-064 arm D",
            "valid_held_folds": list(FOLDS),
            "full_gallery": (
                "all challenge targets present in both held Orion donors; "
                "PerturbMatch and CD4 missingness use frozen fallback"
            ),
            "perturbmatch_fraction": PM_FRACTION,
            "perturbmatch_fallback": (
                "incumbent K562 retained for missing targets and genes"
            ),
            "candidate_covariance_sources": [
                "k562_gwps_pm50_fallback",
                "non-held Orion donor",
            ],
            "candidate_covariance_source_mass": [0.5, 0.5],
            "candidate_covariance_panel_rows_excluded": True,
            "candidate_covariance_panel_coordinates_excluded": True,
            "whitening_gamma": WHITEN_GAMMA,
            "gamma_grid": [WHITEN_GAMMA],
            "cd4_mode": CD4_MODE,
            "cd4_fraction": CD4_FRACTION,
            "cd4_target_support": "all three conditions only",
            "truth": "raw centered, unnormalized, untransformed",
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
        },
        "axis": {
            "official_genes": len(official),
            "challenge_targets": len(panel),
            "gallery_targets": len(targets),
            "target_order": list(targets),
            "target_order_sha256": sha256_lines(targets),
            "panel_order_sha256": sha256_lines(panel),
            "k562_perturbmatch_shared_non_panel_genes": len(shared_columns),
            "k562_perturbmatch_shared_columns_sha256": hashlib.sha256(
                shared_columns.tobytes()
            ).hexdigest(),
        },
        "pm_slot": pm_slot_audit,
        "folds": folds,
        "macro": {
            "baseline_arm_d_pds": baseline_macro,
            "aligned_whitener_arm_d_pds": candidate_macro,
            "delta": macro_delta,
            "bootstrap": bootstrap,
        },
        "stopping_gate": gate,
        "leakage_audit": {
            "held_truth_used_in_candidate_covariance": False,
            "all_300_panel_rows_excluded_from_candidate_covariance": True,
            "all_300_panel_coordinates_excluded_from_both_operators": True,
            "all_300_panel_coordinates_excluded_from_retrieval": True,
            "candidate_count": 1,
            "gamma_retuned": False,
            "pm_fraction_retuned": False,
            "cd4_retuned": False,
            "k562_truth_evaluated": False,
        },
        "input_hashes": {
            "run_py": sha256_file(Path(__file__)),
            "exp064_run_py": sha256_file(EXP064_PATH),
            "exp062_run_py": sha256_file(EXP064.BASE_PATH),
            "panel_csv": sha256_file(args.panel_csv),
            "genes_csv": sha256_file(args.genes_csv),
            "banks": {
                **{name: bank_hashes(bank) for name, bank in banks.items()},
                "perturbmatch": bank_hashes(perturbmatch, include_b_ctrl=False),
                **{bank.name: bank_hashes(bank) for bank in cd4_banks},
            },
        },
        "infrastructure": {
            "provider": "AWS EC2",
            "instance_id": args.instance_id,
            "instance_type": args.instance_type,
            "market": args.market,
            "s3_prefix": args.s3_prefix,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(
        args.output.parent / "fold_summary.csv",
        [
            {
                "held": fold["held"],
                "baseline_arm_d_pds": fold["baseline_arm_d_pds"],
                "aligned_whitener_arm_d_pds": fold["aligned_whitener_arm_d_pds"],
                "delta": fold["delta"],
                "improved_targets": fold["improved_targets"],
                "worsened_targets": fold["worsened_targets"],
                "tied_targets": fold["tied_targets"],
            }
            for fold in folds
        ],
    )
    write_csv(args.output.parent / "per_target_pds.csv", per_target_rows)
    print(json.dumps({"macro": result["macro"], "gate": gate}, indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--baseline-whitener-dir", type=Path, required=True)
    parser.add_argument(
        "--candidate-whitener-dir",
        type=Path,
        default=HERE / "out/whiteners",
    )
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "out/results.json")
    parser.add_argument("--chunk-rows", type=int, default=64)
    parser.add_argument("--max-rss-gib", type=float, default=55.0)
    parser.add_argument("--expected-panel-targets", type=int, default=300)
    parser.add_argument("--expected-full-gallery-targets", type=int, default=289)
    parser.add_argument("--min-pm-shared-genes", type=int, default=7_000)
    parser.add_argument("--min-operator-genes", type=int, default=17_000)
    parser.add_argument("--min-retrieval-genes", type=int, default=17_000)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED)
    parser.add_argument(
        "--expected-hct-whitener-sha256",
        default=EXPECTED_BASELINE_WHITENER_SHA256["orion_hct116"],
    )
    parser.add_argument(
        "--expected-hek-whitener-sha256",
        default=EXPECTED_BASELINE_WHITENER_SHA256["orion_hek293t"],
    )
    parser.add_argument("--instance-id", default="unknown")
    parser.add_argument("--instance-type", default="unknown")
    parser.add_argument("--market", default="unknown")
    parser.add_argument("--s3-prefix", default="unknown")
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing unguarded real-data execution")
    if args.chunk_rows < 1 or args.max_rss_gib <= 0:
        raise SystemExit("chunk rows and memory limit must be positive")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
