#!/usr/bin/env python3
"""Leak-free three-fold PDS/cosine Pareto test for residual de-whitening.

The real-data run is cloud-only. It fits one production-style source-only
union operator at a time and evaluates a frozen convex grid. No count matrix,
cell-eval2 process, challenge result, or submission artifact is touched.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import os
import resource
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

DONOR_NAMES = ("k562_gwps", "orion_hct116", "orion_hek293t")
BLEND_GRID = (0.00, 0.05, 0.10, 0.20, 0.35, 0.50, 1.00)
WHITEN_GAMMA = 0.65
EFFECT_SCALE = 6.817
PDS_SPAN = 0.4525
MSE_BASELINE = 0.989
MSE_REPLICATE = 0.0365
OVERALL_GAP = 0.007


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached "
            f"{limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def sha256_array(array: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(array)
    return hashlib.sha256(memoryview(contiguous).cast("B")).hexdigest()


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column!r}")
        values = tuple(
            value
            for row in reader
            if (value := str(row.get(column, "")).strip())
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column!r}")
    return values


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    pert_pos: dict[str, int]
    delta_path: Path
    b_ctrl_path: Path
    shape: tuple[int, int]

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        path = path.expanduser().resolve()
        genes = read_lines(path / "genes.txt", unique=False)
        perts = read_lines(path / "perts.txt")
        delta_path = path / "delta.npy"
        b_ctrl_path = path / "b_ctrl.npy"
        delta = np.load(delta_path, mmap_mode="r")
        b_ctrl = np.load(b_ctrl_path, mmap_mode="r")
        if not isinstance(delta, np.memmap):
            raise TypeError(f"{name}: delta.npy is not memory mapped")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != "
                f"{(len(perts), len(genes))}"
            )
        if b_ctrl.shape != (len(genes),):
            raise ValueError(
                f"{name}: b_ctrl shape {b_ctrl.shape} != {(len(genes),)}"
            )
        shape = tuple(map(int, delta.shape))
        del delta, b_ctrl
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            b_ctrl_path=b_ctrl_path,
            shape=shape,
        )


@dataclass(frozen=True)
class AxisMap:
    all_source_columns: np.ndarray
    unique_source_columns: np.ndarray
    unique_target_positions: np.ndarray
    duplicate_symbols: tuple[str, ...]


def map_bank_to_axis(
    bank: Bank,
    official_genes: Sequence[str],
) -> AxisMap:
    target_pos = {gene: index for index, gene in enumerate(official_genes)}
    all_source: list[int] = []
    last_source: dict[int, int] = {}
    duplicates: set[str] = set()
    seen: set[int] = set()
    for source_column, gene in enumerate(bank.genes):
        destination = target_pos.get(gene)
        if destination is None:
            continue
        all_source.append(source_column)
        if destination in seen:
            duplicates.add(gene)
        seen.add(destination)
        last_source[destination] = source_column
    if not all_source:
        raise ValueError(f"{bank.name}: no genes map to official axis")
    unique_target = np.asarray(sorted(last_source), dtype=np.int64)
    unique_source = np.asarray(
        [last_source[position] for position in unique_target], dtype=np.int64
    )
    return AxisMap(
        all_source_columns=np.asarray(all_source, dtype=np.int64),
        unique_source_columns=unique_source,
        unique_target_positions=unique_target,
        duplicate_symbols=tuple(sorted(duplicates)),
    )


def stream_column_mean(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    delta = np.load(bank.delta_path, mmap_mode="r")
    total = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        total += np.asarray(
            delta[lo:hi, columns], dtype=np.float64
        ).sum(axis=0, dtype=np.float64)
        enforce_memory(max_rss_gib, f"{bank.name} source centering")
    del delta
    return total / bank.shape[0]


def stream_median_centered_norm(
    bank: Bank,
    columns: np.ndarray,
    mean: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> float:
    delta = np.load(bank.delta_path, mmap_mode="r")
    norms = np.empty(bank.shape[0], dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        centered = np.asarray(
            delta[lo:hi, columns], dtype=np.float64
        ) - mean
        norms[lo:hi] = np.linalg.norm(centered, axis=1)
        enforce_memory(max_rss_gib, f"{bank.name} source scale")
    del delta
    return float(np.median(norms)) or 1.0


@dataclass
class PreparedSource:
    name: str
    effects: np.ndarray
    target_available: np.ndarray
    mean_on_axis: np.ndarray
    raw_median_norm: float
    mapped_gene_count: int
    duplicate_symbols: tuple[str, ...]


def prepare_main_source(
    bank: Bank,
    targets: Sequence[str],
    official_genes: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> PreparedSource:
    """Reproduce the raw, no-EB experiment-031 source preprocessing."""
    mapping = map_bank_to_axis(bank, official_genes)
    mean_all = stream_column_mean(
        bank,
        mapping.all_source_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    raw_median = stream_median_centered_norm(
        bank,
        mapping.all_source_columns,
        mean_all,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    all_column_pos = {
        int(source): index
        for index, source in enumerate(mapping.all_source_columns.tolist())
    }
    unique_mean = np.asarray(
        [
            mean_all[all_column_pos[int(source)]]
            for source in mapping.unique_source_columns
        ],
        dtype=np.float64,
    )
    mean_on_axis = np.zeros(len(official_genes), dtype=np.float64)
    mean_on_axis[mapping.unique_target_positions] = unique_mean

    effects = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    available = np.zeros(len(targets), dtype=bool)
    delta = np.load(bank.delta_path, mmap_mode="r")
    for target_index, target in enumerate(targets):
        row = bank.pert_pos.get(target)
        if row is None:
            continue
        values = np.asarray(
            delta[row, mapping.unique_source_columns], dtype=np.float64
        )
        effects[target_index, mapping.unique_target_positions] = (
            values - unique_mean
        ) * (EFFECT_SCALE / raw_median)
        available[target_index] = True
    del delta
    enforce_memory(max_rss_gib, f"{bank.name} target extraction")
    return PreparedSource(
        name=bank.name,
        effects=effects,
        target_available=available,
        mean_on_axis=mean_on_axis,
        raw_median_norm=raw_median,
        mapped_gene_count=len(mapping.unique_target_positions),
        duplicate_symbols=mapping.duplicate_symbols,
    )


def build_current_baseline(
    sources: Sequence[PreparedSource],
) -> tuple[np.ndarray, np.ndarray]:
    if not sources:
        raise ValueError("at least one source donor is required")
    shape = sources[0].effects.shape
    baseline = np.zeros(shape, dtype=np.float64)
    votes = np.zeros(shape[0], dtype=np.int64)
    for source in sources:
        if source.effects.shape != shape:
            raise ValueError("source effect matrices have different shapes")
        rows = source.target_available
        baseline[rows] += source.effects[rows]
        votes += rows.astype(np.int64)
    covered = votes > 0
    baseline[covered] /= votes[covered, None]
    missing = ~covered
    if missing.any():
        mean_response = sum(
            source.mean_on_axis for source in sources
        ) / len(sources)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        baseline[missing] = EFFECT_SCALE * mean_response / mean_norm
    return baseline, votes


@dataclass(frozen=True)
class Operator:
    evals: np.ndarray
    evecs: np.ndarray
    wcols: np.ndarray
    audit: dict[str, object]


def _last_source_columns(
    bank: Bank,
    official_genes: Sequence[str],
    selected_official_positions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mapping = map_bank_to_axis(bank, official_genes)
    selected_to_axis = {
        int(position): index
        for index, position in enumerate(selected_official_positions.tolist())
    }
    source_columns: list[int] = []
    operator_columns: list[int] = []
    for source, destination in zip(
        mapping.unique_source_columns,
        mapping.unique_target_positions,
    ):
        operator_column = selected_to_axis.get(int(destination))
        if operator_column is not None:
            source_columns.append(int(source))
            operator_columns.append(operator_column)
    return (
        np.asarray(source_columns, dtype=np.int64),
        np.asarray(operator_columns, dtype=np.int64),
    )


def fit_union_operator(
    source_banks: Sequence[Bank],
    prepared_sources: Sequence[PreparedSource],
    official: Sequence[str],
    panel_set: set[str],
    scratch: Path,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> Operator:
    """Fit the exact raw, equal-source production union covariance."""
    source_names = tuple(bank.name for bank in source_banks)
    source_union = set().union(*(set(bank.genes) for bank in source_banks))
    wcols = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in source_union and gene not in panel_set
        ],
        dtype=np.int64,
    )
    operator_genes = tuple(official[index] for index in wcols)
    if panel_set & set(operator_genes):
        raise AssertionError("panel target coordinate entered union operator")

    scratch.mkdir(parents=True, exist_ok=True)
    basis_path = scratch / "source_basis.npy"
    total_rows = sum(bank.shape[0] for bank in source_banks)
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float64,
        shape=(total_rows, len(wcols)),
    )
    source_scale = {
        source.name: source.raw_median_norm for source in prepared_sources
    }
    reference = float(np.median(list(source_scale.values()))) or 1.0
    source_audit: dict[str, object] = {}
    row_offset = 0
    for bank in source_banks:
        source_columns, operator_columns = _last_source_columns(
            bank, official, wcols
        )
        mean = stream_column_mean(
            bank,
            source_columns,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        multiplier = reference / source_scale[bank.name]
        delta = np.load(bank.delta_path, mmap_mode="r")
        try:
            for lo in range(0, bank.shape[0], chunk_rows):
                hi = min(lo + chunk_rows, bank.shape[0])
                destination = basis[row_offset + lo : row_offset + hi]
                destination[:] = 0.0
                destination[:, operator_columns] = (
                    np.asarray(
                        delta[lo:hi, source_columns], dtype=np.float64
                    )
                    - mean
                ) * multiplier
                enforce_memory(max_rss_gib, f"{bank.name} basis")
        finally:
            del delta
        source_audit[bank.name] = {
            "basis_rows": int(bank.shape[0]),
            "mapped_operator_genes": int(len(operator_columns)),
            "raw_median_centered_norm": source_scale[bank.name],
            "basis_relative_scale": multiplier,
            "delta_sha256": sha256_file(bank.delta_path),
            "genes_sha256": sha256_lines(bank.genes),
            "perts_sha256": sha256_lines(bank.perts),
        }
        row_offset += bank.shape[0]
        basis.flush()

    print(
        f"fit {source_names}: basis={basis.shape}, peak={peak_rss_gib():.2f} GiB",
        flush=True,
    )
    covariance = (basis.T @ basis) / float(total_rows)
    del basis
    basis_path.unlink(missing_ok=True)
    gc.collect()
    enforce_memory(max_rss_gib, "covariance")

    evals, evecs = np.linalg.eigh(covariance)
    del covariance
    gc.collect()
    evals = np.maximum(evals.astype(np.float64), float(evals[-1]) * 1e-12)
    evecs = np.asarray(evecs, dtype=np.float32)
    if not np.isfinite(evals).all() or not np.isfinite(evecs[::257, ::257]).all():
        raise FloatingPointError("non-finite union operator")
    enforce_memory(max_rss_gib, "eigendecomposition")

    audit = {
        "source_donors": list(source_names),
        "source_only": True,
        "operator_genes": int(len(wcols)),
        "operator_genes_sha256": sha256_lines(operator_genes),
        "operator_columns_sha256": sha256_array(wcols),
        "excluded_target_coordinates": len(panel_set),
        "panel_overlap_with_operator": len(panel_set & set(operator_genes)),
        "basis_rows": total_rows,
        "basis_dtype": "float64",
        "basis_normalization": (
            "production equal-source renorm: median source norm / donor norm"
        ),
        "evals_sha256": sha256_array(evals),
        "evecs_float32_sha256": sha256_array(evecs),
        "sources": source_audit,
    }
    return Operator(evals=evals, evecs=evecs, wcols=wcols, audit=audit)


def apply_operator(
    matrix: np.ndarray,
    operator: Operator,
) -> tuple[np.ndarray, dict[str, float]]:
    output = np.asarray(matrix, dtype=np.float64).copy()
    block = output[:, operator.wcols]
    before = float(np.median(np.linalg.norm(block, axis=1))) or 1.0
    weights = operator.evals ** (-WHITEN_GAMMA / 2.0)
    transformed = ((block @ operator.evecs) * weights) @ operator.evecs.T
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    scale = before / after
    output[:, operator.wcols] = transformed * scale
    if not np.isfinite(output).all():
        raise FloatingPointError("non-finite whitened prediction")
    return output, {
        "before_operator_median_norm": before,
        "after_operator_median_norm": after,
        "applied_scale": scale,
    }


def load_centered_truth(
    bank: Bank,
    targets: Sequence[str],
    official: Sequence[str],
    retrieval_columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    """Open held-donor values only after predictions have been frozen."""
    source_columns, retrieval_positions = _last_source_columns(
        bank, official, retrieval_columns
    )
    if len(source_columns) != len(retrieval_columns):
        raise ValueError(f"{bank.name}: incomplete truth retrieval mapping")
    if not np.array_equal(
        retrieval_positions, np.arange(len(retrieval_columns))
    ):
        raise ValueError(f"{bank.name}: unordered truth retrieval mapping")
    mean = stream_column_mean(
        bank,
        source_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    rows = np.asarray(
        [bank.pert_pos[target] for target in targets], dtype=np.int64
    )
    delta = np.load(bank.delta_path, mmap_mode="r")
    truth = np.asarray(
        delta[np.ix_(rows, source_columns)], dtype=np.float64
    ) - mean
    del delta
    enforce_memory(max_rss_gib, f"{bank.name} truth")
    return truth


def row_cosines(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    if prediction.shape != truth.shape:
        raise ValueError(
            f"cosine shape mismatch {prediction.shape} vs {truth.shape}"
        )
    numerator = np.sum(prediction * truth, axis=1)
    denominator = np.linalg.norm(prediction, axis=1) * np.linalg.norm(
        truth, axis=1
    )
    return numerator / np.maximum(denominator, 1e-12)


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    prediction = np.asarray(prediction, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    prediction = prediction / np.maximum(
        np.linalg.norm(prediction, axis=1, keepdims=True), 1e-12
    )
    truth = truth / np.maximum(
        np.linalg.norm(truth, axis=1, keepdims=True), 1e-12
    )
    similarities = prediction @ truth.T
    own = np.diag(similarities)[:, None]
    n = similarities.shape[0]
    if similarities.shape != (n, n) or n < 2:
        raise ValueError("PDS requires an aligned full gallery with n >= 2")
    ranks = (similarities > own).sum(axis=1).astype(np.float64)
    ranks += 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


def metric_rows(
    prediction: np.ndarray,
    truth: np.ndarray,
) -> tuple[dict[str, float | int], np.ndarray, np.ndarray]:
    cosines = row_cosines(prediction, truth)
    pds_values = pds_per_target(prediction, truth)
    min_mse = 1.0 - np.square(cosines)
    nonnegative_min_mse = 1.0 - np.square(np.maximum(cosines, 0.0))
    return (
        {
            "pds": float(np.mean(pds_values)),
            "raw_cosine": float(np.mean(cosines)),
            "min_mse_unrestricted": float(np.mean(min_mse)),
            "min_mse_nonnegative_alpha": float(
                np.mean(nonnegative_min_mse)
            ),
            "positive_cosine_targets": int(np.sum(cosines > 0)),
            "negative_cosine_targets": int(np.sum(cosines < 0)),
        },
        pds_values,
        cosines,
    )


def pareto_indices(rows: Sequence[dict[str, float]]) -> list[int]:
    """Return non-dominated rows for PDS/cosine/max and both MSE/min."""
    frontier: list[int] = []
    for i, left in enumerate(rows):
        dominated = False
        for j, right in enumerate(rows):
            if i == j:
                continue
            weak = (
                right["pds"] >= left["pds"]
                and right["raw_cosine"] >= left["raw_cosine"]
                and right["min_mse_unrestricted"]
                <= left["min_mse_unrestricted"]
                and right["min_mse_nonnegative_alpha"]
                <= left["min_mse_nonnegative_alpha"]
            )
            strict = (
                right["pds"] > left["pds"]
                or right["raw_cosine"] > left["raw_cosine"]
                or right["min_mse_unrestricted"]
                < left["min_mse_unrestricted"]
                or right["min_mse_nonnegative_alpha"]
                < left["min_mse_nonnegative_alpha"]
            )
            if weak and strict:
                dominated = True
                break
        if not dominated:
            frontier.append(i)
    return frontier


def score_fold(
    held_name: str,
    banks: dict[str, Bank],
    targets: Sequence[str],
    official: Sequence[str],
    panel_set: set[str],
    scratch_root: Path,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    sources = tuple(name for name in DONOR_NAMES if name != held_name)
    prepared = [
        prepare_main_source(
            banks[name],
            targets,
            official,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        for name in sources
    ]
    raw_prediction, source_votes = build_current_baseline(prepared)
    operator = fit_union_operator(
        [banks[name] for name in sources],
        prepared,
        official,
        panel_set,
        scratch_root / held_name,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    whitened_prediction, renormalization = apply_operator(
        raw_prediction, operator
    )

    # Freeze the complete predeclared grid before opening held-donor values.
    predictions = {
        beta: (1.0 - beta) * whitened_prediction + beta * raw_prediction
        for beta in BLEND_GRID
    }
    held_gene_set = set(banks[held_name].genes)
    retrieval_columns = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in held_gene_set and gene not in panel_set
        ],
        dtype=np.int64,
    )
    retrieval_genes = tuple(official[index] for index in retrieval_columns)
    if panel_set & set(retrieval_genes):
        raise AssertionError("panel target entered retrieval axis")
    off_operator = np.setdiff1d(
        retrieval_columns, operator.wcols, assume_unique=True
    )

    truth = load_centered_truth(
        banks[held_name],
        targets,
        official,
        retrieval_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    fold_rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    for beta in BLEND_GRID:
        prediction = predictions[beta][:, retrieval_columns]
        summary, pds_values, cosines = metric_rows(prediction, truth)
        fold_rows.append({"beta": beta, **summary})
        for target, pds_value, cosine in zip(
            targets, pds_values, cosines
        ):
            per_target_rows.append(
                {
                    "held_donor": held_name,
                    "target_gene": target,
                    "beta_raw_share": beta,
                    "pds": float(pds_value),
                    "raw_cosine": float(cosine),
                    "min_mse_unrestricted": float(1.0 - cosine**2),
                    "min_mse_nonnegative_alpha": float(
                        1.0 - max(cosine, 0.0) ** 2
                    ),
                }
            )
        print(
            f"{held_name} beta={beta:.2f} pds={summary['pds']:.6f} "
            f"cos={summary['raw_cosine']:.6f} "
            f"min_mse+={summary['min_mse_nonnegative_alpha']:.6f}",
            flush=True,
        )

    fold_result = {
        "held_donor": held_name,
        "source_donors": list(sources),
        "targets": len(targets),
        "target_order_sha256": sha256_lines(targets),
        "retrieval_genes": len(retrieval_columns),
        "retrieval_genes_sha256": sha256_lines(retrieval_genes),
        "retrieval_columns_sha256": sha256_array(retrieval_columns),
        "retrieval_genes_outside_operator": int(len(off_operator)),
        "panel_overlap_with_retrieval": len(panel_set & set(retrieval_genes)),
        "source_vote_counts": {
            str(vote): int(np.sum(source_votes == vote))
            for vote in sorted(set(source_votes.tolist()))
        },
        "source_details": [
            {
                "name": source.name,
                "raw_median_centered_norm": source.raw_median_norm,
                "mapped_genes": source.mapped_gene_count,
                "duplicate_symbols": list(source.duplicate_symbols),
            }
            for source in prepared
        ],
        "operator": operator.audit,
        "renormalization": renormalization,
        "metrics": fold_rows,
        "held_delta_sha256": sha256_file(banks[held_name].delta_path),
        "predictions_frozen_before_held_values_opened": True,
        "peak_rss_gib": peak_rss_gib(),
    }
    del (
        operator,
        predictions,
        truth,
        whitened_prediction,
        raw_prediction,
        prepared,
    )
    shutil.rmtree(scratch_root / held_name, ignore_errors=True)
    gc.collect()
    enforce_memory(max_rss_gib, f"finish {held_name}")
    return fold_result, per_target_rows


def optimistic_gap_translation(
    rows: list[dict[str, float]],
) -> list[dict[str, float | bool]]:
    baseline = next(row for row in rows if row["beta"] == 0.0)

    def mse_scaled(value: float) -> float:
        score = (MSE_BASELINE - value) / (MSE_BASELINE - MSE_REPLICATE)
        return float(np.clip(score, 0.0, 1.0))

    baseline_mse_scaled = mse_scaled(
        baseline["min_mse_nonnegative_alpha"]
    )
    output = []
    for row in rows:
        pds_scaled_delta = (row["pds"] - baseline["pds"]) / PDS_SPAN
        mse_scaled_delta = (
            mse_scaled(row["min_mse_nonnegative_alpha"])
            - baseline_mse_scaled
        )
        optimistic_delta = (pds_scaled_delta + mse_scaled_delta) / 6.0
        output.append(
            {
                "beta": row["beta"],
                "pds_scaled_delta_vs_beta0": pds_scaled_delta,
                "mse_scaled_delta_vs_beta0": mse_scaled_delta,
                "optimistic_overall_delta_vs_beta0": optimistic_delta,
                "could_close_0.007_gap": optimistic_delta >= OVERALL_GAP,
            }
        )
    return output


def _write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = args.output_dir.expanduser().resolve()
    if HERE.resolve() not in (output_dir, *output_dir.parents):
        raise ValueError("all outputs must stay inside experiment 063")
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_root = output_dir / "scratch"
    shutil.rmtree(scratch_root, ignore_errors=True)

    panel = read_csv_column(args.panel_csv, "target_gene")
    official = read_csv_column(args.genes_csv, "gene_name")
    panel_set = set(panel)
    if len(panel) != 300 or len(panel_set) != 300:
        raise ValueError(f"expected 300 unique panel targets, got {len(panel_set)}")
    if tuple(args.blend_grid) != BLEND_GRID:
        raise ValueError("blend grid differs from the predeclared frozen grid")
    if args.gamma != WHITEN_GAMMA:
        raise ValueError("gamma differs from the predeclared fixed value")

    paths = {
        "k562_gwps": args.k562,
        "orion_hct116": args.hct116,
        "orion_hek293t": args.hek293t,
    }
    banks = {name: Bank.load(name, path) for name, path in paths.items()}
    targets = tuple(
        target
        for target in panel
        if all(target in banks[name].pert_pos for name in DONOR_NAMES)
    )
    if len(targets) < args.min_shared_targets:
        raise ValueError(
            f"only {len(targets)} targets shared by all donors; "
            f"require {args.min_shared_targets}"
        )

    folds: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    for held_name in DONOR_NAMES:
        fold, rows = score_fold(
            held_name,
            banks,
            targets,
            official,
            panel_set,
            scratch_root,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        folds.append(fold)
        per_target_rows.extend(rows)

    macro_rows: list[dict[str, float]] = []
    for beta in BLEND_GRID:
        selected = [
            next(row for row in fold["metrics"] if row["beta"] == beta)
            for fold in folds
        ]
        macro_rows.append(
            {
                "beta": beta,
                "pds": float(np.mean([row["pds"] for row in selected])),
                "raw_cosine": float(
                    np.mean([row["raw_cosine"] for row in selected])
                ),
                "min_mse_unrestricted": float(
                    np.mean(
                        [row["min_mse_unrestricted"] for row in selected]
                    )
                ),
                "min_mse_nonnegative_alpha": float(
                    np.mean(
                        [
                            row["min_mse_nonnegative_alpha"]
                            for row in selected
                        ]
                    )
                ),
            }
        )
    frontier_indices = pareto_indices(macro_rows)
    frontier_betas = [macro_rows[index]["beta"] for index in frontier_indices]
    gap_translation = optimistic_gap_translation(macro_rows)
    plausible = [
        row["beta"]
        for row in gap_translation
        if row["could_close_0.007_gap"]
    ]

    result: dict[str, object] = {
        "experiment": "063-pds-mse-pareto",
        "status": "completed",
        "protocol": {
            "predeclared_beta_raw_share": list(BLEND_GRID),
            "formula": "(1-beta)*union_whitened + beta*raw_consensus",
            "gamma": WHITEN_GAMMA,
            "source_weighting": "equal across available source donors",
            "source_effects": (
                "raw centered donor effects with production per-donor norm "
                "and global effect-scale normalization"
            ),
            "whitening_axis": "source-donor union on official VCC axis",
            "whitening_renormalization": "production preserve median axis norm",
            "outer_holdouts": list(DONOR_NAMES),
            "targets": len(targets),
            "metric": (
                "full-gallery rank PDS, mean row cosine, and row-wise "
                "theoretical minimum uncorrected normalized MSE"
            ),
            "count_emission": False,
            "cell_eval2": False,
            "challenge_outcomes_used": False,
            "submission": False,
        },
        "leakage_audit": {
            "challenge_target_coordinates_excluded": len(panel_set),
            "challenge_target_order_sha256": sha256_lines(panel),
            "shared_gallery_targets": len(targets),
            "shared_gallery_order_sha256": sha256_lines(targets),
            "held_donor_rows_in_operator": 0,
            "held_values_opened_after_predictions_frozen": all(
                fold["predictions_frozen_before_held_values_opened"]
                for fold in folds
            ),
            "panel_overlap_with_any_operator": max(
                fold["operator"]["panel_overlap_with_operator"]
                for fold in folds
            ),
            "panel_overlap_with_any_retrieval_axis": max(
                fold["panel_overlap_with_retrieval"] for fold in folds
            ),
        },
        "folds": folds,
        "macro": macro_rows,
        "pareto_frontier_betas": frontier_betas,
        "leaderboard_gap_proxy": {
            "current_overall_gap": OVERALL_GAP,
            "anchors": {
                "pds_span": PDS_SPAN,
                "mse_baseline": MSE_BASELINE,
                "mse_replicate": MSE_REPLICATE,
            },
            "assumption": (
                "optimistic: oracle nonnegative per-target alpha, other four "
                "metrics unchanged"
            ),
            "rows": gap_translation,
            "plausible_betas": plausible,
            "any_candidate_could_close_gap": bool(plausible),
        },
        "infrastructure": {
            "provider": "AWS EC2",
            "instance_id": args.instance_id,
            "instance_type": args.instance_type,
            "region": args.region,
            "s3_prefix": args.s3_prefix,
            "process_tree_guard_gib": args.max_rss_gib,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    results_path = output_dir / "results.json"
    results_path.write_text(json.dumps(result, indent=2) + "\n")
    _write_csv(output_dir / "summary.csv", macro_rows)
    _write_csv(output_dir / "per_target_metrics.csv", per_target_rows)
    _write_csv(output_dir / "gap_translation.csv", gap_translation)
    shutil.rmtree(scratch_root, ignore_errors=True)
    print(
        json.dumps(
            {
                "macro": macro_rows,
                "pareto_frontier_betas": frontier_betas,
                "plausible_betas": plausible,
                "peak_rss_gib": peak_rss_gib(),
            },
            indent=2,
        ),
        flush=True,
    )
    return result


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k562", type=Path, required=True)
    ap.add_argument("--hct116", type=Path, required=True)
    ap.add_argument("--hek293t", type=Path, required=True)
    ap.add_argument("--genes-csv", type=Path, required=True)
    ap.add_argument("--panel-csv", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, default=HERE / "out")
    ap.add_argument("--chunk-rows", type=int, default=128)
    ap.add_argument("--max-rss-gib", type=float, default=55.0)
    ap.add_argument("--min-shared-targets", type=int, default=250)
    ap.add_argument(
        "--blend-grid",
        nargs="+",
        type=float,
        default=list(BLEND_GRID),
    )
    ap.add_argument("--gamma", type=float, default=WHITEN_GAMMA)
    ap.add_argument("--instance-id", default="unknown")
    ap.add_argument("--instance-type", default="unknown")
    ap.add_argument("--region", default="unknown")
    ap.add_argument("--s3-prefix", default="unknown")
    return ap


def main() -> None:
    args = parser().parse_args()
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    if not 0 < args.max_rss_gib <= 55:
        raise SystemExit("--max-rss-gib must be in (0, 55]")
    run(args)


if __name__ == "__main__":
    main()
