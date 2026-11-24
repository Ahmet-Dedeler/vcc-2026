"""062 — add an unwhitened CD4 mean3 residual after baseline whitening.

Select one mode/fraction on held HCT116 and held HEK293T, then evaluate that
frozen choice exactly once in the held-K562 target role. Every fold fits its
own gamma-0.65 union whitener from source donors only. All VCC-panel target
coordinates are excluded from fitting, transformation, and retrieval.
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

MAIN_NAMES = ("k562_gwps", "orion_hct116", "orion_hek293t")
DEVELOPMENT_HELD = ("orion_hct116", "orion_hek293t")
FINAL_HELD = "k562_gwps"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
MODES = ("raw", "gene_norm")
FRACTIONS = (0.01, 0.02, 0.05, 0.10, 0.20)
EFFECT_SCALE = 6.817
WHITEN_GAMMA = 0.65


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
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
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: required column {column!r} is missing")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: {column!r} contains duplicate labels")
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
    def load(cls, name: str, path: Path, *, require_cd4_scale: bool = False) -> Bank:
        path = path.expanduser().resolve()
        genes = read_lines(path / "genes.txt", require_unique=False)
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
            raise ValueError(f"{name}: b_ctrl shape {b_ctrl.shape} is invalid")
        if require_cd4_scale:
            metadata = json.loads((path / "meta.json").read_text())
            definition = str(metadata.get("effect_definition", ""))
            if "ln(2)" not in definition or "log2 fold change" not in definition:
                raise ValueError(f"{name}: unrecognized CD4 scale {definition!r}")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            b_ctrl_path=b_ctrl_path,
            shape=tuple(map(int, delta.shape)),
        )


@dataclass(frozen=True)
class AxisMap:
    all_source_columns: np.ndarray
    unique_source_columns: np.ndarray
    unique_target_positions: np.ndarray
    duplicate_symbols: tuple[str, ...]


def map_bank_to_axis(bank: Bank, official_genes: Sequence[str]) -> AxisMap:
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
        total += np.asarray(delta[lo:hi, columns], dtype=np.float64).sum(
            axis=0, dtype=np.float64
        )
        enforce_memory_limit(max_rss_gib, f"{bank.name} centering")
    return total / bank.shape[0]


def centered_statistics(
    bank: Bank,
    columns: np.ndarray,
    mean: np.ndarray,
    *,
    need_gene_scale: bool,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[float, np.ndarray | None, float | None]:
    delta = np.load(bank.delta_path, mmap_mode="r")
    raw_norms = np.empty(bank.shape[0], dtype=np.float64)
    absolute = np.zeros(len(columns), dtype=np.float64) if need_gene_scale else None
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        centered = np.asarray(delta[lo:hi, columns], dtype=np.float64) - mean
        raw_norms[lo:hi] = np.linalg.norm(centered, axis=1)
        if absolute is not None:
            absolute += np.abs(centered).sum(axis=0, dtype=np.float64)
        enforce_memory_limit(max_rss_gib, f"{bank.name} scale")
    raw_median = float(np.median(raw_norms)) or 1.0
    if absolute is None:
        return raw_median, None, None
    mean_abs = absolute / bank.shape[0]
    live = mean_abs[mean_abs > 0]
    floor = float(np.percentile(live, 25.0)) if live.size else 1.0
    gene_scale = np.maximum(mean_abs, floor)
    normalized_norms = np.empty(bank.shape[0], dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        centered = np.asarray(delta[lo:hi, columns], dtype=np.float64) - mean
        normalized_norms[lo:hi] = np.linalg.norm(centered / gene_scale, axis=1)
        enforce_memory_limit(max_rss_gib, f"{bank.name} gene-normalized scale")
    return raw_median, gene_scale, float(np.median(normalized_norms)) or 1.0


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
    panel_positions: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> PreparedSource:
    mapping = map_bank_to_axis(bank, official_genes)
    mean_all = stream_column_mean(
        bank,
        mapping.all_source_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    raw_median, _, _ = centered_statistics(
        bank,
        mapping.all_source_columns,
        mean_all,
        need_gene_scale=False,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    all_position = {
        int(source): index
        for index, source in enumerate(mapping.all_source_columns.tolist())
    }
    unique_mean = np.asarray(
        [mean_all[all_position[int(source)]] for source in mapping.unique_source_columns],
        dtype=np.float64,
    )
    mean_on_axis = np.zeros(len(official_genes), dtype=np.float64)
    mean_on_axis[mapping.unique_target_positions] = unique_mean
    effects = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    available = np.zeros(len(targets), dtype=bool)
    delta = np.load(bank.delta_path, mmap_mode="r")
    for index, target in enumerate(targets):
        row = bank.pert_pos.get(target)
        if row is None:
            continue
        values = np.asarray(
            delta[row, mapping.unique_source_columns], dtype=np.float64
        )
        effects[index, mapping.unique_target_positions] = (
            values - unique_mean
        ) * (EFFECT_SCALE / raw_median)
        available[index] = True
    effects[:, panel_positions] = 0.0
    mean_on_axis[panel_positions] = 0.0
    return PreparedSource(
        name=bank.name,
        effects=effects,
        target_available=available,
        mean_on_axis=mean_on_axis,
        raw_median_norm=raw_median,
        mapped_gene_count=len(mapping.unique_target_positions),
        duplicate_symbols=mapping.duplicate_symbols,
    )


def build_baseline(
    sources: Sequence[PreparedSource],
) -> tuple[np.ndarray, np.ndarray]:
    baseline = np.zeros_like(sources[0].effects)
    votes = np.zeros(baseline.shape[0], dtype=np.int64)
    for source in sources:
        rows = source.target_available
        baseline[rows] += source.effects[rows]
        votes += rows.astype(np.int64)
    covered = votes > 0
    baseline[covered] /= votes[covered, None]
    if (~covered).any():
        mean_response = sum(source.mean_on_axis for source in sources) / len(sources)
        mean_norm = float(np.linalg.norm(mean_response)) or 1.0
        baseline[~covered] = EFFECT_SCALE * mean_response / mean_norm
    return baseline, votes


@dataclass
class CD4Prepared:
    mode: str
    mean_effect: np.ndarray
    gene_votes: np.ndarray
    target_condition_count: np.ndarray
    condition_details: list[dict[str, object]]


def prepare_cd4_mean(
    banks: Sequence[Bank],
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel_targets: set[str],
    mode: str,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> CD4Prepared:
    if mode not in MODES:
        raise ValueError(f"unknown CD4 mode {mode!r}")
    official_pos = {gene: index for index, gene in enumerate(official_genes)}
    value_sum = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    gene_votes = np.zeros(value_sum.shape, dtype=np.uint8)
    condition_count = np.zeros(len(targets), dtype=np.uint8)
    details: list[dict[str, object]] = []
    for bank in banks:
        source_columns: list[int] = []
        target_positions: list[int] = []
        for source_column, gene in enumerate(bank.genes):
            destination = official_pos.get(gene)
            if destination is None or gene in panel_targets:
                continue
            source_columns.append(source_column)
            target_positions.append(destination)
        source_columns_array = np.asarray(source_columns, dtype=np.int64)
        target_positions_array = np.asarray(target_positions, dtype=np.int64)
        if len(np.unique(target_positions_array)) != len(target_positions_array):
            raise ValueError(f"{bank.name}: duplicate compact CD4 symbols")
        mean = stream_column_mean(
            bank,
            source_columns_array,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        raw_median, gene_scale, normalized_median = centered_statistics(
            bank,
            source_columns_array,
            mean,
            need_gene_scale=mode == "gene_norm",
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        if mode == "gene_norm":
            assert gene_scale is not None and normalized_median is not None
            multiplier = EFFECT_SCALE / normalized_median
        else:
            multiplier = EFFECT_SCALE / raw_median
        delta = np.load(bank.delta_path, mmap_mode="r")
        covered = 0
        for index, target in enumerate(targets):
            row = bank.pert_pos.get(target)
            if row is None:
                continue
            centered = (
                np.asarray(delta[row, source_columns_array], dtype=np.float64) - mean
            )
            if mode == "gene_norm":
                centered = centered / gene_scale
            value_sum[index, target_positions_array] += centered * multiplier
            gene_votes[index, target_positions_array] += 1
            condition_count[index] += 1
            covered += 1
        details.append(
            {
                "condition": bank.name,
                "background_perturbations": bank.shape[0],
                "mapped_non_panel_genes": len(source_columns_array),
                "covered_targets": covered,
                "raw_median_centered_norm": raw_median,
                "gene_normalized_median_before_rescale": normalized_median,
            }
        )
    available = gene_votes > 0
    mean_effect = np.zeros_like(value_sum)
    mean_effect[available] = value_sum[available] / gene_votes[available]
    return CD4Prepared(
        mode=mode,
        mean_effect=mean_effect,
        gene_votes=gene_votes,
        target_condition_count=condition_count,
        condition_details=details,
    )


@dataclass(frozen=True)
class Operator:
    evals: np.ndarray
    evecs: np.ndarray
    metadata: dict[str, object]
    path: Path
    metadata_path: Path


def last_source_columns(
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
        mapping.unique_source_columns, mapping.unique_target_positions, strict=True
    ):
        operator_column = selected_to_axis.get(int(destination))
        if operator_column is not None:
            source_columns.append(int(source))
            operator_columns.append(operator_column)
    return (
        np.asarray(source_columns, dtype=np.int64),
        np.asarray(operator_columns, dtype=np.int64),
    )


def fit_union_whitener(
    source_banks: Sequence[Bank],
    prepared_sources: Sequence[PreparedSource],
    official_genes: Sequence[str],
    panel: Sequence[str],
    output_path: Path,
    scratch: Path,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[Operator, np.ndarray]:
    source_names = tuple(bank.name for bank in source_banks)
    panel_set = set(panel)
    donor_union = set().union(*(set(bank.genes) for bank in source_banks))
    wcols = np.asarray(
        [
            index
            for index, gene in enumerate(official_genes)
            if gene in donor_union and gene not in panel_set
        ],
        dtype=np.int64,
    )
    metadata: dict[str, object] = {
        "axis": "union",
        "donors": list(source_names),
        "n_cols": int(len(wcols)),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(panel_set),
        "basis_rows_by_source": {
            bank.name: int(bank.shape[0]) for bank in source_banks
        },
        "basis_normalization": (
            "production equal-source renorm: median source norm / donor norm"
        ),
        "gamma": WHITEN_GAMMA,
    }
    metadata_path = output_path.with_suffix(".json")
    if output_path.exists() and metadata_path.exists():
        actual = json.loads(metadata_path.read_text())
        if actual != metadata:
            raise ValueError(f"{output_path}: stale whitener metadata")
        with np.load(output_path) as archive:
            evals = np.asarray(archive["evals"], dtype=np.float64)
            evecs = np.asarray(archive["evecs"], dtype=np.float32)
        return Operator(evals, evecs, actual, output_path, metadata_path), wcols

    output_path.parent.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    total_rows = sum(bank.shape[0] for bank in source_banks)
    basis_path = scratch / "source_basis.npy"
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float64,
        shape=(total_rows, len(wcols)),
    )
    source_scale = {
        source.name: source.raw_median_norm for source in prepared_sources
    }
    reference_norm = float(np.median(list(source_scale.values()))) or 1.0
    row_offset = 0
    for bank in source_banks:
        source_columns, operator_columns = last_source_columns(
            bank, official_genes, wcols
        )
        mean = stream_column_mean(
            bank,
            source_columns,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        delta = np.load(bank.delta_path, mmap_mode="r")
        multiplier = reference_norm / source_scale[bank.name]
        for lo in range(0, bank.shape[0], chunk_rows):
            hi = min(lo + chunk_rows, bank.shape[0])
            destination = basis[row_offset + lo : row_offset + hi]
            destination[:] = 0.0
            destination[:, operator_columns] = (
                np.asarray(delta[lo:hi, source_columns], dtype=np.float64) - mean
            ) * multiplier
            enforce_memory_limit(max_rss_gib, f"{source_names} basis")
        row_offset += bank.shape[0]
        basis.flush()
    print(
        f"fit held-source whitener donors={source_names} basis={basis.shape} "
        f"peak_rss={peak_rss_gib():.2f} GiB",
        flush=True,
    )
    covariance = (basis.T @ basis) / total_rows
    del basis
    evals, evecs = np.linalg.eigh(covariance)
    del covariance
    evals = np.maximum(evals, float(evals.max()) * 1e-12)
    evecs32 = np.asarray(evecs, dtype=np.float32)
    np.savez(output_path, evals=evals, evecs=evecs32)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    try:
        basis_path.unlink()
    except FileNotFoundError:
        pass
    del evecs
    enforce_memory_limit(max_rss_gib, f"{source_names} whitener fit")
    return Operator(evals, evecs32, metadata, output_path, metadata_path), wcols


def apply_operator(
    matrix: np.ndarray,
    operator: Operator,
    wcols: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    output = np.asarray(matrix, dtype=np.float64).copy()
    before = float(np.median(np.linalg.norm(output[:, wcols], axis=1))) or 1.0
    weights = operator.evals ** (-WHITEN_GAMMA / 2.0)
    transformed = ((output[:, wcols] @ operator.evecs) * weights) @ operator.evecs.T
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    scale = before / after
    output[:, wcols] = transformed * scale
    return output, {
        "before_operator_median_norm": before,
        "after_operator_median_norm": after,
        "applied_scale": scale,
    }


def add_post_whiten_residual(
    baseline_post: np.ndarray,
    cd4: CD4Prepared,
    fraction: float,
    wcols: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    """Add unwhitened CD4, then restore baseline output norm by one scalar."""
    candidate = np.asarray(baseline_post, dtype=np.float64).copy()
    available = cd4.gene_votes > 0
    candidate[available] += fraction * cd4.mean_effect[available]
    baseline_norm = (
        float(np.median(np.linalg.norm(baseline_post[:, wcols], axis=1))) or 1.0
    )
    candidate_norm = (
        float(np.median(np.linalg.norm(candidate[:, wcols], axis=1))) or 1.0
    )
    scale = baseline_norm / candidate_norm
    candidate *= scale
    return candidate, {
        "fraction": fraction,
        "baseline_output_median_norm": baseline_norm,
        "candidate_pre_preserve_median_norm": candidate_norm,
        "post_add_global_preserve_scale": scale,
    }


def raw_truth(prepared: PreparedSource) -> np.ndarray:
    return prepared.effects * (prepared.raw_median_norm / EFFECT_SCALE)


def pds_per_target(
    prediction: np.ndarray,
    truth: np.ndarray,
    columns: np.ndarray,
) -> np.ndarray:
    predicted = np.asarray(prediction[:, columns], dtype=np.float64)
    observed = np.asarray(truth[:, columns], dtype=np.float64)
    predicted /= np.maximum(np.linalg.norm(predicted, axis=1, keepdims=True), 1e-12)
    observed /= np.maximum(np.linalg.norm(observed, axis=1, keepdims=True), 1e-12)
    similarities = predicted @ observed.T
    own = np.diag(similarities)[:, None]
    n = similarities.shape[0]
    ranks = (similarities > own).sum(axis=1)
    ranks = ranks + 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


def choose_candidate(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    grouped: dict[tuple[str, float], list[float]] = {}
    for row in rows:
        key = (str(row["mode"]), float(row["fraction"]))
        grouped.setdefault(key, []).append(float(row["delta_vs_baseline"]))
    if any(len(values) != len(DEVELOPMENT_HELD) for values in grouped.values()):
        raise ValueError("every candidate must have exactly two development folds")
    candidates = []
    for (mode, fraction), values in grouped.items():
        candidates.append(
            {
                "mode": mode,
                "fraction": fraction,
                "worst_fold_delta": float(min(values)),
                "mean_delta": float(np.mean(values)),
                "fold_deltas": values,
            }
        )
    selected = max(
        candidates,
        key=lambda item: (
            item["worst_fold_delta"],
            item["mean_delta"],
            -item["fraction"],
            item["mode"] == "gene_norm",
        ),
    )
    return {
        **selected,
        "selection_rule": (
            "maximize worst development-fold delta; then mean delta; then "
            "smaller fraction; then gene_norm"
        ),
        "candidate_grid": {
            "modes": list(MODES),
            "fractions": list(FRACTIONS),
        },
    }


def evaluate_fold(
    held: str,
    banks: dict[str, Bank],
    prepared: dict[str, PreparedSource],
    cd4_by_mode: dict[str, CD4Prepared],
    official: Sequence[str],
    panel: Sequence[str],
    targets: Sequence[str],
    whitener_path: Path,
    scratch: Path,
    candidates: Sequence[tuple[str, float]],
    *,
    chunk_rows: int,
    max_rss_gib: float,
    min_operator_genes: int,
    min_retrieval_genes: int,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    source_names = tuple(name for name in MAIN_NAMES if name != held)
    source_banks = tuple(banks[name] for name in source_names)
    prepared_sources = tuple(prepared[name] for name in source_names)
    baseline_pre, votes = build_baseline(prepared_sources)
    operator, wcols = fit_union_whitener(
        source_banks,
        prepared_sources,
        official,
        panel,
        whitener_path,
        scratch,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    if len(wcols) < min_operator_genes:
        raise ValueError(f"{held}: only {len(wcols)} operator genes")
    baseline_post, whitening = apply_operator(baseline_pre, operator, wcols)
    held_genes = set(banks[held].genes)
    retrieval_columns = np.asarray(
        [index for index in wcols if official[index] in held_genes],
        dtype=np.int64,
    )
    if len(retrieval_columns) < min_retrieval_genes:
        raise ValueError(f"{held}: only {len(retrieval_columns)} retrieval genes")
    if set(panel) & {official[index] for index in retrieval_columns}:
        raise AssertionError("panel coordinate entered retrieval")
    truth = raw_truth(prepared[held])
    baseline_per_target = pds_per_target(
        baseline_post, truth, retrieval_columns
    )
    baseline_pds = float(np.mean(baseline_per_target))
    candidate_rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    for target, value in zip(targets, baseline_per_target, strict=True):
        per_target_rows.append(
            {
                "held": held,
                "target_gene": target,
                "mode": "baseline",
                "fraction": 0.0,
                "pds": float(value),
                "delta_vs_baseline": 0.0,
            }
        )
    for mode, fraction in candidates:
        candidate, preserve = add_post_whiten_residual(
            baseline_post, cd4_by_mode[mode], fraction, wcols
        )
        per_target = pds_per_target(candidate, truth, retrieval_columns)
        score = float(np.mean(per_target))
        row = {
            "held": held,
            "mode": mode,
            "fraction": fraction,
            "pds": score,
            "baseline_pds": baseline_pds,
            "delta_vs_baseline": score - baseline_pds,
            "post_add_preserve": preserve,
        }
        candidate_rows.append(row)
        for target, value, base_value in zip(
            targets, per_target, baseline_per_target, strict=True
        ):
            per_target_rows.append(
                {
                    "held": held,
                    "target_gene": target,
                    "mode": mode,
                    "fraction": fraction,
                    "pds": float(value),
                    "delta_vs_baseline": float(value - base_value),
                }
            )
    fold = {
        "held": held,
        "source_donors": list(source_names),
        "baseline_pds": baseline_pds,
        "baseline_whitening": whitening,
        "source_vote_counts": {
            str(value): int(np.sum(votes == value))
            for value in sorted(set(votes.tolist()))
        },
        "operator_genes": len(wcols),
        "retrieval_genes": len(retrieval_columns),
        "operator_columns_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "retrieval_columns_sha256": hashlib.sha256(
            retrieval_columns.tobytes()
        ).hexdigest(),
        "whitener_path": str(whitener_path.resolve()),
        "whitener_sha256": sha256_file(whitener_path),
        "whitener_metadata_sha256": sha256_file(
            whitener_path.with_suffix(".json")
        ),
        "candidates": candidate_rows,
    }
    del truth, baseline_pre, baseline_post, operator
    gc.collect()
    shutil.rmtree(scratch, ignore_errors=True)
    enforce_memory_limit(max_rss_gib, f"{held} fold completion")
    return fold, per_target_rows


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    flattened: list[dict[str, object]] = []
    for row in rows:
        flattened.append(
            {
                key: json.dumps(value, sort_keys=True)
                if isinstance(value, (dict, list))
                else value
                for key, value in row.items()
            }
        )
    fieldnames = sorted({key for row in flattened for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(flattened)


def input_hashes(
    banks: dict[str, Bank],
    cd4_banks: Sequence[Bank],
    panel_csv: Path,
    genes_csv: Path,
) -> dict[str, object]:
    output: dict[str, object] = {
        "panel_csv": sha256_file(panel_csv),
        "genes_csv": sha256_file(genes_csv),
        "run_py": sha256_file(Path(__file__)),
        "banks": {},
    }
    for bank in [*banks.values(), *cd4_banks]:
        files = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt")
        if (bank.path / "meta.json").exists():
            files = (*files, "meta.json")
        output["banks"][bank.name] = {
            name: sha256_file(bank.path / name) for name in files
        }
    return output


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = read_csv_column(args.panel_csv, "target_gene")
    official = read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != args.expected_panel_targets:
        raise ValueError(f"expected {args.expected_panel_targets} panel targets")
    official_pos = {gene: index for index, gene in enumerate(official)}
    missing_panel = [gene for gene in panel if gene not in official_pos]
    if missing_panel:
        raise ValueError(f"panel targets absent from official axis: {missing_panel}")
    panel_positions = np.asarray(
        [official_pos[gene] for gene in panel], dtype=np.int64
    )
    panel_set = set(panel)
    banks = {
        "k562_gwps": Bank.load("k562_gwps", args.k562),
        "orion_hct116": Bank.load("orion_hct116", args.hct116),
        "orion_hek293t": Bank.load("orion_hek293t", args.hek293t),
    }
    cd4_banks = tuple(
        Bank.load(name, args.cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    all_banks = [*banks.values(), *cd4_banks]
    targets = tuple(
        target
        for target in panel
        if all(target in bank.pert_pos for bank in all_banks)
    )
    if len(targets) < args.min_shared_targets:
        raise ValueError(f"only {len(targets)} targets shared by all six banks")

    prepared = {
        name: prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name, bank in banks.items()
    }
    cd4_by_mode = {
        mode: prepare_cd4_mean(
            cd4_banks,
            targets,
            official,
            panel_set,
            mode,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for mode in MODES
    }
    coverage = cd4_by_mode["raw"].target_condition_count
    if any(
        not np.array_equal(value.target_condition_count, coverage)
        for value in cd4_by_mode.values()
    ):
        raise AssertionError("CD4 coverage changed between modes")
    if not np.all(coverage == len(CD4_NAMES)):
        raise AssertionError("shared gallery unexpectedly lacks a CD4 condition")

    output_dir = args.output.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    args.whitener_dir.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    grid = [(mode, fraction) for mode in MODES for fraction in FRACTIONS]

    development_folds: list[dict[str, object]] = []
    development_rows: list[dict[str, object]] = []
    development_per_target: list[dict[str, object]] = []
    for held in DEVELOPMENT_HELD:
        fold, per_target = evaluate_fold(
            held,
            banks,
            prepared,
            cd4_by_mode,
            official,
            panel,
            targets,
            args.whitener_dir / f"held_{held}.npz",
            args.scratch / f"held_{held}",
            grid,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
            min_operator_genes=args.min_operator_genes,
            min_retrieval_genes=args.min_retrieval_genes,
        )
        development_folds.append(fold)
        development_rows.extend(fold["candidates"])
        development_per_target.extend(per_target)
        print(
            f"development held={held} baseline={fold['baseline_pds']:.6f} "
            f"peak_rss={peak_rss_gib():.2f} GiB",
            flush=True,
        )

    selected = choose_candidate(development_rows)
    selected_pair = (str(selected["mode"]), float(selected["fraction"]))
    print(f"frozen selection={selected_pair}", flush=True)

    final_fold, final_per_target = evaluate_fold(
        FINAL_HELD,
        banks,
        prepared,
        cd4_by_mode,
        official,
        panel,
        targets,
        args.whitener_dir / f"held_{FINAL_HELD}.npz",
        args.scratch / f"held_{FINAL_HELD}",
        [selected_pair],
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
        min_operator_genes=args.min_operator_genes,
        min_retrieval_genes=args.min_retrieval_genes,
    )
    if len(final_fold["candidates"]) != 1:
        raise AssertionError("held K562 evaluated more than one candidate")
    final_candidate = final_fold["candidates"][0]

    result: dict[str, object] = {
        "experiment": "062-cd4-post-whiten",
        "status": "completed",
        "protocol": {
            "metric": "effect-level rank PDS",
            "main_donors": list(MAIN_NAMES),
            "development_held": list(DEVELOPMENT_HELD),
            "final_held_target_role": FINAL_HELD,
            "gallery_size": len(targets),
            "whitening_gamma": WHITEN_GAMMA,
            "whitening_fit": "source donors only in every fold",
            "baseline_renormalization": (
                "production pre/post operator median norm on source-union axis"
            ),
            "post_whiten_rule": (
                "baseline_post + fraction * unwhitened CD4 mean3 residual; "
                "one global scalar restores baseline output median norm"
            ),
            "truth": "raw centered, unnormalized, untransformed",
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
        },
        "selection": selected,
        "development_folds": development_folds,
        "k562_holdout": {
            **final_fold,
            "selected_mode": selected_pair[0],
            "selected_fraction": selected_pair[1],
            "selected_pds": final_candidate["pds"],
            "selected_delta_vs_baseline": final_candidate["delta_vs_baseline"],
            "number_of_candidates_evaluated": 1,
        },
        "axis_audit": {
            "official_genes": len(official),
            "challenge_targets": len(panel),
            "challenge_targets_excluded_count": len(panel_set),
            "panel_order_sha256": sha256_lines(panel),
            "target_order_sha256": sha256_lines(targets),
        },
        "leakage_audit": {
            "all_panel_coordinates_zeroed_in_main_predictions": True,
            "all_panel_coordinates_excluded_from_cd4": True,
            "all_panel_coordinates_excluded_from_every_operator": True,
            "all_panel_coordinates_excluded_from_every_retrieval": True,
            "truth_transformed": False,
            "truth_gene_normalized": False,
            "k562_target_role_scores_used_for_selection": False,
            "k562_candidate_grid_evaluated": False,
            "k562_selected_candidates_evaluated": 1,
            "k562_effects_used_as_development_source": True,
            "interpretation": (
                "three-line LODO target-role holdout, not a donor-independent "
                "outer split"
            ),
        },
        "cd4": {
            "conditions": list(CD4_NAMES),
            "modes": {
                mode: value.condition_details
                for mode, value in cd4_by_mode.items()
            },
        },
        "input_hashes": input_hashes(
            banks, cd4_banks, args.panel_csv, args.genes_csv
        ),
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
    write_csv(output_dir / "development_summary.csv", development_rows)
    write_csv(output_dir / "development_per_target_pds.csv", development_per_target)
    write_csv(output_dir / "k562_summary.csv", final_fold["candidates"])
    write_csv(output_dir / "k562_per_target_pds.csv", final_per_target)
    shutil.rmtree(args.scratch, ignore_errors=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument(
        "--whitener-dir", type=Path, default=HERE / "whiteners"
    )
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=55.0)
    parser.add_argument("--expected-panel-targets", type=int, default=300)
    parser.add_argument("--min-shared-targets", type=int, default=240)
    parser.add_argument("--min-operator-genes", type=int, default=17_000)
    parser.add_argument("--min-retrieval-genes", type=int, default=7_000)
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
    result = run(args)
    final = result["k562_holdout"]
    print(
        json.dumps(
            {
                "selection": result["selection"],
                "k562_baseline_pds": final["baseline_pds"],
                "k562_selected_pds": final["selected_pds"],
                "k562_delta": final["selected_delta_vs_baseline"],
                "peak_rss_gib": result["peak_rss_gib"],
            },
            indent=2,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
