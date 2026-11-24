#!/usr/bin/env python3
"""Leakage-safe, low-capacity per-target gate over experiment 064 arms.

The experiment deliberately works from compact target-level outcomes and
source-side features.  It never evaluates cell counts, runs cell-eval2, fits a
whitener, downloads data, or uses a held donor's perturbation effects as
features for that donor's outer fold.
"""

from __future__ import annotations

import os

# Keep BLAS deterministic and single-threaded before importing NumPy.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import argparse
import csv
import hashlib
import json
import math
import resource
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.lib import format as npformat


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
FOLDS = ("orion_hct116", "orion_hek293t")
ARMS = ("A", "B", "C", "D")
ALT_ARMS = ("A", "B", "C")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
MAX_RSS_GIB = 0.75
INNER_FOLDS = 5
LABEL_PERMUTATIONS = 500
ROUTING_PERMUTATIONS = 20_000
BOOTSTRAP_REPS = 20_000
SEED = 6_602_026

OUTCOME_PATH = (
    ROOT
    / "experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv"
)
EXP064_JSON = (
    ROOT / "experiments/064-pm-cd4-composition/cloud-output/results.json"
)
EXP058_CSV = ROOT / "experiments/058-perturbmatch/per_target_pds.csv"
EXP058_JSON = ROOT / "experiments/058-perturbmatch/results.json"
EXP062_DEV_CSV = (
    ROOT
    / "experiments/062-cd4-post-whiten/cloud-output/development_per_target_pds.csv"
)

PANEL_CSV = ROOT / "data/challenge/pert_counts.csv"
GENES_CSV = ROOT / "data/challenge/gene_names.csv"
K562_ROOT = ROOT / "data/donors/k562_gwps"
HCT116_ROOT = ROOT / "data/donors/orion_hct116"
HEK293T_ROOT = ROOT / "data/donors/orion_hek293t"
PERTURBMATCH_ROOT = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"


FEATURE_COLUMNS = (
    "source_agreement",
    "pm_incumbent_cosine",
    "pm_companion_advantage",
    "pm_slot_companion_advantage",
    "pm_shift_size",
    "pm_log_norm_ratio",
    "cd4_source_cosine",
    "cd4_condition_agreement",
    "cd4_log_norm",
    "cd4_breadth_log_mean",
    "cd4_quality_mean",
    "held_target_log_control",
    "held_target_control_rank",
    "held_minus_source_log_control",
)

FEATURE_FAMILIES = {
    "pm_core": (
        "source_agreement",
        "pm_incumbent_cosine",
        "pm_companion_advantage",
        "pm_slot_companion_advantage",
        "pm_shift_size",
        "pm_log_norm_ratio",
    ),
    "cd4_core": (
        "source_agreement",
        "cd4_source_cosine",
        "cd4_condition_agreement",
        "cd4_log_norm",
        "cd4_breadth_log_mean",
        "cd4_quality_mean",
    ),
    "context_core": (
        "source_agreement",
        "held_target_log_control",
        "held_target_control_rank",
        "held_minus_source_log_control",
    ),
    "joint_small": (
        "source_agreement",
        "pm_companion_advantage",
        "cd4_source_cosine",
        "cd4_breadth_log_mean",
        "held_target_log_control",
    ),
    "joint_all": FEATURE_COLUMNS,
}


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_rss(stage: str, limit: float = MAX_RSS_GIB) -> None:
    peak = peak_rss_gib()
    if peak >= limit:
        raise MemoryError(f"{stage}: peak RSS {peak:.6f} GiB >= {limit:.6f} GiB")


def sha256_file(path: Path, block_bytes: int = 4 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    rows = read_csv_rows(path)
    values = tuple(row[column].strip() for row in rows if row.get(column, "").strip())
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column}")
    return values


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@dataclass(frozen=True)
class NpyHeader:
    shape: tuple[int, ...]
    fortran_order: bool
    dtype: np.dtype
    offset: int


def read_npy_header(path: Path) -> NpyHeader:
    with path.open("rb") as handle:
        version = npformat.read_magic(handle)
        readers = {
            (1, 0): npformat.read_array_header_1_0,
            (2, 0): npformat.read_array_header_2_0,
            (3, 0): npformat.read_array_header_2_0,
        }
        if version not in readers:
            raise ValueError(f"{path}: unsupported npy version {version}")
        shape, fortran_order, dtype = readers[version](handle)
        return NpyHeader(
            shape=tuple(map(int, shape)),
            fortran_order=bool(fortran_order),
            dtype=np.dtype(dtype),
            offset=handle.tell(),
        )


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    delta_path: Path
    b_ctrl_path: Path | None
    header: NpyHeader

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        path = path.resolve()
        genes = read_lines(path / "genes.txt", require_unique=False)
        perts = read_lines(path / "perts.txt")
        header = read_npy_header(path / "delta.npy")
        if header.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {header.shape} != {(len(perts), len(genes))}"
            )
        b_ctrl = path / "b_ctrl.npy"
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos={gene: index for index, gene in enumerate(genes)},
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=path / "delta.npy",
            b_ctrl_path=b_ctrl if b_ctrl.exists() else None,
            header=header,
        )


def _iter_c_selected(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int = 8,
) -> Iterable[tuple[int, np.ndarray]]:
    """Stream selected columns from a C-order npy without retaining mmap pages."""
    if bank.header.fortran_order:
        raise ValueError(f"{bank.name}: expected C-order matrix")
    n_rows, n_cols = bank.header.shape
    row_bytes = n_cols * bank.header.dtype.itemsize
    with bank.delta_path.open("rb") as handle:
        handle.seek(bank.header.offset)
        for lo in range(0, n_rows, chunk_rows):
            count = min(chunk_rows, n_rows - lo)
            raw = handle.read(count * row_bytes)
            if len(raw) != count * row_bytes:
                raise EOFError(f"{bank.delta_path}: truncated at row {lo}")
            whole = np.frombuffer(raw, dtype=bank.header.dtype).reshape(count, n_cols)
            selected = np.asarray(whole[:, columns], dtype=np.float64)
            yield lo, selected


def _fortran_selected_columns(
    bank: Bank,
    columns: np.ndarray,
) -> Iterable[tuple[int, np.ndarray]]:
    """Stream whole selected columns from a Fortran-order npy."""
    if not bank.header.fortran_order:
        raise ValueError(f"{bank.name}: expected Fortran-order matrix")
    n_rows, _ = bank.header.shape
    column_bytes = n_rows * bank.header.dtype.itemsize
    with bank.delta_path.open("rb") as handle:
        for destination, source_column in enumerate(columns.tolist()):
            handle.seek(bank.header.offset + source_column * column_bytes)
            raw = handle.read(column_bytes)
            if len(raw) != column_bytes:
                raise EOFError(
                    f"{bank.delta_path}: truncated at column {source_column}"
                )
            yield destination, np.asarray(
                np.frombuffer(raw, dtype=bank.header.dtype), dtype=np.float64
            )


def prepare_effects(
    bank: Bank,
    targets: Sequence[str],
    axis: Sequence[str],
    *,
    gene_normalize: bool,
) -> np.ndarray:
    """Return globally centered target rows on a compact common axis.

    C-order matrices are streamed in tiny row blocks through normal reads.
    This intentionally avoids retaining a multi-gigabyte file-backed mmap in
    the process RSS.  PerturbMatch's Fortran-order matrix is streamed by
    contiguous columns.
    """
    columns = np.asarray([bank.gene_pos[gene] for gene in axis], dtype=np.int64)
    target_rows = np.asarray([bank.pert_pos[target] for target in targets], dtype=np.int64)
    n_rows = bank.header.shape[0]
    means = np.zeros(len(columns), dtype=np.float64)

    if not bank.header.fortran_order:
        for _, block in _iter_c_selected(bank, columns):
            means += block.sum(axis=0, dtype=np.float64)
        means /= n_rows
        absolute = np.zeros(len(columns), dtype=np.float64) if gene_normalize else None
        for _, block in _iter_c_selected(bank, columns):
            centered = block - means
            if absolute is not None:
                absolute += np.abs(centered).sum(axis=0, dtype=np.float64)
        target_values = np.empty((len(targets), len(axis)), dtype=np.float64)
        n_cols = bank.header.shape[1]
        row_bytes = n_cols * bank.header.dtype.itemsize
        with bank.delta_path.open("rb") as handle:
            for destination, row in enumerate(target_rows.tolist()):
                handle.seek(bank.header.offset + row * row_bytes)
                raw = handle.read(row_bytes)
                whole = np.frombuffer(raw, dtype=bank.header.dtype, count=n_cols)
                target_values[destination] = np.asarray(
                    whole[columns], dtype=np.float64
                )
        target_values -= means
        if absolute is not None:
            gene_scale = absolute / n_rows
            live = gene_scale[gene_scale > 0]
            floor = float(np.percentile(live, 25.0)) if live.size else 1.0
            gene_scale = np.maximum(gene_scale, floor)
            target_values /= gene_scale
    else:
        target_values = np.empty((len(targets), len(axis)), dtype=np.float64)
        absolute = np.zeros(len(columns), dtype=np.float64) if gene_normalize else None
        for destination, column in _fortran_selected_columns(bank, columns):
            mean = float(column.mean())
            means[destination] = mean
            centered = column - mean
            target_values[:, destination] = centered[target_rows]
            if absolute is not None:
                absolute[destination] = np.abs(centered).mean()
        if absolute is not None:
            live = absolute[absolute > 0]
            floor = float(np.percentile(live, 25.0)) if live.size else 1.0
            target_values /= np.maximum(absolute, floor)

    norms = np.linalg.norm(target_values, axis=1)
    median_norm = float(np.median(norms[norms > 0])) if np.any(norms > 0) else 1.0
    target_values /= max(median_norm, 1e-12)
    enforce_rss(f"prepared {bank.name}")
    return target_values.astype(np.float32)


def control_target_values(bank: Bank, targets: Sequence[str]) -> np.ndarray:
    if bank.b_ctrl_path is None:
        raise ValueError(f"{bank.name}: no control vector")
    control = np.load(bank.b_ctrl_path)
    values = np.asarray(
        [
            control[bank.gene_pos[target]]
            if target in bank.gene_pos
            else float("nan")
            for target in targets
        ],
        dtype=np.float64,
    )
    finite = np.isfinite(values)
    logged = np.full(len(values), np.nan, dtype=np.float64)
    logged[finite] = np.log1p(np.maximum(values[finite], 0.0))
    return logged


def _as_bool(value: str) -> float:
    return float(str(value).strip().lower() in {"1", "true", "yes"})


def rank01(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + 1 + end)
        start = end
    return ranks / max(len(values), 1)


def cd4_metadata_features(
    cd4_banks: dict[str, Bank],
    targets: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    breadth_rows: list[np.ndarray] = []
    quality_rows: list[np.ndarray] = []
    for name in CD4_NAMES:
        rows = read_csv_rows(cd4_banks[name].path / "obs.csv")
        by_target = {row["target_gene"]: row for row in rows}
        ordered = [by_target[target] for target in targets]
        breadth = np.asarray(
            [float(row.get("n_total_de_genes") or 0.0) for row in ordered],
            dtype=np.float64,
        )
        breadth_rows.append(np.log1p(np.maximum(breadth, 0.0)))
        parts = [
            np.asarray([_as_bool(row.get("ontarget_significant", "")) for row in ordered]),
            1.0
            - np.asarray([_as_bool(row.get("low_target_gex", "")) for row in ordered]),
            1.0
            - np.asarray(
                [_as_bool(row.get("distal_offtarget_flag", "")) for row in ordered]
            ),
            1.0
            - np.asarray(
                [_as_bool(row.get("neighboring_gene_KD", "")) for row in ordered]
            ),
        ]
        for column in ("n_guides", "n_cells_target"):
            values = np.asarray(
                [float(row.get(column) or 0.0) for row in ordered], dtype=np.float64
            )
            parts.append(rank01(np.log1p(np.maximum(values, 0.0))))
        quality_rows.append(np.mean(np.vstack(parts), axis=0))
    return np.mean(np.vstack(breadth_rows), axis=0), np.mean(
        np.vstack(quality_rows), axis=0
    )


def row_unit(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    return values / np.maximum(norms, 1e-12)


def row_cos(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.sum(row_unit(left) * row_unit(right), axis=1)


def build_features() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    panel = read_csv_column(PANEL_CSV, "target_gene")
    official = read_csv_column(GENES_CSV, "gene_name")
    panel_set = set(panel)
    banks = {
        "k562_gwps": Bank.load("k562_gwps", K562_ROOT),
        "orion_hct116": Bank.load("orion_hct116", HCT116_ROOT),
        "orion_hek293t": Bank.load("orion_hek293t", HEK293T_ROOT),
        "perturbmatch": Bank.load("perturbmatch", PERTURBMATCH_ROOT),
        **{
            name: Bank.load(name, CD4_ROOT / name)
            for name in CD4_NAMES
        },
    }
    targets = tuple(
        target
        for target in panel
        if all(target in bank.pert_pos for bank in banks.values())
    )
    exp064 = json.loads(EXP064_JSON.read_text())
    if list(targets) != exp064["axis"]["target_order"]:
        raise AssertionError("bank-derived target order does not match experiment 064")
    if sha256_lines(targets) != exp064["axis"]["target_order_sha256"]:
        raise AssertionError("experiment 064 target-order hash mismatch")

    gene_sets = [set(bank.genes) for bank in banks.values()]
    axis = tuple(
        gene
        for gene in official
        if gene not in panel_set and all(gene in genes for genes in gene_sets)
    )
    if len(axis) < 6_000:
        raise ValueError(f"common inference feature axis unexpectedly small: {len(axis)}")

    effects = {
        name: prepare_effects(
            bank,
            targets,
            axis,
            gene_normalize=name in CD4_NAMES,
        )
        for name, bank in banks.items()
    }
    controls = {
        name: control_target_values(banks[name], targets)
        for name in ("k562_gwps", *FOLDS)
    }
    breadth, quality = cd4_metadata_features(banks, targets)

    k562 = effects["k562_gwps"]
    perturbmatch = effects["perturbmatch"]
    cd4_conditions = [effects[name] for name in CD4_NAMES]
    cd4 = np.mean(np.stack(cd4_conditions, axis=0), axis=0)
    cd4_condition_agreement = np.mean(
        np.vstack(
            [
                row_cos(cd4_conditions[0], cd4_conditions[1]),
                row_cos(cd4_conditions[0], cd4_conditions[2]),
                row_cos(cd4_conditions[1], cd4_conditions[2]),
            ]
        ),
        axis=0,
    )
    k_unit = row_unit(k562)
    pm_unit = row_unit(perturbmatch)
    pm_slot = row_unit(k_unit + pm_unit)
    cd4_unit = row_unit(cd4)
    target_index = {target: index for index, target in enumerate(targets)}

    rows: list[dict[str, Any]] = []
    dependencies: dict[str, dict[str, Any]] = {}
    for held in FOLDS:
        companion = next(name for name in FOLDS if name != held)
        companion_values = effects[companion]
        companion_unit = row_unit(companion_values)
        incumbent_source = row_unit(k_unit + companion_unit)
        pm_source = row_unit(pm_slot + companion_unit)
        held_control = controls[held].copy()
        if not np.isfinite(held_control).all():
            held_control[~np.isfinite(held_control)] = np.nanmedian(held_control)
        source_stack = np.vstack([controls["k562_gwps"], controls[companion]])
        with np.errstate(invalid="ignore"):
            source_control = np.nanmean(source_stack, axis=0)
        source_control[~np.isfinite(source_control)] = held_control[
            ~np.isfinite(source_control)
        ]
        held_rank = rank01(held_control)

        values = {
            "source_agreement": np.sum(k_unit * companion_unit, axis=1),
            "pm_incumbent_cosine": np.sum(pm_unit * k_unit, axis=1),
            "pm_companion_advantage": (
                np.sum(pm_unit * companion_unit, axis=1)
                - np.sum(k_unit * companion_unit, axis=1)
            ),
            "pm_slot_companion_advantage": (
                np.sum(pm_slot * companion_unit, axis=1)
                - np.sum(k_unit * companion_unit, axis=1)
            ),
            "pm_shift_size": 1.0 - np.sum(incumbent_source * pm_source, axis=1),
            "pm_log_norm_ratio": np.log(
                (np.linalg.norm(perturbmatch, axis=1) + 1e-12)
                / (np.linalg.norm(k562, axis=1) + 1e-12)
            ),
            "cd4_source_cosine": np.sum(cd4_unit * pm_source, axis=1),
            "cd4_condition_agreement": cd4_condition_agreement,
            "cd4_log_norm": np.log(np.linalg.norm(cd4, axis=1) + 1e-12),
            "cd4_breadth_log_mean": breadth,
            "cd4_quality_mean": quality,
            "held_target_log_control": held_control,
            "held_target_control_rank": held_rank,
            "held_minus_source_log_control": held_control - source_control,
        }
        for column in FEATURE_COLUMNS:
            if not np.isfinite(values[column]).all():
                raise ValueError(f"{held}: nonfinite feature {column}")
        for target in targets:
            index = target_index[target]
            rows.append(
                {
                    "held": held,
                    "target_gene": target,
                    **{column: float(values[column][index]) for column in FEATURE_COLUMNS},
                }
            )
        dependencies[held] = {
            "effect_banks": [
                "k562_gwps",
                "perturbmatch",
                companion,
                *CD4_NAMES,
            ],
            "held_effect_bank_excluded": held,
            "held_input_used": f"{held}/b_ctrl.npy target-gene coordinate only",
        }

    del effects
    enforce_rss("feature extraction complete")
    return rows, {
        "targets": len(targets),
        "target_order_sha256": sha256_lines(targets),
        "feature_axis_genes": len(axis),
        "feature_axis_sha256": sha256_lines(axis),
        "dependencies": dependencies,
        "features": list(FEATURE_COLUMNS),
    }


def load_outcomes(targets: Sequence[str]) -> dict[str, dict[str, dict[str, float]]]:
    rows = read_csv_rows(OUTCOME_PATH)
    result: dict[str, dict[str, dict[str, float]]] = {held: {} for held in FOLDS}
    for row in rows:
        held = row["held"]
        target = row["target_gene"]
        result[held][target] = {
            arm: float(row[f"pds_{arm}"]) for arm in ARMS
        }
    for held in FOLDS:
        if tuple(result[held]) != tuple(targets):
            raise AssertionError(f"{held}: outcome target order mismatch")
    return result


def feature_matrices(
    rows: Sequence[dict[str, Any]],
    targets: Sequence[str],
) -> dict[str, np.ndarray]:
    by_key = {(row["held"], row["target_gene"]): row for row in rows}
    return {
        held: np.asarray(
            [
                [float(by_key[(held, target)][column]) for column in FEATURE_COLUMNS]
                for target in targets
            ],
            dtype=np.float64,
        )
        for held in FOLDS
    }


def outcome_matrices(
    outcomes: dict[str, dict[str, dict[str, float]]],
    targets: Sequence[str],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    pds = {
        held: np.asarray(
            [[outcomes[held][target][arm] for arm in ARMS] for target in targets],
            dtype=np.float64,
        )
        for held in FOLDS
    }
    arm_pos = {arm: index for index, arm in enumerate(ARMS)}
    gains = {
        held: pds[held][:, [arm_pos[arm] for arm in ALT_ARMS]]
        - pds[held][:, [arm_pos["D"]]]
        for held in FOLDS
    }
    return pds, gains


def pearson(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    left = left - left.mean()
    right = right - right.mean()
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator else float("nan")


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    return pearson(rank01(np.asarray(left)), rank01(np.asarray(right)))


@dataclass(frozen=True)
class Candidate:
    kind: str
    name: str
    features: tuple[str, ...] = ()
    alternatives: tuple[str, ...] = ()
    alpha: float = 0.0
    margin: float = 0.0
    stump_feature: str = ""
    stump_alternative: str = ""

    @property
    def complexity(self) -> int:
        if self.kind == "always_d":
            return 0
        if self.kind == "constant":
            return 1
        if self.kind == "stump":
            return 2
        return len(self.features) * max(len(self.alternatives), 1)


def candidate_grid() -> list[Candidate]:
    candidates = [Candidate(kind="always_d", name="always_D")]
    for arm in ALT_ARMS:
        candidates.append(
            Candidate(kind="constant", name=f"always_{arm}", alternatives=(arm,))
        )
    for feature in FEATURE_COLUMNS:
        for arm in ALT_ARMS:
            candidates.append(
                Candidate(
                    kind="stump",
                    name=f"stump:{arm}:{feature}",
                    stump_feature=feature,
                    stump_alternative=arm,
                    alternatives=(arm,),
                )
            )
    ridge_families = (
        ("pm_drop", FEATURE_FAMILIES["pm_core"], ("C",)),
        ("cd4_drop", FEATURE_FAMILIES["cd4_core"], ("B",)),
        ("context_drop_both", FEATURE_FAMILIES["context_core"], ("A",)),
        ("joint_small", FEATURE_FAMILIES["joint_small"], ALT_ARMS),
        ("joint_all", FEATURE_FAMILIES["joint_all"], ALT_ARMS),
    )
    for family, features, alternatives in ridge_families:
        for alpha in (1.0, 10.0, 100.0, 1_000.0):
            for margin in (0.0, 0.0025, 0.005, 0.01):
                candidates.append(
                    Candidate(
                        kind="ridge",
                        name=f"ridge:{family}:a={alpha:g}:m={margin:g}",
                        features=tuple(features),
                        alternatives=tuple(alternatives),
                        alpha=alpha,
                        margin=margin,
                    )
                )
    return candidates


def arm_gain(arms: np.ndarray, gains: np.ndarray) -> np.ndarray:
    result = np.zeros(len(arms), dtype=np.float64)
    for index, arm in enumerate(ALT_ARMS):
        mask = arms == arm
        result[mask] = gains[mask, index]
    return result


def fit_candidate(
    candidate: Candidate,
    x_train: np.ndarray,
    gains_train: np.ndarray,
) -> dict[str, Any]:
    feature_pos = {feature: index for index, feature in enumerate(FEATURE_COLUMNS)}
    if candidate.kind in {"always_d", "constant"}:
        return {}
    if candidate.kind == "stump":
        column = feature_pos[candidate.stump_feature]
        values = x_train[:, column]
        y = gains_train[:, ALT_ARMS.index(candidate.stump_alternative)]
        thresholds = np.unique(np.quantile(values, np.linspace(0.1, 0.9, 9)))
        best = {
            "score": 0.0,
            "threshold": float(np.median(values)),
            "direction": "none",
        }
        for threshold in thresholds:
            for direction in ("low", "high"):
                mask = values <= threshold if direction == "low" else values >= threshold
                if mask.sum() < 15 or (~mask).sum() < 15:
                    continue
                score = float(np.sum(y[mask]) / len(y))
                key = (score, -int(mask.sum()), direction)
                current = (
                    float(best["score"]),
                    -int(np.sum(values <= best["threshold"]))
                    if best["direction"] == "low"
                    else -int(np.sum(values >= best["threshold"]))
                    if best["direction"] == "high"
                    else 0,
                    str(best["direction"]),
                )
                if key > current:
                    best = {
                        "score": score,
                        "threshold": float(threshold),
                        "direction": direction,
                    }
        return best
    if candidate.kind == "ridge":
        columns = np.asarray(
            [feature_pos[feature] for feature in candidate.features], dtype=np.int64
        )
        x = x_train[:, columns]
        mean = x.mean(axis=0)
        scale = x.std(axis=0)
        scale[scale < 1e-12] = 1.0
        z = (x - mean) / scale
        alt_indices = np.asarray(
            [ALT_ARMS.index(arm) for arm in candidate.alternatives], dtype=np.int64
        )
        y = gains_train[:, alt_indices]
        y_mean = y.mean(axis=0)
        gram = z.T @ z + candidate.alpha * np.eye(z.shape[1])
        coefficient = np.linalg.solve(gram, z.T @ (y - y_mean))
        return {
            "columns": columns,
            "mean": mean,
            "scale": scale,
            "y_mean": y_mean,
            "coefficient": coefficient,
        }
    raise ValueError(f"unknown candidate kind {candidate.kind}")


def predict_candidate(
    candidate: Candidate,
    model: dict[str, Any],
    x_test: np.ndarray,
) -> np.ndarray:
    if candidate.kind == "always_d":
        return np.full(len(x_test), "D", dtype="<U1")
    if candidate.kind == "constant":
        return np.full(len(x_test), candidate.alternatives[0], dtype="<U1")
    if candidate.kind == "stump":
        arms = np.full(len(x_test), "D", dtype="<U1")
        if model["direction"] == "none":
            return arms
        column = FEATURE_COLUMNS.index(candidate.stump_feature)
        values = x_test[:, column]
        mask = (
            values <= model["threshold"]
            if model["direction"] == "low"
            else values >= model["threshold"]
        )
        arms[mask] = candidate.stump_alternative
        return arms
    if candidate.kind == "ridge":
        z = (x_test[:, model["columns"]] - model["mean"]) / model["scale"]
        predicted = model["y_mean"] + z @ model["coefficient"]
        best = np.argmax(predicted, axis=1)
        best_gain = predicted[np.arange(len(predicted)), best]
        arms = np.full(len(x_test), "D", dtype="<U1")
        alternatives = np.asarray(candidate.alternatives)
        switch = best_gain > candidate.margin
        arms[switch] = alternatives[best[switch]]
        return arms
    raise ValueError(f"unknown candidate kind {candidate.kind}")


def summarize_model(candidate: Candidate, model: dict[str, Any]) -> dict[str, Any]:
    if candidate.kind in {"always_d", "constant"}:
        return {"kind": candidate.kind}
    if candidate.kind == "stump":
        return {
            "kind": "stump",
            "feature": candidate.stump_feature,
            "alternative": candidate.stump_alternative,
            "threshold": float(model["threshold"]),
            "direction": str(model["direction"]),
        }
    if candidate.kind == "ridge":
        coefficient = np.asarray(model["coefficient"])
        return {
            "kind": "ridge",
            "features": list(candidate.features),
            "alternatives": list(candidate.alternatives),
            "alpha": candidate.alpha,
            "margin": candidate.margin,
            "feature_mean": np.asarray(model["mean"]).tolist(),
            "feature_scale": np.asarray(model["scale"]).tolist(),
            "gain_intercept": np.asarray(model["y_mean"]).tolist(),
            "coefficient": coefficient.tolist(),
        }
    raise ValueError(f"unknown candidate kind {candidate.kind}")


def target_fold_ids(targets: Sequence[str], folds: int = INNER_FOLDS) -> np.ndarray:
    return np.asarray(
        [
            int(hashlib.sha256(target.encode()).hexdigest()[:12], 16) % folds
            for target in targets
        ],
        dtype=np.int64,
    )


def nested_select(
    x_train: np.ndarray,
    gains_train: np.ndarray,
    targets: Sequence[str],
    candidates: Sequence[Candidate],
    *,
    require_nontrivial: bool = False,
) -> tuple[Candidate, dict[str, Any], list[dict[str, Any]], np.ndarray]:
    fold_ids = target_fold_ids(targets)
    evaluations: list[dict[str, Any]] = []
    candidate_oof: dict[str, np.ndarray] = {}
    for candidate in candidates:
        if require_nontrivial and candidate.kind == "always_d":
            continue
        oof_arms = np.full(len(targets), "D", dtype="<U1")
        for inner in range(INNER_FOLDS):
            validation = fold_ids == inner
            training = ~validation
            model = fit_candidate(
                candidate, x_train[training], gains_train[training]
            )
            oof_arms[validation] = predict_candidate(
                candidate, model, x_train[validation]
            )
        realized = arm_gain(oof_arms, gains_train)
        switches = int(np.sum(oof_arms != "D"))
        evaluations.append(
            {
                "candidate": candidate.name,
                "kind": candidate.kind,
                "complexity": candidate.complexity,
                "oof_gain": float(realized.mean()),
                "oof_switches": switches,
            }
        )
        candidate_oof[candidate.name] = oof_arms
    evaluations.sort(
        key=lambda row: (
            -float(row["oof_gain"]),
            int(row["complexity"]),
            int(row["oof_switches"]),
            str(row["candidate"]),
        )
    )
    selected_row = evaluations[0]
    selected = next(
        candidate for candidate in candidates if candidate.name == selected_row["candidate"]
    )
    model = fit_candidate(selected, x_train, gains_train)
    return selected, model, evaluations, candidate_oof[selected.name]


def arm_counts(arms: np.ndarray) -> dict[str, int]:
    return {arm: int(np.sum(arms == arm)) for arm in ARMS}


def bootstrap_cluster(
    per_fold_gain: dict[str, np.ndarray],
    *,
    reps: int,
    seed: int,
) -> dict[str, Any]:
    cluster = np.mean(np.stack([per_fold_gain[held] for held in FOLDS], axis=1), axis=1)
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    for lo in range(0, reps, 1_000):
        hi = min(lo + 1_000, reps)
        sample = rng.integers(0, len(cluster), size=(hi - lo, len(cluster)))
        draws[lo:hi] = cluster[sample].mean(axis=1)
    return {
        "estimate": float(cluster.mean()),
        "ci95": [
            float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)),
        ],
        "positive_probability": float(np.mean(draws > 0)),
        "target_clusters": len(cluster),
        "reps": reps,
        "seed": seed,
    }


def matched_routing_null(
    selected_arms: dict[str, np.ndarray],
    gains: dict[str, np.ndarray],
    *,
    reps: int,
    seed: int,
    observed: float,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    for replicate in range(reps):
        fold_values = []
        for held in FOLDS:
            permuted = selected_arms[held][rng.permutation(len(selected_arms[held]))]
            fold_values.append(float(arm_gain(permuted, gains[held]).mean()))
        draws[replicate] = float(np.mean(fold_values))
    return {
        "method": "permute predicted arm assignments within each outer test fold",
        "reps": reps,
        "seed": seed,
        "mean": float(draws.mean()),
        "p95": float(np.percentile(draws, 95)),
        "p99": float(np.percentile(draws, 99)),
        "one_sided_p": float((1 + np.sum(draws >= observed)) / (reps + 1)),
    }


def label_permutation_null(
    features: dict[str, np.ndarray],
    gains: dict[str, np.ndarray],
    targets: Sequence[str],
    candidates: Sequence[Candidate],
    *,
    reps: int,
    seed: int,
    observed: float,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    for replicate in range(reps):
        fold_values = []
        for test_held in FOLDS:
            train_held = next(held for held in FOLDS if held != test_held)
            permutation = rng.permutation(len(targets))
            permuted_train = gains[train_held][permutation]
            selected, model, _, _ = nested_select(
                features[train_held],
                permuted_train,
                targets,
                candidates,
            )
            arms = predict_candidate(selected, model, features[test_held])
            fold_values.append(float(arm_gain(arms, gains[test_held]).mean()))
        draws[replicate] = float(np.mean(fold_values))
    return {
        "method": (
            "permute training target labels jointly across A/B/C gains and rerun "
            "the complete inner-CV candidate selection in both outer directions"
        ),
        "reps": reps,
        "seed": seed,
        "mean": float(draws.mean()),
        "p95": float(np.percentile(draws, 95)),
        "p99": float(np.percentile(draws, 99)),
        "one_sided_p": float((1 + np.sum(draws >= observed)) / (reps + 1)),
    }


def oracle_and_transfer(
    pds: dict[str, np.ndarray],
    gains: dict[str, np.ndarray],
) -> dict[str, Any]:
    result: dict[str, Any] = {"folds": {}}
    for held in FOLDS:
        fold: dict[str, Any] = {"fixed_D": float(pds[held][:, 3].mean())}
        subsets = {
            "D_or_B": ("D", "B"),
            "D_or_C": ("D", "C"),
            "D_or_A": ("D", "A"),
            "four_arm": ARMS,
        }
        for name, arms in subsets.items():
            indices = np.asarray([ARMS.index(arm) for arm in arms], dtype=np.int64)
            local = np.argmax(pds[held][:, indices], axis=1)
            chosen = np.asarray(arms)[local]
            score = float(pds[held][np.arange(len(chosen)), indices[local]].mean())
            fold[name] = {
                "pds": score,
                "gain_vs_D": score - fold["fixed_D"],
                "arm_counts": arm_counts(chosen),
            }
        result["folds"][held] = fold

    transfer_rows = []
    transfer_gains = []
    for test_held in FOLDS:
        train_held = next(held for held in FOLDS if held != test_held)
        chosen_index = np.argmax(pds[train_held], axis=1)
        chosen = np.asarray(ARMS)[chosen_index]
        test_gain = arm_gain(chosen, gains[test_held])
        transfer_rows.append(
            {
                "train_fold": train_held,
                "test_fold": test_held,
                "gain_vs_D": float(test_gain.mean()),
                "arm_counts": arm_counts(chosen),
            }
        )
        transfer_gains.append(float(test_gain.mean()))
    result["same_target_other_donor_lookup"] = {
        "capacity": "221 target-specific arm choices; diagnostic only",
        "folds": transfer_rows,
        "macro_gain_vs_D": float(np.mean(transfer_gains)),
    }
    return result


def cross_fold_stability(
    pds: dict[str, np.ndarray],
    gains: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    left, right = FOLDS
    rows: list[dict[str, Any]] = []
    names = {
        "A_minus_D": gains[left][:, 0],
        "B_minus_D": gains[left][:, 1],
        "C_minus_D": gains[left][:, 2],
    }
    right_values = {
        "A_minus_D": gains[right][:, 0],
        "B_minus_D": gains[right][:, 1],
        "C_minus_D": gains[right][:, 2],
    }
    for name, values in names.items():
        other = right_values[name]
        rows.append(
            {
                "quantity": name,
                "pearson": pearson(values, other),
                "spearman": spearman(values, other),
                "same_sign_fraction": float(np.mean(np.sign(values) == np.sign(other))),
            }
        )
    for name, left_values, right_values_array in (
        ("B_minus_A", pds[left][:, 1] - pds[left][:, 0], pds[right][:, 1] - pds[right][:, 0]),
        ("C_minus_A", pds[left][:, 2] - pds[left][:, 0], pds[right][:, 2] - pds[right][:, 0]),
        ("D_minus_A", pds[left][:, 3] - pds[left][:, 0], pds[right][:, 3] - pds[right][:, 0]),
    ):
        rows.append(
            {
                "quantity": name,
                "pearson": pearson(left_values, right_values_array),
                "spearman": spearman(left_values, right_values_array),
                "same_sign_fraction": float(
                    np.mean(np.sign(left_values) == np.sign(right_values_array))
                ),
            }
        )
    return rows


def compact_output_consistency(
    outcomes: dict[str, dict[str, dict[str, float]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    order058 = json.loads(EXP058_JSON.read_text())["axis"]["target_order"]
    values058: dict[tuple[str, str, str], float] = {}
    for row in read_csv_rows(EXP058_CSV):
        if row["geometry"] != "source_whitened_g0.65" or row["mode"] != "raw":
            continue
        if row["candidate"] not in {
            "baseline_old_k562_plus_orion",
            "full_ensemble_f0.50",
        }:
            continue
        target = order058[int(row["target_index"])]
        values058[(row["held"], target, row["candidate"])] = float(row["pds"])
    for held in FOLDS:
        shared = [
            target
            for target in outcomes[held]
            if (held, target, "baseline_old_k562_plus_orion") in values058
            and (held, target, "full_ensemble_f0.50") in values058
        ]
        prior = np.asarray(
            [
                values058[(held, target, "full_ensemble_f0.50")]
                - values058[(held, target, "baseline_old_k562_plus_orion")]
                for target in shared
            ]
        )
        current = np.asarray(
            [
                outcomes[held][target]["B"] - outcomes[held][target]["A"]
                for target in shared
            ]
        )
        rows.append(
            {
                "source_experiment": "058",
                "held": held,
                "comparison": "PM gain: exp058 f0.50 vs exp064 B-A",
                "targets": len(shared),
                "prior_mean": float(prior.mean()),
                "current_mean": float(current.mean()),
                "pearson": pearson(prior, current),
                "spearman": spearman(prior, current),
            }
        )

    values062: dict[tuple[str, str, str], float] = {}
    for row in read_csv_rows(EXP062_DEV_CSV):
        if row["mode"] == "baseline":
            key = "baseline"
        elif row["mode"] == "gene_norm" and math.isclose(float(row["fraction"]), 0.05):
            key = "candidate"
        else:
            continue
        values062[(row["held"], row["target_gene"], key)] = float(row["pds"])
    for held in FOLDS:
        shared = [
            target
            for target in outcomes[held]
            if (held, target, "baseline") in values062
            and (held, target, "candidate") in values062
        ]
        prior = np.asarray(
            [
                values062[(held, target, "candidate")]
                - values062[(held, target, "baseline")]
                for target in shared
            ]
        )
        current = np.asarray(
            [
                outcomes[held][target]["C"] - outcomes[held][target]["A"]
                for target in shared
            ]
        )
        rows.append(
            {
                "source_experiment": "062",
                "held": held,
                "comparison": "CD4 gain: exp062 selected vs exp064 C-A",
                "targets": len(shared),
                "prior_mean": float(prior.mean()),
                "current_mean": float(current.mean()),
                "pearson": pearson(prior, current),
                "spearman": spearman(prior, current),
            }
        )
    return rows


def feature_associations(
    features: dict[str, np.ndarray],
    gains: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for held in FOLDS:
        for alternative_index, arm in enumerate(ALT_ARMS):
            y = gains[held][:, alternative_index]
            for feature_index, feature in enumerate(FEATURE_COLUMNS):
                x = features[held][:, feature_index]
                rows.append(
                    {
                        "held": held,
                        "alternative": arm,
                        "gain_definition": f"{arm}_minus_D",
                        "feature": feature,
                        "pearson": pearson(x, y),
                        "spearman": spearman(x, y),
                    }
                )
    return rows


def synthetic_tests() -> dict[str, Any]:
    targets = tuple(f"T{index:03d}" for index in range(100))
    fold_ids = target_fold_ids(targets)
    if set(fold_ids.tolist()) != set(range(INNER_FOLDS)):
        raise AssertionError("inner target folds are not all represented")
    if any(token in feature.lower() for feature in FEATURE_COLUMNS for token in ("pds", "truth")):
        raise AssertionError("an outcome-like field entered the feature list")
    rng = np.random.default_rng(7)
    x = rng.normal(size=(100, len(FEATURE_COLUMNS)))
    gains = np.zeros((100, 3), dtype=np.float64)
    gains[:, 2] = np.where(x[:, 0] > 0, 0.02, -0.02)
    candidate = Candidate(
        kind="stump",
        name="synthetic",
        stump_feature=FEATURE_COLUMNS[0],
        stump_alternative="C",
        alternatives=("C",),
    )
    model = fit_candidate(candidate, x[:80], gains[:80])
    arms = predict_candidate(candidate, model, x[80:])
    if float(arm_gain(arms, gains[80:]).mean()) <= 0:
        raise AssertionError("stump failed to recover a synthetic transferable gate")
    return {
        "inner_folds_nonempty": True,
        "feature_names_exclude_outcome_tokens": True,
        "synthetic_stump_recovers_signal": True,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    tests = synthetic_tests()
    feature_rows, feature_audit = build_features()
    targets = tuple(
        row["target_gene"]
        for row in feature_rows
        if row["held"] == FOLDS[0]
    )
    features = feature_matrices(feature_rows, targets)

    # Outcomes are loaded only after the feature table is complete.
    outcomes = load_outcomes(targets)
    pds, gains = outcome_matrices(outcomes, targets)
    candidates = candidate_grid()

    outer_rows: list[dict[str, Any]] = []
    outer_predictions: list[dict[str, Any]] = []
    inner_rows: list[dict[str, Any]] = []
    selected_arms: dict[str, np.ndarray] = {}
    per_fold_gain: dict[str, np.ndarray] = {}
    selected_specs: dict[str, Candidate] = {}
    forced_rows: list[dict[str, Any]] = []

    for test_held in FOLDS:
        train_held = next(held for held in FOLDS if held != test_held)
        selected, model, evaluations, oof_arms = nested_select(
            features[train_held],
            gains[train_held],
            targets,
            candidates,
        )
        arms = predict_candidate(selected, model, features[test_held])
        realized = arm_gain(arms, gains[test_held])
        selected_arms[test_held] = arms
        per_fold_gain[test_held] = realized
        selected_specs[test_held] = selected
        selected_eval = evaluations[0]
        outer_rows.append(
            {
                "train_fold": train_held,
                "test_fold": test_held,
                "selected_candidate": selected.name,
                "selected_kind": selected.kind,
                "fitted_model": summarize_model(selected, model),
                "inner_oof_gain": float(selected_eval["oof_gain"]),
                "inner_oof_arm_counts": arm_counts(oof_arms),
                "inner_oof_switches": int(selected_eval["oof_switches"]),
                "test_fixed_D_pds": float(pds[test_held][:, 3].mean()),
                "test_selected_pds": float(
                    pds[test_held][:, 3].mean() + realized.mean()
                ),
                "test_gain_vs_D": float(realized.mean()),
                **{f"test_count_{arm}": arm_counts(arms)[arm] for arm in ARMS},
            }
        )
        for rank, evaluation in enumerate(evaluations[:15], start=1):
            inner_rows.append(
                {
                    "train_fold": train_held,
                    "test_fold": test_held,
                    "rank": rank,
                    **evaluation,
                }
            )
        for index, target in enumerate(targets):
            chosen = arms[index]
            outer_predictions.append(
                {
                    "train_fold": train_held,
                    "test_fold": test_held,
                    "target_gene": target,
                    "selected_candidate": selected.name,
                    "predicted_arm": chosen,
                    "pds_D": float(pds[test_held][index, 3]),
                    "pds_selected": float(
                        pds[test_held][index, ARMS.index(chosen)]
                    ),
                    "gain_vs_D": float(realized[index]),
                }
            )

        forced, forced_model, forced_evaluations, _ = nested_select(
            features[train_held],
            gains[train_held],
            targets,
            candidates,
            require_nontrivial=True,
        )
        forced_arms = predict_candidate(forced, forced_model, features[test_held])
        forced_gain = arm_gain(forced_arms, gains[test_held])
        forced_rows.append(
            {
                "train_fold": train_held,
                "test_fold": test_held,
                "selected_candidate": forced.name,
                "inner_oof_gain": float(forced_evaluations[0]["oof_gain"]),
                "test_gain_vs_D": float(forced_gain.mean()),
                "test_arm_counts": arm_counts(forced_arms),
            }
        )

    observed = float(np.mean([per_fold_gain[held].mean() for held in FOLDS]))
    bootstrap = bootstrap_cluster(
        per_fold_gain, reps=args.bootstrap_reps, seed=args.seed
    )
    routing_null = matched_routing_null(
        selected_arms,
        gains,
        reps=args.routing_permutations,
        seed=args.seed + 1,
        observed=observed,
    )
    label_null = label_permutation_null(
        features,
        gains,
        targets,
        candidates,
        reps=args.label_permutations,
        seed=args.seed + 2,
        observed=observed,
    )
    enforce_rss("model and null evaluation complete")

    decision_checks = {
        "macro_gain_at_least_0_003": observed >= 0.003,
        "both_outer_folds_nonnegative": all(
            row["test_gain_vs_D"] >= 0 for row in outer_rows
        ),
        "bootstrap_ci_low_positive": bootstrap["ci95"][0] > 0,
        "full_pipeline_permutation_p_at_most_0_05": label_null["one_sided_p"] <= 0.05,
        "matched_routing_p_at_most_0_05": routing_null["one_sided_p"] <= 0.05,
    }
    decision = "SHIP" if all(decision_checks.values()) else "REJECT"

    feature_rows_for_csv = [
        {
            "held": row["held"],
            "target_gene": row["target_gene"],
            **{feature: f"{float(row[feature]):.12g}" for feature in FEATURE_COLUMNS},
        }
        for row in feature_rows
    ]
    write_csv(HERE / "features.csv", feature_rows_for_csv)
    write_csv(HERE / "outer_predictions.csv", outer_predictions)
    write_csv(HERE / "inner_candidate_ranking.csv", inner_rows)
    stability = cross_fold_stability(pds, gains)
    consistency = compact_output_consistency(outcomes)
    associations = feature_associations(features, gains)
    write_csv(HERE / "cross_fold_stability.csv", stability)
    write_csv(HERE / "compact_output_consistency.csv", consistency)
    write_csv(HERE / "feature_associations.csv", associations)

    exp064 = json.loads(EXP064_JSON.read_text())
    result = {
        "experiment": "066-leakage-safe-target-gate",
        "status": "completed",
        "decision": decision,
        "decision_checks": decision_checks,
        "protocol": {
            "primary_estimand": (
                "macro effect-level PDS gain over fixed experiment-064 arm D"
            ),
            "outer_validation": (
                "two directional donor holdouts: train HCT116/test HEK293T and "
                "train HEK293T/test HCT116"
            ),
            "inner_validation": (
                f"{INNER_FOLDS}-fold deterministic target CV within the outer "
                "training donor; all feature scaling and model fitting are "
                "inner-train-only"
            ),
            "models": (
                "fixed arms, univariate decision stumps, and ridge gain models "
                "with at most 14 features and three alternative-arm outputs"
            ),
            "candidate_count": len(candidates),
            "selection": "highest inner out-of-fold realized PDS gain; D is eligible",
            "ship_rule": decision_checks,
            "no_cell_eval2": True,
            "no_eigendecomposition": True,
            "no_cloud": True,
            "no_download": True,
            "one_python_process": True,
            "rss_limit_gib": MAX_RSS_GIB,
        },
        "feature_audit": feature_audit,
        "leakage_audit": {
            "features_built_before_outcomes_loaded": True,
            "outcome_columns_in_features": False,
            "target_identity_encoded": False,
            "held_perturbation_effects_used_in_own_outer_fold": False,
            "held_control_target_expression_used": True,
            "held_control_is_available_at_challenge_inference": True,
            "source_effects_are_public_training_inputs": True,
            "all_300_panel_coordinates_excluded_from_effect_geometry": True,
            "experiment_058_pds_used_as_feature": False,
            "experiment_062_pds_used_as_feature": False,
            "experiment_064_pds_used_only_as_outer_or_inner_labels": True,
            "compact_prior_outputs_used_only_for_diagnostic_consistency": True,
            "outer_test_outcomes_used_for_tuning": False,
            "permutation_preserves_joint_A_B_C_gain_structure": True,
        },
        "outer_results": outer_rows,
        "macro": {
            "fixed_D_pds": float(np.mean([pds[held][:, 3].mean() for held in FOLDS])),
            "selected_pds": float(
                np.mean(
                    [
                        pds[held][:, 3].mean() + per_fold_gain[held].mean()
                        for held in FOLDS
                    ]
                )
            ),
            "gain_vs_D": observed,
            "bootstrap": bootstrap,
            "matched_routing_null": routing_null,
            "full_pipeline_label_permutation_null": label_null,
        },
        "forced_nontrivial_gate": {
            "folds": forced_rows,
            "macro_gain_vs_D": float(
                np.mean([row["test_gain_vs_D"] for row in forced_rows])
            ),
        },
        "oracle_headroom": oracle_and_transfer(pds, gains),
        "cross_fold_stability": stability,
        "compact_output_consistency": consistency,
        "tests": tests,
        "resource_audit": {
            "peak_self_rss_gib": peak_rss_gib(),
            "limit_gib": MAX_RSS_GIB,
            "within_limit": peak_rss_gib() < MAX_RSS_GIB,
            "large_npy_strategy": (
                "tiny direct row/column streams; no retained multi-GB mmap pages"
            ),
        },
        "input_hashes": {
            "run_py": sha256_file(Path(__file__)),
            "exp064_per_target_pds_csv": sha256_file(OUTCOME_PATH),
            "exp064_results_json": sha256_file(EXP064_JSON),
            "exp058_per_target_pds_csv": sha256_file(EXP058_CSV),
            "exp058_results_json": sha256_file(EXP058_JSON),
            "exp062_development_per_target_pds_csv": sha256_file(EXP062_DEV_CSV),
            "panel_csv": sha256_file(PANEL_CSV),
            "genes_csv": sha256_file(GENES_CSV),
            "exp064_bank_delta_hashes_inherited": exp064["input_hashes"]["banks"],
        },
    }
    (HERE / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "decision": decision,
                "macro": result["macro"],
                "outer_results": outer_rows,
                "forced_nontrivial_gate": result["forced_nontrivial_gate"],
                "peak_rss_gib": result["resource_audit"]["peak_self_rss_gib"],
            },
            indent=2,
        )
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument(
        "--routing-permutations", type=int, default=ROUTING_PERMUTATIONS
    )
    parser.add_argument(
        "--label-permutations", type=int, default=LABEL_PERMUTATIONS
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
