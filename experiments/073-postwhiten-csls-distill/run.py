"""073 — nested CSLS-density distillation after frozen union whitening.

The expensive upstream geometry is fixed to experiment 062's three
gamma=0.65 source-only union operators. This script never refits whitening or
a full covariance eigensystem; the distiller retains experiment 044's small
randomized low-rank SVD. It has two phases:

1. Stream the three main donor banks, construct each outer fold's source-only
   query matrices, apply the matching frozen operator in row chunks, and cache
   only the leakage-safe score coordinates.
2. Distill exact-CSLS candidate-density penalties with strict RPE1/Jurkat
   inner donor LODO, then evaluate the selected query-only transform once on
   the disjoint VCC panel of the held main donor.

All 300 VCC targets and all 600 non-VCC selection targets are excluded from
the effective input/output feature map.  Outer truth is never transformed.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

# Bound BLAS before importing NumPy.  The cloud runner can override these.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

MAIN_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
AUX_DONORS = ("rpe1", "nadig_jurkat")
WHITEN_GAMMA = 0.65
EFFECT_SCALE = 6.817
CACHE_SCHEMA = 1
RESULT_SCHEMA = 1


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS and Linux."""
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached {limit_gib:.3f} GiB"
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
            raise ValueError(f"{path}: required column {column!r} is missing")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column!r}")
    return values


def row_normalize(values: np.ndarray) -> np.ndarray:
    output = np.asarray(values, dtype=np.float64)
    return output / np.maximum(np.linalg.norm(output, axis=1, keepdims=True), 1e-12)


def pds_per_target_from_scores(scores: np.ndarray) -> np.ndarray:
    if scores.ndim != 2 or scores.shape[0] != scores.shape[1]:
        raise ValueError(f"PDS requires a square score matrix, got {scores.shape}")
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    if n < 2:
        raise ValueError("PDS requires at least two gallery members")
    ranks = (scores > own).sum(axis=1)
    ranks = ranks + 0.5 * ((scores == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


def pds_from_scores(scores: np.ndarray) -> float:
    return float(np.mean(pds_per_target_from_scores(scores)))


def candidate_density(scores: np.ndarray, k: int) -> np.ndarray:
    """Mean top-k query similarities for every candidate/gallery column."""
    if k < 1:
        raise ValueError("k must be positive")
    k = min(k, scores.shape[0])
    return np.partition(scores, -k, axis=0)[-k:, :].mean(axis=0)


def exact_csls_scores(q: np.ndarray, truth: np.ndarray, k: int) -> np.ndarray:
    scores = q @ truth.T
    density = candidate_density(scores, k)
    # 2*S-r_q-r_t has the same within-query ranking as S-0.5*r_t.
    return scores - 0.5 * density[None, :]


def stable_seed(names: Iterable[str], extra: int = 0) -> int:
    payload = "|".join(sorted(names)) + f"|{extra}"
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:4], "little")


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    pert_pos: dict[str, int]
    delta_path: Path
    shape: tuple[int, int]

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        path = path.expanduser().resolve()
        genes = read_lines(path / "genes.txt", unique=False)
        perts = read_lines(path / "perts.txt")
        delta_path = path / "delta.npy"
        delta = np.load(delta_path, mmap_mode="r")
        if not isinstance(delta, np.memmap):
            raise TypeError(f"{name}: delta.npy did not memory-map")
        expected = (len(perts), len(genes))
        if delta.shape != expected:
            raise ValueError(f"{name}: delta shape {delta.shape} != {expected}")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            shape=tuple(map(int, delta.shape)),
        )


@dataclass(frozen=True)
class AxisMap:
    all_source_columns: np.ndarray
    unique_source_columns: np.ndarray
    unique_official_positions: np.ndarray
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
        raise ValueError(f"{bank.name}: no genes map to the official axis")
    unique_positions = np.asarray(sorted(last_source), dtype=np.int64)
    unique_sources = np.asarray(
        [last_source[position] for position in unique_positions], dtype=np.int64
    )
    return AxisMap(
        all_source_columns=np.asarray(all_source, dtype=np.int64),
        unique_source_columns=unique_sources,
        unique_official_positions=unique_positions,
        duplicate_symbols=tuple(sorted(duplicates)),
    )


@dataclass(frozen=True)
class BankStats:
    mapping: AxisMap
    unique_mean: np.ndarray
    raw_median_norm: float


def prepare_bank_stats(
    bank: Bank,
    official_genes: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> BankStats:
    """Reproduce experiment 062's full-catalogue centering and norm scale."""
    mapping = map_bank_to_axis(bank, official_genes)
    delta = np.load(bank.delta_path, mmap_mode="r")
    mean_all = np.zeros(len(mapping.all_source_columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        mean_all += np.asarray(
            delta[lo:hi, mapping.all_source_columns], dtype=np.float64
        ).sum(axis=0, dtype=np.float64)
        enforce_memory_limit(max_rss_gib, f"{bank.name} centering")
    mean_all /= bank.shape[0]

    norms = np.empty(bank.shape[0], dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        centered = (
            np.asarray(delta[lo:hi, mapping.all_source_columns], dtype=np.float64)
            - mean_all
        )
        norms[lo:hi] = np.linalg.norm(centered, axis=1)
        enforce_memory_limit(max_rss_gib, f"{bank.name} norm calibration")
    raw_median_norm = float(np.median(norms)) or 1.0

    all_position = {
        int(source): index
        for index, source in enumerate(mapping.all_source_columns.tolist())
    }
    unique_mean = np.asarray(
        [
            mean_all[all_position[int(source)]]
            for source in mapping.unique_source_columns
        ],
        dtype=np.float64,
    )
    del delta, mean_all, norms
    return BankStats(
        mapping=mapping,
        unique_mean=unique_mean,
        raw_median_norm=raw_median_norm,
    )


def extract_scaled_rows(
    bank: Bank,
    stats: BankStats,
    perturbations: Sequence[str],
    operator_columns: np.ndarray,
    official_to_operator: Mapping[int, int],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    """Extract 062-style centered/scaled effects directly on an operator axis."""
    missing = [pert for pert in perturbations if pert not in bank.pert_pos]
    if missing:
        raise ValueError(
            f"{bank.name}: missing {len(missing)} requested perturbations; "
            f"first={missing[:5]}"
        )
    source_columns: list[int] = []
    destination_columns: list[int] = []
    means: list[float] = []
    for source, official, mean in zip(
        stats.mapping.unique_source_columns,
        stats.mapping.unique_official_positions,
        stats.unique_mean,
        strict=True,
    ):
        destination = official_to_operator.get(int(official))
        if destination is None:
            continue
        source_columns.append(int(source))
        destination_columns.append(destination)
        means.append(float(mean))
    if not source_columns:
        raise ValueError(f"{bank.name}: no genes map to the operator axis")

    source_array = np.asarray(source_columns, dtype=np.int64)
    destination_array = np.asarray(destination_columns, dtype=np.int64)
    mean_array = np.asarray(means, dtype=np.float64)
    rows = np.asarray([bank.pert_pos[pert] for pert in perturbations], dtype=np.int64)
    output = np.zeros((len(perturbations), len(operator_columns)), dtype=np.float32)
    delta = np.load(bank.delta_path, mmap_mode="r")
    multiplier = EFFECT_SCALE / stats.raw_median_norm
    for lo in range(0, len(rows), chunk_rows):
        hi = min(lo + chunk_rows, len(rows))
        values = np.asarray(delta[np.ix_(rows[lo:hi], source_array)], dtype=np.float64)
        output[lo:hi, destination_array] = ((values - mean_array) * multiplier).astype(
            np.float32
        )
        enforce_memory_limit(max_rss_gib, f"{bank.name} selected-row extraction")
    return output


def build_source_query(
    banks: Mapping[str, Bank],
    stats: Mapping[str, BankStats],
    source_names: Sequence[str],
    perturbations: Sequence[str],
    operator_columns: np.ndarray,
    blocked_operator_positions: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    official_to_operator = {
        int(position): index for index, position in enumerate(operator_columns.tolist())
    }
    query = np.zeros((len(perturbations), len(operator_columns)), dtype=np.float32)
    for source_name in source_names:
        query += extract_scaled_rows(
            banks[source_name],
            stats[source_name],
            perturbations,
            operator_columns,
            official_to_operator,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
    query /= float(len(source_names))
    query[:, blocked_operator_positions] = 0.0
    return query


@dataclass(frozen=True)
class Operator:
    evals: np.ndarray
    evecs: np.ndarray
    metadata: dict[str, object]


def validate_operator_metadata(
    held: str,
    source_names: Sequence[str],
    metadata: Mapping[str, object],
    operator_columns: np.ndarray,
    official_genes: Sequence[str],
    vcc_panel: set[str],
) -> None:
    if float(metadata.get("gamma", -1.0)) != WHITEN_GAMMA:
        raise ValueError(f"{held}: frozen operator gamma is not {WHITEN_GAMMA}")
    if tuple(metadata.get("donors", [])) != tuple(source_names):
        raise ValueError(
            f"{held}: operator donors {metadata.get('donors')} != {source_names}"
        )
    if int(metadata.get("n_cols", -1)) != len(operator_columns):
        raise ValueError(f"{held}: operator width mismatch")
    if str(metadata.get("wcols_sha256")) != sha256_array(operator_columns):
        raise ValueError(f"{held}: reconstructed operator-axis hash mismatch")
    excluded = set(map(str, metadata.get("excluded_targets", [])))
    if excluded != vcc_panel:
        raise ValueError(f"{held}: operator did not exclude the exact VCC panel")
    if any(official_genes[index] in vcc_panel for index in operator_columns):
        raise AssertionError(f"{held}: VCC coordinate entered operator axis")


def load_frozen_operator(
    held: str,
    operator_dir: Path,
    *,
    expected_npz_sha256: str,
    expected_json_sha256: str,
) -> Operator:
    npz_path = operator_dir / f"held_{held}.npz"
    metadata_path = operator_dir / f"held_{held}.json"
    if sha256_file(npz_path) != expected_npz_sha256:
        raise ValueError(f"{held}: frozen operator NPZ hash mismatch")
    if sha256_file(metadata_path) != expected_json_sha256:
        raise ValueError(f"{held}: frozen operator metadata hash mismatch")
    metadata = json.loads(metadata_path.read_text())
    with np.load(npz_path, allow_pickle=False) as archive:
        evals = np.asarray(archive["evals"], dtype=np.float64)
        # Keep the 1.2 GiB eigenvector matrix float32.
        evecs = np.asarray(archive["evecs"], dtype=np.float32)
    if evecs.shape != (len(evals), len(evals)):
        raise ValueError(f"{held}: malformed frozen eigensystem {evecs.shape}")
    return Operator(evals=evals, evecs=evecs, metadata=metadata)


def apply_frozen_operator(
    query: np.ndarray,
    operator: Operator,
    blocked_operator_positions: np.ndarray,
    selected_operator_positions: np.ndarray,
    output_path: Path,
    *,
    transform_chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, float]:
    """Apply the principal active-coordinate map without refitting.

    Excluded coordinates are zero on input and output.  Selecting only active
    outputs makes this mathematically the W[A,A] principal submatrix action,
    where W is experiment 062's frozen dense operator.
    """
    if query.shape[1] != len(operator.evals):
        raise ValueError("query/operator width mismatch")
    if np.any(query[:, blocked_operator_positions] != 0):
        raise AssertionError("blocked coordinate is nonzero before whitening")
    weights = np.asarray(operator.evals ** (-WHITEN_GAMMA / 2.0), dtype=np.float32)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output = np.lib.format.open_memmap(
        output_path,
        mode="w+",
        dtype=np.float32,
        shape=(query.shape[0], len(selected_operator_positions)),
    )
    before_norms = np.linalg.norm(query, axis=1)
    after_norms = np.empty(query.shape[0], dtype=np.float64)
    for lo in range(0, query.shape[0], transform_chunk_rows):
        hi = min(lo + transform_chunk_rows, query.shape[0])
        block = np.asarray(query[lo:hi], dtype=np.float32)
        latent = block @ operator.evecs
        latent *= weights
        transformed = latent @ operator.evecs.T
        transformed[:, blocked_operator_positions] = 0.0
        after_norms[lo:hi] = np.linalg.norm(transformed, axis=1)
        output[lo:hi] = transformed[:, selected_operator_positions]
        output.flush()
        del block, latent, transformed
        enforce_memory_limit(max_rss_gib, "frozen operator application")
    before = float(np.median(before_norms)) or 1.0
    after = float(np.median(after_norms)) or 1.0
    preserve_scale = before / after
    for lo in range(0, query.shape[0], transform_chunk_rows):
        hi = min(lo + transform_chunk_rows, query.shape[0])
        output[lo:hi] *= preserve_scale
    output.flush()
    del output
    return {
        "input_median_norm_active_map": before,
        "output_median_norm_active_map_before_preserve": after,
        "global_median_norm_preserve_scale": preserve_scale,
    }


def load_manifest(path: Path) -> dict[str, object]:
    manifest = json.loads(path.read_text())
    if int(manifest.get("schema_version", -1)) != 1:
        raise ValueError(f"{path}: unsupported manifest schema")
    return manifest


def verify_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 {actual} != expected {expected}")


def verify_input_artifacts(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
    *,
    include_staging_inputs: bool,
) -> None:
    challenge = manifest["challenge"]
    verify_hash(args.genes_csv, challenge["gene_names.csv"], "gene_names.csv")
    verify_hash(args.panel_csv, challenge["pert_counts.csv"], "pert_counts.csv")

    cache_manifest = manifest["truth_cache"]
    for filename, expected in cache_manifest["files"].items():
        verify_hash(args.truth_cache / filename, expected, f"truth cache/{filename}")

    if include_staging_inputs:
        if args.donor_root is None or args.operator_dir is None:
            raise ValueError("staging requires --donor-root and --operator-dir")
        donor_manifest = manifest["donors"]
        for name in MAIN_DONORS:
            path = args.donor_root / name
            for filename, expected in donor_manifest[name].items():
                verify_hash(path / filename, expected, f"{name}/{filename}")

        operator_manifest = manifest["operators"]
        for held in MAIN_DONORS:
            entry = operator_manifest[held]
            verify_hash(
                args.operator_dir / f"held_{held}.npz",
                entry["npz_sha256"],
                f"operator/{held}.npz",
            )
            verify_hash(
                args.operator_dir / f"held_{held}.json",
                entry["json_sha256"],
                f"operator/{held}.json",
            )


def validate_truth_cache(
    truth_cache: Path,
    official_genes: Sequence[str],
    vcc_panel: Sequence[str],
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], np.ndarray]:
    master_genes = read_lines(truth_cache / "master_genes.txt")
    train_perts = read_lines(truth_cache / "train_perts.txt")
    outer_perts = read_lines(truth_cache / "outer_perts.txt")
    train_mask = np.asarray(
        np.load(truth_cache / "train_mask.npy", allow_pickle=False), dtype=bool
    )
    if len(master_genes) != len(train_mask):
        raise ValueError("truth-cache train_mask length mismatch")
    vcc_set = set(vcc_panel)
    train_set = set(train_perts)
    if train_set & vcc_set:
        raise ValueError("training perturbations overlap the VCC panel")
    if not set(outer_perts) <= vcc_set:
        raise ValueError("outer perturbations are not a VCC-panel subset")
    if set(outer_perts) & train_set:
        raise ValueError("training and outer perturbations overlap")
    if len(train_perts) != 600:
        raise ValueError(f"expected 600 training perturbations, got {len(train_perts)}")
    if len(outer_perts) < 240:
        raise ValueError(f"only {len(outer_perts)} outer perturbations")
    target_positions = np.asarray(
        [gene in train_set for gene in master_genes], dtype=bool
    )
    if np.any(train_mask & target_positions):
        raise ValueError("a training-target coordinate is active in train_mask")
    if not np.any(train_mask):
        raise ValueError("truth-cache train_mask has no active coordinates")
    # Experiment 044 also removed genes absent from either auxiliary teacher.
    # That stricter coverage mask is authenticated by the artifact hash; the
    # invariant needed here is that it never reactivates a training target.
    if any(gene in vcc_set for gene in master_genes):
        raise ValueError("VCC target coordinate entered master score axis")
    official_set = set(official_genes)
    if not set(master_genes) <= official_set:
        raise ValueError("truth-cache master genes are outside the official axis")
    return master_genes, train_perts, outer_perts, train_mask


@dataclass(frozen=True)
class FeatureAxes:
    score_mask: np.ndarray
    correction_mask: np.ndarray
    correction_in_score: np.ndarray
    score_genes: tuple[str, ...]
    correction_genes: tuple[str, ...]


def derive_feature_axes(
    master_genes: Sequence[str],
    train_perts: Sequence[str],
    saved_correction_mask: np.ndarray,
) -> FeatureAxes:
    """Separate full outer retrieval geometry from the teacher-common basis."""
    if len(master_genes) != len(saved_correction_mask):
        raise ValueError("feature-axis mask length mismatch")
    train_set = set(train_perts)
    score_mask = np.asarray(
        [gene not in train_set for gene in master_genes], dtype=bool
    )
    correction_mask = np.asarray(saved_correction_mask, dtype=bool)
    if np.any(correction_mask & ~score_mask):
        raise ValueError("correction axis contains a training-target coordinate")
    score_indices = np.flatnonzero(score_mask)
    correction_indices = np.flatnonzero(correction_mask)
    score_local = {int(master): local for local, master in enumerate(score_indices)}
    correction_in_score = np.asarray(
        [score_local[int(master)] for master in correction_indices], dtype=np.int64
    )
    score_genes = tuple(
        gene for gene, keep in zip(master_genes, score_mask, strict=True) if keep
    )
    correction_genes = tuple(
        gene for gene, keep in zip(master_genes, correction_mask, strict=True) if keep
    )
    if not correction_genes or len(correction_genes) >= len(score_genes):
        raise ValueError(
            "expected a nonempty strict correction subset of the score axis"
        )
    return FeatureAxes(
        score_mask=score_mask,
        correction_mask=correction_mask,
        correction_in_score=correction_in_score,
        score_genes=score_genes,
        correction_genes=correction_genes,
    )


def stage_queries(
    args: argparse.Namespace, manifest: Mapping[str, object]
) -> dict[str, object]:
    """Build query caches without refitting whitening/full covariance geometry."""
    if args.donor_root is None or args.operator_dir is None:
        raise ValueError("staging requires --donor-root and --operator-dir")
    official = read_csv_column(args.genes_csv, "gene_name")
    panel = read_csv_column(args.panel_csv, "target_gene")
    panel_set = set(panel)
    master_genes, train_perts, outer_perts, train_mask = validate_truth_cache(
        args.truth_cache, official, panel
    )
    axes = derive_feature_axes(master_genes, train_perts, train_mask)
    blocked_targets = set(train_perts) | panel_set
    if blocked_targets & set(axes.score_genes):
        raise AssertionError("a panel target entered the outer score axis")

    banks = {name: Bank.load(name, args.donor_root / name) for name in MAIN_DONORS}
    stats: dict[str, BankStats] = {}
    for name in MAIN_DONORS:
        stats[name] = prepare_bank_stats(
            banks[name],
            official,
            chunk_rows=args.source_chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        print(
            f"prepared {name}: rows={banks[name].shape[0]} "
            f"mapped={len(stats[name].mapping.unique_official_positions)} "
            f"median_norm={stats[name].raw_median_norm:.6f} "
            f"peak_rss={peak_rss_gib():.2f} GiB",
            flush=True,
        )

    args.query_cache.mkdir(parents=True, exist_ok=True)
    fold_metadata: dict[str, object] = {}
    official_pos = {gene: index for index, gene in enumerate(official)}
    operator_manifest = manifest["operators"]
    for held in MAIN_DONORS:
        source_names = tuple(name for name in MAIN_DONORS if name != held)
        donor_union = set().union(*(set(banks[name].genes) for name in source_names))
        operator_columns = np.asarray(
            [
                index
                for index, gene in enumerate(official)
                if gene in donor_union and gene not in panel_set
            ],
            dtype=np.int64,
        )
        metadata_path = args.operator_dir / f"held_{held}.json"
        operator_metadata = json.loads(metadata_path.read_text())
        validate_operator_metadata(
            held,
            source_names,
            operator_metadata,
            operator_columns,
            official,
            panel_set,
        )
        official_to_operator = {
            int(position): index
            for index, position in enumerate(operator_columns.tolist())
        }
        score_official = np.asarray(
            [official_pos[gene] for gene in axes.score_genes], dtype=np.int64
        )
        missing_score = [
            gene
            for gene, position in zip(axes.score_genes, score_official, strict=True)
            if int(position) not in official_to_operator
        ]
        if missing_score:
            raise ValueError(
                f"{held}: {len(missing_score)} score genes absent from operator; "
                f"first={missing_score[:5]}"
            )
        selected_operator_positions = np.asarray(
            [official_to_operator[int(position)] for position in score_official],
            dtype=np.int64,
        )
        blocked_operator_positions = np.asarray(
            [
                index
                for index, official_index in enumerate(operator_columns)
                if official[official_index] in blocked_targets
            ],
            dtype=np.int64,
        )
        blocked_symbols = {
            official[int(operator_columns[index])]
            for index in blocked_operator_positions
        }
        expected_blocked = set(train_perts) & {
            official[int(index)] for index in operator_columns
        }
        if blocked_symbols != expected_blocked:
            raise AssertionError(f"{held}: blocked operator-coordinate mismatch")
        if np.intersect1d(selected_operator_positions, blocked_operator_positions).size:
            raise AssertionError(f"{held}: selected and blocked coordinates overlap")

        entry = operator_manifest[held]
        operator = load_frozen_operator(
            held,
            args.operator_dir,
            expected_npz_sha256=entry["npz_sha256"],
            expected_json_sha256=entry["json_sha256"],
        )
        if len(operator.evals) != len(operator_columns):
            raise ValueError(f"{held}: operator array width mismatch")

        fold_dir = args.query_cache / held
        fold_dir.mkdir(parents=True, exist_ok=True)
        panel_metadata: dict[str, object] = {}
        for label, perturbations in (
            ("train", train_perts),
            ("outer", outer_perts),
        ):
            query = build_source_query(
                banks,
                stats,
                source_names,
                perturbations,
                operator_columns,
                blocked_operator_positions,
                chunk_rows=args.source_chunk_rows,
                max_rss_gib=args.max_rss_gib,
            )
            output_path = fold_dir / f"{label}_query.npy"
            transform = apply_frozen_operator(
                query,
                operator,
                blocked_operator_positions,
                selected_operator_positions,
                output_path,
                transform_chunk_rows=args.transform_chunk_rows,
                max_rss_gib=args.max_rss_gib,
            )
            panel_metadata[label] = {
                "rows": len(perturbations),
                "columns": len(axes.score_genes),
                "query_sha256": sha256_file(output_path),
                "transform": transform,
            }
            del query
            gc.collect()
            print(
                f"staged held={held} panel={label} "
                f"shape=({len(perturbations)},{len(axes.score_genes)}) "
                f"peak_rss={peak_rss_gib():.2f} GiB",
                flush=True,
            )

        fold_metadata[held] = {
            "source_donors": list(source_names),
            "operator_columns": len(operator_columns),
            "operator_columns_sha256": sha256_array(operator_columns),
            "blocked_training_target_coordinates": len(blocked_operator_positions),
            "blocked_training_target_symbols_sha256": sha256_lines(
                sorted(blocked_symbols)
            ),
            "selected_score_coordinates": len(selected_operator_positions),
            "selected_score_symbols_sha256": sha256_lines(axes.score_genes),
            "correction_coordinates": len(axes.correction_genes),
            "correction_symbols_sha256": sha256_lines(axes.correction_genes),
            "operator_npz_sha256": entry["npz_sha256"],
            "operator_json_sha256": entry["json_sha256"],
            "panels": panel_metadata,
        }
        del operator
        gc.collect()
        enforce_memory_limit(args.max_rss_gib, f"{held} stage completion")

    metadata: dict[str, object] = {
        "schema_version": CACHE_SCHEMA,
        "experiment": "073-postwhiten-csls-distill",
        "status": "query_cache_complete",
        "gamma": WHITEN_GAMMA,
        "whitening_or_full_covariance_eigendecompositions": 0,
        "operator_action": (
            "reuse frozen experiment-062 eigensystems; zero both-panel target "
            "coordinates on input and retain the full leakage-safe score outputs"
        ),
        "official_genes": len(official),
        "vcc_targets": len(panel),
        "training_targets": len(train_perts),
        "outer_targets": len(outer_perts),
        "train_outer_overlap": len(set(train_perts) & set(outer_perts)),
        "master_score_genes": len(master_genes),
        "outer_score_genes": len(axes.score_genes),
        "correction_genes": len(axes.correction_genes),
        "score_only_genes": len(axes.score_genes) - len(axes.correction_genes),
        "vcc_target_coordinate_overlap": len(panel_set & set(axes.score_genes)),
        "training_target_coordinate_overlap": len(
            set(train_perts) & set(axes.score_genes)
        ),
        "official_order_sha256": sha256_lines(official),
        "vcc_panel_order_sha256": sha256_lines(panel),
        "train_order_sha256": sha256_lines(train_perts),
        "outer_order_sha256": sha256_lines(outer_perts),
        "score_gene_order_sha256": sha256_lines(axes.score_genes),
        "correction_gene_order_sha256": sha256_lines(axes.correction_genes),
        "correction_in_score_sha256": sha256_array(axes.correction_in_score),
        "folds": fold_metadata,
        "peak_rss_gib": peak_rss_gib(),
    }
    metadata_path = args.query_cache / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


@dataclass(frozen=True)
class Candidate:
    method: str
    k: int
    rank: int
    ridge: float
    eta: float

    @property
    def key(self) -> str:
        return (
            f"{self.method}:k={self.k}:r={self.rank}:"
            f"ridge={self.ridge:g}:eta={self.eta:g}"
        )


@dataclass(frozen=True)
class Fitted:
    basis: np.ndarray
    bias: np.ndarray | None = None
    bands: np.ndarray | None = None


def build_candidates(args: argparse.Namespace) -> list[Candidate]:
    return [
        Candidate(method, k, rank, ridge, eta)
        for method in args.methods
        for k in args.ks
        for rank in args.ranks
        for ridge in args.ridges
        for eta in args.etas
    ]


def randomized_right_basis(
    matrix: np.ndarray,
    rank: int,
    seed: int,
    *,
    oversample: int = 8,
    power_iters: int = 1,
) -> np.ndarray:
    matrix32 = np.asarray(matrix, dtype=np.float32)
    width = min(rank + oversample, min(matrix32.shape))
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((matrix32.shape[1], width), dtype=np.float32)
    projected = matrix32 @ omega
    for _ in range(power_iters):
        q, _ = np.linalg.qr(projected, mode="reduced")
        projected = matrix32 @ (matrix32.T @ q)
    q, _ = np.linalg.qr(projected, mode="reduced")
    small = q.T @ matrix32
    _, _, vt = np.linalg.svd(small, full_matrices=False)
    return vt[:rank].T.astype(np.float64)


def standardized_ridge(
    gram: np.ndarray,
    rhs: np.ndarray,
    observations: int,
    ridge: float,
) -> np.ndarray:
    gram = gram / max(observations, 1)
    rhs = rhs / max(observations, 1)
    scale = np.sqrt(np.maximum(np.diag(gram), 1e-12))
    standardized_gram = gram / scale[:, None] / scale[None, :]
    standardized_rhs = rhs / scale
    coefficient = np.linalg.solve(
        standardized_gram + ridge * np.eye(len(scale)),
        standardized_rhs,
    )
    return coefficient / scale


def fit_bundle(
    query: np.ndarray,
    truths: Mapping[str, np.ndarray],
    candidates: Sequence[Candidate],
    *,
    seed_names: Iterable[str],
) -> dict[str, Fitted]:
    if not truths:
        raise ValueError("at least one teacher donor is required")
    max_rank = max(candidate.rank for candidate in candidates)
    normalized_truths = {name: row_normalize(value) for name, value in truths.items()}
    basis = randomized_right_basis(
        np.vstack(list(normalized_truths.values())),
        max_rank,
        stable_seed(seed_names, extra=max_rank),
    )
    normalized_query = row_normalize(query)
    fitted: dict[str, Fitted] = {}
    for k in sorted({candidate.k for candidate in candidates}):
        projections: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        labels: list[np.ndarray] = []
        for truth in normalized_truths.values():
            scores = normalized_query @ truth.T
            label = -0.5 * candidate_density(scores, k)
            label -= label.mean()
            projections.append((normalized_query @ basis, truth @ basis, label))
            labels.append(label)

        for rank in sorted({candidate.rank for candidate in candidates}):
            local_basis = basis[:, :rank]
            bias_features = np.vstack(
                [truth_projection[:, :rank] for _, truth_projection, _ in projections]
            )
            bias_labels = np.concatenate(labels)
            bias_gram = bias_features.T @ bias_features
            bias_rhs = bias_features.T @ bias_labels

            bands_gram = np.zeros((rank, rank), dtype=np.float64)
            bands_rhs = np.zeros(rank, dtype=np.float64)
            combined_gram = np.zeros((2 * rank, 2 * rank), dtype=np.float64)
            combined_rhs = np.zeros(2 * rank, dtype=np.float64)
            pair_observations = 0
            for query_projection0, truth_projection0, label in projections:
                query_projection = query_projection0[:, :rank]
                truth_projection = truth_projection0[:, :rank]
                query_gram = query_projection.T @ query_projection
                truth_gram = truth_projection.T @ truth_projection
                query_sum = query_projection.sum(axis=0)
                truth_label = truth_projection.T @ label
                bands_gram += query_gram * truth_gram
                bands_rhs += query_sum * truth_label
                pair_observations += len(query_projection) * len(truth_projection)

                combined_gram[:rank, :rank] += len(query_projection) * truth_gram
                cross = truth_gram * query_sum[None, :]
                combined_gram[:rank, rank:] += cross
                combined_gram[rank:, :rank] += cross.T
                combined_gram[rank:, rank:] += query_gram * truth_gram
                combined_rhs[:rank] += len(query_projection) * truth_label
                combined_rhs[rank:] += query_sum * truth_label

            for ridge in sorted({candidate.ridge for candidate in candidates}):
                bias = standardized_ridge(
                    bias_gram, bias_rhs, len(bias_features), ridge
                )
                bands = standardized_ridge(
                    bands_gram, bands_rhs, pair_observations, ridge
                )
                combined = standardized_ridge(
                    combined_gram, combined_rhs, pair_observations, ridge
                )
                for method in sorted({candidate.method for candidate in candidates}):
                    for eta in sorted({candidate.eta for candidate in candidates}):
                        candidate = Candidate(method, k, rank, ridge, eta)
                        if candidate not in candidates:
                            continue
                        if method == "bias":
                            fitted[candidate.key] = Fitted(basis=local_basis, bias=bias)
                        elif method == "bands":
                            fitted[candidate.key] = Fitted(
                                basis=local_basis, bands=bands
                            )
                        elif method == "combo":
                            fitted[candidate.key] = Fitted(
                                basis=local_basis,
                                bias=combined[:rank],
                                bands=combined[rank:],
                            )
                        else:
                            raise ValueError(f"unknown method {method!r}")
    return fitted


def fitted_score_matrix(
    query_score: np.ndarray,
    truth_score: np.ndarray,
    fitted: Fitted,
    eta: float,
    correction_in_score: np.ndarray | None = None,
) -> np.ndarray:
    normalized_query = row_normalize(query_score)
    normalized_truth = row_normalize(truth_score)
    scores = normalized_query @ normalized_truth.T
    if correction_in_score is None:
        correction_query = normalized_query
        correction_truth = normalized_truth
    else:
        correction_query = normalized_query[:, correction_in_score]
        correction_truth = normalized_truth[:, correction_in_score]
    if correction_query.shape[1] != fitted.basis.shape[0]:
        raise ValueError(
            "correction feature width does not match the fitted distillation basis"
        )
    query_projection = correction_query @ fitted.basis
    truth_projection = correction_truth @ fitted.basis
    if fitted.bias is not None:
        scores += eta * (truth_projection @ fitted.bias)[None, :]
    if fitted.bands is not None:
        scores += eta * ((query_projection * fitted.bands) @ truth_projection.T)
    return scores


def evaluate_bundle(
    query_score: np.ndarray,
    truth_score: np.ndarray,
    fitted: Mapping[str, Fitted],
    candidates: Sequence[Candidate],
    correction_in_score: np.ndarray | None = None,
) -> tuple[float, dict[str, float], dict[int, float]]:
    normalized_query = row_normalize(query_score)
    normalized_truth = row_normalize(truth_score)
    baseline_scores = normalized_query @ normalized_truth.T
    baseline = pds_from_scores(baseline_scores)
    by_key = {candidate.key: candidate for candidate in candidates}
    scores = {
        key: pds_from_scores(
            fitted_score_matrix(
                query_score,
                truth_score,
                model,
                by_key[key].eta,
                correction_in_score,
            )
        )
        for key, model in fitted.items()
    }
    oracle = {
        k: pds_from_scores(exact_csls_scores(normalized_query, normalized_truth, k))
        for k in sorted({candidate.k for candidate in candidates})
    }
    return baseline, scores, oracle


def choose_candidate(
    mean_scores: Mapping[str, float],
    baseline: float,
    candidates: Sequence[Candidate],
    min_gain: float,
) -> str:
    by_key = {candidate.key: candidate for candidate in candidates}

    def complexity(key: str) -> tuple[float | int | str, ...]:
        candidate = by_key[key]
        return (
            candidate.rank,
            -candidate.ridge,
            abs(candidate.eta - 1.0),
            candidate.k,
            candidate.method,
        )

    best_score = max(mean_scores.values())
    near = [key for key, score in mean_scores.items() if score >= best_score - 5e-4]
    selected = min(near, key=complexity)
    return selected if mean_scores[selected] >= baseline + min_gain else "baseline"


def load_active_truth(
    truth_cache: Path,
    donor: str,
    panel: str,
    active_mask: np.ndarray,
) -> np.ndarray:
    path = truth_cache / f"{donor}_raw_{panel}.npy"
    values = np.load(path, mmap_mode="r")
    if values.shape[1] != len(active_mask):
        raise ValueError(f"{path}: feature width mismatch")
    return np.asarray(values[:, active_mask], dtype=np.float64)


def validate_query_cache(
    query_cache: Path,
    train_count: int,
    outer_count: int,
    score_count: int,
    *,
    train_order_sha256: str,
    outer_order_sha256: str,
    score_gene_order_sha256: str,
    correction_gene_order_sha256: str,
    correction_in_score_sha256: str,
    correction_count: int,
    operator_manifest: Mapping[str, object],
) -> dict[str, object]:
    metadata = json.loads((query_cache / "metadata.json").read_text())
    if int(metadata.get("schema_version", -1)) != CACHE_SCHEMA:
        raise ValueError("unsupported query-cache schema")
    expected = {
        "training_targets": train_count,
        "outer_targets": outer_count,
        "outer_score_genes": score_count,
        "correction_genes": correction_count,
        "score_only_genes": score_count - correction_count,
        "train_outer_overlap": 0,
        "vcc_target_coordinate_overlap": 0,
        "training_target_coordinate_overlap": 0,
        "whitening_or_full_covariance_eigendecompositions": 0,
        "train_order_sha256": train_order_sha256,
        "outer_order_sha256": outer_order_sha256,
        "score_gene_order_sha256": score_gene_order_sha256,
        "correction_gene_order_sha256": correction_gene_order_sha256,
        "correction_in_score_sha256": correction_in_score_sha256,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"query-cache metadata {key}={metadata.get(key)!r}, expected {value!r}"
            )
    for held in MAIN_DONORS:
        fold = metadata["folds"][held]
        operator_entry = operator_manifest[held]
        if fold.get("operator_npz_sha256") != operator_entry["npz_sha256"]:
            raise ValueError(f"{held}: staged operator NPZ hash is not frozen")
        if fold.get("operator_json_sha256") != operator_entry["json_sha256"]:
            raise ValueError(f"{held}: staged operator metadata hash is not frozen")
        if tuple(fold.get("source_donors", [])) != tuple(
            operator_entry["source_donors"]
        ):
            raise ValueError(f"{held}: staged source donors do not match manifest")
        for panel, rows in (("train", train_count), ("outer", outer_count)):
            path = query_cache / held / f"{panel}_query.npy"
            values = np.load(path, mmap_mode="r")
            if values.shape != (rows, score_count):
                raise ValueError(f"{path}: unexpected shape {values.shape}")
            if sha256_file(path) != fold["panels"][panel]["query_sha256"]:
                raise ValueError(f"{path}: cache hash mismatch")
    return metadata


def clustered_target_bootstrap(
    deltas_by_fold: np.ndarray,
    *,
    replicates: int,
    seed: int,
    chunk: int = 1000,
) -> dict[str, float | int]:
    if deltas_by_fold.ndim != 2:
        raise ValueError("bootstrap input must be fold x target")
    rng = np.random.default_rng(seed)
    draws = np.empty(replicates, dtype=np.float64)
    n_targets = deltas_by_fold.shape[1]
    for lo in range(0, replicates, chunk):
        hi = min(lo + chunk, replicates)
        indices = rng.integers(0, n_targets, size=(hi - lo, n_targets))
        sampled = deltas_by_fold[:, indices]
        draws[lo:hi] = sampled.mean(axis=(0, 2))
    return {
        "replicates": replicates,
        "seed": seed,
        "mean": float(np.mean(draws)),
        "lower_95": float(np.quantile(draws, 0.025)),
        "upper_95": float(np.quantile(draws, 0.975)),
        "probability_positive": float(np.mean(draws > 0)),
    }


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_distillation(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
) -> dict[str, object]:
    official = read_csv_column(args.genes_csv, "gene_name")
    panel = read_csv_column(args.panel_csv, "target_gene")
    master_genes, train_perts, outer_perts, train_mask = validate_truth_cache(
        args.truth_cache, official, panel
    )
    axes = derive_feature_axes(master_genes, train_perts, train_mask)
    query_metadata = validate_query_cache(
        args.query_cache,
        len(train_perts),
        len(outer_perts),
        len(axes.score_genes),
        train_order_sha256=sha256_lines(train_perts),
        outer_order_sha256=sha256_lines(outer_perts),
        score_gene_order_sha256=sha256_lines(axes.score_genes),
        correction_gene_order_sha256=sha256_lines(axes.correction_genes),
        correction_in_score_sha256=sha256_array(axes.correction_in_score),
        correction_count=len(axes.correction_genes),
        operator_manifest=manifest["operators"],
    )
    candidates = build_candidates(args)
    candidate_by_key = {candidate.key: candidate for candidate in candidates}

    aux_truths = {
        name: load_active_truth(args.truth_cache, name, "train", axes.correction_mask)
        for name in AUX_DONORS
    }
    outer_truths = {
        name: load_active_truth(args.truth_cache, name, "outer", axes.score_mask)
        for name in MAIN_DONORS
    }
    folds: dict[str, object] = {}
    per_target_rows: list[dict[str, object]] = []
    selected_target_deltas: list[np.ndarray] = []

    for outer in MAIN_DONORS:
        train_query = np.asarray(
            np.load(args.query_cache / outer / "train_query.npy", mmap_mode="r"),
            dtype=np.float64,
        )
        outer_query = np.asarray(
            np.load(args.query_cache / outer / "outer_query.npy", mmap_mode="r"),
            dtype=np.float64,
        )
        train_query_correction = train_query[:, axes.correction_in_score]
        inner_scores = {candidate.key: [] for candidate in candidates}
        inner_baselines: list[float] = []
        inner_detail: dict[str, object] = {}
        for valid in AUX_DONORS:
            fit_names = tuple(name for name in AUX_DONORS if name != valid)
            fitted = fit_bundle(
                train_query_correction,
                {name: aux_truths[name] for name in fit_names},
                candidates,
                seed_names=("inner", outer, *fit_names),
            )
            baseline, scores, oracle = evaluate_bundle(
                train_query_correction, aux_truths[valid], fitted, candidates
            )
            inner_baselines.append(baseline)
            for key, score in scores.items():
                inner_scores[key].append(score)
            inner_detail[valid] = {
                "fit_teacher_donors": list(fit_names),
                "validation_teacher_donor": valid,
                "baseline": baseline,
                "oracle_by_k": oracle,
                "best_candidate_score": max(scores.values()),
            }
            del fitted
            gc.collect()

        mean_inner = {
            key: float(np.mean(values)) for key, values in inner_scores.items()
        }
        inner_baseline = float(np.mean(inner_baselines))
        selected = choose_candidate(
            mean_inner, inner_baseline, candidates, args.min_inner_gain
        )

        final_bundle = fit_bundle(
            train_query_correction,
            aux_truths,
            candidates,
            seed_names=("final", outer, *AUX_DONORS),
        )
        outer_truth = outer_truths[outer]
        outer_baseline, outer_scores, outer_oracle = evaluate_bundle(
            outer_query,
            outer_truth,
            final_bundle,
            candidates,
            axes.correction_in_score,
        )
        selected_score = (
            outer_baseline if selected == "baseline" else outer_scores[selected]
        )
        baseline_matrix = row_normalize(outer_query) @ row_normalize(outer_truth).T
        baseline_per_target = pds_per_target_from_scores(baseline_matrix)
        if selected == "baseline":
            selected_per_target = baseline_per_target.copy()
        else:
            selected_per_target = pds_per_target_from_scores(
                fitted_score_matrix(
                    outer_query,
                    outer_truth,
                    final_bundle[selected],
                    candidate_by_key[selected].eta,
                    axes.correction_in_score,
                )
            )
        delta_per_target = selected_per_target - baseline_per_target
        selected_target_deltas.append(delta_per_target)
        for target, baseline_value, selected_value in zip(
            outer_perts,
            baseline_per_target,
            selected_per_target,
            strict=True,
        ):
            per_target_rows.append(
                {
                    "outer_donor": outer,
                    "target_gene": target,
                    "selected_config": selected,
                    "baseline_pds": float(baseline_value),
                    "selected_pds": float(selected_value),
                    "delta": float(selected_value - baseline_value),
                }
            )

        fixed_order = sorted(
            outer_scores,
            key=lambda key: outer_scores[key],
            reverse=True,
        )
        folds[outer] = {
            "source_donors": [name for name in MAIN_DONORS if name != outer],
            "inner_teacher_donors": list(AUX_DONORS),
            "inner": {
                "baseline": inner_baseline,
                "selected": selected,
                "selected_score": (
                    inner_baseline if selected == "baseline" else mean_inner[selected]
                ),
                "selected_delta": (
                    0.0
                    if selected == "baseline"
                    else mean_inner[selected] - inner_baseline
                ),
                "detail": inner_detail,
            },
            "outer": {
                "baseline": outer_baseline,
                "selected": selected_score,
                "delta": selected_score - outer_baseline,
                "oracle_csls_by_k": outer_oracle,
                "oracle_csls_best": max(outer_oracle.values()),
                "oracle_headroom": max(outer_oracle.values()) - outer_baseline,
                "fixed_grid_optimistic_top": [
                    {
                        "config": key,
                        "pds": outer_scores[key],
                        "delta": outer_scores[key] - outer_baseline,
                    }
                    for key in fixed_order[:20]
                ],
            },
        }
        print(
            f"outer={outer} inner_selected={selected} "
            f"baseline={outer_baseline:.6f} selected={selected_score:.6f} "
            f"delta={selected_score - outer_baseline:+.6f} "
            f"oracle_headroom={max(outer_oracle.values()) - outer_baseline:+.6f}",
            flush=True,
        )
        del (
            train_query,
            train_query_correction,
            outer_query,
            final_bundle,
            baseline_matrix,
        )
        gc.collect()
        enforce_memory_limit(args.max_rss_gib, f"{outer} distillation completion")

    macro_baseline = float(
        np.mean([folds[name]["outer"]["baseline"] for name in MAIN_DONORS])
    )
    macro_selected = float(
        np.mean([folds[name]["outer"]["selected"] for name in MAIN_DONORS])
    )
    macro_oracle = float(
        np.mean([folds[name]["outer"]["oracle_csls_best"] for name in MAIN_DONORS])
    )
    bootstrap = clustered_target_bootstrap(
        np.vstack(selected_target_deltas),
        replicates=args.bootstrap_replicates,
        seed=args.bootstrap_seed,
    )
    fold_deltas = [folds[name]["outer"]["delta"] for name in MAIN_DONORS]
    promotion = {
        "minimum_required_macro_delta": 0.003,
        "all_outer_folds_nonnegative": all(delta >= 0 for delta in fold_deltas),
        "macro_delta_at_least_0_003": macro_selected - macro_baseline >= 0.003,
        "bootstrap_lower_bound_positive": bootstrap["lower_95"] > 0,
    }
    promotion["passed"] = all(
        value
        for key, value in promotion.items()
        if key != "minimum_required_macro_delta"
    )

    result: dict[str, object] = {
        "schema_version": RESULT_SCHEMA,
        "experiment": "073-postwhiten-csls-distill",
        "status": "completed",
        "protocol": {
            "outer_donors": list(MAIN_DONORS),
            "inner_teacher_donors": list(AUX_DONORS),
            "training_perturbations": len(train_perts),
            "outer_vcc_perturbations": len(outer_perts),
            "train_outer_overlap": len(set(train_perts) & set(outer_perts)),
            "outer_score_genes": len(axes.score_genes),
            "correction_genes": len(axes.correction_genes),
            "score_only_genes": len(axes.score_genes) - len(axes.correction_genes),
            "gamma": WHITEN_GAMMA,
            "operator_source": "frozen experiment-062 fold-specific eigensystems",
            "operator_refits": 0,
            "whitening_or_full_covariance_eigendecompositions": 0,
            "randomized_low_rank_distillation_svd": True,
            "query": (
                "equal two-main-donor 062-style raw-centered/norm-scaled "
                "consensus, then frozen source-only union whitening"
            ),
            "truth": "raw centered, untransformed",
            "selection": (
                "within each outer fold, fit on one auxiliary teacher and "
                "validate on the other in both directions; select by mean "
                "inner PDS with experiment-044 complexity tie-break and "
                f"minimum gain {args.min_inner_gain:g}"
            ),
            "outer_evaluation": (
                "one selected transform on disjoint VCC-panel truth from the "
                "held main donor; baseline and exact CSLS use the full outer "
                "score axis while the learned correction is embedded only on "
                "the teacher-common correction subset"
            ),
        },
        "leakage_audit": {
            "outer_donor_used_in_query": False,
            "outer_donor_used_in_operator_fit": False,
            "outer_donor_used_in_distillation_fit": False,
            "outer_donor_used_in_hyperparameter_selection": False,
            "inner_validation_teacher_used_in_inner_fit": False,
            "training_and_outer_perturbations_disjoint": True,
            "vcc_target_coordinates_in_effective_map": False,
            "training_target_coordinates_in_effective_map": False,
            "outer_truth_transformed": False,
            "outer_scores_used_to_select_candidate": False,
        },
        "query_cache": query_metadata,
        "folds": folds,
        "nested_macro": {
            "baseline": macro_baseline,
            "selected": macro_selected,
            "delta": macro_selected - macro_baseline,
            "oracle_csls": macro_oracle,
            "oracle_headroom": macro_oracle - macro_baseline,
            "bootstrap": bootstrap,
        },
        "promotion_gate": promotion,
        "candidate_grid": {
            "methods": list(args.methods),
            "ks": list(args.ks),
            "ranks": list(args.ranks),
            "ridges": list(args.ridges),
            "etas": list(args.etas),
            "count": len(candidates),
        },
        "artifact_manifest_sha256": sha256_file(args.artifact_manifest),
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(args.per_target_output, per_target_rows)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-root", type=Path)
    parser.add_argument("--truth-cache", type=Path, required=True)
    parser.add_argument("--operator-dir", type=Path)
    parser.add_argument("--artifact-manifest", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--query-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument(
        "--per-target-output",
        type=Path,
        default=HERE / "per_target_pds.csv",
    )
    parser.add_argument(
        "--phase",
        choices=("all", "stage", "distill"),
        default="all",
    )
    parser.add_argument("--source-chunk-rows", type=int, default=128)
    parser.add_argument("--transform-chunk-rows", type=int, default=16)
    parser.add_argument("--max-rss-gib", type=float, default=48.0)
    parser.add_argument("--methods", nargs="+", default=["bias", "bands", "combo"])
    parser.add_argument("--ks", nargs="+", type=int, default=[5, 10, 20])
    parser.add_argument("--ranks", nargs="+", type=int, default=[8, 16, 32])
    parser.add_argument("--ridges", nargs="+", type=float, default=[0.01, 0.1, 1.0])
    parser.add_argument("--etas", nargs="+", type=float, default=[0.5, 1.0])
    parser.add_argument("--min-inner-gain", type=float, default=0.001)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--bootstrap-seed", type=int, default=7302026)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    if args.source_chunk_rows < 1 or args.transform_chunk_rows < 1:
        raise SystemExit("chunk sizes must be positive")
    if args.max_rss_gib <= 0:
        raise SystemExit("--max-rss-gib must be positive")
    if not set(args.methods) <= {"bias", "bands", "combo"}:
        raise SystemExit("methods must be bias, bands, and/or combo")
    manifest = load_manifest(args.artifact_manifest)
    verify_input_artifacts(
        args,
        manifest,
        include_staging_inputs=args.phase in {"all", "stage"},
    )
    if args.phase in {"all", "stage"}:
        stage_queries(args, manifest)
    if args.phase in {"all", "distill"}:
        result = run_distillation(args, manifest)
        print(json.dumps(result["nested_macro"], indent=2), flush=True)
        print(json.dumps(result["promotion_gate"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
