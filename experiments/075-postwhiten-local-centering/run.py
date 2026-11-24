"""075 — reciprocal local centering of post-whitened query vectors.

The transform is gallery-independent and query-only.  For each row of a
post-whitened query matrix it finds reciprocal cosine-nearest query neighbors,
subtracts a fraction of their unit-vector centroid, renormalizes, and restores
the original row norm:

    u_i = q_i / ||q_i||
    m_i = mean({u_j: i and j are reciprocal k-nearest neighbors})
    q'_i = ||q_i|| * (u_i - beta*m_i) / ||u_i - beta*m_i||

No hidden reference effect, target label, or challenge score enters the
transform.  Hyperparameters are selected only on the 600 perturbations in the
experiment-044 cache, which are disjoint from the 2026 VCC panel.  Evaluation
uses the frozen post-whitening query cache produced by experiment 073.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")

import numpy as np


HERE = Path(__file__).resolve().parent
MAIN_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
AUX_DONORS = ("rpe1", "nadig_jurkat")
EXPECTED_UPSTREAM_EXPERIMENT = "073-postwhiten-csls-distill"
EXPECTED_GAMMA = 0.65
QUERY_CACHE_SCHEMA = 1
RESULT_SCHEMA = 1
TIE_WINDOW = 5e-4


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


def sha256_json(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


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
    array = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    return np.divide(array, np.maximum(norms, 1e-12))


def pds_per_target_from_scores(scores: np.ndarray) -> np.ndarray:
    if scores.ndim != 2 or scores.shape[0] != scores.shape[1]:
        raise ValueError(f"PDS requires a square score matrix, got {scores.shape}")
    if scores.shape[0] < 2:
        raise ValueError("PDS requires at least two gallery members")
    own = np.diag(scores)[:, None]
    ranks = (scores > own).sum(axis=1)
    ranks = ranks + 0.5 * ((scores == own).sum(axis=1) - 1)
    return 1.0 - ranks / (scores.shape[1] - 1)


def pds_from_scores(scores: np.ndarray) -> float:
    return float(np.mean(pds_per_target_from_scores(scores)))


@dataclass(frozen=True, order=True)
class Candidate:
    k: int
    beta: float

    @property
    def key(self) -> str:
        return f"reciprocal-local-center:k={self.k}:beta={self.beta:g}"


@dataclass(frozen=True)
class QueryGeometry:
    unit: np.ndarray
    row_norms: np.ndarray
    neighbor_order: np.ndarray


@dataclass(frozen=True)
class LocalCentroid:
    values: np.ndarray
    mutual_degree: np.ndarray
    fallback_rows: int


@dataclass(frozen=True)
class FeatureAxes:
    score_mask: np.ndarray
    correction_mask: np.ndarray
    correction_in_score: np.ndarray
    score_genes: tuple[str, ...]
    correction_genes: tuple[str, ...]


def build_candidates(ks: Sequence[int], betas: Sequence[float]) -> list[Candidate]:
    candidates = sorted({Candidate(int(k), float(beta)) for k in ks for beta in betas})
    if not candidates:
        raise ValueError("candidate grid is empty")
    if any(candidate.k < 1 for candidate in candidates):
        raise ValueError("all k values must be positive")
    if any(not (0.0 < candidate.beta <= 1.0) for candidate in candidates):
        raise ValueError("all beta values must be in (0, 1]")
    return candidates


def prepare_query_geometry(query: np.ndarray, max_k: int) -> QueryGeometry:
    """Normalize rows and compute deterministic cosine-neighbor order."""
    query64 = np.asarray(query, dtype=np.float64)
    if query64.ndim != 2 or query64.shape[0] < 3:
        raise ValueError("query must be a 2-D matrix with at least three rows")
    if max_k < 1 or max_k >= query64.shape[0]:
        raise ValueError(f"max_k must be in [1, {query64.shape[0] - 1}]")
    if not np.isfinite(query64).all():
        raise ValueError("query contains non-finite values")
    row_norms = np.linalg.norm(query64, axis=1)
    unit = np.divide(query64, np.maximum(row_norms[:, None], 1e-12))
    similarities = unit @ unit.T
    np.fill_diagonal(similarities, -np.inf)
    # Stable full sorting is cheap for 600 rows and makes ties deterministic.
    neighbor_order = np.argsort(-similarities, axis=1, kind="stable")[:, :max_k]
    return QueryGeometry(
        unit=unit,
        row_norms=row_norms,
        neighbor_order=neighbor_order,
    )


def reciprocal_local_centroid(
    geometry: QueryGeometry,
    k: int,
    *,
    minimum_mutual: int = 2,
) -> LocalCentroid:
    """Return reciprocal-kNN centroids, with a deterministic kNN fallback.

    Reciprocal neighbors reduce the influence of query hubs.  If a row has
    fewer than ``minimum_mutual`` reciprocal neighbors, its ordinary top-k
    set is used so that the transform remains defined for isolated rows.
    """
    if k < 1 or k > geometry.neighbor_order.shape[1]:
        raise ValueError("k is outside the prepared neighbor range")
    n_rows = geometry.unit.shape[0]
    top = geometry.neighbor_order[:, :k]
    directed = np.zeros((n_rows, n_rows), dtype=bool)
    directed[np.arange(n_rows)[:, None], top] = True
    reciprocal = directed & directed.T
    centroid = np.empty_like(geometry.unit)
    degree = reciprocal.sum(axis=1).astype(np.int64)
    fallback_rows = 0
    for row in range(n_rows):
        neighbors = np.flatnonzero(reciprocal[row])
        if len(neighbors) < minimum_mutual:
            neighbors = top[row]
            fallback_rows += 1
        if row in neighbors:
            raise AssertionError("self-neighbor entered local centroid")
        centroid[row] = geometry.unit[neighbors].mean(axis=0)
    return LocalCentroid(
        values=centroid,
        mutual_degree=degree,
        fallback_rows=fallback_rows,
    )


def centered_unit_rows(
    geometry: QueryGeometry,
    centroid: LocalCentroid,
    beta: float,
) -> np.ndarray:
    adjusted = geometry.unit - float(beta) * centroid.values
    norms = np.linalg.norm(adjusted, axis=1, keepdims=True)
    if np.any(norms <= 1e-12):
        raise ValueError("local centering collapsed a query row")
    return adjusted / norms


def apply_reciprocal_local_centering(
    query: np.ndarray,
    *,
    k: int,
    beta: float,
    minimum_mutual: int = 2,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Apply the deployable transform and preserve every input row norm."""
    geometry = prepare_query_geometry(query, k)
    centroid = reciprocal_local_centroid(geometry, k, minimum_mutual=minimum_mutual)
    transformed_unit = centered_unit_rows(geometry, centroid, beta)
    transformed = transformed_unit * geometry.row_norms[:, None]
    zero_rows = geometry.row_norms <= 1e-12
    transformed[zero_rows] = 0.0
    diagnostics: dict[str, float | int] = {
        "rows": int(query.shape[0]),
        "genes": int(query.shape[1]),
        "k": int(k),
        "beta": float(beta),
        "minimum_mutual": int(minimum_mutual),
        "fallback_rows": int(centroid.fallback_rows),
        "mean_mutual_degree": float(np.mean(centroid.mutual_degree)),
        "median_mutual_degree": float(np.median(centroid.mutual_degree)),
        "median_centroid_norm": float(
            np.median(np.linalg.norm(centroid.values, axis=1))
        ),
        "max_row_norm_error": float(
            np.max(np.abs(np.linalg.norm(transformed, axis=1) - geometry.row_norms))
        ),
    }
    return transformed, diagnostics


def candidate_score_matrices(
    query: np.ndarray,
    truth: np.ndarray,
    candidates: Sequence[Candidate],
    *,
    minimum_mutual: int,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[int, dict[str, float | int]]]:
    """Score every candidate without materializing one transformed query per beta."""
    if query.shape != truth.shape:
        raise ValueError(f"query/truth shape mismatch: {query.shape} != {truth.shape}")
    max_k = max(candidate.k for candidate in candidates)
    geometry = prepare_query_geometry(query, max_k)
    truth_unit = row_normalize(truth)
    baseline = geometry.unit @ truth_unit.T
    by_key: dict[str, np.ndarray] = {}
    diagnostics: dict[int, dict[str, float | int]] = {}
    for k in sorted({candidate.k for candidate in candidates}):
        centroid = reciprocal_local_centroid(geometry, k, minimum_mutual=minimum_mutual)
        centroid_truth = centroid.values @ truth_unit.T
        dot = np.sum(geometry.unit * centroid.values, axis=1)
        centroid_norm_sq = np.sum(centroid.values * centroid.values, axis=1)
        diagnostics[k] = {
            "fallback_rows": int(centroid.fallback_rows),
            "mean_mutual_degree": float(np.mean(centroid.mutual_degree)),
            "median_mutual_degree": float(np.median(centroid.mutual_degree)),
            "median_centroid_norm": float(np.median(np.sqrt(centroid_norm_sq))),
        }
        for candidate in candidates:
            if candidate.k != k:
                continue
            denominator = np.sqrt(
                np.maximum(
                    1.0
                    - 2.0 * candidate.beta * dot
                    + candidate.beta**2 * centroid_norm_sq,
                    1e-24,
                )
            )
            by_key[candidate.key] = (
                baseline - candidate.beta * centroid_truth
            ) / denominator[:, None]
            by_key[candidate.key][geometry.row_norms <= 1e-12] = 0.0
    return baseline, by_key, diagnostics


def score_candidates(
    query: np.ndarray,
    truth: np.ndarray,
    candidates: Sequence[Candidate],
    *,
    minimum_mutual: int,
) -> dict[str, object]:
    baseline_scores, candidate_scores, diagnostics = candidate_score_matrices(
        query,
        truth,
        candidates,
        minimum_mutual=minimum_mutual,
    )
    baseline = pds_from_scores(baseline_scores)
    return {
        "baseline": baseline,
        "scores": {
            key: pds_from_scores(scores) for key, scores in candidate_scores.items()
        },
        "deltas": {
            key: pds_from_scores(scores) - baseline
            for key, scores in candidate_scores.items()
        },
        "neighbor_diagnostics": diagnostics,
    }


def select_robust_candidate(
    cells: Sequence[Mapping[str, object]],
    candidates: Sequence[Candidate],
    *,
    minimum_mean_gain: float,
) -> tuple[str, dict[str, object]]:
    """Choose using only candidates nonnegative in every supplied validation cell."""
    if not cells:
        raise ValueError("selection requires at least one validation cell")
    by_key = {candidate.key: candidate for candidate in candidates}
    summaries: dict[str, dict[str, float | list[float]]] = {}
    eligible: list[str] = []
    for candidate in candidates:
        deltas = [float(cell["deltas"][candidate.key]) for cell in cells]
        summary: dict[str, float | list[float]] = {
            "deltas": deltas,
            "minimum_delta": float(min(deltas)),
            "mean_delta": float(np.mean(deltas)),
        }
        summaries[candidate.key] = summary
        if min(deltas) >= 0.0 and np.mean(deltas) >= minimum_mean_gain:
            eligible.append(candidate.key)

    if not eligible:
        return "identity", {
            "selected": "identity",
            "minimum_mean_gain": minimum_mean_gain,
            "eligible_count": 0,
            "candidate_summaries": summaries,
        }

    best_minimum = max(float(summaries[key]["minimum_delta"]) for key in eligible)
    near_minimum = [
        key
        for key in eligible
        if float(summaries[key]["minimum_delta"]) >= best_minimum - TIE_WINDOW
    ]
    best_mean = max(float(summaries[key]["mean_delta"]) for key in near_minimum)
    near_mean = [
        key
        for key in near_minimum
        if float(summaries[key]["mean_delta"]) >= best_mean - TIE_WINDOW
    ]
    # Prefer the smaller intervention, then the broader neighborhood.
    selected = min(
        near_mean,
        key=lambda key: (by_key[key].beta, -by_key[key].k, key),
    )
    return selected, {
        "selected": selected,
        "minimum_mean_gain": minimum_mean_gain,
        "eligible_count": len(eligible),
        "best_minimum_delta": best_minimum,
        "best_mean_delta_within_tie_window": best_mean,
        "candidate_summaries": summaries,
    }


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
    target_count = deltas_by_fold.shape[1]
    for lo in range(0, replicates, chunk):
        hi = min(lo + chunk, replicates)
        indices = rng.integers(0, target_count, size=(hi - lo, target_count))
        draws[lo:hi] = deltas_by_fold[:, indices].mean(axis=(0, 2))
    return {
        "replicates": int(replicates),
        "seed": int(seed),
        "mean": float(np.mean(draws)),
        "lower_95": float(np.quantile(draws, 0.025)),
        "upper_95": float(np.quantile(draws, 0.975)),
        "probability_positive": float(np.mean(draws > 0.0)),
    }


def verify_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 {actual} != expected {expected}")


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
    train_set = set(train_perts)
    vcc_set = set(vcc_panel)
    if len(train_perts) != 600:
        raise ValueError(f"expected 600 training perturbations, got {len(train_perts)}")
    if len(outer_perts) < 240:
        raise ValueError(f"only {len(outer_perts)} outer perturbations")
    if train_set & vcc_set:
        raise ValueError("training perturbations overlap the VCC panel")
    if not set(outer_perts) <= vcc_set:
        raise ValueError("outer perturbations are not a VCC-panel subset")
    if set(outer_perts) & train_set:
        raise ValueError("training and outer perturbations overlap")
    if any(gene in vcc_set for gene in master_genes):
        raise ValueError("a VCC target coordinate entered the master score axis")
    if np.any(
        train_mask
        & np.asarray([gene in train_set for gene in master_genes], dtype=bool)
    ):
        raise ValueError("a training-target coordinate is active in train_mask")
    if not np.any(train_mask):
        raise ValueError("truth-cache train_mask has no active coordinates")
    if not set(master_genes) <= set(official_genes):
        raise ValueError("truth-cache master genes are outside the official axis")
    return master_genes, train_perts, outer_perts, train_mask


def derive_feature_axes(
    master_genes: Sequence[str],
    train_perts: Sequence[str],
    train_mask: np.ndarray,
) -> FeatureAxes:
    train_set = set(train_perts)
    score_mask = np.asarray([gene not in train_set for gene in master_genes])
    correction_mask = np.asarray(train_mask, dtype=bool)
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
        raise ValueError("expected a nonempty strict correction subset")
    return FeatureAxes(
        score_mask=score_mask,
        correction_mask=correction_mask,
        correction_in_score=correction_in_score,
        score_genes=score_genes,
        correction_genes=correction_genes,
    )


def verify_inputs(
    args: argparse.Namespace,
    upstream_manifest: Mapping[str, object],
) -> None:
    challenge = upstream_manifest["challenge"]
    verify_hash(args.genes_csv, challenge["gene_names.csv"], "gene_names.csv")
    verify_hash(args.panel_csv, challenge["pert_counts.csv"], "pert_counts.csv")
    for filename, expected in upstream_manifest["truth_cache"]["files"].items():
        verify_hash(args.truth_cache / filename, expected, f"truth cache/{filename}")


def validate_query_cache(
    query_cache: Path,
    *,
    train_count: int,
    outer_count: int,
    axes: FeatureAxes,
    train_order_sha256: str,
    outer_order_sha256: str,
    upstream_manifest: Mapping[str, object],
) -> dict[str, object]:
    metadata_path = query_cache / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    expected = {
        "schema_version": QUERY_CACHE_SCHEMA,
        "experiment": EXPECTED_UPSTREAM_EXPERIMENT,
        "gamma": EXPECTED_GAMMA,
        "training_targets": train_count,
        "outer_targets": outer_count,
        "train_outer_overlap": 0,
        "outer_score_genes": len(axes.score_genes),
        "correction_genes": len(axes.correction_genes),
        "vcc_target_coordinate_overlap": 0,
        "training_target_coordinate_overlap": 0,
        "whitening_or_full_covariance_eigendecompositions": 0,
        "train_order_sha256": train_order_sha256,
        "outer_order_sha256": outer_order_sha256,
        "score_gene_order_sha256": sha256_lines(axes.score_genes),
        "correction_gene_order_sha256": sha256_lines(axes.correction_genes),
        "correction_in_score_sha256": sha256_array(axes.correction_in_score),
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"query-cache metadata {key}={metadata.get(key)!r}, expected {value!r}"
            )
    for held in MAIN_DONORS:
        fold = metadata["folds"][held]
        operator = upstream_manifest["operators"][held]
        expected_sources = tuple(name for name in MAIN_DONORS if name != held)
        if tuple(fold.get("source_donors", ())) != expected_sources:
            raise ValueError(f"{held}: query sources are not strict outer LODO")
        if tuple(operator.get("source_donors", ())) != expected_sources:
            raise ValueError(f"{held}: frozen operator sources are not strict LODO")
        if fold.get("operator_npz_sha256") != operator["npz_sha256"]:
            raise ValueError(f"{held}: frozen operator NPZ hash mismatch")
        if fold.get("operator_json_sha256") != operator["json_sha256"]:
            raise ValueError(f"{held}: frozen operator metadata hash mismatch")
        for panel, rows in (("train", train_count), ("outer", outer_count)):
            path = query_cache / held / f"{panel}_query.npy"
            values = np.load(path, mmap_mode="r")
            if values.shape != (rows, len(axes.score_genes)):
                raise ValueError(f"{path}: unexpected shape {values.shape}")
            if sha256_file(path) != fold["panels"][panel]["query_sha256"]:
                raise ValueError(f"{path}: cache hash mismatch")
    return metadata


def load_truth(
    truth_cache: Path,
    donor: str,
    panel: str,
    mask: np.ndarray,
) -> np.ndarray:
    path = truth_cache / f"{donor}_raw_{panel}.npy"
    values = np.load(path, mmap_mode="r")
    if values.shape[1] != len(mask):
        raise ValueError(f"{path}: feature width mismatch")
    return np.asarray(values[:, mask], dtype=np.float64)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_experiment(
    args: argparse.Namespace,
    upstream_manifest: Mapping[str, object],
) -> dict[str, object]:
    official = read_csv_column(args.genes_csv, "gene_name")
    vcc_panel = read_csv_column(args.panel_csv, "target_gene")
    master_genes, train_perts, outer_perts, train_mask = validate_truth_cache(
        args.truth_cache, official, vcc_panel
    )
    axes = derive_feature_axes(master_genes, train_perts, train_mask)
    query_metadata = validate_query_cache(
        args.query_cache,
        train_count=len(train_perts),
        outer_count=len(outer_perts),
        axes=axes,
        train_order_sha256=sha256_lines(train_perts),
        outer_order_sha256=sha256_lines(outer_perts),
        upstream_manifest=upstream_manifest,
    )
    candidates = build_candidates(args.ks, args.betas)
    candidate_by_key = {candidate.key: candidate for candidate in candidates}

    # Phase 1: freeze all selectors using only the disjoint non-VCC panel.
    auxiliary_truth = {
        donor: load_truth(args.truth_cache, donor, "train", axes.correction_mask)
        for donor in AUX_DONORS
    }
    inner_cells_by_role: dict[str, list[dict[str, object]]] = {}
    inner_detail: dict[str, object] = {}
    per_role_selected: dict[str, str] = {}
    all_inner_cells: list[dict[str, object]] = []
    for outer_role in MAIN_DONORS:
        train_query_full = np.asarray(
            np.load(
                args.query_cache / outer_role / "train_query.npy",
                mmap_mode="r",
            ),
            dtype=np.float64,
        )
        train_query = train_query_full[:, axes.correction_in_score]
        role_cells: list[dict[str, object]] = []
        teacher_detail: dict[str, object] = {}
        for teacher in AUX_DONORS:
            scored = score_candidates(
                train_query,
                auxiliary_truth[teacher],
                candidates,
                minimum_mutual=args.minimum_mutual,
            )
            cell = {
                "outer_role": outer_role,
                "teacher": teacher,
                "baseline": scored["baseline"],
                "scores": scored["scores"],
                "deltas": scored["deltas"],
            }
            role_cells.append(cell)
            all_inner_cells.append(cell)
            teacher_detail[teacher] = scored
        selected, selection = select_robust_candidate(
            role_cells,
            candidates,
            minimum_mean_gain=args.minimum_inner_gain,
        )
        inner_cells_by_role[outer_role] = role_cells
        per_role_selected[outer_role] = selected
        inner_detail[outer_role] = {
            "source_donors": [donor for donor in MAIN_DONORS if donor != outer_role],
            "teachers": teacher_detail,
            "selection": selection,
        }
        del train_query, train_query_full
        enforce_memory_limit(args.max_rss_gib, f"{outer_role} inner scoring")

    role_nested: dict[str, object] = {}
    role_nested_validation_deltas: list[float] = []
    for held_role in MAIN_DONORS:
        fit_cells = [
            cell
            for role, cells in inner_cells_by_role.items()
            if role != held_role
            for cell in cells
        ]
        selected, selection = select_robust_candidate(
            fit_cells,
            candidates,
            minimum_mean_gain=args.minimum_inner_gain,
        )
        held_cells = inner_cells_by_role[held_role]
        held_deltas = [
            0.0 if selected == "identity" else float(cell["deltas"][selected])
            for cell in held_cells
        ]
        role_nested_validation_deltas.extend(held_deltas)
        role_nested[held_role] = {
            "selected_without_role": selected,
            "fit_roles": [role for role in MAIN_DONORS if role != held_role],
            "validation_role": held_role,
            "validation_teacher_deltas": {
                teacher: delta
                for teacher, delta in zip(AUX_DONORS, held_deltas, strict=True)
            },
            "validation_mean_delta": float(np.mean(held_deltas)),
            "selection": selection,
        }

    deployment_selected, deployment_selection = select_robust_candidate(
        all_inner_cells,
        candidates,
        minimum_mean_gain=args.minimum_inner_gain,
    )
    selection_freeze = {
        "deployment_selected": deployment_selected,
        "per_role_selected": per_role_selected,
        "role_nested_selected": {
            role: value["selected_without_role"] for role, value in role_nested.items()
        },
        "training_order_sha256": sha256_lines(train_perts),
        "candidate_grid": [candidate.key for candidate in candidates],
    }
    selection_freeze_sha256 = sha256_json(selection_freeze)

    # Phase 2: only now open the held main-donor VCC-panel truths.
    outer_folds: dict[str, object] = {}
    per_target_rows: list[dict[str, object]] = []
    deployment_target_deltas: list[np.ndarray] = []
    for outer_role in MAIN_DONORS:
        outer_query = np.asarray(
            np.load(
                args.query_cache / outer_role / "outer_query.npy",
                mmap_mode="r",
            ),
            dtype=np.float64,
        )
        outer_truth = load_truth(args.truth_cache, outer_role, "outer", axes.score_mask)
        baseline_scores = row_normalize(outer_query) @ row_normalize(outer_truth).T
        baseline_per_target = pds_per_target_from_scores(baseline_scores)
        baseline = float(np.mean(baseline_per_target))

        keys = {
            key
            for key in (deployment_selected, per_role_selected[outer_role])
            if key != "identity"
        }
        needed = [candidate_by_key[key] for key in sorted(keys)]
        score_by_key: dict[str, np.ndarray] = {}
        diagnostics: dict[int, dict[str, float | int]] = {}
        if needed:
            _, score_by_key, diagnostics = candidate_score_matrices(
                outer_query,
                outer_truth,
                needed,
                minimum_mutual=args.minimum_mutual,
            )

        def per_target(key: str) -> np.ndarray:
            return (
                baseline_per_target.copy()
                if key == "identity"
                else pds_per_target_from_scores(score_by_key[key])
            )

        deployment_per_target = per_target(deployment_selected)
        nested_selected = per_role_selected[outer_role]
        nested_per_target = per_target(nested_selected)
        deployment_delta = deployment_per_target - baseline_per_target
        deployment_target_deltas.append(deployment_delta)
        for target, base, deployed, nested in zip(
            outer_perts,
            baseline_per_target,
            deployment_per_target,
            nested_per_target,
            strict=True,
        ):
            per_target_rows.append(
                {
                    "outer_donor": outer_role,
                    "target_gene": target,
                    "deployment_config": deployment_selected,
                    "nested_config": nested_selected,
                    "baseline_pds": float(base),
                    "deployment_pds": float(deployed),
                    "deployment_delta": float(deployed - base),
                    "nested_pds": float(nested),
                    "nested_delta": float(nested - base),
                }
            )
        outer_folds[outer_role] = {
            "source_donors": [donor for donor in MAIN_DONORS if donor != outer_role],
            "baseline": baseline,
            "deployment_config": deployment_selected,
            "deployment_pds": float(np.mean(deployment_per_target)),
            "deployment_delta": float(np.mean(deployment_delta)),
            "nested_config": nested_selected,
            "nested_pds": float(np.mean(nested_per_target)),
            "nested_delta": float(np.mean(nested_per_target - baseline_per_target)),
            "neighbor_diagnostics": diagnostics,
        }
        del outer_query, outer_truth, baseline_scores
        enforce_memory_limit(args.max_rss_gib, f"{outer_role} outer scoring")

    deployment_deltas = [
        float(outer_folds[role]["deployment_delta"]) for role in MAIN_DONORS
    ]
    deployment_macro_delta = float(np.mean(deployment_deltas))
    deployment_macro_baseline = float(
        np.mean([outer_folds[role]["baseline"] for role in MAIN_DONORS])
    )
    deployment_macro_pds = deployment_macro_baseline + deployment_macro_delta
    nested_macro_delta = float(
        np.mean([outer_folds[role]["nested_delta"] for role in MAIN_DONORS])
    )
    bootstrap = clustered_target_bootstrap(
        np.vstack(deployment_target_deltas),
        replicates=args.bootstrap_replicates,
        seed=args.bootstrap_seed,
    )
    role_nested_macro = float(np.mean(role_nested_validation_deltas))
    promotion_gate = {
        "predeclared": True,
        "deployment_candidate_nonidentity": deployment_selected != "identity",
        "source_role_nested_all_teacher_cells_nonnegative": all(
            delta >= 0.0 for delta in role_nested_validation_deltas
        ),
        "source_role_nested_macro_at_least_0_001": role_nested_macro >= 0.001,
        "outer_all_three_deltas_nonnegative": all(
            delta >= 0.0 for delta in deployment_deltas
        ),
        "outer_macro_delta_at_least_0_003": deployment_macro_delta >= 0.003,
        "outer_cluster_bootstrap_lower_95_positive": bootstrap["lower_95"] > 0.0,
    }
    promotion_gate["passed"] = all(
        value for key, value in promotion_gate.items() if key != "predeclared"
    )

    result: dict[str, object] = {
        "schema_version": RESULT_SCHEMA,
        "experiment": "075-postwhiten-local-centering",
        "status": "completed",
        "method": {
            "name": "reciprocal-kNN local query centering",
            "formula": (
                "q_i' = ||q_i|| * normalize(normalize(q_i) - beta * "
                "mean(reciprocal_kNN(normalized queries)))"
            ),
            "query_only": True,
            "hidden_gallery_used": False,
            "target_labels_used_by_transform": False,
            "row_norm_preserved": True,
            "minimum_mutual_neighbors_before_knn_fallback": args.minimum_mutual,
        },
        "protocol": {
            "outer_donors": list(MAIN_DONORS),
            "auxiliary_truth_donors": list(AUX_DONORS),
            "training_perturbations": len(train_perts),
            "outer_vcc_perturbations": len(outer_perts),
            "train_outer_overlap": len(set(train_perts) & set(outer_perts)),
            "inner_correction_genes": len(axes.correction_genes),
            "outer_score_genes": len(axes.score_genes),
            "upstream_query_cache": EXPECTED_UPSTREAM_EXPERIMENT,
            "upstream_gamma": EXPECTED_GAMMA,
            "outer_truth_load_order": (
                "outer truth opened only after per-role and deployment "
                "hyperparameters were frozen"
            ),
        },
        "leakage_audit": {
            "training_and_vcc_perturbations_disjoint": True,
            "vcc_target_coordinates_in_query_cache": False,
            "training_target_coordinates_in_query_cache": False,
            "outer_donor_used_in_its_query_or_whitener": False,
            "outer_truth_used_in_selection": False,
            "hidden_reference_gallery_used_by_transform": False,
            "challenge_context_controls_used": False,
        },
        "candidate_grid": {
            "ks": sorted({candidate.k for candidate in candidates}),
            "betas": sorted({candidate.beta for candidate in candidates}),
            "count": len(candidates),
        },
        "inner": {
            "per_role": inner_detail,
            "role_leave_one_out": role_nested,
            "role_nested_validation_macro_delta": role_nested_macro,
            "deployment_selection": deployment_selection,
            "selection_freeze": selection_freeze,
            "selection_freeze_sha256": selection_freeze_sha256,
        },
        "outer": {
            "folds": outer_folds,
            "deployment": {
                "selected": deployment_selected,
                "baseline_macro_pds": deployment_macro_baseline,
                "macro_pds": deployment_macro_pds,
                "macro_delta": deployment_macro_delta,
                "fold_deltas": deployment_deltas,
                "bootstrap": bootstrap,
            },
            "per_role_nested_macro_delta": nested_macro_delta,
        },
        "promotion_gate": promotion_gate,
        "query_cache_metadata_sha256": sha256_file(args.query_cache / "metadata.json"),
        "upstream_manifest_sha256": sha256_file(args.upstream_manifest),
        "query_cache_summary": {
            "experiment": query_metadata["experiment"],
            "gamma": query_metadata["gamma"],
            "outer_score_genes": query_metadata["outer_score_genes"],
            "correction_genes": query_metadata["correction_genes"],
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(args.per_target_output, per_target_rows)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-cache", type=Path, required=True)
    parser.add_argument("--query-cache", type=Path, required=True)
    parser.add_argument("--upstream-manifest", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument(
        "--per-target-output",
        type=Path,
        default=HERE / "per_target_pds.csv",
    )
    parser.add_argument("--ks", nargs="+", type=int, default=[8, 16, 32, 64])
    parser.add_argument(
        "--betas",
        nargs="+",
        type=float,
        default=[0.10, 0.25, 0.50, 0.75],
    )
    parser.add_argument("--minimum-mutual", type=int, default=2)
    parser.add_argument("--minimum-inner-gain", type=float, default=0.001)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--bootstrap-seed", type=int, default=7502026)
    parser.add_argument("--max-rss-gib", type=float, default=2.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    if args.minimum_mutual < 1:
        raise SystemExit("--minimum-mutual must be positive")
    if args.minimum_inner_gain < 0:
        raise SystemExit("--minimum-inner-gain must be nonnegative")
    if args.bootstrap_replicates < 100:
        raise SystemExit("--bootstrap-replicates must be at least 100")
    if args.max_rss_gib <= 0:
        raise SystemExit("--max-rss-gib must be positive")
    upstream_manifest = json.loads(args.upstream_manifest.read_text())
    if int(upstream_manifest.get("schema_version", -1)) != 1:
        raise SystemExit("unsupported upstream artifact manifest")
    verify_inputs(args, upstream_manifest)
    result = run_experiment(args, upstream_manifest)
    print(json.dumps(result["outer"]["deployment"], indent=2), flush=True)
    print(json.dumps(result["promotion_gate"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
