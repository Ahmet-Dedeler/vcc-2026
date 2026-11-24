#!/usr/bin/env python3
"""Leakage-safe listwise low-rank adapter over frozen arm-D predictions.

The trainable map is a near-identity residual restricted to a source-derived
feature basis:

    z' = normalize(z + beta * B U V^T B^T z)

where B is fitted only from non-panel fitting targets, rank(UV^T) is at most
16, and beta is fixed by the candidate grid.  The objective is a
temperature-scaled listwise softmax: the positive candidate is the same target
in another donor and every other target in that donor is a negative.

Held VCC-panel truths are not opened until target splits, hyperparameters,
weights for the non-panel pre-gate, and the label-permutation null are frozen.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import resource
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
DEVELOPMENT_FOLDS = ("orion_hct116", "orion_hek293t")
KOLF_FOLD = "kolf2"
FINAL_FOLDS = (*DEVELOPMENT_FOLDS, KOLF_FOLD)
PRIMARY_FOLDS = DEVELOPMENT_FOLDS
OUTER_DONORS = DEVELOPMENT_FOLDS
TEACHER_DONORS = ("rpe1", "nadig_jurkat")
IDENTITY = "identity"
RESULT_SCHEMA = 5
CACHE_SCHEMA = 5
ANGLE_LIMIT_DEGREES = 10.0
EXPECTED_FINAL_SCORE_GENES = 17_806


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


def sha256_json(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def stable_seed(*parts: object) -> int:
    text = "|".join(map(str, parts))
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4], "little")


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_gene_lines(path: Path) -> tuple[str, ...]:
    """Read a gene axis where duplicate symbols are valid source columns."""
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values:
        raise ValueError(f"{path}: empty gene axis")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column!r}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate values in {column!r}")
    return values


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def row_normalize(values: np.ndarray) -> np.ndarray:
    values64 = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(values64, axis=1, keepdims=True)
    return np.divide(values64, np.maximum(norms, 1e-12))


def pds_per_query_from_scores(
    scores: np.ndarray,
    positive_columns: np.ndarray,
) -> np.ndarray:
    if scores.ndim != 2:
        raise ValueError(f"PDS requires a score matrix, got {scores.shape}")
    if scores.shape[1] < 2:
        raise ValueError("PDS requires at least two gallery members")
    positives = np.asarray(positive_columns, dtype=np.int64)
    if positives.shape != (scores.shape[0],):
        raise ValueError("one positive-gallery index is required per query")
    if np.any(positives < 0) or np.any(positives >= scores.shape[1]):
        raise ValueError("positive-gallery index is out of bounds")
    own = scores[np.arange(scores.shape[0]), positives][:, None]
    ranks = (scores > own).sum(axis=1)
    ranks = ranks + 0.5 * ((scores == own).sum(axis=1) - 1)
    return 1.0 - ranks / (scores.shape[1] - 1)


def pds_per_target_from_scores(scores: np.ndarray) -> np.ndarray:
    if scores.ndim != 2 or scores.shape[0] != scores.shape[1]:
        raise ValueError(f"PDS requires a square score matrix, got {scores.shape}")
    return pds_per_query_from_scores(scores, np.arange(scores.shape[0], dtype=np.int64))


def pds_from_scores(scores: np.ndarray) -> float:
    return float(np.mean(pds_per_target_from_scores(scores)))


def primary_fold_mean(values: Mapping[str, float]) -> float:
    missing = [fold for fold in PRIMARY_FOLDS if fold not in values]
    if missing:
        raise ValueError(f"missing primary folds: {missing}")
    return float(np.mean([float(values[fold]) for fold in PRIMARY_FOLDS]))


@dataclass(frozen=True, order=True)
class Candidate:
    rank: int
    beta: float
    temperature: float
    regularization: float

    @property
    def key(self) -> str:
        return (
            f"listwise:r={self.rank}:beta={self.beta:g}:"
            f"temp={self.temperature:g}:reg={self.regularization:g}"
        )


@dataclass(frozen=True)
class Adapter:
    candidate: Candidate
    basis: np.ndarray
    u: np.ndarray
    v: np.ndarray
    fit_audit: Mapping[str, object]


@dataclass(frozen=True)
class FeatureAxes:
    train_score_mask: np.ndarray
    panel_score_official: np.ndarray
    correction_mask: np.ndarray
    correction_in_train_score: np.ndarray
    correction_in_panel_score: np.ndarray
    train_score_genes: tuple[str, ...]
    panel_score_genes: tuple[str, ...]
    correction_genes: tuple[str, ...]


@dataclass(frozen=True)
class TrainingPair:
    base_scores: np.ndarray
    query_projection: np.ndarray
    truth_projection: np.ndarray
    labels: np.ndarray


@dataclass(frozen=True)
class FitCache:
    basis: np.ndarray
    query_projection: np.ndarray
    pairs: Mapping[str, TrainingPair]
    rows_sha256: str
    correction_sha256: str
    basis_seed: int


@dataclass(frozen=True)
class PermutationPlan:
    replicate: int
    root_seed: int
    global_labels: Mapping[str, np.ndarray]
    labels: Mapping[str, Mapping[str, np.ndarray]]
    priority_sha256: Mapping[str, str]
    plan_sha256: str


def build_candidates(
    ranks: Sequence[int],
    betas: Sequence[float],
    temperatures: Sequence[float],
    regularizations: Sequence[float],
) -> list[Candidate]:
    candidates = sorted(
        {
            Candidate(int(rank), float(beta), float(temperature), float(regularization))
            for rank in ranks
            for beta in betas
            for temperature in temperatures
            for regularization in regularizations
        }
    )
    if not candidates:
        raise ValueError("candidate grid is empty")
    if any(candidate.rank not in (4, 8, 16) for candidate in candidates):
        raise ValueError("adapter ranks must be drawn from {4,8,16}")
    if any(candidate.rank > 16 for candidate in candidates):
        raise ValueError("adapter rank exceeds the predeclared limit")
    if any(not (0.0 < candidate.beta <= 0.20) for candidate in candidates):
        raise ValueError("beta must be in (0,0.20]")
    if any(candidate.temperature <= 0 for candidate in candidates):
        raise ValueError("temperature must be positive")
    if any(candidate.regularization <= 0 for candidate in candidates):
        raise ValueError("regularization must be positive")
    return candidates


def deterministic_target_split(
    targets: Sequence[str],
    *,
    seed: int,
    fit_fraction: float = 0.60,
    tune_fraction: float = 0.20,
) -> dict[str, np.ndarray]:
    if len(targets) < 15:
        raise ValueError("target split requires at least 15 perturbations")
    if not (0 < fit_fraction < 1 and 0 < tune_fraction < 1):
        raise ValueError("split fractions must lie in (0,1)")
    if fit_fraction + tune_fraction >= 1:
        raise ValueError("fit+tune fractions must leave a test partition")
    ordered = sorted(
        range(len(targets)),
        key=lambda index: hashlib.sha256(f"{seed}:{targets[index]}".encode()).digest(),
    )
    n_fit = int(round(len(targets) * fit_fraction))
    n_tune = int(round(len(targets) * tune_fraction))
    split = {
        "fit": np.asarray(sorted(ordered[:n_fit]), dtype=np.int64),
        "tune": np.asarray(sorted(ordered[n_fit : n_fit + n_tune]), dtype=np.int64),
        "test": np.asarray(sorted(ordered[n_fit + n_tune :]), dtype=np.int64),
    }
    sets = [set(indices.tolist()) for indices in split.values()]
    if any(left & right for i, left in enumerate(sets) for right in sets[i + 1 :]):
        raise AssertionError("target partitions overlap")
    if set().union(*sets) != set(range(len(targets))):
        raise AssertionError("target partitions are not exhaustive")
    return split


def randomized_right_basis(
    matrix: np.ndarray,
    rank: int,
    seed: int,
    *,
    oversample: int = 8,
    power_iterations: int = 1,
) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    width = min(rank + oversample, min(values.shape))
    if width < rank:
        raise ValueError(f"basis rank {rank} exceeds matrix shape {values.shape}")
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((values.shape[1], width))
    projected = values @ omega
    for _ in range(power_iterations):
        q, _ = np.linalg.qr(projected, mode="reduced")
        projected = values @ (values.T @ q)
    q, _ = np.linalg.qr(projected, mode="reduced")
    small = q.T @ values
    _, _, vt = np.linalg.svd(small, full_matrices=False)
    basis = vt[:rank].T
    gram = basis.T @ basis
    np.testing.assert_allclose(gram, np.eye(rank), atol=1e-8)
    return basis


def _softmax_cross_entropy_gradient(
    pair: TrainingPair,
    u: np.ndarray,
    v: np.ndarray,
    *,
    beta: float,
    temperature: float,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return listwise loss and gradients without materializing full gene vectors."""
    x = pair.query_projection
    t = pair.truth_projection
    latent = x @ v
    residual = latent @ u.T
    residual_scores = residual @ t.T
    cross = np.sum(x * residual, axis=1)
    residual_norm_sq = np.sum(residual * residual, axis=1)
    denominator = np.sqrt(
        np.maximum(1.0 + 2.0 * beta * cross + beta**2 * residual_norm_sq, 1e-18)
    )
    numerators = pair.base_scores + beta * residual_scores
    scores = numerators / denominator[:, None]
    logits = scores / temperature
    logits -= logits.max(axis=1, keepdims=True)
    exp_logits = np.exp(logits)
    probabilities = exp_logits / exp_logits.sum(axis=1, keepdims=True)
    n = len(pair.labels)
    loss = float(
        np.mean(np.log(exp_logits.sum(axis=1)) - logits[np.arange(n), pair.labels])
    )
    gradient_scores = probabilities
    gradient_scores[np.arange(n), pair.labels] -= 1.0
    gradient_scores /= n * temperature

    direct = (gradient_scores @ t) * (beta / denominator[:, None])
    row_weight = np.sum(gradient_scores * numerators, axis=1)
    norm_gradient = (row_weight / np.maximum(denominator**3, 1e-18))[:, None] * (
        beta * x + beta**2 * residual
    )
    gradient_residual = direct - norm_gradient
    gradient_u = gradient_residual.T @ latent
    gradient_v = x.T @ (gradient_residual @ u)
    return loss, gradient_u, gradient_v


def permute_labels_within_norm_deciles(
    truth: np.ndarray,
    *,
    seed: int,
    deciles: int = 10,
) -> np.ndarray:
    priorities = np.random.default_rng(seed).random(len(truth))
    return permutation_labels_from_priorities(
        truth,
        priorities=priorities,
        deciles=deciles,
    )


def permutation_labels_from_priorities(
    truth: np.ndarray,
    *,
    priorities: np.ndarray,
    deciles: int = 10,
) -> np.ndarray:
    """Build a stratified permutation from one frozen per-target random ordering."""
    norms = np.linalg.norm(np.asarray(truth, dtype=np.float64), axis=1)
    priorities = np.asarray(priorities, dtype=np.float64)
    if priorities.shape != norms.shape:
        raise ValueError("permutation priorities and truth rows differ")
    if not np.isfinite(priorities).all():
        raise ValueError("permutation priorities contain non-finite values")
    order = np.argsort(norms, kind="stable")
    labels = np.arange(len(norms), dtype=np.int64)
    for bucket in np.array_split(order, min(deciles, len(order))):
        labels[bucket] = bucket[np.argsort(priorities[bucket], kind="stable")]
    if sorted(labels.tolist()) != list(range(len(labels))):
        raise AssertionError("norm-stratified permutation is not bijective")
    return labels


def build_permutation_plan(
    truths: Mapping[str, np.ndarray],
    splits: Mapping[str, np.ndarray],
    *,
    replicate: int,
    seed: int,
) -> PermutationPlan:
    """Pre-generate one coherent null world for every later fit stage."""
    if replicate < 0:
        raise ValueError("permutation replicate must be nonnegative")
    required = {
        "fit": np.asarray(splits["fit"], dtype=np.int64),
        "fit_tune": np.sort(
            np.concatenate(
                [
                    np.asarray(splits["fit"], dtype=np.int64),
                    np.asarray(splits["tune"], dtype=np.int64),
                ]
            )
        ),
        "all": np.sort(
            np.concatenate(
                [
                    np.asarray(splits["fit"], dtype=np.int64),
                    np.asarray(splits["tune"], dtype=np.int64),
                    np.asarray(splits["test"], dtype=np.int64),
                ]
            )
        ),
    }
    labels: dict[str, dict[str, np.ndarray]] = {stage: {} for stage in required}
    global_labels: dict[str, np.ndarray] = {}
    priority_hashes: dict[str, str] = {}
    label_hashes: dict[str, dict[str, str]] = {stage: {} for stage in required}
    root_seed = stable_seed(seed, "coherent-null", replicate)
    for teacher, truth0 in sorted(truths.items()):
        truth = np.asarray(truth0, dtype=np.float64)
        priorities = np.random.default_rng(
            stable_seed(root_seed, teacher, "target-priorities")
        ).random(len(truth))
        priority_hashes[teacher] = sha256_array(priorities)
        global_mapping = permutation_labels_from_priorities(
            truth,
            priorities=priorities,
        )
        global_mapping.setflags(write=False)
        global_labels[teacher] = global_mapping
        for stage, rows in required.items():
            stage_labels = restrict_permutation(global_mapping, rows)
            stage_labels.setflags(write=False)
            labels[stage][teacher] = stage_labels
            label_hashes[stage][teacher] = sha256_array(stage_labels)
    audit = {
        "replicate": replicate,
        "root_seed": root_seed,
        "priority_sha256": priority_hashes,
        "global_label_sha256": {
            teacher: sha256_array(value) for teacher, value in global_labels.items()
        },
        "label_sha256": label_hashes,
        "stage_row_sha256": {
            stage: sha256_array(rows) for stage, rows in required.items()
        },
    }
    return PermutationPlan(
        replicate=replicate,
        root_seed=root_seed,
        global_labels=global_labels,
        labels=labels,
        priority_sha256=priority_hashes,
        plan_sha256=sha256_json(audit),
    )


def restrict_permutation(
    global_labels: np.ndarray,
    rows: np.ndarray,
) -> np.ndarray:
    """Restrict permutation cycles to a subset without drawing new randomness."""
    mapping = np.asarray(global_labels, dtype=np.int64)
    rows = np.asarray(rows, dtype=np.int64)
    if sorted(mapping.tolist()) != list(range(len(mapping))):
        raise ValueError("global permutation is not bijective")
    selected = set(map(int, rows))
    local = {int(row): index for index, row in enumerate(rows)}
    result = np.empty(len(rows), dtype=np.int64)
    for index, row0 in enumerate(rows):
        row = int(mapping[int(row0)])
        while row not in selected:
            row = int(mapping[row])
        result[index] = local[row]
    if sorted(result.tolist()) != list(range(len(rows))):
        raise AssertionError("restricted permutation is not bijective")
    return result


def save_permutation_plans(
    path: Path,
    plans: Sequence[PermutationPlan],
) -> dict[str, object]:
    """Freeze exact per-stage labels before any model selection starts."""
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {
        f"r{plan.replicate:04d}__{stage}__{teacher}": np.asarray(labels)
        for plan in plans
        for stage, by_teacher in plan.labels.items()
        for teacher, labels in by_teacher.items()
    }
    arrays.update(
        {
            f"r{plan.replicate:04d}__global__{teacher}": np.asarray(labels)
            for plan in plans
            for teacher, labels in plan.global_labels.items()
        }
    )
    np.savez(path, **arrays)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "replicates": len(plans),
        "plans": [
            {
                "replicate": plan.replicate,
                "root_seed": plan.root_seed,
                "plan_sha256": plan.plan_sha256,
                "priority_sha256": dict(plan.priority_sha256),
                "label_sha256": {
                    stage: {
                        teacher: sha256_array(labels)
                        for teacher, labels in by_teacher.items()
                    }
                    for stage, by_teacher in plan.labels.items()
                },
            }
            for plan in plans
        ],
    }


def prepare_fit_cache(
    query: np.ndarray,
    truths: Mapping[str, np.ndarray],
    correction_in_score: np.ndarray,
    rows: np.ndarray,
    *,
    basis_dim: int,
    basis_seed: int,
) -> FitCache:
    """Cache candidate-invariant normalization, basis, projections, and scores."""
    if not truths:
        raise ValueError("adapter fit requires at least one teacher donor")
    rows = np.asarray(rows, dtype=np.int64)
    q = row_normalize(np.asarray(query[rows], dtype=np.float64))
    normalized_truths = {
        name: row_normalize(np.asarray(values[rows], dtype=np.float64))
        for name, values in truths.items()
    }
    if np.any(np.linalg.norm(q, axis=1) < 0.99):
        raise ValueError("zero query row entered adapter fit")
    if any(
        np.any(np.linalg.norm(value, axis=1) < 0.99)
        for value in normalized_truths.values()
    ):
        raise ValueError("zero truth row entered adapter fit")
    correction = np.asarray(correction_in_score, dtype=np.int64)
    basis_rank = min(
        int(basis_dim),
        len(correction),
        len(rows) * (1 + len(normalized_truths)),
    )
    basis_source = np.vstack(
        [
            q[:, correction],
            *(truth[:, correction] for truth in normalized_truths.values()),
        ]
    )
    basis = randomized_right_basis(basis_source, basis_rank, basis_seed)
    query_projection = q[:, correction] @ basis
    pairs = {
        teacher: TrainingPair(
            base_scores=q @ truth.T,
            query_projection=query_projection,
            truth_projection=truth[:, correction] @ basis,
            labels=np.arange(len(rows), dtype=np.int64),
        )
        for teacher, truth in normalized_truths.items()
    }
    return FitCache(
        basis=basis,
        query_projection=query_projection,
        pairs=pairs,
        rows_sha256=sha256_array(rows),
        correction_sha256=sha256_array(correction),
        basis_seed=int(basis_seed),
    )


def fit_adapter(
    query: np.ndarray,
    truths: Mapping[str, np.ndarray],
    correction_in_score: np.ndarray,
    candidate: Candidate,
    rows: np.ndarray,
    *,
    basis_dim: int,
    epochs: int,
    learning_rate: float,
    gradient_clip: float,
    seed: int,
    basis_seed: int | None = None,
    label_permutations: Mapping[str, np.ndarray] | None = None,
    prepared_cache: FitCache | None = None,
) -> Adapter:
    correction = np.asarray(correction_in_score, dtype=np.int64)
    rows = np.asarray(rows, dtype=np.int64)
    cache = prepared_cache or prepare_fit_cache(
        query,
        truths,
        correction,
        rows,
        basis_dim=basis_dim,
        basis_seed=seed if basis_seed is None else basis_seed,
    )
    if cache.rows_sha256 != sha256_array(rows):
        raise ValueError("prepared fit cache row contract changed")
    if cache.correction_sha256 != sha256_array(correction):
        raise ValueError("prepared fit cache correction axis changed")
    if tuple(cache.pairs) != tuple(truths):
        raise ValueError("prepared fit cache teacher order changed")
    basis = cache.basis
    basis_rank = basis.shape[1]
    if basis_rank < candidate.rank:
        raise ValueError("feature basis is smaller than adapter rank")
    pairs: list[TrainingPair] = []
    permutation_hashes: dict[str, str] = {}
    for teacher, cached_pair in cache.pairs.items():
        labels = (
            np.asarray(label_permutations[teacher], dtype=np.int64)
            if label_permutations is not None
            else np.arange(len(rows), dtype=np.int64)
        )
        if labels.shape != (len(rows),):
            raise ValueError(f"{teacher}: frozen permutation label shape changed")
        if sorted(labels.tolist()) != list(range(len(rows))):
            raise ValueError(f"{teacher}: frozen permutation is not bijective")
        permutation_hashes[teacher] = sha256_array(labels)
        pairs.append(
            TrainingPair(
                base_scores=cached_pair.base_scores,
                query_projection=cached_pair.query_projection,
                truth_projection=cached_pair.truth_projection,
                labels=labels,
            )
        )

    rng = np.random.default_rng(seed)
    initial = rng.standard_normal((basis_rank, candidate.rank))
    v, _ = np.linalg.qr(initial, mode="reduced")
    u = np.zeros((basis_rank, candidate.rank), dtype=np.float64)
    first_u = np.zeros_like(u)
    second_u = np.zeros_like(u)
    first_v = np.zeros_like(v)
    second_v = np.zeros_like(v)
    beta1, beta2 = 0.9, 0.999
    losses: list[float] = []
    gradient_clipped_steps = 0
    for step in range(1, epochs + 1):
        loss = 0.0
        gradient_u = np.zeros_like(u)
        gradient_v = np.zeros_like(v)
        for pair in pairs:
            pair_loss, pair_gu, pair_gv = _softmax_cross_entropy_gradient(
                pair,
                u,
                v,
                beta=candidate.beta,
                temperature=candidate.temperature,
            )
            loss += pair_loss
            gradient_u += pair_gu
            gradient_v += pair_gv
        loss /= len(pairs)
        gradient_u /= len(pairs)
        gradient_v /= len(pairs)

        gram_u = u.T @ u
        gram_v = v.T @ v
        residual_frobenius_sq = float(np.sum(gram_u * gram_v))
        loss += 0.5 * candidate.regularization * residual_frobenius_sq
        gradient_u += candidate.regularization * (u @ gram_v)
        gradient_v += candidate.regularization * (v @ gram_u)

        gradient_norm = math.sqrt(
            float(np.sum(gradient_u * gradient_u) + np.sum(gradient_v * gradient_v))
        )
        if gradient_norm > gradient_clip:
            gradient_clipped_steps += 1
            scale = gradient_clip / gradient_norm
            gradient_u *= scale
            gradient_v *= scale
        first_u = beta1 * first_u + (1 - beta1) * gradient_u
        second_u = beta2 * second_u + (1 - beta2) * (gradient_u * gradient_u)
        first_v = beta1 * first_v + (1 - beta1) * gradient_v
        second_v = beta2 * second_v + (1 - beta2) * (gradient_v * gradient_v)
        correction1 = 1 - beta1**step
        correction2 = 1 - beta2**step
        u -= (
            learning_rate
            * (first_u / correction1)
            / (np.sqrt(second_u / correction2) + 1e-8)
        )
        v -= (
            learning_rate
            * (first_v / correction1)
            / (np.sqrt(second_v / correction2) + 1e-8)
        )
        if (
            not np.isfinite(loss)
            or not np.isfinite(u).all()
            or not np.isfinite(v).all()
        ):
            raise FloatingPointError(f"{candidate.key}: non-finite optimizer state")
        losses.append(loss)

    product_frobenius = math.sqrt(float(np.sum((u.T @ u) * (v.T @ v))))
    return Adapter(
        candidate=candidate,
        basis=basis,
        u=u,
        v=v,
        fit_audit={
            "rows": int(len(rows)),
            "teachers": list(truths),
            "basis_dim": int(basis_rank),
            "epochs": int(epochs),
            "learning_rate": float(learning_rate),
            "gradient_clip": float(gradient_clip),
            "gradient_clipped_steps": int(gradient_clipped_steps),
            "first_loss": float(losses[0]),
            "last_loss": float(losses[-1]),
            "minimum_loss": float(min(losses)),
            "residual_product_frobenius": product_frobenius,
            "regularization_objective": "0.5 * lambda * ||U V^T||_F^2",
            "factor_penalty": False,
            "labels_permuted_within_norm_deciles": label_permutations is not None,
            "label_permutation_sha256": permutation_hashes,
            "prepared_cache": True,
            "prepared_cache_rows_sha256": cache.rows_sha256,
            "prepared_cache_correction_sha256": cache.correction_sha256,
            "prepared_cache_basis_seed": cache.basis_seed,
            "prepared_cache_basis_sha256": sha256_array(cache.basis),
        },
    )


def fit_adapter_batch(
    candidates: Sequence[Candidate],
    cache: FitCache,
    *,
    epochs: int,
    learning_rate: float,
    gradient_clip: float,
    seed: int,
    label_permutations: Mapping[str, np.ndarray] | None = None,
) -> list[Adapter]:
    """Fit same-rank candidates together while preserving independent Adam state."""
    candidates = list(candidates)
    if not candidates:
        return []
    ranks = {candidate.rank for candidate in candidates}
    if len(ranks) != 1:
        raise ValueError("batched candidates must share one adapter rank")
    rank = candidates[0].rank
    basis_rank = cache.basis.shape[1]
    if basis_rank < rank:
        raise ValueError("feature basis is smaller than adapter rank")
    pairs: list[TrainingPair] = []
    permutation_hashes: dict[str, str] = {}
    for teacher, cached_pair in cache.pairs.items():
        labels = (
            np.asarray(label_permutations[teacher], dtype=np.int64)
            if label_permutations is not None
            else np.asarray(cached_pair.labels, dtype=np.int64)
        )
        if labels.shape != cached_pair.labels.shape:
            raise ValueError(f"{teacher}: frozen permutation label shape changed")
        if sorted(labels.tolist()) != list(range(len(labels))):
            raise ValueError(f"{teacher}: frozen permutation is not bijective")
        permutation_hashes[teacher] = sha256_array(labels)
        pairs.append(
            TrainingPair(
                base_scores=cached_pair.base_scores,
                query_projection=cached_pair.query_projection,
                truth_projection=cached_pair.truth_projection,
                labels=labels,
            )
        )

    batch = len(candidates)
    rng = np.random.default_rng(seed)
    initial = rng.standard_normal((basis_rank, rank))
    initial_v, _ = np.linalg.qr(initial, mode="reduced")
    v = np.broadcast_to(initial_v, (batch, basis_rank, rank)).copy()
    u = np.zeros_like(v)
    first_u = np.zeros_like(u)
    second_u = np.zeros_like(u)
    first_v = np.zeros_like(v)
    second_v = np.zeros_like(v)
    betas = np.asarray([candidate.beta for candidate in candidates], dtype=np.float64)
    temperatures = np.asarray(
        [candidate.temperature for candidate in candidates], dtype=np.float64
    )
    regularizations = np.asarray(
        [candidate.regularization for candidate in candidates], dtype=np.float64
    )
    first_losses: np.ndarray | None = None
    minimum_losses = np.full(batch, np.inf, dtype=np.float64)
    last_losses = np.full(batch, np.nan, dtype=np.float64)
    gradient_clipped_steps = np.zeros(batch, dtype=np.int64)
    beta1, beta2 = 0.9, 0.999

    for step in range(1, epochs + 1):
        gradient_u = np.zeros_like(u)
        gradient_v = np.zeros_like(v)
        losses = np.zeros(batch, dtype=np.float64)
        x = pairs[0].query_projection
        latent = np.einsum("nk,bkr->bnr", x, v, optimize=True)
        residual = np.einsum("bnr,bkr->bnk", latent, u, optimize=True)
        cross = np.einsum("nk,bnk->bn", x, residual, optimize=True)
        residual_norm_sq = np.einsum("bnk,bnk->bn", residual, residual, optimize=True)
        denominator = np.sqrt(
            np.maximum(
                1.0
                + 2.0 * betas[:, None] * cross
                + betas[:, None] ** 2 * residual_norm_sq,
                1e-18,
            )
        )
        for pair in pairs:
            residual_scores = np.einsum(
                "bnk,mk->bnm",
                residual,
                pair.truth_projection,
                optimize=True,
            )
            numerators = (
                pair.base_scores[None, :, :] + betas[:, None, None] * residual_scores
            )
            scores = numerators / denominator[:, :, None]
            logits = scores / temperatures[:, None, None]
            logits -= logits.max(axis=2, keepdims=True)
            exp_logits = np.exp(logits)
            probabilities = exp_logits / exp_logits.sum(axis=2, keepdims=True)
            n = len(pair.labels)
            losses += np.mean(
                np.log(exp_logits.sum(axis=2)) - logits[:, np.arange(n), pair.labels],
                axis=1,
            )
            probabilities[:, np.arange(n), pair.labels] -= 1.0
            gradient_scores = probabilities / (n * temperatures[:, None, None])
            direct = np.einsum(
                "bnm,mk->bnk",
                gradient_scores,
                pair.truth_projection,
                optimize=True,
            )
            direct *= betas[:, None, None] / denominator[:, :, None]
            row_weight = np.sum(gradient_scores * numerators, axis=2)
            norm_gradient = (row_weight / np.maximum(denominator**3, 1e-18))[
                :, :, None
            ] * (
                betas[:, None, None] * x[None, :, :]
                + betas[:, None, None] ** 2 * residual
            )
            gradient_residual = direct - norm_gradient
            gradient_u += np.einsum(
                "bnk,bnr->bkr",
                gradient_residual,
                latent,
                optimize=True,
            )
            gradient_v += np.einsum(
                "nk,bnr->bkr",
                x,
                np.einsum(
                    "bnk,bkr->bnr",
                    gradient_residual,
                    u,
                    optimize=True,
                ),
                optimize=True,
            )
        losses /= len(pairs)
        gradient_u /= len(pairs)
        gradient_v /= len(pairs)

        gram_u = np.einsum("bkr,bks->brs", u, u, optimize=True)
        gram_v = np.einsum("bkr,bks->brs", v, v, optimize=True)
        residual_frobenius_sq = np.sum(gram_u * gram_v, axis=(1, 2))
        losses += 0.5 * regularizations * residual_frobenius_sq
        gradient_u += regularizations[:, None, None] * np.matmul(u, gram_v)
        gradient_v += regularizations[:, None, None] * np.matmul(v, gram_u)

        gradient_norm = np.sqrt(
            np.sum(gradient_u * gradient_u, axis=(1, 2))
            + np.sum(gradient_v * gradient_v, axis=(1, 2))
        )
        gradient_clipped_steps += gradient_norm > gradient_clip
        scale = np.minimum(
            1.0,
            gradient_clip / np.maximum(gradient_norm, 1e-30),
        )
        gradient_u *= scale[:, None, None]
        gradient_v *= scale[:, None, None]
        first_u = beta1 * first_u + (1 - beta1) * gradient_u
        second_u = beta2 * second_u + (1 - beta2) * (gradient_u * gradient_u)
        first_v = beta1 * first_v + (1 - beta1) * gradient_v
        second_v = beta2 * second_v + (1 - beta2) * (gradient_v * gradient_v)
        correction1 = 1 - beta1**step
        correction2 = 1 - beta2**step
        u -= (
            learning_rate
            * (first_u / correction1)
            / (np.sqrt(second_u / correction2) + 1e-8)
        )
        v -= (
            learning_rate
            * (first_v / correction1)
            / (np.sqrt(second_v / correction2) + 1e-8)
        )
        if (
            not np.isfinite(losses).all()
            or not np.isfinite(u).all()
            or not np.isfinite(v).all()
        ):
            raise FloatingPointError("batched candidate optimizer became non-finite")
        if first_losses is None:
            first_losses = losses.copy()
        minimum_losses = np.minimum(minimum_losses, losses)
        last_losses = losses

    assert first_losses is not None
    product_frobenius = np.sqrt(
        np.sum(
            np.einsum("bkr,bks->brs", u, u, optimize=True)
            * np.einsum("bkr,bks->brs", v, v, optimize=True),
            axis=(1, 2),
        )
    )
    return [
        Adapter(
            candidate=candidate,
            basis=cache.basis,
            u=u[index],
            v=v[index],
            fit_audit={
                "rows": int(pairs[0].base_scores.shape[0]),
                "teachers": list(cache.pairs),
                "basis_dim": int(basis_rank),
                "epochs": int(epochs),
                "learning_rate": float(learning_rate),
                "gradient_clip": float(gradient_clip),
                "gradient_clipped_steps": int(gradient_clipped_steps[index]),
                "first_loss": float(first_losses[index]),
                "last_loss": float(last_losses[index]),
                "minimum_loss": float(minimum_losses[index]),
                "residual_product_frobenius": float(product_frobenius[index]),
                "regularization_objective": "0.5 * lambda * ||U V^T||_F^2",
                "factor_penalty": False,
                "labels_permuted_within_norm_deciles": (label_permutations is not None),
                "label_permutation_sha256": permutation_hashes,
                "prepared_cache": True,
                "prepared_cache_rows_sha256": cache.rows_sha256,
                "prepared_cache_correction_sha256": cache.correction_sha256,
                "prepared_cache_basis_seed": cache.basis_seed,
                "prepared_cache_basis_sha256": sha256_array(cache.basis),
                "same_rank_batch_size": batch,
            },
        )
        for index, candidate in enumerate(candidates)
    ]


def adapter_scores_and_angles(
    query: np.ndarray,
    truth: np.ndarray,
    correction_in_score: np.ndarray,
    adapter: Adapter | None,
) -> tuple[np.ndarray, np.ndarray]:
    q = row_normalize(query)
    t = row_normalize(truth)
    baseline = q @ t.T
    if adapter is None:
        return baseline, np.zeros(len(q), dtype=np.float64)
    candidate = adapter.candidate
    x = q[:, correction_in_score] @ adapter.basis
    truth_projection = t[:, correction_in_score] @ adapter.basis
    latent = x @ adapter.v
    residual = latent @ adapter.u.T
    numerator = baseline + candidate.beta * (residual @ truth_projection.T)
    cross = np.sum(x * residual, axis=1)
    residual_norm_sq = np.sum(residual * residual, axis=1)
    denominator = np.sqrt(
        np.maximum(
            1.0 + 2.0 * candidate.beta * cross + candidate.beta**2 * residual_norm_sq,
            1e-18,
        )
    )
    scores = numerator / denominator[:, None]
    cosine_to_original = np.clip(
        (1.0 + candidate.beta * cross) / denominator,
        -1.0,
        1.0,
    )
    angles = np.degrees(np.arccos(cosine_to_original))
    return scores, angles


def transform_query(
    query: np.ndarray,
    correction_in_score: np.ndarray,
    adapter: Adapter | None,
) -> np.ndarray:
    """Materialize the exact normalized adapter output for deployment/tests."""
    q = row_normalize(query)
    if adapter is None:
        return q
    x = q[:, correction_in_score] @ adapter.basis
    residual = (x @ adapter.v) @ adapter.u.T
    output = q.copy()
    output[:, correction_in_score] += (
        adapter.candidate.beta * residual @ adapter.basis.T
    )
    return row_normalize(output)


def evaluate_adapter(
    query: np.ndarray,
    truth_gallery: np.ndarray,
    correction_in_query: np.ndarray,
    adapter: Adapter | None,
    query_rows: np.ndarray,
    positive_columns: np.ndarray,
    *,
    query_score_columns: np.ndarray | None = None,
    truth_score_columns: np.ndarray | None = None,
) -> dict[str, object]:
    rows = np.asarray(query_rows, dtype=np.int64)
    positives = np.asarray(positive_columns, dtype=np.int64)
    if query_score_columns is None:
        query_score_columns = np.arange(query.shape[1], dtype=np.int64)
    if truth_score_columns is None:
        truth_score_columns = np.arange(truth_gallery.shape[1], dtype=np.int64)
    query_columns = np.asarray(query_score_columns, dtype=np.int64)
    truth_columns = np.asarray(truth_score_columns, dtype=np.int64)
    if len(query_columns) != len(truth_columns):
        raise ValueError("query and truth score axes have different widths")

    original_full = row_normalize(np.asarray(query, dtype=np.float64))
    adapted_full = transform_query(query, correction_in_query, adapter)
    original = row_normalize(original_full[np.ix_(rows, query_columns)])
    transformed = row_normalize(adapted_full[np.ix_(rows, query_columns)])
    truth = row_normalize(np.asarray(truth_gallery[:, truth_columns], dtype=np.float64))
    baseline_scores = original @ truth.T
    adapted_scores = transformed @ truth.T
    baseline = pds_per_query_from_scores(baseline_scores, positives)
    adapted = pds_per_query_from_scores(adapted_scores, positives)
    cosines = np.clip(np.sum(original * transformed, axis=1), -1.0, 1.0)
    angles = np.degrees(np.arccos(cosines))
    return {
        "baseline": float(np.mean(baseline)),
        "adapted": float(np.mean(adapted)),
        "delta": float(np.mean(adapted - baseline)),
        "queries": int(len(rows)),
        "gallery": int(len(truth)),
        "score_genes": int(len(query_columns)),
        "per_target_baseline": baseline,
        "per_target_adapted": adapted,
        "per_target_delta": adapted - baseline,
        "angles_degrees": angles,
        "median_angle_degrees": float(np.median(angles)),
        "max_angle_degrees": float(np.max(angles)),
    }


def choose_candidate(
    cells_by_candidate: Mapping[str, Sequence[Mapping[str, object]]],
    candidates: Sequence[Candidate],
) -> tuple[str, dict[str, object]]:
    by_key = {candidate.key: candidate for candidate in candidates}
    summaries: dict[str, object] = {}
    eligible: list[str] = []
    for candidate in candidates:
        cells = cells_by_candidate[candidate.key]
        deltas = [float(cell["delta"]) for cell in cells]
        angles = [float(cell["median_angle_degrees"]) for cell in cells]
        summary = {
            "cell_deltas": deltas,
            "cell_median_angles_degrees": angles,
            "minimum_delta": float(min(deltas)),
            "macro_delta": float(np.mean(deltas)),
            "maximum_cell_median_angle_degrees": float(max(angles)),
        }
        summaries[candidate.key] = summary
        if (
            min(deltas) >= 0.0
            and np.mean(deltas) > 0.0
            and max(angles) < ANGLE_LIMIT_DEGREES
        ):
            eligible.append(candidate.key)
    if not eligible:
        return IDENTITY, {
            "selected": IDENTITY,
            "eligible_count": 0,
            "candidate_summaries": summaries,
            "rule": (
                "identity unless a candidate is nonnegative in every inner "
                "donor-holdout cell and stays below the 10-degree angle guard"
            ),
        }

    best = max(float(summaries[key]["macro_delta"]) for key in eligible)
    exact_best = [
        key for key in eligible if float(summaries[key]["macro_delta"]) == best
    ]
    selected = min(
        exact_best,
        key=lambda key: (
            by_key[key].rank,
            by_key[key].beta,
            -by_key[key].regularization,
            abs(by_key[key].temperature - 0.07),
            key,
        ),
    )
    return selected, {
        "selected": selected,
        "eligible_count": len(eligible),
        "best_macro_delta": best,
        "exact_tie_count": len(exact_best),
        "candidate_summaries": summaries,
        "rule": (
            "exactly maximize macro tune PDS among candidates nonnegative in "
            "every outer-role x held-teacher cell and under 10 degrees; only "
            "an exact floating-point tie uses the deterministic lower-complexity "
            "ordering"
        ),
    }


def source_effect_clusters(
    queries: Sequence[np.ndarray],
    *,
    reciprocal_k: int,
    cosine_threshold: float,
) -> tuple[np.ndarray, dict[str, object]]:
    if not queries:
        raise ValueError("source clustering requires at least one query matrix")
    units = [row_normalize(query) for query in queries]
    consensus = row_normalize(np.mean(units, axis=0))
    similarities = consensus @ consensus.T
    np.fill_diagonal(similarities, -np.inf)
    k = min(reciprocal_k, len(consensus) - 1)
    neighbors = np.argsort(-similarities, axis=1, kind="stable")[:, :k]
    directed = np.zeros((len(consensus), len(consensus)), dtype=bool)
    directed[np.arange(len(consensus))[:, None], neighbors] = True
    reciprocal = directed & directed.T & (similarities >= cosine_threshold)

    parent = np.arange(len(consensus), dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[max(root_left, root_right)] = min(root_left, root_right)

    rows, columns = np.where(np.triu(reciprocal, k=1))
    for left, right in zip(rows.tolist(), columns.tolist(), strict=True):
        union(left, right)
    roots = [find(index) for index in range(len(consensus))]
    root_order = {root: index for index, root in enumerate(sorted(set(roots)))}
    labels = np.asarray([root_order[root] for root in roots], dtype=np.int64)
    sizes = np.bincount(labels)
    return labels, {
        "method": "source-only reciprocal-kNN connected components",
        "reciprocal_k": int(k),
        "cosine_threshold": float(cosine_threshold),
        "clusters": int(len(sizes)),
        "singleton_clusters": int(np.sum(sizes == 1)),
        "largest_cluster": int(np.max(sizes)),
        "cluster_sizes_sha256": sha256_array(sizes),
    }


def every_fold_angle_below_limit(
    fold_median_angles: Mapping[str, float],
    *,
    limit_degrees: float = ANGLE_LIMIT_DEGREES,
) -> bool:
    return bool(fold_median_angles) and all(
        float(angle) < limit_degrees for angle in fold_median_angles.values()
    )


def load_target_clusters(
    path: Path,
    targets: Sequence[str],
) -> tuple[np.ndarray, dict[str, object]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"target_gene", "cluster"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        mapping: dict[str, str] = {}
        for row in reader:
            target = str(row["target_gene"]).strip()
            cluster = str(row["cluster"]).strip()
            if not target or not cluster:
                continue
            if target in mapping:
                raise ValueError(f"{path}: duplicate target {target}")
            mapping[target] = cluster
    missing = [target for target in targets if target not in mapping]
    extras = sorted(set(mapping) - set(targets))
    if missing or extras:
        raise ValueError(
            f"{path}: cluster coverage mismatch; missing={missing[:5]}, "
            f"extras={extras[:5]}"
        )
    names = sorted(set(mapping.values()))
    cluster_index = {name: index for index, name in enumerate(names)}
    labels = np.asarray(
        [cluster_index[mapping[target]] for target in targets], dtype=np.int64
    )
    sizes = np.bincount(labels)
    return labels, {
        "method": "provided gene-family/pathway clusters",
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "clusters": int(len(names)),
        "singleton_clusters": int(np.sum(sizes == 1)),
        "largest_cluster": int(np.max(sizes)),
        "cluster_names_sha256": sha256_lines(names),
    }


def fold_specific_cluster_bootstrap(
    deltas_by_fold: Mapping[str, np.ndarray],
    cluster_labels_by_fold: Mapping[str, np.ndarray],
    *,
    replicates: int,
    seed: int,
) -> dict[str, object]:
    folds = tuple(deltas_by_fold)
    if not folds or set(folds) != set(cluster_labels_by_fold):
        raise ValueError("bootstrap folds and cluster-label folds differ")
    members_by_fold: dict[str, list[np.ndarray]] = {}
    for fold in folds:
        values = np.asarray(deltas_by_fold[fold], dtype=np.float64)
        labels = np.asarray(cluster_labels_by_fold[fold], dtype=np.int64)
        if values.ndim != 1 or values.shape != labels.shape:
            raise ValueError(f"{fold}: deltas and cluster labels do not align")
        clusters = sorted(set(map(int, labels.tolist())))
        if not clusters:
            raise ValueError(f"{fold}: no bootstrap clusters")
        members_by_fold[fold] = [
            np.flatnonzero(labels == cluster) for cluster in clusters
        ]
    rng = np.random.default_rng(seed)
    draws = np.empty(replicates, dtype=np.float64)
    for replicate in range(replicates):
        fold_means = []
        for fold in folds:
            members = members_by_fold[fold]
            sampled = rng.integers(0, len(members), size=len(members))
            indices = np.concatenate([members[index] for index in sampled])
            fold_means.append(float(np.mean(deltas_by_fold[fold][indices])))
        draws[replicate] = float(np.mean(fold_means))
    return {
        "method": (
            "fold-specific source-query target-cluster bootstrap; independently "
            "resample each fold's clusters, then average fold means"
        ),
        "replicates": int(replicates),
        "seed": int(seed),
        "clusters_by_fold": {fold: int(len(members_by_fold[fold])) for fold in folds},
        "point_estimate": float(
            np.mean([np.mean(deltas_by_fold[fold]) for fold in folds])
        ),
        "lower_95": float(np.quantile(draws, 0.025)),
        "upper_95": float(np.quantile(draws, 0.975)),
        "probability_positive": float(np.mean(draws > 0.0)),
    }


def fold_specific_cluster_dependence(
    deltas_by_fold: Mapping[str, np.ndarray],
    cluster_labels_by_fold: Mapping[str, np.ndarray],
) -> dict[str, object]:
    if set(deltas_by_fold) != set(cluster_labels_by_fold):
        raise ValueError("dependence folds and cluster-label folds differ")
    per_fold: dict[str, object] = {}
    for fold, raw_values in deltas_by_fold.items():
        values = np.asarray(raw_values, dtype=np.float64)
        labels = np.asarray(cluster_labels_by_fold[fold], dtype=np.int64)
        if values.shape != labels.shape:
            raise ValueError(f"{fold}: dependence inputs do not align")
        clusters = sorted(set(map(int, labels.tolist())))
        cluster_sums = np.asarray(
            [float(np.sum(values[labels == cluster])) for cluster in clusters]
        )
        positive = np.maximum(cluster_sums, 0.0)
        positive_total = float(np.sum(positive))
        largest_share = (
            float(np.max(positive) / positive_total) if positive_total > 0 else 1.0
        )
        leave_one_out = [
            float(np.mean(values[labels != cluster])) for cluster in clusters
        ]
        passes = bool(
            len(clusters) >= 2 and largest_share <= 0.50 and min(leave_one_out) >= 0.0
        )
        per_fold[fold] = {
            "clusters": len(clusters),
            "largest_positive_cluster_share": largest_share,
            "leave_one_cluster_out_minimum_delta": (
                float(min(leave_one_out)) if leave_one_out else None
            ),
            "leave_one_cluster_out_maximum_delta": (
                float(max(leave_one_out)) if leave_one_out else None
            ),
            "passes_no_obvious_single_family_dependence": passes,
        }
    return {
        "method": "fold-specific diagnostics on source-only query clusters",
        "folds": per_fold,
        "passes_no_obvious_single_family_dependence": all(
            bool(value["passes_no_obvious_single_family_dependence"])
            for value in per_fold.values()
        ),
    }


def validate_truth_cache(
    truth_cache: Path,
    official: Sequence[str],
    panel: Sequence[str],
    donor_root: Path,
    manifest: Mapping[str, object],
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], FeatureAxes]:
    master = read_lines(truth_cache / "master_genes.txt")
    train_targets = read_lines(truth_cache / "train_perts.txt")
    outer_targets = read_lines(truth_cache / "outer_perts.txt")
    correction_mask = np.asarray(
        np.load(truth_cache / "train_mask.npy", allow_pickle=False), dtype=bool
    )
    if len(train_targets) != 600:
        raise ValueError(f"expected 600 non-panel targets, got {len(train_targets)}")
    if set(train_targets) & set(panel):
        raise ValueError("training targets overlap VCC panel")
    if not set(outer_targets) <= set(panel):
        raise ValueError("outer targets are not a VCC-panel subset")
    if set(train_targets) & set(outer_targets):
        raise ValueError("training and panel targets overlap")
    if len(master) != len(correction_mask):
        raise ValueError("truth-cache correction mask width mismatch")
    if set(master) & set(panel):
        raise ValueError("VCC coordinate entered master feature axis")
    if not set(master) <= set(official):
        raise ValueError("truth-cache feature lies outside official genes")
    train_set = set(train_targets)
    train_score_mask = np.asarray(
        [gene not in train_set for gene in master], dtype=bool
    )
    panel_set = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    donor_axes = []
    for donor in DEVELOPMENT_FOLDS:
        donor_genes = set(read_gene_lines(donor_root / donor / "genes.txt"))
        donor_axes.append(
            tuple(
                gene
                for gene in official
                if gene in donor_genes and gene not in panel_set
            )
        )
    panel_score_genes = donor_axes[0]
    if any(axis != panel_score_genes for axis in donor_axes[1:]):
        raise ValueError("development held-donor score axes are not identical")
    panel_score_official = np.asarray(
        [official_pos[gene] for gene in panel_score_genes], dtype=np.int64
    )
    final_axis = manifest.get("final_score_axis", {})
    expected_count = int(final_axis.get("gene_count", -1))
    if (
        len(official) == 18_533
        and len(panel) == 300
        and expected_count != EXPECTED_FINAL_SCORE_GENES
    ):
        raise ValueError(
            "manifest final score-axis count is not the established 17,806"
        )
    if len(panel_score_genes) != expected_count:
        raise ValueError(
            f"final score axis {len(panel_score_genes)} != frozen {expected_count}"
        )
    if sha256_lines(panel_score_genes) != final_axis.get("gene_order_sha256"):
        raise ValueError("final score-gene order/hash changed")
    if sha256_array(panel_score_official) != final_axis.get(
        "official_positions_sha256"
    ):
        raise ValueError("final score official-position order/hash changed")
    if np.any(
        correction_mask & np.asarray([gene in train_set for gene in master], dtype=bool)
    ):
        raise ValueError("training target coordinate entered correction axis")
    train_score_indices = np.flatnonzero(train_score_mask)
    correction_indices = np.flatnonzero(correction_mask)
    train_score_position = {
        int(value): index for index, value in enumerate(train_score_indices)
    }
    panel_score_position = {gene: index for index, gene in enumerate(panel_score_genes)}
    correction_in_train_score = np.asarray(
        [train_score_position[int(value)] for value in correction_indices],
        dtype=np.int64,
    )
    correction_in_panel_score = np.asarray(
        [panel_score_position[master[int(value)]] for value in correction_indices],
        dtype=np.int64,
    )
    train_score_genes = tuple(
        gene for gene, keep in zip(master, train_score_mask, strict=True) if keep
    )
    correction_genes = tuple(
        gene for gene, keep in zip(master, correction_mask, strict=True) if keep
    )
    if set(panel_score_genes) & set(panel):
        raise AssertionError("VCC coordinate entered final panel score axis")
    if set(train_score_genes) & train_set:
        raise AssertionError("non-panel target coordinate entered training score axis")
    if not (set(panel_score_genes) & train_set):
        raise ValueError(
            "final panel score axis unexpectedly excludes all training targets"
        )
    return (
        master,
        train_targets,
        outer_targets,
        FeatureAxes(
            train_score_mask=train_score_mask,
            panel_score_official=panel_score_official,
            correction_mask=correction_mask,
            correction_in_train_score=correction_in_train_score,
            correction_in_panel_score=correction_in_panel_score,
            train_score_genes=train_score_genes,
            panel_score_genes=panel_score_genes,
            correction_genes=correction_genes,
        ),
    )


def verify_manifest_inputs(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
    *,
    include_outer_truth: bool,
) -> None:
    expected_protocol_contracts = {
        "optimizer_regularization": "0.5 * lambda * ||U V^T||_F^2",
        "factor_penalty": False,
        "selection_objective": "exact_macro_tune_pds",
        "tie_window": None,
        "exact_tie_break": [
            "lower_rank",
            "lower_beta",
            "stronger_regularization",
            "temperature_nearest_0.07",
            "candidate_key",
        ],
        "batch_equivalence_regression": {
            "ranks": [4, 8, 16],
            "teachers": 2,
            "permuted_labels": True,
            "forced_gradient_clipping": True,
        },
        "source_clustering_query": "panel_rows_retained_coordinate_safe",
    }
    if manifest.get("protocol_contracts") != expected_protocol_contracts:
        raise ValueError("artifact manifest protocol contracts changed")
    challenge = manifest["challenge"]
    for path, name in (
        (args.genes_csv, "gene_names.csv"),
        (args.panel_csv, "pert_counts.csv"),
    ):
        actual = sha256_file(path)
        expected = str(challenge[name])
        if actual != expected:
            raise ValueError(f"{name}: SHA-256 {actual} != {expected}")
    for filename, expected0 in manifest["truth_cache"]["files"].items():
        if filename.endswith("_raw_outer.npy") and not include_outer_truth:
            continue
        path = args.truth_cache / filename
        expected = str(expected0)
        actual = sha256_file(path)
        if actual != expected:
            raise ValueError(f"truth-cache/{filename}: SHA-256 {actual} != {expected}")
    if include_outer_truth:
        for donor in DEVELOPMENT_FOLDS:
            for filename, expected0 in manifest["donors"][donor].items():
                path = args.donor_root / donor / filename
                expected = str(expected0)
                actual = sha256_file(path)
                if actual != expected:
                    raise ValueError(
                        f"{donor}/{filename}: SHA-256 {actual} != {expected}"
                    )
    kolf = manifest.get("kolf_confirmation", {})
    if include_outer_truth and bool(kolf.get("supported", False)):
        if args.kolf_root is None:
            raise ValueError("KOLF confirmation is supported but --kolf-root is absent")
        for filename, expected0 in kolf["files"].items():
            path = args.kolf_root / filename
            expected = str(expected0)
            actual = sha256_file(path)
            if actual != expected:
                raise ValueError(f"kolf/{filename}: SHA-256 {actual} != {expected}")


def validate_query_cache(
    query_cache: Path,
    axes: FeatureAxes,
    train_targets: Sequence[str],
    outer_targets: Sequence[str],
    panel: Sequence[str],
    manifest: Mapping[str, object],
) -> dict[str, object]:
    metadata_path = query_cache / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    panel_official_by_gene = dict(
        zip(
            axes.panel_score_genes,
            map(int, axes.panel_score_official),
            strict=True,
        )
    )
    train_score_official = np.asarray(
        [panel_official_by_gene[gene] for gene in axes.train_score_genes],
        dtype=np.int64,
    )
    expected = {
        "schema_version": CACHE_SCHEMA,
        "experiment": "080-listwise-metric-adapter",
        "arm": "D",
        "development_folds": list(DEVELOPMENT_FOLDS),
        "final_folds": list(FINAL_FOLDS),
        "training_targets": len(train_targets),
        "outer_targets": len(outer_targets),
        "training_outer_overlap": 0,
        "vcc_targets": len(panel),
        "train_score_genes": len(axes.train_score_genes),
        "cluster_score_genes": len(axes.train_score_genes),
        "panel_score_genes": len(axes.panel_score_genes),
        "correction_genes": len(axes.correction_genes),
        "train_score_gene_order_sha256": sha256_lines(axes.train_score_genes),
        "cluster_score_gene_order_sha256": sha256_lines(axes.train_score_genes),
        "cluster_score_official_positions_sha256": sha256_array(train_score_official),
        "panel_score_gene_order_sha256": sha256_lines(axes.panel_score_genes),
        "panel_score_official_positions_sha256": sha256_array(
            axes.panel_score_official
        ),
        "correction_gene_order_sha256": sha256_lines(axes.correction_genes),
        "correction_in_train_score_sha256": sha256_array(
            axes.correction_in_train_score
        ),
        "correction_in_panel_score_sha256": sha256_array(
            axes.correction_in_panel_score
        ),
        "train_order_sha256": sha256_lines(train_targets),
        "outer_order_sha256": sha256_lines(outer_targets),
        "panel_order_sha256": sha256_lines(panel),
        "final_score_axis": {
            "gene_count": len(axes.panel_score_genes),
            "gene_order_sha256": sha256_lines(axes.panel_score_genes),
            "official_positions_sha256": sha256_array(axes.panel_score_official),
            "development_donor_axes_equal": True,
            "excludes_exactly_vcc_coordinates_within_held_support": True,
        },
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"query-cache metadata {key}={metadata.get(key)!r}, expected {value!r}"
            )
    for held in FINAL_FOLDS:
        fold = metadata["folds"][held]
        operator = manifest["operators"][held]
        if fold.get("held_donor") != held:
            raise ValueError(f"{held}: cache held-donor label mismatch")
        source_donors = tuple(map(str, fold.get("source_donors", [])))
        if held in set(source_donors):
            raise ValueError(f"{held}: held donor entered query sources")
        if source_donors != tuple(map(str, operator["source_donors"])):
            raise ValueError(f"{held}: cache/operator source-donor mismatch")
        if fold.get("operator_npz_sha256") != operator["npz_sha256"]:
            raise ValueError(f"{held}: frozen operator NPZ hash mismatch")
        if fold.get("operator_json_sha256") != operator["json_sha256"]:
            raise ValueError(f"{held}: frozen operator metadata hash mismatch")
        if fold.get("held_donor_used_in_prediction") is not False:
            raise ValueError(f"{held}: held-donor query audit failed")
        if fold.get("held_donor_used_in_operator") is not False:
            raise ValueError(f"{held}: held-donor operator audit failed")
        expected_panels = (
            ("train", len(train_targets), len(axes.train_score_genes)),
            (
                "cluster",
                len(panel) if held == KOLF_FOLD else len(outer_targets),
                len(axes.train_score_genes),
            ),
            (
                "panel",
                len(panel) if held == KOLF_FOLD else len(outer_targets),
                len(axes.panel_score_genes),
            ),
        )
        for label, rows, columns in expected_panels:
            path = query_cache / held / f"{label}_query.npy"
            values = np.load(path, mmap_mode="r", allow_pickle=False)
            if values.shape != (rows, columns):
                raise ValueError(f"{path}: unexpected shape {values.shape}")
            expected_hash = fold["panels"][label]["query_sha256"]
            if sha256_file(path) != expected_hash:
                raise ValueError(f"{path}: cache hash mismatch")
            panel_audit = fold["panels"][label]
            if label == "panel":
                required_axis = {
                    "score_gene_count": len(axes.panel_score_genes),
                    "score_gene_order_sha256": sha256_lines(axes.panel_score_genes),
                    "score_official_positions_sha256": sha256_array(
                        axes.panel_score_official
                    ),
                    "saved_query_is_complete_scored_axis": True,
                }
                for key, value in required_axis.items():
                    if panel_audit.get(key) != value:
                        raise ValueError(
                            f"{held}: panel query axis audit {key} changed"
                        )
            elif label == "cluster":
                required_axis = {
                    "score_gene_count": len(axes.train_score_genes),
                    "score_gene_order_sha256": sha256_lines(axes.train_score_genes),
                    "saved_query_is_complete_cluster_axis": True,
                }
                for key, value in required_axis.items():
                    if panel_audit.get(key) != value:
                        raise ValueError(
                            f"{held}: cluster query axis audit {key} changed"
                        )
        restriction = fold["panels"]["train"]["operator_restriction"]
        if restriction.get("all_nonpanel_target_coordinates_excluded") is not True:
            raise ValueError(
                f"{held}: a non-panel target entered the training operator"
            )
        if restriction.get("nonpanel_target_overlap_with_active_operator") != 0:
            raise ValueError(f"{held}: non-panel/operator overlap is nonzero")
        if restriction.get("dense_zero_equivalence_proved") is not True:
            raise ValueError(f"{held}: dense-zero equivalence was not proved")
        train_support = fold["panels"]["train"].get("support", {})
        cluster_support = fold["panels"]["cluster"].get("support", {})
        panel_support = fold["panels"]["panel"].get("support", {})
        if train_support.get("used_for_adapter_fit_or_selection") is not True:
            raise ValueError(f"{held}: training feature-role audit changed")
        if panel_support.get("inference_only_fixed_baseline") is not True:
            raise ValueError(f"{held}: panel baseline-role audit changed")
        if panel_support.get("used_for_adapter_fit_or_selection") is not False:
            raise ValueError(f"{held}: panel baseline entered a learned statistic")
        if cluster_support.get("used_for_source_clustering") is not True:
            raise ValueError(f"{held}: source clustering query role changed")
        if cluster_support.get("used_for_adapter_fit_or_selection") is not False:
            raise ValueError(f"{held}: cluster query entered adapter fitting/selection")
        cluster_restriction = fold["panels"]["cluster"]["operator_restriction"]
        if (
            cluster_restriction.get(
                "all_vcc_and_historical_target_coordinates_excluded"
            )
            is not True
        ):
            raise ValueError(f"{held}: cluster query retained a target coordinate")
        if held in DEVELOPMENT_FOLDS:
            reproduction = fold.get("exp064_reproduction", {})
            if reproduction.get("passed") is not True:
                raise ValueError(f"{held}: exp064 arm-D reproduction did not pass")
            if (
                reproduction.get("independent_from_staged_panel_preparation")
                is not True
            ):
                raise ValueError(
                    f"{held}: exp064 comparison was not independently reconstructed"
                )
            if float(reproduction.get("max_abs_error", math.inf)) > float(
                manifest["exp064_reference"]["required_max_abs_error"]
            ):
                raise ValueError(f"{held}: exp064 reconstruction tolerance failed")
            if reproduction.get("saved_query_score_axis_exact") is not True:
                raise ValueError(
                    f"{held}: exp064 audit did not cover saved panel_query"
                )
            if reproduction.get("retrieval_genes") != len(axes.panel_score_genes):
                raise ValueError(f"{held}: exp064 retrieval-axis width changed")
            if reproduction.get("retrieval_gene_order_sha256") != sha256_lines(
                axes.panel_score_genes
            ):
                raise ValueError(f"{held}: exp064 retrieval gene order changed")
    return metadata


def load_truth(
    truth_cache: Path,
    donor: str,
    panel: str,
    score_mask: np.ndarray,
) -> np.ndarray:
    path = truth_cache / f"{donor}_raw_{panel}.npy"
    values = np.load(path, mmap_mode="r", allow_pickle=False)
    if values.shape[1] != len(score_mask):
        raise ValueError(f"{path}: feature width mismatch")
    selected = np.asarray(values[:, score_mask], dtype=np.float64)
    if not np.isfinite(selected).all():
        raise ValueError(f"{path}: non-finite truth value")
    return selected


def load_compact_panel_truth(
    donor_root: Path,
    donor: str,
    targets: Sequence[str],
    score_genes: Sequence[str],
    *,
    chunk_rows: int,
) -> tuple[np.ndarray, dict[str, object]]:
    """Load a held donor on the complete frozen score axis after authentication."""
    root = donor_root / donor
    genes = read_gene_lines(root / "genes.txt")
    perts = read_lines(root / "perts.txt")
    gene_last = {gene: index for index, gene in enumerate(genes)}
    missing_genes = [gene for gene in score_genes if gene not in gene_last]
    if missing_genes:
        raise ValueError(f"{donor}: {len(missing_genes)} frozen score genes are absent")
    pert_pos = {target: index for index, target in enumerate(perts)}
    missing_targets = [target for target in targets if target not in pert_pos]
    if missing_targets:
        raise ValueError(f"{donor}: {len(missing_targets)} held targets are absent")
    columns = np.asarray([gene_last[gene] for gene in score_genes], dtype=np.int64)
    rows = np.asarray([pert_pos[target] for target in targets], dtype=np.int64)
    delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape != (len(perts), len(genes)):
        raise ValueError(f"{donor}: compact delta shape does not match its axes")
    center = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, len(perts), chunk_rows):
        hi = min(lo + chunk_rows, len(perts))
        center += np.asarray(delta[lo:hi, columns], dtype=np.float64).sum(
            axis=0,
            dtype=np.float64,
        )
    center /= len(perts)
    values = np.asarray(delta[np.ix_(rows, columns)], dtype=np.float64) - center
    mmap = getattr(delta, "_mmap", None)
    if mmap is not None:
        mmap.close()
    if not np.isfinite(values).all():
        raise ValueError(f"{donor}: held truth contains non-finite values")
    return values, {
        "targets": len(targets),
        "score_genes": len(score_genes),
        "target_order_sha256": sha256_lines(targets),
        "score_gene_order_sha256": sha256_lines(score_genes),
        "representation": "raw donor delta centered over every donor perturbation",
        "used_for_selection": False,
        "used_for_fitting": False,
    }


def summarize_evaluation(result: Mapping[str, object]) -> dict[str, float]:
    return {
        "baseline": float(result["baseline"]),
        "adapted": float(result["adapted"]),
        "delta": float(result["delta"]),
        "queries": int(result["queries"]),
        "gallery": int(result["gallery"]),
        "score_genes": int(result["score_genes"]),
        "median_angle_degrees": float(result["median_angle_degrees"]),
        "max_angle_degrees": float(result["max_angle_degrees"]),
    }


def prepare_all_fit_caches(
    *,
    train_queries: Mapping[str, np.ndarray],
    auxiliary_truth: Mapping[str, np.ndarray],
    correction_in_train_score: np.ndarray,
    splits: Mapping[str, np.ndarray],
    args: argparse.Namespace,
) -> tuple[
    dict[tuple[str, str], FitCache],
    dict[str, FitCache],
    dict[str, FitCache],
    dict[str, object],
]:
    """Build the nine candidate- and permutation-invariant fit contexts once."""
    fit_rows = np.asarray(splits["fit"], dtype=np.int64)
    fit_tune_rows = np.sort(
        np.concatenate(
            [
                np.asarray(splits["fit"], dtype=np.int64),
                np.asarray(splits["tune"], dtype=np.int64),
            ]
        )
    )
    all_rows = np.sort(
        np.concatenate(
            [
                np.asarray(splits["fit"], dtype=np.int64),
                np.asarray(splits["tune"], dtype=np.int64),
                np.asarray(splits["test"], dtype=np.int64),
            ]
        )
    )
    inner: dict[tuple[str, str], FitCache] = {}
    for outer in DEVELOPMENT_FOLDS:
        for held_teacher in TEACHER_DONORS:
            fit_teachers = {
                donor: auxiliary_truth[donor]
                for donor in TEACHER_DONORS
                if donor != held_teacher
            }
            inner[(outer, held_teacher)] = prepare_fit_cache(
                train_queries[outer],
                fit_teachers,
                correction_in_train_score,
                fit_rows,
                basis_dim=args.basis_dim,
                basis_seed=stable_seed(
                    args.optimizer_seed,
                    "inner-basis",
                    outer,
                    held_teacher,
                ),
            )
    pretest = {
        outer: prepare_fit_cache(
            train_queries[outer],
            auxiliary_truth,
            correction_in_train_score,
            fit_tune_rows,
            basis_dim=args.basis_dim,
            basis_seed=stable_seed(
                args.optimizer_seed,
                "pretest-basis",
                outer,
            ),
        )
        for outer in DEVELOPMENT_FOLDS
    }
    final = {
        fold: prepare_fit_cache(
            train_queries[fold],
            auxiliary_truth,
            correction_in_train_score,
            all_rows,
            basis_dim=args.basis_dim,
            basis_seed=stable_seed(args.optimizer_seed, "final-basis", fold),
        )
        for fold in FINAL_FOLDS
    }
    all_caches = [*inner.values(), *pretest.values(), *final.values()]
    audit = {
        "contexts": len(all_caches),
        "expected_contexts": 9,
        "candidate_invariant_values": [
            "row normalization",
            "randomized SVD basis",
            "query projection",
            "truth projection",
            "base score matrix",
        ],
        "basis_sha256": [sha256_array(cache.basis) for cache in all_caches],
        "rows_sha256": [cache.rows_sha256 for cache in all_caches],
        "semantics_changed": False,
    }
    if audit["contexts"] != audit["expected_contexts"]:
        raise AssertionError("unexpected fit-cache context count")
    return inner, pretest, final, audit


def complete_candidate_selection(
    *,
    train_queries: Mapping[str, np.ndarray],
    auxiliary_truth: Mapping[str, np.ndarray],
    correction_in_train_score: np.ndarray,
    splits: Mapping[str, np.ndarray],
    candidates: Sequence[Candidate],
    args: argparse.Namespace,
    permutation_replicate: int | None,
    permutation_plan: PermutationPlan | None = None,
    fit_caches: Mapping[tuple[str, str], FitCache] | None = None,
) -> dict[str, object]:
    """Repeat the complete grid search for real labels or one permutation null."""
    if (permutation_replicate is None) != (permutation_plan is None):
        raise ValueError("null replicate and frozen permutation plan must agree")
    if (
        permutation_plan is not None
        and permutation_plan.replicate != permutation_replicate
    ):
        raise ValueError("wrong frozen permutation plan supplied")
    cells_by_candidate: dict[str, list[dict[str, object]]] = {
        candidate.key: [] for candidate in candidates
    }
    cells: list[dict[str, object]] = []
    for outer in DEVELOPMENT_FOLDS:
        for held_teacher in TEACHER_DONORS:
            fit_teachers = {
                donor: auxiliary_truth[donor]
                for donor in TEACHER_DONORS
                if donor != held_teacher
            }
            cache = (
                fit_caches[(outer, held_teacher)]
                if fit_caches is not None
                else prepare_fit_cache(
                    train_queries[outer],
                    fit_teachers,
                    correction_in_train_score,
                    splits["fit"],
                    basis_dim=args.basis_dim,
                    basis_seed=stable_seed(
                        args.optimizer_seed,
                        "inner-basis",
                        outer,
                        held_teacher,
                    ),
                )
            )
            by_rank = {
                rank: [candidate for candidate in candidates if candidate.rank == rank]
                for rank in sorted({candidate.rank for candidate in candidates})
            }
            fitted: list[tuple[Candidate, Adapter]] = []
            for rank, rank_candidates in by_rank.items():
                models = fit_adapter_batch(
                    rank_candidates,
                    cache,
                    epochs=args.epochs,
                    learning_rate=args.learning_rate,
                    gradient_clip=args.gradient_clip,
                    seed=stable_seed(
                        args.optimizer_seed,
                        "inner-init",
                        outer,
                        held_teacher,
                        rank,
                    ),
                    label_permutations=(
                        None
                        if permutation_plan is None
                        else permutation_plan.labels["fit"]
                    ),
                )
                fitted.extend(zip(rank_candidates, models, strict=True))
            for candidate, model in fitted:
                evaluated = evaluate_adapter(
                    train_queries[outer],
                    auxiliary_truth[held_teacher],
                    correction_in_train_score,
                    model,
                    splits["tune"],
                    splits["tune"],
                )
                cell = {
                    "outer_donor": outer,
                    "held_teacher": held_teacher,
                    "fit_teachers": list(fit_teachers),
                    "candidate": candidate.key,
                    "labels_permuted": permutation_replicate is not None,
                    "permutation_replicate": permutation_replicate,
                    **summarize_evaluation(evaluated),
                    "fit_last_loss": model.fit_audit["last_loss"],
                    "residual_product_frobenius": model.fit_audit[
                        "residual_product_frobenius"
                    ],
                    "label_permutation_sha256": model.fit_audit[
                        "label_permutation_sha256"
                    ],
                }
                cells_by_candidate[candidate.key].append(cell)
                cells.append(cell)
        enforce_memory_limit(args.max_rss_gib, f"selection {outer}")
    selected_key, selection = choose_candidate(cells_by_candidate, candidates)
    by_key = {candidate.key: candidate for candidate in candidates}
    return {
        "selected": selected_key,
        "candidate": None if selected_key == IDENTITY else by_key[selected_key],
        "selection": selection,
        "cells": cells,
        "permutation_replicate": permutation_replicate,
        "permutation_plan_sha256": (
            None if permutation_plan is None else permutation_plan.plan_sha256
        ),
        "complete_grid_repeated": True,
        "candidate_invariant_fit_cache_reused": True,
        "same_rank_candidates_batched": True,
        "candidate_count_including_identity": len(candidates) + 1,
    }


def evaluate_selected_pretest(
    *,
    selected: Mapping[str, object],
    train_queries: Mapping[str, np.ndarray],
    auxiliary_truth: Mapping[str, np.ndarray],
    correction_in_train_score: np.ndarray,
    splits: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    permutation_replicate: int | None,
    permutation_plan: PermutationPlan | None = None,
    fit_caches: Mapping[str, FitCache] | None = None,
) -> dict[str, object]:
    if (permutation_replicate is None) != (permutation_plan is None):
        raise ValueError("null replicate and frozen permutation plan must agree")
    if (
        permutation_plan is not None
        and permutation_plan.replicate != permutation_replicate
    ):
        raise ValueError("wrong frozen permutation plan supplied")
    candidate = selected["candidate"]
    fit_tune = np.sort(np.concatenate([splits["fit"], splits["tune"]]))
    cells: list[dict[str, object]] = []
    fold_angles: dict[str, float] = {}
    for outer in DEVELOPMENT_FOLDS:
        model = None
        if candidate is not None:
            assert isinstance(candidate, Candidate)
            cache = (
                fit_caches[outer]
                if fit_caches is not None
                else prepare_fit_cache(
                    train_queries[outer],
                    auxiliary_truth,
                    correction_in_train_score,
                    fit_tune,
                    basis_dim=args.basis_dim,
                    basis_seed=stable_seed(
                        args.optimizer_seed,
                        "pretest-basis",
                        outer,
                    ),
                )
            )
            model = fit_adapter(
                train_queries[outer],
                auxiliary_truth,
                correction_in_train_score,
                candidate,
                fit_tune,
                basis_dim=args.basis_dim,
                epochs=args.epochs,
                learning_rate=args.learning_rate,
                gradient_clip=args.gradient_clip,
                seed=stable_seed(args.optimizer_seed, "pretest", outer),
                basis_seed=stable_seed(args.optimizer_seed, "pretest-basis", outer),
                label_permutations=(
                    None
                    if permutation_plan is None
                    else permutation_plan.labels["fit_tune"]
                ),
                prepared_cache=cache,
            )
        fold_cell_angles = []
        for teacher in TEACHER_DONORS:
            evaluated = evaluate_adapter(
                train_queries[outer],
                auxiliary_truth[teacher],
                correction_in_train_score,
                model,
                splits["test"],
                splits["test"],
            )
            fold_cell_angles.append(float(evaluated["median_angle_degrees"]))
            cells.append(
                {
                    "outer_donor": outer,
                    "teacher": teacher,
                    "selected": selected["selected"],
                    "labels_permuted": permutation_replicate is not None,
                    "permutation_replicate": permutation_replicate,
                    "label_permutation_sha256": (
                        {}
                        if model is None
                        else model.fit_audit["label_permutation_sha256"]
                    ),
                    **summarize_evaluation(evaluated),
                }
            )
        fold_angles[outer] = float(max(fold_cell_angles))
    return {
        "cells": cells,
        "macro_delta": float(np.mean([float(cell["delta"]) for cell in cells])),
        "fold_median_angles_degrees": fold_angles,
        "all_fold_angles_below_10_degrees": every_fold_angle_below_limit(fold_angles),
        "permutation_plan_sha256": (
            None if permutation_plan is None else permutation_plan.plan_sha256
        ),
    }


def save_adapter(path: Path, adapter: Adapter) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        basis=np.asarray(adapter.basis, dtype=np.float64),
        u=np.asarray(adapter.u, dtype=np.float64),
        v=np.asarray(adapter.v, dtype=np.float64),
    )
    return {
        "identity": False,
        "path": str(path),
        "sha256": sha256_file(path),
        "candidate": asdict(adapter.candidate),
        "basis_sha256": sha256_array(adapter.basis),
        "u_sha256": sha256_array(adapter.u),
        "v_sha256": sha256_array(adapter.v),
        "fit_audit": dict(adapter.fit_audit),
    }


def load_frozen_adapter(record: Mapping[str, object]) -> Adapter:
    path = Path(str(record["path"]))
    if sha256_file(path) != record["sha256"]:
        raise ValueError(f"{path}: frozen adapter file hash changed")
    with np.load(path, allow_pickle=False) as archive:
        basis = np.asarray(archive["basis"], dtype=np.float64)
        u = np.asarray(archive["u"], dtype=np.float64)
        v = np.asarray(archive["v"], dtype=np.float64)
    for values, key in (
        (basis, "basis_sha256"),
        (u, "u_sha256"),
        (v, "v_sha256"),
    ):
        if sha256_array(values) != record[key]:
            raise ValueError(f"{path}: frozen adapter {key} changed")
    return Adapter(
        candidate=Candidate(**record["candidate"]),
        basis=basis,
        u=u,
        v=v,
        fit_audit=dict(record["fit_audit"]),
    )


def fit_and_freeze_final_models(
    *,
    root: Path,
    label: str,
    selected: Mapping[str, object],
    train_queries: Mapping[str, np.ndarray],
    auxiliary_truth: Mapping[str, np.ndarray],
    correction_in_train_score: np.ndarray,
    args: argparse.Namespace,
    permutation_replicate: int | None,
    permutation_plan: PermutationPlan | None = None,
    fit_caches: Mapping[str, FitCache] | None = None,
) -> dict[str, object]:
    if (permutation_replicate is None) != (permutation_plan is None):
        raise ValueError("null replicate and frozen permutation plan must agree")
    if (
        permutation_plan is not None
        and permutation_plan.replicate != permutation_replicate
    ):
        raise ValueError("wrong frozen permutation plan supplied")
    candidate = selected["candidate"]
    records: dict[str, object] = {}
    rows = np.arange(next(iter(train_queries.values())).shape[0], dtype=np.int64)
    for fold in FINAL_FOLDS:
        if candidate is None:
            records[fold] = {
                "identity": True,
                "candidate": None,
                "path": None,
                "permutation_plan_sha256": (
                    None if permutation_plan is None else permutation_plan.plan_sha256
                ),
            }
            continue
        assert isinstance(candidate, Candidate)
        cache = (
            fit_caches[fold]
            if fit_caches is not None
            else prepare_fit_cache(
                train_queries[fold],
                auxiliary_truth,
                correction_in_train_score,
                rows,
                basis_dim=args.basis_dim,
                basis_seed=stable_seed(args.optimizer_seed, "final-basis", fold),
            )
        )
        model = fit_adapter(
            train_queries[fold],
            auxiliary_truth,
            correction_in_train_score,
            candidate,
            rows,
            basis_dim=args.basis_dim,
            epochs=args.epochs,
            learning_rate=args.learning_rate,
            gradient_clip=args.gradient_clip,
            seed=stable_seed(args.optimizer_seed, "final", fold),
            basis_seed=stable_seed(args.optimizer_seed, "final-basis", fold),
            label_permutations=(
                None if permutation_plan is None else permutation_plan.labels["all"]
            ),
            prepared_cache=cache,
        )
        records[fold] = {
            **save_adapter(root / label / f"{fold}.npz", model),
            "permutation_plan_sha256": (
                None if permutation_plan is None else permutation_plan.plan_sha256
            ),
        }
        del model
        enforce_memory_limit(args.max_rss_gib, f"freeze {label}/{fold}")
    return records


def load_kolf_truth(
    root: Path,
    panel: Sequence[str],
    panel_score_genes: Sequence[str],
    manifest: Mapping[str, object],
) -> dict[str, object]:
    targets = read_lines(root / "perts.txt")
    genes = read_gene_lines(root / "genes.txt")
    panel_pos = {target: index for index, target in enumerate(panel)}
    shared_targets = tuple(target for target in panel if target in set(targets))
    if tuple(targets) != shared_targets:
        raise ValueError("KOLF compact target order is not the VCC-panel subset order")
    query_rows = np.asarray(
        [panel_pos[target] for target in shared_targets], dtype=np.int64
    )
    gene_pos = {gene: index for index, gene in enumerate(genes)}
    score_genes = tuple(gene for gene in panel_score_genes if gene in gene_pos)
    final_axis = manifest.get("final_score_axis", {})
    if score_genes != tuple(panel_score_genes):
        raise ValueError("KOLF score axis differs from the frozen final score axis")
    if len(score_genes) != int(final_axis.get("gene_count", -1)):
        raise ValueError("KOLF score-gene count changed")
    if sha256_lines(score_genes) != final_axis.get("gene_order_sha256"):
        raise ValueError("KOLF score-gene order/hash changed")
    query_columns = np.arange(len(score_genes), dtype=np.int64)
    truth_columns = np.asarray([gene_pos[gene] for gene in score_genes], dtype=np.int64)
    delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape != (len(targets), len(genes)):
        raise ValueError("KOLF delta shape does not match its axes")
    raw = np.asarray(delta[:, truth_columns], dtype=np.float64)
    values = raw - raw.mean(axis=0, keepdims=True)
    if not np.isfinite(values).all():
        raise ValueError("KOLF truth contains non-finite values")
    return {
        "values": values,
        "targets": shared_targets,
        "query_rows": query_rows,
        "query_columns": query_columns,
        "score_genes": score_genes,
        "audit": {
            "targets": len(shared_targets),
            "score_genes": len(score_genes),
            "target_order_sha256": sha256_lines(shared_targets),
            "score_gene_order_sha256": sha256_lines(score_genes),
            "representation": (
                "raw decoded quantized z-score, centered per gene across all "
                "compact KOLF panel rows"
            ),
            "used_for_selection": False,
            "used_for_fitting": False,
        },
    }


def run_experiment(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
) -> dict[str, object]:
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite nonempty output directory {args.output_dir}"
        )
    official = read_csv_column(args.genes_csv, "gene_name")
    panel = read_csv_column(args.panel_csv, "target_gene")
    _, train_targets, outer_targets, axes = validate_truth_cache(
        args.truth_cache,
        official,
        panel,
        args.donor_root,
        manifest,
    )
    query_metadata = validate_query_cache(
        args.query_cache,
        axes,
        train_targets,
        outer_targets,
        panel,
        manifest,
    )
    kolf_protocol = manifest.get("kolf_confirmation", {})
    if not bool(kolf_protocol.get("supported", False)):
        raise ValueError(
            "KOLF confirmation is not supported by the frozen manifest; "
            "promotion must fail closed"
        )
    splits = deterministic_target_split(train_targets, seed=args.split_seed)
    candidates = build_candidates(
        args.ranks,
        args.betas,
        args.temperatures,
        args.regularizations,
    )
    auxiliary_truth = {
        donor: load_truth(args.truth_cache, donor, "train", axes.train_score_mask)
        for donor in TEACHER_DONORS
    }
    train_queries = {
        held: np.asarray(
            np.load(
                args.query_cache / held / "train_query.npy",
                mmap_mode="r",
                allow_pickle=False,
            ),
            dtype=np.float64,
        )
        for held in FINAL_FOLDS
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inner_fit_caches, pretest_fit_caches, final_fit_caches, fit_cache_audit = (
        prepare_all_fit_caches(
            train_queries=train_queries,
            auxiliary_truth=auxiliary_truth,
            correction_in_train_score=axes.correction_in_train_score,
            splits=splits,
            args=args,
        )
    )
    permutation_plans = [
        build_permutation_plan(
            auxiliary_truth,
            splits,
            replicate=replicate,
            seed=args.null_seed,
        )
        for replicate in range(args.null_replicates)
    ]
    permutation_plan_freeze = save_permutation_plans(
        args.output_dir / "permutation_plans.npz",
        permutation_plans,
    )
    enforce_memory_limit(args.max_rss_gib, "training inputs")

    real_selection = complete_candidate_selection(
        train_queries=train_queries,
        auxiliary_truth=auxiliary_truth,
        correction_in_train_score=axes.correction_in_train_score,
        splits=splits,
        candidates=candidates,
        args=args,
        permutation_replicate=None,
        permutation_plan=None,
        fit_caches=inner_fit_caches,
    )
    real_pretest = evaluate_selected_pretest(
        selected=real_selection,
        train_queries=train_queries,
        auxiliary_truth=auxiliary_truth,
        correction_in_train_score=axes.correction_in_train_score,
        splits=splits,
        args=args,
        permutation_replicate=None,
        permutation_plan=None,
        fit_caches=pretest_fit_caches,
    )

    # Every null replicate repeats the complete inner grid selection before
    # its selected null candidate is refit on fit+tune and tested.
    null_runs: list[dict[str, object]] = []
    for permutation_plan in permutation_plans:
        replicate = permutation_plan.replicate
        null_selection = complete_candidate_selection(
            train_queries=train_queries,
            auxiliary_truth=auxiliary_truth,
            correction_in_train_score=axes.correction_in_train_score,
            splits=splits,
            candidates=candidates,
            args=args,
            permutation_replicate=replicate,
            permutation_plan=permutation_plan,
            fit_caches=inner_fit_caches,
        )
        null_pretest = evaluate_selected_pretest(
            selected=null_selection,
            train_queries=train_queries,
            auxiliary_truth=auxiliary_truth,
            correction_in_train_score=axes.correction_in_train_score,
            splits=splits,
            args=args,
            permutation_replicate=replicate,
            permutation_plan=permutation_plan,
            fit_caches=pretest_fit_caches,
        )
        null_runs.append(
            {
                "replicate": replicate,
                "selection": null_selection,
                "pretest": null_pretest,
            }
        )
        enforce_memory_limit(args.max_rss_gib, f"complete null selection {replicate}")

    selected_key = str(real_selection["selected"])
    selected_candidate = real_selection["candidate"]
    pretest_cells = real_pretest["cells"]
    pretest_deltas = [float(cell["delta"]) for cell in pretest_cells]
    pretest_macro = float(real_pretest["macro_delta"])
    pretest_null = [float(run["pretest"]["macro_delta"]) for run in null_runs]
    pretest_null_p95 = float(np.quantile(pretest_null, 0.95))
    inner_candidate_cells = real_selection["selection"]["candidate_summaries"].get(
        selected_key, {}
    )
    inner_deltas = list(inner_candidate_cells.get("cell_deltas", [0.0]))
    pre_gate = {
        "selected_nonidentity": selected_candidate is not None,
        "all_inner_holdout_deltas_nonnegative": bool(
            selected_candidate is not None and min(map(float, inner_deltas)) >= 0.0
        ),
        "all_nonpanel_test_cells_nonnegative": min(pretest_deltas) >= 0.0,
        "nonpanel_test_macro_positive": pretest_macro > 0.0,
        "nonpanel_test_gain_above_complete_selection_null_p95": (
            selected_candidate is not None and pretest_macro > pretest_null_p95
        ),
        "angular_median_below_10_degrees_in_every_development_fold": bool(
            real_pretest["all_fold_angles_below_10_degrees"]
        ),
    }
    pre_gate["passed"] = all(pre_gate.values())

    def selection_for_json(value: Mapping[str, object]) -> dict[str, object]:
        candidate = value["candidate"]
        return {
            "selected": value["selected"],
            "candidate": None if candidate is None else asdict(candidate),
            "selection": value["selection"],
            "complete_grid_repeated": value["complete_grid_repeated"],
            "candidate_count_including_identity": value[
                "candidate_count_including_identity"
            ],
            "permutation_replicate": value["permutation_replicate"],
            "permutation_plan_sha256": value["permutation_plan_sha256"],
            "candidate_invariant_fit_cache_reused": value[
                "candidate_invariant_fit_cache_reused"
            ],
            "same_rank_candidates_batched": value["same_rank_candidates_batched"],
        }

    selection_freeze: dict[str, object] = {
        "selected": selected_key,
        "candidate": None if selected_candidate is None else asdict(selected_candidate),
        "candidate_grid": {
            "identity": True,
            "ranks": list(args.ranks),
            "betas": list(args.betas),
            "temperatures": list(args.temperatures),
            "regularizations": list(args.regularizations),
            "count_including_identity": len(candidates) + 1,
        },
        "target_split": {
            name: {
                "count": int(len(indices)),
                "targets_sha256": sha256_lines(
                    train_targets[index] for index in indices.tolist()
                ),
            }
            for name, indices in splits.items()
        },
        "real_selection": selection_for_json(real_selection),
        "complete_selection_nulls": [
            {
                "replicate": run["replicate"],
                "selection": selection_for_json(run["selection"]),
                "pretest_macro_delta": run["pretest"]["macro_delta"],
            }
            for run in null_runs
        ],
        "permutation_plan_freeze": permutation_plan_freeze,
        "fit_cache_audit": fit_cache_audit,
        "pretest_gate": pre_gate,
        "final_models_fitted": False,
        "final_model_files": {},
        "fold_cluster_files": {},
        "outer_truth_opened": False,
    }
    selection_path = args.output_dir / "selection_freeze.json"

    result: dict[str, object] = {
        "schema_version": RESULT_SCHEMA,
        "experiment": "080-listwise-metric-adapter",
        "status": "pre_gate_passed" if pre_gate["passed"] else "pre_gate_failed",
        "method": {
            "formula": "normalize(z + beta * B U V^T B^T z)",
            "objective": (
                "temperature-scaled listwise softmax with same-target "
                "cross-donor positive and all other targets as negatives"
            ),
            "basis": (
                "randomized right singular basis fitted only on the current "
                "non-panel fitting partition and permitted teacher donors"
            ),
            "rank_limit": 16,
            "regularization": "0.5 * lambda * ||U V^T||_F^2",
            "identity_candidate": True,
        },
        "protocol": {
            "outer_donors": list(OUTER_DONORS),
            "development_folds": list(DEVELOPMENT_FOLDS),
            "required_final_folds": list(FINAL_FOLDS),
            "teacher_donors": list(TEACHER_DONORS),
            "invalid_k562_outer_reason": (
                "PerturbMatch is a reanalysis of the K562 assay"
            ),
            "arm": "frozen experiment-064 arm D",
            "baseline_reproduction": (
                "saved panel queries use exp064's original all-mapped-coordinate "
                "normalization followed by VCC-coordinate blocking and must match "
                "an independent raw-bank reconstruction within 1e-6"
            ),
            "learned_feature_separation": (
                "all adapter bases, weights, selection metrics, nulls, and source "
                "clusters use retained-coordinate-safe artifacts: train_query for "
                "fitting/selection and panel-row cluster_query for clustering; both "
                "exclude both target panels during centering, norm calibration, "
                "and whitener restriction"
            ),
            "training_targets": len(train_targets),
            "outer_targets": len(outer_targets),
            "target_split_seed": args.split_seed,
            "fit_targets": int(len(splits["fit"])),
            "tune_targets": int(len(splits["tune"])),
            "test_targets": int(len(splits["test"])),
            "outer_truth_open_policy": (
                "verify and open only after the non-panel pre-gate, exact final "
                "real/null weights, and fold-specific source-query clusters are frozen"
            ),
            "tune_and_test_gallery": "complete 600-target non-panel truth gallery",
            "null_plan": (
                "one pre-generated per-replicate target-priority realization "
                "induces and freezes fit, fit+tune, and all-target permutations"
            ),
            "fit_cache": (
                "nine candidate-invariant normalization/basis/projection/base-score "
                "contexts reused without changing epochs or optimizer state"
            ),
            "candidate_batching": (
                "the 36 candidates at each rank share batched linear algebra "
                "while retaining separate parameters and Adam state"
            ),
            "candidate_selection": (
                "exact macro tune PDS maximization among eligible candidates; "
                "only exact floating-point ties use rank, beta, regularization, "
                "temperature, then key as a deterministic ordering"
            ),
            "source_clustering": (
                "separately staged panel-row query using retained-coordinate-only "
                "centering, norm calibration, and principal whitener restriction; "
                "all VCC and historical target coordinates remain excluded"
            ),
            "kolf_confirmation": "required and supported by pinned artifacts",
        },
        "leakage_audit": {
            "vcc_target_coordinates_removed": True,
            "training_target_coordinates_removed_only_from_nonpanel_protocol": True,
            "final_panel_score_excludes_exactly_vcc_target_coordinates": True,
            "historical_target_coordinates_retained_for_final_panel": True,
            "dense_whitener_principal_restriction_used_for_nonpanel_protocol": True,
            "exact_exp064_baseline_used_only_for_fixed_panel_inference": True,
            "panel_baseline_used_for_adapter_fit_or_selection": False,
            "learned_statistics_exclude_all_target_coordinates": True,
            "fit_tune_test_disjoint": True,
            "outer_donor_absent_from_query_and_operator": True,
            "inner_held_teacher_absent_from_adapter_fit": True,
            "vcc_panel_used_for_feature_or_hyperparameter_selection": False,
            "outer_truth_opened_before_selection": False,
            "outer_truth_opened_before_pre_gate": False,
            "final_weights_fitted_after_panel_truth_authentication": False,
            "k562_truth_used": False,
        },
        "selection_freeze_sha256": None,
        "permutation_plan_freeze": permutation_plan_freeze,
        "fit_cache_audit": fit_cache_audit,
        "selection": real_selection["selection"],
        "inner_cells": real_selection["cells"],
        "pretest": {
            "cells": pretest_cells,
            "macro_delta": pretest_macro,
            "null_replicates": len(pretest_null),
            "null_deltas": pretest_null,
            "null_p95": pretest_null_p95,
            "null_complete_grid_repetitions": len(null_runs),
            "fold_median_angles_degrees": real_pretest["fold_median_angles_degrees"],
            "gate": pre_gate,
        },
        "outer": {
            "evaluated": False,
            "reason": (
                "non-panel pre-gate failed; held VCC-panel truth was not opened"
                if not pre_gate["passed"]
                else "pending"
            ),
        },
        "query_cache_metadata_sha256": sha256_file(args.query_cache / "metadata.json"),
        "query_cache_summary": {
            "arm": query_metadata["arm"],
            "development_folds": query_metadata["development_folds"],
            "final_folds": query_metadata["final_folds"],
            "train_score_genes": query_metadata["train_score_genes"],
            "panel_score_genes": query_metadata["panel_score_genes"],
            "cluster_score_genes": query_metadata["cluster_score_genes"],
            "correction_genes": query_metadata["correction_genes"],
        },
        "peak_rss_gib": peak_rss_gib(),
    }

    if not pre_gate["passed"]:
        selection_path.write_text(json.dumps(selection_freeze, indent=2) + "\n")
        result["selection_freeze_sha256"] = sha256_file(selection_path)
        (args.output_dir / "results.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result

    # Fit and serialize every real/null final model before authenticating any
    # panel truth. Post-authentication code only reloads these exact bytes.
    frozen_root = args.output_dir / "frozen_models"
    real_models = fit_and_freeze_final_models(
        root=frozen_root,
        label="real",
        selected=real_selection,
        train_queries=train_queries,
        auxiliary_truth=auxiliary_truth,
        correction_in_train_score=axes.correction_in_train_score,
        args=args,
        permutation_replicate=None,
        permutation_plan=None,
        fit_caches=final_fit_caches,
    )
    null_model_records = []
    for run in null_runs:
        replicate = int(run["replicate"])
        records = fit_and_freeze_final_models(
            root=frozen_root,
            label=f"null-{replicate:04d}",
            selected=run["selection"],
            train_queries=train_queries,
            auxiliary_truth=auxiliary_truth,
            correction_in_train_score=axes.correction_in_train_score,
            args=args,
            permutation_replicate=replicate,
            permutation_plan=permutation_plans[replicate],
            fit_caches=final_fit_caches,
        )
        null_model_records.append(
            {
                "replicate": replicate,
                "selected": run["selection"]["selected"],
                "candidate": (
                    None
                    if run["selection"]["candidate"] is None
                    else asdict(run["selection"]["candidate"])
                ),
                "models": records,
            }
        )

    cluster_queries = {
        held: np.asarray(
            np.load(
                args.query_cache / held / "cluster_query.npy",
                mmap_mode="r",
                allow_pickle=False,
            ),
            dtype=np.float64,
        )
        for held in FINAL_FOLDS
    }
    panel_queries = {
        held: np.asarray(
            np.load(
                args.query_cache / held / "panel_query.npy",
                mmap_mode="r",
                allow_pickle=False,
            ),
            dtype=np.float64,
        )
        for held in FINAL_FOLDS
    }
    cluster_root = args.output_dir / "frozen_clusters"
    cluster_root.mkdir(parents=True, exist_ok=True)
    cluster_labels_full: dict[str, np.ndarray] = {}
    cluster_audits: dict[str, object] = {}
    cluster_records: dict[str, object] = {}
    for fold in FINAL_FOLDS:
        labels, audit = source_effect_clusters(
            [cluster_queries[fold]],
            reciprocal_k=args.cluster_reciprocal_k,
            cosine_threshold=args.cluster_cosine_threshold,
        )
        path = cluster_root / f"{fold}.npy"
        np.save(path, labels, allow_pickle=False)
        cluster_labels_full[fold] = labels
        cluster_audits[fold] = {
            **audit,
            "fold": fold,
            "source_donors": query_metadata["folds"][fold]["source_donors"],
            "held_donor_present_in_cluster_source": False,
        }
        cluster_records[fold] = {
            "path": str(path),
            "sha256": sha256_file(path),
            "labels_sha256": sha256_array(labels),
            "query_sha256": query_metadata["folds"][fold]["panels"]["cluster"][
                "query_sha256"
            ],
            "audit": cluster_audits[fold],
        }

    selection_freeze["final_models_fitted"] = True
    selection_freeze["final_model_files"] = {
        "real": real_models,
        "null": null_model_records,
    }
    selection_freeze["fold_cluster_files"] = cluster_records
    selection_path.write_text(json.dumps(selection_freeze, indent=2) + "\n")
    selection_hash = sha256_file(selection_path)
    result["selection_freeze_sha256"] = selection_hash

    # The only panel-truth authentication boundary is below. No fit_adapter
    # call is permitted after this point.
    verify_manifest_inputs(args, manifest, include_outer_truth=True)
    outer_truth_with_audit = {
        held: load_compact_panel_truth(
            args.donor_root,
            held,
            outer_targets,
            axes.panel_score_genes,
            chunk_rows=args.truth_chunk_rows,
        )
        for held in DEVELOPMENT_FOLDS
    }
    outer_truth = {
        held: values for held, (values, _audit) in outer_truth_with_audit.items()
    }
    outer_truth_audit = {
        held: audit for held, (_values, audit) in outer_truth_with_audit.items()
    }
    kolf = load_kolf_truth(
        args.kolf_root,
        panel,
        axes.panel_score_genes,
        manifest,
    )
    if len(kolf["targets"]) != int(kolf_protocol["expected_targets"]):
        raise ValueError("KOLF target coverage differs from the frozen protocol")
    fold_specs: dict[str, dict[str, object]] = {
        held: {
            "query": panel_queries[held],
            "truth": outer_truth[held],
            "query_rows": np.arange(len(outer_targets), dtype=np.int64),
            "positive_columns": np.arange(len(outer_targets), dtype=np.int64),
            "query_columns": np.arange(len(axes.panel_score_genes), dtype=np.int64),
            "truth_columns": np.arange(len(axes.panel_score_genes), dtype=np.int64),
            "targets": outer_targets,
            "clusters": cluster_labels_full[held],
        }
        for held in DEVELOPMENT_FOLDS
    }
    fold_specs[KOLF_FOLD] = {
        "query": panel_queries[KOLF_FOLD],
        "truth": kolf["values"],
        "query_rows": kolf["query_rows"],
        "positive_columns": np.arange(len(kolf["targets"]), dtype=np.int64),
        "query_columns": kolf["query_columns"],
        "truth_columns": np.arange(len(kolf["score_genes"]), dtype=np.int64),
        "targets": kolf["targets"],
        "clusters": cluster_labels_full[KOLF_FOLD][kolf["query_rows"]],
    }

    outer_fold_results: dict[str, object] = {}
    per_target_rows: list[dict[str, object]] = []
    outer_target_deltas: dict[str, np.ndarray] = {}
    outer_angles: dict[str, np.ndarray] = {}
    for held in FINAL_FOLDS:
        model_record = real_models[held]
        model = None if model_record["identity"] else load_frozen_adapter(model_record)
        spec = fold_specs[held]
        evaluated = evaluate_adapter(
            spec["query"],
            spec["truth"],
            axes.correction_in_panel_score,
            model,
            spec["query_rows"],
            spec["positive_columns"],
            query_score_columns=spec["query_columns"],
            truth_score_columns=spec["truth_columns"],
        )
        outer_target_deltas[held] = evaluated["per_target_delta"]
        outer_angles[held] = evaluated["angles_degrees"]
        outer_fold_results[held] = {
            **summarize_evaluation(evaluated),
            "frozen_model": model_record,
            "source_donors": query_metadata["folds"][held]["source_donors"],
            "target_cluster_audit": cluster_audits[held],
            "held_truth": (
                kolf["audit"] if held == KOLF_FOLD else outer_truth_audit[held]
            ),
        }
        for index, target in enumerate(spec["targets"]):
            per_target_rows.append(
                {
                    "outer_donor": held,
                    "target_gene": target,
                    "source_cluster": int(spec["clusters"][index]),
                    "baseline_pds": float(evaluated["per_target_baseline"][index]),
                    "adapted_pds": float(evaluated["per_target_adapted"][index]),
                    "delta": float(evaluated["per_target_delta"][index]),
                    "angle_degrees": float(evaluated["angles_degrees"][index]),
                }
            )

    fold_deltas = {
        held: float(outer_fold_results[held]["delta"]) for held in FINAL_FOLDS
    }
    macro_delta = primary_fold_mean(fold_deltas)
    bootstrap = fold_specific_cluster_bootstrap(
        {held: outer_target_deltas[held] for held in PRIMARY_FOLDS},
        {
            held: np.asarray(fold_specs[held]["clusters"], dtype=np.int64)
            for held in PRIMARY_FOLDS
        },
        replicates=args.bootstrap_replicates,
        seed=args.bootstrap_seed,
    )
    dependence = fold_specific_cluster_dependence(
        {held: outer_target_deltas[held] for held in PRIMARY_FOLDS},
        {
            held: np.asarray(fold_specs[held]["clusters"], dtype=np.int64)
            for held in PRIMARY_FOLDS
        },
    )

    outer_null: list[float] = []
    outer_null_details: list[dict[str, object]] = []
    for null_record in null_model_records:
        null_fold_deltas: dict[str, float] = {}
        for held in FINAL_FOLDS:
            model_record = null_record["models"][held]
            null_model = (
                None if model_record["identity"] else load_frozen_adapter(model_record)
            )
            spec = fold_specs[held]
            evaluated = evaluate_adapter(
                spec["query"],
                spec["truth"],
                axes.correction_in_panel_score,
                null_model,
                spec["query_rows"],
                spec["positive_columns"],
                query_score_columns=spec["query_columns"],
                truth_score_columns=spec["truth_columns"],
            )
            null_fold_deltas[held] = float(evaluated["delta"])
        null_macro = primary_fold_mean(null_fold_deltas)
        outer_null.append(null_macro)
        outer_null_details.append(
            {
                "replicate": null_record["replicate"],
                "selected": null_record["selected"],
                "fold_deltas": null_fold_deltas,
                "macro_delta": null_macro,
            }
        )
        enforce_memory_limit(args.max_rss_gib, f"outer null {null_record['replicate']}")

    null_p95 = float(np.quantile(outer_null, 0.95))
    fold_median_angles = {
        held: float(np.median(outer_angles[held])) for held in FINAL_FOLDS
    }
    input_nonzero = all(
        np.all(
            np.linalg.norm(
                np.asarray(fold_specs[held]["query"])[
                    np.ix_(
                        np.asarray(fold_specs[held]["query_rows"], dtype=np.int64),
                        np.asarray(fold_specs[held]["query_columns"], dtype=np.int64),
                    )
                ],
                axis=1,
            )
            > 1e-12
        )
        for held in FINAL_FOLDS
    )
    output_nonzero = all(
        np.all(
            np.linalg.norm(
                transform_query(
                    np.asarray(fold_specs[held]["query"]),
                    axes.correction_in_panel_score,
                    (
                        None
                        if real_models[held]["identity"]
                        else load_frozen_adapter(real_models[held])
                    ),
                )[
                    np.ix_(
                        np.asarray(fold_specs[held]["query_rows"], dtype=np.int64),
                        np.asarray(fold_specs[held]["query_columns"], dtype=np.int64),
                    )
                ],
                axis=1,
            )
            > 1e-12
        )
        for held in FINAL_FOLDS
    )
    promotion = {
        "predeclared": True,
        "nonpanel_pre_gate_passed": True,
        "selected_nonidentity": selected_candidate is not None,
        "adapter_rank_at_most_16": selected_candidate.rank <= 16,
        "all_primary_outer_folds_nonnegative": min(
            fold_deltas[held] for held in PRIMARY_FOLDS
        )
        >= 0.0,
        "kolf_independent_confirmation_included": KOLF_FOLD in fold_deltas,
        "kolf_independent_confirmation_nonnegative": (fold_deltas[KOLF_FOLD] >= 0.0),
        "outer_macro_raw_pds_delta_at_least_0_003": macro_delta >= 0.003,
        "cluster_bootstrap_lower_95_above_zero": bootstrap["lower_95"] > 0.0,
        "observed_gain_above_permutation_null_p95": macro_delta > null_p95,
        "angular_median_below_10_degrees_in_every_fold": every_fold_angle_below_limit(
            fold_median_angles
        ),
        "target_coverage_unchanged": input_nonzero and output_nonzero,
        "no_obvious_single_target_cluster_dependence": dependence[
            "passes_no_obvious_single_family_dependence"
        ],
    }
    promotion["passed"] = all(
        bool(value) for key, value in promotion.items() if key != "predeclared"
    )
    result["status"] = "completed"
    result["leakage_audit"]["outer_truth_opened_after_pre_gate"] = True
    result["outer"] = {
        "evaluated": True,
        "selection_freeze_sha256_before_truth_open": selection_hash,
        "folds": outer_fold_results,
        "primary_folds": list(PRIMARY_FOLDS),
        "confirmation_fold": KOLF_FOLD,
        "macro_baseline_pds": float(
            np.mean([outer_fold_results[held]["baseline"] for held in PRIMARY_FOLDS])
        ),
        "macro_adapted_pds": float(
            np.mean([outer_fold_results[held]["adapted"] for held in PRIMARY_FOLDS])
        ),
        "macro_delta": macro_delta,
        "fold_deltas": fold_deltas,
        "fold_median_angular_displacement_degrees": fold_median_angles,
        "target_clusters": cluster_audits,
        "cluster_bootstrap": bootstrap,
        "cluster_dependence": dependence,
        "kolf_truth": kolf["audit"],
        "permutation_null": {
            "method": (
                "repeat the complete candidate grid and selection under each "
                "response-norm-decile label permutation, freeze that null's "
                "selected final weights before panel truth, then evaluate"
            ),
            "replicates": len(outer_null),
            "macro_deltas": outer_null,
            "details": outer_null_details,
            "p95": null_p95,
        },
        "promotion_gate": promotion,
    }
    result["peak_rss_gib"] = peak_rss_gib()
    (args.output_dir / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    write_csv(args.output_dir / "per_target_pds.csv", per_target_rows)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-root", type=Path, required=True)
    parser.add_argument("--truth-cache", type=Path, required=True)
    parser.add_argument("--kolf-root", type=Path, required=True)
    parser.add_argument("--query-cache", type=Path, required=True)
    parser.add_argument("--artifact-manifest", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=HERE / "output")
    parser.add_argument("--ranks", nargs="+", type=int, default=[4, 8, 16])
    parser.add_argument(
        "--betas", nargs="+", type=float, default=[0.02, 0.05, 0.10, 0.20]
    )
    parser.add_argument(
        "--temperatures", nargs="+", type=float, default=[0.03, 0.07, 0.15]
    )
    parser.add_argument(
        "--regularizations",
        nargs="+",
        type=float,
        default=[1e-4, 1e-3, 1e-2],
    )
    parser.add_argument("--basis-dim", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--learning-rate", type=float, default=0.03)
    parser.add_argument("--gradient-clip", type=float, default=5.0)
    parser.add_argument("--split-seed", type=int, default=8002026)
    parser.add_argument("--optimizer-seed", type=int, default=8012026)
    parser.add_argument("--null-seed", type=int, default=8022026)
    parser.add_argument("--null-replicates", type=int, default=100)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--bootstrap-seed", type=int, default=8032026)
    parser.add_argument("--cluster-reciprocal-k", type=int, default=5)
    parser.add_argument("--cluster-cosine-threshold", type=float, default=0.80)
    parser.add_argument("--truth-chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=16.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    if args.basis_dim < max(args.ranks):
        raise SystemExit("--basis-dim must be at least the maximum adapter rank")
    if (
        args.epochs < 1
        or args.learning_rate <= 0
        or args.gradient_clip <= 0
        or args.truth_chunk_rows < 1
    ):
        raise SystemExit("optimizer controls must be positive")
    if args.null_replicates < 20:
        raise SystemExit("--null-replicates must be at least 20")
    if args.bootstrap_replicates < 100:
        raise SystemExit("--bootstrap-replicates must be at least 100")
    if args.max_rss_gib <= 0:
        raise SystemExit("--max-rss-gib must be positive")
    manifest = json.loads(args.artifact_manifest.read_text())
    if int(manifest.get("schema_version", -1)) != 2:
        raise SystemExit("unsupported artifact manifest")
    # Deliberately defer hashing/opening *_raw_outer.npy until the pre-gate.
    verify_manifest_inputs(args, manifest, include_outer_truth=False)
    result = run_experiment(args, manifest)
    print(json.dumps(result["pretest"]["gate"], indent=2), flush=True)
    if result["outer"]["evaluated"]:
        print(json.dumps(result["outer"]["promotion_gate"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
