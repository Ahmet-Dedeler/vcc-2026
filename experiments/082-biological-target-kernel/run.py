#!/usr/bin/env python3
"""Leakage-safe ESM2 target-kernel residual model around frozen arm D.

The model predicts only a low-rank residual:

    residual(d, p) = truth(d, p) - arm_D(excluding d, p)
    prediction(p)  = arm_D(p) + beta * residual_hat(ESM2(p))

All feature standardization, residual PCA, ridge fitting, and hyperparameter
selection use deterministic non-panel partitions. Held VCC-panel truth is
opened only after the real non-panel pre-gate and both repeated-selection null
gates pass.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
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
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STAGE_PATH = HERE / "stage_cache.py"
DEVELOPMENT_FOLDS = ("orion_hct116", "orion_hek293t")
KOLF_FOLD = "kolf2"
FINAL_FOLDS = (*DEVELOPMENT_FOLDS, KOLF_FOLD)
OUTER_DONORS = DEVELOPMENT_FOLDS
TEACHER_DONORS = ("rpe1", "nadig_jurkat")
IDENTITY = "identity"
CACHE_SCHEMA = 3
RESULT_SCHEMA = 1
TIE_WINDOW = 5e-4
ANGLE_LIMIT_DEGREES = 10.0
EXPECTED_EMBEDDING_ROWS = 19_790
EXPECTED_EMBEDDING_DIM = 5_120
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


def stable_seed(*parts: object) -> int:
    text = "|".join(map(str, parts))
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4], "little")


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values:
        raise ValueError(f"{path}: axis is empty")
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: axis contains duplicates")
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
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: {column!r} is empty or contains duplicates")
    return values


def row_normalize(values: np.ndarray) -> np.ndarray:
    values64 = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(values64, axis=1, keepdims=True)
    return np.divide(values64, np.maximum(norms, 1e-12))


def pds_per_target_from_scores(
    scores: np.ndarray,
    positive_columns: np.ndarray,
) -> np.ndarray:
    if scores.ndim != 2 or scores.shape[1] < 2:
        raise ValueError(f"PDS requires query x gallery scores, got {scores.shape}")
    positives = np.asarray(positive_columns, dtype=np.int64)
    if positives.shape != (scores.shape[0],):
        raise ValueError("positive columns do not align with query rows")
    if np.any(positives < 0) or np.any(positives >= scores.shape[1]):
        raise ValueError("positive gallery column is out of range")
    own = scores[np.arange(len(scores)), positives][:, None]
    ranks = (scores > own).sum(axis=1)
    ties = (scores == own).sum(axis=1) - 1
    return 1.0 - (ranks + 0.5 * ties) / (scores.shape[1] - 1)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@dataclass(frozen=True, order=True)
class Candidate:
    rank: int
    ridge: float
    beta: float

    @property
    def key(self) -> str:
        return f"esm2-ridge:rank={self.rank}:ridge={self.ridge:g}:beta={self.beta:g}"


@dataclass(frozen=True)
class EmbeddingTable:
    targets: tuple[str, ...]
    values: np.ndarray
    available: np.ndarray
    feature_names: tuple[str, ...]
    index_column: str
    audit: Mapping[str, object]


@dataclass(frozen=True)
class FeatureAxes:
    master_genes: tuple[str, ...]
    dev_score_mask: np.ndarray
    dev_score_genes: tuple[str, ...]
    panel_score_genes: tuple[str, ...]
    panel_score_official: np.ndarray
    dev_to_panel_positions: np.ndarray
    donor_axis_audits: Mapping[str, object]


@dataclass(frozen=True)
class ResidualBasis:
    residual_mean: np.ndarray
    components: np.ndarray
    coefficients: np.ndarray
    fit_rows: np.ndarray
    audit: Mapping[str, object]


@dataclass(frozen=True)
class ResidualFamily:
    basis: ResidualBasis
    embedding_mean: np.ndarray
    embedding_scale: np.ndarray
    fit_embedding: np.ndarray
    dual_by_ridge: Mapping[float, np.ndarray]
    audit: Mapping[str, object]


@dataclass(frozen=True)
class FrozenPredictor:
    candidate: Candidate
    embedding_mean: np.ndarray
    embedding_scale: np.ndarray
    primal_weight: np.ndarray
    residual_mean: np.ndarray
    components: np.ndarray
    audit: Mapping[str, object]


@dataclass(frozen=True)
class SelectionResult:
    selected_key: str
    selected: Candidate | None
    audit: Mapping[str, object]
    cells: tuple[Mapping[str, object], ...]


def build_candidates(
    ranks: Sequence[int],
    ridges: Sequence[float],
    betas: Sequence[float],
) -> list[Candidate]:
    candidates = sorted(
        {
            Candidate(int(rank), float(ridge), float(beta))
            for rank in ranks
            for ridge in ridges
            for beta in betas
        }
    )
    if not candidates:
        raise ValueError("candidate grid is empty")
    if any(candidate.rank <= 0 or candidate.rank > 64 for candidate in candidates):
        raise ValueError("residual PCA rank must lie in [1,64]")
    if any(candidate.ridge <= 0 for candidate in candidates):
        raise ValueError("ridge strengths must be positive")
    if any(not (0.0 < candidate.beta <= 0.20) for candidate in candidates):
        raise ValueError("non-identity beta must lie in (0,0.20]")
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


def verify_embedding_artifact(
    path: Path,
    manifest: Mapping[str, object],
) -> dict[str, object]:
    expected = manifest["embedding"]
    actual_size = path.stat().st_size
    if actual_size != int(expected["size_bytes"]):
        raise ValueError(
            f"embedding size {actual_size} != expected {expected['size_bytes']}"
        )
    actual_sha256 = sha256_file(path)
    if actual_sha256 != str(expected["sha256"]):
        raise ValueError(
            f"embedding SHA-256 {actual_sha256} != expected {expected['sha256']}"
        )
    return {
        "path": str(path.resolve()),
        "size_bytes": actual_size,
        "sha256": actual_sha256,
        "repository": expected["repository"],
        "revision": expected["revision"],
        "filename": expected["filename"],
    }


def _parquet_schema(path: Path):
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as error:
        raise RuntimeError("pyarrow is required to read ESM2 parquet") from error
    parquet = pq.ParquetFile(path)
    schema = parquet.schema_arrow
    numeric = sorted(
        (
            name
            for name, field in zip(schema.names, schema, strict=True)
            if pa.types.is_floating(field.type) and name.isdigit()
        ),
        key=int,
    )
    other = [name for name in schema.names if name not in set(numeric)]
    if len(other) != 1:
        raise ValueError(
            f"expected one target-index column, found {other[:5]} ({len(other)} total)"
        )
    if [int(name) for name in numeric] != list(range(len(numeric))):
        raise ValueError("numeric embedding columns are not contiguous from zero")
    return parquet, tuple(numeric), other[0]


def load_target_embeddings(
    path: Path,
    targets: Sequence[str],
    *,
    expected_rows: int | None = EXPECTED_EMBEDDING_ROWS,
    expected_dim: int | None = EXPECTED_EMBEDDING_DIM,
    column_chunk: int = 128,
) -> EmbeddingTable:
    """Read only requested rows while scanning numeric Parquet columns in chunks."""
    parquet, feature_names, index_column = _parquet_schema(path)
    if expected_rows is not None and parquet.metadata.num_rows != expected_rows:
        raise ValueError(
            f"embedding rows {parquet.metadata.num_rows} != expected {expected_rows}"
        )
    if expected_dim is not None and len(feature_names) != expected_dim:
        raise ValueError(
            f"embedding width {len(feature_names)} != expected {expected_dim}"
        )
    index_table = parquet.read(columns=[index_column], use_threads=False)
    identifiers = tuple(map(str, index_table.column(0).to_pylist()))
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("embedding target index contains duplicates")
    positions = {target: index for index, target in enumerate(identifiers)}
    requested_positions = np.asarray(
        [positions.get(target, -1) for target in targets], dtype=np.int64
    )
    available = requested_positions >= 0
    values = np.full((len(targets), len(feature_names)), np.nan, dtype=np.float32)
    destination = np.flatnonzero(available)
    source_rows = requested_positions[available]
    for start in range(0, len(feature_names), column_chunk):
        names = feature_names[start : start + column_chunk]
        table = parquet.read(columns=list(names), use_threads=True)
        block = np.column_stack(
            [
                np.asarray(table.column(index).to_numpy(zero_copy_only=False))
                for index in range(len(names))
            ]
        )
        values[np.ix_(destination, np.arange(start, start + len(names)))] = block[
            source_rows
        ]
    if available.any() and not np.isfinite(values[available]).all():
        raise ValueError("available embedding rows contain non-finite values")
    if np.any(np.isfinite(values[~available])):
        raise AssertionError("missing embedding rows were not preserved as NaN")
    return EmbeddingTable(
        targets=tuple(targets),
        values=values,
        available=available,
        feature_names=feature_names,
        index_column=index_column,
        audit={
            "parquet_rows": int(parquet.metadata.num_rows),
            "parquet_row_groups": int(parquet.metadata.num_row_groups),
            "embedding_dimensions": len(feature_names),
            "index_column": index_column,
            "requested_targets": len(targets),
            "available_targets": int(np.sum(available)),
            "missing_targets": [
                target
                for target, present in zip(targets, available, strict=True)
                if not present
            ],
            "target_order_sha256": sha256_lines(targets),
            "feature_order_sha256": sha256_lines(feature_names),
            "available_mask_sha256": sha256_array(available),
        },
    )


def _stable_component_signs(components: np.ndarray) -> np.ndarray:
    result = np.asarray(components, dtype=np.float64).copy()
    for row in result:
        pivot = int(np.argmax(np.abs(row)))
        if row[pivot] < 0:
            row *= -1
    return result


def residual_pca(
    residuals: np.ndarray,
    max_rank: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(residuals, dtype=np.float64)
    if values.ndim != 2 or len(values) < max_rank + 1:
        raise ValueError(
            f"residual PCA rank {max_rank} needs at least {max_rank + 1} rows"
        )
    mean = values.mean(axis=0)
    centered = values - mean
    gram = centered @ centered.T
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = np.maximum(eigenvalues[order], 0.0)
    eigenvectors = eigenvectors[:, order]
    tolerance = max(float(eigenvalues[0]), 1.0) * 1e-12
    usable = np.flatnonzero(eigenvalues > tolerance)
    if len(usable) < max_rank:
        raise ValueError(
            f"residual matrix has only {len(usable)} stable components; need {max_rank}"
        )
    selected = usable[:max_rank]
    components = (eigenvectors[:, selected].T @ centered) / np.sqrt(
        eigenvalues[selected]
    )[:, None]
    components = _stable_component_signs(components)
    coefficients = centered @ components.T
    np.testing.assert_allclose(
        components @ components.T,
        np.eye(max_rank),
        atol=2e-6,
        rtol=2e-6,
    )
    explained = eigenvalues[selected] / max(float(np.sum(eigenvalues)), 1e-18)
    return mean, components, coefficients, explained


def fit_residual_basis(
    query: np.ndarray,
    truths: Mapping[str, np.ndarray],
    available: np.ndarray,
    rows: np.ndarray,
    *,
    max_rank: int,
) -> ResidualBasis:
    if not truths:
        raise ValueError("residual fit requires at least one teacher donor")
    rows = np.asarray(rows, dtype=np.int64)
    fit_rows = rows[np.asarray(available[rows], dtype=bool)]
    if len(fit_rows) < max_rank + 1:
        raise ValueError(
            f"only {len(fit_rows)} embedded fit targets for rank {max_rank}"
        )
    base = np.asarray(query[fit_rows], dtype=np.float64)
    residuals = np.mean(
        np.stack(
            [
                np.asarray(truth[fit_rows], dtype=np.float64) - base
                for truth in truths.values()
            ],
            axis=0,
        ),
        axis=0,
    )
    residual_mean, components, coefficients, explained = residual_pca(
        residuals, max_rank
    )
    return ResidualBasis(
        residual_mean=residual_mean,
        components=components,
        coefficients=coefficients,
        fit_rows=fit_rows,
        audit={
            "requested_fit_rows": int(len(rows)),
            "embedded_fit_rows": int(len(fit_rows)),
            "teachers": list(truths),
            "max_residual_rank": int(max_rank),
            "explained_variance_sum": float(np.sum(explained)),
            "fit_rows_sha256": sha256_array(fit_rows),
        },
    )


def fit_residual_family(
    embeddings: np.ndarray,
    available: np.ndarray,
    basis: ResidualBasis,
    *,
    ridges: Sequence[float],
) -> ResidualFamily:
    fit_rows = basis.fit_rows
    raw_x = np.asarray(embeddings[fit_rows], dtype=np.float64)
    embedding_mean = raw_x.mean(axis=0)
    embedding_scale = raw_x.std(axis=0)
    constant = embedding_scale < 1e-8
    embedding_scale[constant] = 1.0
    fit_embedding = (raw_x - embedding_mean) / embedding_scale
    fit_embedding /= np.sqrt(fit_embedding.shape[1])
    kernel = fit_embedding @ fit_embedding.T
    identity = np.eye(len(fit_embedding), dtype=np.float64)
    dual_by_ridge = {
        float(ridge): np.linalg.solve(
            kernel + float(ridge) * identity,
            basis.coefficients,
        )
        for ridge in sorted(set(map(float, ridges)))
    }
    return ResidualFamily(
        basis=basis,
        embedding_mean=embedding_mean,
        embedding_scale=embedding_scale,
        fit_embedding=fit_embedding,
        dual_by_ridge=dual_by_ridge,
        audit={
            **basis.audit,
            "embedding_dimensions": int(raw_x.shape[1]),
            "constant_embedding_features": int(np.sum(constant)),
            "ridge_strengths": sorted(set(map(float, ridges))),
        },
    )


def predict_coefficients(
    family: ResidualFamily,
    embeddings: np.ndarray,
    available: np.ndarray,
    ridge: float,
) -> np.ndarray:
    values = np.zeros(
        (len(embeddings), family.basis.components.shape[0]), dtype=np.float64
    )
    rows = np.flatnonzero(available)
    if len(rows):
        standardized = (
            np.asarray(embeddings[rows], dtype=np.float64) - family.embedding_mean
        ) / family.embedding_scale
        standardized /= np.sqrt(standardized.shape[1])
        kernel = standardized @ family.fit_embedding.T
        values[rows] = kernel @ family.dual_by_ridge[float(ridge)]
    return values


def explicit_candidate_query(
    query: np.ndarray,
    embeddings: np.ndarray,
    available: np.ndarray,
    family: ResidualFamily | FrozenPredictor,
    candidate: Candidate,
    residual_positions: np.ndarray | None = None,
) -> np.ndarray:
    if isinstance(family, FrozenPredictor):
        coefficients = frozen_predict_coefficients(family, embeddings, available)
        mean = family.residual_mean
        components = family.components
    else:
        coefficients = predict_coefficients(
            family, embeddings, available, candidate.ridge
        )[:, : candidate.rank]
        mean = family.basis.residual_mean
        components = family.basis.components[: candidate.rank]
    residual = mean + coefficients @ components
    residual[~available] = 0.0
    output = np.asarray(query, dtype=np.float64).copy()
    positions = (
        np.arange(output.shape[1], dtype=np.int64)
        if residual_positions is None
        else np.asarray(residual_positions, dtype=np.int64)
    )
    if len(positions) != residual.shape[1]:
        raise ValueError("residual positions do not match residual width")
    output[:, positions] += candidate.beta * residual
    return output


def freeze_predictor(
    family: ResidualFamily,
    candidate: Candidate,
) -> FrozenPredictor:
    dual = family.dual_by_ridge[float(candidate.ridge)][:, : candidate.rank]
    primal = family.fit_embedding.T @ dual
    predictor = FrozenPredictor(
        candidate=candidate,
        embedding_mean=np.asarray(family.embedding_mean, dtype=np.float32),
        embedding_scale=np.asarray(family.embedding_scale, dtype=np.float32),
        primal_weight=np.asarray(primal, dtype=np.float32),
        residual_mean=np.asarray(family.basis.residual_mean, dtype=np.float32),
        components=np.asarray(
            family.basis.components[: candidate.rank], dtype=np.float32
        ),
        audit={
            **family.audit,
            "candidate": asdict(candidate),
            "representation": (
                "fit-only embedding mean/scale, primal ridge weights, "
                "residual mean, and residual PCA components"
            ),
        },
    )
    return predictor


def frozen_predict_coefficients(
    predictor: FrozenPredictor,
    embeddings: np.ndarray,
    available: np.ndarray,
) -> np.ndarray:
    values = np.zeros((len(embeddings), predictor.candidate.rank), dtype=np.float64)
    rows = np.flatnonzero(available)
    if len(rows):
        standardized = (
            np.asarray(embeddings[rows], dtype=np.float64) - predictor.embedding_mean
        ) / predictor.embedding_scale
        standardized /= np.sqrt(standardized.shape[1])
        values[rows] = standardized @ predictor.primal_weight
    return values


def predictor_sha256(predictor: FrozenPredictor | None) -> str:
    if predictor is None:
        return hashlib.sha256(IDENTITY.encode()).hexdigest()
    digest = hashlib.sha256()
    digest.update(
        json.dumps(
            asdict(predictor.candidate), sort_keys=True, separators=(",", ":")
        ).encode()
    )
    for array in (
        predictor.embedding_mean,
        predictor.embedding_scale,
        predictor.primal_weight,
        predictor.residual_mean,
        predictor.components,
    ):
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(np.asarray(contiguous.shape, dtype="<i8").tobytes())
        digest.update(contiguous.view(np.uint8))
    return digest.hexdigest()


def candidate_scores(
    query: np.ndarray,
    truth: np.ndarray,
    embeddings: np.ndarray,
    available: np.ndarray,
    family: ResidualFamily | FrozenPredictor | None,
    candidate: Candidate | None,
    residual_positions: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q = np.asarray(query, dtype=np.float64)
    truth_unit = row_normalize(truth)
    q_norm_sq = np.sum(q * q, axis=1)
    base_numerator = q @ truth_unit.T
    base_scores = base_numerator / np.sqrt(np.maximum(q_norm_sq, 1e-18))[:, None]
    if candidate is None:
        return (
            base_scores,
            np.sqrt(np.maximum(q_norm_sq, 0.0)),
            np.zeros(len(q), dtype=np.float64),
        )
    if family is None:
        raise ValueError("non-identity candidate requires a fitted residual family")
    if isinstance(family, FrozenPredictor):
        if family.candidate != candidate:
            raise ValueError("frozen predictor/candidate mismatch")
        coefficients = frozen_predict_coefficients(family, embeddings, available)
        components = np.asarray(family.components, dtype=np.float64)
        mean = np.asarray(family.residual_mean, dtype=np.float64)
    else:
        coefficients = predict_coefficients(
            family, embeddings, available, candidate.ridge
        )[:, : candidate.rank]
        components = family.basis.components[: candidate.rank]
        mean = family.basis.residual_mean
    positions = (
        np.arange(q.shape[1], dtype=np.int64)
        if residual_positions is None
        else np.asarray(residual_positions, dtype=np.int64)
    )
    if len(positions) != len(mean):
        raise ValueError("residual positions do not match fitted residual axis")
    q_residual_axis = q[:, positions]
    truth_residual_axis = truth_unit[:, positions]
    residual_dot_truth = mean @ truth_residual_axis.T + coefficients @ (
        components @ truth_residual_axis.T
    )
    q_dot_residual = q_residual_axis @ mean + np.sum(
        (q_residual_axis @ components.T) * coefficients, axis=1
    )
    mean_component = components @ mean
    residual_norm_sq = (
        float(mean @ mean)
        + 2.0 * (coefficients @ mean_component)
        + np.sum(coefficients * coefficients, axis=1)
    )
    inactive = ~np.asarray(available, dtype=bool)
    residual_dot_truth[inactive] = 0.0
    q_dot_residual[inactive] = 0.0
    residual_norm_sq[inactive] = 0.0
    numerator = base_numerator + candidate.beta * residual_dot_truth
    norm_sq = (
        q_norm_sq
        + 2.0 * candidate.beta * q_dot_residual
        + candidate.beta**2 * residual_norm_sq
    )
    norms = np.sqrt(np.maximum(norm_sq, 1e-18))
    cosine_to_baseline = np.clip(
        (q_norm_sq + candidate.beta * q_dot_residual)
        / np.maximum(np.sqrt(q_norm_sq) * norms, 1e-18),
        -1.0,
        1.0,
    )
    angles = np.degrees(np.arccos(cosine_to_baseline))
    return numerator / norms[:, None], norms, angles


def evaluate_candidate(
    query: np.ndarray,
    truth_gallery: np.ndarray,
    embeddings: np.ndarray,
    available: np.ndarray,
    family: ResidualFamily | FrozenPredictor | None,
    candidate: Candidate | None,
    query_rows: np.ndarray,
    positive_columns: np.ndarray,
    residual_positions: np.ndarray | None = None,
) -> dict[str, object]:
    rows = np.asarray(query_rows, dtype=np.int64)
    positives = np.asarray(positive_columns, dtype=np.int64)
    local_query = np.asarray(query[rows], dtype=np.float64)
    local_embedding = np.asarray(embeddings[rows])
    local_available = np.asarray(available[rows], dtype=bool)
    baseline_scores, baseline_norms, _ = candidate_scores(
        local_query,
        truth_gallery,
        local_embedding,
        local_available,
        None,
        None,
        residual_positions,
    )
    adapted_scores, adapted_norms, angles = candidate_scores(
        local_query,
        truth_gallery,
        local_embedding,
        local_available,
        family,
        candidate,
        residual_positions,
    )
    baseline = pds_per_target_from_scores(baseline_scores, positives)
    adapted = pds_per_target_from_scores(adapted_scores, positives)
    delta = adapted - baseline
    return {
        "baseline": float(np.mean(baseline)),
        "adapted": float(np.mean(adapted)),
        "delta": float(np.mean(delta)),
        "per_target_baseline": baseline,
        "per_target_adapted": adapted,
        "per_target_delta": delta,
        "baseline_norms": baseline_norms,
        "adapted_norms": adapted_norms,
        "angles_degrees": angles,
        "median_angle_degrees": float(np.median(angles)),
        "max_angle_degrees": float(np.max(angles)),
        "embedding_available": local_available,
        "gallery_size": int(len(truth_gallery)),
    }


def summarize_evaluation(result: Mapping[str, object]) -> dict[str, object]:
    return {
        "baseline": float(result["baseline"]),
        "adapted": float(result["adapted"]),
        "delta": float(result["delta"]),
        "targets": int(len(result["per_target_delta"])),
        "gallery_size": int(result["gallery_size"]),
        "embedding_available": int(np.sum(result["embedding_available"])),
        "median_angle_degrees": float(result["median_angle_degrees"]),
        "max_angle_degrees": float(result["max_angle_degrees"]),
    }


def evaluate_candidate_grid(
    query: np.ndarray,
    truth_gallery: np.ndarray,
    embeddings: np.ndarray,
    available: np.ndarray,
    family: ResidualFamily,
    candidates: Sequence[Candidate],
    query_rows: np.ndarray,
    positive_columns: np.ndarray,
) -> dict[str, dict[str, object]]:
    """Evaluate a whole rank/ridge/beta grid with one complete-gallery baseline."""
    rows = np.asarray(query_rows, dtype=np.int64)
    positives = np.asarray(positive_columns, dtype=np.int64)
    q = np.asarray(query[rows], dtype=np.float64)
    truth_unit = row_normalize(truth_gallery)
    local_embedding = np.asarray(embeddings[rows])
    local_available = np.asarray(available[rows], dtype=bool)
    q_norm_sq = np.sum(q * q, axis=1)
    base_numerator = q @ truth_unit.T
    base_scores = base_numerator / np.sqrt(np.maximum(q_norm_sq, 1e-18))[:, None]
    baseline = pds_per_target_from_scores(base_scores, positives)
    output: dict[str, dict[str, object]] = {}
    by_ridge: dict[float, list[Candidate]] = {}
    for candidate in candidates:
        by_ridge.setdefault(candidate.ridge, []).append(candidate)
    for ridge, ridge_candidates in by_ridge.items():
        all_coefficients = predict_coefficients(
            family, local_embedding, local_available, ridge
        )
        for rank in sorted({candidate.rank for candidate in ridge_candidates}):
            coefficients = all_coefficients[:, :rank]
            components = family.basis.components[:rank]
            mean = family.basis.residual_mean
            residual_dot_truth = mean @ truth_unit.T + coefficients @ (
                components @ truth_unit.T
            )
            q_dot_residual = q @ mean + np.sum(
                (q @ components.T) * coefficients, axis=1
            )
            mean_component = components @ mean
            residual_norm_sq = (
                float(mean @ mean)
                + 2.0 * (coefficients @ mean_component)
                + np.sum(coefficients * coefficients, axis=1)
            )
            inactive = ~local_available
            residual_dot_truth[inactive] = 0.0
            q_dot_residual[inactive] = 0.0
            residual_norm_sq[inactive] = 0.0
            for candidate in (
                value for value in ridge_candidates if value.rank == rank
            ):
                numerator = base_numerator + candidate.beta * residual_dot_truth
                norm_sq = (
                    q_norm_sq
                    + 2.0 * candidate.beta * q_dot_residual
                    + candidate.beta**2 * residual_norm_sq
                )
                norms = np.sqrt(np.maximum(norm_sq, 1e-18))
                scores = numerator / norms[:, None]
                adapted = pds_per_target_from_scores(scores, positives)
                angles = np.degrees(
                    np.arccos(
                        np.clip(
                            (q_norm_sq + candidate.beta * q_dot_residual)
                            / np.maximum(np.sqrt(q_norm_sq) * norms, 1e-18),
                            -1.0,
                            1.0,
                        )
                    )
                )
                output[candidate.key] = {
                    "baseline": float(np.mean(baseline)),
                    "adapted": float(np.mean(adapted)),
                    "delta": float(np.mean(adapted - baseline)),
                    "per_target_baseline": baseline,
                    "per_target_adapted": adapted,
                    "per_target_delta": adapted - baseline,
                    "baseline_norms": np.sqrt(np.maximum(q_norm_sq, 0.0)),
                    "adapted_norms": norms,
                    "angles_degrees": angles,
                    "median_angle_degrees": float(np.median(angles)),
                    "max_angle_degrees": float(np.max(angles)),
                    "embedding_available": local_available,
                    "gallery_size": int(len(truth_gallery)),
                }
    if set(output) != {candidate.key for candidate in candidates}:
        raise AssertionError("candidate grid evaluation is incomplete")
    return output


def choose_candidate(
    cells_by_candidate: Mapping[str, Sequence[Mapping[str, object]]],
    candidates: Sequence[Candidate],
) -> tuple[str, dict[str, object]]:
    by_key = {candidate.key: candidate for candidate in candidates}
    summaries: dict[str, object] = {}
    eligible: list[str] = []
    for candidate in candidates:
        deltas = [float(cell["delta"]) for cell in cells_by_candidate[candidate.key]]
        summaries[candidate.key] = {
            "cell_deltas": deltas,
            "cell_median_angles_degrees": [
                float(cell["median_angle_degrees"])
                for cell in cells_by_candidate[candidate.key]
            ],
            "minimum_delta": float(min(deltas)),
            "macro_delta": float(np.mean(deltas)),
            "maximum_cell_median_angle_degrees": float(
                max(
                    float(cell["median_angle_degrees"])
                    for cell in cells_by_candidate[candidate.key]
                )
            ),
        }
        if (
            min(deltas) >= 0.0
            and np.mean(deltas) > 0.0
            and summaries[candidate.key]["maximum_cell_median_angle_degrees"]
            < ANGLE_LIMIT_DEGREES
        ):
            eligible.append(candidate.key)
    if not eligible:
        return IDENTITY, {
            "selected": IDENTITY,
            "eligible_count": 0,
            "candidate_summaries": summaries,
            "rule": (
                "identity unless every inner cell is nonnegative and every "
                "cell median angular displacement is below 10 degrees"
            ),
        }
    best = max(float(summaries[key]["macro_delta"]) for key in eligible)
    near = [
        key
        for key in eligible
        if float(summaries[key]["macro_delta"]) >= best - TIE_WINDOW
    ]
    selected = min(
        near,
        key=lambda key: (
            by_key[key].rank,
            by_key[key].beta,
            -by_key[key].ridge,
            key,
        ),
    )
    return selected, {
        "selected": selected,
        "eligible_count": len(eligible),
        "best_macro_delta": best,
        "tie_window": TIE_WINDOW,
        "candidate_summaries": summaries,
        "rule": (
            "maximize macro tune PDS among candidates nonnegative in every "
            "outer-role x held-teacher cell and below 10 degrees in every "
            "cell; within 5e-4 prefer lower rank, lower beta, then stronger ridge"
        ),
    }


def build_inner_residual_bases(
    train_queries: Mapping[str, np.ndarray],
    teacher_truth: Mapping[str, np.ndarray],
    available: np.ndarray,
    fit_rows: np.ndarray,
    *,
    max_rank: int,
) -> dict[tuple[str, str], ResidualBasis]:
    result: dict[tuple[str, str], ResidualBasis] = {}
    for outer in OUTER_DONORS:
        for held_teacher in TEACHER_DONORS:
            fit_teachers = {
                donor: teacher_truth[donor]
                for donor in TEACHER_DONORS
                if donor != held_teacher
            }
            result[(outer, held_teacher)] = fit_residual_basis(
                train_queries[outer],
                fit_teachers,
                available,
                fit_rows,
                max_rank=max_rank,
            )
    return result


def run_selection(
    train_queries: Mapping[str, np.ndarray],
    teacher_truth: Mapping[str, np.ndarray],
    embeddings: np.ndarray,
    available: np.ndarray,
    splits: Mapping[str, np.ndarray],
    candidates: Sequence[Candidate],
    *,
    max_rank: int,
    ridges: Sequence[float],
    max_rss_gib: float,
    residual_bases: Mapping[tuple[str, str], ResidualBasis] | None = None,
) -> SelectionResult:
    cells_by_candidate: dict[str, list[dict[str, object]]] = {
        candidate.key: [] for candidate in candidates
    }
    all_cells: list[dict[str, object]] = []
    for outer in OUTER_DONORS:
        for held_teacher in TEACHER_DONORS:
            fit_teachers = tuple(
                donor for donor in TEACHER_DONORS if donor != held_teacher
            )
            basis = (
                residual_bases[(outer, held_teacher)]
                if residual_bases is not None
                else fit_residual_basis(
                    train_queries[outer],
                    {donor: teacher_truth[donor] for donor in fit_teachers},
                    available,
                    splits["fit"],
                    max_rank=max_rank,
                )
            )
            family = fit_residual_family(
                embeddings,
                available,
                basis,
                ridges=ridges,
            )
            grid = evaluate_candidate_grid(
                train_queries[outer],
                teacher_truth[held_teacher],
                embeddings,
                available,
                family,
                candidates,
                splits["tune"],
                splits["tune"],
            )
            for candidate in candidates:
                evaluated = grid[candidate.key]
                cell = {
                    "outer_donor": outer,
                    "held_teacher": held_teacher,
                    "fit_teachers": list(fit_teachers),
                    "candidate": candidate.key,
                    **summarize_evaluation(evaluated),
                }
                cells_by_candidate[candidate.key].append(cell)
                all_cells.append(cell)
            enforce_memory_limit(max_rss_gib, f"selection {outer}/{held_teacher}")
    selected_key, audit = choose_candidate(cells_by_candidate, candidates)
    by_key = {candidate.key: candidate for candidate in candidates}
    return SelectionResult(
        selected_key=selected_key,
        selected=None if selected_key == IDENTITY else by_key[selected_key],
        audit=audit,
        cells=tuple(all_cells),
    )


def evaluate_refit(
    train_queries: Mapping[str, np.ndarray],
    teacher_truth: Mapping[str, np.ndarray],
    eval_queries: Mapping[str, np.ndarray],
    eval_truth: Mapping[str, np.ndarray],
    train_embeddings: np.ndarray,
    train_available: np.ndarray,
    eval_embeddings: np.ndarray,
    eval_available: np.ndarray,
    selected: Candidate | None,
    fit_rows: np.ndarray,
    eval_rows: np.ndarray,
    *,
    max_rank: int,
    ridges: Sequence[float],
    residual_bases: Mapping[str, ResidualBasis] | None = None,
    residual_positions: np.ndarray | None = None,
) -> tuple[dict[str, object], dict[str, ResidualFamily | None]]:
    evaluations: dict[str, object] = {}
    families: dict[str, ResidualFamily | None] = {}
    for outer in OUTER_DONORS:
        basis = (
            None
            if selected is None
            else (
                residual_bases[outer]
                if residual_bases is not None
                else fit_residual_basis(
                    train_queries[outer],
                    teacher_truth,
                    train_available,
                    fit_rows,
                    max_rank=max_rank,
                )
            )
        )
        family = (
            None
            if selected is None
            else fit_residual_family(
                train_embeddings,
                train_available,
                basis,
                ridges=ridges,
            )
        )
        families[outer] = family
        target_truths = (
            TEACHER_DONORS
            if eval_queries is train_queries and eval_truth is teacher_truth
            else (outer,)
        )
        for truth_name in target_truths:
            key = f"{outer}|{truth_name}"
            evaluations[key] = evaluate_candidate(
                eval_queries[outer],
                eval_truth[truth_name],
                eval_embeddings,
                eval_available,
                family,
                selected,
                eval_rows,
                eval_rows,
                residual_positions,
            )
    return evaluations, families


def make_source_strata(
    queries: Mapping[str, np.ndarray],
    source_counts: Mapping[str, np.ndarray],
    partition_labels: np.ndarray,
    *,
    norm_bins: int = 10,
) -> np.ndarray:
    names = tuple(queries)
    n = len(next(iter(queries.values())))
    if any(len(values) != n for values in queries.values()):
        raise ValueError("query lengths do not align for null strata")
    if any(len(values) != n for values in source_counts.values()):
        raise ValueError("source-count lengths do not align for null strata")
    norms = np.mean(
        np.stack([np.linalg.norm(np.asarray(queries[name]), axis=1) for name in names]),
        axis=0,
    )
    labels = np.empty(n, dtype=object)
    coverage_signatures = [
        tuple(int(source_counts[name][index]) for name in names) for index in range(n)
    ]
    groups: dict[tuple[str, tuple[int, ...]], list[int]] = {}
    for index, (partition, coverage) in enumerate(
        zip(partition_labels.tolist(), coverage_signatures, strict=True)
    ):
        groups.setdefault((str(partition), coverage), []).append(index)
    for (partition, coverage), members0 in sorted(groups.items()):
        members = np.asarray(members0, dtype=np.int64)
        ordered = members[np.argsort(norms[members], kind="stable")]
        for decile, bucket in enumerate(
            np.array_split(ordered, min(norm_bins, len(ordered)))
        ):
            for index in bucket:
                labels[int(index)] = (
                    f"{partition}|coverage={','.join(map(str, coverage))}|"
                    f"norm_bin={decile}"
                )
    if any(value is None for value in labels):
        raise AssertionError("null stratum assignment is incomplete")
    return np.asarray(labels, dtype=str)


def permute_embeddings_within_strata(
    embeddings: np.ndarray,
    available: np.ndarray,
    strata: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, dict[str, object]]:
    values = np.asarray(embeddings).copy()
    rng = np.random.default_rng(seed)
    mapping = np.arange(len(values), dtype=np.int64)
    for label in sorted(set(strata.tolist())):
        rows = np.flatnonzero((strata == label) & available)
        if len(rows) > 1:
            mapping[rows] = rng.permutation(rows)
    values[available] = embeddings[mapping[available]]
    values[~available] = np.nan
    if not np.array_equal(np.isfinite(values).all(axis=1), available):
        raise AssertionError("embedding permutation changed the availability mask")
    moved = int(np.sum(mapping[available] != np.flatnonzero(available)))
    return values, {
        "method": (
            "permute ESM2 rows within target split, consensus arm-D "
            "response-norm decile, and two-fold source-coverage signature"
        ),
        "seed": int(seed),
        "strata": int(len(set(strata.tolist()))),
        "available_rows": int(np.sum(available)),
        "moved_rows": moved,
        "mapping_sha256": sha256_array(mapping),
    }


def gaussian_embeddings_like(
    embeddings: np.ndarray,
    available: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, dict[str, object]]:
    values = np.full(np.shape(embeddings), np.nan, dtype=np.float32)
    rng = np.random.default_rng(seed)
    values[available] = rng.standard_normal(
        (int(np.sum(available)), embeddings.shape[1]), dtype=np.float32
    )
    return values, {
        "method": (
            "equal-dimensional iid standard Gaussian target embedding with "
            "the exact ESM2 availability mask"
        ),
        "seed": int(seed),
        "dimensions": int(embeddings.shape[1]),
        "available_rows": int(np.sum(available)),
        "values_sha256": sha256_array(values[available]),
    }


def make_null_embedding(
    kind: str,
    embeddings: np.ndarray,
    available: np.ndarray,
    strata: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, dict[str, object]]:
    if kind == "permuted":
        return permute_embeddings_within_strata(
            embeddings, available, strata, seed=seed
        )
    if kind == "gaussian":
        return gaussian_embeddings_like(embeddings, available, seed=seed)
    raise ValueError(f"unknown null kind {kind!r}")


def stratum_summaries(
    per_target_delta: np.ndarray,
    source_count: np.ndarray,
    *,
    complete_count: int,
) -> dict[str, dict[str, object]]:
    definitions = {
        "complete_source": source_count == complete_count,
        "partial_source": source_count < complete_count,
        "direct_row": source_count >= 1,
        "missing_source": source_count == 0,
    }
    result: dict[str, dict[str, object]] = {}
    for label, mask in definitions.items():
        result[label] = {
            "count": int(np.sum(mask)),
            "delta": (
                None
                if not np.any(mask)
                else float(np.mean(np.asarray(per_target_delta)[mask]))
            ),
        }
    return result


def aggregate_stratum_delta(
    fold_evaluations: Mapping[str, Mapping[str, object]],
    source_counts: Mapping[str, np.ndarray],
    complete_counts: Mapping[str, int],
    label: str,
) -> tuple[float | None, int]:
    values: list[np.ndarray] = []
    for held in fold_evaluations:
        count = source_counts[held]
        complete_count = int(complete_counts[held])
        if label == "complete_source":
            mask = count == complete_count
        elif label == "partial_source":
            mask = count < complete_count
        elif label == "direct_row":
            mask = count >= 1
        elif label == "missing_source":
            mask = count == 0
        else:
            raise ValueError(label)
        if np.any(mask):
            values.append(np.asarray(fold_evaluations[held]["per_target_delta"])[mask])
    if not values:
        return None, 0
    concatenated = np.concatenate(values)
    return float(np.mean(concatenated)), int(len(concatenated))


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
    legacy_mask = np.asarray(
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
    if len(master) != len(legacy_mask):
        raise ValueError("truth-cache legacy mask width mismatch")
    if not set(master) <= set(official):
        raise ValueError("truth-cache feature lies outside official genes")
    if set(master) & set(panel):
        raise ValueError("VCC coordinate entered final panel score axis")
    train_set = set(train_targets)
    dev_mask = legacy_mask
    dev_genes = tuple(gene for gene, keep in zip(master, dev_mask, strict=True) if keep)
    teacher_gene_sets = {
        donor: set(
            read_lines(
                donor_root / donor / "genes.txt",
                require_unique=False,
            )
        )
        for donor in TEACHER_DONORS
    }
    expected_dev_genes = (
        set(master)
        .intersection(*(teacher_gene_sets[donor] for donor in TEACHER_DONORS))
        .difference(train_set)
    )
    if set(dev_genes) != expected_dev_genes:
        raise ValueError(
            "teacher-common score mask does not equal the measured-gene "
            "intersection after target-coordinate exclusion"
        )
    panel_set = set(panel)
    official_pos = {gene: index for index, gene in enumerate(official)}
    if len(official_pos) != len(official):
        raise ValueError("official gene axis contains duplicate symbols")
    donor_axes = []
    donor_axis_audits: dict[str, object] = {}
    for donor in DEVELOPMENT_FOLDS:
        donor_gene_rows = read_lines(
            donor_root / donor / "genes.txt",
            require_unique=False,
        )
        donor_gene_lookup, donor_axis_audit = last_occurrence_lookup(donor_gene_rows)
        donor_axis_audits[donor] = donor_axis_audit
        donor_axes.append(
            tuple(
                gene
                for gene in official
                if gene in donor_gene_lookup and gene not in panel_set
            )
        )
    panel_genes = donor_axes[0]
    if any(axis != panel_genes for axis in donor_axes[1:]):
        raise ValueError("development held-donor score axes are not identical")
    panel_score_official = np.asarray(
        [official_pos[gene] for gene in panel_genes], dtype=np.int64
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
    if len(panel_genes) != expected_count:
        raise ValueError(
            f"final score axis {len(panel_genes)} != frozen {expected_count}"
        )
    if sha256_lines(panel_genes) != final_axis.get("gene_order_sha256"):
        raise ValueError("final score-gene order/hash changed")
    if sha256_array(panel_score_official) != final_axis.get(
        "official_positions_sha256"
    ):
        raise ValueError("final score official-position order/hash changed")
    if np.any(
        legacy_mask & np.asarray([gene in train_set for gene in master], dtype=bool)
    ):
        raise ValueError("training target coordinate entered legacy truth mask")
    panel_position = {gene: index for index, gene in enumerate(panel_genes)}
    missing_dev = [gene for gene in dev_genes if gene not in panel_position]
    if missing_dev:
        raise ValueError(
            f"{len(missing_dev)} development score genes are absent from the "
            "frozen final score axis"
        )
    dev_to_panel = np.asarray(
        [panel_position[gene] for gene in dev_genes], dtype=np.int64
    )
    if set(dev_genes) & (set(panel) | train_set):
        raise AssertionError("held target coordinate entered development score axis")
    historical_coordinates = set(panel_genes) & train_set
    expected_historical = train_set & set(official_pos)
    if historical_coordinates != expected_historical:
        raise AssertionError(
            "final panel axis did not restore every available development coordinate"
        )
    return (
        master,
        train_targets,
        outer_targets,
        FeatureAxes(
            master_genes=master,
            dev_score_mask=dev_mask,
            dev_score_genes=dev_genes,
            panel_score_genes=panel_genes,
            panel_score_official=panel_score_official,
            dev_to_panel_positions=dev_to_panel,
            donor_axis_audits=donor_axis_audits,
        ),
    )


def source_effect_clusters(
    query: np.ndarray,
    *,
    reciprocal_k: int,
    cosine_threshold: float,
) -> tuple[np.ndarray, dict[str, object]]:
    """Cluster one outer fold from that fold's source-only arm-D queries."""
    unit = row_normalize(query)
    similarities = unit @ unit.T
    np.fill_diagonal(similarities, -np.inf)
    k = min(reciprocal_k, len(unit) - 1)
    neighbors = np.argsort(-similarities, axis=1, kind="stable")[:, :k]
    directed = np.zeros((len(unit), len(unit)), dtype=bool)
    directed[np.arange(len(unit))[:, None], neighbors] = True
    reciprocal = directed & directed.T & (similarities >= cosine_threshold)
    parent = np.arange(len(unit), dtype=np.int64)

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
    roots = [find(index) for index in range(len(unit))]
    root_order = {root: index for index, root in enumerate(sorted(set(roots)))}
    labels = np.asarray([root_order[root] for root in roots], dtype=np.int64)
    sizes = np.bincount(labels)
    return labels, {
        "method": ("outer-fold source-only arm-D reciprocal-kNN connected components"),
        "reciprocal_k": int(k),
        "cosine_threshold": float(cosine_threshold),
        "clusters": int(len(sizes)),
        "singleton_clusters": int(np.sum(sizes == 1)),
        "largest_cluster": int(np.max(sizes)),
        "cluster_sizes_sha256": sha256_array(sizes),
    }


def outer_fold_cluster_bootstrap(
    deltas_by_fold: Mapping[str, np.ndarray],
    cluster_labels_by_fold: Mapping[str, np.ndarray],
    *,
    replicates: int,
    seed: int,
) -> dict[str, object]:
    folds = tuple(deltas_by_fold)
    if not folds or set(folds) != set(cluster_labels_by_fold):
        raise ValueError("bootstrap folds and cluster-label folds differ")
    rng = np.random.default_rng(seed)
    draws = np.empty(replicates, dtype=np.float64)
    fold_members: dict[str, list[np.ndarray]] = {}
    for held in folds:
        values = np.asarray(deltas_by_fold[held])
        labels = np.asarray(cluster_labels_by_fold[held], dtype=np.int64)
        if values.shape != labels.shape:
            raise ValueError(f"{held}: cluster labels do not align")
        fold_members[held] = [
            np.flatnonzero(labels == cluster)
            for cluster in sorted(set(labels.tolist()))
        ]
    for replicate in range(replicates):
        fold_means = []
        for held in folds:
            members = fold_members[held]
            sampled = rng.integers(0, len(members), size=len(members))
            rows = np.concatenate([members[index] for index in sampled])
            fold_means.append(float(np.mean(deltas_by_fold[held][rows])))
        draws[replicate] = float(np.mean(fold_means))
    point = float(np.mean([np.mean(deltas_by_fold[held]) for held in folds]))
    return {
        "method": (
            "independent outer-fold source-only cluster bootstrap; sample "
            "clusters within each fold, then average fold means"
        ),
        "replicates": int(replicates),
        "seed": int(seed),
        "point_estimate": point,
        "lower_95": float(np.quantile(draws, 0.025)),
        "upper_95": float(np.quantile(draws, 0.975)),
        "probability_positive": float(np.mean(draws > 0.0)),
        "clusters_by_fold": {held: len(fold_members[held]) for held in folds},
    }


def cluster_dependence_diagnostics(
    deltas: np.ndarray,
    labels: np.ndarray,
) -> dict[str, object]:
    clusters = sorted(set(map(int, labels.tolist())))
    cluster_sums = np.asarray(
        [float(np.sum(deltas[labels == cluster])) for cluster in clusters]
    )
    positive = np.maximum(cluster_sums, 0.0)
    positive_total = float(np.sum(positive))
    largest_share = (
        float(np.max(positive) / positive_total) if positive_total > 0 else 1.0
    )
    if len(clusters) < 2:
        return {
            "largest_positive_cluster_share": largest_share,
            "leave_one_cluster_out_minimum_delta": None,
            "passes_no_obvious_single_family_dependence": False,
        }
    leave_one_out = [float(np.mean(deltas[labels != cluster])) for cluster in clusters]
    return {
        "largest_positive_cluster_share": largest_share,
        "leave_one_cluster_out_minimum_delta": float(min(leave_one_out)),
        "leave_one_cluster_out_maximum_delta": float(max(leave_one_out)),
        "passes_no_obvious_single_family_dependence": bool(
            largest_share <= 0.50 and min(leave_one_out) >= 0.0
        ),
    }


def verify_manifest_inputs(
    args: argparse.Namespace,
    manifest: Mapping[str, object],
    *,
    include_outer_truth: bool,
) -> None:
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
        actual = sha256_file(args.truth_cache / filename)
        expected = str(expected0)
        if actual != expected:
            raise ValueError(f"truth-cache/{filename}: SHA-256 {actual} != {expected}")
    for donor in TEACHER_DONORS:
        for filename, expected0 in manifest["teacher_truth_donors"][donor].items():
            path = args.donor_root / donor / filename
            expected = str(expected0)
            actual = sha256_file(path)
            if actual != expected:
                raise ValueError(
                    f"teacher truth {donor}/{filename}: SHA-256 {actual} != {expected}"
                )
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
    verify_embedding_artifact(args.embedding_parquet, manifest)


def validate_query_cache(
    query_cache: Path,
    axes,
    train_targets: Sequence[str],
    outer_targets: Sequence[str],
    panel: Sequence[str],
    manifest: Mapping[str, object],
) -> tuple[dict[str, object], dict[str, tuple[str, ...]]]:
    metadata = json.loads((query_cache / "metadata.json").read_text())
    expected = {
        "schema_version": CACHE_SCHEMA,
        "experiment": "082-biological-target-kernel",
        "arm": "D",
        "development_folds": list(DEVELOPMENT_FOLDS),
        "final_folds": list(FINAL_FOLDS),
        "valid_outer_donors": list(OUTER_DONORS),
        "training_targets": len(train_targets),
        "outer_targets": len(outer_targets),
        "training_outer_overlap": 0,
        "vcc_targets": len(panel),
        "master_genes": len(axes.master_genes),
        "dev_score_genes": len(axes.dev_score_genes),
        "panel_score_genes": len(axes.panel_score_genes),
        "dev_score_gene_order_sha256": sha256_lines(axes.dev_score_genes),
        "panel_score_gene_order_sha256": sha256_lines(axes.panel_score_genes),
        "panel_score_official_positions_sha256": sha256_array(
            axes.panel_score_official
        ),
        "dev_to_panel_positions_sha256": sha256_array(axes.dev_to_panel_positions),
        "train_order_sha256": sha256_lines(train_targets),
        "outer_order_sha256": sha256_lines(outer_targets),
        "panel_order_sha256": sha256_lines(panel),
        "final_score_axis": {
            "gene_count": len(axes.panel_score_genes),
            "gene_order_sha256": sha256_lines(axes.panel_score_genes),
            "official_positions_sha256": sha256_array(axes.panel_score_official),
            "development_donor_axes_equal": True,
            "excludes_exactly_vcc_coordinates_within_held_support": True,
            "donor_gene_axes": axes.donor_axis_audits,
        },
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"query-cache metadata {key}={metadata.get(key)!r}, expected {value!r}"
            )
    implementation = metadata.get("implementation", {})
    if (
        implementation.get("exp062_run_sha256")
        != manifest["dependencies"]["exp062_run_sha256"]
    ):
        raise ValueError("experiment-062 dependency hash mismatch")
    if (
        implementation.get("exp064_run_sha256")
        != manifest["dependencies"]["exp064_run_sha256"]
    ):
        raise ValueError("experiment-064 dependency hash mismatch")
    if implementation.get("exact_arm_d_formula") is not True:
        raise ValueError("exact arm-D formula audit failed")
    contract = metadata.get("coordinate_contract", {})
    if (
        contract.get("development_blocked", {}).get("before_dense_whitening")
        is not True
    ):
        raise ValueError("development coordinates were not blocked pre-whitening")
    if (
        contract.get("final_panel_blocked", {}).get(
            "score_axis_excludes_exactly_vcc_targets"
        )
        is not True
    ):
        raise ValueError("final panel score-axis contract failed")
    mapping_path = query_cache / "dev_to_panel_positions.npy"
    mapping = np.load(mapping_path, allow_pickle=False)
    if not np.array_equal(mapping, axes.dev_to_panel_positions):
        raise ValueError("dev-to-panel axis mapping mismatch")
    fold_panel_targets: dict[str, tuple[str, ...]] = {}
    for held in FINAL_FOLDS:
        fold = metadata["folds"][held]
        operator = manifest["operators"][held]
        target_path = query_cache / held / "panel_targets.txt"
        saved_targets = read_lines(target_path)
        if held == KOLF_FOLD:
            expected_targets = tuple(panel)
            if saved_targets != expected_targets:
                raise ValueError("KOLF source query must cover all 300 VCC targets")
        else:
            gallery = manifest["held_panel_galleries"][held]
            if int(gallery["targets"]) != len(saved_targets):
                raise ValueError(f"{held}: complete panel target count changed")
            if str(gallery["target_order_sha256"]) != sha256_lines(saved_targets):
                raise ValueError(f"{held}: complete panel target order changed")
            if not set(saved_targets) <= set(panel):
                raise ValueError(f"{held}: non-panel target entered final gallery")
        if int(fold.get("panel_targets", -1)) != len(saved_targets):
            raise ValueError(f"{held}: panel target count audit changed")
        if fold.get("panel_target_order_sha256") != sha256_lines(saved_targets):
            raise ValueError(f"{held}: panel target order audit changed")
        if fold.get("panel_targets_file_sha256") != sha256_file(target_path):
            raise ValueError(f"{held}: panel target file hash changed")
        fold_panel_targets[held] = saved_targets
        if held in set(map(str, fold["source_donors"])):
            raise ValueError(f"{held}: held donor entered arm-D sources")
        if tuple(fold["source_donors"]) != tuple(operator["source_donors"]):
            raise ValueError(f"{held}: source donor mismatch")
        if fold["operator_npz_sha256"] != operator["npz_sha256"]:
            raise ValueError(f"{held}: operator NPZ hash mismatch")
        if fold["operator_json_sha256"] != operator["json_sha256"]:
            raise ValueError(f"{held}: operator JSON hash mismatch")
        if fold.get("held_donor_used_in_prediction") is not False:
            raise ValueError(f"{held}: held donor entered prediction")
        if fold.get("held_donor_used_in_operator") is not False:
            raise ValueError(f"{held}: held donor entered operator")
        if fold.get("held_donor_delta_loaded_in_own_fold") is not False:
            raise ValueError(f"{held}: held donor delta was loaded in its own fold")
        if held in DEVELOPMENT_FOLDS and fold.get("complete_held_panel") is not True:
            raise ValueError(f"{held}: final query is not marked as the complete panel")
        for label, count in (
            ("train", len(train_targets)),
            ("panel", len(saved_targets)),
        ):
            query_path = query_cache / held / f"{label}_query.npy"
            query = np.load(query_path, mmap_mode="r", allow_pickle=False)
            width = (
                len(axes.dev_score_genes)
                if label == "train"
                else len(axes.panel_score_genes)
            )
            if query.shape != (count, width):
                raise ValueError(f"{query_path}: unexpected shape {query.shape}")
            if sha256_file(query_path) != fold["panels"][label]["query_sha256"]:
                raise ValueError(f"{query_path}: query hash mismatch")
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
                    if fold["panels"][label].get(key) != value:
                        raise ValueError(
                            f"{held}: panel query axis audit {key} changed"
                        )
            count_path = query_cache / held / f"{label}_source_count.npy"
            source_count = np.load(count_path, allow_pickle=False)
            if source_count.shape != (count,):
                raise ValueError(f"{count_path}: unexpected shape")
            if not np.all(np.isin(source_count, range(len(fold["source_donors"]) + 1))):
                raise ValueError(f"{count_path}: invalid source count")
            if sha256_file(count_path) != fold["panels"][label]["source_count_sha256"]:
                raise ValueError(f"{count_path}: source-count hash mismatch")
        restriction = fold["panels"]["train"]["operator_restriction"]
        if restriction.get("all_nonpanel_target_coordinates_excluded") is not True:
            raise ValueError(f"{held}: non-panel target entered training operator")
        if restriction.get("nonpanel_target_overlap_with_active_operator") != 0:
            raise ValueError(f"{held}: non-panel/operator overlap is nonzero")
        if restriction.get("dense_zero_equivalence_proved") is not True:
            raise ValueError(f"{held}: dense-zero equivalence was not proved")
        if held in DEVELOPMENT_FOLDS:
            reproduction = fold.get("exp064_reproduction", {})
            if reproduction.get("passed") is not True:
                raise ValueError(f"{held}: exp064 arm-D reproduction did not pass")
            if reproduction.get("saved_query_score_axis_exact") is not True:
                raise ValueError(f"{held}: arm-D audit did not cover saved panel_query")
            if reproduction.get("retrieval_genes") != len(axes.panel_score_genes):
                raise ValueError(f"{held}: arm-D retrieval-axis width changed")
            if reproduction.get("retrieval_gene_order_sha256") != sha256_lines(
                axes.panel_score_genes
            ):
                raise ValueError(f"{held}: arm-D retrieval gene order changed")
    expected_counts = {held: len(fold_panel_targets[held]) for held in FINAL_FOLDS}
    expected_hashes = {
        held: sha256_lines(fold_panel_targets[held]) for held in FINAL_FOLDS
    }
    if metadata.get("panel_query_target_counts") != expected_counts:
        raise ValueError("fold-specific panel query target counts changed")
    if metadata.get("panel_query_target_order_sha256") != expected_hashes:
        raise ValueError("fold-specific panel query target hashes changed")
    return metadata, fold_panel_targets


def load_truth(
    donor_root: Path,
    donor: str,
    targets: Sequence[str],
    score_genes: Sequence[str],
    *,
    chunk_rows: int,
) -> np.ndarray:
    values, _audit = load_compact_panel_truth(
        donor_root,
        donor,
        targets,
        score_genes,
        chunk_rows=chunk_rows,
    )
    return values


def load_compact_panel_truth(
    donor_root: Path,
    donor: str,
    targets: Sequence[str],
    score_genes: Sequence[str],
    *,
    chunk_rows: int,
    complete_panel: Sequence[str] | None = None,
    require_complete_genes: bool = True,
) -> tuple[np.ndarray, dict[str, object]]:
    """Load held donor truth on the complete frozen score axis."""
    root = donor_root / donor
    genes = read_lines(root / "genes.txt", require_unique=False)
    perts = read_lines(root / "perts.txt")
    gene_last, gene_axis_audit = last_occurrence_lookup(genes)
    missing_genes = [gene for gene in score_genes if gene not in gene_last]
    if require_complete_genes and missing_genes:
        raise ValueError(f"{donor}: {len(missing_genes)} frozen score genes are absent")
    pert_pos = {target: index for index, target in enumerate(perts)}
    missing_targets = [target for target in targets if target not in pert_pos]
    if missing_targets:
        raise ValueError(f"{donor}: {len(missing_targets)} held targets are absent")
    if complete_panel is not None:
        complete_targets = tuple(
            target for target in complete_panel if target in pert_pos
        )
        if tuple(targets) != complete_targets:
            raise ValueError(
                f"{donor}: final target list is not the complete held-panel gallery"
            )
    available_gene_rows = np.asarray(
        [index for index, gene in enumerate(score_genes) if gene in gene_last],
        dtype=np.int64,
    )
    columns = np.asarray(
        [gene_last[score_genes[index]] for index in available_gene_rows],
        dtype=np.int64,
    )
    rows = np.asarray([pert_pos[target] for target in targets], dtype=np.int64)
    delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape != (len(perts), len(genes)):
        raise ValueError(f"{donor}: compact delta shape does not match its axes")
    if chunk_rows <= 0:
        raise ValueError("chunk_rows must be positive")
    values = np.zeros((len(rows), len(score_genes)), dtype=np.float64)
    for lo in range(0, len(rows), chunk_rows):
        hi = min(lo + chunk_rows, len(rows))
        values[np.ix_(np.arange(lo, hi), available_gene_rows)] = np.asarray(
            delta[np.ix_(rows[lo:hi], columns)],
            dtype=np.float64,
        )
    mmap = getattr(delta, "_mmap", None)
    if mmap is not None:
        mmap.close()
    if not np.isfinite(values).all():
        raise ValueError(f"{donor}: held truth contains non-finite values")
    return values, {
        "targets": len(targets),
        "score_genes": len(score_genes),
        "available_score_genes": int(len(available_gene_rows)),
        "missing_score_genes": int(len(missing_genes)),
        "missing_score_genes_sha256": sha256_lines(sorted(missing_genes)),
        "target_order_sha256": sha256_lines(targets),
        "score_gene_order_sha256": sha256_lines(score_genes),
        "representation": "raw control-referenced donor delta, uncentered",
        "complete_held_panel": complete_panel is not None,
        "gene_axis": gene_axis_audit,
        "selected_source_columns_sha256": sha256_array(columns),
        "used_for_selection": False,
        "used_for_fitting": False,
    }


def partition_labels(
    length: int,
    splits: Mapping[str, np.ndarray],
) -> np.ndarray:
    labels = np.empty(length, dtype=object)
    labels[:] = None
    for name, rows in splits.items():
        labels[np.asarray(rows, dtype=np.int64)] = name
    if any(value is None for value in labels):
        raise AssertionError("target split labels are incomplete")
    return np.asarray(labels, dtype=str)


def null_pretest_distribution(
    kinds: Sequence[str],
    replicates: int,
    *,
    train_embeddings: np.ndarray,
    train_available: np.ndarray,
    train_strata: np.ndarray,
    train_queries: Mapping[str, np.ndarray],
    teacher_truth: Mapping[str, np.ndarray],
    splits: Mapping[str, np.ndarray],
    candidates: Sequence[Candidate],
    max_rank: int,
    ridges: Sequence[float],
    null_seed: int,
    max_rss_gib: float,
    inner_residual_bases: Mapping[tuple[str, str], ResidualBasis],
    pretest_residual_bases: Mapping[str, ResidualBasis],
) -> tuple[dict[str, list[dict[str, object]]], dict[tuple[str, int], str]]:
    records: dict[str, list[dict[str, object]]] = {kind: [] for kind in kinds}
    selected: dict[tuple[str, int], str] = {}
    fit_tune = np.sort(np.concatenate([splits["fit"], splits["tune"]]))
    for kind in kinds:
        for replicate in range(replicates):
            seed = stable_seed(null_seed, kind, replicate, "train")
            null_embedding, null_audit = make_null_embedding(
                kind,
                train_embeddings,
                train_available,
                train_strata,
                seed=seed,
            )
            selection = run_selection(
                train_queries,
                teacher_truth,
                null_embedding,
                train_available,
                splits,
                candidates,
                max_rank=max_rank,
                ridges=ridges,
                max_rss_gib=max_rss_gib,
                residual_bases=inner_residual_bases,
            )
            evaluated, _ = evaluate_refit(
                train_queries,
                teacher_truth,
                train_queries,
                teacher_truth,
                null_embedding,
                train_available,
                null_embedding,
                train_available,
                selection.selected,
                fit_tune,
                splits["test"],
                max_rank=max_rank,
                ridges=ridges,
                residual_bases=pretest_residual_bases,
            )
            deltas = [float(value["delta"]) for value in evaluated.values()]
            selected[(kind, replicate)] = selection.selected_key
            records[kind].append(
                {
                    "replicate": replicate,
                    "selected": selection.selected_key,
                    "macro_test_delta": float(np.mean(deltas)),
                    "minimum_test_cell_delta": float(min(deltas)),
                    "embedding_audit": null_audit,
                    "selection": selection.audit,
                    "complete_grid_repeated": True,
                    "candidate_count_including_identity": len(candidates) + 1,
                }
            )
            enforce_memory_limit(max_rss_gib, f"{kind} null pretest {replicate}")
    return records, selected


def build_refit_residual_bases(
    train_queries: Mapping[str, np.ndarray],
    teacher_truth: Mapping[str, np.ndarray],
    available: np.ndarray,
    fit_rows: np.ndarray,
    *,
    max_rank: int,
    folds: Sequence[str] = OUTER_DONORS,
) -> dict[str, ResidualBasis]:
    return {
        fold: fit_residual_basis(
            train_queries[fold],
            teacher_truth,
            available,
            fit_rows,
            max_rank=max_rank,
        )
        for fold in folds
    }


def build_frozen_predictors(
    embeddings: np.ndarray,
    available: np.ndarray,
    selected: Candidate | None,
    bases: Mapping[str, ResidualBasis],
    *,
    ridges: Sequence[float],
) -> dict[str, FrozenPredictor | None]:
    if selected is None:
        return {fold: None for fold in bases}
    return {
        fold: freeze_predictor(
            fit_residual_family(
                embeddings,
                available,
                bases[fold],
                ridges=[selected.ridge],
            ),
            selected,
        )
        for fold in bases
    }


def save_frozen_predictor(
    path: Path,
    predictor: FrozenPredictor | None,
) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if predictor is None:
        marker = path.with_suffix(".identity.json")
        marker.write_text(
            json.dumps(
                {"candidate": IDENTITY, "sha256": predictor_sha256(None)},
                indent=2,
            )
            + "\n"
        )
        return {
            "candidate": IDENTITY,
            "path": str(marker),
            "file_sha256": sha256_file(marker),
            "weights_sha256": predictor_sha256(None),
        }
    np.savez_compressed(
        path,
        candidate=np.asarray(
            [
                predictor.candidate.rank,
                predictor.candidate.ridge,
                predictor.candidate.beta,
            ],
            dtype=np.float64,
        ),
        embedding_mean=predictor.embedding_mean,
        embedding_scale=predictor.embedding_scale,
        primal_weight=predictor.primal_weight,
        residual_mean=predictor.residual_mean,
        components=predictor.components,
    )
    return {
        "candidate": asdict(predictor.candidate),
        "path": str(path),
        "file_sha256": sha256_file(path),
        "weights_sha256": predictor_sha256(predictor),
    }


def load_frozen_predictor(
    path: Path,
    expected: Mapping[str, object],
) -> FrozenPredictor | None:
    if expected["candidate"] == IDENTITY:
        marker = path.with_suffix(".identity.json")
        if sha256_file(marker) != expected["file_sha256"]:
            raise ValueError("identity freeze marker hash mismatch")
        return None
    if sha256_file(path) != expected["file_sha256"]:
        raise ValueError(f"{path}: frozen predictor file hash mismatch")
    with np.load(path, allow_pickle=False) as archive:
        raw_candidate = np.asarray(archive["candidate"], dtype=np.float64)
        candidate = Candidate(
            int(raw_candidate[0]),
            float(raw_candidate[1]),
            float(raw_candidate[2]),
        )
        predictor = FrozenPredictor(
            candidate=candidate,
            embedding_mean=np.asarray(archive["embedding_mean"], dtype=np.float32),
            embedding_scale=np.asarray(archive["embedding_scale"], dtype=np.float32),
            primal_weight=np.asarray(archive["primal_weight"], dtype=np.float32),
            residual_mean=np.asarray(archive["residual_mean"], dtype=np.float32),
            components=np.asarray(archive["components"], dtype=np.float32),
            audit={"loaded_from_freeze": str(path)},
        )
    if predictor_sha256(predictor) != expected["weights_sha256"]:
        raise ValueError(f"{path}: frozen predictor weight hash mismatch")
    return predictor


def load_kolf_truth(
    root: Path,
    panel: Sequence[str],
    panel_score_genes: Sequence[str],
) -> dict[str, object]:
    targets = read_lines(root / "perts.txt")
    genes = read_lines(root / "genes.txt")
    panel_pos = {target: index for index, target in enumerate(panel)}
    shared_targets = tuple(target for target in panel if target in set(targets))
    if tuple(targets) != shared_targets:
        raise ValueError("KOLF compact target order is not the VCC-panel subset order")
    query_rows = np.asarray(
        [panel_pos[target] for target in shared_targets], dtype=np.int64
    )
    gene_pos = {gene: index for index, gene in enumerate(genes)}
    query_gene_pos = {gene: index for index, gene in enumerate(panel_score_genes)}
    score_genes = tuple(gene for gene in panel_score_genes if gene in gene_pos)
    query_columns = np.asarray(
        [query_gene_pos[gene] for gene in score_genes], dtype=np.int64
    )
    truth_columns = np.asarray([gene_pos[gene] for gene in score_genes], dtype=np.int64)
    delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape != (len(targets), len(genes)):
        raise ValueError("KOLF delta shape does not match its axes")
    values = np.asarray(delta[:, truth_columns], dtype=np.float64)
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
            "representation": "raw decoded control-referenced quantized z-score, uncentered",
            "used_for_selection": False,
            "used_for_fitting": False,
        },
    }


def evaluate_frozen_fold(
    predictor: FrozenPredictor | None,
    query: np.ndarray,
    truth: np.ndarray,
    embeddings: np.ndarray,
    available: np.ndarray,
    query_rows: np.ndarray,
    positive_columns: np.ndarray,
    query_columns: np.ndarray,
    truth_columns: np.ndarray,
    *,
    residual_positions: np.ndarray,
) -> dict[str, object]:
    rows = np.asarray(query_rows, dtype=np.int64)
    qcols = np.asarray(query_columns, dtype=np.int64)
    tcols = np.asarray(truth_columns, dtype=np.int64)
    local_query_full = np.asarray(query[rows], dtype=np.float64)
    local_embeddings = np.asarray(embeddings[rows])
    local_available = np.asarray(available[rows], dtype=bool)
    baseline = local_query_full[:, qcols]
    if predictor is None:
        adapted = baseline.copy()
    else:
        adapted_full = explicit_candidate_query(
            local_query_full,
            local_embeddings,
            local_available,
            predictor,
            predictor.candidate,
            residual_positions,
        )
        adapted = adapted_full[:, qcols]
    truth_scored = np.asarray(truth[:, tcols], dtype=np.float64)
    baseline_scores = row_normalize(baseline) @ row_normalize(truth_scored).T
    adapted_scores = row_normalize(adapted) @ row_normalize(truth_scored).T
    baseline_pds = pds_per_target_from_scores(
        baseline_scores, np.asarray(positive_columns, dtype=np.int64)
    )
    adapted_pds = pds_per_target_from_scores(
        adapted_scores, np.asarray(positive_columns, dtype=np.int64)
    )
    angles = np.degrees(
        np.arccos(
            np.clip(
                np.sum(row_normalize(baseline) * row_normalize(adapted), axis=1),
                -1.0,
                1.0,
            )
        )
    )
    return {
        "baseline": float(np.mean(baseline_pds)),
        "adapted": float(np.mean(adapted_pds)),
        "delta": float(np.mean(adapted_pds - baseline_pds)),
        "per_target_baseline": baseline_pds,
        "per_target_adapted": adapted_pds,
        "per_target_delta": adapted_pds - baseline_pds,
        "baseline_norms": np.linalg.norm(baseline, axis=1),
        "adapted_norms": np.linalg.norm(adapted, axis=1),
        "angles_degrees": angles,
        "median_angle_degrees": float(np.median(angles)),
        "max_angle_degrees": float(np.max(angles)),
        "embedding_available": local_available,
        "gallery_size": int(len(truth_scored)),
    }


def evaluate_frozen_predictors(
    predictors: Mapping[str, FrozenPredictor | None],
    queries: Mapping[str, np.ndarray],
    truths: Mapping[str, np.ndarray],
    embeddings: np.ndarray,
    available: np.ndarray,
    rows: np.ndarray,
    *,
    residual_positions: np.ndarray,
) -> dict[str, dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    for held in OUTER_DONORS:
        predictor = predictors[held]
        candidate = None if predictor is None else predictor.candidate
        output[held] = evaluate_candidate(
            queries[held],
            truths[held],
            embeddings,
            available,
            predictor,
            candidate,
            rows,
            rows,
            residual_positions,
        )
    return output


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
    query_metadata, fold_panel_targets = validate_query_cache(
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
    all_embedding_targets = (*train_targets, *panel)
    embedding_table = load_target_embeddings(
        args.embedding_parquet,
        all_embedding_targets,
        expected_rows=args.expected_embedding_rows,
        expected_dim=args.expected_embedding_dim,
        column_chunk=args.embedding_column_chunk,
    )
    split_at = len(train_targets)
    train_embeddings = embedding_table.values[:split_at]
    train_available = embedding_table.available[:split_at]
    panel_embeddings = embedding_table.values[split_at:]
    panel_available = embedding_table.available[split_at:]
    panel_pos = {target: index for index, target in enumerate(panel)}
    fold_panel_rows = {
        held: np.asarray(
            [panel_pos[target] for target in fold_panel_targets[held]],
            dtype=np.int64,
        )
        for held in FINAL_FOLDS
    }
    fold_panel_embeddings = {
        held: panel_embeddings[rows] for held, rows in fold_panel_rows.items()
    }
    fold_panel_available = {
        held: panel_available[rows] for held, rows in fold_panel_rows.items()
    }
    splits = deterministic_target_split(train_targets, seed=args.split_seed)
    candidates = build_candidates(args.ranks, args.ridges, args.betas)
    max_rank = max(args.ranks)
    teacher_truth = {
        donor: load_truth(
            args.donor_root,
            donor,
            train_targets,
            axes.dev_score_genes,
            chunk_rows=args.truth_chunk_rows,
        )
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
    train_source_counts = {
        held: np.asarray(
            np.load(
                args.query_cache / held / "train_source_count.npy",
                allow_pickle=False,
            ),
            dtype=np.uint8,
        )
        for held in FINAL_FOLDS
    }
    train_strata = make_source_strata(
        {held: train_queries[held] for held in OUTER_DONORS},
        {held: train_source_counts[held] for held in OUTER_DONORS},
        partition_labels(len(train_targets), splits),
        norm_bins=args.null_norm_bins,
    )
    enforce_memory_limit(args.max_rss_gib, "inputs")

    inner_residual_bases = build_inner_residual_bases(
        train_queries,
        teacher_truth,
        train_available,
        splits["fit"],
        max_rank=max_rank,
    )
    real_selection = run_selection(
        train_queries,
        teacher_truth,
        train_embeddings,
        train_available,
        splits,
        candidates,
        max_rank=max_rank,
        ridges=args.ridges,
        max_rss_gib=args.max_rss_gib,
        residual_bases=inner_residual_bases,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selection_freeze: dict[str, object] = {
        "selected": real_selection.selected_key,
        "candidate": (
            None if real_selection.selected is None else asdict(real_selection.selected)
        ),
        "candidate_grid": {
            "identity": True,
            "ranks": list(args.ranks),
            "ridges": list(args.ridges),
            "betas": list(args.betas),
            "count_including_identity": len(candidates) + 1,
        },
        "target_split": {
            name: {
                "count": int(len(rows)),
                "targets_sha256": sha256_lines(
                    train_targets[index] for index in rows.tolist()
                ),
            }
            for name, rows in splits.items()
        },
        "selection": real_selection.audit,
        "complete_selection_nulls": {},
        "pretest_gate": None,
        "outer_truth_opened": False,
    }
    selection_path = args.output_dir / "selection_freeze.json"

    fit_tune = np.sort(np.concatenate([splits["fit"], splits["tune"]]))
    pretest_residual_bases = build_refit_residual_bases(
        train_queries,
        teacher_truth,
        train_available,
        fit_tune,
        max_rank=max_rank,
    )
    real_pretest, _ = evaluate_refit(
        train_queries,
        teacher_truth,
        train_queries,
        teacher_truth,
        train_embeddings,
        train_available,
        train_embeddings,
        train_available,
        real_selection.selected,
        fit_tune,
        splits["test"],
        max_rank=max_rank,
        ridges=args.ridges,
        residual_bases=pretest_residual_bases,
        residual_positions=np.arange(len(axes.dev_score_genes), dtype=np.int64),
    )
    pretest_cells = [
        {
            "outer_donor": key.split("|", 1)[0],
            "teacher": key.split("|", 1)[1],
            **summarize_evaluation(value),
        }
        for key, value in real_pretest.items()
    ]
    null_records, null_selected = null_pretest_distribution(
        ("permuted", "gaussian"),
        args.null_replicates,
        train_embeddings=train_embeddings,
        train_available=train_available,
        train_strata=train_strata,
        train_queries=train_queries,
        teacher_truth=teacher_truth,
        splits=splits,
        candidates=candidates,
        max_rank=max_rank,
        ridges=args.ridges,
        null_seed=args.null_seed,
        max_rss_gib=args.max_rss_gib,
        inner_residual_bases=inner_residual_bases,
        pretest_residual_bases=pretest_residual_bases,
    )
    pretest_macro = float(np.mean([cell["delta"] for cell in pretest_cells]))
    pretest_minimum = float(min(cell["delta"] for cell in pretest_cells))
    pretest_fold_angles = {
        outer: float(
            max(
                cell["median_angle_degrees"]
                for cell in pretest_cells
                if cell["outer_donor"] == outer
            )
        )
        for outer in OUTER_DONORS
    }
    null_p95 = {
        kind: float(
            np.quantile([record["macro_test_delta"] for record in records], 0.95)
        )
        for kind, records in null_records.items()
    }
    selected_cells = (
        []
        if real_selection.selected is None
        else [
            cell
            for cell in real_selection.cells
            if cell["candidate"] == real_selection.selected_key
        ]
    )
    pre_gate = {
        "selected_nonidentity": real_selection.selected is not None,
        "all_inner_holdout_deltas_nonnegative": bool(
            selected_cells
            and min(float(cell["delta"]) for cell in selected_cells) >= 0.0
        ),
        "all_nonpanel_test_cells_nonnegative": pretest_minimum >= 0.0,
        "every_outer_fold_pretest_angle_below_10_degrees": all(
            value < ANGLE_LIMIT_DEGREES for value in pretest_fold_angles.values()
        ),
        "nonpanel_test_macro_positive": pretest_macro > 0.0,
        "nonpanel_test_gain_above_permuted_selection_null_p95": (
            pretest_macro > null_p95["permuted"]
        ),
        "nonpanel_test_gain_above_gaussian_selection_null_p95": (
            pretest_macro > null_p95["gaussian"]
        ),
    }
    pre_gate["passed"] = all(pre_gate.values())
    selection_freeze["complete_selection_nulls"] = null_records
    selection_freeze["pretest_gate"] = pre_gate
    selection_path.write_text(json.dumps(selection_freeze, indent=2) + "\n")
    selection_hash = sha256_file(selection_path)

    result: dict[str, object] = {
        "schema_version": RESULT_SCHEMA,
        "experiment": "082-biological-target-kernel",
        "status": "pre_gate_passed" if pre_gate["passed"] else "pre_gate_failed",
        "method": {
            "formula": "arm_D + beta * ESM2_kernel_residual",
            "target_kernel": "linear kernel after fit-only feature standardization",
            "kernel_scaling": "X X^T / embedding_dimension",
            "residual_basis": "fit-only centered PCA",
            "identity_candidate": True,
            "missing_embedding_behavior": "exact arm-D identity",
        },
        "protocol": {
            "outer_donors": list(OUTER_DONORS),
            "development_folds": list(DEVELOPMENT_FOLDS),
            "required_final_folds": list(FINAL_FOLDS),
            "teacher_donors": list(TEACHER_DONORS),
            "invalid_k562_outer_reason": (
                "PerturbMatch is a reanalysis of the K562 assay"
            ),
            "training_targets": len(train_targets),
            "outer_targets": len(outer_targets),
            "split_seed": args.split_seed,
            "fit_targets": int(len(splits["fit"])),
            "tune_targets": int(len(splits["tune"])),
            "test_targets": int(len(splits["test"])),
            "complete_gallery_evaluation": True,
            "development_gallery_size": len(train_targets),
            "development_query_aggregation": (
                "fit/tune/test choose query rows only; every row ranks against "
                "the complete 600-target truth gallery"
            ),
            "outer_truth_open_policy": (
                "select and load role-specific held-panel rows only after "
                "selection freeze, both null selection runs, and exact final "
                "model-weight freeze; donor files may be authenticated and "
                "used as sources in the opposite outer fold"
            ),
            "independent_kolf": {
                "required_for_promotion": True,
                "query": "source-only three-donor frozen production arm D",
                "fit": "KOLF-specific residual weights fit on 600 non-panel rows only",
                "gallery": "complete pinned 273-target KOLF panel",
                "primary_macro_bootstrap_and_nulls": False,
            },
        },
        "leakage_audit": {
            "vcc_target_coordinates_removed": True,
            "development_target_coordinates_removed_before_dev_whitening": True,
            "development_target_coordinates_excluded_from_dev_score_axis": True,
            "development_target_coordinates_retained_in_final_panel_score_axis": True,
            "final_panel_excludes_exactly_vcc_target_coordinates": True,
            "fit_tune_test_disjoint": True,
            "embedding_standardization_fit_only": True,
            "residual_pca_fit_only": True,
            "ridge_fit_only": True,
            "outer_donor_absent_from_query_and_operator": True,
            "inner_held_teacher_absent_from_residual_fit": True,
            "vcc_panel_used_for_selection": False,
            "outer_truth_opened_before_selection": False,
            "outer_truth_opened_before_pre_gate": False,
            "nulls_repeat_selection": True,
            "development_queries_use_complete_600_target_gallery": True,
            "kolf_truth_used_for_selection_or_fitting": False,
            "kolf_excluded_from_primary_macro_bootstrap_and_null_p95": True,
            "final_panel_null_strata_are_fold_specific_and_source_only": True,
            "all_final_real_and_null_weights_serialized_before_panel_truth": True,
        },
        "embedding": {
            **embedding_table.audit,
            **verify_embedding_artifact(args.embedding_parquet, manifest),
            "train_available": int(np.sum(train_available)),
            "fold_panel_available": {
                held: int(np.sum(fold_panel_available[held])) for held in FINAL_FOLDS
            },
            "panel_available": int(np.sum(panel_available)),
        },
        "selection_freeze_sha256": selection_hash,
        "selection": real_selection.audit,
        "inner_cells": list(real_selection.cells),
        "pretest": {
            "cells": pretest_cells,
            "macro_delta": pretest_macro,
            "minimum_cell_delta": pretest_minimum,
            "nulls": {
                kind: {
                    "selection_repeated": True,
                    "replicates": len(records),
                    "p95": null_p95[kind],
                    "records": records,
                }
                for kind, records in null_records.items()
            },
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
            "outer_donors": query_metadata["valid_outer_donors"],
            "final_folds": query_metadata["final_folds"],
            "dev_score_genes": query_metadata["dev_score_genes"],
            "panel_score_genes": query_metadata["panel_score_genes"],
            "panel_query_target_counts": query_metadata["panel_query_target_counts"],
            "panel_query_target_order_sha256": query_metadata[
                "panel_query_target_order_sha256"
            ],
            "coordinate_contract": query_metadata["coordinate_contract"],
            "implementation": query_metadata["implementation"],
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    if not pre_gate["passed"]:
        (args.output_dir / "results.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        return result

    # Freeze exact final real and null model weights before selecting or
    # loading any role-specific held VCC-panel rows for evaluation. HCT116
    # and HEK293T files can already have served as source donors in the
    # opposite outer fold. No model fit occurs below this block.
    all_train_rows = np.arange(len(train_targets), dtype=np.int64)
    final_residual_bases = build_refit_residual_bases(
        train_queries,
        teacher_truth,
        train_available,
        all_train_rows,
        max_rank=max_rank,
        folds=FINAL_FOLDS,
    )
    real_predictors = build_frozen_predictors(
        train_embeddings,
        train_available,
        real_selection.selected,
        final_residual_bases,
        ridges=args.ridges,
    )
    real_freeze = {
        held: save_frozen_predictor(
            args.output_dir / "frozen_models" / "real" / f"{held}.npz",
            real_predictors[held],
        )
        for held in FINAL_FOLDS
    }
    null_freeze_records: dict[str, list[dict[str, object]]] = {
        "permuted": [],
        "gaussian": [],
    }
    by_key = {candidate.key: candidate for candidate in candidates}
    for kind in ("permuted", "gaussian"):
        for replicate in range(args.null_replicates):
            train_null, _ = make_null_embedding(
                kind,
                train_embeddings,
                train_available,
                train_strata,
                seed=stable_seed(args.null_seed, kind, replicate, "train"),
            )
            selected_key = null_selected[(kind, replicate)]
            selected_candidate = (
                None if selected_key == IDENTITY else by_key[selected_key]
            )
            frozen = build_frozen_predictors(
                train_null,
                train_available,
                selected_candidate,
                final_residual_bases,
                ridges=args.ridges,
            )
            model_records = {
                held: save_frozen_predictor(
                    args.output_dir
                    / "frozen_models"
                    / f"null-{kind}-{replicate:04d}"
                    / f"{held}.npz",
                    frozen[held],
                )
                for held in FINAL_FOLDS
            }
            null_freeze_records[kind].append(
                {
                    "replicate": replicate,
                    "selected": selected_key,
                    "models": model_records,
                }
            )
            del frozen, train_null
            enforce_memory_limit(
                args.max_rss_gib, f"{kind} final weight freeze {replicate}"
            )
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
    panel_source_counts = {
        held: np.asarray(
            np.load(
                args.query_cache / held / "panel_source_count.npy",
                allow_pickle=False,
            ),
            dtype=np.uint8,
        )
        for held in FINAL_FOLDS
    }
    cluster_root = args.output_dir / "frozen_clusters"
    cluster_root.mkdir(parents=True, exist_ok=True)
    cluster_labels_full: dict[str, np.ndarray] = {}
    cluster_audit_by_fold: dict[str, object] = {}
    cluster_records: dict[str, object] = {}
    for held in FINAL_FOLDS:
        labels, audit = source_effect_clusters(
            panel_queries[held],
            reciprocal_k=args.cluster_reciprocal_k,
            cosine_threshold=args.cluster_cosine_threshold,
        )
        path = cluster_root / f"{held}.npy"
        np.save(path, labels, allow_pickle=False)
        cluster_labels_full[held] = labels
        cluster_audit_by_fold[held] = {
            **audit,
            "fold": held,
            "source_donors": query_metadata["folds"][held]["source_donors"],
            "held_donor_present_in_cluster_source": False,
        }
        cluster_records[held] = {
            "path": str(path),
            "file_sha256": sha256_file(path),
            "labels_sha256": sha256_array(labels),
            "query_sha256": query_metadata["folds"][held]["panels"]["panel"][
                "query_sha256"
            ],
            "audit": cluster_audit_by_fold[held],
        }

    final_freeze = {
        "frozen_before_role_specific_outer_truth_evaluation": True,
        "no_post_panel_refit": True,
        "selection_freeze_sha256": selection_hash,
        "real": real_freeze,
        "nulls": null_freeze_records,
        "source_query_clusters": cluster_records,
    }
    final_freeze_path = args.output_dir / "final_model_freeze.json"
    final_freeze_path.write_text(json.dumps(final_freeze, indent=2) + "\n")
    final_freeze_hash = sha256_file(final_freeze_path)

    del real_predictors, final_residual_bases

    # Only exact persisted real/null weights may cross this boundary.
    verify_manifest_inputs(args, manifest, include_outer_truth=True)
    outer_truth_with_audit = {
        held: load_compact_panel_truth(
            args.donor_root,
            held,
            fold_panel_targets[held],
            axes.panel_score_genes,
            chunk_rows=args.truth_chunk_rows,
            complete_panel=panel,
        )
        for held in DEVELOPMENT_FOLDS
    }
    outer_truth = {
        held: values for held, (values, _audit) in outer_truth_with_audit.items()
    }
    outer_truth_audit = {
        held: audit for held, (_values, audit) in outer_truth_with_audit.items()
    }
    for held in DEVELOPMENT_FOLDS:
        gallery = manifest["held_panel_galleries"][held]
        if int(gallery["targets"]) != int(outer_truth_audit[held]["targets"]):
            raise ValueError(f"{held}: authenticated gallery count changed")
        if (
            str(gallery["target_order_sha256"])
            != outer_truth_audit[held]["target_order_sha256"]
        ):
            raise ValueError(f"{held}: authenticated gallery order changed")
    kolf = load_kolf_truth(args.kolf_root, panel, axes.panel_score_genes)
    if len(kolf["targets"]) != int(kolf_protocol["expected_targets"]):
        raise ValueError("KOLF target coverage differs from the frozen protocol")
    kolf_gallery = manifest["held_panel_galleries"][KOLF_FOLD]
    if int(kolf_gallery["targets"]) != len(kolf["targets"]):
        raise ValueError("KOLF complete-gallery count changed")
    if str(kolf_gallery["target_order_sha256"]) != sha256_lines(kolf["targets"]):
        raise ValueError("KOLF complete-gallery target order changed")

    all_panel_columns = np.arange(len(axes.panel_score_genes), dtype=np.int64)
    fold_specs: dict[str, dict[str, object]] = {
        held: {
            "query": panel_queries[held],
            "truth": outer_truth[held],
            "embeddings": fold_panel_embeddings[held],
            "available": fold_panel_available[held],
            "query_rows": np.arange(len(fold_panel_targets[held]), dtype=np.int64),
            "positive_columns": np.arange(
                len(fold_panel_targets[held]), dtype=np.int64
            ),
            "query_columns": all_panel_columns,
            "truth_columns": all_panel_columns,
            "targets": fold_panel_targets[held],
            "source_count": panel_source_counts[held],
            "complete_count": len(query_metadata["folds"][held]["source_donors"]),
            "clusters": cluster_labels_full[held],
        }
        for held in DEVELOPMENT_FOLDS
    }
    fold_specs[KOLF_FOLD] = {
        "query": panel_queries[KOLF_FOLD],
        "truth": kolf["values"],
        "embeddings": panel_embeddings,
        "available": panel_available,
        "query_rows": kolf["query_rows"],
        "positive_columns": np.arange(len(kolf["targets"]), dtype=np.int64),
        "query_columns": kolf["query_columns"],
        "truth_columns": np.arange(len(kolf["score_genes"]), dtype=np.int64),
        "targets": kolf["targets"],
        "source_count": panel_source_counts[KOLF_FOLD][kolf["query_rows"]],
        "complete_count": len(query_metadata["folds"][KOLF_FOLD]["source_donors"]),
        "clusters": cluster_labels_full[KOLF_FOLD][kolf["query_rows"]],
    }

    loaded_real_predictors = {
        held: load_frozen_predictor(
            args.output_dir / "frozen_models" / "real" / f"{held}.npz",
            real_freeze[held],
        )
        for held in FINAL_FOLDS
    }
    real_outer: dict[str, dict[str, object]] = {}
    for held in FINAL_FOLDS:
        spec = fold_specs[held]
        real_outer[held] = evaluate_frozen_fold(
            loaded_real_predictors[held],
            spec["query"],
            spec["truth"],
            spec["embeddings"],
            spec["available"],
            spec["query_rows"],
            spec["positive_columns"],
            spec["query_columns"],
            spec["truth_columns"],
            residual_positions=axes.dev_to_panel_positions,
        )
    deltas_by_fold = {
        held: np.asarray(real_outer[held]["per_target_delta"]) for held in FINAL_FOLDS
    }
    primary_deltas_by_fold = {held: deltas_by_fold[held] for held in OUTER_DONORS}
    macro_delta = float(
        np.mean([np.mean(primary_deltas_by_fold[held]) for held in OUTER_DONORS])
    )
    bootstrap = outer_fold_cluster_bootstrap(
        primary_deltas_by_fold,
        {
            held: np.asarray(fold_specs[held]["clusters"], dtype=np.int64)
            for held in OUTER_DONORS
        },
        replicates=args.bootstrap_replicates,
        seed=args.bootstrap_seed,
    )
    dependence_by_fold = {
        held: cluster_dependence_diagnostics(
            deltas_by_fold[held],
            np.asarray(fold_specs[held]["clusters"], dtype=np.int64),
        )
        for held in FINAL_FOLDS
    }
    panel_null_inputs: dict[str, dict[str, object]] = {}
    for held in FINAL_FOLDS:
        embeddings = fold_panel_embeddings[held]
        available = fold_panel_available[held]
        fold_partition = np.asarray(["panel"] * len(embeddings))
        strata = make_source_strata(
            {held: panel_queries[held]},
            {held: panel_source_counts[held]},
            fold_partition,
            norm_bins=args.null_norm_bins,
        )
        panel_null_inputs[held] = {
            "embeddings": embeddings,
            "available": available,
            "strata": strata,
            "strata_sha256": sha256_lines(strata.tolist()),
        }
    outer_null_records: dict[str, list[dict[str, object]]] = {
        "permuted": [],
        "gaussian": [],
    }
    for kind in ("permuted", "gaussian"):
        for replicate in range(args.null_replicates):
            record = null_freeze_records[kind][replicate]
            fold_deltas: dict[str, float] = {}
            fold_embedding_audits: dict[str, object] = {}
            for held in FINAL_FOLDS:
                null_input = panel_null_inputs[held]
                panel_null, panel_null_audit = make_null_embedding(
                    kind,
                    np.asarray(null_input["embeddings"]),
                    np.asarray(null_input["available"], dtype=bool),
                    np.asarray(null_input["strata"], dtype=str),
                    seed=stable_seed(
                        args.null_seed,
                        kind,
                        replicate,
                        "panel",
                        held,
                    ),
                )
                model_record = record["models"][held]
                null_predictor = load_frozen_predictor(
                    args.output_dir
                    / "frozen_models"
                    / f"null-{kind}-{replicate:04d}"
                    / f"{held}.npz",
                    model_record,
                )
                spec = fold_specs[held]
                evaluated = evaluate_frozen_fold(
                    null_predictor,
                    spec["query"],
                    spec["truth"],
                    panel_null,
                    np.asarray(null_input["available"], dtype=bool),
                    spec["query_rows"],
                    spec["positive_columns"],
                    spec["query_columns"],
                    spec["truth_columns"],
                    residual_positions=axes.dev_to_panel_positions,
                )
                fold_deltas[held] = float(evaluated["delta"])
                fold_embedding_audits[held] = {
                    **panel_null_audit,
                    "source_fold": held,
                    "source_only_strata": True,
                    "strata_sha256": null_input["strata_sha256"],
                }
            outer_null_records[kind].append(
                {
                    "replicate": replicate,
                    "selected": null_selected[(kind, replicate)],
                    "fold_deltas": fold_deltas,
                    "primary_macro_delta": float(
                        np.mean([fold_deltas[held] for held in OUTER_DONORS])
                    ),
                    "kolf_delta": fold_deltas[KOLF_FOLD],
                    "embedding_audit_by_fold": fold_embedding_audits,
                }
            )
            enforce_memory_limit(args.max_rss_gib, f"{kind} outer null {replicate}")
    outer_null_p95 = {
        kind: float(np.quantile([row["primary_macro_delta"] for row in records], 0.95))
        for kind, records in outer_null_records.items()
    }
    kolf_null_p95 = {
        kind: float(np.quantile([row["kolf_delta"] for row in records], 0.95))
        for kind, records in outer_null_records.items()
    }
    fold_deltas = {held: float(real_outer[held]["delta"]) for held in FINAL_FOLDS}
    primary_fold_deltas = {held: fold_deltas[held] for held in OUTER_DONORS}
    baseline_nonzero = all(
        np.all(np.asarray(real_outer[held]["baseline_norms"]) > 1e-12)
        for held in FINAL_FOLDS
    )
    candidate_nonzero = all(
        np.all(np.asarray(real_outer[held]["adapted_norms"]) > 1e-12)
        for held in FINAL_FOLDS
    )
    fold_gates: dict[str, object] = {}
    for held in FINAL_FOLDS:
        source_count = np.asarray(fold_specs[held]["source_count"], dtype=np.uint8)
        complete_count = int(fold_specs[held]["complete_count"])
        strata = stratum_summaries(
            np.asarray(real_outer[held]["per_target_delta"]),
            source_count,
            complete_count=complete_count,
        )
        complete = strata["complete_source"]
        partial = strata["partial_source"]
        fold_gates[held] = {
            "delta_nonnegative": float(real_outer[held]["delta"]) >= 0.0,
            "median_angular_displacement_below_10_degrees": (
                float(real_outer[held]["median_angle_degrees"]) < ANGLE_LIMIT_DEGREES
            ),
            "complete_source_present_and_nonnegative": (
                int(complete["count"]) > 0
                and complete["delta"] is not None
                and float(complete["delta"]) >= 0.0
            ),
            "partial_source_present_and_nonnegative": (
                int(partial["count"]) > 0
                and partial["delta"] is not None
                and float(partial["delta"]) >= 0.0
            ),
            "no_obvious_single_cluster_dependence": dependence_by_fold[held][
                "passes_no_obvious_single_family_dependence"
            ],
        }
        fold_gates[held]["passed"] = all(fold_gates[held].values())
    promotion = {
        "predeclared": True,
        "nonpanel_pre_gate_passed": True,
        "selected_nonidentity": real_selection.selected is not None,
        "all_valid_outer_folds_nonnegative": min(primary_fold_deltas.values()) >= 0.0,
        "kolf_independent_confirmation_included": KOLF_FOLD in fold_deltas,
        "kolf_independent_confirmation_nonnegative": fold_deltas[KOLF_FOLD] >= 0.0,
        "all_primary_fold_angular_and_strata_gates_pass": all(
            bool(fold_gates[held]["passed"]) for held in OUTER_DONORS
        ),
        "kolf_independent_confirmation_gate_passed": bool(
            fold_gates[KOLF_FOLD]["passed"]
        ),
        "outer_macro_raw_pds_delta_at_least_0_003": macro_delta >= 0.003,
        "cluster_bootstrap_lower_95_above_zero": bootstrap["lower_95"] > 0.0,
        "observed_gain_above_permuted_selection_null_p95": (
            macro_delta > outer_null_p95["permuted"]
        ),
        "observed_gain_above_gaussian_selection_null_p95": (
            macro_delta > outer_null_p95["gaussian"]
        ),
        "target_coverage_unchanged": baseline_nonzero and candidate_nonzero,
        "exact_final_weights_frozen_before_panel_truth": True,
    }
    promotion["passed"] = all(
        bool(value) for key, value in promotion.items() if key != "predeclared"
    )
    per_target_rows: list[dict[str, object]] = []
    fold_summaries: dict[str, object] = {}
    source_counts_by_fold: dict[str, np.ndarray] = {}
    complete_counts = {
        held: int(fold_specs[held]["complete_count"]) for held in FINAL_FOLDS
    }
    for held in FINAL_FOLDS:
        evaluated = real_outer[held]
        source_count_values = np.asarray(
            fold_specs[held]["source_count"], dtype=np.uint8
        )
        source_counts_by_fold[held] = source_count_values
        strata = stratum_summaries(
            np.asarray(evaluated["per_target_delta"]),
            source_count_values,
            complete_count=complete_counts[held],
        )
        fold_summaries[held] = {
            **summarize_evaluation(evaluated),
            "source_donors": query_metadata["folds"][held]["source_donors"],
            "source_strata": strata,
            "held_truth": (
                kolf["audit"] if held == KOLF_FOLD else outer_truth_audit[held]
            ),
        }
        targets = tuple(fold_specs[held]["targets"])
        available = np.asarray(fold_specs[held]["available"], dtype=bool)[
            np.asarray(fold_specs[held]["query_rows"], dtype=np.int64)
        ]
        clusters = np.asarray(fold_specs[held]["clusters"], dtype=np.int64)
        for index, target in enumerate(targets):
            source_count = int(source_count_values[index])
            per_target_rows.append(
                {
                    "outer_donor": held,
                    "target_gene": target,
                    "source_cluster": int(clusters[index]),
                    "source_count": source_count,
                    "source_stratum": (
                        "complete_source"
                        if source_count == complete_counts[held]
                        else "partial_source"
                    ),
                    "direct_row": source_count >= 1,
                    "missing_source": source_count == 0,
                    "embedding_available": bool(available[index]),
                    "baseline_pds": float(evaluated["per_target_baseline"][index]),
                    "adapted_pds": float(evaluated["per_target_adapted"][index]),
                    "delta": float(evaluated["per_target_delta"][index]),
                }
            )
    result["status"] = "completed"
    result["leakage_audit"]["outer_truth_opened_after_pre_gate"] = True
    result["outer"] = {
        "evaluated": True,
        "complete_gallery": True,
        "selection_freeze_sha256_before_truth_open": selection_hash,
        "final_model_freeze_sha256_before_truth_open": final_freeze_hash,
        "no_post_panel_refit": True,
        "folds": fold_summaries,
        "macro_baseline_pds": float(
            np.mean([real_outer[held]["baseline"] for held in OUTER_DONORS])
        ),
        "macro_adapted_pds": float(
            np.mean([real_outer[held]["adapted"] for held in OUTER_DONORS])
        ),
        "macro_delta": macro_delta,
        "primary_folds": list(OUTER_DONORS),
        "primary_fold_deltas": primary_fold_deltas,
        "fold_deltas": fold_deltas,
        "aggregate_source_strata": {
            label: dict(
                zip(
                    ("delta", "count"),
                    aggregate_stratum_delta(
                        {held: real_outer[held] for held in OUTER_DONORS},
                        {held: source_counts_by_fold[held] for held in OUTER_DONORS},
                        {held: complete_counts[held] for held in OUTER_DONORS},
                        label,
                    ),
                    strict=True,
                )
            )
            for label in (
                "complete_source",
                "partial_source",
                "direct_row",
                "missing_source",
            )
        },
        "fold_gates": fold_gates,
        "target_clusters": cluster_audit_by_fold,
        "cluster_bootstrap": bootstrap,
        "cluster_dependence": {held: dependence_by_fold[held] for held in OUTER_DONORS},
        "kolf_confirmation": {
            "included_in_primary_macro": False,
            "included_in_primary_bootstrap": False,
            "included_in_primary_null_p95": False,
            "delta": fold_deltas[KOLF_FOLD],
            "fold": fold_summaries[KOLF_FOLD],
            "gate": fold_gates[KOLF_FOLD],
            "cluster_dependence": dependence_by_fold[KOLF_FOLD],
            "truth": kolf["audit"],
            "null_p95": kolf_null_p95,
        },
        "nulls": {
            kind: {
                "selection_repeated": True,
                "final_weights_frozen_before_panel_truth": True,
                "primary_folds": list(OUTER_DONORS),
                "kolf_reported_separately": True,
                "fold_specific_source_only_strata": True,
                "replicates": len(records),
                "p95": outer_null_p95[kind],
                "records": records,
            }
            for kind, records in outer_null_records.items()
        },
        "promotion_gate": promotion,
    }
    result["peak_rss_gib"] = peak_rss_gib()
    (args.output_dir / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    write_csv(args.output_dir / "per_target_pds.csv", per_target_rows)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-cache", type=Path, required=True)
    parser.add_argument("--donor-root", type=Path, required=True)
    parser.add_argument("--kolf-root", type=Path, required=True)
    parser.add_argument("--query-cache", type=Path, required=True)
    parser.add_argument("--embedding-parquet", type=Path, required=True)
    parser.add_argument("--artifact-manifest", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=HERE / "output")
    parser.add_argument("--ranks", nargs="+", type=int, default=[16, 32, 64])
    parser.add_argument(
        "--ridges",
        nargs="+",
        type=float,
        default=[1e-3, 1e-2, 1e-1, 1.0, 10.0],
    )
    parser.add_argument(
        "--betas", nargs="+", type=float, default=[0.02, 0.05, 0.10, 0.20]
    )
    parser.add_argument("--split-seed", type=int, default=8202026)
    parser.add_argument("--null-seed", type=int, default=8212026)
    parser.add_argument("--null-replicates", type=int, default=100)
    parser.add_argument("--null-norm-bins", type=int, default=10)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--bootstrap-seed", type=int, default=8222026)
    parser.add_argument("--cluster-reciprocal-k", type=int, default=5)
    parser.add_argument("--cluster-cosine-threshold", type=float, default=0.80)
    parser.add_argument("--expected-embedding-rows", type=int, default=19_790)
    parser.add_argument("--expected-embedding-dim", type=int, default=5_120)
    parser.add_argument("--embedding-column-chunk", type=int, default=128)
    parser.add_argument("--truth-chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=24.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    if sorted(args.ranks) != [16, 32, 64]:
        raise SystemExit("production ranks must be exactly 16 32 64")
    if sorted(args.ridges) != [1e-3, 1e-2, 1e-1, 1.0, 10.0]:
        raise SystemExit("production ridge grid must be exactly 1e-3..10")
    if sorted(args.betas) != [0.02, 0.05, 0.10, 0.20]:
        raise SystemExit("production non-identity betas must be exactly 0.02..0.20")
    if args.null_replicates < 20:
        raise SystemExit("--null-replicates must be at least 20")
    if args.bootstrap_replicates < 100:
        raise SystemExit("--bootstrap-replicates must be at least 100")
    if args.null_norm_bins < 2:
        raise SystemExit("--null-norm-bins must be at least 2")
    if (
        args.embedding_column_chunk < 1
        or args.truth_chunk_rows < 1
        or args.max_rss_gib <= 0
    ):
        raise SystemExit("chunk and memory controls must be positive")
    manifest = json.loads(args.artifact_manifest.read_text())
    if int(manifest.get("schema_version", -1)) != 1:
        raise SystemExit("unsupported artifact manifest")
    # Deliberately defer hashing/opening *_raw_outer.npy until the full
    # repeated-null non-panel pre-gate succeeds.
    verify_manifest_inputs(args, manifest, include_outer_truth=False)
    result = run_experiment(args, manifest)
    print(json.dumps(result["pretest"]["gate"], indent=2), flush=True)
    if result["outer"]["evaluated"]:
        print(json.dumps(result["outer"]["promotion_gate"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
