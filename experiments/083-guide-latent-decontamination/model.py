#!/usr/bin/env python3
"""Pure numerical primitives for guide-latent decontamination.

The real-data runner keeps I/O and frozen arm-D reconstruction separate from
this module so that every learning rule can be exercised on tiny synthetic
arrays.  No function here knows VCC-panel outcomes.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


SPLIT_SEED = 83_2026
FIT_FRACTION = 0.60
TUNE_FRACTION = 0.20
EPS = 1e-12


def stable_seed(*parts: object) -> int:
    payload = "\0".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def deterministic_target_split(
    targets: Iterable[str],
    *,
    seed: int = SPLIT_SEED,
    fit_fraction: float = FIT_FRACTION,
    tune_fraction: float = TUNE_FRACTION,
) -> dict[str, tuple[str, ...]]:
    """Hash-order targets into exact, disjoint fit/tune/test partitions."""
    unique = sorted(set(map(str, targets)))
    if len(unique) < 5:
        raise ValueError("need at least five targets for fit/tune/test")
    if not (0 < fit_fraction < 1 and 0 < tune_fraction < 1):
        raise ValueError("split fractions must be in (0, 1)")
    if fit_fraction + tune_fraction >= 1:
        raise ValueError("fit+tune fraction must be below one")
    ordered = sorted(
        unique,
        key=lambda value: (
            hashlib.sha256(f"{seed}\0{value}".encode()).digest(),
            value,
        ),
    )
    n = len(ordered)
    n_fit = max(1, int(math.floor(n * fit_fraction)))
    n_tune = max(1, int(math.floor(n * tune_fraction)))
    if n_fit + n_tune >= n:
        n_tune = 1
        n_fit = n - 2
    return {
        "fit": tuple(sorted(ordered[:n_fit])),
        "tune": tuple(sorted(ordered[n_fit : n_fit + n_tune])),
        "test": tuple(sorted(ordered[n_fit + n_tune :])),
    }


def assert_split_contract(
    splits: Mapping[str, Sequence[str]], panel_targets: Iterable[str]
) -> None:
    expected = {"fit", "tune", "test"}
    if set(splits) != expected:
        raise ValueError(f"split keys {set(splits)} != {expected}")
    sets = {name: set(map(str, values)) for name, values in splits.items()}
    if any(sets[a] & sets[b] for a in expected for b in expected if a < b):
        raise ValueError("fit/tune/test targets overlap")
    panel = set(map(str, panel_targets))
    leaked = panel & set().union(*sets.values())
    if leaked:
        raise ValueError(f"VCC-panel targets entered model split: {sorted(leaked)[:5]}")


def excluded_coordinate_mask(
    genes: Sequence[str], panel_targets: Iterable[str]
) -> np.ndarray:
    panel = set(map(str, panel_targets))
    mask = np.asarray([str(gene) not in panel for gene in genes], dtype=bool)
    if not np.any(mask):
        raise ValueError("all response coordinates were excluded")
    return mask


def target_rows(targets: Sequence[str]) -> dict[str, np.ndarray]:
    grouped: dict[str, list[int]] = {}
    for index, target in enumerate(map(str, targets)):
        grouped.setdefault(target, []).append(index)
    return {
        target: np.asarray(indices, dtype=np.int64)
        for target, indices in grouped.items()
    }


def identifiable_targets(
    targets: Sequence[str],
    identifiable: np.ndarray,
    *,
    minimum_guides: int = 2,
) -> tuple[str, ...]:
    identifiable = np.asarray(identifiable, dtype=bool)
    if len(targets) != len(identifiable):
        raise ValueError("target/identifiable length mismatch")
    grouped = target_rows(targets)
    return tuple(
        sorted(
            target
            for target, rows in grouped.items()
            if int(np.sum(identifiable[rows])) >= minimum_guides
        )
    )


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("non-positive weight sum")
    return np.sum(values * weights[:, None], axis=0) / total


def robust_center(
    values: np.ndarray,
    cell_counts: np.ndarray,
    *,
    weight_power: float,
    huber_c: float,
    iterations: int = 4,
) -> np.ndarray:
    """Robust vector center without deleting any guide observation."""
    values = np.asarray(values, dtype=np.float64)
    counts = np.asarray(cell_counts, dtype=np.float64)
    if values.ndim != 2 or counts.shape != (len(values),):
        raise ValueError("invalid robust-center shapes")
    if len(values) == 0:
        raise ValueError("cannot center zero guides")
    if np.any(counts <= 0) or not np.isfinite(values).all():
        raise ValueError("non-positive counts or non-finite guide effects")
    base = np.power(counts / np.median(counts), weight_power)
    center = _weighted_mean(values, base)
    for _ in range(iterations):
        distances = np.linalg.norm(values - center, axis=1)
        positive = distances[distances > EPS]
        if not len(positive):
            break
        scale = float(np.median(positive))
        robust = np.minimum(1.0, huber_c * max(scale, EPS) / np.maximum(distances, EPS))
        center = _weighted_mean(values, base * robust)
    return center


def _center_values_by_labels(
    values: np.ndarray,
    labels: Sequence[str],
) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    local_targets = np.asarray(labels, dtype=object)
    if len(values) != len(local_targets):
        raise ValueError("value/label row mismatch")
    output = np.empty_like(values, dtype=np.float64)
    for target in sorted(set(map(str, local_targets.tolist()))):
        local = np.flatnonzero(local_targets == target)
        block = values[local]
        output[local] = block - block.mean(axis=0, keepdims=True)
    return output


def randomized_right_basis(
    matrix: np.ndarray,
    rank: int,
    *,
    seed: int,
    oversample: int = 8,
    power_iterations: int = 1,
) -> np.ndarray:
    """Return rank orthonormal right singular vectors without a full SVD."""
    values = np.asarray(matrix, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError("matrix must be two-dimensional")
    max_rank = min(values.shape)
    if not 0 < rank <= max_rank:
        raise ValueError(f"rank {rank} outside 1..{max_rank}")
    width = min(max_rank, rank + max(0, oversample))
    rng = np.random.default_rng(seed)
    omega = rng.normal(size=(values.shape[1], width))
    projected = values @ omega
    for _ in range(power_iterations):
        projected = values @ (values.T @ projected)
    q, _ = np.linalg.qr(projected, mode="reduced")
    small = q.T @ values
    _, _, vt = np.linalg.svd(small, full_matrices=False)
    basis = np.asarray(vt[:rank], dtype=np.float64)
    gram = basis @ basis.T
    np.testing.assert_allclose(gram, np.eye(rank), atol=1e-8, rtol=1e-8)
    return basis


@dataclass(frozen=True, order=True)
class Candidate:
    rank: int
    ridge: float
    eta: float
    weight_power: float
    huber_c: float

    @property
    def is_identity(self) -> bool:
        return self.rank == 0 or self.eta == 0


@dataclass(frozen=True)
class GuideLatentModel:
    candidate: Candidate
    feature_mean: np.ndarray
    feature_scale: np.ndarray
    feature_keep: np.ndarray
    response_mask: np.ndarray
    basis: np.ndarray
    coefficients: np.ndarray

    def standardized_features(self, features: np.ndarray) -> np.ndarray:
        values = np.asarray(features, dtype=np.float64)
        kept = values[:, self.feature_keep]
        return (kept - self.feature_mean) / self.feature_scale

    def nuisance(self, features: np.ndarray) -> np.ndarray:
        if self.candidate.is_identity:
            return np.zeros(
                (len(features), int(np.sum(self.response_mask))), dtype=np.float64
            )
        return self.standardized_features(features) @ self.coefficients @ self.basis


def candidate_grid() -> tuple[Candidate, ...]:
    values = [Candidate(0, 0.0, 0.0, 0.0, 2.0)]
    for rank in (8, 16, 32):
        for ridge in (0.1, 1.0, 10.0, 100.0):
            for eta in (0.25, 0.50, 0.75, 1.0):
                for weight_power in (0.0, 0.5):
                    values.append(Candidate(rank, ridge, eta, weight_power, 2.0))
    return tuple(values)


def fit_latent_model(
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    identifiable: np.ndarray,
    fit_targets: Iterable[str],
    response_mask: np.ndarray,
    candidate: Candidate,
    *,
    seed: int = SPLIT_SEED,
    fixed_basis: np.ndarray | None = None,
) -> GuideLatentModel:
    """Fit the low-rank guide nuisance map on fit targets only."""
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    identifiable = np.asarray(identifiable, dtype=bool)
    response_mask = np.asarray(response_mask, dtype=bool)
    if effects.ndim != 2 or features.ndim != 2:
        raise ValueError("effects/features must be matrices")
    if len(effects) != len(features) or len(effects) != len(guide_targets):
        raise ValueError("guide arrays have inconsistent rows")
    if effects.shape[1] != len(response_mask):
        raise ValueError("response mask width mismatch")
    fit_set = set(map(str, fit_targets))
    eligible_targets = (
        set(identifiable_targets(guide_targets, identifiable, minimum_guides=2))
        & fit_set
    )
    rows = np.asarray(
        [
            i
            for i, target in enumerate(map(str, guide_targets))
            if identifiable[i] and target in eligible_targets
        ],
        dtype=np.int64,
    )
    if len(rows) < 4:
        raise ValueError("too few identifiable fit-guide rows")

    fit_features = np.asarray(features[rows], dtype=np.float64)
    feature_mean_all = fit_features.mean(axis=0)
    feature_scale_all = fit_features.std(axis=0)
    feature_keep = np.isfinite(feature_scale_all) & (feature_scale_all > 1e-8)
    if not np.any(feature_keep):
        raise ValueError("no varying guide features")
    feature_mean = feature_mean_all[feature_keep]
    feature_scale = feature_scale_all[feature_keep]
    standardized = (fit_features[:, feature_keep] - feature_mean) / feature_scale

    row_labels = tuple(str(guide_targets[index]) for index in rows)
    centered_x = _center_values_by_labels(standardized, row_labels)
    response = np.asarray(
        effects[np.ix_(rows, np.flatnonzero(response_mask))], dtype=np.float64
    )
    centered_y = _center_values_by_labels(response, row_labels)

    if candidate.is_identity:
        return GuideLatentModel(
            candidate=candidate,
            feature_mean=feature_mean,
            feature_scale=feature_scale,
            feature_keep=feature_keep,
            response_mask=response_mask,
            basis=np.zeros((0, response.shape[1]), dtype=np.float64),
            coefficients=np.zeros((int(feature_keep.sum()), 0), dtype=np.float64),
        )

    if fixed_basis is None:
        basis = randomized_right_basis(
            centered_y,
            candidate.rank,
            seed=stable_seed(seed, "basis", candidate.rank),
        )
    else:
        basis = np.asarray(fixed_basis[: candidate.rank], dtype=np.float64)
        if basis.shape != (candidate.rank, response.shape[1]):
            raise ValueError("fixed basis shape mismatch")
    scores = centered_y @ basis.T
    gram = centered_x.T @ centered_x
    gram.flat[:: len(gram) + 1] += candidate.ridge
    coefficients = np.linalg.solve(gram, centered_x.T @ scores)
    return GuideLatentModel(
        candidate=candidate,
        feature_mean=feature_mean,
        feature_scale=feature_scale,
        feature_keep=feature_keep,
        response_mask=response_mask,
        basis=basis,
        coefficients=coefficients,
    )


def _held_guide_index(target: str, rows: np.ndarray, seed: int) -> int:
    return int(rows[stable_seed(seed, "held-guide", target) % len(rows)])


def held_guide_predictions(
    model: GuideLatentModel,
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    cell_counts: np.ndarray,
    identifiable: np.ndarray,
    evaluation_targets: Iterable[str],
    *,
    seed: int = SPLIT_SEED,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    """Predict one deterministic held guide per target.

    Returns truth, identity prediction, candidate prediction, and target order.
    """
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    counts = np.asarray(cell_counts)
    identifiable = np.asarray(identifiable, dtype=bool)
    grouped = target_rows(guide_targets)
    truth_rows: list[np.ndarray] = []
    baseline_rows: list[np.ndarray] = []
    candidate_rows: list[np.ndarray] = []
    labels: list[str] = []
    mask = model.response_mask
    for target in sorted(set(map(str, evaluation_targets))):
        rows = grouped.get(target)
        if rows is None:
            continue
        rows = rows[identifiable[rows]]
        if len(rows) < 2:
            continue
        held = _held_guide_index(target, rows, seed)
        train = rows[rows != held]
        baseline = robust_center(
            effects[np.ix_(train, np.flatnonzero(mask))],
            counts[train],
            weight_power=model.candidate.weight_power,
            huber_c=model.candidate.huber_c,
        )
        nuisance = model.nuisance(features[np.concatenate((train, [held]))])
        train_nuisance = nuisance[:-1]
        held_nuisance = nuisance[-1]
        theta = robust_center(
            effects[np.ix_(train, np.flatnonzero(mask))] - train_nuisance,
            counts[train],
            weight_power=model.candidate.weight_power,
            huber_c=model.candidate.huber_c,
        )
        modeled = theta + held_nuisance
        candidate = baseline + model.candidate.eta * (modeled - baseline)
        truth_rows.append(np.asarray(effects[held, mask], dtype=np.float64))
        baseline_rows.append(baseline)
        candidate_rows.append(candidate)
        labels.append(target)
    if not labels:
        width = int(np.sum(mask))
        empty = np.empty((0, width), dtype=np.float64)
        return empty, empty.copy(), empty.copy(), ()
    return (
        np.stack(truth_rows),
        np.stack(baseline_rows),
        np.stack(candidate_rows),
        tuple(labels),
    )


def row_cosines(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    left = np.asarray(a, dtype=np.float64)
    right = np.asarray(b, dtype=np.float64)
    numerator = np.sum(left * right, axis=1)
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return numerator / np.maximum(denominator, EPS)


def guide_prediction_metrics(
    truth: np.ndarray, baseline: np.ndarray, candidate: np.ndarray
) -> dict[str, float]:
    if len(truth) == 0:
        return {
            "targets": 0,
            "baseline_cosine": float("nan"),
            "candidate_cosine": float("nan"),
            "cosine_delta": float("nan"),
            "baseline_nmse": float("nan"),
            "candidate_nmse": float("nan"),
            "nmse_delta": float("nan"),
        }
    baseline_cos = row_cosines(truth, baseline)
    candidate_cos = row_cosines(truth, candidate)
    denominator = np.maximum(np.sum(truth * truth, axis=1), EPS)
    baseline_nmse = np.sum((baseline - truth) ** 2, axis=1) / denominator
    candidate_nmse = np.sum((candidate - truth) ** 2, axis=1) / denominator
    return {
        "targets": int(len(truth)),
        "baseline_cosine": float(baseline_cos.mean()),
        "candidate_cosine": float(candidate_cos.mean()),
        "cosine_delta": float((candidate_cos - baseline_cos).mean()),
        "baseline_nmse": float(baseline_nmse.mean()),
        "candidate_nmse": float(candidate_nmse.mean()),
        "nmse_delta": float((candidate_nmse - baseline_nmse).mean()),
    }


def infer_cleaned_target(
    model: GuideLatentModel,
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    cell_counts: np.ndarray,
    identifiable: np.ndarray,
    pooled_target_effect: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    """Infer one cleaned target row while preserving the exact fallback."""
    effects = np.asarray(guide_effects)
    features = np.asarray(guide_features)
    counts = np.asarray(cell_counts)
    identifiable = np.asarray(identifiable, dtype=bool)
    pooled = np.asarray(pooled_target_effect, dtype=np.float64)
    usable = np.flatnonzero(identifiable)
    if model.candidate.is_identity or len(usable) < 2:
        return pooled.copy(), {
            "modified": False,
            "reason": "identity"
            if model.candidate.is_identity
            else "fewer_than_two_identifiable_guides",
            "guide_count": int(len(effects)),
            "identifiable_guides": int(len(usable)),
            "correction_norm": 0.0,
        }
    nuisance = model.nuisance(features[usable])
    adjusted = (
        np.asarray(
            effects[np.ix_(usable, np.flatnonzero(model.response_mask))],
            dtype=np.float64,
        )
        - nuisance
    )
    theta = robust_center(
        adjusted,
        counts[usable],
        weight_power=model.candidate.weight_power,
        huber_c=model.candidate.huber_c,
    )
    cleaned = pooled.copy()
    correction = theta - pooled[model.response_mask]
    cleaned[model.response_mask] += model.candidate.eta * correction
    return cleaned, {
        "modified": True,
        "reason": "latent_decontamination",
        "guide_count": int(len(effects)),
        "identifiable_guides": int(len(usable)),
        "correction_norm": float(np.linalg.norm(model.candidate.eta * correction)),
        "raw_norm": float(np.linalg.norm(pooled[model.response_mask])),
        "norm_ratio": float(
            np.linalg.norm(cleaned[model.response_mask])
            / max(np.linalg.norm(pooled[model.response_mask]), EPS)
        ),
    }


def pds_per_target(predictions: np.ndarray, truth: np.ndarray) -> np.ndarray:
    pred = np.asarray(predictions, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    if pred.shape != observed.shape or pred.ndim != 2:
        raise ValueError("PDS matrices must have the same 2D shape")
    pred /= np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), EPS)
    observed /= np.maximum(np.linalg.norm(observed, axis=1, keepdims=True), EPS)
    scores = pred @ observed.T
    own = np.diag(scores)
    ahead = np.sum(scores > own[:, None], axis=1)
    ties = np.sum(scores == own[:, None], axis=1) - 1
    return 1.0 - (ahead + 0.5 * ties) / max(len(pred) - 1, 1)


def source_effect_clusters(
    source_queries: np.ndarray,
    *,
    k: int = 8,
) -> np.ndarray:
    """Reciprocal-kNN connected components using source-derived queries only."""
    values = np.asarray(source_queries, dtype=np.float64)
    if values.ndim != 2 or len(values) < 2:
        raise ValueError("need at least two source queries")
    normalized = values / np.maximum(np.linalg.norm(values, axis=1, keepdims=True), EPS)
    similarities = normalized @ normalized.T
    np.fill_diagonal(similarities, -np.inf)
    width = min(max(1, k), len(values) - 1)
    neighbors = np.argpartition(-similarities, width - 1, axis=1)[:, :width]
    parent = np.arange(len(values), dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    neighbor_sets = [set(map(int, row.tolist())) for row in neighbors]
    for left, row in enumerate(neighbors):
        for right in map(int, row.tolist()):
            if left in neighbor_sets[right]:
                union(left, right)
    roots = [find(i) for i in range(len(values))]
    remap = {root: index for index, root in enumerate(sorted(set(roots)))}
    return np.asarray([remap[root] for root in roots], dtype=np.int64)


def clustered_bootstrap(
    target_deltas: np.ndarray,
    cluster_ids: np.ndarray,
    *,
    repetitions: int = 20_000,
    seed: int = SPLIT_SEED,
) -> dict[str, float]:
    deltas = np.asarray(target_deltas, dtype=np.float64)
    clusters = np.asarray(cluster_ids, dtype=np.int64)
    if deltas.shape != clusters.shape or len(deltas) == 0:
        raise ValueError("invalid clustered bootstrap inputs")
    unique = np.unique(clusters)
    members = [np.flatnonzero(clusters == cluster) for cluster in unique]
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=np.float64)
    for i in range(repetitions):
        sampled = rng.integers(0, len(members), size=len(members))
        rows = np.concatenate([members[index] for index in sampled])
        draws[i] = float(np.mean(deltas[rows]))
    return {
        "repetitions": int(repetitions),
        "clusters": int(len(unique)),
        "mean": float(np.mean(deltas)),
        "lower_95": float(np.percentile(draws, 2.5)),
        "upper_95": float(np.percentile(draws, 97.5)),
    }


def stratified_feature_permutation(
    features: np.ndarray,
    guide_targets: Sequence[str],
    target_norms: Mapping[str, float],
    *,
    seed: int,
    bundle: bool,
) -> np.ndarray:
    """Permute guide features within guide-count and response-norm strata.

    ``bundle=True`` moves complete same-size guide bundles between targets,
    implementing the guide-to-target assignment null. ``bundle=False``
    independently permutes rows inside each stratum, implementing the
    seed-feature null.
    """
    values = np.asarray(features)
    grouped = target_rows(guide_targets)
    ordered_targets = sorted(grouped)
    norms = np.asarray([float(target_norms[target]) for target in ordered_targets])
    if len(norms) < 2:
        return values.copy()
    quantiles = np.quantile(norms, np.linspace(0, 1, 11))
    norm_bin = np.searchsorted(quantiles[1:-1], norms, side="right")
    strata: dict[tuple[int, int], list[str]] = {}
    for target, decile in zip(ordered_targets, norm_bin, strict=True):
        strata.setdefault((len(grouped[target]), int(decile)), []).append(target)
    rng = np.random.default_rng(seed)
    output = values.copy()
    for targets_in_stratum in strata.values():
        if len(targets_in_stratum) < 2:
            continue
        if bundle:
            donors = list(targets_in_stratum)
            rng.shuffle(donors)
            for recipient, donor in zip(targets_in_stratum, donors, strict=True):
                output[grouped[recipient]] = values[grouped[donor]]
        else:
            rows = np.concatenate([grouped[target] for target in targets_in_stratum])
            shuffled = rows.copy()
            rng.shuffle(shuffled)
            output[rows] = values[shuffled]
    return output


def arm_d_k562_adjustment(
    raw_adjustment: np.ndarray,
    *,
    effect_scale: float,
    raw_median_norm: float,
    pm_shared_mask: np.ndarray,
    pm_available: bool,
    companion_available: bool,
) -> np.ndarray:
    """Map a raw K562-row correction into frozen arm-D pre-whitening space.

    The operation is continuous and never removes a source row.  It mirrors
    arm D's 50:50 K562/PerturbMatch slot and subsequent equal-source average.
    """
    raw = np.asarray(raw_adjustment, dtype=np.float64)
    shared = np.asarray(pm_shared_mask, dtype=bool)
    if raw.shape != shared.shape:
        raise ValueError("adjustment/shared-mask shape mismatch")
    if not pm_available:
        return np.zeros_like(raw)
    scaled = raw * (effect_scale / max(float(raw_median_norm), EPS))
    scaled[shared] *= 0.5
    source_votes = 1 + int(companion_available)
    return scaled / source_votes


def shrink_control(original: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    """Direction-preserving norm match used to reject shrink-only gains."""
    original = np.asarray(original, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    scale = np.linalg.norm(candidate) / max(np.linalg.norm(original), EPS)
    return original * scale
