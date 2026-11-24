#!/usr/bin/env python3
"""Pure numerical utilities for experiment 110.

This module has no file-system readers and no knowledge of held-truth paths.
Callers must construct arrays from separately frozen, provenance-checked
artifacts. Target and gene identifiers are used only as grouping or axis
metadata; they can never enter the selector feature matrix.
"""

from __future__ import annotations

import hashlib
import itertools
from dataclasses import dataclass, replace
from typing import Iterable, Sequence

import numpy as np
from sklearn.linear_model import LogisticRegression


FEATURE_NAMES = (
    "a14_abs_lfc",
    "state_abs_lfc",
    "signed_log_magnitude_ratio",
    "absolute_log_magnitude_ratio",
    "a14_cross_donor_agreement",
    "control_log_cpm",
    "target_level_a14_state_cosine",
    "state_control_bank_sign_agreement",
    "state_control_bank_variance",
)
DEVELOPMENT_CONTEXTS = ("rpe1", "jurkat", "hepg2")
RIDGE_LAMBDAS = (0.01, 0.1, 1.0, 10.0)
PROBABILITY_THRESHOLDS = (0.6, 0.7, 0.8)
STATE_BLEND_FRACTIONS = (0.25, 0.5, 0.75, 1.0)
MAXIMUM_COSINE_DISTANCES = (0.001, 0.003, 0.005)
CONTROL_BANKS = 4
NULL_REPETITIONS = 100
GRID_ATOL = 1e-12


class SelectionError(ValueError):
    """Raised when the frozen selection contract cannot identify one winner."""


@dataclass(frozen=True, order=True)
class CandidateConfig:
    """One member of the preregistered experiment-110 grid."""

    ridge_lambda: float
    probability_threshold: float
    state_blend_fraction: float
    maximum_cosine_distance: float

    def validate(self) -> None:
        if self.ridge_lambda not in RIDGE_LAMBDAS:
            raise ValueError(f"unregistered ridge lambda: {self.ridge_lambda}")
        if self.probability_threshold not in PROBABILITY_THRESHOLDS:
            raise ValueError(
                f"unregistered probability threshold: {self.probability_threshold}"
            )
        if self.state_blend_fraction not in STATE_BLEND_FRACTIONS:
            raise ValueError(f"unregistered STATE blend: {self.state_blend_fraction}")
        if self.maximum_cosine_distance not in MAXIMUM_COSINE_DISTANCES:
            raise ValueError(
                f"unregistered cosine-distance cap: {self.maximum_cosine_distance}"
            )

    def as_array(self) -> np.ndarray:
        return np.asarray(
            (
                self.ridge_lambda,
                self.probability_threshold,
                self.state_blend_fraction,
                self.maximum_cosine_distance,
            ),
            dtype=np.float64,
        )


FIXED_GRID = tuple(
    CandidateConfig(*values)
    for values in itertools.product(
        RIDGE_LAMBDAS,
        PROBABILITY_THRESHOLDS,
        STATE_BLEND_FRACTIONS,
        MAXIMUM_COSINE_DISTANCES,
    )
)


def fixed_grid_array() -> np.ndarray:
    """Return the exact preregistered grid in deterministic product order."""

    return np.vstack([config.as_array() for config in FIXED_GRID])


def sha256_array(array: np.ndarray) -> str:
    """Hash an array with its dtype and shape bound into the digest."""

    values = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(values.dtype.str.encode())
    digest.update(np.asarray(values.shape, dtype=np.int64).tobytes())
    digest.update(values.tobytes())
    return digest.hexdigest()


def stable_seed(*parts: object) -> int:
    payload = "\0".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _float_matrix(name: str, values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError(f"{name} must be two-dimensional")
    if not np.isfinite(array).all():
        raise ValueError(f"{name} must be finite")
    return array


def _bool_matrix(name: str, values: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    array = np.asarray(values, dtype=bool)
    if array.shape != shape:
        raise ValueError(f"{name} shape {array.shape} != {shape}")
    return array


def _gene_mask(name: str, values: np.ndarray, genes: int) -> np.ndarray:
    array = np.asarray(values, dtype=bool)
    if array.shape != (genes,):
        raise ValueError(f"{name} shape {array.shape} != {(genes,)}")
    return array


def disagreement_mask(
    incumbent_lfc: np.ndarray,
    state_lfc: np.ndarray,
    *,
    mutable_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Return finite coordinates where STATE and a14 have different signs."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    state = _float_matrix("state_lfc", state_lfc)
    if state.shape != incumbent.shape:
        raise ValueError("incumbent and STATE shapes differ")
    mask = np.sign(incumbent) != np.sign(state)
    if mutable_mask is not None:
        mask &= _bool_matrix("mutable_mask", mutable_mask, incumbent.shape)
    return mask


def target_level_cosine(
    incumbent_lfc: np.ndarray,
    state_lfc: np.ndarray,
    geometry_mask: np.ndarray,
) -> np.ndarray:
    """Cosine similarity for each target on the prediction-only gene geometry."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    state = _float_matrix("state_lfc", state_lfc)
    if state.shape != incumbent.shape:
        raise ValueError("incumbent and STATE shapes differ")
    use = _gene_mask("geometry_mask", geometry_mask, incumbent.shape[1])
    if not np.any(use):
        raise ValueError("geometry_mask selects no genes")
    left = incumbent[:, use]
    right = state[:, use]
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    out = np.zeros(len(incumbent), dtype=np.float64)
    valid = denominator > 0
    out[valid] = np.einsum("ij,ij->i", left[valid], right[valid]) / denominator[valid]
    np.clip(out, -1.0, 1.0, out=out)
    return out


@dataclass(frozen=True)
class FeatureBatch:
    """Permitted features and their coordinate map."""

    values: np.ndarray
    rows: np.ndarray
    cols: np.ndarray
    state_mean: np.ndarray
    disagreement: np.ndarray


def build_permitted_feature_matrix(
    incumbent_lfc: np.ndarray,
    state_bank_lfc: np.ndarray,
    a14_cross_donor_agreement: np.ndarray,
    control_log_cpm: np.ndarray,
    geometry_mask: np.ndarray,
    *,
    pair_mask: np.ndarray | None = None,
    ratio_epsilon: float = 1e-8,
) -> FeatureBatch:
    """Build exactly the nine preregistered numeric features.

    No target symbols, gene symbols, integer IDs, context IDs, embeddings, or
    arbitrary caller-provided feature columns are accepted by this API.
    """

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    banks = np.asarray(state_bank_lfc, dtype=np.float64)
    if banks.ndim != 3 or banks.shape[0] != CONTROL_BANKS:
        raise ValueError(
            f"state_bank_lfc must have shape ({CONTROL_BANKS}, targets, genes)"
        )
    if banks.shape[1:] != incumbent.shape or not np.isfinite(banks).all():
        raise ValueError("STATE bank predictions are misaligned or non-finite")
    agreement = _float_matrix("a14_cross_donor_agreement", a14_cross_donor_agreement)
    if agreement.shape != incumbent.shape:
        raise ValueError("cross-donor agreement shape differs from effects")
    control = np.asarray(control_log_cpm, dtype=np.float64)
    if control.shape != (incumbent.shape[1],) or not np.isfinite(control).all():
        raise ValueError("control_log_cpm must be one finite value per gene")
    if ratio_epsilon <= 0 or not np.isfinite(ratio_epsilon):
        raise ValueError("ratio_epsilon must be positive and finite")

    state = banks.mean(axis=0)
    if pair_mask is None:
        pairs = disagreement_mask(incumbent, state)
    else:
        pairs = _bool_matrix("pair_mask", pair_mask, incumbent.shape)
    rows, cols = np.nonzero(pairs)

    absolute_incumbent = np.abs(incumbent)
    absolute_state = np.abs(state)
    log_ratio = np.log(
        (absolute_state + ratio_epsilon) / (absolute_incumbent + ratio_epsilon)
    )
    signed_log_ratio = np.sign(incumbent * state) * log_ratio
    bank_sign_agreement = np.mean(
        np.sign(banks) == np.sign(state)[None, :, :],
        axis=0,
        dtype=np.float64,
    )
    bank_variance = np.var(banks, axis=0, ddof=0)
    cosine = target_level_cosine(incumbent, state, geometry_mask)

    columns = (
        absolute_incumbent[rows, cols],
        absolute_state[rows, cols],
        signed_log_ratio[rows, cols],
        np.abs(log_ratio[rows, cols]),
        agreement[rows, cols],
        control[cols],
        cosine[rows],
        bank_sign_agreement[rows, cols],
        bank_variance[rows, cols],
    )
    features = np.column_stack(columns).astype(np.float64, copy=False)
    if features.shape != (len(rows), len(FEATURE_NAMES)):
        raise AssertionError("permitted feature construction changed width")
    if not np.isfinite(features).all():
        raise ValueError("permitted feature matrix contains non-finite values")
    return FeatureBatch(features, rows, cols, state, pairs)


def target_equal_sample_weights(target_blocks: Sequence[object]) -> np.ndarray:
    """Give each distinct target block exactly one unit of total weight."""

    blocks = np.asarray(target_blocks)
    if blocks.ndim != 1 or len(blocks) == 0:
        raise ValueError("target_blocks must be a nonempty one-dimensional array")
    keys = np.asarray([str(value) for value in blocks], dtype="U")
    unique, inverse, counts = np.unique(keys, return_inverse=True, return_counts=True)
    if len(unique) == 0 or np.any(counts <= 0):
        raise AssertionError("invalid target block counts")
    weights = 1.0 / counts[inverse].astype(np.float64)
    totals = np.bincount(inverse, weights=weights)
    if not np.allclose(
        totals, np.ones(len(unique), dtype=np.float64), rtol=0, atol=1e-15
    ):
        raise AssertionError("target-equal weights failed exact block normalization")
    return weights


@dataclass(frozen=True)
class Standardizer:
    mean: np.ndarray
    scale: np.ndarray

    def transform(self, features: np.ndarray) -> np.ndarray:
        values = np.asarray(features, dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != len(self.mean):
            raise ValueError("feature matrix is incompatible with standardizer")
        if not np.isfinite(values).all():
            raise ValueError("feature matrix must be finite")
        return (values - self.mean[None, :]) / self.scale[None, :]


def fit_standardizer(
    training_features: np.ndarray,
    sample_weight: np.ndarray,
) -> Standardizer:
    """Fit weighted means/scales using the training rows only."""

    features = np.asarray(training_features, dtype=np.float64)
    weights = np.asarray(sample_weight, dtype=np.float64)
    if features.ndim != 2 or features.shape[1] != len(FEATURE_NAMES):
        raise ValueError(f"training features must have width {len(FEATURE_NAMES)}")
    if weights.shape != (len(features),):
        raise ValueError("sample weights do not align with training features")
    if (
        len(features) == 0
        or not np.isfinite(features).all()
        or not np.isfinite(weights).all()
        or np.any(weights <= 0)
    ):
        raise ValueError("training features and weights must be finite and positive")
    total = float(weights.sum())
    mean = np.einsum("i,ij->j", weights, features) / total
    variance = np.einsum("i,ij->j", weights, np.square(features - mean)) / total
    scale = np.sqrt(np.maximum(variance, 0.0))
    scale[scale == 0] = 1.0
    return Standardizer(mean.astype(np.float64), scale.astype(np.float64))


@dataclass(frozen=True)
class RidgeLogisticModel:
    """Frozen train-only normalization and weighted ridge logistic parameters."""

    ridge_lambda: float
    standardizer: Standardizer
    coefficient: np.ndarray
    intercept: float

    def predict_probability(self, features: np.ndarray) -> np.ndarray:
        standardized = self.standardizer.transform(features)
        logits = standardized @ self.coefficient + self.intercept
        probabilities = np.empty_like(logits, dtype=np.float64)
        positive = logits >= 0
        probabilities[positive] = 1.0 / (1.0 + np.exp(-logits[positive]))
        exp_logits = np.exp(logits[~positive])
        probabilities[~positive] = exp_logits / (1.0 + exp_logits)
        return probabilities


def fit_weighted_ridge_logistic(
    training_features: np.ndarray,
    labels: np.ndarray,
    target_blocks: Sequence[object],
    ridge_lambda: float,
) -> RidgeLogisticModel:
    """Fit deterministic L2 logistic regression with target-equal weights."""

    if ridge_lambda not in RIDGE_LAMBDAS:
        raise ValueError(f"ridge lambda {ridge_lambda} is outside the frozen grid")
    features = np.asarray(training_features, dtype=np.float64)
    outcomes = np.asarray(labels, dtype=bool)
    if outcomes.shape != (len(features),):
        raise ValueError("labels do not align with training features")
    if len(np.unique(outcomes)) != 2:
        raise ValueError("weighted logistic regression requires both classes")
    weights = target_equal_sample_weights(target_blocks)
    standardizer = fit_standardizer(features, weights)
    standardized = standardizer.transform(features)
    estimator = LogisticRegression(
        C=1.0 / ridge_lambda,
        fit_intercept=True,
        solver="lbfgs",
        max_iter=2_000,
        tol=1e-10,
        random_state=stable_seed("experiment-110", "weighted-ridge") % (2**32),
    )
    estimator.fit(standardized, outcomes.astype(np.int8), sample_weight=weights)
    if estimator.n_iter_[0] >= estimator.max_iter:
        raise RuntimeError("weighted ridge logistic regression did not converge")
    coefficient = np.asarray(estimator.coef_[0], dtype=np.float64)
    intercept = float(estimator.intercept_[0])
    if coefficient.shape != (len(FEATURE_NAMES),):
        raise AssertionError("unexpected selector coefficient width")
    if not np.isfinite(coefficient).all() or not np.isfinite(intercept):
        raise ValueError("selector fit produced non-finite parameters")
    return RidgeLogisticModel(ridge_lambda, standardizer, coefficient, intercept)


def sparse_state_replacement(
    incumbent_lfc: np.ndarray,
    state_lfc: np.ndarray,
    replacement_mask: np.ndarray,
    blend_fraction: float,
) -> np.ndarray:
    """Apply the preregistered sparse STATE replacement formula."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    state = _float_matrix("state_lfc", state_lfc)
    if state.shape != incumbent.shape:
        raise ValueError("incumbent and STATE shapes differ")
    replace_mask = _bool_matrix("replacement_mask", replacement_mask, incumbent.shape)
    if blend_fraction not in STATE_BLEND_FRACTIONS:
        raise ValueError("blend fraction is outside the frozen grid")
    out = incumbent.copy()
    out[replace_mask] += blend_fraction * (
        state[replace_mask] - incumbent[replace_mask]
    )
    return out


def cosine_distance_rows(
    incumbent_lfc: np.ndarray,
    candidate_lfc: np.ndarray,
    geometry_mask: np.ndarray,
) -> np.ndarray:
    """Compute row-level cosine distance with a fail-closed zero-row convention."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    candidate = _float_matrix("candidate_lfc", candidate_lfc)
    if candidate.shape != incumbent.shape:
        raise ValueError("incumbent and candidate shapes differ")
    use = _gene_mask("geometry_mask", geometry_mask, incumbent.shape[1])
    left = incumbent[:, use]
    right = candidate[:, use]
    left_norm = np.linalg.norm(left, axis=1)
    right_norm = np.linalg.norm(right, axis=1)
    out = np.ones(len(left), dtype=np.float64)
    both_zero = (left_norm == 0) & (right_norm == 0)
    out[both_zero] = 0.0
    valid = (left_norm > 0) & (right_norm > 0)
    cosine = np.einsum("ij,ij->i", left[valid], right[valid]) / (
        left_norm[valid] * right_norm[valid]
    )
    out[valid] = 1.0 - np.clip(cosine, -1.0, 1.0)
    return out


def apply_row_cosine_distance_cap(
    incumbent_lfc: np.ndarray,
    candidate_lfc: np.ndarray,
    geometry_mask: np.ndarray,
    maximum_cosine_distance: float,
    *,
    iterations: int = 80,
) -> tuple[np.ndarray, np.ndarray]:
    """Scale each sparse candidate delta until its row meets the angular cap.

    Scaling the candidate delta, rather than rotating the whole row, preserves
    exact sparsity: coordinates untouched by the selector remain bit-identical
    to the incumbent.
    """

    if maximum_cosine_distance not in MAXIMUM_COSINE_DISTANCES:
        raise ValueError("cosine-distance cap is outside the frozen grid")
    if iterations < 40:
        raise ValueError("at least 40 bisection iterations are required")
    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    candidate = _float_matrix("candidate_lfc", candidate_lfc)
    if candidate.shape != incumbent.shape:
        raise ValueError("incumbent and candidate shapes differ")
    use = _gene_mask("geometry_mask", geometry_mask, incumbent.shape[1])
    out = candidate.copy()
    scales = np.ones(len(out), dtype=np.float64)
    distances = cosine_distance_rows(incumbent, candidate, use)
    for row in np.flatnonzero(distances > maximum_cosine_distance):
        base = incumbent[row]
        delta = candidate[row] - base
        if not np.any(delta):
            scales[row] = 0.0
            out[row] = base
            continue
        if np.linalg.norm(base[use]) == 0:
            scales[row] = 0.0
            out[row] = base
            continue
        low = 0.0
        high = 1.0
        for _ in range(iterations):
            middle = 0.5 * (low + high)
            probe = base + middle * delta
            distance = cosine_distance_rows(base[None, :], probe[None, :], use)[0]
            if distance <= maximum_cosine_distance:
                low = middle
            else:
                high = middle
        scales[row] = low
        out[row] = base + low * delta
    final_distance = cosine_distance_rows(incumbent, out, use)
    if np.any(final_distance > maximum_cosine_distance + 5e-13):
        raise AssertionError("row-level cosine cap failed")
    return out, scales


def apply_row_cap_scales(
    incumbent_lfc: np.ndarray,
    uncapped_candidate_lfc: np.ndarray,
    row_scales: np.ndarray,
) -> np.ndarray:
    """Reconstruct a frozen capped candidate from its row delta scales."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    candidate = _float_matrix("uncapped_candidate_lfc", uncapped_candidate_lfc)
    scales = np.asarray(row_scales, dtype=np.float64)
    if candidate.shape != incumbent.shape or scales.shape != (len(incumbent),):
        raise ValueError("row-cap reconstruction shapes differ")
    if not np.isfinite(scales).all() or np.any((scales < 0) | (scales > 1)):
        raise ValueError("row cap scales must be finite and in [0, 1]")
    return incumbent + scales[:, None] * (candidate - incumbent)


def enforce_coordinate_contract(
    candidate_lfc: np.ndarray,
    incumbent_lfc: np.ndarray,
    mutable_mask: np.ndarray,
    forbidden_gene_mask: np.ndarray,
) -> np.ndarray:
    """Enforce exact a14 fallback and exact zero on forbidden coordinates."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    candidate = _float_matrix("candidate_lfc", candidate_lfc)
    if candidate.shape != incumbent.shape:
        raise ValueError("incumbent and candidate shapes differ")
    mutable = _bool_matrix("mutable_mask", mutable_mask, incumbent.shape)
    forbidden = _gene_mask(
        "forbidden_gene_mask", forbidden_gene_mask, incumbent.shape[1]
    )
    if np.any(mutable[:, forbidden]):
        raise ValueError("forbidden coordinates cannot be mutable")
    if np.any(incumbent[:, forbidden] != 0):
        raise ValueError("incumbent is nonzero on a forbidden coordinate")
    out = candidate.copy()
    out[~mutable] = incumbent[~mutable]
    out[:, forbidden] = 0.0
    if not np.array_equal(out[~mutable], incumbent[~mutable]):
        raise AssertionError("unsupported fallback is not exact")
    if np.any(out[:, forbidden] != 0):
        raise AssertionError("forbidden coordinate is not exactly zero")
    return out


def build_candidate(
    context: "DevelopmentContext",
    probabilities: np.ndarray,
    config: CandidateConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build one sparse, capped, coordinate-safe candidate."""

    config.validate()
    q = np.asarray(probabilities, dtype=np.float64)
    if q.shape != (len(context.features),):
        raise ValueError("probabilities do not align with disagreement pairs")
    if not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError("selector probabilities must lie in [0, 1]")
    replace_pairs = q >= config.probability_threshold
    replacement_mask = np.zeros_like(context.mutable_mask, dtype=bool)
    replacement_mask[
        context.pair_rows[replace_pairs], context.pair_cols[replace_pairs]
    ] = True
    if np.any(replacement_mask & ~context.mutable_mask):
        raise AssertionError("selector requested a forbidden replacement")
    uncapped = sparse_state_replacement(
        context.incumbent,
        context.state,
        replacement_mask,
        config.state_blend_fraction,
    )
    capped, row_scales = apply_row_cosine_distance_cap(
        context.incumbent,
        uncapped,
        context.geometry_mask,
        config.maximum_cosine_distance,
    )
    candidate = enforce_coordinate_contract(
        capped,
        context.incumbent,
        context.mutable_mask,
        context.forbidden_gene_mask,
    )
    return candidate, replacement_mask, row_scales


def directional_credit_per_target(
    predicted_lfc: np.ndarray,
    reference_lfc: np.ndarray,
    direction_mask: np.ndarray,
) -> np.ndarray:
    """Compute the experiment-110 effect-level directional-credit proxy."""

    predicted = _float_matrix("predicted_lfc", predicted_lfc)
    reference = np.asarray(reference_lfc, dtype=np.float64)
    if reference.ndim != 2:
        raise ValueError("reference_lfc must be two-dimensional")
    if reference.shape != predicted.shape:
        raise ValueError("predicted and reference shapes differ")
    mask = _bool_matrix("direction_mask", direction_mask, predicted.shape)
    if not np.isfinite(reference[mask]).all():
        raise ValueError("direction mask includes undefined reference LFC")
    values = np.full(len(predicted), np.nan, dtype=np.float64)
    for row in range(len(predicted)):
        use = mask[row]
        if not np.any(use):
            continue
        denominator = float(np.abs(reference[row, use]).sum())
        if denominator <= 0:
            continue
        values[row] = (
            float(np.sign(reference[row, use]) @ predicted[row, use]) / denominator
        )
    return values


def directional_credit_delta(
    candidate_lfc: np.ndarray,
    context: "DevelopmentContext",
) -> np.ndarray:
    return directional_credit_per_target(
        candidate_lfc, context.reference_lfc, context.direction_mask
    ) - directional_credit_per_target(
        context.incumbent, context.reference_lfc, context.direction_mask
    )


@dataclass(frozen=True)
class DevelopmentContext:
    """One context's source-development arrays.

    `target_blocks`, `pair_rows`, and `pair_cols` are used for weighting,
    resampling, and reconstruction only. They are never passed to the model.
    """

    name: str
    features: np.ndarray
    labels: np.ndarray
    target_blocks: np.ndarray
    pair_rows: np.ndarray
    pair_cols: np.ndarray
    incumbent: np.ndarray
    state: np.ndarray
    reference_lfc: np.ndarray
    direction_mask: np.ndarray
    mutable_mask: np.ndarray
    geometry_mask: np.ndarray
    forbidden_gene_mask: np.ndarray
    labels_are_permuted: bool = False

    def validate(self) -> None:
        if self.name not in DEVELOPMENT_CONTEXTS:
            raise ValueError(f"unexpected development context: {self.name}")
        incumbent = _float_matrix("incumbent", self.incumbent)
        state = _float_matrix("state", self.state)
        reference = np.asarray(self.reference_lfc, dtype=np.float64)
        if reference.ndim != 2:
            raise ValueError("reference_lfc must be two-dimensional")
        if state.shape != incumbent.shape or reference.shape != incumbent.shape:
            raise ValueError("context effect matrices are misaligned")
        features = np.asarray(self.features, dtype=np.float64)
        if features.ndim != 2 or features.shape[1] != len(FEATURE_NAMES):
            raise ValueError("context feature matrix has a forbidden width")
        if not np.isfinite(features).all():
            raise ValueError("context feature matrix must be finite")
        labels = np.asarray(self.labels, dtype=bool)
        rows = np.asarray(self.pair_rows, dtype=np.int64)
        cols = np.asarray(self.pair_cols, dtype=np.int64)
        blocks = np.asarray(self.target_blocks)
        expected_vector = (len(features),)
        for name, vector in (
            ("labels", labels),
            ("target_blocks", blocks),
            ("pair_rows", rows),
            ("pair_cols", cols),
        ):
            if vector.shape != expected_vector:
                raise ValueError(f"{name} does not align with feature rows")
        if len(features) == 0:
            raise ValueError("context has no eligible disagreement pairs")
        if (
            np.any(rows < 0)
            or np.any(rows >= incumbent.shape[0])
            or np.any(cols < 0)
            or np.any(cols >= incumbent.shape[1])
        ):
            raise ValueError("pair coordinates are outside the effect matrix")
        linear = rows * incumbent.shape[1] + cols
        if len(np.unique(linear)) != len(linear):
            raise ValueError("pair coordinates contain duplicates")
        mutable = _bool_matrix("mutable_mask", self.mutable_mask, incumbent.shape)
        direction = _bool_matrix("direction_mask", self.direction_mask, incumbent.shape)
        geometry = _gene_mask("geometry_mask", self.geometry_mask, incumbent.shape[1])
        forbidden = _gene_mask(
            "forbidden_gene_mask",
            self.forbidden_gene_mask,
            incumbent.shape[1],
        )
        if np.any(geometry & forbidden):
            raise ValueError("forbidden genes cannot enter cosine geometry")
        if np.any(mutable[:, forbidden]):
            raise ValueError("forbidden genes cannot be mutable")
        if np.any(incumbent[:, forbidden] != 0):
            raise ValueError("incumbent must be zero on forbidden genes")
        if not np.all(mutable[rows, cols]) or not np.all(direction[rows, cols]):
            raise ValueError("training pairs must be mutable directional coordinates")
        if np.any(np.sign(incumbent[rows, cols]) == np.sign(state[rows, cols])):
            raise ValueError("training population contains a sign-agreement pair")
        if not self.labels_are_permuted:
            expected_labels = np.sign(state[rows, cols]) == np.sign(
                reference[rows, cols]
            )
            if not np.array_equal(labels, expected_labels):
                raise ValueError("selector labels do not match STATE/reference signs")
        if not np.all(np.isfinite(reference[rows, cols])):
            raise ValueError("training population contains undefined reference LFC")
        if np.any(reference[rows, cols] == 0):
            raise ValueError("training population contains zero reference LFC")


def make_development_context(
    *,
    name: str,
    incumbent_lfc: np.ndarray,
    state_bank_lfc: np.ndarray,
    a14_cross_donor_agreement: np.ndarray,
    control_log_cpm: np.ndarray,
    reference_lfc: np.ndarray,
    reference_significant: np.ndarray,
    state_supported_gene_mask: np.ndarray,
    forbidden_gene_mask: np.ndarray,
    own_target_coordinate_mask: np.ndarray,
    target_blocks_by_row: Sequence[object],
) -> DevelopmentContext:
    """Construct a validated development context from already-frozen arrays."""

    incumbent = _float_matrix("incumbent_lfc", incumbent_lfc)
    reference = np.asarray(reference_lfc, dtype=np.float64)
    if reference.ndim != 2:
        raise ValueError("reference_lfc must be two-dimensional")
    if reference.shape != incumbent.shape:
        raise ValueError("reference and incumbent shapes differ")
    significant = _bool_matrix(
        "reference_significant", reference_significant, incumbent.shape
    )
    own_target = _bool_matrix(
        "own_target_coordinate_mask",
        own_target_coordinate_mask,
        incumbent.shape,
    )
    supported = _gene_mask(
        "state_supported_gene_mask",
        state_supported_gene_mask,
        incumbent.shape[1],
    )
    forbidden = _gene_mask(
        "forbidden_gene_mask", forbidden_gene_mask, incumbent.shape[1]
    )
    if np.any(incumbent[:, forbidden] != 0):
        raise ValueError("incumbent must already be zero on forbidden coordinates")
    row_blocks = np.asarray([str(value) for value in target_blocks_by_row], dtype="U")
    if row_blocks.shape != (incumbent.shape[0],):
        raise ValueError("target_blocks_by_row must provide one block per target row")
    if len(set(row_blocks)) != len(row_blocks):
        raise ValueError("target_blocks_by_row must be unique within a context")
    mutable = np.broadcast_to(supported & ~forbidden, incumbent.shape).copy()
    direction = significant & np.isfinite(reference) & (reference != 0) & ~own_target
    geometry = ~forbidden
    banks = np.asarray(state_bank_lfc, dtype=np.float64)
    if banks.shape != (CONTROL_BANKS, *incumbent.shape):
        raise ValueError("STATE bank prediction shape differs from contract")
    state = banks.mean(axis=0)
    pairs = disagreement_mask(incumbent, state, mutable_mask=mutable) & direction
    feature_batch = build_permitted_feature_matrix(
        incumbent,
        banks,
        a14_cross_donor_agreement,
        control_log_cpm,
        geometry,
        pair_mask=pairs,
    )
    labels = np.sign(state[feature_batch.rows, feature_batch.cols]) == np.sign(
        reference[feature_batch.rows, feature_batch.cols]
    )
    context = DevelopmentContext(
        name=name,
        features=feature_batch.values,
        labels=labels,
        target_blocks=row_blocks[feature_batch.rows],
        pair_rows=feature_batch.rows,
        pair_cols=feature_batch.cols,
        incumbent=incumbent,
        state=state,
        reference_lfc=reference,
        direction_mask=direction,
        mutable_mask=mutable,
        geometry_mask=geometry,
        forbidden_gene_mask=forbidden,
    )
    context.validate()
    return context


@dataclass(frozen=True)
class GridScores:
    contexts: tuple[str, ...]
    grid: tuple[CandidateConfig, ...]
    rotation_scores: np.ndarray
    rotation_winner_indices: np.ndarray
    rotation_unique_winner: np.ndarray


@dataclass(frozen=True)
class LocoSelection:
    grid_scores: GridScores
    selected_index: int
    selected_config: CandidateConfig
    final_model: RidgeLogisticModel


def _canonical_contexts(
    contexts: Sequence[DevelopmentContext],
) -> tuple[DevelopmentContext, ...]:
    if len(contexts) != len(DEVELOPMENT_CONTEXTS):
        raise ValueError("selection requires exactly three development contexts")
    lookup: dict[str, DevelopmentContext] = {}
    for context in contexts:
        context.validate()
        if context.name in lookup:
            raise ValueError(f"duplicate development context: {context.name}")
        lookup[context.name] = context
    if set(lookup) != set(DEVELOPMENT_CONTEXTS):
        raise ValueError(
            f"development contexts {sorted(lookup)} != {list(DEVELOPMENT_CONTEXTS)}"
        )
    return tuple(lookup[name] for name in DEVELOPMENT_CONTEXTS)


def _pool_contexts(
    contexts: Sequence[DevelopmentContext],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    features = np.concatenate([context.features for context in contexts], axis=0)
    labels = np.concatenate([context.labels for context in contexts], axis=0)
    blocks = np.concatenate(
        [np.asarray(context.target_blocks, dtype="U") for context in contexts]
    )
    return features, labels, blocks


def _unique_winner(scores: np.ndarray) -> tuple[int, bool]:
    values = np.asarray(scores, dtype=np.float64)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError("grid scores must be one-dimensional and finite")
    maximum = float(values.max())
    winners = np.flatnonzero(np.isclose(values, maximum, rtol=0, atol=GRID_ATOL))
    return int(winners[0]), len(winners) == 1


def score_loco_grid(
    contexts: Sequence[DevelopmentContext],
    *,
    grid: Sequence[CandidateConfig] = FIXED_GRID,
) -> GridScores:
    """Score the fixed grid under all leave-one-context-out rotations."""

    ordered = _canonical_contexts(contexts)
    candidates = tuple(grid)
    if candidates != FIXED_GRID:
        raise ValueError("grid differs from preregistered_protocol.json")
    scores = np.empty((len(ordered), len(candidates)), dtype=np.float64)
    winners = np.empty(len(ordered), dtype=np.int64)
    unique = np.empty(len(ordered), dtype=bool)

    for held_index, held in enumerate(ordered):
        training = tuple(
            context for index, context in enumerate(ordered) if index != held_index
        )
        train_features, train_labels, train_blocks = _pool_contexts(training)
        models = {
            ridge: fit_weighted_ridge_logistic(
                train_features, train_labels, train_blocks, ridge
            )
            for ridge in RIDGE_LAMBDAS
        }
        probabilities = {
            ridge: model.predict_probability(held.features)
            for ridge, model in models.items()
        }
        for config_index, config in enumerate(candidates):
            candidate, _, _ = build_candidate(
                held, probabilities[config.ridge_lambda], config
            )
            per_target = directional_credit_delta(candidate, held)
            if not np.isfinite(per_target).any():
                raise ValueError(f"{held.name}: no finite directional-credit targets")
            scores[held_index, config_index] = float(np.nanmean(per_target))
        winners[held_index], unique[held_index] = _unique_winner(scores[held_index])

    return GridScores(
        tuple(context.name for context in ordered),
        candidates,
        scores,
        winners,
        unique,
    )


def select_loco_configuration(
    contexts: Sequence[DevelopmentContext],
) -> LocoSelection:
    """Require one unique, identical winner in all three rotations."""

    ordered = _canonical_contexts(contexts)
    grid_scores = score_loco_grid(ordered)
    if not np.all(grid_scores.rotation_unique_winner):
        ambiguous = [
            grid_scores.contexts[index]
            for index in np.flatnonzero(~grid_scores.rotation_unique_winner)
        ]
        raise SelectionError(f"LOCO rotation has a tied winner: {ambiguous}")
    if len(set(map(int, grid_scores.rotation_winner_indices))) != 1:
        detail = {
            context: int(winner)
            for context, winner in zip(
                grid_scores.contexts,
                grid_scores.rotation_winner_indices,
                strict=True,
            )
        }
        raise SelectionError(f"LOCO rotations disagree on the winner: {detail}")
    selected_index = int(grid_scores.rotation_winner_indices[0])
    selected = grid_scores.grid[selected_index]
    features, labels, blocks = _pool_contexts(ordered)
    final_model = fit_weighted_ridge_logistic(
        features, labels, blocks, selected.ridge_lambda
    )
    return LocoSelection(grid_scores, selected_index, selected, final_model)


def _ordered_block_indices(target_blocks: Sequence[object]) -> list[np.ndarray]:
    blocks = np.asarray([str(value) for value in target_blocks], dtype="U")
    if blocks.ndim != 1 or len(blocks) == 0:
        raise ValueError("target_blocks must be a nonempty vector")
    ordered: dict[str, list[int]] = {}
    for index, block in enumerate(blocks):
        ordered.setdefault(str(block), []).append(index)
    return [np.asarray(indices, dtype=np.int64) for indices in ordered.values()]


def permute_labels_within_target(
    labels: np.ndarray,
    target_blocks: Sequence[object],
    repetition: int,
    *,
    base_seed: int = 20_260_830,
) -> np.ndarray:
    """Deterministically permute labels only within each target block."""

    if not 0 <= repetition < NULL_REPETITIONS:
        raise ValueError(f"repetition must be in [0, {NULL_REPETITIONS - 1}]")
    outcomes = np.asarray(labels, dtype=bool)
    if outcomes.ndim != 1 or len(outcomes) != len(target_blocks):
        raise ValueError("labels and target blocks must be aligned vectors")
    permuted = outcomes.copy()
    for block_number, indices in enumerate(_ordered_block_indices(target_blocks)):
        rng = np.random.default_rng(
            stable_seed(
                "experiment-110",
                "within-target-label-permutation",
                base_seed,
                repetition,
                block_number,
            )
        )
        permuted[indices] = outcomes[indices[rng.permutation(len(indices))]]
        if int(permuted[indices].sum()) != int(outcomes[indices].sum()):
            raise AssertionError("within-target permutation changed class counts")
    return permuted


def iter_within_target_label_permutations(
    labels: np.ndarray,
    target_blocks: Sequence[object],
    *,
    repetitions: int = NULL_REPETITIONS,
    base_seed: int = 20_260_830,
) -> Iterable[np.ndarray]:
    """Yield the exact deterministic permutation sequence without materializing it."""

    if repetitions != NULL_REPETITIONS:
        raise ValueError(f"the preregistered null requires {NULL_REPETITIONS} repeats")
    for repetition in range(repetitions):
        yield permute_labels_within_target(
            labels,
            target_blocks,
            repetition,
            base_seed=base_seed,
        )


@dataclass(frozen=True)
class FullSelectionNull:
    """Models and selection diagnostics from 100 complete null reruns."""

    selected_config_indices: np.ndarray
    selected_configs: np.ndarray
    selected_macro_scores: np.ndarray
    rotation_scores: np.ndarray
    rotation_winner_indices: np.ndarray
    rotation_unique_winner: np.ndarray
    common_rotation_winner: np.ndarray
    label_sha256: np.ndarray
    scaler_mean: np.ndarray
    scaler_scale: np.ndarray
    coefficient: np.ndarray
    intercept: np.ndarray


def run_full_selection_null(
    contexts: Sequence[DevelopmentContext],
    *,
    repetitions: int = NULL_REPETITIONS,
    base_seed: int = 20_260_830,
) -> FullSelectionNull:
    """Repeat label permutation, all LOCO fits, grid scoring, and final refit."""

    if repetitions != NULL_REPETITIONS:
        raise ValueError(f"the preregistered null requires {NULL_REPETITIONS} repeats")
    ordered = _canonical_contexts(contexts)
    selected_indices = np.empty(repetitions, dtype=np.int64)
    selected_configs = np.empty((repetitions, 4), dtype=np.float64)
    selected_scores = np.empty(repetitions, dtype=np.float64)
    all_rotation_scores = np.empty(
        (repetitions, len(DEVELOPMENT_CONTEXTS), len(FIXED_GRID)),
        dtype=np.float64,
    )
    rotation_winners = np.empty(
        (repetitions, len(DEVELOPMENT_CONTEXTS)), dtype=np.int64
    )
    rotation_unique = np.empty((repetitions, len(DEVELOPMENT_CONTEXTS)), dtype=bool)
    common_winner = np.empty(repetitions, dtype=bool)
    label_hashes = np.empty((repetitions, len(DEVELOPMENT_CONTEXTS)), dtype="U64")
    means = np.empty((repetitions, len(FEATURE_NAMES)), dtype=np.float64)
    scales = np.empty_like(means)
    coefficients = np.empty_like(means)
    intercepts = np.empty(repetitions, dtype=np.float64)

    for repetition in range(repetitions):
        permuted_contexts = []
        for context_index, context in enumerate(ordered):
            labels = permute_labels_within_target(
                context.labels,
                context.target_blocks,
                repetition,
                base_seed=stable_seed(base_seed, context.name),
            )
            label_hashes[repetition, context_index] = sha256_array(
                labels.astype(np.uint8)
            )
            permuted_contexts.append(
                replace(context, labels=labels, labels_are_permuted=True)
            )

        grid_scores = score_loco_grid(permuted_contexts)
        all_rotation_scores[repetition] = grid_scores.rotation_scores
        rotation_winners[repetition] = grid_scores.rotation_winner_indices
        rotation_unique[repetition] = grid_scores.rotation_unique_winner
        common_winner[repetition] = bool(
            np.all(grid_scores.rotation_unique_winner)
            and len(set(map(int, grid_scores.rotation_winner_indices))) == 1
        )
        macro_scores = grid_scores.rotation_scores.mean(axis=0)
        selected_index, _ = _unique_winner(macro_scores)
        selected_config = FIXED_GRID[selected_index]
        selected_indices[repetition] = selected_index
        selected_configs[repetition] = selected_config.as_array()
        selected_scores[repetition] = macro_scores[selected_index]

        features, labels, blocks = _pool_contexts(permuted_contexts)
        model = fit_weighted_ridge_logistic(
            features,
            labels,
            blocks,
            selected_config.ridge_lambda,
        )
        means[repetition] = model.standardizer.mean
        scales[repetition] = model.standardizer.scale
        coefficients[repetition] = model.coefficient
        intercepts[repetition] = model.intercept

    return FullSelectionNull(
        selected_indices,
        selected_configs,
        selected_scores,
        all_rotation_scores,
        rotation_winners,
        rotation_unique,
        common_winner,
        label_hashes,
        means,
        scales,
        coefficients,
        intercepts,
    )


def target_block_bootstrap(
    values: np.ndarray,
    target_blocks: Sequence[object],
    *,
    repetitions: int = 2_000,
    base_seed: int = 20_260_830,
    label: str = "directional-credit",
) -> dict[str, np.ndarray | float | int]:
    """Bootstrap target-block means, never individual target-gene pairs."""

    observations = np.asarray(values, dtype=np.float64)
    blocks = np.asarray([str(value) for value in target_blocks], dtype="U")
    if observations.ndim != 1 or blocks.shape != observations.shape:
        raise ValueError("values and target blocks must be aligned vectors")
    if repetitions <= 0:
        raise ValueError("bootstrap repetitions must be positive")
    finite = np.isfinite(observations)
    observations = observations[finite]
    blocks = blocks[finite]
    if len(observations) == 0:
        raise ValueError("bootstrap has no finite observations")
    block_values = np.asarray(
        [
            observations[blocks == block].mean()
            for block in dict.fromkeys(map(str, blocks))
        ],
        dtype=np.float64,
    )
    rng = np.random.default_rng(
        stable_seed("experiment-110", "target-block-bootstrap", base_seed, label)
    )
    draws = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        draws[index] = block_values[
            rng.integers(0, len(block_values), len(block_values))
        ].mean()
    return {
        "estimate": float(block_values.mean()),
        "lower_95": float(np.percentile(draws, 2.5)),
        "upper_95": float(np.percentile(draws, 97.5)),
        "blocks": int(len(block_values)),
        "draws": draws,
    }
