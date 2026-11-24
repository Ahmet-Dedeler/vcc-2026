"""Pure geometry and scoring primitives for Experiment 137."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np


CAP_DEGREES = (0.0, 1.5, 3.0, 6.0)
DEGENERATE_TOLERANCE = 1e-12


class GeometryError(ValueError):
    """Raised when a candidate violates the frozen geometry contract."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise GeometryError(message)


def _matrix(values: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(values)
    require(array.ndim == 2, f"{name} must be a matrix")
    require(np.issubdtype(array.dtype, np.floating), f"{name} must be floating")
    require(np.all(np.isfinite(array)), f"{name} must be finite")
    return array


def spherical_tangent_candidate(
    incumbent: np.ndarray,
    source: np.ndarray,
    transferable: np.ndarray,
    coordinate_supported: np.ndarray,
    cap_degrees: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotate supported incumbent subvectors toward source by at most one cap.

    Unsupported rows and columns remain bit-identical. The zero-cap branch is
    returned before any floating-point geometry is evaluated.
    """

    base = _matrix(incumbent, "incumbent")
    proposal = _matrix(source, "source")
    require(proposal.shape == base.shape, "source and incumbent shapes differ")
    rows = np.asarray(transferable)
    columns = np.asarray(coordinate_supported)
    require(rows.dtype == np.bool_, "transferable must be boolean")
    require(columns.dtype == np.bool_, "coordinate_supported must be boolean")
    require(rows.shape == (base.shape[0],), "transferable shape differs")
    require(columns.shape == (base.shape[1],), "coordinate support shape differs")
    cap = float(cap_degrees)
    require(cap in CAP_DEGREES, "angular cap is outside the frozen grid")

    candidate = base.copy()
    effective = np.zeros(base.shape[0], dtype=bool)
    realized = np.zeros(base.shape[0], dtype=np.float64)
    if cap == 0.0:
        return candidate, effective, realized

    cap_radians = math.radians(cap)
    for row in np.flatnonzero(rows):
        b = np.asarray(base[row, columns], dtype=np.float64)
        s = np.asarray(proposal[row, columns], dtype=np.float64)
        b_norm = float(np.linalg.norm(b))
        s_norm = float(np.linalg.norm(s))
        if (
            not math.isfinite(b_norm)
            or not math.isfinite(s_norm)
            or b_norm <= DEGENERATE_TOLERANCE
            or s_norm <= DEGENERATE_TOLERANCE
        ):
            continue
        b_unit = b / b_norm
        s_unit = s / s_norm
        cosine = float(np.clip(np.dot(b_unit, s_unit), -1.0, 1.0))
        tangent = s_unit - cosine * b_unit
        tangent_norm = float(np.linalg.norm(tangent))
        if not math.isfinite(tangent_norm) or tangent_norm <= DEGENERATE_TOLERANCE:
            continue
        phi = math.atan2(tangent_norm, cosine)
        theta = min(phi, cap_radians)
        tangent_unit = tangent / tangent_norm
        rotated = b_norm * (
            math.cos(theta) * b_unit + math.sin(theta) * tangent_unit
        )
        candidate[row, columns] = rotated.astype(base.dtype, copy=False)
        effective[row] = True
        realized[row] = math.degrees(theta)

    assert_exact_fallback(base, candidate, effective, columns)
    return candidate, effective, realized


def assert_exact_fallback(
    incumbent: np.ndarray,
    candidate: np.ndarray,
    effective: np.ndarray,
    coordinate_supported: np.ndarray,
) -> None:
    base = np.asarray(incumbent)
    changed = np.asarray(candidate)
    rows = np.asarray(effective, dtype=bool)
    columns = np.asarray(coordinate_supported, dtype=bool)
    require(changed.shape == base.shape, "candidate and incumbent shapes differ")
    require(rows.shape == (base.shape[0],), "effective rows differ")
    require(columns.shape == (base.shape[1],), "supported columns differ")
    require(
        np.array_equal(
            np.ascontiguousarray(changed[~rows]).view(np.uint8),
            np.ascontiguousarray(base[~rows]).view(np.uint8),
        ),
        "an ineffective target changed",
    )
    require(
        np.array_equal(
            np.ascontiguousarray(changed[:, ~columns]).view(np.uint8),
            np.ascontiguousarray(base[:, ~columns]).view(np.uint8),
        ),
        "an unsupported coordinate changed",
    )


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Official cosine-rank PDS with exact midranks."""

    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    require(
        predicted.ndim == 2 and predicted.shape == observed.shape,
        "PDS arrays differ",
    )
    require(predicted.shape[0] >= 2, "PDS needs at least two targets")
    require(
        np.all(np.isfinite(predicted)) and np.all(np.isfinite(observed)),
        "PDS arrays must be finite",
    )
    pred_norm = np.linalg.norm(predicted, axis=1)
    truth_norm = np.linalg.norm(observed, axis=1)
    pred_unit = np.divide(
        predicted,
        pred_norm[:, None],
        out=np.zeros_like(predicted),
        where=pred_norm[:, None] > 0,
    )
    truth_unit = np.divide(
        observed,
        truth_norm[:, None],
        out=np.zeros_like(observed),
        where=truth_norm[:, None] > 0,
    )
    similarity = pred_unit @ truth_unit.T
    own = np.diag(similarity)[:, None]
    rank = (similarity > own).sum(axis=1).astype(np.float64)
    rank += 0.5 * ((similarity == own).sum(axis=1) - 1)
    result = 1.0 - rank / (predicted.shape[0] - 1)
    result[pred_norm == 0] = 0.5
    return result


def signed_credit_per_target(
    prediction: np.ndarray, truth: np.ndarray
) -> np.ndarray:
    """Target-equal directional credit without truth-derived filtering."""

    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    require(
        predicted.ndim == 2 and predicted.shape == observed.shape,
        "signed-credit arrays differ",
    )
    require(
        np.all(np.isfinite(predicted)) and np.all(np.isfinite(observed)),
        "signed-credit arrays must be finite",
    )
    numerator = np.sum(np.sign(observed) * predicted, axis=1)
    denominator = np.maximum(np.sum(np.abs(predicted), axis=1), 1e-300)
    return numerator / denominator


def select_cap(
    fold_scores: dict[float, Sequence[np.ndarray]],
    caps: Sequence[float] = CAP_DEGREES,
) -> tuple[tuple[float, ...], np.ndarray]:
    """Five-fold nested selection with the smallest-cap tie break."""

    grid = tuple(sorted(map(float, caps)))
    require(grid == CAP_DEGREES, "cap grid changed")
    require(set(fold_scores) == set(grid), "cap score grid differs")
    fold_count = len(fold_scores[grid[0]])
    require(fold_count == 5, "nested selection requires five folds")
    means = np.empty((len(grid), fold_count), dtype=np.float64)
    for cap_index, cap in enumerate(grid):
        require(len(fold_scores[cap]) == fold_count, "fold counts differ")
        for fold, values in enumerate(fold_scores[cap]):
            array = np.asarray(values, dtype=np.float64)
            require(array.ndim == 1 and array.size == 63, "fold size changed")
            require(np.all(np.isfinite(array)), "fold scores are non-finite")
            means[cap_index, fold] = float(array.mean())
    selected: list[float] = []
    held = np.empty(fold_count, dtype=np.float64)
    for fold in range(fold_count):
        objective = np.delete(means, fold, axis=1).mean(axis=1)
        chosen = int(np.flatnonzero(objective == objective.max())[0])
        selected.append(grid[chosen])
        held[fold] = means[chosen, fold]
    return tuple(selected), held


def final_cap(
    fold_scores: dict[float, Sequence[np.ndarray]],
    caps: Sequence[float] = CAP_DEGREES,
) -> float:
    grid = tuple(sorted(map(float, caps)))
    require(grid == CAP_DEGREES, "cap grid changed")
    require(set(fold_scores) == set(grid), "cap score grid differs")
    objective: list[float] = []
    for cap in grid:
        folds = fold_scores[cap]
        require(len(folds) == 5, "final selection requires five folds")
        means: list[float] = []
        for values in folds:
            array = np.asarray(values, dtype=np.float64)
            require(array.ndim == 1 and array.size == 63, "fold size changed")
            require(np.all(np.isfinite(array)), "fold scores are non-finite")
            means.append(float(array.mean()))
        objective.append(float(np.mean(means)))
    return grid[objective.index(max(objective))]
