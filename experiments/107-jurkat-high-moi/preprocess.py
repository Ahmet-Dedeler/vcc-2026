#!/usr/bin/env python3
"""Frozen numerical primitives for experiment 107 preprocessing.

This module deliberately contains only deterministic, testable pieces.  The
real GSE249595 transcriptome coordinates must not be parsed until
``test_preprocess.py`` and ``test_protocol.py`` pass.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import BinaryIO

import numpy as np
from scipy.special import gammaln
from scipy.stats import binom

try:
    from numba import njit
except ImportError:  # pragma: no cover - the project environment includes numba.
    njit = None


EMPTYDROPS_BASE_SEED = 1_072_026
GUIDE_HYPOTHESES = 83_401
GUIDE_FDR = 0.05


class PreprocessError(RuntimeError):
    """The frozen preprocessing contract cannot be satisfied."""


@dataclass(frozen=True)
class BarcodeRankResult:
    """The two thresholds returned by the pinned barcode-rank algorithm."""

    knee: float
    inflection: float


@dataclass(frozen=True)
class MonteCarloResult:
    """Monte Carlo p-values and lower-bound flags in input order."""

    p_value: np.ndarray
    limited: np.ndarray
    exceedances: np.ndarray


@dataclass(frozen=True)
class MatrixMarketSummary:
    """Compact first-pass statistics for one transcriptome matrix."""

    rows: int
    columns: int
    declared_nnz: int
    observed_nnz: int
    coordinate_order: str
    totals: np.ndarray
    detected_features: np.ndarray
    mitochondrial_counts: np.ndarray
    feature_totals: np.ndarray


def _matrix_market_header(stream: BinaryIO) -> tuple[int, int, int]:
    banner = stream.readline().decode("ascii").strip().split()
    if banner != [
        "%%MatrixMarket",
        "matrix",
        "coordinate",
        "integer",
        "general",
    ]:
        raise PreprocessError("expected an integer general coordinate matrix")
    while True:
        raw = stream.readline()
        if not raw:
            raise PreprocessError("Matrix Market dimensions are missing")
        if not raw.startswith(b"%"):
            break
    fields = raw.split()
    if len(fields) != 3:
        raise PreprocessError("Matrix Market dimensions must have three fields")
    rows, columns, nnz = (int(field) for field in fields)
    if rows <= 0 or columns <= 0 or nnz < 0:
        raise PreprocessError("Matrix Market dimensions are invalid")
    return rows, columns, nnz


def _iter_matrix_market_coordinates(
    stream: BinaryIO,
) -> tuple[tuple[int, int, int], Sequence[tuple[int, int, int]]]:
    """Yield validated zero-based coordinates from a Matrix Market stream."""

    rows, columns, declared_nnz = _matrix_market_header(stream)

    def iterator() -> Sequence[tuple[int, int, int]]:
        observed = 0
        for line_number, raw in enumerate(stream, start=1):
            fields = raw.split()
            if len(fields) != 3:
                raise PreprocessError(
                    f"coordinate line {line_number} does not have three fields"
                )
            row, column, value = (int(field) for field in fields)
            if not 1 <= row <= rows or not 1 <= column <= columns:
                raise PreprocessError(
                    f"coordinate line {line_number} is outside matrix dimensions"
                )
            if value <= 0:
                raise PreprocessError(
                    f"coordinate line {line_number} is not a positive integer"
                )
            observed += 1
            if observed > declared_nnz:
                raise PreprocessError("matrix has more entries than declared")
            yield row - 1, column - 1, value
        if observed != declared_nnz:
            raise PreprocessError(
                f"matrix has {observed} entries but declares {declared_nnz}"
            )

    return (rows, columns, declared_nnz), iterator()


def scan_transcriptome_summary(
    stream: BinaryIO,
    *,
    mitochondrial_features: np.ndarray | None = None,
) -> MatrixMarketSummary:
    """Scan one raw transcriptome coordinate stream without materializing it."""

    (rows, columns, declared_nnz), coordinates = _iter_matrix_market_coordinates(stream)
    mitochondrial = np.zeros(rows, dtype=bool)
    if mitochondrial_features is not None:
        indices = np.asarray(mitochondrial_features, dtype=np.int64)
        if np.any(indices < 0) or np.any(indices >= rows):
            raise PreprocessError("mitochondrial feature index is out of bounds")
        mitochondrial[indices] = True

    totals = np.zeros(columns, dtype=np.int64)
    detected = np.zeros(columns, dtype=np.uint16)
    mitochondrial_counts = np.zeros(columns, dtype=np.int64)
    feature_totals = np.zeros(rows, dtype=np.int64)
    previous_row = previous_column = -1
    row_major = column_major = True
    observed_nnz = 0

    for row, column, value in coordinates:
        if observed_nnz:
            row_major &= (row, column) > (previous_row, previous_column)
            column_major &= (column, row) > (previous_column, previous_row)
        previous_row, previous_column = row, column
        observed_nnz += 1
        totals[column] += value
        detected[column] += 1
        feature_totals[row] += value
        if mitochondrial[row]:
            mitochondrial_counts[column] += value

    if not (row_major or column_major):
        raise PreprocessError("matrix coordinates are not in a strict canonical order")
    coordinate_order = "column_major" if column_major else "row_major"
    return MatrixMarketSummary(
        rows=rows,
        columns=columns,
        declared_nnz=declared_nnz,
        observed_nnz=observed_nnz,
        coordinate_order=coordinate_order,
        totals=totals,
        detected_features=detected,
        mitochondrial_counts=mitochondrial_counts,
        feature_totals=feature_totals,
    )


def ambient_feature_counts(
    stream: BinaryIO,
    *,
    totals: np.ndarray,
    lower: int = 100,
) -> np.ndarray:
    """Sum feature counts over the exact EmptyDrops ambient barcode set."""

    (rows, columns, _), coordinates = _iter_matrix_market_coordinates(stream)
    barcode_totals = np.asarray(totals, dtype=np.int64)
    if barcode_totals.shape != (columns,):
        raise PreprocessError("barcode totals do not match matrix columns")
    ambient = barcode_totals <= lower
    output = np.zeros(rows, dtype=np.int64)
    for row, column, value in coordinates:
        if ambient[column]:
            output[row] += value
    return output


def observed_data_log_probabilities(
    stream: BinaryIO,
    *,
    totals: np.ndarray,
    ambient_probability: np.ndarray,
    tested: np.ndarray,
) -> np.ndarray:
    """Compute DropletUtils' sparse data term for selected barcodes."""

    (rows, columns, _), coordinates = _iter_matrix_market_coordinates(stream)
    barcode_totals = np.asarray(totals, dtype=np.int64)
    probability = np.asarray(ambient_probability, dtype=np.float64)
    selected = np.asarray(tested, dtype=bool)
    if barcode_totals.shape != (columns,) or selected.shape != (columns,):
        raise PreprocessError("barcode vectors do not match matrix columns")
    if probability.shape != (rows,):
        raise PreprocessError("ambient probabilities do not match matrix rows")
    if np.any(probability < 0) or not np.isclose(
        probability.sum(),
        1.0,
        rtol=0,
        atol=1e-12,
    ):
        raise PreprocessError(
            "ambient probabilities must be nonnegative and sum to one"
        )

    output = np.zeros(columns, dtype=np.float64)
    log_probability = np.full(probability.shape, -np.inf, dtype=np.float64)
    positive = probability > 0
    log_probability[positive] = np.log(probability[positive])
    for row, column, value in coordinates:
        if selected[column]:
            if not positive[row]:
                raise PreprocessError(
                    "positive observed count has zero ambient probability"
                )
            output[column] += value * log_probability[row] - gammaln(value + 1)
    if np.any(barcode_totals[selected] <= 0):
        raise PreprocessError("tested barcodes must have positive totals")
    return output[selected]


def _interpolate_on_curve(
    target_distance: np.ndarray,
    cumulative_distance: np.ndarray,
    stepwise_distance: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Port DropletUtils ``.interpolate_on_curve`` with left-open intervals."""

    index = np.searchsorted(cumulative_distance, target_distance, side="left") - 1
    if np.any(index < 0) or np.any(index + 1 >= cumulative_distance.size):
        raise PreprocessError("curve interpolation target is outside the curve")
    index_p1 = index + 1
    proportion = (target_distance - cumulative_distance[index]) / stepwise_distance[
        index_p1
    ]
    return (
        x[index] + proportion * (x[index_p1] - x[index]),
        y[index] + proportion * (y[index_p1] - y[index]),
    )


def barcode_ranks_thresholds(
    totals: np.ndarray,
    *,
    lower: int = 100,
    exclude_from: int = 50,
    window: float = 1.0,
    gradient_threshold: float = -1.0,
) -> BarcodeRankResult:
    """Port the knee/inflection calculation in DropletUtils 1.31.2."""

    values = np.asarray(totals)
    if values.ndim != 1 or values.size == 0:
        raise PreprocessError("totals must be a nonempty one-dimensional array")
    if np.any(values < 0) or np.any(values != np.floor(values)):
        raise PreprocessError("totals must contain nonnegative integers")

    descending = np.sort(values.astype(np.int64, copy=False))[::-1]
    change = np.r_[True, descending[1:] != descending[:-1]]
    starts = np.flatnonzero(change)
    lengths = np.diff(np.r_[starts, descending.size])
    run_totals = descending[starts].astype(np.float64)
    run_ranks = np.cumsum(lengths) - (lengths - 1) / 2

    keep = (run_totals > lower) & (run_ranks > exclude_from)
    if np.count_nonzero(keep) < 2:
        raise PreprocessError(
            "insufficient unique points for computing knee/inflection points"
        )

    y = np.log10(run_totals[keep])
    x = np.log10(run_ranks[keep])
    stepwise = np.r_[0.0, np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)]
    cumulative = np.cumsum(stepwise)
    rhs = cumulative + window
    scan = rhs <= cumulative[-1]

    if not np.any(scan):
        threshold = float(10 ** y[-1])
        return BarcodeRankResult(knee=threshold, inflection=threshold)

    left_x = x[scan]
    left_y = y[scan]
    right_x, right_y = _interpolate_on_curve(
        rhs[scan],
        cumulative,
        stepwise,
        x,
        y,
    )
    window_gap = np.sqrt((left_x - right_x) ** 2 + (left_y - right_y) ** 2)
    window_gradient = (right_y - left_y) / (right_x - left_x)
    window_intercept = right_y - window_gradient * right_x
    mid_x, mid_y = _interpolate_on_curve(
        cumulative[scan] + window / 2,
        cumulative,
        stepwise,
        x,
        y,
    )
    mid_above = mid_y > window_gradient * mid_x + window_intercept

    elbows = np.flatnonzero((~mid_above) & (window_gradient < gradient_threshold))
    if elbows.size == 0:
        inflection_window = int(np.argmin(window_gradient))
        knee_candidates = np.flatnonzero(mid_above)
    else:
        first_elbow = int(elbows[0])
        knee_candidates = np.flatnonzero(mid_above[: first_elbow + 1])
        before_first_elbow = int(
            np.searchsorted(left_x, mid_x[first_elbow], side="right")
        )
        inflection_window = int(np.argmin(window_gradient[:before_first_elbow]))

    if knee_candidates.size:
        knee_window = int(knee_candidates[np.argmin(window_gap[knee_candidates])])
    else:
        knee_window = inflection_window

    return BarcodeRankResult(
        knee=float(10 ** mid_y[knee_window]),
        inflection=float(10 ** mid_y[inflection_window]),
    )


def simple_good_turing(counts: np.ndarray, *, confidence: float = 1.96) -> np.ndarray:
    """Port edgeR's simple Good-Turing estimator for one count vector."""

    observed = np.asarray(counts)
    if observed.ndim != 1 or observed.size == 0:
        raise PreprocessError("Good-Turing counts must be a nonempty vector")
    if np.any(observed < 0) or np.any(observed != np.floor(observed)):
        raise PreprocessError("Good-Turing counts must be nonnegative integers")
    observed = observed.astype(np.int64, copy=False)

    positive = observed[observed > 0]
    if positive.size == 0:
        raise PreprocessError("no counts available to estimate the ambient profile")
    frequencies = np.bincount(positive)
    r = np.flatnonzero(frequencies).astype(np.int64)
    n = frequencies[r].astype(np.float64)
    n_zero = int(np.count_nonzero(observed == 0))

    big_n = float(np.dot(r, n))
    log_r = np.log(r.astype(np.float64))
    previous = np.r_[0, r[:-1]]
    following_gap = np.r_[r[1:] - previous[:-1], 2 * (r[-1] - previous[-1])]
    z = 2 * n / following_gap
    centered_x = log_r - log_r.mean()
    centered_y = np.log(z) - np.log(z).mean()
    denominator = float(np.dot(centered_x, centered_x))
    if denominator == 0:
        raise PreprocessError("Good-Turing regression has zero frequency variance")
    slope = float(np.dot(centered_x, centered_y) / denominator)

    p_zero = float(n[0] / big_n) if r[0] == 1 else 0.0
    adjusted = np.empty(r.size, dtype=np.float64)
    use_smoothed = False
    for index, current in enumerate(r):
        next_observed = int(current + 1)
        smoothed = next_observed * np.exp(
            slope * (np.log(next_observed) - log_r[index])
        )
        if index == r.size - 1 or r[index + 1] != next_observed:
            use_smoothed = True
        if not use_smoothed:
            empirical = next_observed * n[index + 1] / n[index]
            tolerance = (
                confidence * empirical * np.sqrt(1.0 / n[index + 1] + 1.0 / n[index])
            )
            if abs(empirical - smoothed) <= tolerance:
                use_smoothed = True
            else:
                adjusted[index] = empirical
        if use_smoothed:
            adjusted[index] = smoothed

    adjusted *= (1.0 - p_zero) / float(np.dot(adjusted, n))
    output = np.empty(observed.size, dtype=np.float64)
    if n_zero:
        output[observed == 0] = p_zero / n_zero
    positions = np.searchsorted(r, observed[observed > 0])
    output[observed > 0] = adjusted[positions]
    return output


def safe_good_turing_proportions(counts: np.ndarray) -> np.ndarray:
    """Port DropletUtils ``.safe_good_turing`` including zero protection."""

    observed = np.asarray(counts, dtype=np.int64)
    proportions = simple_good_turing(observed)
    still_zero = proportions <= 0
    if np.any(still_zero):
        total = int(observed.sum())
        if total <= 0:
            raise PreprocessError("no counts available to estimate the ambient profile")
        pseudo_probability = 1.0 / total
        proportions[still_zero] = pseudo_probability / int(still_zero.sum())
        proportions[~still_zero] *= 1.0 - pseudo_probability
    return proportions


def multinomial_log_probability(
    feature_indices: np.ndarray,
    counts: np.ndarray,
    total: int,
    ambient_probability: np.ndarray,
) -> float:
    """Return the pinned multinomial log probability for one sparse barcode."""

    indices = np.asarray(feature_indices, dtype=np.int64)
    values = np.asarray(counts, dtype=np.int64)
    probabilities = np.asarray(ambient_probability, dtype=np.float64)
    if indices.shape != values.shape or indices.ndim != 1:
        raise PreprocessError(
            "sparse feature indices and counts must be matching vectors"
        )
    if total != int(values.sum()):
        raise PreprocessError("sparse counts do not sum to the supplied total")
    selected = probabilities[indices]
    if np.any(selected <= 0):
        raise PreprocessError("positive observed count has zero ambient probability")
    return float(
        gammaln(total + 1)
        + np.dot(values, np.log(selected))
        - np.sum(gammaln(values + 1))
    )


def multinomial_data_log_probability(
    feature_indices: np.ndarray,
    counts: np.ndarray,
    ambient_probability: np.ndarray,
) -> float:
    """Return the total-independent term compared by DropletUtils."""

    values = np.asarray(counts, dtype=np.int64)
    return multinomial_log_probability(
        feature_indices,
        values,
        int(values.sum()),
        ambient_probability,
    ) - float(gammaln(int(values.sum()) + 1))


def _simulated_prefix_log_probabilities(
    choices: np.ndarray,
    log_probability: np.ndarray,
    checkpoints: np.ndarray,
) -> np.ndarray:
    """Evaluate one multinomial draw path at ascending total-count checkpoints."""

    selected = np.asarray(choices, dtype=np.int64)
    log_p = np.asarray(log_probability, dtype=np.float64)
    stops = np.asarray(checkpoints, dtype=np.int64)
    if stops.ndim != 1 or np.any(np.diff(stops) <= 0):
        raise PreprocessError("Monte Carlo checkpoints must be strictly increasing")
    if stops.size == 0 or stops[-1] != selected.size:
        raise PreprocessError("last Monte Carlo checkpoint must equal draw count")

    feature_counts = np.zeros(log_p.size, dtype=np.int64)
    output = np.empty(stops.size, dtype=np.float64)
    current = 0.0
    checkpoint_index = 0
    for draw_index, feature in enumerate(selected, start=1):
        feature_counts[feature] += 1
        current += log_p[feature] - np.log(feature_counts[feature])
        if draw_index == stops[checkpoint_index]:
            output[checkpoint_index] = current
            checkpoint_index += 1
            if checkpoint_index == stops.size:
                break
    return output


def _simulated_prefix_log_probabilities_batch_python(
    choices: np.ndarray,
    log_probability: np.ndarray,
    checkpoints: np.ndarray,
) -> np.ndarray:
    return np.stack(
        [
            _simulated_prefix_log_probabilities(row, log_probability, checkpoints)
            for row in choices
        ],
        axis=0,
    )


if njit is not None:

    @njit(cache=True)
    def _simulated_prefix_log_probabilities_batch_numba(
        choices: np.ndarray,
        log_probability: np.ndarray,
        checkpoints: np.ndarray,
    ) -> np.ndarray:
        """Compiled prefix likelihoods for a block of independent draw paths."""

        block_size, draw_count = choices.shape
        output = np.empty((block_size, checkpoints.size), dtype=np.float64)
        for block_index in range(block_size):
            feature_counts = np.zeros(log_probability.size, dtype=np.int64)
            current = 0.0
            checkpoint_index = 0
            for draw_index in range(draw_count):
                feature = choices[block_index, draw_index]
                feature_counts[feature] += 1
                current += log_probability[feature] - np.log(feature_counts[feature])
                if draw_index + 1 == checkpoints[checkpoint_index]:
                    output[block_index, checkpoint_index] = current
                    checkpoint_index += 1
                    if checkpoint_index == checkpoints.size:
                        break
        return output


def simulated_prefix_log_probabilities_batch(
    choices: np.ndarray,
    log_probability: np.ndarray,
    checkpoints: np.ndarray,
) -> np.ndarray:
    """Evaluate a block of multinomial paths at ascending checkpoints."""

    selected = np.asarray(choices)
    if selected.ndim != 2:
        raise PreprocessError("Monte Carlo choices must be a two-dimensional block")
    log_p = np.asarray(log_probability, dtype=np.float64)
    stops = np.asarray(checkpoints, dtype=np.int64)
    if stops.ndim != 1 or np.any(np.diff(stops) <= 0):
        raise PreprocessError("Monte Carlo checkpoints must be strictly increasing")
    if stops.size == 0 or stops[-1] != selected.shape[1]:
        raise PreprocessError("last Monte Carlo checkpoint must equal draw count")
    if np.any(selected < 0) or np.any(selected >= log_p.size):
        raise PreprocessError("Monte Carlo choice is outside the feature axis")
    if njit is None:
        return _simulated_prefix_log_probabilities_batch_python(
            selected,
            log_p,
            stops,
        )
    return _simulated_prefix_log_probabilities_batch_numba(
        selected.astype(np.int32, copy=False),
        log_p,
        stops,
    )


def emptydrops_monte_carlo(
    totals: np.ndarray,
    observed_data_log_probability: np.ndarray,
    ambient_probability: np.ndarray,
    *,
    channel_number: int,
    iterations: int = 10_000,
    base_seed: int = EMPTYDROPS_BASE_SEED,
    iteration_block_size: int = 64,
) -> MonteCarloResult:
    """Compute frozen EmptyDrops multinomial Monte Carlo p-values.

    Each iteration has an independent seed, so changing ``iteration_block_size``
    cannot change the result.  The block argument controls scheduling only.
    """

    barcode_totals = np.asarray(totals, dtype=np.int64)
    observed = np.asarray(observed_data_log_probability, dtype=np.float64)
    ambient = np.asarray(ambient_probability, dtype=np.float64)
    if barcode_totals.ndim != 1 or observed.shape != barcode_totals.shape:
        raise PreprocessError("Monte Carlo totals and probabilities must match")
    if np.any(barcode_totals <= 0) or np.any(~np.isfinite(observed)):
        raise PreprocessError(
            "tested barcodes need positive totals and finite probabilities"
        )
    if ambient.ndim != 1 or np.any(ambient < 0) or not np.any(ambient > 0):
        raise PreprocessError(
            "ambient probabilities must be a nonnegative nonzero vector"
        )
    if not np.isclose(ambient.sum(), 1.0, rtol=0, atol=1e-12):
        raise PreprocessError("ambient probabilities must sum to one")
    if iterations <= 0 or iteration_block_size <= 0 or channel_number <= 0:
        raise PreprocessError(
            "Monte Carlo iteration, block, and channel values must be positive"
        )

    order = np.lexsort((observed, barcode_totals))
    sorted_totals = barcode_totals[order]
    sorted_observed = observed[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_totals)) + 1]
    ends = np.r_[starts[1:], sorted_totals.size]
    unique_totals = sorted_totals[starts]
    events = np.zeros(sorted_totals.size, dtype=np.int64)
    log_ambient = np.full(ambient.shape, -np.inf, dtype=np.float64)
    positive = ambient > 0
    log_ambient[positive] = np.log(ambient[positive])
    maximum_total = int(unique_totals[-1])

    for block_start in range(0, iterations, iteration_block_size):
        block_end = min(iterations, block_start + iteration_block_size)
        choices = np.empty(
            (block_end - block_start, maximum_total),
            dtype=np.int32,
        )
        for block_index, iteration in enumerate(range(block_start, block_end)):
            seed = np.random.SeedSequence([base_seed, channel_number, iteration])
            generator = np.random.Generator(np.random.PCG64(seed))
            choices[block_index] = generator.choice(
                ambient.size,
                size=maximum_total,
                p=ambient,
            )
        simulated = simulated_prefix_log_probabilities_batch(
            choices,
            log_ambient,
            unique_totals,
        )
        for group_index, (start, end) in enumerate(zip(starts, ends, strict=True)):
            insertions = np.searchsorted(
                sorted_observed[start:end],
                simulated[:, group_index],
                side="left",
            )
            valid = insertions[insertions < end - start]
            if valid.size:
                events[start:end] += np.bincount(
                    valid,
                    minlength=end - start,
                )

    sorted_exceedances = np.empty_like(events)
    for start, end in zip(starts, ends, strict=True):
        sorted_exceedances[start:end] = np.cumsum(events[start:end])
    exceedances = np.empty_like(sorted_exceedances)
    exceedances[order] = sorted_exceedances
    return MonteCarloResult(
        p_value=(exceedances + 1) / (iterations + 1),
        limited=exceedances == 0,
        exceedances=exceedances,
    )


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """R-compatible BH adjustment for a finite one-dimensional vector."""

    values = np.asarray(p_values, dtype=np.float64)
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise PreprocessError("BH p-values must be a finite one-dimensional vector")
    if np.any((values < 0) | (values > 1)):
        raise PreprocessError("BH p-values must lie in [0, 1]")
    count = values.size
    order = np.argsort(values, kind="stable")
    ranked = values[order] * count / np.arange(1, count + 1)
    adjusted_sorted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted


def emptydrops_fdr(
    totals: np.ndarray,
    p_values: np.ndarray,
    *,
    lower: int,
    retain: float,
) -> np.ndarray:
    """Apply DropletUtils retention and BH behavior to Monte Carlo p-values."""

    barcode_totals = np.asarray(totals, dtype=np.int64)
    probabilities = np.asarray(p_values, dtype=np.float64)
    if probabilities.shape != barcode_totals.shape:
        raise PreprocessError("EmptyDrops totals and p-values must match")
    output = np.full(probabilities.shape, np.nan, dtype=np.float64)
    tested = barcode_totals > lower
    if np.any(~np.isfinite(probabilities[tested])):
        raise PreprocessError("tested EmptyDrops p-values must be finite")
    temporary = probabilities[tested].copy()
    temporary[barcode_totals[tested] >= retain] = 0.0
    output[tested] = benjamini_hochberg(temporary)
    return output


def sparse_guide_qvalues(
    guide_indices: np.ndarray,
    observed_umi: np.ndarray,
    *,
    total_guide_umi: int,
    channel_probability: np.ndarray,
    hypothesis_count: int = GUIDE_HYPOTHESES,
) -> np.ndarray:
    """BH-adjust nonzero guide observations while accounting for all guides."""

    indices = np.asarray(guide_indices, dtype=np.int64)
    counts = np.asarray(observed_umi, dtype=np.int64)
    probability = np.asarray(channel_probability, dtype=np.float64)
    if indices.shape != counts.shape or indices.ndim != 1:
        raise PreprocessError("guide indices and counts must be matching vectors")
    if probability.shape != (hypothesis_count,):
        raise PreprocessError(
            "channel probability axis does not match hypothesis count"
        )
    if np.unique(indices).size != indices.size:
        raise PreprocessError("guide indices must be unique within a cell")
    if np.any(indices < 0) or np.any(indices >= hypothesis_count):
        raise PreprocessError("guide index is outside the canonical guide axis")
    if np.any(counts <= 0) or total_guide_umi != int(counts.sum()):
        raise PreprocessError("guide counts must be positive and sum to the cell total")

    selected_probability = probability[indices]
    if np.any(selected_probability <= 0):
        raise PreprocessError("positive guide observation has zero channel probability")
    p_values = binom.sf(counts - 1, total_guide_umi, selected_probability)

    order = np.lexsort((indices, p_values))
    ranked = p_values[order] * hypothesis_count / np.arange(1, indices.size + 1)
    adjusted_sorted = np.minimum.accumulate(
        np.r_[ranked, np.ones(hypothesis_count - indices.size)][::-1]
    )[::-1][: indices.size]
    adjusted = np.empty(indices.size, dtype=np.float64)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return adjusted


def called_guides(
    guide_indices: np.ndarray,
    observed_umi: np.ndarray,
    *,
    total_guide_umi: int,
    channel_probability: np.ndarray,
    hypothesis_count: int = GUIDE_HYPOTHESES,
    fdr: float = GUIDE_FDR,
) -> np.ndarray:
    """Return canonical guide indices passing the frozen guide-call rule."""

    q_values = sparse_guide_qvalues(
        guide_indices,
        observed_umi,
        total_guide_umi=total_guide_umi,
        channel_probability=channel_probability,
        hypothesis_count=hypothesis_count,
    )
    return np.asarray(guide_indices, dtype=np.int64)[q_values <= fdr]


def run_demultiplex_rounds(
    cell_ids: Sequence[str],
    classify_round: Callable[[tuple[str, ...], int], dict[str, str]],
    *,
    maximum_rounds: int = 3,
) -> dict[str, str]:
    """Execute the frozen deMULTIplex negative-peeling state machine."""

    active = tuple(cell_ids)
    if len(set(active)) != len(active):
        raise PreprocessError("cell identifiers must be unique")
    permanent_negative: set[str] = set()
    current: dict[str, str] = {}

    for round_number in range(1, maximum_rounds + 1):
        current = classify_round(active, round_number)
        if set(current) != set(active):
            raise PreprocessError(
                "round classifier did not return exactly the active cells"
            )
        new_negative = {cell for cell, label in current.items() if label == "Negative"}
        permanent_negative.update(new_negative)
        if not new_negative:
            break
        active = tuple(cell for cell in active if cell not in new_negative)
        if not active:
            current = {}
            break

    final = {cell: "Negative" for cell in permanent_negative}
    final.update(current)
    if set(final) != set(cell_ids):
        raise PreprocessError("deMULTIplex state machine lost or duplicated cells")
    return final


def within_group_demean(values: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Demean every column independently within each integer group."""

    matrix = np.asarray(values, dtype=np.float64)
    labels = np.asarray(groups)
    if matrix.ndim != 2 or labels.shape != (matrix.shape[0],):
        raise PreprocessError("group labels must match the matrix rows")
    output = matrix.copy()
    for label in np.unique(labels):
        selected = labels == label
        output[selected] -= output[selected].mean(axis=0, keepdims=True)
    return output


def residualize_training_and_held(
    nuisance_training: np.ndarray,
    response_training: np.ndarray,
    target_training: np.ndarray,
    nuisance_held: np.ndarray,
    response_held: np.ndarray,
    target_held: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fit nuisance slopes on training data and apply them to both outcomes."""

    nuisance_training = np.asarray(nuisance_training, dtype=np.float64)
    nuisance_held = np.asarray(nuisance_held, dtype=np.float64)
    response_training = np.asarray(response_training, dtype=np.float64)
    response_held = np.asarray(response_held, dtype=np.float64)
    target_training = np.asarray(target_training, dtype=np.float64)
    target_held = np.asarray(target_held, dtype=np.float64)
    response_slope = np.linalg.pinv(nuisance_training) @ response_training
    target_slope = np.linalg.pinv(nuisance_training) @ target_training
    return (
        response_training - nuisance_training @ response_slope,
        target_training - nuisance_training @ target_slope,
        response_held - nuisance_held @ response_slope,
        target_held - nuisance_held @ target_slope,
    )


def identifiable_projection(
    residualized_target_design: np.ndarray,
    *,
    relative_tolerance: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the identifiable row-space projector and singular values."""

    design = np.asarray(residualized_target_design, dtype=np.float64)
    if design.ndim != 2 or design.shape[1] == 0:
        raise PreprocessError("target design must be a nonempty matrix")
    _, singular_values, right = np.linalg.svd(design, full_matrices=False)
    if singular_values.size == 0 or singular_values[0] <= 0:
        return np.zeros((design.shape[1], design.shape[1])), singular_values
    keep = singular_values / singular_values[0] > relative_tolerance
    identifiable = right[keep]
    return identifiable.T @ identifiable, singular_values


def project_identifiable_effect(
    residualized_target_design: np.ndarray,
    coefficients_binary_units: np.ndarray,
    *,
    relative_tolerance: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Project fitted coefficients onto the design-identifiable row space."""

    projection, singular_values = identifiable_projection(
        residualized_target_design,
        relative_tolerance=relative_tolerance,
    )
    return projection @ np.asarray(coefficients_binary_units), singular_values


def masked_blend_exact_fallback(
    incumbent: np.ndarray,
    source: np.ndarray,
    supported: np.ndarray,
    fraction: float,
) -> np.ndarray:
    """Blend supported entries while preserving all other values bit-exactly."""

    incumbent_array = np.asarray(incumbent)
    source_array = np.asarray(source)
    mask = np.asarray(supported, dtype=bool)
    if (
        incumbent_array.shape != source_array.shape
        or mask.shape != incumbent_array.shape
    ):
        raise PreprocessError(
            "blend arrays and support mask must have identical shapes"
        )
    output = incumbent_array.copy()
    output[mask] = (1.0 - fraction) * incumbent_array[mask] + fraction * source_array[
        mask
    ]
    return output
