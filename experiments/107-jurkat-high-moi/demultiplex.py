#!/usr/bin/env python3
"""Exact Python port of the frozen deMULTIplex hash-classification core.

The implementation follows ``deMULTIplex`` 1.0.2's
``MULTIseq.Classification.Suite.R`` at commit
``9883367c16d374962955e8eeb9fab720a64d7a57``.  It deliberately contains no
GSE249595 loading code; callers provide one channel's synthetic or already
extracted raw hash-count matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


Q_VALUES = np.arange(1, 100, 2, dtype=np.float64) / 100.0
KDE_GRID_SIZE = 401
EXTREMA_GRID_SIZE = 100
MAXIMUM_ROUNDS = 3
NEGATIVE = "Negative"
DOUBLET = "Doublet"


class DemultiplexError(RuntimeError):
    """The frozen deMULTIplex contract cannot produce a defined result."""


@dataclass(frozen=True)
class HashThresholdModel:
    """KDE evaluation and selected low/high extrema for one hash."""

    evaluation_points: np.ndarray
    density: np.ndarray
    extrema: np.ndarray
    low_extreme: int
    high_extreme: int

    def threshold(self, q: float) -> float:
        """Return the type-7 quantile between the two selected extrema."""

        return type7_quantile(
            np.array(
                [
                    self.evaluation_points[self.high_extreme],
                    self.evaluation_points[self.low_extreme],
                ],
                dtype=np.float64,
            ),
            q,
        )


@dataclass(frozen=True)
class RoundResult:
    """Complete quantile sweep and selected calls for one peeling round."""

    round_number: int
    active_cell_ids: tuple[str, ...]
    normalized_counts: np.ndarray
    hash_models: tuple[HashThresholdModel, ...]
    thresholds: np.ndarray
    calls_by_q: np.ndarray
    singlet_fraction: np.ndarray
    q_local_maxima: np.ndarray
    selected_q_index: int
    selected_q: float
    selected_calls: tuple[str, ...]


@dataclass(frozen=True)
class DemultiplexResult:
    """Final calls in original cell order plus every executed round."""

    cell_ids: tuple[str, ...]
    calls: tuple[str, ...]
    rounds: tuple[RoundResult, ...]

    def as_dict(self) -> dict[str, str]:
        """Return final calls keyed by the original cell identifiers."""

        return dict(zip(self.cell_ids, self.calls, strict=True))


def _as_count_matrix(raw_counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(raw_counts)
    if counts.ndim != 2 or counts.shape[0] == 0 or counts.shape[1] == 0:
        raise DemultiplexError("hash counts must be a nonempty two-dimensional matrix")
    if not np.issubdtype(counts.dtype, np.number):
        raise DemultiplexError("hash counts must be numeric")
    counts = np.asarray(counts, dtype=np.float64)
    if np.any(~np.isfinite(counts)) or np.any(counts < 0):
        raise DemultiplexError("hash counts must be finite and nonnegative")
    return counts


def _validate_names(
    names: Sequence[str],
    expected_length: int,
    description: str,
) -> tuple[str, ...]:
    result = tuple(str(name) for name in names)
    if len(result) != expected_length:
        raise DemultiplexError(
            f"{description} count does not match the hash-count matrix"
        )
    if any(not name for name in result) or len(set(result)) != len(result):
        raise DemultiplexError(f"{description} must be nonempty and unique")
    return result


def normalize_hash_counts(raw_counts: np.ndarray) -> np.ndarray:
    """Apply per-hash log2, nonfinite-to-zero, then column centering."""

    counts = _as_count_matrix(raw_counts)
    with np.errstate(divide="ignore", invalid="ignore"):
        normalized = np.log2(counts)
    normalized[~np.isfinite(normalized)] = 0.0
    normalized -= normalized.mean(axis=0, keepdims=True)
    return normalized


def type7_quantile(values: np.ndarray, probability: float) -> float:
    """Reproduce R's scalar ``quantile(..., type=7)``."""

    data = np.asarray(values, dtype=np.float64)
    if data.ndim != 1 or data.size == 0 or np.any(~np.isfinite(data)):
        raise DemultiplexError("type-7 quantile requires finite one-dimensional data")
    q = float(probability)
    if not np.isfinite(q) or q < 0.0 or q > 1.0:
        raise DemultiplexError("quantile probability must be in [0, 1]")
    ordered = np.sort(data)
    index = (ordered.size - 1) * q
    lower = int(np.floor(index))
    upper = int(np.ceil(index))
    if lower == upper:
        return float(ordered[lower])
    fraction = index - lower
    return float((1.0 - fraction) * ordered[lower] + fraction * ordered[upper])


def r_seq_length_out(start: float, stop: float, length: int) -> np.ndarray:
    """Match ``seq(from=start, to=stop, length.out=length)`` for doubles."""

    if length < 0:
        raise DemultiplexError("sequence length cannot be negative")
    if length == 0:
        return np.empty(0, dtype=np.float64)
    if length == 1:
        return np.array([start], dtype=np.float64)
    if length == 2:
        return np.array([start, stop], dtype=np.float64)
    start = float(start)
    stop = float(stop)
    if start == stop:
        return np.full(length, start, dtype=np.float64)
    output = np.empty(length, dtype=np.float64)
    output[0] = start
    output[-1] = stop
    step = (stop - start) / (length - 1)
    output[1:-1] = start + np.arange(1, length - 1, dtype=np.float64) * step
    return output


def local_maxima(values: Sequence[float]) -> np.ndarray:
    """Return zero-based indices from deMULTIplex's exact ``localMaxima``."""

    data = np.asarray(values)
    if data.ndim != 1 or data.size < 2:
        raise DemultiplexError(
            "local maxima require at least two one-dimensional values"
        )
    if np.any(~np.isfinite(data)):
        raise DemultiplexError("local maxima require finite values")

    rising = (
        np.diff(
            np.concatenate(
                (
                    np.array([-float(np.iinfo(np.int32).max)], dtype=np.float64),
                    np.asarray(data, dtype=np.float64),
                )
            )
        )
        > 0
    )
    run_ends = np.concatenate(
        (np.flatnonzero(rising[1:] != rising[:-1]) + 1, [rising.size])
    )
    maxima = run_ends[::2] - 1
    if data.size >= 2 and data[0] == data[1]:
        maxima = maxima[1:]
    return np.asarray(maxima, dtype=np.int64)


def _linear_bin(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Port KernSmooth's linear binning for points inside an ordered grid."""

    position = (values - grid[0]) * (grid.size - 1) / (grid[-1] - grid[0])
    left = np.floor(position).astype(np.int64)
    fraction = position - left
    if np.any(left < 0) or np.any(left >= grid.size):
        raise DemultiplexError("KDE values fell outside the binning grid")

    counts = np.zeros(grid.size, dtype=np.float64)
    right_edge = left == grid.size - 1
    np.add.at(counts, left[right_edge], 1.0)
    interior = ~right_edge
    np.add.at(counts, left[interior], 1.0 - fraction[interior])
    np.add.at(counts, left[interior] + 1, fraction[interior])
    return counts


def bkde_normal(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce ``KernSmooth::bkde(x, kernel="normal")`` defaults."""

    data = np.asarray(values, dtype=np.float64)
    if data.ndim != 1 or data.size < 2 or np.any(~np.isfinite(data)):
        raise DemultiplexError("KDE requires at least two finite values")

    sample_variance = np.var(data, ddof=1)
    bandwidth = (
        (1.0 / (4.0 * np.pi)) ** 0.1
        * (243.0 / (35.0 * data.size)) ** 0.2
        * np.sqrt(sample_variance)
    )
    if not np.isfinite(bandwidth) or bandwidth <= 0.0:
        raise DemultiplexError("KDE default bandwidth is undefined")

    grid = r_seq_length_out(
        float(np.min(data) - 4.0 * bandwidth),
        float(np.max(data) + 4.0 * bandwidth),
        KDE_GRID_SIZE,
    )
    grid_counts = _linear_bin(data, grid)
    delta = (grid[-1] - grid[0]) / (bandwidth * (KDE_GRID_SIZE - 1))
    kernel_radius = min(int(np.floor(4.0 / delta)), KDE_GRID_SIZE)
    offsets = np.arange(kernel_radius + 1, dtype=np.float64)
    kernel = (
        np.exp(-0.5 * (offsets * delta) ** 2)
        / np.sqrt(2.0 * np.pi)
        / (data.size * bandwidth)
    )

    fft_size = 2 ** int(
        np.ceil(np.log(KDE_GRID_SIZE + kernel_radius + 1) / np.log(2.0))
    )
    padding = fft_size - 2 * kernel_radius - 1
    if padding < 0:
        raise DemultiplexError("KDE convolution grid is invalid")
    kernel = np.concatenate(
        (kernel, np.zeros(padding, dtype=np.float64), kernel[1:][::-1])
    )
    grid_counts = np.pad(grid_counts, (0, fft_size - KDE_GRID_SIZE))
    total = np.sum(kernel) * (grid[-1] - grid[0]) / (KDE_GRID_SIZE - 1) * data.size
    density = np.fft.ifft(np.fft.fft(kernel / total) * np.fft.fft(grid_counts)).real[
        :KDE_GRID_SIZE
    ]
    return grid, density


def fit_hash_threshold_model(values: np.ndarray) -> HashThresholdModel:
    """Fit one hash KDE and select the frozen low/high extrema."""

    data = np.asarray(values, dtype=np.float64)
    kde_x, kde_y = bkde_normal(data)
    evaluation_points = r_seq_length_out(
        type7_quantile(data, 0.001),
        type7_quantile(data, 0.999),
        EXTREMA_GRID_SIZE,
    )
    density = np.interp(evaluation_points, kde_x, kde_y)
    extrema = local_maxima(density)
    if extrema.size <= 1:
        raise DemultiplexError("hash KDE has fewer than two local maxima")

    low_extreme = int(extrema[np.argmax(density[extrema])])
    high_extreme = int(extrema[-1])
    if low_extreme == high_extreme:
        raise DemultiplexError("hash KDE low and high extrema are identical")
    return HashThresholdModel(
        evaluation_points=evaluation_points,
        density=density,
        extrema=extrema,
        low_extreme=low_extreme,
        high_extreme=high_extreme,
    )


def classify_normalized(
    normalized_counts: np.ndarray,
    hash_labels: Sequence[str],
    thresholds: np.ndarray,
) -> tuple[str, ...]:
    """Classify zero, one, or multiple ``>=`` threshold hits."""

    normalized = np.asarray(normalized_counts, dtype=np.float64)
    labels = _validate_names(hash_labels, normalized.shape[1], "hash labels")
    threshold = np.asarray(thresholds, dtype=np.float64)
    if normalized.ndim != 2 or threshold.shape != (normalized.shape[1],):
        raise DemultiplexError("thresholds must match the normalized hash columns")
    if np.any(~np.isfinite(normalized)) or np.any(~np.isfinite(threshold)):
        raise DemultiplexError("classification values and thresholds must be finite")

    hits = normalized >= threshold
    hit_count = hits.sum(axis=1)
    calls = np.full(normalized.shape[0], NEGATIVE, dtype=object)
    calls[hit_count > 1] = DOUBLET
    singlets = hit_count == 1
    if np.any(singlets):
        calls[singlets] = np.asarray(labels, dtype=object)[
            np.argmax(hits[singlets], axis=1)
        ]
    return tuple(str(call) for call in calls)


def select_quantile(
    singlet_fraction: np.ndarray,
    q_values: np.ndarray = Q_VALUES,
) -> tuple[int, np.ndarray]:
    """Choose the first best singlet-fraction local maximum."""

    fractions = np.asarray(singlet_fraction, dtype=np.float64)
    q = np.asarray(q_values, dtype=np.float64)
    if fractions.ndim != 1 or fractions.shape != q.shape:
        raise DemultiplexError("singlet fractions must match the q grid")
    maxima = local_maxima(fractions)
    if maxima.size == 0:
        raise DemultiplexError("singlet-fraction sweep has no local maximum")
    selected = int(maxima[np.argmax(fractions[maxima])])
    return selected, maxima


def classify_round(
    raw_counts: np.ndarray,
    cell_ids: Sequence[str],
    hash_labels: Sequence[str],
    *,
    round_number: int = 1,
    q_values: np.ndarray = Q_VALUES,
) -> RoundResult:
    """Run one complete frozen q sweep on the current active cell set."""

    counts = _as_count_matrix(raw_counts)
    cells = _validate_names(cell_ids, counts.shape[0], "cell identifiers")
    labels = _validate_names(hash_labels, counts.shape[1], "hash labels")
    q = np.asarray(q_values, dtype=np.float64)
    if (
        q.ndim != 1
        or q.size == 0
        or np.any(~np.isfinite(q))
        or np.any(q < 0.0)
        or np.any(q > 1.0)
        or np.any(np.diff(q) <= 0.0)
    ):
        raise DemultiplexError(
            "q values must be finite, unique, ascending, and in [0, 1]"
        )

    normalized = normalize_hash_counts(counts)
    models = tuple(
        fit_hash_threshold_model(normalized[:, column])
        for column in range(normalized.shape[1])
    )
    thresholds = np.empty((q.size, normalized.shape[1]), dtype=np.float64)
    calls_by_q = np.empty((q.size, normalized.shape[0]), dtype=object)
    singlet_fraction = np.empty(q.size, dtype=np.float64)

    for q_index, q_value in enumerate(q):
        thresholds[q_index] = [model.threshold(float(q_value)) for model in models]
        calls = classify_normalized(normalized, labels, thresholds[q_index])
        calls_by_q[q_index] = calls
        call_array = np.asarray(calls, dtype=object)
        singlet_fraction[q_index] = np.mean(
            (call_array != NEGATIVE) & (call_array != DOUBLET)
        )

    selected_q_index, q_local_maxima = select_quantile(singlet_fraction, q)
    return RoundResult(
        round_number=round_number,
        active_cell_ids=cells,
        normalized_counts=normalized,
        hash_models=models,
        thresholds=thresholds,
        calls_by_q=calls_by_q,
        singlet_fraction=singlet_fraction,
        q_local_maxima=q_local_maxima,
        selected_q_index=selected_q_index,
        selected_q=float(q[selected_q_index]),
        selected_calls=tuple(
            str(call) for call in calls_by_q[selected_q_index].tolist()
        ),
    )


def demultiplex(
    raw_counts: np.ndarray,
    cell_ids: Sequence[str],
    hash_labels: Sequence[str],
    *,
    maximum_rounds: int = MAXIMUM_ROUNDS,
) -> DemultiplexResult:
    """Run exactly the frozen negative-peeling state machine."""

    counts = _as_count_matrix(raw_counts)
    cells = _validate_names(cell_ids, counts.shape[0], "cell identifiers")
    labels = _validate_names(hash_labels, counts.shape[1], "hash labels")
    if maximum_rounds != MAXIMUM_ROUNDS:
        raise DemultiplexError("the frozen protocol requires exactly three rounds")

    active = np.arange(counts.shape[0], dtype=np.int64)
    final = np.full(counts.shape[0], None, dtype=object)
    rounds: list[RoundResult] = []

    for round_number in range(1, MAXIMUM_ROUNDS + 1):
        result = classify_round(
            counts[active],
            tuple(cells[index] for index in active),
            labels,
            round_number=round_number,
        )
        rounds.append(result)
        selected = np.asarray(result.selected_calls, dtype=object)
        negative = selected == NEGATIVE
        final[active[negative]] = NEGATIVE

        if not np.any(negative):
            final[active] = selected
            break

        survivors = active[~negative]
        if survivors.size == 0:
            break
        if round_number == MAXIMUM_ROUNDS:
            final[survivors] = selected[~negative]
            break
        active = survivors

    if any(call is None for call in final):
        raise DemultiplexError("negative peeling lost one or more cells")
    return DemultiplexResult(
        cell_ids=cells,
        calls=tuple(str(call) for call in final),
        rounds=tuple(rounds),
    )
