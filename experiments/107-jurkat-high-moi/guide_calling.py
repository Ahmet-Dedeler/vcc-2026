#!/usr/bin/env python3
"""Sparse, streaming guide calling for experiment 107.

The raw guide matrices have 83,401 rows and millions of columns per channel.
This module therefore keeps only one dense vector on the guide axis and
represents retained hash-singlet cells as CSR-like batches.  Callers must make
two passes over the retained cells:

1. estimate channel-wide guide probabilities from every retained hash singlet;
2. call each retained cell against those frozen probabilities.

No transcriptome values are read here.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

import numpy as np
from scipy.stats import binom


GUIDE_HYPOTHESES = 83_401
GUIDE_FDR = 0.05
NON_TARGETING_GUIDES = 250
NO_UNIQUE_TARGET = -1


class GuideCallingError(RuntimeError):
    """The frozen experiment-107 guide-calling contract was violated."""


def _integer_vector(values: np.ndarray, *, name: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 1:
        raise GuideCallingError(f"{name} must be one-dimensional")
    if not np.issubdtype(array.dtype, np.integer) or np.issubdtype(
        array.dtype, np.bool_
    ):
        raise GuideCallingError(f"{name} must contain integers")
    if np.issubdtype(array.dtype, np.unsignedinteger) and array.size:
        if int(array.max()) > np.iinfo(np.int64).max:
            raise GuideCallingError(f"{name} exceeds int64 range")
    return array.astype(np.int64, copy=False)


@dataclass(frozen=True)
class SparseGuideBatch:
    """CSR-like guide counts for retained hash-singlet cells only.

    ``cell_indices`` contains canonical channel-column indices. ``indptr``
    partitions ``guide_indices`` and ``observed_umi`` by cell. Empty retained
    cells are represented by equal adjacent ``indptr`` entries.
    """

    cell_indices: np.ndarray
    indptr: np.ndarray
    guide_indices: np.ndarray
    observed_umi: np.ndarray

    def __post_init__(self) -> None:
        cells = _integer_vector(self.cell_indices, name="cell_indices")
        pointers = _integer_vector(self.indptr, name="indptr")
        guides = _integer_vector(self.guide_indices, name="guide_indices")
        counts = _integer_vector(self.observed_umi, name="observed_umi")

        if pointers.size != cells.size + 1:
            raise GuideCallingError(
                "indptr length must equal retained cell count plus one"
            )
        if pointers.size == 0 or pointers[0] != 0:
            raise GuideCallingError("indptr must start at zero")
        if np.any(np.diff(pointers) < 0) or pointers[-1] != guides.size:
            raise GuideCallingError("indptr does not partition sparse guide entries")
        if guides.shape != counts.shape:
            raise GuideCallingError(
                "guide_indices and observed_umi must have matching shapes"
            )
        if np.any(cells < 0) or np.any(np.diff(cells) <= 0):
            raise GuideCallingError(
                "cell_indices must be nonnegative and strictly increasing"
            )
        if np.any(counts <= 0):
            raise GuideCallingError("stored guide observations must be positive")

        for start, stop in zip(pointers[:-1], pointers[1:], strict=True):
            row = guides[start:stop]
            if np.unique(row).size != row.size:
                raise GuideCallingError("guide indices must be unique within a cell")

        object.__setattr__(self, "cell_indices", cells)
        object.__setattr__(self, "indptr", pointers)
        object.__setattr__(self, "guide_indices", guides)
        object.__setattr__(self, "observed_umi", counts)

    @property
    def cell_count(self) -> int:
        return int(self.cell_indices.size)

    def iter_cells(self) -> Iterator[tuple[int, np.ndarray, np.ndarray]]:
        """Yield ``(cell_index, guide_indices, observed_umi)`` without copies."""

        for row, cell_index in enumerate(self.cell_indices):
            start = int(self.indptr[row])
            stop = int(self.indptr[row + 1])
            yield (
                int(cell_index),
                self.guide_indices[start:stop],
                self.observed_umi[start:stop],
            )


@dataclass(frozen=True)
class ChannelGuideNull:
    """Channel-wide null estimated from all retained hash-singlet cells."""

    hypothesis_count: int
    retained_cell_count: int
    total_guide_umi: int
    guide_umi: np.ndarray
    probability: np.ndarray


@dataclass(frozen=True)
class CellGuideCall:
    """Sparse guide calls and downstream design summaries for one cell."""

    cell_index: int
    total_guide_umi: int
    called_guide_indices: np.ndarray
    called_observed_umi: np.ndarray
    called_q_values: np.ndarray
    total_called_guide_burden: int
    non_targeting_guide_burden: int
    target_indices: np.ndarray


def _validated_probability(
    channel_probability: np.ndarray,
    *,
    hypothesis_count: int,
) -> np.ndarray:
    if hypothesis_count <= 0:
        raise GuideCallingError("hypothesis_count must be positive")
    probability = np.asarray(channel_probability, dtype=np.float64)
    if probability.shape != (hypothesis_count,):
        raise GuideCallingError(
            "channel probability axis does not match hypothesis_count"
        )
    if np.any(~np.isfinite(probability)) or np.any(
        (probability < 0) | (probability > 1)
    ):
        raise GuideCallingError("channel probabilities must be finite and in [0, 1]")
    if not np.isclose(probability.sum(), 1.0, rtol=0, atol=1e-12):
        raise GuideCallingError("channel probabilities must sum to one")
    return probability


def _validated_annotations(
    guide_to_target: np.ndarray,
    non_targeting_guides: np.ndarray,
    *,
    hypothesis_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    target_map = _integer_vector(guide_to_target, name="guide_to_target")
    ntc = np.asarray(non_targeting_guides)
    if target_map.shape != (hypothesis_count,):
        raise GuideCallingError("guide_to_target axis does not match hypothesis_count")
    if np.any(target_map < NO_UNIQUE_TARGET):
        raise GuideCallingError(
            "guide_to_target values must be -1 or nonnegative target indices"
        )
    if ntc.shape != (hypothesis_count,) or ntc.dtype != np.bool_:
        raise GuideCallingError(
            "non_targeting_guides must be a boolean canonical guide-axis vector"
        )
    if (
        hypothesis_count == GUIDE_HYPOTHESES
        and np.count_nonzero(ntc) != NON_TARGETING_GUIDES
    ):
        raise GuideCallingError(
            "the canonical guide axis must mark exactly 250 non-targeting records"
        )
    return target_map, ntc


def _validated_fdr(fdr: float) -> float:
    try:
        threshold = float(fdr)
    except (TypeError, ValueError) as error:
        raise GuideCallingError("fdr must lie in [0, 1]") from error
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise GuideCallingError("fdr must lie in [0, 1]")
    return threshold


def fit_channel_guide_null(
    retained_hash_singlet_batches: Iterable[SparseGuideBatch],
    *,
    hypothesis_count: int = GUIDE_HYPOTHESES,
) -> ChannelGuideNull:
    """Estimate guide probabilities from every retained hash-singlet cell.

    The iterable is consumed once and never materialized. At most one
    ``hypothesis_count``-length integer vector is accumulated.
    """

    if hypothesis_count <= 0:
        raise GuideCallingError("hypothesis_count must be positive")

    guide_umi = np.zeros(hypothesis_count, dtype=np.int64)
    retained_cell_count = 0
    total_guide_umi = 0
    previous_cell_index = -1

    for batch in retained_hash_singlet_batches:
        if not isinstance(batch, SparseGuideBatch):
            raise GuideCallingError(
                "retained_hash_singlet_batches must contain SparseGuideBatch values"
            )
        if batch.cell_count:
            if int(batch.cell_indices[0]) <= previous_cell_index:
                raise GuideCallingError(
                    "retained cell indices must be globally strictly increasing"
                )
            previous_cell_index = int(batch.cell_indices[-1])
        if np.any(batch.guide_indices < 0) or np.any(
            batch.guide_indices >= hypothesis_count
        ):
            raise GuideCallingError("guide index is outside the canonical guide axis")

        np.add.at(guide_umi, batch.guide_indices, batch.observed_umi)
        retained_cell_count += batch.cell_count
        total_guide_umi += int(batch.observed_umi.sum(dtype=np.int64))

    if retained_cell_count == 0:
        raise GuideCallingError("channel has no retained hash-singlet cells")
    if total_guide_umi == 0:
        raise GuideCallingError(
            "retained hash-singlet cells contain no guide observations"
        )
    if int(guide_umi.sum(dtype=np.int64)) != total_guide_umi:
        raise GuideCallingError("channel guide UMI accumulation overflowed")

    probability = guide_umi.astype(np.float64) / total_guide_umi
    return ChannelGuideNull(
        hypothesis_count=hypothesis_count,
        retained_cell_count=retained_cell_count,
        total_guide_umi=total_guide_umi,
        guide_umi=guide_umi,
        probability=probability,
    )


def _sparse_observed_q_values_validated(
    guide_indices: np.ndarray,
    observed_umi: np.ndarray,
    *,
    total_guide_umi: int,
    channel_probability: np.ndarray,
    hypothesis_count: int,
) -> np.ndarray:
    guides = _integer_vector(guide_indices, name="guide_indices")
    counts = _integer_vector(observed_umi, name="observed_umi")

    if guides.shape != counts.shape:
        raise GuideCallingError(
            "guide_indices and observed_umi must have matching shapes"
        )
    if guides.size > hypothesis_count:
        raise GuideCallingError("cell has more observations than guide hypotheses")
    if np.unique(guides).size != guides.size:
        raise GuideCallingError("guide indices must be unique within a cell")
    if np.any(guides < 0) or np.any(guides >= hypothesis_count):
        raise GuideCallingError("guide index is outside the canonical guide axis")
    if np.any(counts <= 0):
        raise GuideCallingError("stored guide observations must be positive")
    if total_guide_umi < 0 or total_guide_umi != int(counts.sum(dtype=np.int64)):
        raise GuideCallingError("observed guide counts must sum to total_guide_umi")
    if guides.size == 0:
        return np.empty(0, dtype=np.float64)

    selected_probability = channel_probability[guides]
    if np.any(selected_probability == 0):
        raise GuideCallingError(
            "positive guide observation has zero channel probability"
        )

    p_values = np.asarray(
        binom.sf(counts - 1, total_guide_umi, selected_probability),
        dtype=np.float64,
    )
    if np.any(~np.isfinite(p_values)):
        raise GuideCallingError("binomial guide p-value is non-finite")

    # All omitted hypotheses have p=1. Their minimum reverse-BH contribution
    # is exactly one at rank hypothesis_count, so one is the tail sentinel.
    order = np.lexsort((guides, p_values))
    ranks = np.arange(1, guides.size + 1, dtype=np.float64)
    ranked = p_values[order] * hypothesis_count / ranks
    adjusted_sorted = np.minimum.accumulate(np.minimum(ranked, 1.0)[::-1])[::-1]
    adjusted = np.empty(guides.size, dtype=np.float64)
    adjusted[order] = adjusted_sorted
    return adjusted


def sparse_observed_q_values(
    guide_indices: np.ndarray,
    observed_umi: np.ndarray,
    *,
    total_guide_umi: int,
    channel_probability: np.ndarray,
    hypothesis_count: int = GUIDE_HYPOTHESES,
) -> np.ndarray:
    """Return BH q-values for observed guides without allocating a dense cell.

    Omitted guides have observation zero and therefore p-value one. Sorting is
    by p-value and then canonical guide index, which is exactly a stable sort
    of a dense canonical-axis p-value vector. Returned q-values follow the
    input order.
    """

    probability = _validated_probability(
        channel_probability,
        hypothesis_count=hypothesis_count,
    )
    return _sparse_observed_q_values_validated(
        guide_indices,
        observed_umi,
        total_guide_umi=total_guide_umi,
        channel_probability=probability,
        hypothesis_count=hypothesis_count,
    )


def _call_sparse_cell_validated(
    cell_index: int,
    guide_indices: np.ndarray,
    observed_umi: np.ndarray,
    *,
    channel_probability: np.ndarray,
    guide_to_target: np.ndarray,
    non_targeting_guides: np.ndarray,
    hypothesis_count: int,
    fdr: float,
) -> CellGuideCall:
    if not isinstance(cell_index, (int, np.integer)) or int(cell_index) < 0:
        raise GuideCallingError("cell_index must be a nonnegative integer")

    guides = _integer_vector(guide_indices, name="guide_indices")
    counts = _integer_vector(observed_umi, name="observed_umi")
    total = int(counts.sum(dtype=np.int64))
    q_values = _sparse_observed_q_values_validated(
        guides,
        counts,
        total_guide_umi=total,
        channel_probability=channel_probability,
        hypothesis_count=hypothesis_count,
    )
    selected = q_values <= fdr
    called_guides = guides[selected]
    called_counts = counts[selected]
    called_q = q_values[selected]

    canonical_order = np.argsort(called_guides, kind="stable")
    called_guides = called_guides[canonical_order]
    called_counts = called_counts[canonical_order]
    called_q = called_q[canonical_order]

    eligible_targets = guide_to_target[called_guides]
    target_indices = np.unique(eligible_targets[eligible_targets != NO_UNIQUE_TARGET])
    total_burden = int(called_guides.size)
    ntc_burden = int(np.count_nonzero(non_targeting_guides[called_guides]))
    return CellGuideCall(
        cell_index=int(cell_index),
        total_guide_umi=total,
        called_guide_indices=called_guides,
        called_observed_umi=called_counts,
        called_q_values=called_q,
        total_called_guide_burden=total_burden,
        non_targeting_guide_burden=ntc_burden,
        target_indices=target_indices,
    )


def call_sparse_cell(
    cell_index: int,
    guide_indices: np.ndarray,
    observed_umi: np.ndarray,
    *,
    channel_probability: np.ndarray,
    guide_to_target: np.ndarray,
    non_targeting_guides: np.ndarray,
    hypothesis_count: int = GUIDE_HYPOTHESES,
    fdr: float = GUIDE_FDR,
) -> CellGuideCall:
    """Call one sparse cell and collapse only unique single-target guides.

    ``guide_to_target`` must contain ``-1`` for every non-targeting, unmapped,
    or composite guide record. Only a record already proven to map to exactly
    one source target may contain a nonnegative target index.
    """

    probability = _validated_probability(
        channel_probability,
        hypothesis_count=hypothesis_count,
    )
    target_map, ntc = _validated_annotations(
        guide_to_target,
        non_targeting_guides,
        hypothesis_count=hypothesis_count,
    )
    threshold = _validated_fdr(fdr)
    return _call_sparse_cell_validated(
        cell_index,
        guide_indices,
        observed_umi,
        channel_probability=probability,
        guide_to_target=target_map,
        non_targeting_guides=ntc,
        hypothesis_count=hypothesis_count,
        fdr=threshold,
    )


def iter_channel_guide_calls(
    retained_hash_singlet_batches: Iterable[SparseGuideBatch],
    *,
    channel_null: ChannelGuideNull,
    guide_to_target: np.ndarray,
    non_targeting_guides: np.ndarray,
    fdr: float = GUIDE_FDR,
) -> Iterator[CellGuideCall]:
    """Stream calls for a second pass over retained hash-singlet cells."""

    if not isinstance(channel_null, ChannelGuideNull):
        raise GuideCallingError("channel_null must be a ChannelGuideNull")

    probability = _validated_probability(
        channel_null.probability,
        hypothesis_count=channel_null.hypothesis_count,
    )
    guide_umi = _integer_vector(channel_null.guide_umi, name="channel_null.guide_umi")
    if guide_umi.shape != (channel_null.hypothesis_count,) or np.any(guide_umi < 0):
        raise GuideCallingError(
            "channel_null guide UMI axis must be nonnegative and canonical"
        )
    if int(guide_umi.sum(dtype=np.int64)) != channel_null.total_guide_umi:
        raise GuideCallingError("channel_null guide UMI totals are inconsistent")
    if channel_null.retained_cell_count <= 0 or channel_null.total_guide_umi <= 0:
        raise GuideCallingError("channel_null counts must be positive")
    target_map, ntc = _validated_annotations(
        guide_to_target,
        non_targeting_guides,
        hypothesis_count=channel_null.hypothesis_count,
    )
    threshold = _validated_fdr(fdr)

    previous_cell_index = -1
    observed_cells = 0
    observed_umi = 0
    observed_guide_umi = np.zeros(channel_null.hypothesis_count, dtype=np.int64)

    for batch in retained_hash_singlet_batches:
        if not isinstance(batch, SparseGuideBatch):
            raise GuideCallingError(
                "retained_hash_singlet_batches must contain SparseGuideBatch values"
            )
        if batch.cell_count:
            if int(batch.cell_indices[0]) <= previous_cell_index:
                raise GuideCallingError(
                    "retained cell indices must be globally strictly increasing"
                )
            previous_cell_index = int(batch.cell_indices[-1])
        if np.any(batch.guide_indices < 0) or np.any(
            batch.guide_indices >= channel_null.hypothesis_count
        ):
            raise GuideCallingError("guide index is outside the canonical guide axis")

        np.add.at(observed_guide_umi, batch.guide_indices, batch.observed_umi)
        observed_cells += batch.cell_count
        observed_umi += int(batch.observed_umi.sum(dtype=np.int64))
        for cell_index, guides, counts in batch.iter_cells():
            yield _call_sparse_cell_validated(
                cell_index,
                guides,
                counts,
                channel_probability=probability,
                guide_to_target=target_map,
                non_targeting_guides=ntc,
                hypothesis_count=channel_null.hypothesis_count,
                fdr=threshold,
            )

    if observed_cells != channel_null.retained_cell_count:
        raise GuideCallingError(
            "guide-calling pass retained-cell count differs from null-fit pass"
        )
    if observed_umi != channel_null.total_guide_umi or not np.array_equal(
        observed_guide_umi, guide_umi
    ):
        raise GuideCallingError(
            "guide-calling pass observations differ from channel null-fit pass"
        )
