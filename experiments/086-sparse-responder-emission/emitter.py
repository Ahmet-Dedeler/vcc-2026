#!/usr/bin/env python3
"""Exact-margin sparse-responder count emission.

The construction keeps a deterministic subset of emitted cells as close as
possible to real control templates. The complementary responder subset absorbs
the exact residual needed to reproduce a pre-existing integer pseudobulk.

This module is intentionally block-local: one perturbation (normally
400 x 18,533) is materialized at a time.
"""

from __future__ import annotations

import hashlib
import math
import sys
from dataclasses import asdict, dataclass

import numpy as np
import scipy.sparse as sp


BULK_TARGET_SUM = 50_000.0
CPM_TARGET_SUM = 1_000_000.0
MAX_COUNTS_PER_CELL = 1_000_000


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    return hashlib.sha256(array.view(np.uint8)).hexdigest()


def stable_seed(*parts: object) -> int:
    text = "|".join(map(str, parts))
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def _integer_csr(values) -> sp.csr_matrix:
    matrix = sp.csr_matrix(values)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if matrix.ndim != 2 or min(matrix.shape) < 1:
        raise ValueError(f"control template must be non-empty 2-D, got {matrix.shape}")
    if np.any(matrix.data < 0) or not np.isfinite(matrix.data).all():
        raise ValueError("control template contains negative or non-finite counts")
    rounded = np.rint(matrix.data)
    if not np.array_equal(matrix.data, rounded):
        raise ValueError("control template must contain integer raw counts")
    matrix.data = rounded.astype(np.int64)
    return matrix


def _integer_totals(values: np.ndarray, genes: int) -> np.ndarray:
    totals = np.asarray(values)
    if totals.shape != (genes,):
        raise ValueError(
            f"column_totals must have shape ({genes},), got {totals.shape}"
        )
    if np.any(totals < 0) or not np.isfinite(totals).all():
        raise ValueError("column_totals must be finite and non-negative")
    rounded = np.rint(totals)
    if not np.array_equal(totals, rounded):
        raise ValueError("column_totals must be integer-valued")
    return rounded.astype(np.int64)


def deterministic_responder_indices(
    cells: int,
    responder_fraction: float,
    *,
    key: str,
    seed: int,
) -> np.ndarray:
    """Choose one stable cell subset without relying on NumPy RNG behavior."""
    if cells < 2:
        raise ValueError("sparse response requires at least two cells")
    if not (0.0 < responder_fraction < 1.0):
        raise ValueError("responder_fraction must be strictly between 0 and 1")
    responders = int(math.ceil(cells * float(responder_fraction)))
    responders = min(max(responders, 1), cells - 1)
    priorities = [
        hashlib.sha256(f"{seed}|{key}|responder|{row}".encode()).digest()
        for row in range(cells)
    ]
    order = sorted(range(cells), key=lambda row: (priorities[row], row))
    return np.asarray(sorted(order[:responders]), dtype=np.int64)


def _proportional_cap(
    counts: np.ndarray,
    keep_total: int,
    *,
    tie_offset: int,
) -> np.ndarray:
    """Reduce one integer column to an exact total without creating support."""
    values = np.asarray(counts, dtype=np.int64)
    current = int(values.sum())
    keep_total = int(keep_total)
    if keep_total < 0 or keep_total > current:
        raise ValueError("cap total lies outside [0, current]")
    if keep_total == current:
        return values.copy()
    if keep_total == 0:
        return np.zeros_like(values)
    numerators = values * keep_total
    out = numerators // current
    left = keep_total - int(out.sum())
    if left:
        remainder = numerators % current
        eligible = np.flatnonzero(out < values)
        if left > len(eligible):
            raise RuntimeError("Hamilton cap has insufficient eligible rows")
        tie = (eligible - int(tie_offset)) % max(len(values), 1)
        order = np.lexsort((tie, -remainder[eligible]))
        out[eligible[order[:left]]] += 1
    if int(out.sum()) != keep_total or np.any(out > values) or np.any(out < 0):
        raise RuntimeError("proportional cap violated an integer invariant")
    return out


def _allocate_under_row_deficits(
    total: int,
    deficits: np.ndarray,
    preferred_rows: np.ndarray,
    *,
    tie_offset: int,
) -> np.ndarray:
    """Allocate one positive-gene increment without changing any row margin."""
    total = int(total)
    caps = np.asarray(deficits, dtype=np.int64)
    if total < 0 or total > int(caps.sum()):
        raise ValueError("positive increment exceeds available row deficits")
    out = np.zeros_like(caps)
    left = total
    preferred = np.asarray(preferred_rows, dtype=np.int64)
    fallback = np.setdiff1d(
        np.arange(len(caps), dtype=np.int64), preferred, assume_unique=True
    )
    for rows in (preferred, fallback):
        if left == 0:
            break
        available = int(caps[rows].sum())
        if available <= 0:
            continue
        take = min(left, available)
        local_caps = caps[rows]
        # Proportional-to-deficit Hamilton allocation. It uses exactly `take`
        # units, never exceeds a row's missing-count capacity, and minimizes
        # avoidable row-composition distortion at this step.
        numerators = local_caps * take
        local = numerators // available
        remainder = take - int(local.sum())
        if remainder:
            eligible = np.flatnonzero(local < local_caps)
            fractions = numerators % available
            tie = (eligible - int(tie_offset)) % max(len(rows), 1)
            order = np.lexsort((tie, -fractions[eligible]))
            local[eligible[order[:remainder]]] += 1
        out[rows] += local
        left -= take
    if left or int(out.sum()) != total or np.any(out > caps):
        raise RuntimeError("row-deficit allocation failed")
    return out


@dataclass(frozen=True)
class SparseResponderDiagnostics:
    cells: int
    genes: int
    responder_count: int
    nonresponder_count: int
    requested_responder_fraction: float
    actual_responder_fraction: float
    responder_indices_sha256: str
    input_column_totals_sha256: str
    output_column_totals_sha256: str
    grand_total: int
    max_abs_column_error: int
    output_nnz: int
    row_sum_min: int
    row_sum_median: float
    row_sum_max: int
    row_sum_cv: float
    maximum_row_share: float
    baseline_row_margins_preserved: bool
    baseline_column_totals_sha256: str
    nonresponder_template_counts: int
    nonresponder_retained_counts: int
    nonresponder_clipped_counts: int
    nonresponder_clipped_mass_fraction: float
    nonresponder_clipped_genes: int
    scored_nonresponder_template_counts: int
    scored_nonresponder_clipped_counts: int
    scored_nonresponder_clipped_mass_fraction: float
    scored_nonresponder_clipped_genes: int
    unchanged_nonresponder_rows: int
    unchanged_nonresponder_row_fraction: float
    nonresponder_absolute_change_counts: int
    scored_nonresponder_absolute_change_counts: int
    absolute_change_mass_in_responders_fraction: float
    positive_effect_genes: int
    negative_effect_genes: int
    scored_negative_effect_genes: int
    representative_subset_screening_fraction_all: float
    representative_subset_screening_fraction_scored: float
    selected_subset_exact_control_feasible_all: bool
    selected_subset_exact_control_feasible_scored: bool

    def to_dict(self) -> dict[str, int | float | bool | str]:
        return asdict(self)


def _minimum_responder_fraction(
    target: np.ndarray,
    control: np.ndarray,
    mask: np.ndarray,
) -> float:
    active = mask & (control > 0) & (target < control)
    if not np.any(active):
        return 0.0
    required = 1.0 - target[active].astype(np.float64) / control[active]
    return float(np.max(required))


def emit_sparse_responders(
    control_baseline,
    *,
    column_totals: np.ndarray,
    responder_fraction: float,
    key: str,
    seed: int = 0,
    control_reference_totals: np.ndarray | None = None,
    scored_gene_mask: np.ndarray | None = None,
    max_counts_per_cell: int = MAX_COUNTS_PER_CELL,
    minimum_counts_per_cell: int = 1,
) -> tuple[sp.csr_matrix, SparseResponderDiagnostics]:
    """Emit one exact-margin perturbation block.

    ``control_baseline`` must already be an integer control-like realization at
    the same grand total as ``column_totals``. Negative gene deltas are removed
    from responders first. Positive deltas are then allocated into the exact row
    deficits those removals created. Thus every row margin stays fixed.

    If responders do not carry enough counts for a negative gene, the
    mathematically unavoidable remainder is removed from nonresponders and the
    corresponding positive mass later returns to nonresponders. That leakage is
    measured explicitly.

    The output is deterministic for fixed inputs, ``key``, and ``seed``.
    """
    template = _integer_csr(control_baseline)
    cells, genes = template.shape
    totals = _integer_totals(column_totals, genes)
    if int(totals.sum()) <= 0:
        raise ValueError("column_totals sum to zero")
    baseline_columns = np.asarray(template.sum(axis=0), dtype=np.int64).ravel()
    if int(baseline_columns.sum()) != int(totals.sum()):
        raise ValueError(
            "control baseline and target columns must have the same grand total"
        )
    if max_counts_per_cell < 1 or minimum_counts_per_cell < 1:
        raise ValueError("row-count limits must be positive")
    if minimum_counts_per_cell > max_counts_per_cell:
        raise ValueError("minimum row count exceeds maximum")

    responders = deterministic_responder_indices(
        cells, responder_fraction, key=key, seed=seed
    )
    responder_mask = np.zeros(cells, dtype=bool)
    responder_mask[responders] = True
    nonresponders = np.flatnonzero(~responder_mask)
    actual_fraction = len(responders) / cells

    dense_template = template.toarray().astype(np.int64, copy=False)
    output = dense_template.copy()
    delta = totals - baseline_columns
    negative_genes = np.flatnonzero(delta < 0)
    positive_genes = np.flatnonzero(delta > 0)
    row_deficits = np.zeros(cells, dtype=np.int64)
    nonresponder_removed_by_gene = np.zeros(genes, dtype=np.int64)

    for gene in negative_genes:
        need = int(-delta[gene])
        responder_available = int(output[responders, gene].sum())
        take = min(need, responder_available)
        if take:
            kept = _proportional_cap(
                output[responders, gene],
                responder_available - take,
                tie_offset=stable_seed(seed, key, "remove-responder", int(gene))
                % len(responders),
            )
            removed = output[responders, gene] - kept
            output[responders, gene] = kept
            row_deficits[responders] += removed
        left = need - take
        if left:
            nonresponder_available = int(output[nonresponders, gene].sum())
            if left > nonresponder_available:
                raise RuntimeError(
                    "target margin is negative after exhausting every cell"
                )
            kept = _proportional_cap(
                output[nonresponders, gene],
                nonresponder_available - left,
                tie_offset=stable_seed(seed, key, "remove-nonresponder", int(gene))
                % len(nonresponders),
            )
            removed = output[nonresponders, gene] - kept
            output[nonresponders, gene] = kept
            row_deficits[nonresponders] += removed
            nonresponder_removed_by_gene[gene] = int(removed.sum())

    for gene in positive_genes:
        addition = int(delta[gene])
        allocated = _allocate_under_row_deficits(
            addition,
            row_deficits,
            responders,
            tie_offset=stable_seed(seed, key, "add", int(gene)) % cells,
        )
        output[:, gene] += allocated
        row_deficits -= allocated
    if np.any(row_deficits):
        raise RuntimeError("positive deltas did not exactly fill removed row counts")

    if int(output.max(initial=0)) > np.iinfo(np.int32).max:
        raise ValueError("an emitted count exceeds int32 capacity")
    output = output.astype(np.int32)
    if np.any(output < 0):
        raise RuntimeError("negative emitted count")
    got_columns = output.sum(axis=0, dtype=np.int64)
    if not np.array_equal(got_columns, totals):
        raise RuntimeError("exact column-margin invariant failed")
    row_sums = output.sum(axis=1, dtype=np.int64)
    baseline_rows = dense_template.sum(axis=1, dtype=np.int64)
    if not np.array_equal(row_sums, baseline_rows):
        raise RuntimeError("baseline row margins changed")
    if np.any(row_sums < minimum_counts_per_cell):
        raise ValueError(
            f"emitted row below minimum {minimum_counts_per_cell}: "
            f"{int(row_sums.min())}"
        )
    if np.any(row_sums > max_counts_per_cell):
        raise ValueError(
            f"emitted row exceeds cap {max_counts_per_cell}: {int(row_sums.max())}"
        )

    if scored_gene_mask is None:
        scored = np.ones(genes, dtype=bool)
    else:
        scored = np.asarray(scored_gene_mask, dtype=bool)
        if scored.shape != (genes,):
            raise ValueError("scored_gene_mask has wrong shape")
    if control_reference_totals is None:
        reference = baseline_columns
    else:
        reference = _integer_totals(control_reference_totals, genes)

    nonresponder_before = dense_template[nonresponders]
    nonresponder_after = output[nonresponders].astype(np.int64)
    clipped_by_gene = nonresponder_removed_by_gene
    template_mass = int(nonresponder_before.sum())
    clipped_mass = int(nonresponder_removed_by_gene.sum())
    scored_template_mass = int(nonresponder_before[:, scored].sum())
    scored_clipped_mass = int(nonresponder_removed_by_gene[scored].sum())
    nonresponder_absolute = np.abs(nonresponder_after - nonresponder_before)
    scored_nonresponder_absolute = nonresponder_absolute[:, scored]
    unchanged_rows = int(np.sum(np.all(nonresponder_absolute == 0, axis=1)))

    absolute_change = np.abs(output.astype(np.int64) - dense_template)
    change_total = int(absolute_change.sum())
    responder_change = int(absolute_change[responders].sum())
    negative = (totals < reference) & (reference > 0)
    min_all = _minimum_responder_fraction(totals, reference, np.ones(genes, dtype=bool))
    min_scored = _minimum_responder_fraction(totals, reference, scored)

    result = sp.csr_matrix(output)
    result.eliminate_zeros()
    diagnostics = SparseResponderDiagnostics(
        cells=cells,
        genes=genes,
        responder_count=len(responders),
        nonresponder_count=len(nonresponders),
        requested_responder_fraction=float(responder_fraction),
        actual_responder_fraction=float(actual_fraction),
        responder_indices_sha256=sha256_array(responders),
        input_column_totals_sha256=sha256_array(totals),
        output_column_totals_sha256=sha256_array(got_columns),
        grand_total=int(totals.sum()),
        max_abs_column_error=int(np.max(np.abs(got_columns - totals), initial=0)),
        output_nnz=int(result.nnz),
        row_sum_min=int(row_sums.min()),
        row_sum_median=float(np.median(row_sums)),
        row_sum_max=int(row_sums.max()),
        row_sum_cv=float(
            row_sums.std(ddof=1) / row_sums.mean() if row_sums.mean() else 0.0
        ),
        maximum_row_share=float(row_sums.max() / totals.sum()),
        baseline_row_margins_preserved=True,
        baseline_column_totals_sha256=sha256_array(baseline_columns),
        nonresponder_template_counts=template_mass,
        nonresponder_retained_counts=int(nonresponder_before.sum() - clipped_mass),
        nonresponder_clipped_counts=clipped_mass,
        nonresponder_clipped_mass_fraction=float(
            clipped_mass / template_mass if template_mass else 0.0
        ),
        nonresponder_clipped_genes=int(np.count_nonzero(clipped_by_gene)),
        scored_nonresponder_template_counts=scored_template_mass,
        scored_nonresponder_clipped_counts=scored_clipped_mass,
        scored_nonresponder_clipped_mass_fraction=float(
            scored_clipped_mass / scored_template_mass if scored_template_mass else 0.0
        ),
        scored_nonresponder_clipped_genes=int(
            np.count_nonzero(clipped_by_gene[scored])
        ),
        unchanged_nonresponder_rows=unchanged_rows,
        unchanged_nonresponder_row_fraction=float(
            unchanged_rows / len(nonresponders) if len(nonresponders) else 0.0
        ),
        nonresponder_absolute_change_counts=int(nonresponder_absolute.sum()),
        scored_nonresponder_absolute_change_counts=int(
            scored_nonresponder_absolute.sum()
        ),
        absolute_change_mass_in_responders_fraction=float(
            responder_change / change_total if change_total else 1.0
        ),
        positive_effect_genes=int(len(positive_genes)),
        negative_effect_genes=int(np.count_nonzero(negative)),
        scored_negative_effect_genes=int(np.count_nonzero(negative & scored)),
        representative_subset_screening_fraction_all=min_all,
        representative_subset_screening_fraction_scored=min_scored,
        selected_subset_exact_control_feasible_all=bool(
            np.sum(nonresponder_absolute) == 0
        ),
        selected_subset_exact_control_feasible_scored=bool(
            np.sum(scored_nonresponder_absolute) == 0
        ),
    )
    return result, diagnostics


def scorer_pseudobulk(counts, target_sum: float = BULK_TARGET_SUM) -> np.ndarray:
    matrix = sp.csr_matrix(counts)
    totals = np.asarray(matrix.sum(axis=0), dtype=np.float64).ravel()
    grand = float(totals.sum())
    if grand <= 0:
        raise ValueError("count matrix has zero grand total")
    return np.log1p(float(target_sum) * totals / grand)


def arithmetic_mean_cpm(counts) -> np.ndarray:
    matrix = sp.csr_matrix(counts).astype(np.float64)
    libraries = np.asarray(matrix.sum(axis=1), dtype=np.float64).ravel()
    if np.any(libraries <= 0):
        raise ValueError("all rows need positive library sizes")
    row_of_nnz = np.repeat(np.arange(matrix.shape[0]), np.diff(matrix.indptr))
    matrix.data *= CPM_TARGET_SUM / libraries[row_of_nnz]
    return np.asarray(matrix.mean(axis=0), dtype=np.float64).ravel()


def jackknife_correction(
    counts,
    *,
    target_sum: float = BULK_TARGET_SUM,
    chunk_genes: int = 256,
) -> float:
    """Cell-eval2-style delete-one-cell correction, streamed over genes."""
    matrix = sp.csr_matrix(counts)
    n, genes = matrix.shape
    if n < 2:
        return 0.0
    libraries = np.asarray(matrix.sum(axis=1), dtype=np.float64).ravel()
    grand = float(libraries.sum())
    denominators = grand - libraries
    if grand <= 0 or np.any(denominators <= 0):
        raise ValueError("jackknife requires positive leave-one-out totals")
    columns = np.asarray(matrix.sum(axis=0), dtype=np.float64).ravel()
    total = 0.0
    for start in range(0, genes, chunk_genes):
        stop = min(start + chunk_genes, genes)
        block = matrix[:, start:stop].toarray().astype(np.float64, copy=False)
        leave_one = np.log1p(
            float(target_sum)
            * (columns[start:stop][None, :] - block)
            / denominators[:, None]
        )
        centered = leave_one - leave_one.mean(axis=0, keepdims=True)
        total += float(np.sum(centered * centered))
    return float((n - 1) / n * total)


def process_peak_rss_gib() -> float:
    import resource

    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2
