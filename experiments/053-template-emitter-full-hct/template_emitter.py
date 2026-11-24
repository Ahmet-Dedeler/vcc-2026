"""Self-contained fixed-margin control-template emitter for experiment 053.

The output has exact requested row and column margins. Counts are initially
allocated within each gene according to that gene's support and abundance in
real control cells. A bounded repair then moves counts within columns until
every requested library size is exact.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class TemplateEmitterDiagnostics:
    grand_total: int
    row_l1_before_repair: int
    max_row_error_before_repair: int
    moved_counts: int
    template_nnz: int
    output_nnz: int
    retained_template_nnz: int
    removed_template_nnz: int
    introduced_nnz: int

    def to_dict(self) -> dict[str, int]:
        return asdict(self)


def as_integer_csr(matrix) -> sp.csr_matrix:
    """Return canonical int64 CSR after validating raw non-negative counts."""
    x = sp.csr_matrix(matrix)
    x.sum_duplicates()
    x.eliminate_zeros()
    if x.ndim != 2 or x.shape[0] == 0 or x.shape[1] == 0:
        raise ValueError(f"matrix must be non-empty and 2-D, got {x.shape}")
    if not np.all(np.isfinite(x.data)) or np.any(x.data < 0):
        raise ValueError("counts must be finite and non-negative")
    rounded = np.rint(x.data)
    if not np.array_equal(x.data, rounded):
        raise ValueError("counts must be integer-valued")
    x.data = rounded.astype(np.int64)
    x.sort_indices()
    return x


def largest_remainder(
    total: int,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Hamilton-apportion ``total`` integer units over non-negative weights."""
    total = int(total)
    w = np.asarray(weights, dtype=np.float64)
    if total < 0:
        raise ValueError("total must be non-negative")
    if w.ndim != 1 or np.any(w < 0) or not np.all(np.isfinite(w)):
        raise ValueError("weights must be a finite non-negative vector")
    if total == 0:
        return np.zeros(w.size, dtype=np.int64)
    weight_sum = float(w.sum())
    if weight_sum <= 0:
        raise ValueError("positive total requires a positive weight")

    ideal = total * w / weight_sum
    result = np.floor(ideal).astype(np.int64)
    remainder = total - int(result.sum())
    if remainder:
        ties = rng.random(w.size)
        order = np.lexsort((ties, -(ideal - result)))
        result[order[:remainder]] += 1
    return result


def integer_margin(
    weights: np.ndarray,
    total: int,
    *,
    seed: int,
) -> np.ndarray:
    """Create an exact integer margin proportional to ``weights``."""
    return largest_remainder(total, weights, np.random.default_rng(seed))


def scorer_pseudobulk(counts, target_sum: float = 50_000.0) -> np.ndarray:
    """Compute the scorer's bulk log-normalized expression vector."""
    x = sp.csr_matrix(counts)
    columns = np.asarray(x.sum(axis=0), dtype=np.float64).ravel()
    total = float(columns.sum())
    if total <= 0:
        return np.zeros(x.shape[1], dtype=np.float64)
    return np.log1p(float(target_sum) * columns / total)


def _spread_under_caps(
    total: int,
    capacities: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Spread integer units approximately evenly under integer capacities."""
    caps = np.asarray(capacities, dtype=np.int64)
    if total < 0 or total > int(caps.sum()):
        raise ValueError("requested transfer exceeds available capacity")
    result = np.zeros(caps.size, dtype=np.int64)
    left = int(total)
    while left:
        active = np.flatnonzero(result < caps)
        if active.size == 0:
            raise RuntimeError("capacity accounting failed")
        tranche = max(1, left // active.size)
        addition = np.minimum(caps[active] - result[active], tranche)
        used = int(addition.sum())
        if used > left:
            selected = rng.permutation(active.size)[:left]
            addition[:] = 0
            addition[selected] = 1
            used = left
        result[active] += addition
        left -= used
    return result


def _move_with_priority(
    counts: np.ndarray,
    template_support: np.ndarray,
    donor: int,
    recipient: int,
    amount: int,
    rng: np.random.Generator,
) -> int:
    """Move counts donor-to-recipient within genes, preserving columns."""
    left = int(amount)
    if left <= 0:
        return 0

    priorities = (
        lambda: (counts[donor] > 1) & (counts[recipient] > 0),
        lambda: (counts[donor] > 1) & template_support[recipient],
        lambda: (counts[donor] > 0) & (counts[recipient] > 0),
        lambda: (counts[donor] > 0) & template_support[recipient],
        lambda: counts[donor] > 1,
        lambda: counts[donor] > 0,
    )
    for make_mask in priorities:
        if left == 0:
            break
        columns = np.flatnonzero(make_mask())
        if columns.size == 0:
            continue
        keep = 1 if np.all(counts[donor, columns] > 1) else 0
        capacities = counts[donor, columns].astype(np.int64) - keep
        available = int(capacities.sum())
        if available <= 0:
            continue
        take = min(left, available)
        moved = _spread_under_caps(take, capacities, rng)
        nonzero = moved > 0
        columns = columns[nonzero]
        moved = moved[nonzero].astype(counts.dtype, copy=False)
        counts[donor, columns] -= moved
        counts[recipient, columns] += moved
        left -= take

    if left:
        raise RuntimeError(
            f"could not move {left} of {amount} counts from row {donor} to {recipient}"
        )
    return int(amount)


def _repair_row_margins(
    counts: np.ndarray,
    desired_rows: np.ndarray,
    template_support: np.ndarray,
    rng: np.random.Generator,
) -> int:
    """Repair row margins using only within-column transfers."""
    current = counts.sum(axis=1, dtype=np.int64)
    delta = current - desired_rows
    if int(delta.sum()) != 0:
        raise RuntimeError("row and column grand totals disagree")

    donors = np.flatnonzero(delta > 0)
    recipients = np.flatnonzero(delta < 0)
    inverse_depth = 1.0 / desired_rows.astype(np.float64)
    donors = donors[np.argsort(inverse_depth[donors], kind="stable")]
    recipients = recipients[np.argsort(inverse_depth[recipients], kind="stable")]

    donor_pos = 0
    recipient_pos = 0
    moved_total = 0
    while donor_pos < donors.size and recipient_pos < recipients.size:
        donor = int(donors[donor_pos])
        recipient = int(recipients[recipient_pos])
        amount = min(int(delta[donor]), int(-delta[recipient]))
        moved_total += _move_with_priority(
            counts,
            template_support,
            donor,
            recipient,
            amount,
            rng,
        )
        delta[donor] -= amount
        delta[recipient] += amount
        if delta[donor] == 0:
            donor_pos += 1
        if delta[recipient] == 0:
            recipient_pos += 1

    if np.any(delta):
        raise RuntimeError(f"row-margin repair left residuals: {delta}")
    return moved_total


def emit_template_fixed_margins(
    template,
    *,
    column_totals: np.ndarray,
    row_margins: np.ndarray,
    seed: int,
) -> tuple[sp.csr_matrix, TemplateEmitterDiagnostics]:
    """Emit one sparse count block with exact rows and columns."""
    source = as_integer_csr(template)
    n_cells, n_genes = source.shape
    rng = np.random.default_rng(seed)

    rows = np.asarray(row_margins)
    columns = np.asarray(column_totals)
    if (
        rows.shape != (n_cells,)
        or np.any(rows <= 0)
        or not np.all(np.isfinite(rows))
        or not np.array_equal(rows, np.rint(rows))
    ):
        raise ValueError(
            f"row_margins must be positive integer-valued length {n_cells}"
        )
    if (
        columns.shape != (n_genes,)
        or np.any(columns < 0)
        or not np.all(np.isfinite(columns))
        or not np.array_equal(columns, np.rint(columns))
    ):
        raise ValueError(
            f"column_totals must be non-negative integer-valued length {n_genes}"
        )
    rows = np.rint(rows).astype(np.int64)
    columns = np.rint(columns).astype(np.int64)
    if int(rows.sum()) != int(columns.sum()):
        raise ValueError("row and column margins have different grand totals")
    if np.any(rows > 1_000_000):
        raise ValueError("a requested row margin exceeds 1,000,000 counts")
    if int(columns.max(initial=0)) > np.iinfo(np.int32).max:
        raise ValueError("a column total exceeds int32 capacity")

    output = np.zeros((n_cells, n_genes), dtype=np.int32)
    source_csc = source.tocsc()
    for gene in np.flatnonzero(columns):
        lo, hi = source_csc.indptr[gene : gene + 2]
        support_rows = source_csc.indices[lo:hi]
        weights = source_csc.data[lo:hi]
        if support_rows.size == 0:
            support_rows = np.arange(n_cells, dtype=np.int64)
            weights = rows
        allocation = largest_remainder(int(columns[gene]), weights, rng)
        nonzero = allocation > 0
        output[support_rows[nonzero], gene] = allocation[nonzero].astype(np.int32)

    before = output.sum(axis=1, dtype=np.int64) - rows
    template_support = source.astype(bool).toarray()
    moved = _repair_row_margins(output, rows, template_support, rng)

    got_rows = output.sum(axis=1, dtype=np.int64)
    got_columns = output.sum(axis=0, dtype=np.int64)
    if not np.array_equal(got_rows, rows):
        raise RuntimeError("internal error: row margins are not exact")
    if not np.array_equal(got_columns, columns):
        raise RuntimeError("internal error: column margins are not exact")
    if np.any(output < 0):
        raise RuntimeError("internal error: negative count")

    output_support = output > 0
    retained = int(np.count_nonzero(output_support & template_support))
    removed = int(np.count_nonzero(template_support & ~output_support))
    introduced = int(np.count_nonzero(output_support & ~template_support))
    result = sp.csr_matrix(output)
    result.eliminate_zeros()
    result.sort_indices()

    diagnostics = TemplateEmitterDiagnostics(
        grand_total=int(columns.sum()),
        row_l1_before_repair=int(np.abs(before).sum()),
        max_row_error_before_repair=int(np.abs(before).max(initial=0)),
        moved_counts=int(moved),
        template_nnz=int(source.nnz),
        output_nnz=int(result.nnz),
        retained_template_nnz=retained,
        removed_template_nnz=removed,
        introduced_nnz=introduced,
    )
    return result, diagnostics
