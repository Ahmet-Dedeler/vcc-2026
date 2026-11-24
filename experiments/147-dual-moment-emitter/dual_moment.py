"""Exact-margin emitter with independently targeted pooled and arithmetic CPM moments."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import scipy.sparse as sp


CPM_SUM = 1_000_000.0
BULK_SUM = 50_000.0
MAX_COUNTS_PER_CELL = 1_000_000


@dataclass(frozen=True)
class DualMomentDiagnostics:
    grand_total: int
    low_total: int
    high_total: int
    low_cells: int
    high_cells: int
    depth_contrast: float
    blend_lambda: float
    quasi_spread: float
    projected_low_gene_mass_l1: float
    clipped_gene_fraction: float
    desired_cpm_relative_l1: float
    feasible_cpm_relative_l1: float
    max_row_total: int
    output_nnz: int

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


def as_integer_csr(matrix) -> sp.csr_matrix:
    """Return canonical int64 CSR after validating raw counts."""
    x = sp.csr_matrix(matrix)
    x.sum_duplicates()
    x.eliminate_zeros()
    x.sort_indices()
    if x.ndim != 2 or not x.shape[0] or not x.shape[1]:
        raise ValueError(f"counts must be a non-empty matrix, got {x.shape}")
    if not np.all(np.isfinite(x.data)) or np.any(x.data < 0):
        raise ValueError("counts must be finite and non-negative")
    rounded = np.rint(x.data)
    if not np.array_equal(x.data, rounded):
        raise ValueError("counts must be integer-valued")
    x.data = rounded.astype(np.int64)
    return x


def scorer_pseudobulk(counts) -> np.ndarray:
    """Exact pooled profile used by VCC PDS and the MSE plug-in term."""
    x = sp.csr_matrix(counts)
    columns = np.asarray(x.sum(axis=0), dtype=np.float64).ravel()
    total = float(columns.sum())
    if total <= 0:
        raise ValueError("count block has zero total")
    return np.log1p(BULK_SUM * columns / total)


def arithmetic_mean_cpm(counts) -> np.ndarray:
    """Arithmetic mean of per-cell CPM, used by the DE fold-change path."""
    x = sp.csr_matrix(counts).astype(np.float64)
    rows = np.asarray(x.sum(axis=1), dtype=np.float64).ravel()
    if np.any(rows <= 0):
        raise ValueError("all cells must have positive depth")
    row_index = np.repeat(np.arange(x.shape[0]), np.diff(x.indptr))
    scaled = x.copy()
    scaled.data *= CPM_SUM / rows[row_index]
    return np.asarray(scaled.mean(axis=0), dtype=np.float64).ravel()


def blend_cpm(base: np.ndarray, target: np.ndarray, blend_lambda: float) -> np.ndarray:
    """Convexly blend two positive-composition arithmetic CPM profiles."""
    if not 0.0 <= blend_lambda <= 1.0:
        raise ValueError("blend_lambda must be in [0, 1]")
    left = np.asarray(base, dtype=np.float64)
    right = np.asarray(target, dtype=np.float64)
    if (
        left.shape != right.shape
        or left.ndim != 1
        or np.any(left < 0)
        or np.any(right < 0)
        or not np.all(np.isfinite(left))
        or not np.all(np.isfinite(right))
    ):
        raise ValueError("CPM profiles must be same-shape finite non-negative vectors")
    if float(left.sum()) <= 0 or float(right.sum()) <= 0:
        raise ValueError("CPM profiles must have positive mass")
    left = CPM_SUM * left / left.sum()
    right = CPM_SUM * right / right.sum()
    mixed = (1.0 - blend_lambda) * left + blend_lambda * right
    return CPM_SUM * mixed / mixed.sum()


def _balanced_partition(total: int, parts: int, rng: np.random.Generator) -> np.ndarray:
    """Positive integer margins differing by at most one."""
    if parts <= 0 or total < parts:
        raise ValueError("cannot split total into positive requested parts")
    base, remainder = divmod(int(total), int(parts))
    result = np.full(parts, base, dtype=np.int64)
    if remainder:
        result[rng.permutation(parts)[:remainder]] += 1
    return result


def two_depth_row_margins(
    grand_total: int,
    n_cells: int,
    depth_contrast: float,
    *,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Create exact balanced low/high row margins at a fixed depth contrast."""
    if n_cells < 2 or n_cells % 2:
        raise ValueError("two-depth emission requires a positive even cell count")
    if not 0.0 < depth_contrast < 1.0:
        raise ValueError("depth_contrast must lie strictly between zero and one")
    if grand_total <= n_cells:
        raise ValueError("grand total is too small for positive row margins")

    low_cells = n_cells // 2
    high_cells = n_cells - low_cells
    low_fraction = low_cells * (1.0 - depth_contrast) / n_cells
    low_total = int(np.rint(grand_total * low_fraction))
    low_total = min(max(low_total, low_cells), grand_total - high_cells)
    high_total = grand_total - low_total

    rng = np.random.default_rng(seed)
    low = _balanced_partition(low_total, low_cells, rng)
    high = _balanced_partition(high_total, high_cells, rng)
    rows = np.concatenate((low, high))
    mask = np.zeros(n_cells, dtype=bool)
    mask[:low_cells] = True
    if int(rows.sum()) != grand_total:
        raise RuntimeError("internal row-margin accounting error")
    if int(rows.max(initial=0)) > MAX_COUNTS_PER_CELL:
        raise ValueError("a requested row exceeds the VCC count cap")
    return rows, mask


def _project_box_sum(
    raw: np.ndarray,
    caps: np.ndarray,
    total: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Euclidean projection onto 0<=x<=caps, sum(x)=total, then exact rounding."""
    values = np.asarray(raw, dtype=np.float64)
    upper = np.asarray(caps, dtype=np.int64)
    if values.shape != upper.shape or values.ndim != 1:
        raise ValueError("raw and caps must be same-shape vectors")
    if np.any(upper < 0) or not np.all(np.isfinite(values)):
        raise ValueError("invalid projection inputs")
    if total < 0 or total > int(upper.sum()):
        raise ValueError("projection total lies outside box capacity")

    lo = float(np.min(values - upper, initial=0.0)) - 1.0
    hi = float(np.max(values, initial=0.0)) + 1.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        current = np.clip(values - mid, 0.0, upper)
        if float(current.sum()) > total:
            lo = mid
        else:
            hi = mid
    projected = np.clip(values - hi, 0.0, upper)

    integer = np.floor(projected).astype(np.int64)
    remaining = int(total - integer.sum())
    if remaining:
        eligible = np.flatnonzero(integer < upper)
        if remaining > eligible.size:
            raise RuntimeError("bounded rounding lacks enough eligible coordinates")
        fractions = projected[eligible] - integer[eligible]
        tie = rng.random(eligible.size)
        order = np.lexsort((tie, -fractions))
        integer[eligible[order[:remaining]]] += 1

    if int(integer.sum()) != total or np.any(integer < 0) or np.any(integer > upper):
        raise RuntimeError("bounded integer projection failed")
    return projected, integer


def solve_low_group_totals(
    column_totals: np.ndarray,
    desired_cpm: np.ndarray,
    row_margins: np.ndarray,
    low_mask: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve and project low-depth per-gene totals for the desired arithmetic CPM."""
    columns = np.asarray(column_totals, dtype=np.int64)
    desired = np.asarray(desired_cpm, dtype=np.float64)
    rows = np.asarray(row_margins, dtype=np.int64)
    low = np.asarray(low_mask, dtype=bool)
    if columns.shape != desired.shape or columns.ndim != 1:
        raise ValueError("column_totals and desired_cpm must be same-shape vectors")
    if rows.ndim != 1 or low.shape != rows.shape or not np.any(low) or np.all(low):
        raise ValueError("invalid row margins or low-group mask")
    if int(columns.sum()) != int(rows.sum()):
        raise ValueError("row and column grand totals differ")

    desired = CPM_SUM * desired / desired.sum()
    n_cells = rows.size
    low_total = int(rows[low].sum())
    high_total = int(rows[~low].sum())
    r_low = low_total / int(low.sum())
    r_high = high_total / int((~low).sum())
    denominator = 1.0 / r_low - 1.0 / r_high
    if denominator <= 0:
        raise ValueError("low group must have smaller average depth than high group")
    raw = (n_cells * desired / CPM_SUM - columns / r_high) / denominator
    projected, integer = _project_box_sum(
        raw,
        columns,
        low_total,
        np.random.default_rng(seed),
    )
    feasible = CPM_SUM / n_cells * (
        integer / r_low + (columns - integer) / r_high
    )
    feasible = CPM_SUM * feasible / feasible.sum()
    return raw, integer, feasible


def _spread_under_caps(
    total: int,
    capacities: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Spread integer transfers across available genes under exact capacities."""
    caps = np.asarray(capacities, dtype=np.int64)
    if total < 0 or total > int(caps.sum()):
        raise ValueError("requested transfer exceeds available capacity")
    out = np.zeros(caps.size, dtype=np.int64)
    left = int(total)
    while left:
        active = np.flatnonzero(out < caps)
        if not active.size:
            raise RuntimeError("transfer capacity accounting failed")
        tranche = max(1, left // active.size)
        add = np.minimum(caps[active] - out[active], tranche)
        used = int(add.sum())
        if used > left:
            chosen = rng.permutation(active.size)[:left]
            add[:] = 0
            add[chosen] = 1
            used = left
        out[active] += add
        left -= used
    return out


def _repair_rows(
    block: np.ndarray,
    desired_rows: np.ndarray,
    rng: np.random.Generator,
) -> None:
    """Repair row totals using within-gene transfers, preserving every column."""
    delta = block.sum(axis=1, dtype=np.int64) - desired_rows
    if int(delta.sum()) != 0:
        raise RuntimeError("row and column totals disagree before repair")
    donors = list(map(int, np.flatnonzero(delta > 0)))
    recipients = list(map(int, np.flatnonzero(delta < 0)))
    di = ri = 0
    while di < len(donors) and ri < len(recipients):
        donor = donors[di]
        recipient = recipients[ri]
        amount = min(int(delta[donor]), int(-delta[recipient]))
        columns = np.flatnonzero(block[donor] > 0)
        if not columns.size:
            raise RuntimeError("positive row residual has no movable counts")
        moved = _spread_under_caps(amount, block[donor, columns], rng)
        use = moved > 0
        columns = columns[use]
        moved = moved[use].astype(block.dtype, copy=False)
        block[donor, columns] -= moved
        block[recipient, columns] += moved
        delta[donor] -= amount
        delta[recipient] += amount
        if delta[donor] == 0:
            di += 1
        if delta[recipient] == 0:
            ri += 1
    if np.any(delta):
        raise RuntimeError("row repair left nonzero residuals")


def _allocate_group(
    column_totals: np.ndarray,
    row_margins: np.ndarray,
    *,
    quasi_spread: float,
    seed: int,
) -> np.ndarray:
    """Allocate exact group columns with quasi-Poisson spread and exact rows."""
    columns = np.asarray(column_totals, dtype=np.int64)
    rows = np.asarray(row_margins, dtype=np.int64)
    if not 0.0 <= quasi_spread <= 1.0:
        raise ValueError("quasi_spread must be in [0, 1]")
    if np.any(columns < 0) or np.any(rows <= 0):
        raise ValueError("invalid group margins")
    if int(columns.sum()) != int(rows.sum()):
        raise ValueError("group row and column totals differ")

    n_cells = rows.size
    genes = columns.size
    rng = np.random.default_rng(seed)
    random_totals = np.floor(columns * quasi_spread).astype(np.int64)
    even_totals = columns - random_totals
    block = np.zeros((n_cells, genes), dtype=np.int32)

    nonzero = np.flatnonzero(even_totals)
    if nonzero.size:
        base, remainder = np.divmod(even_totals[nonzero], n_cells)
        offsets = np.concatenate(([0], np.cumsum(remainder)[:-1])) % n_cells
        positions = (
            np.arange(n_cells, dtype=np.int64)[:, None] - offsets[None, :]
        ) % n_cells
        block[:, nonzero] = (
            base[None, :] + (positions < remainder[None, :])
        ).astype(np.int32)

    probabilities = rows.astype(np.float64) / rows.sum()
    for gene in np.flatnonzero(random_totals):
        block[:, gene] += rng.multinomial(
            int(random_totals[gene]),
            probabilities,
        ).astype(np.int32)

    _repair_rows(block, rows, rng)
    if not np.array_equal(block.sum(axis=0, dtype=np.int64), columns):
        raise RuntimeError("group allocation changed column totals")
    if not np.array_equal(block.sum(axis=1, dtype=np.int64), rows):
        raise RuntimeError("group allocation missed row totals")
    if np.any(block < 0):
        raise RuntimeError("group allocation produced negative counts")
    return block


def emit_dual_moment(
    incumbent_counts,
    direction_counts,
    *,
    blend_lambda: float,
    depth_contrast: float = 0.5,
    quasi_spread: float = 0.25,
    seed: int = 14700,
) -> tuple[sp.csr_matrix, DualMomentDiagnostics]:
    """Preserve incumbent pooled counts while targeting another arithmetic CPM profile."""
    incumbent = as_integer_csr(incumbent_counts)
    direction = as_integer_csr(direction_counts)
    if incumbent.shape != direction.shape:
        raise ValueError("incumbent and direction blocks must have identical shape")
    n_cells, _ = incumbent.shape
    if n_cells % 2:
        raise ValueError("dual-moment emitter requires an even cell count")

    columns = np.asarray(incumbent.sum(axis=0), dtype=np.int64).ravel()
    grand_total = int(columns.sum())
    rows, low_mask = two_depth_row_margins(
        grand_total,
        n_cells,
        depth_contrast,
        seed=seed + 1,
    )
    base_cpm = arithmetic_mean_cpm(incumbent)
    direction_cpm = arithmetic_mean_cpm(direction)
    desired_cpm = blend_cpm(base_cpm, direction_cpm, blend_lambda)
    raw_low, low_columns, feasible_cpm = solve_low_group_totals(
        columns,
        desired_cpm,
        rows,
        low_mask,
        seed=seed + 2,
    )
    high_columns = columns - low_columns

    low_block = _allocate_group(
        low_columns,
        rows[low_mask],
        quasi_spread=quasi_spread,
        seed=seed + 3,
    )
    high_block = _allocate_group(
        high_columns,
        rows[~low_mask],
        quasi_spread=quasi_spread,
        seed=seed + 4,
    )
    dense = np.vstack((low_block, high_block))
    emitted = sp.csr_matrix(dense)
    emitted.eliminate_zeros()
    emitted.sort_indices()

    got_rows = np.asarray(emitted.sum(axis=1), dtype=np.int64).ravel()
    got_columns = np.asarray(emitted.sum(axis=0), dtype=np.int64).ravel()
    if not np.array_equal(got_rows, rows):
        raise RuntimeError("final row margins are not exact")
    if not np.array_equal(got_columns, columns):
        raise RuntimeError("final column margins are not exact")
    if not np.array_equal(scorer_pseudobulk(emitted), scorer_pseudobulk(incumbent)):
        raise RuntimeError("final scorer pseudobulk changed")
    if emitted.data.size and (
        np.any(emitted.data < 0)
        or not np.array_equal(emitted.data, np.rint(emitted.data))
    ):
        raise RuntimeError("final count matrix is not non-negative integer")
    if int(got_rows.max(initial=0)) > MAX_COUNTS_PER_CELL:
        raise RuntimeError("final row exceeds the VCC count cap")

    achieved_cpm = arithmetic_mean_cpm(emitted)
    desired_denominator = float(np.abs(desired_cpm).sum())
    feasible_denominator = float(np.abs(feasible_cpm).sum())
    raw_scale = max(1.0, float(np.abs(raw_low).sum()))
    diagnostics = DualMomentDiagnostics(
        grand_total=grand_total,
        low_total=int(rows[low_mask].sum()),
        high_total=int(rows[~low_mask].sum()),
        low_cells=int(low_mask.sum()),
        high_cells=int((~low_mask).sum()),
        depth_contrast=float(depth_contrast),
        blend_lambda=float(blend_lambda),
        quasi_spread=float(quasi_spread),
        projected_low_gene_mass_l1=float(
            np.abs(low_columns - raw_low).sum() / raw_scale
        ),
        clipped_gene_fraction=float(
            np.mean((raw_low < 0.0) | (raw_low > columns))
        ),
        desired_cpm_relative_l1=float(
            np.abs(achieved_cpm - desired_cpm).sum() / desired_denominator
        ),
        feasible_cpm_relative_l1=float(
            np.abs(achieved_cpm - feasible_cpm).sum() / feasible_denominator
        ),
        max_row_total=int(got_rows.max(initial=0)),
        output_nnz=int(emitted.nnz),
    )
    return emitted, diagnostics
