"""Fixed-margin emission from real control-cell templates.

This is an isolated experiment, not production submission code.

The VCC scorer sees several summaries of the same emitted cells:

* PDS and MSE's plug-in distance see each perturbation's raw column totals,
  normalized once. MSE additionally sees a cell-layout jackknife correction.
* The DE metrics see per-cell CPM ranks and arithmetic per-cell CPM means.

The emitter therefore fixes:

1. row margins to selected control-cell library sizes; and
2. column margins to an integer representation of the model's target profile.

Within those margins, each target gene's counts are first allocated across cells
in proportion to that gene's observed control-template counts. A small repair
then moves counts between rows *within the same gene* until every row margin is
exact. The repair pairs similar-depth donor/recipient rows first, minimizing the
per-count change in arithmetic mean CPM.

The returned matrix is integer CSR. Its scorer pseudobulk is exactly determined
by ``column_totals``; no random sampling can move it.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import scipy.sparse as sp


BULK_TARGET_SUM = 50_000.0
CPM_TARGET_SUM = 1_000_000.0


@dataclass(frozen=True)
class TemplateEmitterDiagnostics:
    """Small diagnostics that can be recorded beside an emitted block."""

    grand_total: int
    row_l1_before_repair: int
    max_row_error_before_repair: int
    moved_counts: int
    template_nnz: int
    output_nnz: int
    retained_template_nnz: int
    removed_template_nnz: int
    introduced_nnz: int
    template_cpm_relative_l1: float

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


def _as_integer_csr(template) -> sp.csr_matrix:
    """Validate raw non-negative integer counts and return canonical int64 CSR."""
    x = sp.csr_matrix(template)
    x.sum_duplicates()
    x.eliminate_zeros()
    if x.ndim != 2 or x.shape[0] == 0 or x.shape[1] == 0:
        raise ValueError(f"template must be a non-empty 2-D matrix, got {x.shape}")
    if not np.all(np.isfinite(x.data)) or np.any(x.data < 0):
        raise ValueError("template counts must be finite and non-negative")
    rounded = np.rint(x.data)
    if not np.array_equal(x.data, rounded):
        raise ValueError("template must contain raw integer counts")
    x.data = rounded.astype(np.int64)
    return x


def _largest_remainder(
    total: int,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Hamilton apportionment: non-negative integers summing exactly to ``total``."""
    total = int(total)
    w = np.asarray(weights, dtype=np.float64)
    if total < 0:
        raise ValueError("total must be non-negative")
    if w.ndim != 1 or np.any(w < 0) or not np.all(np.isfinite(w)):
        raise ValueError("weights must be a finite non-negative vector")
    if total == 0:
        return np.zeros(w.size, dtype=np.int64)
    sw = float(w.sum())
    if sw <= 0:
        raise ValueError("positive total requires at least one positive weight")

    ideal = total * (w / sw)
    out = np.floor(ideal).astype(np.int64)
    left = total - int(out.sum())
    if left:
        # Random tie-breaking is reproducible and prevents one fixed cell index from
        # receiving every equal-remainder count across thousands of genes.
        tie = rng.random(w.size)
        order = np.lexsort((tie, -(ideal - out)))
        out[order[:left]] += 1
    return out


def integer_column_totals(
    profile: np.ndarray,
    grand_total: int,
    *,
    seed: int = 0,
) -> np.ndarray:
    """Closest Hamilton-rounded integer profile at the requested grand total.

    Exact equality to an arbitrary floating profile is generally impossible with
    integer counts. This function defines the intended *integer* pseudobulk. Each
    column differs from its real-valued quota by less than one count.
    """
    p = np.asarray(profile, dtype=np.float64)
    if p.ndim != 1 or np.any(p < 0) or not np.all(np.isfinite(p)):
        raise ValueError("profile must be a finite non-negative vector")
    if float(p.sum()) <= 0:
        raise ValueError("profile sums to zero")
    return _largest_remainder(int(grand_total), p, np.random.default_rng(seed))


def scorer_pseudobulk(counts, target_sum: float = BULK_TARGET_SUM) -> np.ndarray:
    """The exact bulk-log-normalization used by cell-eval2 for PDS/MSE."""
    x = sp.csr_matrix(counts)
    cols = np.asarray(x.sum(axis=0), dtype=np.float64).ravel()
    total = float(cols.sum())
    if total <= 0:
        return np.zeros(x.shape[1], dtype=np.float64)
    return np.log1p(float(target_sum) * cols / total)


def arithmetic_mean_cpm(counts) -> np.ndarray:
    """Arithmetic mean of per-cell CPM, used for VCC DE log2 fold changes."""
    x = sp.csr_matrix(counts).astype(np.float64)
    libs = np.asarray(x.sum(axis=1), dtype=np.float64).ravel()
    if np.any(libs <= 0):
        raise ValueError("all emitted cells must have positive library size")
    row_of_nnz = np.repeat(np.arange(x.shape[0]), np.diff(x.indptr))
    scaled = x.copy()
    scaled.data *= CPM_TARGET_SUM / libs[row_of_nnz]
    return np.asarray(scaled.mean(axis=0), dtype=np.float64).ravel()


def _spread_under_caps(
    total: int,
    capacities: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Spread ``total`` integer units as evenly as possible under integer caps."""
    caps = np.asarray(capacities, dtype=np.int64)
    if total < 0 or total > int(caps.sum()):
        raise ValueError("requested transfer exceeds available capacity")
    out = np.zeros(caps.size, dtype=np.int64)
    left = int(total)
    while left:
        active = np.flatnonzero(out < caps)
        if active.size == 0:
            raise RuntimeError("capacity accounting failed")
        # Give every active gene a common tranche. The random rotation only breaks
        # the final partial-tranche tie and is reproducible.
        tranche = max(1, left // active.size)
        add = np.minimum(caps[active] - out[active], tranche)
        used = int(add.sum())
        if used > left:
            take = rng.permutation(active.size)[:left]
            add[:] = 0
            add[take] = 1
            used = left
        out[active] += add
        left -= used
    return out


def _move_with_priority(
    x: np.ndarray,
    template_support: np.ndarray,
    donor: int,
    recipient: int,
    amount: int,
    rng: np.random.Generator,
) -> int:
    """Move counts donor->recipient within columns, preserving column totals."""
    left = int(amount)
    if left <= 0:
        return 0

    # Earlier passes preserve more of the observed zero pattern. Recomputing masks
    # after every pass is intentional because prior transfers change capacities.
    priorities = (
        lambda: (x[donor] > 1) & (x[recipient] > 0),
        lambda: (x[donor] > 1) & template_support[recipient],
        lambda: (x[donor] > 0) & (x[recipient] > 0),
        lambda: (x[donor] > 0) & template_support[recipient],
        lambda: x[donor] > 1,
        lambda: x[donor] > 0,  # guaranteed fallback
    )

    for make_mask in priorities:
        if left == 0:
            break
        cols = np.flatnonzero(make_mask())
        if cols.size == 0:
            continue
        keep = 1 if np.all(x[donor, cols] > 1) else 0
        caps = x[donor, cols].astype(np.int64) - keep
        available = int(caps.sum())
        if available <= 0:
            continue
        take = min(left, available)
        moved = _spread_under_caps(take, caps, rng)
        nz = moved > 0
        cols = cols[nz]
        moved = moved[nz].astype(x.dtype, copy=False)
        x[donor, cols] -= moved
        x[recipient, cols] += moved
        left -= take

    if left:
        raise RuntimeError(
            f"could not move {left} of {amount} counts from row {donor} to {recipient}"
        )
    return int(amount)


def _repair_row_margins(
    x: np.ndarray,
    desired_rows: np.ndarray,
    template_support: np.ndarray,
    rng: np.random.Generator,
) -> int:
    """Make row margins exact via within-column moves.

    Donors and recipients are matched in sorted inverse-library-size order. For
    cost ``|1/r_i - 1/r_j|`` this monotone matching minimizes the total absolute
    arithmetic-CPM disturbance before gene/support constraints are considered.
    """
    current = x.sum(axis=1, dtype=np.int64)
    delta = current - desired_rows
    if int(delta.sum()) != 0:
        raise RuntimeError("row and column grand totals disagree")

    donors = np.flatnonzero(delta > 0)
    recipients = np.flatnonzero(delta < 0)
    inv_depth = 1.0 / desired_rows.astype(np.float64)
    donors = donors[np.argsort(inv_depth[donors], kind="stable")]
    recipients = recipients[np.argsort(inv_depth[recipients], kind="stable")]

    di = ri = 0
    moved_total = 0
    while di < donors.size and ri < recipients.size:
        donor = int(donors[di])
        recipient = int(recipients[ri])
        amount = min(int(delta[donor]), int(-delta[recipient]))
        moved_total += _move_with_priority(
            x, template_support, donor, recipient, amount, rng
        )
        delta[donor] -= amount
        delta[recipient] += amount
        if delta[donor] == 0:
            di += 1
        if delta[recipient] == 0:
            ri += 1

    if np.any(delta):
        raise RuntimeError(f"row-margin repair left residuals: {delta}")
    return moved_total


def template_scaled_cpm_target(
    template: sp.csr_matrix,
    column_totals: np.ndarray,
    row_margins: np.ndarray,
    *,
    template_library_sizes: np.ndarray | None = None,
) -> np.ndarray:
    """Gene-wise CPM target preserving control depth/expression correlation."""
    if template_library_sizes is None:
        base_rows = np.asarray(template.sum(axis=1), dtype=np.float64).ravel()
    else:
        base_rows = np.asarray(template_library_sizes, dtype=np.float64)
        if (
            base_rows.shape != (template.shape[0],)
            or np.any(base_rows <= 0)
            or not np.all(np.isfinite(base_rows))
        ):
            raise ValueError(
                "template_library_sizes must be a positive finite vector with "
                "one value per template row"
            )
    base_cols = np.asarray(template.sum(axis=0), dtype=np.float64).ravel()
    row_of_nnz = np.repeat(np.arange(template.shape[0]), np.diff(template.indptr))
    scaled = template.astype(np.float64, copy=True)
    scaled.data *= CPM_TARGET_SUM / base_rows[row_of_nnz]
    base_mean_cpm = np.asarray(scaled.mean(axis=0), dtype=np.float64).ravel()

    out = np.zeros(template.shape[1], dtype=np.float64)
    observed = base_cols > 0
    out[observed] = (
        base_mean_cpm[observed] * column_totals[observed] / base_cols[observed]
    )
    absent = ~observed & (column_totals > 0)
    out[absent] = CPM_TARGET_SUM * column_totals[absent] / float(row_margins.sum())
    return out


def emit_template_fixed_margins(
    template,
    *,
    profile: np.ndarray | None = None,
    column_totals: np.ndarray | None = None,
    row_margins: np.ndarray | None = None,
    template_library_sizes: np.ndarray | None = None,
    seed: int = 0,
) -> tuple[sp.csr_matrix, TemplateEmitterDiagnostics]:
    """Emit an integer block with exact row and column margins.

    Parameters
    ----------
    template
        Selected control cells, raw counts, shape ``(n_cells, n_genes)``.
    profile
        Model relative-expression profile. It is Hamilton-rounded at the grand
        total implied by ``row_margins``.
    column_totals
        Exact intended integer pseudobulk. Use this instead of ``profile`` to
        preserve an existing emitter's pseudobulk bit-for-bit.
    row_margins
        Exact desired library sizes. Defaults to the template's row sums. When
        exact ``column_totals`` have a different grand total and row margins were
        not supplied, the template library sizes are Hamilton-rescaled.
    template_library_sizes
        Optional original full-axis control library sizes. This matters when the
        template has been aligned to a smaller prediction gene axis: use these
        values for row-margin defaults and depth-aware CPM diagnostics rather
        than the smaller aligned row sums.
    seed
        Reproducible tie-breaking only; it cannot change either margin.
    """
    if (profile is None) == (column_totals is None):
        raise ValueError("provide exactly one of profile or column_totals")

    z = _as_integer_csr(template)
    n, genes = z.shape
    rng = np.random.default_rng(seed)
    aligned_template_rows = np.asarray(z.sum(axis=1), dtype=np.int64).ravel()
    if np.any(aligned_template_rows <= 0):
        raise ValueError("template rows must all have positive library size")
    if template_library_sizes is None:
        template_rows = aligned_template_rows
    else:
        template_rows = np.asarray(template_library_sizes)
        if (
            template_rows.shape != (n,)
            or np.any(template_rows <= 0)
            or not np.all(np.isfinite(template_rows))
            or not np.array_equal(template_rows, np.rint(template_rows))
        ):
            raise ValueError(
                "template_library_sizes must be a positive integer-valued "
                f"length-{n} vector"
            )
        template_rows = np.rint(template_rows).astype(np.int64)

    explicit_rows = row_margins is not None
    if explicit_rows:
        rows = np.asarray(row_margins)
        if rows.shape != (n,) or np.any(rows <= 0) or not np.all(np.isfinite(rows)):
            raise ValueError(f"row_margins must be a positive length-{n} vector")
        if not np.array_equal(rows, np.rint(rows)):
            raise ValueError("row_margins must be integer-valued")
        rows = np.rint(rows).astype(np.int64)
    else:
        rows = template_rows.copy()

    if column_totals is None:
        p = np.asarray(profile, dtype=np.float64)
        if p.shape != (genes,):
            raise ValueError(f"profile must have shape ({genes},), got {p.shape}")
        cols = integer_column_totals(p, int(rows.sum()), seed=seed + 1)
    else:
        cols = np.asarray(column_totals)
        if cols.shape != (genes,) or np.any(cols < 0) or not np.all(np.isfinite(cols)):
            raise ValueError(f"column_totals must be a non-negative length-{genes} vector")
        if not np.array_equal(cols, np.rint(cols)):
            raise ValueError("column_totals must be integer-valued")
        cols = np.rint(cols).astype(np.int64)
        if not explicit_rows and int(rows.sum()) != int(cols.sum()):
            rows = _largest_remainder(int(cols.sum()), template_rows, rng)

    if int(rows.sum()) != int(cols.sum()):
        raise ValueError(
            f"row margins sum to {rows.sum()}, column margins sum to {cols.sum()}"
        )
    if np.any(rows > 1_000_000):
        raise ValueError("a requested row margin exceeds the VCC 1,000,000-count cap")
    if int(cols.max(initial=0)) > np.iinfo(np.int32).max:
        raise ValueError("a column total exceeds int32 capacity")

    # Column-wise template projection. This exactly preserves every column
    # margin and preserves the template support wherever that column was seen.
    # Dense int32 is only ~30 MiB at 400 x 18,533 and makes exact repair cheap.
    out = np.zeros((n, genes), dtype=np.int32)
    csc = z.tocsc()
    for gene in np.flatnonzero(cols):
        lo, hi = csc.indptr[gene : gene + 2]
        support_rows = csc.indices[lo:hi]
        weights = csc.data[lo:hi]
        if support_rows.size == 0:
            support_rows = np.arange(n, dtype=np.int64)
            weights = rows
        alloc = _largest_remainder(int(cols[gene]), weights, rng)
        nz = alloc > 0
        out[support_rows[nz], gene] = alloc[nz].astype(np.int32)

    before = out.sum(axis=1, dtype=np.int64) - rows
    template_support = z.astype(bool).toarray()
    moved = _repair_row_margins(out, rows, template_support, rng)

    # Hard invariants: fail here, close to the construction.
    got_rows = out.sum(axis=1, dtype=np.int64)
    got_cols = out.sum(axis=0, dtype=np.int64)
    if not np.array_equal(got_rows, rows):
        raise RuntimeError("internal error: row margins are not exact")
    if not np.array_equal(got_cols, cols):
        raise RuntimeError("internal error: column margins are not exact")
    if np.any(out < 0):
        raise RuntimeError("internal error: negative count")

    output_support = out > 0
    retained = int(np.count_nonzero(output_support & template_support))
    removed = int(np.count_nonzero(template_support & ~output_support))
    introduced = int(np.count_nonzero(output_support & ~template_support))

    result = sp.csr_matrix(out)
    result.eliminate_zeros()
    target_cpm = template_scaled_cpm_target(
        z,
        cols,
        rows,
        template_library_sizes=template_rows,
    )
    actual_cpm = arithmetic_mean_cpm(result)
    denom = float(np.abs(target_cpm).sum())
    cpm_rel_l1 = (
        float(np.abs(actual_cpm - target_cpm).sum() / denom) if denom > 0 else 0.0
    )
    diagnostics = TemplateEmitterDiagnostics(
        grand_total=int(cols.sum()),
        row_l1_before_repair=int(np.abs(before).sum()),
        max_row_error_before_repair=int(np.abs(before).max(initial=0)),
        moved_counts=int(moved),
        template_nnz=int(z.nnz),
        output_nnz=int(result.nnz),
        retained_template_nnz=retained,
        removed_template_nnz=removed,
        introduced_nnz=introduced,
        template_cpm_relative_l1=cpm_rel_l1,
    )
    return result, diagnostics
