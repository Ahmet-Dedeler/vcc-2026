"""Source-population scale estimation, with no recipient truth access."""

from contextlib import contextmanager

import numpy as np


def gene_scale(x):
    x = np.asarray(x, float)
    if x.ndim != 2 or x.shape[0] < 2 or not x.shape[1] or not np.isfinite(x).all():
        raise ValueError("expected finite target by gene matrix")
    center = x.mean(0)
    mean_abs = np.abs(x-center).mean(0)
    positive = mean_abs[mean_abs > 0]
    if positive.size != mean_abs.size:
        raise ValueError("degenerate background response column")
    floor = float(np.percentile(positive, 25))
    return np.maximum(mean_abs, floor), {"n": len(x), "floor": floor}


def mixed_scale(original, columns, replacement, eligible_columns):
    """Replace explicitly indexed source columns only; preserve the rest exactly."""
    original = np.asarray(original, float)
    columns, eligible_columns = np.asarray(columns), np.asarray(eligible_columns)
    replacement = np.asarray(replacement, float)
    if (original.shape != columns.shape or eligible_columns.shape != replacement.shape or
            len(set(columns.tolist())) != len(columns) or
            len(set(eligible_columns.tolist())) != len(eligible_columns) or
            not np.isfinite(replacement).all() or np.any(replacement <= 0)):
        raise ValueError("invalid replacement scale or axes")
    lookup = dict(zip(eligible_columns.tolist(), replacement.tolist()))
    out = original.copy()
    for i, column in enumerate(columns):
        if column in lookup:
            out[i] = lookup[column]
    if not np.isfinite(out).all() or np.any(out <= 0):
        raise ValueError("nonpositive combined scale")
    return out


def median_normalized_norm(delta, columns, mean, scale, chunk_rows=128):
    values = np.empty(len(delta), float)
    for lo in range(0, len(delta), chunk_rows):
        block = np.asarray(delta[lo:lo+chunk_rows, columns], float)-mean
        values[lo:lo+len(block)] = np.linalg.norm(block/scale, axis=1)
    return float(np.median(values)) or 1.


def permute_scales(scales, seed=20260909):
    scales = np.asarray(scales, float)
    if scales.ndim != 2:
        raise ValueError("condition by gene matrix required")
    rng = np.random.default_rng(seed)
    return np.array([rng.permutation(row) for row in scales])


@contextmanager
def scale_override(module, replacements, eligible_columns, audits):
    """Patch the imported function only during preparation, then restore it."""
    original = module.centered_statistics

    def wrapped(bank, columns, mean, *, need_gene_scale, chunk_rows, max_rss_gib):
        raw, scale, median = original(
            bank, columns, mean, need_gene_scale=need_gene_scale,
            chunk_rows=chunk_rows, max_rss_gib=max_rss_gib)
        if replacements is None:
            return raw, scale, median
        if not need_gene_scale or bank.name not in replacements:
            raise ValueError("unexpected source or normalization mode")
        changed = mixed_scale(scale, columns, replacements[bank.name], eligible_columns)
        delta = np.load(bank.delta_path, mmap_mode="r", allow_pickle=False)
        new_median = median_normalized_norm(delta, columns, mean, changed, chunk_rows)
        ratios = changed/scale
        audits[bank.name] = {
            "changed_coordinates": int(np.count_nonzero(changed != scale)),
            "scale_ratio_min_median_max": list(map(float, [ratios.min(), np.median(ratios), ratios.max()])),
            "original_normalized_median": median, "new_normalized_median": new_median}
        return raw, changed, new_median

    module.centered_statistics = wrapped
    try:
        yield
    finally:
        module.centered_statistics = original
