"""Small pure functions for the frozen coordinate-median construction."""

import numpy as np


def eligible_mask(observed, operator_mask):
    observed = np.asarray(observed, dtype=bool)
    if observed.ndim != 3 or observed.shape[0] != 3:
        raise ValueError("Exactly three donor views are required")
    return observed.all(axis=0) & np.asarray(operator_mask, dtype=bool)[None, :]


def restore_subvector(base, candidate, eligible):
    """Restore only eligible norms, leaving all other entries bit-identical."""
    base, candidate = np.asarray(base, float), np.asarray(candidate, float)
    eligible = np.asarray(eligible, bool)
    if base.shape != candidate.shape or base.shape != eligible.shape:
        raise ValueError("axis/shape mismatch")
    if not np.isfinite(base).all() or not np.isfinite(candidate).all():
        raise ValueError("nonfinite effect")
    result = base.copy()
    for i in range(len(base)):
        mask = eligible[i]
        nb, nc = np.linalg.norm(base[i, mask]), np.linalg.norm(candidate[i, mask])
        if nb > 0 and nc > 0:
            result[i, mask] = candidate[i, mask] * (nb / nc)
    return result


def interpolate(base, alternate, eligible, dose):
    if not 0 <= dose <= 1:
        raise ValueError("dose outside [0,1]")
    if dose == 0:
        return np.asarray(base).copy()
    candidate = np.asarray(base, float) + dose * (
        np.asarray(alternate, float) - np.asarray(base, float))
    return restore_subvector(base, candidate, eligible)


def arms(base, sources, observed, operator_mask, dose=.25, seed=20260909):
    sources = np.asarray(sources, float)
    if sources.shape != np.asarray(observed).shape:
        raise ValueError("source and observation shapes differ")
    if sources.shape[1:] != np.asarray(base).shape:
        raise ValueError("source and baseline shapes differ")
    if not np.isfinite(sources).all():
        raise ValueError("nonfinite source")
    eligible = eligible_mask(observed, operator_mask)
    median = np.median(sources, axis=0)
    robust = interpolate(base, median, eligible, dose)
    magnitude = robust.copy()
    nonzero = np.asarray(base) != 0
    magnitude[nonzero] = np.sign(base[nonzero]) * np.abs(robust[nonzero])
    power = restore_subvector(base, np.sign(base) * np.square(base), eligible)
    power = interpolate(base, power, eligible, dose)
    rng = np.random.default_rng(seed)
    shuffled = np.asarray(base, float).copy()
    correction = robust - base
    for i in range(len(base)):
        mask = eligible[i]
        shuffled[i, mask] += rng.permutation(correction[i, mask])
    shuffled = restore_subvector(base, shuffled, eligible)
    return {
        "cap15": np.asarray(base).copy(), "median25": robust,
        "magnitude25": magnitude, "power25": power, "shuffled25": shuffled,
    }, eligible


def geometry(candidate, baseline, mask=None):
    a, b = np.asarray(candidate, float), np.asarray(baseline, float)
    if mask is not None:
        a, b = a[:, mask], b[:, mask]
    nb = np.linalg.norm(b, axis=1)
    na = np.linalg.norm(a, axis=1)
    ratios = np.divide(na, nb, out=np.ones_like(na), where=nb > 0)
    return {
        "relative_change_frobenius": float(np.linalg.norm(a-b) / max(np.linalg.norm(b), 1e-300)),
        "row_norm_ratio_min_median_max": list(map(float, [ratios.min(), np.median(ratios), ratios.max()])),
        "sign_changes_nonzero_baseline": int(((np.sign(a) != np.sign(b)) & (b != 0)).sum()),
        "nonzero_baseline_entries": int(np.count_nonzero(b)),
    }
