"""Frozen posterior inference and audit utilities, independent of H1 truth."""

import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "328-joint-cd4-estimation"))
import mixture328 as M


def extractor_delta(logfc):
    values = np.asarray(logfc, dtype=np.float32)
    # Match the extractor's Python-float scalar, not np.float64. NumPy 2's
    # scalar-promotion rules otherwise change rounding before float32 storage.
    return np.where(np.isfinite(values), values * math.log(2.), 0.).astype(np.float32)


def validate_fit(fit, *, diagonal=False):
    if not fit.get("converged", False):
        raise ValueError("nonconverged source fit")
    if fit["components"] != M.components(diagonal=diagonal):
        raise ValueError("source component grid changed")
    weights = np.asarray(fit["weights"], float)
    if (weights.shape != (len(fit["components"]),) or
            not np.isfinite(weights).all() or np.any(weights < 0) or
            abs(weights.sum()-1) > 1e-10):
        raise ValueError("invalid source mixture weights")


def replace_deltas(original, logfc, se, present, gene_eligible, fit, *,
                   diagonal=False, batch_size=8192):
    """All arrays use target/gene/condition order; absence never supplies evidence."""
    original = np.asarray(original)
    logfc, se = np.asarray(logfc), np.asarray(se)
    present, gene_eligible = np.asarray(present, bool), np.asarray(gene_eligible, bool)
    if (original.ndim != 3 or original.shape[-1] != 3 or original.dtype != np.float32 or
            logfc.shape != original.shape or se.shape != original.shape or
            present.shape != (original.shape[0], 3) or
            gene_eligible.shape != (original.shape[1],)):
        raise ValueError("target/gene/condition shape or dtype mismatch")
    if not np.isfinite(original).all() or batch_size < 1:
        raise ValueError("invalid original deltas or batch size")
    validate_fit(fit, diagonal=diagonal)
    eligible = (present.all(1)[:, None] & gene_eligible[None, :] &
                np.all(np.isfinite(logfc) & np.isfinite(se) & (se > 0), axis=2))
    out = original.copy()
    flat = np.flatnonzero(eligible.ravel())
    dst = out.reshape(-1, 3)
    y, s = logfc.reshape(-1, 3), se.reshape(-1, 3)
    for lo in range(0, len(flat), batch_size):
        idx = flat[lo:lo+batch_size]
        posterior = M.posterior(y[idx], s[idx], fit["components"], fit["weights"])
        if not np.isfinite(posterior).all():
            raise ArithmeticError("nonfinite posterior")
        dst[idx] = extractor_delta(posterior)
    if not np.array_equal(out[~eligible], original[~eligible]):
        raise AssertionError("raw fallback changed")
    return out, eligible


def geometry(candidate, baseline, mask=None):
    a, b = np.asarray(candidate, float), np.asarray(baseline, float)
    if a.shape != b.shape or a.ndim != 2 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("invalid geometry arrays")
    if mask is not None:
        a, b = a[:, mask], b[:, mask]
    na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
    ratios = np.divide(na, nb, out=np.ones_like(na), where=nb > 0)
    return {
        "relative_change_frobenius": float(np.linalg.norm(a-b) / max(np.linalg.norm(b), 1e-300)),
        "row_norm_ratio_min_median_max": list(map(float, [ratios.min(), np.median(ratios), ratios.max()])),
        "sign_changes_nonzero_baseline": int(((np.sign(a) != np.sign(b)) & (b != 0)).sum()),
        "nonzero_baseline_entries": int(np.count_nonzero(b)),
        "zero_baseline_rows": int(np.sum(nb == 0)),
        "nonzero_candidate_at_zero_baseline_rows": int(np.sum((nb == 0) & (na > 0))),
    }
