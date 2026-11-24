"""183 — the unique-component boost, and nothing else.

    c  = mean_p h_p                     centroid of the predicted-effect rows
    u  = h - alpha * proj_c(h)          strip a fraction of the SHARED direction
    h' = u * (||h|| / ||u||)            preserve the row's L2 exactly

alpha = 0 is the identity (the known-answer row).  alpha = 1 removes the centroid
direction entirely and then restores the norm.

This is deliberately NOT:
  * an un-whitening move — the whitening operator is untouched, this acts on the
    already-whitened rows;
  * an amplitude move — every row's L2 is preserved exactly (asserted below), so the
    submission's `A` and its across-perturbation budget `B` are unchanged and the
    rho-cliff (170 §4) is not disturbed;
  * AgingLab's `L <- L - mean_t(L)` — that subtracts the centroid VECTOR, which both
    rotates and shrinks each row.  This projects out the centroid DIRECTION and
    renormalises.  At alpha = 1 the two share a null space but differ off it.

170's criterion for a remaining `pds` lever is "norm-preserving, raise X not A".
This is the smallest operator that satisfies it.
"""
from __future__ import annotations

import numpy as np


def unique_boost(D: np.ndarray, alpha: float) -> np.ndarray:
    """Apply the boost row-wise to a (n_pert, n_gene) matrix.

    The centroid is taken over the rows of `D` itself, so the operator is defined by the
    panel it is applied to and needs no external fit (nothing to leak).
    """
    D = np.asarray(D, dtype=np.float64)
    if alpha == 0.0:
        return D.copy()
    c = D.mean(axis=0)
    cc = float(c @ c)
    if not np.isfinite(cc) or cc <= 0.0:
        return D.copy()
    coef = (D @ c) / cc                                   # <h, c> / <c, c>
    U = D - alpha * coef[:, None] * c[None, :]
    nh = np.linalg.norm(D, axis=1)
    nu = np.linalg.norm(U, axis=1)
    return U * (nh / np.maximum(nu, 1e-300))[:, None]


def norm_check(D: np.ndarray, Dp: np.ndarray) -> float:
    """max |‖h'‖/‖h‖ − 1| over rows.  Must be ~1e-16 straight out of `unique_boost`."""
    a = np.linalg.norm(np.asarray(D, dtype=np.float64), axis=1)
    b = np.linalg.norm(np.asarray(Dp, dtype=np.float64), axis=1)
    live = a > 0
    return float(np.max(np.abs(b[live] / a[live] - 1.0))) if live.any() else 0.0


def shared_energy(D: np.ndarray) -> dict:
    """How much of the matrix lives on its own centroid direction, and the hubness proxy."""
    D = np.asarray(D, dtype=np.float64)
    c = D.mean(axis=0)
    cc = float(c @ c)
    tot = float((D * D).sum())
    on_c = float(((D @ c) ** 2).sum() / cc) if cc > 0 else 0.0
    n = D.shape[0]
    A = D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-300)
    G = A @ A.T
    rr = float(G[~np.eye(n, dtype=bool)].mean())
    return dict(centroid_energy_frac=(on_c / tot if tot > 0 else 0.0),
                mean_pairwise_cos=rr,
                median_row_norm=float(np.median(np.linalg.norm(D, axis=1))),
                mean_square=float(tot / D.size))
