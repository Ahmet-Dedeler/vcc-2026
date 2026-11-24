#!/usr/bin/env python3
"""380 — SVD deflation of the FINAL cap15 design's target-shared directions.

Pure numpy construction, no I/O beyond what the caller passes in. See PLAN.md for the exact
recipe and the record check that establishes this is a new object (not 156/175/204/209/246/249).
"""
from __future__ import annotations

import numpy as np

RAND_SEED = 20260910   # 380's own seed, distinct from 218's 20726 and 246's 20260904


def _restore_row_norms(original: np.ndarray, deflated: np.ndarray) -> np.ndarray:
    orig_norm = np.linalg.norm(original, axis=1, keepdims=True)
    new_norm = np.linalg.norm(deflated, axis=1, keepdims=True)
    safe = np.where(new_norm > 0, new_norm, 1.0)
    scale = np.where(new_norm > 0, orig_norm / safe, 1.0)
    return deflated * scale


def svd_energy(X: np.ndarray, centered: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Economy SVD of X (or its across-target-centred form). Returns (Xc, U, S, Vt)."""
    mu = X.mean(axis=0, keepdims=True) if centered else np.zeros((1, X.shape[1]))
    Xc = X - mu
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Xc, U, S, Vt, mu


def deflate_topk(X: np.ndarray, k: int, *, centered: bool) -> dict:
    """Remove the top-k target-shared singular directions of X, norm-preserving per row.

    Returns dict with keys: effects (n,g), energy_removed_frac, centered, k, singular_values
    (first min(20,rank) of them, for the record).
    """
    if k == 0:
        return dict(effects=X.copy(), energy_removed_frac=0.0, centered=centered, k=0,
                    singular_values=None, total_energy=float((X * X).sum()))
    Xc, U, S, Vt, mu = svd_energy(X, centered)
    total_energy = float((S ** 2).sum())
    removed_energy = float((S[:k] ** 2).sum())
    proj = (U[:, :k] * S[:k]) @ Vt[:k]
    Xc_def = Xc - proj
    Xk = Xc_def + mu if centered else Xc_def
    out = _restore_row_norms(X, Xk)
    return dict(effects=out, energy_removed_frac=removed_energy / total_energy if total_energy > 0 else 0.0,
                centered=centered, k=int(k), singular_values=S[:20].tolist(),
                total_energy=total_energy)


def deflate_random(X: np.ndarray, k: int, *, centered: bool, seed: int = RAND_SEED) -> dict:
    """Control: remove a random k-dim subspace of the (centred) ROW space -- 246's `rand` mode
    (QR of a seeded Gaussian n_targets x k matrix), not energy-matched (see PLAN.md for why an
    exact match is not achievable at k=3 by chance, and why that gap is itself informative)."""
    Xc, U, S, Vt, mu = svd_energy(X, centered)
    total_energy = float((S ** 2).sum())
    rng = np.random.default_rng(seed + 1000 * k)
    n = X.shape[0]
    Uk, _ = np.linalg.qr(rng.standard_normal((n, k)))
    proj = Uk @ (Uk.T @ Xc)
    removed_energy = float((proj * proj).sum())
    Xc_def = Xc - proj
    Xk = Xc_def + mu if centered else Xc_def
    out = _restore_row_norms(X, Xk)
    return dict(effects=out, energy_removed_frac=removed_energy / total_energy if total_energy > 0 else 0.0,
                centered=centered, k=int(k), mode="random", seed=int(seed + 1000 * k),
                total_energy=total_energy)


def deflate_shuffled(X: np.ndarray, k: int, *, centered: bool, seed: int = RAND_SEED) -> dict:
    """Control: the REAL top-k right-singular vectors of X, with a single fixed random
    derangement applied to their GENE-AXIS coordinates, then removed from the true
    (unpermuted) X. NOTE on a design mistake caught and fixed before any heavy run: learning
    the basis on a gene-permuted COPY of X and mapping it back via the inverse permutation
    does NOT scramble anything -- a column permutation is an orthogonal change of basis, so
    that construction recovers X's own true top-k subspace exactly (verified empirically:
    the two constructions' `energy_removed_frac` agreed to float noise in a standalone check).
    The correct null keeps the k learned singular-vector PATTERNS (their exact coefficient
    values/shape/energy) but reassigns each coefficient to a random different gene before
    subtracting from the real (correctly-indexed) X -- same shape, wrong identity."""
    Xc, U, S, Vt, mu = svd_energy(X, centered)
    total_energy = float((S ** 2).sum())
    g = X.shape[1]
    rng = np.random.default_rng(seed + 2_000_003 * k)
    perm = rng.permutation(g)
    while np.any(perm == np.arange(g)):
        perm = rng.permutation(g)
    Vk_scrambled = Vt[:k][:, perm]           # same k row-vectors, gene coordinates reassigned
    proj = (Xc @ Vk_scrambled.T) @ Vk_scrambled
    removed_energy = float((proj * proj).sum())
    Xc_def = Xc - proj
    Xk = Xc_def + mu if centered else Xc_def
    out = _restore_row_norms(X, Xk)
    return dict(effects=out, energy_removed_frac=removed_energy / total_energy if total_energy > 0 else 0.0,
                centered=centered, k=int(k), mode="shuffled", seed=int(seed + 2_000_003 * k),
                permutation_is_derangement=bool(np.all(perm != np.arange(g))),
                real_topk_energy_frac=float((S[:k] ** 2).sum() / total_energy) if total_energy > 0 else 0.0,
                total_energy=total_energy)
