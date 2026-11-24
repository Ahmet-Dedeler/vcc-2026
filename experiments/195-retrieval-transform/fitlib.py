"""195 — fitting a retrieval operator to the ranking objective, in reduced coordinates.

## The family

The production operator is one member of a hand-chosen one-parameter family,
`W = (C + lam I)^(-gamma/2)` at `gamma = 0.65`.  We keep it as the initialisation and hang a
learnable correction off it, exactly so that the zero-parameter row IS the incumbent:

    W(M) = W_prod (I + E_k M E_k^T),      M = U V^T + diag(a),   U, V in R^{k x r}

`E_k` is an orthonormal k-dimensional subspace built from donor data only.  With `M = 0` the
operator is bit-identical to production's, so the learning curve starts at the incumbent's
score by construction rather than by luck.

## Why this is cheap

`pds` and the listwise loss depend on `W(M)` only through cosines between whitened queries and
raw gallery vectors, and `E_k` is orthonormal, so the whole 7,760-dimensional problem collapses
to k dimensions plus two per-batch constants.  Writing `Z = P W_prod`, `Zp = Z E_k`,
`Tb = T E_k`, and `G = Zp (I + M)`:

    <Z'_i, t_j>  =  (<Z_i, t_j> - <Zp_i, Tb_j>)  +  <G_i, Tb_j>
    ||Z'_i||^2   =  (||Z_i||^2 - ||Zp_i||^2)     +  ||G_i||^2

Both parenthesised terms are independent of `M`, so they are precomputed once and every
gradient step costs O(B^2 k) instead of O(B G).  That is what makes a rank sweep, a learning
curve, and a permutation null affordable on one machine.

## Objective

In-batch listwise softmax cross-entropy.  The batch is the gallery, so a batch of 256 makes
each training problem a 256-way retrieval -- deliberately matched to the 259-way panel problem
the metric actually poses.  The reported curve is in `pds` units (mean fraction of the gallery
ranked worse than the true answer), not in loss units, because 4.2's whole lesson is that the
proxy and the rank can move in opposite directions.
"""

from __future__ import annotations

import numpy as np


def make_batches(n: int, batch: int, rng: np.random.Generator) -> list[np.ndarray]:
    """Random equal-size partition; the tail short batch is dropped so gallery size is fixed."""
    perm = rng.permutation(n)
    nb = n // batch
    return [perm[i * batch:(i + 1) * batch] for i in range(nb)]


class BatchStats:
    """The `M`-independent constants for one batch (see module docstring)."""

    __slots__ = ("R", "Zp", "Tb", "c", "tn", "B")

    def __init__(self, Z: np.ndarray, T: np.ndarray, Ek: np.ndarray):
        Z = np.asarray(Z, np.float64)
        T = np.asarray(T, np.float64)
        self.Zp = Z @ Ek
        self.Tb = T @ Ek
        self.R = Z @ T.T - self.Zp @ self.Tb.T
        self.c = np.einsum("ij,ij->i", Z, Z) - np.einsum("ij,ij->i", self.Zp, self.Zp)
        self.tn = np.linalg.norm(T, axis=1)
        self.B = Z.shape[0]


def batch_forward(bs: BatchStats, M: np.ndarray, tau: float,
                  labels: np.ndarray | None = None):
    """Loss, in-batch pds, and dL/dM for one batch."""
    k = M.shape[0]
    G = bs.Zp + bs.Zp @ M                       # (B, k)
    num = bs.R + G @ bs.Tb.T                    # (B, B)
    q = np.sqrt(bs.c + np.einsum("ij,ij->i", G, G))
    cos = num / (q[:, None] * bs.tn[None, :])
    lab = np.arange(bs.B) if labels is None else labels

    s = cos / tau
    s = s - s.max(1, keepdims=True)
    e = np.exp(s)
    p = e / e.sum(1, keepdims=True)
    loss = float(-np.log(np.maximum(p[np.arange(bs.B), lab], 1e-300)).mean())

    sd = cos[np.arange(bs.B), lab][:, None]
    worse = (cos < sd).sum(1) + 0.5 * ((cos == sd).sum(1) - 1)
    pds = float((worse / (bs.B - 1)).mean())

    gp = p.copy()
    gp[np.arange(bs.B), lab] -= 1.0
    g = gp / (bs.B * tau)                        # dL/dcos
    dG = (g / bs.tn[None, :]) @ bs.Tb / q[:, None]
    dG -= (G / (q ** 2)[:, None]) * np.einsum("ij,ij->i", g, cos)[:, None]
    dM = bs.Zp.T @ dG
    assert dM.shape == (k, k)
    return loss, pds, dM


class Adapter:
    """`M = U V^T + diag(a)`, Adam, with a Frobenius penalty on the product (as in 080)."""

    def __init__(self, k: int, rank: int, diag: bool, seed: int = 0, init: float = 1e-3):
        rng = np.random.default_rng(seed)
        self.rank = rank
        self.U = rng.normal(0, init, (k, rank)) if rank else np.zeros((k, 0))
        self.V = np.zeros((k, rank)) if rank else np.zeros((k, 0))
        self.a = np.zeros(k) if diag else None
        self.k = k
        self._m = {}
        self._v = {}
        self._t = 0

    def M(self) -> np.ndarray:
        M = self.U @ self.V.T if self.rank else np.zeros((self.k, self.k))
        if self.a is not None:
            M = M + np.diag(self.a)
        return M

    def step(self, dM: np.ndarray, lr: float, reg: float) -> None:
        grads = {}
        Mp = self.U @ self.V.T if self.rank else None
        if self.rank:
            grads["U"] = dM @ self.V + reg * (Mp @ self.V)
            grads["V"] = dM.T @ self.U + reg * (Mp.T @ self.U)
        if self.a is not None:
            grads["a"] = np.diag(dM) + reg * self.a
        self._t += 1
        b1, b2, eps = 0.9, 0.999, 1e-8
        for name, gr in grads.items():
            self._m[name] = b1 * self._m.get(name, 0.0) + (1 - b1) * gr
            self._v[name] = b2 * self._v.get(name, 0.0) + (1 - b2) * gr ** 2
            mh = self._m[name] / (1 - b1 ** self._t)
            vh = self._v[name] / (1 - b2 ** self._t)
            setattr(self, name, getattr(self, name) - lr * mh / (np.sqrt(vh) + eps))

    def n_params(self) -> int:
        return (2 * self.k * self.rank if self.rank else 0) + (self.k if self.a is not None else 0)


def apply_adapter(Z: np.ndarray, Ek: np.ndarray, M: np.ndarray) -> np.ndarray:
    """`Z (I + E_k M E_k^T)` -- the whitened query after the learned correction."""
    if M is None:
        return Z
    return Z + ((Z @ Ek) @ M) @ Ek.T
