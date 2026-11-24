"""Spectral whitening of predicted perturbation effects (experiment 017).

## The mechanism

Perturbation effects are wildly anisotropic. A handful of shared transcriptional programmes
-- growth arrest, ribosome biogenesis, stress response -- account for most of the variance,
so two unrelated knockdowns that both slow the cell down look similar under a cosine. `pds`
asks a retrieval question ("is this prediction nearer its own true effect than to any of the
other 299?"), and retrieval in an anisotropic space suffers from *hubness*: whichever true
effects sit near the dominant programme are near everything, and discrimination collapses.

The fix is the standard one from retrieval and from word-embedding post-processing: rescale
the space so the dominant directions stop dominating.

    W = (C + lam*I)^(-gamma/2),     C = B^T B / n

`B` is a matrix of donor effect vectors, so `C` is the covariance of perturbation responses.
gamma = 0 leaves the prediction alone, gamma = 1 is standard whitening, gamma = 2 is a full
Mahalanobis inverse. Predictions are mapped `P -> P W`.

## Why this is not the low-rank truncation that failed

Experiment 004 truncated the donor SVD and lost *monotonically* -- rank 5 gave cosine 0.039,
rank 800 gave 0.052, full rank 0.073. The tail directions carry signal, and truncation
deletes them. Whitening does the opposite: it **amplifies** exactly those directions, and it
is invertible, so nothing is discarded. The failure of truncation is the reason to expect
whitening to work, not a reason to avoid it.

## Where the basis comes from matters

Estimating a 7,543-dimensional covariance from the 300 panel predictions alone is badly
conditioned, and the resulting operator has no interior optimum -- it just keeps helping as
gamma grows, which is the signature of an artefact rather than a mechanism. Built instead
from every perturbation the donor lines measured (~12,000), the covariance is full rank and
the sweep has a clean interior peak. Use the donor basis.

The operator never sees the target line. It is a function of donor screens only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


class Whitener:
    """`P @ W` with `W = (C + lam I)^(-gamma/2)`, C from a donor effect basis."""

    def __init__(self, evals: np.ndarray, evecs: np.ndarray):
        self.evals = evals            # (G,) ascending or descending, non-negative
        self.evecs = evecs            # (G, G) orthonormal columns

    @classmethod
    def fit(cls, basis: np.ndarray, floor_rel: float = 1e-12) -> "Whitener":
        """Eigendecompose the effect covariance C = basis^T basis / n.

        `eigh` on the G x G Gram is far cheaper than an SVD of the tall basis when the basis
        has more rows than columns, which is the regime we are in (12k perturbations,
        7.5k genes).
        """
        B = np.asarray(basis, dtype=np.float64)
        C = (B.T @ B) / B.shape[0]
        ev, V = np.linalg.eigh(C)
        ev = np.maximum(ev, ev.max() * floor_rel)
        return cls(ev, V)

    def apply(self, P: np.ndarray, gamma: float, lam: float = 0.0,
              renorm: bool = True) -> np.ndarray:
        """Whiten `P` (n x G), optionally restoring its original median row norm.

        `renorm` exists so that the downstream `alpha` keeps its meaning: whitening changes
        the overall scale by orders of magnitude, and alpha is calibrated against a fixed
        reference effect magnitude.
        """
        if gamma == 0.0:
            return P
        P = np.asarray(P, dtype=np.float64)
        w = (self.evals + lam) ** (-gamma / 2.0)
        out = ((P @ self.evecs) * w) @ self.evecs.T
        if renorm:
            before = float(np.median(np.linalg.norm(P, axis=1)))
            after = float(np.median(np.linalg.norm(out, axis=1))) or 1.0
            out *= before / after
        return out

    def save(self, path: str | Path) -> None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        np.savez(p, evals=self.evals.astype(np.float64),
                 evecs=self.evecs.astype(np.float32))

    @classmethod
    def load(cls, path: str | Path) -> "Whitener":
        z = np.load(path)
        return cls(z["evals"], z["evecs"].astype(np.float64))
