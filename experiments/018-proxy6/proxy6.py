"""018 — a fast six-metric proxy on the LODO folds.

## Why this exists

Every DE metric so far has needed a count-level `cell-eval2` run: hours per configuration,
and the one fast count-level harness we built (`h1-fast`) demonstrably lied about `jac`.
That has meant `nmae`, `fid`, `reach` and `jac` were effectively untunable, which is a
problem because the board says `nmae` is the second-largest reachable gap we have
(Vivai +0.120 scaled, us -0.020).

The observation that makes a proxy possible: **five of the six metrics are functionals of the
predicted and true effect vectors**, not of the counts, once two things are stipulated:

  - which genes the reference calls DE  -> proxy: the top `K` by |true effect|
  - which genes *we* call DE            -> proxy: the top `K_pred` by |predicted effect|,
    which is what the Wilcoxon ordering reduces to under a minimum-variance emitter

Both stipulations are approximations, and they are why this harness ranks configurations
rather than predicting leaderboard values. `pds` alone is exact.

`mse` is omitted: it is provably clamped to 0 at any cosine we can reach (log, 2026-08-21).
"""

from __future__ import annotations

import numpy as np


def pds_cosine(P: np.ndarray, T: np.ndarray) -> float:
    Pn = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    Tn = T / (np.linalg.norm(T, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - Pn @ Tn.T
    n = D.shape[1]
    sd = np.diag(D)[:, None]
    return float((1.0 - ((D < sd).sum(1) + 0.5 * ((D == sd).sum(1) - 1)) / (n - 1)).mean())


def topk_mask(A: np.ndarray, k: int) -> np.ndarray:
    """Boolean mask of the k largest |A| per row."""
    k = min(k, A.shape[1])
    idx = np.argpartition(-np.abs(A), k - 1, axis=1)[:, :k]
    m = np.zeros(A.shape, bool)
    np.put_along_axis(m, idx, True, axis=1)
    return m


def nmae_curve(P: np.ndarray, T: np.ndarray, R: np.ndarray,
               alphas: np.ndarray) -> np.ndarray:
    """nmae vs the global magnitude alpha, per-perturbation then averaged.

    Piecewise linear and convex in alpha, so a grid resolves the optimum cleanly.
    """
    den = np.maximum(np.where(R, np.abs(T), 0.0).sum(1), 1e-12)
    return np.array([float((np.where(R, np.abs(a * P - T), 0.0).sum(1) / den).mean())
                     for a in alphas])


def de_metrics(P: np.ndarray, T: np.ndarray, R: np.ndarray, k_pred: int,
               conf: np.ndarray | None = None) -> dict:
    """fid / jac / reach given the reference set `R` and our top-`k_pred` calls.

    `conf` orders the reference's DE genes for `reach`; it defaults to |P|, which is what
    the emitted Wilcoxon p-values reduce to under a minimum-variance emitter.
    """
    C = topk_mask(P, k_pred)
    agree = np.sign(P) == np.sign(T)
    k_real = float(R[0].sum())

    n_pred = C.sum(1)
    fid = float(((C & agree).sum(1) / np.maximum(n_pred, 1)
                 * np.minimum(1.0, n_pred / k_real)).mean())
    jac = float(((C & R).sum(1) / np.maximum((C | R).sum(1), 1)).mean())

    if conf is None:
        conf = np.abs(P)
    reach = []
    for i in range(P.shape[0]):
        g = np.flatnonzero(R[i])
        order = g[np.argsort(-conf[i, g], kind="stable")]
        purity = np.cumsum(agree[i, order]) / np.arange(1, len(order) + 1)
        good = np.flatnonzero(purity >= 0.9)
        reach.append((good[-1] + 1) / len(order) if good.size else 0.0)
    return dict(fid=fid, jac=jac, reach=float(np.mean(reach)))


# Official 2026 anchors (AGENTS.md section 4), midpoints of the A/B/C ranges. Used only to
# put proxy raw values on a comparable footing -- never to predict a leaderboard score.
ANCHORS = {
    "pds":   (0.500, 0.9525),
    "nmae":  (1.0013, 0.400),
    "fid":   (0.5135, 0.8135),
    "reach": (0.072, 0.968),
    "jac":   (0.029, 0.399),
}


def scale(name: str, raw: float) -> float:
    b, r = ANCHORS[name]
    return (raw - b) / (r - b)
