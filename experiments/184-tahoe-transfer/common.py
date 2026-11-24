"""184 -- shared loading and the centring convention.

Centring: R(d, L) = delta(d, L) - mean over the drugs measured on the SAME (line, plate) of
delta(., L).  Two things are removed and both must be:

  * the line's mean response across drugs -- the metric's zero (159 sec.5.1, 161 sec.1); a
    predictor credited for it is credited for a hub.
  * the plate-matched DMSO control's own noise, which enters every delta on that (line, plate)
    with the same sign and is exactly the shape of the 86.5% "cell-pool blob" that made 159's
    brief misread the Orion harnesses.

`global` centring (across all of a line's drugs, pooling plates) is kept as a robustness variant.
"""
from __future__ import annotations

import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
SEED = 20260903


def load(centring: str = "plate"):
    z = np.load(os.path.join(WORK, "tensor.npz"), allow_pickle=True)
    D = z["delta"].astype(np.float32)                 # (T, L, G), NaN where unmeasured
    mask = z["group"]                                 # (T, L) bitmask over plates 1..14
    genes = np.asarray([str(g) for g in z["genes"]])
    lines = np.asarray([str(l) for l in z["lines"]])
    treats = np.asarray([str(t) for t in z["treatments"]])
    basal = z["basal_line"].astype(np.float64)        # (L, G) fraction of control library
    R = centre(D, mask, centring)
    return R, D, genes, lines, treats, basal, mask


def centre(D, mask, how="plate"):
    """Return the centred residual, NaN-safe.

    `plate`  subtract, from each measurement, the mean of the (line, plate) blobs of the plates
             it averages -- each blob estimated from that plate's ~94 drugs in that line.  Never
             annihilates a row, unlike centring inside the measurement's own plate-set group.
    `global` subtract the line's mean over all its drugs.
    """
    R = D.astype(np.float32).copy()
    T, L, G = D.shape
    if how == "none":
        return R
    if how == "global":
        for l in range(L):
            R[:, l, :] -= np.nanmean(R[:, l, :], axis=0)
        return R
    PM = np.load(os.path.join(WORK, "platemean.npy"))
    for l in range(L):
        for code in np.unique(mask[:, l]):
            if code <= 0:
                continue
            ps = [p for p in range(1, 15) if code >> p & 1]
            R[mask[:, l] == code, l, :] -= PM[ps, l].mean(0)
    return R


def unit(A, axis=-1):
    n = np.linalg.norm(A, axis=axis, keepdims=True)
    return A / np.maximum(n, 1e-30)


def rowcos(A, B):
    """Per-row cosine between two (n, G) stacks."""
    return (unit(A) * unit(B)).sum(-1)


def sign_acc(pred, truth, top=200):
    """Per-row sign accuracy over each row's own top-|truth| genes, macro-averaged (rule 8).

    Never pools (row, gene) pairs: every row contributes one number, then they are averaged.
    """
    out = np.empty(len(pred))
    for i in range(len(pred)):
        t = truth[i]
        k = min(top, t.shape[0]) if top else t.shape[0]
        idx = np.argpartition(-np.abs(t), k - 1)[:k] if k < t.shape[0] else np.arange(t.shape[0])
        s = np.sign(t[idx]) * np.sign(pred[i][idx])
        ok = s > 0
        tie = s == 0
        out[i] = (ok.sum() + 0.5 * tie.sum()) / len(idx)
    return out
