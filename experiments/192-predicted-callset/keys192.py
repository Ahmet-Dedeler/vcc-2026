"""192 — shared definitions: the build-time feature block and the candidate call-set keys.

EVERYTHING in `buildable_keys` is computable before the reference is seen:

  * `H`          our predicted effect in the scorer's profile space `log1p(TS*p)`  (we ship it)
  * `cpm`        the control's pooled CPM, measured exactly from the 18,400 NTC cells
  * `ref_lin`    the control's MEAN PER-CELL CPM -- the quantity the DE test's lfc divides by
  * `nz`         the control's per-gene detection rate (fraction of control cells non-zero)
  * `cv2`        the control's per-gene squared coefficient of variation of per-cell CPM
  * `auc`        the EXACT expected Mann-Whitney AUC of the control distribution shifted
                 multiplicatively by our own predicted fold change -- a power proxy that,
                 unlike `|lfc|*sqrt(cpm)`, knows about zero-inflation and discreteness

`oracle_keys` use the reference's own `lfc` and exist only to bound the family (188's 0.7219).
"""
from __future__ import annotations

import os

import numpy as np

TS = 50_000.0
EPS = 1e-9


def h_to_lfc(H, cpm):
    """Invert h = log1p(TS*1e-6*cpm*2^lfc) - log1p(TS*1e-6*cpm).  Exact."""
    base = TS * 1e-6 * np.asarray(cpm, dtype=np.float64)
    num = np.expm1(np.asarray(H, dtype=np.float64) + np.log1p(base)[None, :])
    return np.log2(np.maximum(num, 1e-300) / np.maximum(base, 1e-300)[None, :])


def lfc_to_h(L, cpm):
    base = TS * 1e-6 * np.asarray(cpm, dtype=np.float64)
    return np.log1p(base[None, :] * np.exp2(np.nan_to_num(L))) - np.log1p(base)[None, :]


# ----------------------------------------------------------------------------- control block
def gene_stats(cs_path, cache):
    """Per-gene control statistics from the sorted log1p(CPM) control block (genes x cells).

    Returns nz (detection rate), cv2, mean per-cell CPM, and a 257-point quantile grid of the
    control's log1p(CPM) per gene (for the exact shift-AUC below).
    """
    if os.path.exists(cache):
        z = np.load(cache)
        return {k: z[k] for k in z.files}
    cs = np.load(cs_path, mmap_mode="r")
    G, m = cs.shape
    nz = np.empty(G)
    cv2 = np.empty(G)
    mu = np.empty(G)
    qi = np.linspace(0, m - 1, 129).round().astype(np.int64)
    grid = np.empty((G, qi.size))
    for a in range(0, G, 512):
        b = min(a + 512, G)
        blk = np.asarray(cs[a:b], dtype=np.float64)
        lin = np.expm1(blk)
        nz[a:b] = (blk > 0).mean(1)
        mu[a:b] = lin.mean(1)
        cv2[a:b] = lin.var(1) / np.maximum(lin.mean(1) ** 2, 1e-30)
        grid[a:b] = blk[:, qi]
    out = dict(nz=nz, cv2=cv2, mu=mu, grid=grid)
    np.savez_compressed(cache, **out)
    return out


def shift_auc_table(cs_path, grid, shifts, cache):
    """A_g(s) = P(2^s * X > Y) + 0.5 P(2^s * X = Y), X,Y iid from gene g's control CPM
    distribution; `shifts` are LOG2 fold changes.

    X is the 257-point quantile grid; Y is the full 18,400-cell sorted block.  Computed once
    per (gene, shift) and interpolated in log2(c) afterwards.  This is the Mann-Whitney AUC a
    perturbation would produce if it multiplied every control cell's expression by `c`; it is
    the exact power object the reference's Wilcoxon reads, and it needs no prediction.
    """
    if os.path.exists(cache):
        return np.load(cache)["A"]
    cs = np.load(cs_path, mmap_mode="r")
    G, m = cs.shape
    A = np.empty((G, shifts.size), np.float32)
    for a in range(0, G, 256):
        b = min(a + 256, G)
        Y = np.asarray(cs[a:b], dtype=np.float64)                    # sorted, (nb, m)
        X = np.expm1(grid[a:b])                                      # (nb, 257) linear CPM
        for j, sft in enumerate(shifts):
            Xs = np.log1p(X * (2.0 ** sft))                          # (nb, 129)
            lo = np.empty(Xs.shape)
            hi = np.empty(Xs.shape)
            for r in range(b - a):
                lo[r] = np.searchsorted(Y[r], Xs[r], side="left")
                hi[r] = np.searchsorted(Y[r], Xs[r], side="right")
            A[a:b, j] = ((lo + hi) * 0.5 / m).mean(1)
    np.savez_compressed(cache, A=A, shifts=shifts)
    return A


def auc_key(lfc, A, shifts):
    """|AUC - 0.5| for each (perturbation, gene), interpolating the per-gene AUC table."""
    x = np.clip(np.nan_to_num(lfc), shifts[0], shifts[-1])
    out = np.empty(x.shape, np.float64)
    for g in range(x.shape[1]):
        out[:, g] = np.interp(x[:, g], shifts, A[g])
    return np.abs(out - 0.5)
