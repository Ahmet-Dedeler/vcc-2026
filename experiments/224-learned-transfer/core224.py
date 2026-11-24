"""224 — shared primitives: surfaces, the LOLO pair builder, and every rung's estimator.

`ladder224.py` (Phase 1) and `falsify224.py` (Phase 2) both import from here so that a
falsifier is always the SAME code path as the arm it falsifies, with one input changed.

THE MODEL-SELECTION RULE, and why it is the whole ballgame
---------------------------------------------------------
A first version selected (rank, lambda) for the learned operator by holding out one inner
TRAINING LINE while keeping the same targets.  That is contaminated: the held fold's input
row for target `t` is the mean of `k-2` of the same `k-1` lines that produced the fit rows
for the same `t`, so a near-interpolating ridge (N < G) scores highly on it by memorising
the target.  It selected full rank at lambda 0.01 and then read **-0.010 on the held line
against a +0.188 donor mean.**  Selection here therefore holds out a **disjoint set of
TARGETS as well as a line**, which is the regime the model is deployed in.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

RANKS = (16, 32, 64, 128, 256, 512, 1024, 0)      # 0 = full
LAMS = (1e-3, 1e-2, 1e-1, 1.0, 3.0, 10.0, 30.0, 100.0)
DIAG_SHRINK = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)


def unit(A: np.ndarray) -> np.ndarray:
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def macro_cos(P: np.ndarray, T: np.ndarray) -> float:
    return float((unit(P) * unit(T)).sum(1).mean())


def per_target_cos(P: np.ndarray, T: np.ndarray) -> np.ndarray:
    return (unit(P) * unit(T)).sum(1)


def decoy_cos(P: np.ndarray, T: np.ndarray) -> float:
    """Mean cosine of each prediction against every OTHER target's truth.

    159 §3/§4's attenuation-invariant instrument.  A gain that raises this as much as it
    raises the same-target cosine is amplitude, not transfer — that is exactly the shape
    that killed 184 (arm +0.0665 / control +0.0667 / wrong-recipient null +0.0663)."""
    S = unit(P) @ unit(T).T
    n = S.shape[0]
    return float((S.sum() - np.trace(S)) / (n * (n - 1)))


def pds(P: np.ndarray, T: np.ndarray) -> float:
    A, B = unit(P.astype(np.float64)), unit(T.astype(np.float64))
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S)
    beaten = (S < own[:, None]).sum(1)
    ties = (S == own[:, None]).sum(1) - 1
    return float(np.mean((beaten + 0.5 * ties) / (n - 1)))


def load_surface(name: str):
    """One TABLE per LINE (first listed wins); k562_essential is dropped as a second K562
    library.  Returns C (L,T,G) float32 mean-removed, basal (L,G), names, groups, ..."""
    z = np.load(WORK / f"surface-{name}.npz", allow_pickle=True)
    members = [str(x) for x in z["members"]]
    groups = [str(x) for x in z["groups"]]
    keep, seen = [], set()
    for i, g in enumerate(groups):
        if g not in seen:
            seen.add(g)
            keep.append(i)
    return (np.asarray(z["C"], dtype=np.float32)[keep],
            np.asarray(z["basal"], dtype=np.float32)[keep],
            [members[i] for i in keep], [groups[i] for i in keep],
            [str(t) for t in z["targets"]], [str(g) for g in z["genes"]],
            [members[i] for i in range(len(members)) if i not in keep])


def donor_mean(C: np.ndarray, lines, rows) -> np.ndarray:
    """Accumulated in place: building the list-of-arrays and then np.mean over the stack
    held 3 copies at once and took `crosspanel224.py` to 9.4 GiB."""
    acc = None
    for x in lines:
        u = unit(C[x][rows].astype(np.float64))
        acc = u if acc is None else np.add(acc, u, out=acc)
    return unit(acc / len(lines))


def inner_pairs(C, D, rows, perm=None):
    """Leave-one-training-LINE-out pairs inside the donor pool.
    `perm` deranges the TARGET correspondence of the labels (Phase-2 falsifier 2)."""
    P, Y, F = [], [], []
    for j, B in enumerate(D):
        rest = [x for x in D if x != B]
        P.append(donor_mean(C, rest, rows))
        lab = unit(C[B][rows].astype(np.float64))
        Y.append(lab[perm] if perm is not None else lab)
        F.append(np.full(len(rows), j))
    return np.concatenate(P), np.concatenate(Y), np.concatenate(F)


# ------------------------------------------------------------------- rung 3: the operator
def ridge_svd(Xp: np.ndarray, Y: np.ndarray):
    """SVD once; every (rank, lambda) is a rescaling.  M (G x G) is never materialised."""
    U, s, Vt = np.linalg.svd(Xp, full_matrices=False)
    UtY = U.T @ Y

    def apply(Pnew: np.ndarray, rank: int, lam: float) -> np.ndarray:
        r = len(s) if rank in (0, None) else min(rank, len(s))
        return (Pnew @ Vt[:r].T * (s[:r] / (s[:r] ** 2 + lam))) @ UtY[:r]
    return apply, s


def select_operator(C, D, tsel, tval, hold, perm_sel=None, perm_val=None):
    """(rank, lambda) chosen on a held-out LINE **and** a disjoint target block."""
    Ps, Ys, Fs = inner_pairs(C, D, tsel, perm_sel)
    Pv, Yv, Fv = inner_pairs(C, D, tval, perm_val)
    m, mv = Fs != hold, Fv == hold
    apply, _ = ridge_svd(Ps[m], Ys[m])
    best, bs = (128, 10.0), -np.inf
    grid = []
    for r in RANKS:
        for lam in LAMS:
            v = macro_cos(apply(Pv[mv], r, lam), Yv[mv])
            grid.append((r, lam, v))
            if v > bs:
                bs, best = v, (r, lam)
    return best, float(bs), grid


# ------------------------------------------------------ rung 4: per-gene diagonal operators
def diag_fit(Xp: np.ndarray, Y: np.ndarray, bins: np.ndarray, n_bins: int, shrink: float):
    """g_b = <P,Y>_b / <P,P>_b shrunk toward 1.  `bins` = each gene's group index.
    n_bins == n_genes gives the recipient-INDEPENDENT per-gene diagonal (rung 4's null)."""
    num = np.bincount(bins, weights=(Xp * Y).sum(0), minlength=n_bins)
    den = np.bincount(bins, weights=(Xp * Xp).sum(0), minlength=n_bins)
    if shrink <= 0:
        return num / np.maximum(den, 1e-300)
    pri = np.median(den[den > 0]) if (den > 0).any() else 1.0
    return (num + shrink * pri) / np.maximum(den + shrink * pri, 1e-300)


def basal_bins(b: np.ndarray, n_bins: int) -> np.ndarray:
    q = np.quantile(b, np.linspace(0, 1, n_bins + 1)[1:-1])
    return np.clip(np.searchsorted(q, b), 0, n_bins - 1).astype(np.int64)
