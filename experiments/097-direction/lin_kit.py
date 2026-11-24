"""Shared machinery for the Ahlmann-Eltze linear-transfer probe (exp097, `lin_*` arms).

Everything is reported through the SAME two dimensionless shape quantities `exact.py` defines:

    kappa = X / sqrt(B * sumD)        B = cell_eval2.metrics.delta._across_pert_budget
    eps   = (A - B) / B

with `MSE_min = 1 + [Kc + claim*eps - 2*sqrt(claim*sumD)*kappa] / sumD` after the global
amplitude has been re-optimised onto the rho cliff.  kappa is scale-free in the prediction, so
no arm here ever tunes an amplitude.

Two extras beyond `exact.py`:

* `shape()` accepts a ROW SUBSET, so a supervised arm can be scored on a held-out half of the
  perturbations with the scorer's own weights, budget and denominator restricted to those rows.
  `sumD_S = sum_{p in S} w_p (r_p - C_p - C_ctrl)`, exactly the rows the scorer would emit.
* `donor_block()` memory-maps a donor `delta.npy` and slices to the harness gene axis BEFORE
  materialising anything, because the two Orion donors are 2.8 GB dense each.

`_across_pert_budget` is invariant to whatever sits at each row's own target coordinate (the
coordinate is dropped from its sums and only enters through a shift that cancels exactly), so
zeroing dh[k, cols[k]] before the einsums -- which the scorer's distance requires -- does not
disturb B.
"""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/src")

from decompose import terms  # noqa: E402
from msekit import group_stats  # noqa: E402

from cell_eval2.metrics.delta import _across_pert_budget  # noqa: E402

HARNESS = {
    "h1-test": (
        f"{ROOT}/data/harness/h1-test.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        "h1_2025",          # BANNED donor: this harness IS that dataset
    ),
    "orion": (
        f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
        f"{ROOT}/experiments/095-clean-harness/out/runs/"
        f"g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        "orion_hct116",     # BANNED donor: this harness IS that cell line
    ),
}

DONOR_DIR = f"{ROOT}/data/donors"


# --------------------------------------------------------------------------------------- harness

def prep(name: str):
    """`terms(...)` for one harness, with the on-target coordinates zeroed and targets attached."""
    real, pred, banned = HARNESS[name]
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    assert np.array_equal(rg, pg)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0
        t["d"][k, c] = 0.0
    t["targets"] = np.array([str(x) for x in pp[t["keep"]]])
    t["genes"] = rg
    t["banned"] = banned
    t["name"] = name
    # per-row denominator contribution, so a row subset can be scored exactly
    t["drow"] = (t["r"] - t["Cp"] - t["Cctrl"]) * t["w"]
    return t


def shape_full(DH, t, rows=None):
    """Every scored quantity for a candidate effect matrix, optionally on a row subset."""
    P = t["dh"].shape[0]
    idx = np.arange(P) if rows is None else np.asarray(rows, np.intp)
    kap, eps, X, A, B, sumD = shape(DH, t, rows)
    claim = float((t["capped"] * t["w"])[idx].sum())
    Kc = float((t["Cctrl"] * t["w"])[idx].sum())
    m, s = mse_min(X, A, B, sumD, claim, Kc)
    return dict(kappa=kap, eps=eps, X=X, A=A, B=B, sumD=sumD, claim=claim, Kc=Kc,
                cos=X / np.sqrt(max(A, 1e-300) * sumD), mse_min=m, s_star=s, n=idx.size)


def shape(DH, t, rows=None):
    """`(kappa, eps, X, A, B, sumD)` for a candidate effect matrix, optionally on a row subset.

    `DH` is always full-panel `[P, G]`; `rows` selects which perturbations are scored.
    """
    P = t["dh"].shape[0]
    idx = np.arange(P) if rows is None else np.asarray(rows, np.intp)
    w = t["w"][idx]
    D = t["d"][idx]
    H = DH[idx]
    cols = {}
    for new, old in enumerate(idx):
        if int(old) in t["cols"]:
            cols[new] = t["cols"][int(old)]
    H = H.copy()
    for k, c in cols.items():
        H[k, c] = 0.0
    X = float((np.einsum("ij,ij->i", H, D) * w).sum())
    A = float((np.einsum("ij,ij->i", H, H) * w).sum())
    B = float(_across_pert_budget(H, cols, w))
    sumD = float(t["drow"][idx].sum())
    if B <= 1e-12 * max(A, 1e-300) or sumD <= 0:
        # A perturbation-constant prediction has B == 0 EXACTLY: `_across_pert_budget` measures
        # across-perturbation spread, so a broadcast vector is granted no correction at all
        # (rho = B/claim = 0).  kappa and eps are both singular there; `mse_min` still is not.
        return np.nan, np.nan, X, A, B, sumD
    return X / np.sqrt(B * sumD), (A - B) / B, X, A, B, sumD


def mse_min(X, A, B, sumD, claim, Kc):
    """`min_s [ s^2 A - 2 s X + Kc - min(claim, s^2 B) ] / sumD + 1`, evaluated EXACTLY.

    `exact.py` asserts the minimum sits on the rho cliff `s_c = sqrt(claim/B)`; that is true only
    while `X/A <= s_c <= X/(A-B)`.  Both branches are checked here so a degenerate arm (`B = 0`,
    i.e. a constant predictor, whose cliff is at infinity) is still scored honestly.
    """
    def f(s):
        return 1.0 + (s * s * A - 2 * s * X + Kc - min(claim, s * s * B)) / sumD
    cands = [0.0]
    if A > 0:
        cands.append(X / A)                                   # above the cliff
    if A - B > 0:
        cands.append(X / (A - B))                             # below the cliff
    if B > 0:
        cands.append(np.sqrt(claim / B))                      # the cliff itself
    cands = [c for c in cands if np.isfinite(c) and c >= 0]
    vals = [f(c) for c in cands]
    i = int(np.argmin(vals))
    return float(vals[i]), float(cands[i])


# ---------------------------------------------------------------------------------------- donors

_CACHE: dict = {}


def donor_block(donor: str, genes: np.ndarray, *, decompose: bool = True):
    """Donor effects projected onto `genes` (the harness gene axis), memory-lean.

    Returns `dict(S[Pd, G], m[G], perts, index, hit[G])` where `S` is the gene-specific part
    (`delta - panel mean`) when `decompose`, else `delta` itself, laid out on the TARGET gene
    axis with unmatched target genes left at exactly 0.  `hit` marks the target columns the
    donor actually covers.

    The panel mean is taken column-wise, so restricting columns first is identical to
    `Donor.decompose()` followed by the same restriction.
    """
    key = (donor, genes.size, str(genes[0]), str(genes[-1]), decompose)
    if key in _CACHE:
        return _CACHE[key]
    d = f"{DONOR_DIR}/{donor}"
    dgenes = np.asarray(open(f"{d}/genes.txt").read().split("\n"))
    dperts = np.asarray(open(f"{d}/perts.txt").read().split("\n"))
    A = np.load(f"{d}/delta.npy", mmap_mode="r")
    assert A.shape == (dperts.size, dgenes.size), (A.shape, dperts.size, dgenes.size)
    gi = {str(g): i for i, g in enumerate(dgenes)}
    src = np.full(genes.size, -1, np.int64)
    for j, g in enumerate(genes):
        k = gi.get(str(g))
        if k is not None:
            src[j] = k
    hit = src >= 0
    take = src[hit]
    order = np.argsort(take)                      # mmap fancy-index wants ascending columns
    S = np.zeros((dperts.size, genes.size), np.float32)
    cols_target = np.nonzero(hit)[0][order]
    tk = take[order]
    for lo in range(0, dperts.size, 2048):        # chunk rows: the Orion blocks are 2.8 GB
        hi = min(lo + 2048, dperts.size)
        S[lo:hi, cols_target] = np.asarray(A[lo:hi])[:, tk]
    del A
    m = S.mean(axis=0)
    if decompose:
        S -= m[None, :]
        S[:, ~hit] = 0.0                          # keep uncovered genes at exactly zero
    out = dict(S=S, m=m, perts=dperts, index={str(p): i for i, p in enumerate(dperts)},
               hit=hit, n_genes_hit=int(hit.sum()))
    _CACHE[key] = out
    return out


def embed_pca(S: np.ndarray, kmax: int, seed: int = 0):
    """Rank-`kmax` PCA of a donor effect matrix `S[Pd, G]` (rows already centred by `donor_block`).

    Returns `(U*Sv [Pd, kmax], V [G, kmax], evr[kmax], total_var)`.  Randomised SVD so the
    17,977 x 18,533 Orion blocks stay tractable.
    """
    from sklearn.utils.extmath import randomized_svd
    U, sv, Vt = randomized_svd(S, n_components=kmax, n_iter=7, random_state=seed)
    tot = float(np.einsum("ij,ij->", S, S))
    evr = (sv ** 2) / max(tot, 1e-30)
    return (U * sv).astype(np.float64), Vt.T.astype(np.float64), evr, tot


def ridge_fit(E: np.ndarray, Y: np.ndarray, lam: float, intercept: bool):
    """Ridge `W = argmin ||E W - Y||^2 + lam ||W||^2`, with an optional un-penalised intercept.

    `E [n, k]`, `Y [n, G]`.  Returns `(W [k, G], c [G])`.  The intercept is the control-centred
    constant vector -- and since `d_p` is already control-centred, that constant IS the panel
    mean response, i.e. exactly the baseline Ahlmann-Eltze's linear model has to beat.
    """
    E = np.asarray(E, np.float64)
    Y = np.asarray(Y, np.float64)
    if intercept:
        ec, yc = E.mean(0), Y.mean(0)
        Ec, Yc = E - ec, Y - yc
    else:
        ec = np.zeros(E.shape[1])
        yc = np.zeros(Y.shape[1])
        Ec, Yc = E, Y
    k = E.shape[1]
    G = Ec.T @ Ec + lam * np.eye(k)
    W = np.linalg.solve(G, Ec.T @ Yc)
    c = yc - ec @ W
    return W, c


# --------------------------------------------------------- mapping a shape onto the board panel

# exact.py's board-anchored constants (MSE-MECHANISM.md (a), a12's emission policy at the #247
# cap).  `kappa0` is the incumbent's own deconfounded cosine on the challenge panel, so an arm
# measured as `M = kappa_arm / kappa_incumbent` on a harness maps to `kappa0 * M` there.
CH = dict(claim=0.24653, Kc=0.005410, sumD=0.228572, kappa0=0.07402, B=0.280503)
B_ANCHOR, SPAN = 0.987863, 0.950163


def challenge_raw(M, eps):
    c, sD = CH["claim"], CH["sumD"]
    return 1 + (CH["Kc"] + c * eps - 2 * np.sqrt(c * sD) * CH["kappa0"] * M) / sD


def challenge_score(M, eps):
    return max(0.0, min(1.0, (B_ANCHOR - challenge_raw(M, eps)) / SPAN))


def r2(Y, P):
    """Fraction of the total sum of squares of `Y` explained by `P` (can go negative)."""
    ss = float(np.einsum("ij,ij->", Y, Y))
    re = float(np.einsum("ij,ij->", Y - P, Y - P))
    return 1.0 - re / max(ss, 1e-30)
