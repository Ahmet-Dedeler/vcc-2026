"""189 — shared loaders and the per-gene cross-line reliability statistic.

WHAT THIS FILE IS
-----------------
The definition of `r_g`, the per-gene cross-line sign-agreement rate, plus the donor-table
plumbing needed to compute it.  Everything else in this experiment imports from here so that
one definition of the statistic is used in the pilot reproduction, in the LOLO build, in the
three priced uses and in every falsifier.

THE STATISTIC
-------------
For an ordered pair of donor lines (A, B) sharing a target set T and a gene axis G:

    1. each matrix is centred ACROSS TARGETS   (A <- A - mean_t A)   -- production removes the
       shared programme, and a target-independent addition is the metric's zero (161 §1).
    2. per gene, per line, a gate at that line's OWN median |effect| over the shared targets:
           gA[t,g] = |A[t,g]| > median_t |A[.,g]|
       so exactly half of each line's entries pass, by construction, independent of the gene's
       overall responsiveness.  Genes are therefore NOT selected by magnitude.
    3. supported pairs are the entries where BOTH lines pass their own gate.
    4. r_g = fraction of supported pairs with sign(A) == sign(B).

Under independence the support is ~25% of |T|, and `r_g` is a binomial proportion with a known
`n_g` -- which is what makes the empirical-Bayes shrinkage in `eb_shrink` legitimate.

DONOR TABLES.  `b = log1p(TS * P / sum P)`, `delta_p = b_p - b_ctrl`, TS = 5e4 (161 common.py).
`delta` is memory-mapped; orion is 2.7 GB per line.
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")

# table id -> (path, cell LINE group, human label).  Line groups matter, not tables:
# k562_gwps and k562_essential are two libraries of ONE line (161 common.py).
TABLES = {
    "k562_gwps":      ("data/donors/k562_gwps",      "K562",    "K562 (Replogle gwps)"),
    "k562_essential": ("data/donors/k562_essential", "K562",    "K562 (Replogle essential)"),
    "rpe1":           ("data/donors/rpe1",           "RPE1",    "RPE1"),
    "nadig_jurkat":   ("data/donors/nadig_jurkat",   "JURKAT",  "Jurkat"),
    "orion_hct116":   ("data/donors/orion_hct116",   "HCT116",  "HCT116 (X-Atlas/Orion)"),
    "orion_hek293t":  ("data/donors/orion_hek293t",  "HEK293T", "HEK293T (X-Atlas/Orion)"),
    "h1_2025":        ("data/donors/h1_2025",        "H1",      "H1 hESC (VCC2025) -- BANNED"),
}

# AGENTS.md trap #1: h1_2025 covers 85/85 of h1-test's targets.  Never opened here.
BANNED = {"h1_2025"}


def panel_targets() -> set[str]:
    with open(os.path.join(REPO, "data/challenge/pert_counts.csv")) as fh:
        return {r["target_gene"] for r in csv.DictReader(fh)}


def load_table(tid: str) -> dict | None:
    assert tid not in BANNED, f"{tid} is BANNED (oracle leak on h1-test)"
    p = os.path.join(REPO, TABLES[tid][0])
    if not os.path.exists(os.path.join(p, "delta.npy")):
        return None
    perts = [l.strip() for l in open(os.path.join(p, "perts.txt")) if l.strip()]
    genes = [l.strip() for l in open(os.path.join(p, "genes.txt")) if l.strip()]
    delta = np.load(os.path.join(p, "delta.npy"), mmap_mode="r")
    b = np.load(os.path.join(p, "b_ctrl.npy"))
    assert delta.shape == (len(perts), len(genes))
    meta = {}
    mp = os.path.join(p, "meta.json")
    if os.path.exists(mp):
        meta = json.load(open(mp))
    return dict(id=tid, line=TABLES[tid][1], label=TABLES[tid][2], perts=perts, genes=genes,
                delta=delta, b=b, meta=meta)


def submatrix(tab: dict, targets, genes, dtype=np.float32) -> np.ndarray:
    """(targets x genes) block from a memory-mapped delta, rows read in ascending order."""
    pi = {p: i for i, p in enumerate(tab["perts"])}
    gi = {g: i for i, g in enumerate(tab["genes"])}
    ridx = np.asarray([pi[t] for t in targets])
    cidx = np.asarray([gi[g] for g in genes])
    out = np.empty((len(targets), len(genes)), dtype=dtype)
    for k in np.argsort(ridx):
        out[k] = np.nan_to_num(np.asarray(tab["delta"][ridx[k]])[cidx], copy=False)
    return out


def bvec(tab: dict, genes) -> np.ndarray:
    gi = {g: i for i, g in enumerate(tab["genes"])}
    return tab["b"][np.asarray([gi[g] for g in genes])].astype(np.float64)


def centre(A: np.ndarray) -> np.ndarray:
    """Remove the target-independent mean -- production's across-target centring."""
    return A - A.mean(axis=0, keepdims=True)


def gate_mask(A: np.ndarray) -> np.ndarray:
    """Per-gene gate at that line's own median |effect| over the shared targets."""
    a = np.abs(A)
    return a > np.median(a, axis=0, keepdims=True)


def pair_agreement(A: np.ndarray, B: np.ndarray, gA=None, gB=None):
    """(k_g, n_g) -- agreeing supported pairs and supported pairs, per gene.

    A, B must already be centred and share (targets x genes).
    """
    if gA is None:
        gA = gate_mask(A)
    if gB is None:
        gB = gate_mask(B)
    sup = gA & gB
    agree = sup & ((A > 0) == (B > 0))
    return agree.sum(axis=0).astype(np.float64), sup.sum(axis=0).astype(np.float64)


def eb_shrink(k: np.ndarray, n: np.ndarray):
    """Empirical-Bayes posterior mean of a binomial proportion under a Beta(a, b) prior.

    `r_g` is a binomial proportion with a KNOWN `n_g`, so the honest estimate is the posterior
    mean, not the raw fraction.  (a, b) are fitted by maximum marginal (beta-binomial)
    likelihood over the genes -- the standard EB estimator, no moment approximation.

    Returns (post, a, b, mom) where `mom` is the method-of-moments starting point, kept so the
    fit can be sanity-checked against it.
    """
    from scipy.optimize import minimize
    from scipy.special import betaln

    ok = n > 0
    kk, nn = k[ok].astype(np.float64), n[ok].astype(np.float64)
    p = kk / nn
    mu = float(np.sum(kk) / np.sum(nn))
    var_obs = float(np.average((p - mu) ** 2, weights=nn))
    var_samp = float(np.average(mu * (1 - mu) / nn, weights=nn))
    var_true = max(var_obs - var_samp, 1e-8)
    nu0 = max(mu * (1 - mu) / var_true - 1.0, 1e-2)
    mom = (mu * nu0, (1 - mu) * nu0)

    def nll(theta):
        a, b = np.exp(theta)
        return -float(np.sum(betaln(kk + a, nn - kk + b) - betaln(a, b)))

    res = minimize(nll, np.log(np.maximum(mom, 1e-3)), method="Nelder-Mead",
                   options=dict(xatol=1e-8, fatol=1e-8, maxiter=4000))
    a, b = np.exp(res.x)
    post = (k + a) / (n + a + b)
    return post, float(a), float(b), mom


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    from scipy.stats import rankdata
    rx, ry = rankdata(x), rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else 0.0


def partial_rank(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Residual of rank(y) after least-squares removal of rank(X) columns (with intercept)."""
    from scipy.stats import rankdata
    ry = rankdata(y).astype(np.float64)
    RX = np.column_stack([rankdata(X[:, j]) for j in range(X.shape[1])]).astype(np.float64)
    RX = np.column_stack([np.ones(len(ry)), RX])
    beta, *_ = np.linalg.lstsq(RX, ry, rcond=None)
    return ry - RX @ beta


def decile_spread(key: np.ndarray, value: np.ndarray, w: np.ndarray | None = None):
    """Top vs bottom decile of `key`, reading the support-weighted mean of `value`."""
    n = len(key)
    order = np.argsort(key, kind="stable")
    d = max(1, n // 10)
    lo, hi = order[:d], order[-d:]
    if w is None:
        a, b = float(value[hi].mean()), float(value[lo].mean())
    else:
        a = float(np.sum(w[hi] * value[hi]) / np.sum(w[hi]))
        b = float(np.sum(w[lo] * value[lo]) / np.sum(w[lo]))
    return a, b, a - b
