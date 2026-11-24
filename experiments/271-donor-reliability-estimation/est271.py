"""271 -- the NOISE-AWARE DONOR ESTIMATORS, frozen before any fold was measured (sha256 recorded in every sweep's meta).

Notation.  For a PLOLO fold with held line Y and donors j = 1..J (J = 2 on the three cancer folds), `row_j(t)` is donor
j's production un-whitened row for target t on the official axis (E062.prepare_main_source: (delta - bank mean) x
EFFECT_SCALE / raw_median; the K562 row carries the PerturbMatch 0.5 blend on the shared columns, exactly as the fold
cache's `pre` was built), `m(t)` the equal-weight mean over the carrying donors (= cd4svd5's base row), `A_j(t)` and
`B_j(t)` donor j's rows built from the two BATCH-PARITY halves of its own cells (267's half-banks, each half centred
on its own bank-wide mean, converted to the production row's units by the real bank's raw_median), and

    r_j(t)  = cos(A_j(t), B_j(t))                       the split-half reliability of the row at half depth
    R_j(t)  = 2 r / (1 + r)                              Spearman-Brown to the full depth (the row is the union of the halves)
    clipped to [R_MIN, R_MAX] = [0.005, 0.99] wherever it enters a weight (a reliability <= 0 means "no measurable
    response", not "negative information").

The estimators.  Every combination happens BEFORE whitening, at the one place a per-target, per-donor weight would
enter production (`build_v9_pre_whiten_effects` / `MultiDonorModel.effect(..., weights)`), and every combined row is
NORM-MATCHED to ||m(t)|| (263's `:nm`: the direction-only reading, so the global preserve and the residual sizing see
the incumbent's row norms and the arm is shape-only).  Targets carried by ONE donor keep the cached row (a weight over
one donor is 1).  Downstream is cd4svd5 verbatim (lib263: whitening through `plolo.Operator.apply`, CD4 deflated k = 5
re-solved at f = 0.34, the H1-pair union re-solved at g = 0.28, one add, one preserve).

  (a) relw       reliability-weighted mean:      w_j(t) ∝ R_j(t) / (1 - R_j(t)), normalised over the carrying donors.
                 (the BLUE of J measurements of one signal with noise variances sigma_j^2 = Var(s)(1 - R_j)/R_j)
      relw_eb    the same with R computed from EMPIRICAL-BAYES-shrunk half reliabilities
                     r~_j(t) = rbar_j + kappa_j (r_j(t) - rbar_j),  kappa_j = tau_j^2 / (tau_j^2 + s_j^2),
                     s_j^2 = Var_t cos(A_j(t), B_j(pi(t)))  (a fixed derangement pi: the sampling variance of ONE
                     split-half cosine under no target-specific signal),  tau_j^2 = max(Var_t r_j(t) - s_j^2, 0);
                 knob-free: kappa_j is estimated per donor per fold and reported.
      relw_glob  w_j ∝ Rbar_j / (1 - Rbar_j) with the DONOR-level mean reliability (constant over targets): the
                 scalar-reliability arm (019 / 142 / VCC_DONOR_WEIGHTS class), the baseline (a) must beat.
      relw_wht   (a) with r_j(t) measured AFTER the fold's operator (whitened halves, cosine on the operator columns).
  (b) ebrow      empirical-Bayes row shrinkage toward the consensus, then the average:
                     row~_j(t) = m(t) + R_j(t) (row_j(t) - m(t));   pre(t) = mean_j row~_j(t)
                 which is the weight vector w_j(t) = (1 + R_j(t) - Rbar(t)) / J_t, Rbar(t) = mean over the carrying
                 donors (exact identity, so it goes through the same `combine`).
      ebrow_eb   with the EB-shrunk reliabilities.
  (c) wiener_own per-gene Wiener shrinkage of each donor row with the row's OWN between-half variance (a diagonal
                 filter): with x = row_j(t, g) and d = A_j(t, g) - B_j(t, g),  sigma^2 = d^2 / 4  (Var((A+B)/2) =
                 Var(A - B) / 4 for independent halves; one degree of freedom), and the positive-part Wiener factor
                     w(t, g) = max(0, 1 - sigma^2 / x^2) = max(0, 1 - d^2 / (4 x^2)),   x' = w x;
                 then the equal-weight mean of the filtered rows, norm-matched.
      wiener_gene the same filter with the noise variance POOLED over the donor's targets per gene and scaled by the
                 target's cell count (Poisson: sigma^2 ∝ 1/n_t):  c_g = mean_t [ n_t d_tg^2 / 4 ],  sigma_tg^2 = c_g / n_t,
                 signal variance pooled per gene s_g^2 = max(mean_t x_tg^2 - mean_t sigma_tg^2, 0),
                     w(t, g) = s_g^2 / (s_g^2 + sigma_tg^2).
  (d) a13        the BASELINE: a13's per-gene donor-AGREEMENT damping (scripts/make_submission_pm_cd4.py
                 `agreement_shrink_vector`, gamma 0.25, rank scale 2000, toward zero), re-implemented: per target the
                 key  agree_g = (prod_j |row_j(g)|)^(1/J_t) if every carrying donor is nonzero and of one sign at g,
                 else 0;  alpha_g = exp(-0.25 * midrank_g / 2000) with the DESCENDING midrank over the whole official
                 axis (ties share a midrank; a row with < 2 distinct key values keeps alpha = 1); applied to the FINAL
                 composed effect row (where production applies it, at render: gating is a per-gene mask, so the order
                 does not matter).  `a13:nm` re-normalises each damped row to its undamped norm (shape-only, the
                 comparable reading); `a13:lit` is the literal amplitude-changing arm.
  (e) cd4svd5    equal weights = the incumbent (known answer at 0.0e+00).

Controls (frozen with the estimators):  `:shuf` = the donor's reliability vector (or its half-difference rows) permuted
across targets within the donor (seed SHUF_SEED); `relw:rev` = w_j ∝ (1 - R_j) / R_j (reversed weights); 263's random
pick (seed 0) and 263's whitened oracle (the ceiling) are re-run as reference rows.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

R_MIN, R_MAX = 0.005, 0.99
GAMMA_A13, RANK_SCALE_A13 = 0.25, 2000.0
SHUF_SEED = 20260906
DERANGE_SEED = 20260906


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------ reliabilities
def row_cos(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


def spearman_brown(r: np.ndarray) -> np.ndarray:
    """Two halves -> the full row (the union of the halves): R = 2r / (1 + r)."""
    r = np.asarray(r, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(np.isfinite(r), 2.0 * r / (1.0 + np.clip(r, -0.999, None)), np.nan)


def clip_r(R: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(R, dtype=np.float64), R_MIN, R_MAX)


def derangement(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    while True:
        p = rng.permutation(n)
        if n < 2 or not np.any(p == np.arange(n)):
            return p


def eb_shrink_half_reliability(r: np.ndarray, A: np.ndarray, B: np.ndarray, ok: np.ndarray, seed: int = DERANGE_SEED) -> dict:
    """r~ = rbar + kappa (r - rbar) over the targets `ok`; s^2 from a fixed derangement of the B rows (no
    target-specific signal), tau^2 = max(Var r - s^2, 0), kappa = tau^2 / (tau^2 + s^2)."""
    idx = np.flatnonzero(ok)
    r_ok = r[idx]
    pi = derangement(idx.size, seed)
    null = row_cos(A[idx], B[idx[pi]])
    s2 = float(np.nanvar(null, ddof=1))
    v = float(np.nanvar(r_ok, ddof=1))
    tau2 = max(v - s2, 0.0)
    kappa = tau2 / (tau2 + s2) if (tau2 + s2) > 0 else 0.0
    rbar = float(np.nanmean(r_ok))
    out = np.full(r.shape, np.nan)
    out[idx] = rbar + kappa * (r_ok - rbar)
    return dict(r_shrunk=out, rbar=rbar, var_r=v, s2_null=s2, tau2=tau2, kappa=float(kappa), null_mean=float(np.nanmean(null)),
                n=int(idx.size))


# ------------------------------------------------------------------------------------------ weights (a), (b)
def _fill(Rj: np.ndarray, avail_j: np.ndarray) -> np.ndarray:
    """A carrying donor with an undefined reliability gets the donor's mean over its defined targets."""
    R = Rj.copy()
    bad = avail_j & ~np.isfinite(R)
    if bad.any():
        R[bad] = np.nanmean(R[avail_j & np.isfinite(R)])
    return R


def weights(R_full: np.ndarray, avail: np.ndarray, mode: str) -> np.ndarray:
    """R_full: (J, n) Spearman-Brown reliabilities (NaN allowed); avail: (J, n) bool.  Returns W (n, J): NaN rows
    where fewer than two donors carry the target (the cached mean row is kept), else the per-target weights.
      mode 'relw'  : w_j ∝ R_j / (1 - R_j)           (a)
           'rev'   : w_j ∝ (1 - R_j) / R_j           reversed control
           'ebrow' : w_j = (1 + R_j - mean_k R_k) / J_t   (b), the exact weight form of the row shrinkage
           'glob'  : w_j ∝ Rbar_j / (1 - Rbar_j), Rbar_j = the donor's mean R over its carrying targets
    """
    J, n = avail.shape
    W = np.full((n, J), np.nan)
    R = np.stack([_fill(R_full[j], avail[j]) for j in range(J)])
    Rc = clip_r(R)
    if mode == "glob":
        rbar = np.array([float(np.nanmean(R[j][avail[j]])) for j in range(J)])
        rbar = clip_r(rbar)
    for t in np.flatnonzero(avail.sum(0) >= 2):
        js = np.flatnonzero(avail[:, t])
        if mode == "relw":
            w = Rc[js, t] / (1.0 - Rc[js, t])
        elif mode == "rev":
            w = (1.0 - Rc[js, t]) / Rc[js, t]
        elif mode == "ebrow":
            w = (1.0 + Rc[js, t] - Rc[js, t].mean()) / js.size
        elif mode == "glob":
            w = rbar[js] / (1.0 - rbar[js])
        else:
            raise ValueError(mode)
        W[t] = 0.0
        W[t, js] = w / w.sum()
    return W


def shuffle_within_donor(R_half: np.ndarray, avail: np.ndarray, seed: int = SHUF_SEED) -> np.ndarray:
    """The `:shuf` control: each donor's reliability vector permuted across ITS carrying targets."""
    rng = np.random.default_rng(seed)
    out = R_half.copy()
    for j in range(R_half.shape[0]):
        idx = np.flatnonzero(avail[j])
        out[j, idx] = R_half[j, idx][rng.permutation(idx.size)]
    return out


# ------------------------------------------------------------------------------------------ (c) diagonal filters
def wiener_own(x: np.ndarray, dAB: np.ndarray) -> np.ndarray:
    """w = max(0, 1 - d^2 / (4 x^2)) per coordinate; x' = w x.  x == 0 stays 0."""
    x2 = x * x
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(x2 > 0, 1.0 - (dAB * dAB) / (4.0 * x2), 0.0)
    return np.clip(w, 0.0, 1.0) * x


def wiener_gene(x: np.ndarray, dAB: np.ndarray, ncells: np.ndarray, avail_j: np.ndarray) -> tuple[np.ndarray, dict]:
    """Gene-pooled noise variance scaled by cell count, gene-pooled signal variance; w = s_g^2 / (s_g^2 + sigma_tg^2).
    x, dAB: (n, G) for one donor; ncells: (n,) the target's full cell count in the donor; avail_j: (n,) bool."""
    idx = np.flatnonzero(avail_j & (ncells > 0))
    v = (dAB[idx] ** 2) / 4.0                                        # (n_ok, G): one-df estimates of sigma_tg^2
    nt = ncells[idx].astype(np.float64)[:, None]
    c_g = (v * nt).mean(0)                                           # sigma_tg^2 = c_g / n_t
    sig2 = c_g[None, :] / nt
    s2_g = np.maximum((x[idx] ** 2).mean(0) - sig2.mean(0), 0.0)
    out = np.zeros_like(x)
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.where(s2_g[None, :] + sig2 > 0, s2_g[None, :] / (s2_g[None, :] + sig2), 0.0)
    out[idx] = w * x[idx]
    audit = dict(n_targets=int(idx.size), median_w=float(np.median(w)), mean_w=float(w.mean()),
                 share_genes_signal_zero=float((s2_g == 0).mean()))
    return out, audit


def combine_filtered(st: dict, rows_f: list[np.ndarray], *, norm_match: bool = True) -> np.ndarray:
    """lib263.combine's arithmetic for FILTERED rows: the equal-weight mean over the carrying donors, norm-matched to
    the cached mean row, float32-cast; single-donor targets are combined too (the filter applies to every row)."""
    pre = st["pre"].copy()
    avail = st["avail"]
    J, n = avail.shape
    acc = np.zeros_like(pre)
    cnt = np.zeros(n, dtype=np.int64)
    for j in range(J):
        m = avail[j]
        acc[m] += rows_f[j][m]
        cnt += m.astype(np.int64)
    cov = cnt > 0
    acc[cov] /= cnt[cov, None]
    if norm_match:
        na = np.linalg.norm(acc, axis=1)
        nm = np.linalg.norm(pre, axis=1)
        ok = cov & (na > 0) & (nm > 0)
        acc[ok] *= (nm[ok] / na[ok])[:, None]
    pre[cov] = acc[cov].astype(np.float32).astype(np.float64)
    return pre


# ------------------------------------------------------------------------------------------ (d) a13's damping
def descending_midrank(v: np.ndarray) -> np.ndarray:
    """0-based descending rank with ties sharing their mean position (the builder's `_descending_midrank`)."""
    order = np.argsort(-v, kind="stable")
    sv = v[order]
    n = v.size
    new_group = np.r_[True, sv[1:] != sv[:-1]]
    gid = np.cumsum(new_group) - 1
    first = np.flatnonzero(new_group)
    size = np.bincount(gid)
    mid = first + 0.5 * (size - 1)
    ranks = np.empty(n, dtype=np.float64)
    ranks[order] = mid[gid]
    return ranks


def a13_alpha(rows: list[np.ndarray], avail: np.ndarray, gamma: float = GAMMA_A13, rank_scale: float = RANK_SCALE_A13) -> tuple[np.ndarray, dict]:
    """Per-target per-gene damping factors from the cross-donor agreement key (see the module docstring)."""
    J, n = avail.shape
    G = rows[0].shape[1]
    alpha = np.ones((n, G))
    n_active = 0
    graded = []
    for t in range(n):
        js = np.flatnonzero(avail[:, t])
        if js.size < 2:
            continue
        X = np.stack([rows[j][t] for j in js])
        nz = np.all(X != 0, axis=0)
        same = np.all(np.sign(X) == np.sign(X[0])[None, :], axis=0)
        key = np.where(nz & same, np.exp(np.log(np.abs(np.where(nz, X, 1.0))).mean(0)), 0.0)
        if np.unique(key).size < 2:
            continue
        alpha[t] = np.exp(-gamma * descending_midrank(key) / rank_scale)
        n_active += 1
        graded.append(float((key > 0).mean()))
    return alpha, dict(n_active=n_active, median_graded_fraction=float(np.median(graded)) if graded else float("nan"),
                       gamma=gamma, rank_scale=rank_scale)


def apply_alpha(final: np.ndarray, alpha: np.ndarray, *, norm_match: bool) -> np.ndarray:
    out = final * alpha
    if norm_match:
        n0 = np.linalg.norm(final, axis=1)
        n1 = np.linalg.norm(out, axis=1)
        ok = (n0 > 0) & (n1 > 0)
        out[ok] *= (n0[ok] / n1[ok])[:, None]
    return out
