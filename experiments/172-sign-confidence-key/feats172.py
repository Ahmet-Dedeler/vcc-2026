"""172 — the feature map.  ONE code path for training folds and for deployment.

`c(t, g)` must be a function of measurable scalars only: **no gene identity, no target
identity**, so that it transfers to targets never seen in training and to the unseen D/E/F
panel.  That is competitor `AgingLab`'s published shape ("a universal function of four
measurable scalars per (target, gene) — expression, saturation, consensus magnitude,
cross-source agreement — with no gene identity and no target identity"), extended here with
the donor-panel and cross-donor quantities this repo already has.

The single most important property of this file is that `build()` is called with exactly the
same code for a LODO fold (recipient = a held-out donor line) and for h1-test (recipient = the
h1-test reference).  A feature that could only be computed on one of the two would be a
training/serving skew, and it is the failure mode that kills transfer.

WHAT IS DELIBERATELY NOT IN HERE, AND WHY
-----------------------------------------
* **The whitened effect `H`.**  Production whitens with `union_pruned`, an operator fitted on
  donor perturbations.  There is no leak-free per-fold analogue (the operator would have to be
  refitted without the held-out line on that line's own gene axis), and AGENTS.md §4.2's
  operational corollary says the DE ordering signal is the **un-whitened** cross-donor
  agreement anyway.  So every feature here is un-whitened, and `|H|` is kept as the BASELINE
  KEY the learned key must beat, not as an input.  §"whitening" of RESULTS.md measures how much
  that costs.
* **Recipient control-cell detection rate and dispersion.**  Buildable for h1-test but not
  symmetrically for the donor lines without re-deriving three raw h5ads, and experiment 161
  measured that exactly these two scalars make a cross-line learned model *worse*
  (worst-line dpds -0.00097 -> +0.00382 when they are dropped).  Stated as unmeasured.

SCALE.  a15 ships `beta = 0.0`, so the transferred effect is `alpha * global * mean_d(renorm_d *
s_d[t])` for one build-wide scalar `alpha * global`.  Every magnitude feature here is therefore
either a within-target rank, a within-target log-ratio, or centred on the recipient's own
median -- none of them can carry the build's arbitrary amplitude into the model.

TIES (097 §10).  Every rank in this file is a MIDRANK.  The agreement vector is exactly zero
wherever a donor is silent or the donors disagree, which is a majority of most rows; ordinal
ranks would hand that block an arbitrary ramp ordered by position in the gene axis.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata

from donors172 import DonorHandle

EPS = 1e-12

FEATURES = [
    "rank_cons", "rank_geo", "rank_agree", "rank_z",
    "log_cons_rel", "log_geo_rel", "log_z_rel",
    "concord", "n_don", "rel_dis", "ratio_mm",
    "log_cpm_rec", "log_cpm_don", "log_marg", "log_sd", "resp",
    "sat",
    "t_cos", "t_lognorm", "t_concfrac", "t_ndon",
]


def _midrank01(v: np.ndarray) -> np.ndarray:
    """Midrank of `v` mapped to (0, 1]; 1 = largest.  Ties share their average rank."""
    if v.size == 0:
        return v
    return rankdata(v, method="average") / v.size


class Recipient:
    """Everything the feature map needs about the line the prediction is being made FOR."""

    def __init__(self, name: str, genes: np.ndarray, b_ctrl: np.ndarray,
                 targets: np.ndarray, live: np.ndarray):
        self.name = name
        self.genes = np.asarray([str(g) for g in genes])
        self.b_ctrl = np.asarray(b_ctrl, dtype=np.float64)
        self.targets = np.asarray([str(t) for t in targets])
        self.live = np.asarray(live, dtype=bool)     # CPM-gated: the coordinates that can move
        e = np.expm1(self.b_ctrl)
        self.cpm = 1e6 * e / max(e.sum(), EPS)


def build(rec: Recipient, source_names: tuple[str, ...], handles=None,
          verbose: bool = True) -> dict:
    """Feature tensors of shape (n_targets, n_genes) plus the consensus effect itself."""
    handles = handles or {n: DonorHandle(n) for n in source_names}
    T, G = len(rec.targets), len(rec.genes)
    live = rec.live

    S, cov, carried, renorm_in = {}, {}, {}, {}
    marg_abs, marg_sd, resp, cpm_don = {}, {}, {}, {}
    for n in source_names:
        h = handles[n]
        tcols, dcols = h.project(rec.genes)
        S[n], carried[n] = h.rows(rec.targets, rec.genes, tcols, dcols)
        c = np.zeros(G, bool)
        c[tcols] = True
        cov[n] = c
        renorm_in[n] = h.scale
        for dst, src in ((marg_abs, h.marg_abs), (marg_sd, h.marg_sd),
                         (resp, h.resp_frac), (cpm_don, h.cpm)):
            v = np.zeros(G)
            v[tcols] = src[dcols]
            dst[n] = v
        if verbose:
            print(f"    {n:<16} covers {c.sum():>6}/{G} of the recipient axis, "
                  f"carries {carried[n].sum():>5}/{T} targets")

    # MultiDonorModel's per-donor renormalisation, verbatim (multi_donor.py:98-101).
    ref = float(np.median(list(renorm_in.values())))
    renorm = {n: ref / (renorm_in[n] or 1.0) for n in source_names}
    for n in source_names:
        S[n] *= np.float32(renorm[n])
        marg_abs[n] *= renorm[n]
        marg_sd[n] *= renorm[n]

    # ---- per-(t, g) presence: donor carries the target AND covers the gene -------------
    P = {n: carried[n][:, None] & cov[n][None, :] for n in source_names}
    n_don = np.zeros((T, G), np.float32)
    for n in source_names:
        n_don += P[n]

    A = np.stack([np.where(P[n], S[n], 0.0) for n in source_names])          # (D, T, G)
    absA = np.abs(A)
    present = np.stack([P[n] for n in source_names])
    npres = np.maximum(n_don, 1.0)

    cons = A.sum(0) / npres                                    # production's uniform consensus
    with np.errstate(divide="ignore", invalid="ignore"):
        geo = np.exp(np.where(present, np.log(absA + 1e-300), 0.0).sum(0) / npres)
    geo = np.where(n_don > 0, geo, 0.0).astype(np.float32)
    sgn = np.where(present, np.sign(A), 0.0)
    unanimous = (np.abs(sgn.sum(0)) == n_don) & (n_don >= 2) & np.all(
        np.where(present, absA > 0, True), axis=0)
    concord = unanimous.astype(np.float32)
    agree = (geo * concord).astype(np.float32)                 # 013 / 095's best buildable key

    mx = absA.max(0)
    mn = np.where(n_don >= 2, np.where(present, absA, np.inf).min(0), 0.0)
    ratio_mm = np.where(n_don >= 2, mn / np.maximum(mx, EPS), 0.0).astype(np.float32)
    if len(source_names) == 2:
        a, b = (np.where(P[n], S[n], 0.0) for n in source_names)
        rel_dis = np.where(n_don >= 2, np.abs(a - b) / np.maximum(np.abs(a) + np.abs(b), EPS),
                           1.0).astype(np.float32)
    else:                                                       # generic: spread / sum
        rel_dis = np.where(n_don >= 2, 1.0 - ratio_mm, 1.0).astype(np.float32)

    # ---- gene-level donor-panel statistics, averaged over the donors that COVER g -------
    cov_stack = np.stack([cov[n] for n in source_names]).astype(np.float64)
    ncov = np.maximum(cov_stack.sum(0), 1.0)
    g_marg = np.stack([marg_abs[n] for n in source_names]).sum(0) / ncov
    g_sd = np.stack([marg_sd[n] for n in source_names]).sum(0) / ncov
    g_resp = np.stack([resp[n] for n in source_names]).sum(0) / ncov
    g_cpmd = np.stack([np.log10(cpm_don[n] + 1.0) for n in source_names]).sum(0) / ncov

    zcons = np.abs(cons) / np.maximum(g_sd[None, :], 1e-6)

    # ---- per-target scalars -----------------------------------------------------------
    t_ndon = np.zeros(T, np.float32)
    for n in source_names:
        t_ndon += carried[n]
    t_cos = np.zeros(T, np.float32)
    t_concfrac = np.zeros(T, np.float32)
    if len(source_names) == 2:
        u, v = (S[n] for n in source_names)
        both = cov[source_names[0]] & cov[source_names[1]] & live
        uu, vv = u[:, both].astype(np.float64), v[:, both].astype(np.float64)
        den = np.linalg.norm(uu, axis=1) * np.linalg.norm(vv, axis=1)
        t_cos = np.where(den > EPS, (uu * vv).sum(1) / np.maximum(den, EPS), 0.0)
        t_concfrac = (np.sign(uu) == np.sign(vv)).mean(1)
        bad = (t_ndon < 2)
        t_cos[bad] = 0.0
        t_concfrac[bad] = 0.0
    t_norm = np.linalg.norm(cons[:, live], axis=1)
    t_lognorm = np.log10(t_norm + 1e-9)
    t_lognorm = t_lognorm - np.median(t_lognorm)       # recipient-centred: no label needed

    # ---- assemble ---------------------------------------------------------------------
    def per_target_rank(M):
        R = np.zeros((T, G), np.float32)
        for i in range(T):
            R[i, live] = _midrank01(M[i, live].astype(np.float64))
        return R

    def per_target_logrel(M):
        R = np.zeros((T, G), np.float32)
        for i in range(T):
            v = np.abs(M[i, live].astype(np.float64))
            med = np.median(v[v > 0]) if np.any(v > 0) else 1.0
            R[i, live] = np.log10((v + 1e-9) / max(med, 1e-9))
        return R

    F = {}
    F["rank_cons"] = per_target_rank(np.abs(cons))
    F["rank_geo"] = per_target_rank(geo)
    F["rank_agree"] = per_target_rank(agree)
    F["rank_z"] = per_target_rank(zcons)
    F["log_cons_rel"] = per_target_logrel(cons)
    F["log_geo_rel"] = per_target_logrel(geo)
    F["log_z_rel"] = per_target_logrel(zcons)
    F["concord"] = concord
    F["n_don"] = n_don
    F["rel_dis"] = rel_dis
    F["ratio_mm"] = ratio_mm
    F["log_cpm_rec"] = np.broadcast_to(
        np.log10(rec.cpm + 1.0).astype(np.float32), (T, G))
    F["log_cpm_don"] = np.broadcast_to(g_cpmd.astype(np.float32), (T, G))
    lm = np.log10(g_marg + 1e-6)
    ls = np.log10(g_sd + 1e-6)
    F["log_marg"] = np.broadcast_to((lm - np.median(lm[live])).astype(np.float32), (T, G))
    F["log_sd"] = np.broadcast_to((ls - np.median(ls[live])).astype(np.float32), (T, G))
    F["resp"] = np.broadcast_to(g_resp.astype(np.float32), (T, G))
    F["sat"] = np.log1p(np.abs(cons) / (np.abs(rec.b_ctrl)[None, :] + 0.5)).astype(np.float32)
    F["t_cos"] = np.broadcast_to(t_cos.astype(np.float32)[:, None], (T, G))
    F["t_lognorm"] = np.broadcast_to(t_lognorm.astype(np.float32)[:, None], (T, G))
    F["t_concfrac"] = np.broadcast_to(t_concfrac.astype(np.float32)[:, None], (T, G))
    F["t_ndon"] = np.broadcast_to(t_ndon.astype(np.float32)[:, None], (T, G))

    return dict(F=F, cons=cons.astype(np.float32), geo=geo, agree=agree,
                n_don=n_don, live=live, carried_any=(t_ndon > 0))


def stack(F: dict, rows_t: np.ndarray, rows_g: np.ndarray) -> np.ndarray:
    """(n_rows, n_features) in the fixed FEATURES order."""
    return np.column_stack([np.asarray(F[k])[rows_t, rows_g] for k in FEATURES]).astype(
        np.float32)
