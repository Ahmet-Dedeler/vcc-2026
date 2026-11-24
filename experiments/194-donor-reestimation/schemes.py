"""194 — the estimator variants, expressed as PER-CELL WEIGHTS.

The incumbent donor estimator is

    S_T = sum_{c in T} x_c ;   delta_T = log1p(TS * S_T / |S_T|) - b_ctrl

i.e. a *pooled count sum*.  The key structural observation of this experiment is that
almost every standard alternative in the Perturb-seq literature is the SAME functional
with a per-cell weight in front of the sum:

    pooled           w_c = 1                       (incumbent)
    per-cell mean    w_c = 1 / depth_c             (mean of relative profiles, not of counts)
    cell QC filter   w_c in {0, 1}
    batch matching   w_c = f_ctrl(batch_c) / f_T(batch_c)
    knockdown        w_c = some decreasing function of the target's own residual expression
    guide averaging  w_c = 1 / (n_guides_T * n_cells_of_guide(c))

so K variants can be accumulated in ONE streaming pass over the cells by carrying K
weighted sums.  Only genuinely non-linear estimators (median/trimmed mean over guides,
MODZ) fall outside this, and section 2 of RESULTS.md records why they are unavailable
here anyway.

`w` is only ever used inside a ratio that is renormalised to TS, so the scale of `w`
within a target is irrelevant; only its shape across the target's cells matters.
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammaln

# knockdown model: a CRISPRi guide either fails (residual expression ~ 1x control) or
# works (the challenge's own criterion is >80% repression, i.e. residual ~ 0.2x).
KD_RESIDUAL = 0.2
KD_PRIOR = 0.8
# a per-cell knockdown readout is only meaningful if the target gene is expressed enough
# in the control that a single cell carries information.  Below this expected-count
# threshold the readout is dropout noise (which correlates with library depth, so
# weighting on it would import a depth confound), and every kd scheme falls back to w=1.
KD_MIN_EXPECTED = 0.5


def poisson_logpmf(k: np.ndarray, lam: np.ndarray) -> np.ndarray:
    lam = np.maximum(lam, 1e-12)
    return k * np.log(lam) - lam - gammaln(k + 1.0)


def kd_readout(obs: np.ndarray, expected: np.ndarray) -> np.ndarray:
    """Shrunken residual expression r_hat = (obs + 1) / (expected + 1).

    The +1 is one *control-scale* pseudo-count: when `expected` is large the estimate is
    the raw ratio, when it is small r_hat -> 1 and the cell looks un-knocked-down for
    everyone, so the weight becomes flat instead of becoming noise.
    """
    return (obs + 1.0) / (expected + 1.0)


def kd_posterior(obs: np.ndarray, expected: np.ndarray) -> np.ndarray:
    """P(cell is knocked down | its own target count), two-component Poisson mixture.

    Self-regularising: where `expected` is small both components have almost the same
    likelihood and the posterior collapses to the prior, i.e. a flat weight.
    """
    l_kd = poisson_logpmf(obs, expected * KD_RESIDUAL) + np.log(KD_PRIOR)
    l_no = poisson_logpmf(obs, expected) + np.log1p(-KD_PRIOR)
    m = np.maximum(l_kd, l_no)
    return np.exp(l_kd - m) / (np.exp(l_kd - m) + np.exp(l_no - m))


class WeightBook:
    """Per-cell weights for every scheme, computed once from the metadata table.

    Parameters
    ----------
    depth        (n,) library size per cell
    pct_mt       (n,) mitochondrial percentage
    n_genes      (n,) genes detected
    batch        (n,) batch index
    tcode        (n,) target index, -1 for cells not accumulated, `ctrl_code` for control
    obs_target   (n,) the cell's own count of its target gene, NaN where unmeasurable
    exp_target   (n,) expected count of that gene under the control profile, at this
                      cell's depth; NaN where unmeasurable
    """

    def __init__(self, *, depth, pct_mt, n_genes, batch, tcode, ctrl_code,
                 obs_target, exp_target, seed=20260903):
        self.n = len(depth)
        self.depth = np.asarray(depth, dtype=np.float64)
        self.pct_mt = np.asarray(pct_mt, dtype=np.float64)
        self.n_genes = np.asarray(n_genes, dtype=np.float64)
        self.batch = np.asarray(batch, dtype=np.int64)
        self.tcode = np.asarray(tcode, dtype=np.int64)
        self.ctrl_code = int(ctrl_code)
        self.obs = np.asarray(obs_target, dtype=np.float64)
        self.exp = np.asarray(exp_target, dtype=np.float64)
        self.rng = np.random.default_rng(seed)
        self.live = self.tcode >= 0
        self.is_ctrl = self.tcode == self.ctrl_code
        self.audit: dict[str, dict] = {}
        self._build()

    # ------------------------------------------------------------------ helpers
    def _norm_within_target(self, w: np.ndarray) -> np.ndarray:
        """Rescale so every target's weights average 1.  Pure bookkeeping (the profile is
        renormalised anyway) but it keeps the effective-cell-count audit interpretable and
        stops a target with a degenerate weight vector from underflowing."""
        w = np.where(self.live, w, 0.0)
        ng = int(self.tcode.max()) + 1
        s = np.bincount(self.tcode[self.live], weights=w[self.live], minlength=ng)
        n = np.bincount(self.tcode[self.live], minlength=ng).astype(np.float64)
        scale = np.where(s > 0, n / np.maximum(s, 1e-300), 0.0)
        # a target whose weights all vanished falls back to uniform rather than to zero
        dead = s <= 0
        out = np.zeros_like(w)
        out[self.live] = w[self.live] * scale[self.tcode[self.live]]
        if dead.any():
            m = self.live & dead[self.tcode]
            out[m] = 1.0
        return out

    def _shuffle_within(self, dest: np.ndarray, src: np.ndarray,
                        idx: np.ndarray, key: np.ndarray) -> None:
        """dest[block] <- a random permutation of src[block], for each block of `key`."""
        o = np.argsort(key, kind="stable")
        idx_s, key_s = idx[o], key[o]
        uk = np.unique(key_s)
        st = np.searchsorted(key_s, uk, side="left")
        en = np.searchsorted(key_s, uk, side="right")
        for a, b in zip(st, en):
            blk = idx_s[a:b]
            dest[blk] = src[self.rng.permutation(blk)]

    def _kd_applicable(self) -> np.ndarray:
        """Cells of targets whose own gene is expressed enough for a per-cell readout.

        Applicability is decided PER TARGET, on the target's median expected count, not
        per cell.  A per-cell threshold on `exp = p_ctrl[g] * depth_c` would leave the
        shallow cells of an applicable target un-weighted and the deep ones down-weighted,
        which silently builds a library-depth preference into the weighting — the exact
        confound that broke `data/harness/orion-hct116.h5ad`.  The first version of this
        file did that and the depth audit caught it (`kd_perm` shifted mean depth to
        0.947x, i.e. as much as `kd_soft` itself).
        """
        cand = self.live & ~self.is_ctrl & np.isfinite(self.exp)
        ng = int(self.tcode.max()) + 1
        med = np.zeros(ng)
        idx = np.flatnonzero(cand)
        o = np.argsort(self.tcode[idx], kind="stable")
        idx_s = idx[o]
        t_s = self.tcode[idx_s]
        ut = np.unique(t_s)
        st = np.searchsorted(t_s, ut, side="left")
        en = np.searchsorted(t_s, ut, side="right")
        for t, a, b in zip(ut, st, en):
            med[t] = np.median(self.exp[idx_s[a:b]])
        ok_target = med >= KD_MIN_EXPECTED
        self._kd_median_expected = med
        return cand & ok_target[self.tcode]

    # -------------------------------------------------------------------- build
    def _build(self) -> None:
        w: dict[str, np.ndarray] = {}
        live = self.live

        w["v0"] = np.where(live, 1.0, 0.0)

        # --- candidate 5a: relative-profile mean instead of pooled counts -------
        w["cellmean"] = np.where(live, 1.0 / np.maximum(self.depth, 1.0), 0.0)

        # --- candidate 5b: cell-quality filtering ------------------------------
        d_lo, d_hi = np.percentile(self.depth[live], [2.0, 99.5])
        mt_hi = np.percentile(self.pct_mt[live], 97.5)
        g_lo = np.percentile(self.n_genes[live], 2.0)
        keep = live & (self.depth >= d_lo) & (self.depth <= d_hi) & \
            (self.pct_mt <= mt_hi) & (self.n_genes >= g_lo)
        w["qc"] = keep.astype(np.float64)
        self.audit["qc"] = dict(depth_lo=float(d_lo), depth_hi=float(d_hi),
                                mt_hi=float(mt_hi), genes_lo=float(g_lo),
                                kept=float(keep[live].mean()))

        # --- candidate 4: batch-composition matching to the control ------------
        nb = int(self.batch.max()) + 1
        ng = int(self.tcode.max()) + 1
        cb = np.bincount(self.batch[self.is_ctrl], minlength=nb).astype(np.float64)
        f_ctrl = cb / max(cb.sum(), 1.0)
        flat = self.tcode[live] * nb + self.batch[live]
        tb = np.bincount(flat, minlength=ng * nb).astype(np.float64).reshape(ng, nb)
        tb_frac = tb / np.maximum(tb.sum(1, keepdims=True), 1.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(tb_frac > 0, f_ctrl[None, :] / tb_frac, 0.0)
        bw = np.zeros(self.n)
        bw[live] = ratio[self.tcode[live], self.batch[live]]
        bw[self.is_ctrl] = 1.0
        # cap so a target with one cell in a batch the control barely used cannot be
        # dominated by that cell
        bw = np.minimum(bw, 20.0)
        w["batchmatch"] = self._norm_within_target(bw)
        self.audit["batchmatch"] = dict(
            n_batches=nb, max_weight=float(w["batchmatch"][live].max()),
            control_batches_used=int((cb > 0).sum()))

        # --- candidate 1: knockdown weighting / filtering -----------------------
        app = self._kd_applicable()
        r = np.ones(self.n)
        r[app] = kd_readout(self.obs[app], self.exp[app])
        post = np.full(self.n, KD_PRIOR)
        post[app] = kd_posterior(self.obs[app], self.exp[app])

        soft = np.where(live, 1.0, 0.0)
        soft[app] = np.clip(1.0 - r[app], 0.0, 1.0)
        w["kd_soft"] = self._norm_within_target(soft)

        pw = np.where(live, 1.0, 0.0)
        pw[app] = post[app]
        w["kd_post"] = self._norm_within_target(pw)

        # Hard gates are placed on the POSTERIOR, not on the shrunken ratio.  A cell that
        # is genuinely repressed usually carries obs = 0, and (0+1)/(exp+1) is ~0.5 at
        # exp ~ 1, so a ratio gate at 0.2 discards essentially every cell.  The posterior
        # gate is the challenge's own ">80% knockdown" criterion read probabilistically.
        for tag, thr in (("kd_gate50", 0.5), ("kd_gate80", 0.8)):
            h = np.where(live, 1.0, 0.0)
            h[app] = (post[app] >= thr).astype(np.float64)
            w[tag] = self._norm_within_target(h)

        # --- falsifier 1: weights permuted within target -----------------------
        idx = np.flatnonzero(app)
        perm = pw.copy()
        self._shuffle_within(perm, pw, idx, self.tcode[idx])
        w["kd_perm"] = self._norm_within_target(perm)

        # --- falsifier 2: permuted within (target, depth decile) ---------------
        # A cell-level knockdown readout is confounded with library depth (a shallow cell
        # shows fewer target counts whether or not the guide worked).  Permuting inside a
        # depth stratum keeps the weight/depth coupling and destroys only the
        # knockdown-specific part, so `kd_soft` must beat THIS, not just `kd_perm`.
        dec = np.zeros(self.n, dtype=np.int64)
        qs = np.percentile(self.depth[live], np.arange(10, 100, 10))
        dec[live] = np.searchsorted(qs, self.depth[live])
        permd = pw.copy()
        self._shuffle_within(permd, pw, idx, self.tcode[idx] * 10 + dec[idx])
        w["kd_perm_depth"] = self._norm_within_target(permd)

        # --- the composite -------------------------------------------------------
        w["combo"] = self._norm_within_target(w["kd_post"] * w["qc"] * w["batchmatch"])

        self.audit["kd"] = dict(
            applicable_cells=float(app.sum() / max(live.sum() - self.is_ctrl.sum(), 1)),
            applicable_targets=int(len(np.unique(self.tcode[app]))),
            median_expected=float(np.median(self.exp[app])) if app.any() else 0.0,
            mean_soft_weight=float(soft[app].mean()) if app.any() else 0.0,
            mean_posterior=float(post[app].mean()) if app.any() else 0.0,
            gate50_kept=float((post[app] >= 0.5).mean()) if app.any() else 0.0,
            gate80_kept=float((post[app] >= 0.8).mean()) if app.any() else 0.0,
            median_residual=float(np.median(r[app])) if app.any() else 0.0,
        )
        # depth confound audit: does a scheme systematically prefer shallow cells?
        base_depth = float(np.average(self.depth[live & ~self.is_ctrl]))
        for k, v in w.items():
            m = live & ~self.is_ctrl
            ww = v[m]
            self.audit.setdefault("depth_shift", {})[k] = float(
                np.average(self.depth[m], weights=ww) / base_depth) if ww.sum() > 0 else 0.0
        self.w = w

    def names(self) -> list[str]:
        return list(self.w)
