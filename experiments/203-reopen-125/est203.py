"""203 — 125's hierarchical panel posterior, generalised to the PRODUCTION gene axis.

125 was fitted and scored on a 300x300 block: panel targets x panel GENES as responses.
Production zeroes every one of those 300 response coordinates (`062.prepare_main_source`
line 287 and `162.run_arms.effects_for`), and experiment 127 measured that they carry
2.283% of reference DE genes.  So the frozen 125 candidate cannot be scored on `pds` at
all -- it lives entirely inside production's null space.  To price it we must re-fit the
same estimator on the axis production actually ships, which is what this module does.

Everything else is verbatim 125:
  * per (source, target) row scale = median |effect| over that source's valid columns;
  * per-response-gene paired moments over the fit rows (latent = max(cov,0), noise =
    max(var - cov, 0)), EB-shrunk to the df-weighted global with pseudo-count k;
  * a variance floor at floor_fraction x global total;
  * Gaussian posterior mean, times 122's output scale (median contributing row scale);
  * posterior sign confidence |2 Phi(mean/sd) - 1|.

The fit rows honour 125's own frozen 5-fold target split: the moments used to predict a
target never see that target's fold.  The split is read from
`125/work/preregistered-freeze-v1/target-folds.npy`, so no new randomisation is introduced.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
F125 = ROOT / "experiments/125-hierarchical-panel-posterior/work/preregistered-freeze-v1"

SETTINGS = {
    "eb_balanced":     dict(pseudo=40.0, floor_frac=0.05, prior_scale=1.0),
    "eb_conservative": dict(pseudo=80.0, floor_frac=0.10, prior_scale=0.5),
    "eb_low_floor":    dict(pseudo=40.0, floor_frac=0.01, prior_scale=1.0),
    "eb_broad_prior":  dict(pseudo=40.0, floor_frac=0.05, prior_scale=2.0),  # 125's pick, 14/15
}


@dataclass
class Panel:
    """The two production sources for one PLOLO fold, on the official gene axis."""
    effects: np.ndarray        # (2, T, G) float64
    valid: np.ndarray          # (2, T, G) bool
    normalized: np.ndarray     # (2, T, G) float64
    row_scale: np.ndarray      # (2, T)    float64


def prepare(effects, available=None, support=None, valid_mask=None) -> Panel:
    """125 `transform_training`, on the production axis.  Row scale = median |effect|."""
    F, T, G = effects.shape
    eff = np.asarray(effects, dtype=np.float64)
    if valid_mask is not None:
        valid = np.asarray(valid_mask, dtype=bool).copy()
    else:
        valid = (available[:, :, None] & support[:, None, :]
                 & np.isfinite(eff) & (eff != 0.0))
    scale = np.zeros((F, T), dtype=np.float64)
    norm = np.zeros_like(eff)
    for f in range(F):
        for t in range(T):
            v = np.abs(eff[f, t, valid[f, t]])
            if v.size:
                scale[f, t] = float(np.median(v))
        usable = valid[f] & (scale[f, :, None] > 0.0)
        np.divide(eff[f], scale[f, :, None], out=norm[f], where=usable)
        valid[f] = usable
    return Panel(effects=eff, valid=valid, normalized=norm, row_scale=scale)


def output_scale(p: Panel):
    """122's output scale: per coordinate, median row scale among contributing families."""
    F, T, G = p.effects.shape
    s = np.broadcast_to(p.row_scale[:, :, None], p.valid.shape)
    masked = np.where(p.valid, s, np.nan)
    with np.errstate(invalid="ignore"):
        res = np.nanmedian(masked, axis=0)
    res = np.where(np.isfinite(res), res, 0.0)
    return res, p.valid.sum(axis=0).astype(np.int16)


def moments(p: Panel, fit_rows: np.ndarray):
    """125 `fit_raw_moments`, per response gene, generalised to F sources.

    For F = 2 this is bit-identical to 125: latent = max(cov_01, 0), noise_f =
    max(var_f - cov_01, 0) using the UNCLIPPED covariance, prior mean = mean of the
    family means.  For F = 3 (the kolf2 fold) `cov` is the mean of the three pairwise
    covariances, which is the standard one-way random-effects estimator and reduces to
    the F = 2 case exactly.  `test_est203.py` asserts that reduction.
    """
    F = p.effects.shape[0]
    pair = fit_rows[:, None].copy()
    for f in range(F):
        pair = pair & p.valid[f]
    n = pair.sum(axis=0).astype(np.float64)
    G = p.effects.shape[2]
    z = np.zeros(G)
    s = [np.where(pair, p.normalized[f], 0.0).sum(0) for f in range(F)]
    mean = [np.divide(s[f], n, out=z.copy(), where=n > 0) for f in range(F)]
    df = np.maximum(n - 1.0, 0.0)

    def centred_prod(a, b, sa, sb):
        return np.divide(
            np.where(pair, p.normalized[a] * p.normalized[b], 0.0).sum(0)
            - np.divide(sa * sb, n, out=z.copy(), where=n > 0),
            df, out=z.copy(), where=df > 0)

    var = [centred_prod(f, f, s[f], s[f]) for f in range(F)]
    pairs = [(a, b) for a in range(F) for b in range(a + 1, F)]
    cov = sum(centred_prod(a, b, s[a], s[b]) for a, b in pairs) / len(pairs)
    stable = df >= 2.0
    lat = np.where(stable, np.maximum(cov, 0.0), 0.0)
    noise = [np.where(stable, np.maximum(var[f] - cov, 0.0), 0.0) for f in range(F)]
    w = np.where(stable, df, 0.0)
    ws = float(w.sum())
    assert ws > 0.0, "no paired degrees of freedom"

    def g(v):
        return max(float(np.sum(w * v) / ws), 0.0)

    gl = g(lat)
    gn = [g(v) for v in noise]
    return dict(n=n, prior_mean_raw=sum(mean) / F, latent=lat, noise=noise,
                g_latent=gl, g_noise=gn,
                g_total=max(gl + float(np.mean(gn)), 1e-8))


def parameters(m, setting):
    k = setting["pseudo"]
    n = m["n"]
    mw = n / (n + k)
    cw = np.maximum(n - 1.0, 0.0) / (np.maximum(n - 1.0, 0.0) + k)
    floor = max(float(setting["floor_frac"] * m["g_total"]), 1e-10)
    lat = (cw * m["latent"] + (1 - cw) * m["g_latent"]) * setting["prior_scale"]
    noise = [np.maximum(cw * v + (1 - cw) * gv, floor)
             for v, gv in zip(m["noise"], m["g_noise"])]
    return dict(prior_mean=mw * m["prior_mean_raw"], latent=np.maximum(lat, floor),
                noise=noise, floor=floor)


def posterior(p: Panel, par):
    """125 `posterior_prediction`, over every observed source."""
    T, G = p.effects.shape[1], p.effects.shape[2]
    prec = np.broadcast_to(1.0 / par["latent"], (T, G)).copy()
    num = np.broadcast_to(par["prior_mean"] / par["latent"], (T, G)).copy()
    for f in range(p.effects.shape[0]):
        w = 1.0 / par["noise"][f]
        prec += p.valid[f] * w[None, :]
        num += p.valid[f] * p.normalized[f] * w[None, :]
    pvar = 1.0 / prec
    pmean = pvar * num
    scale, support = output_scale(p)
    allowed = support > 0
    effect = np.where(allowed, pmean * scale, 0.0)
    zz = pmean / np.sqrt(pvar)
    conf = np.where(allowed, np.abs(2.0 * ndtr(zz) - 1.0), 0.0)
    return effect, conf, support


def equal_mean_agreement(p: Panel):
    """122's comparator, verbatim (`125/run.py:565-601`), on the production axis."""
    support = p.valid.sum(axis=0).astype(np.int16)
    loc = np.divide(np.where(p.valid, p.normalized, 0.0).sum(0), support,
                    out=np.zeros(p.normalized.shape[1:]), where=support > 0)
    nzc = (p.valid & (p.normalized != 0.0)).sum(0)
    ssum = np.where(p.valid & (p.normalized != 0.0), np.sign(p.normalized), 0.0).sum(0)
    agree = np.divide(np.abs(ssum), nzc, out=np.zeros_like(loc), where=nzc > 0)
    scale, _ = output_scale(p)
    effect = np.where(support > 0, loc * scale, 0.0)
    conf = np.abs(loc) * (support / p.effects.shape[0]) * agree
    return effect, conf, support


def fold_assignment(targets) -> np.ndarray:
    """125's frozen 5-fold target split, restricted to this fold's targets."""
    panel = [line.strip() for line in (F125 / "panel-targets.txt").read_text().splitlines()
             if line.strip()]
    folds = np.load(F125 / "target-folds.npy")
    pos = {g: i for i, g in enumerate(panel)}
    return np.asarray([int(folds[pos[t]]) for t in targets], dtype=np.int64)


def nested_posterior(p: Panel, targets, setting):
    """Predict every target from moments fitted WITHOUT that target's frozen 125 fold."""
    fa = fold_assignment(targets)
    T = p.effects.shape[1]
    eff = np.zeros(p.effects.shape[1:], dtype=np.float64)
    conf = np.zeros_like(eff)
    sup = np.zeros(eff.shape, dtype=np.int16)
    for k in np.unique(fa):
        test = fa == k
        par = parameters(moments(p, ~test), setting)
        e, c, s = posterior(p, par)
        eff[test], conf[test], sup[test] = e[test], c[test], s[test]
    return eff, conf, sup
