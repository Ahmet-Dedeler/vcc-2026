"""Canonical covariance-mixture empirical Bayes in three-condition space."""

from itertools import product

import numpy as np
from scipy.special import logsumexp

TAUS = (.05, .1, .2, .4, .8, 1.6, 3.2)


def components(diagonal=False):
    result = [{"kind": "zero", "tau": 0.}]
    result += [{"kind": "iso", "tau": t} for t in TAUS]
    if not diagonal:
        for values in product((-1, 0, 1), repeat=3):
            nonzero = [v for v in values if v]
            if not nonzero or nonzero[0] < 0:
                continue
            u = np.asarray(values, dtype=float)
            u /= np.linalg.norm(u)
            result += [{"kind": "rank1", "tau": t, "u": u.tolist()} for t in TAUS]
    return result


def validate(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float)
    if y.shape != se.shape or y.ndim != 2 or y.shape[1] != 3 or not len(y):
        raise ValueError("expected nonempty N by 3 effect and SE arrays")
    if not np.isfinite(y).all() or not np.isfinite(se).all() or np.any(se <= 0):
        raise ValueError("effects and SE must be finite and SE strictly positive")
    return y, se


def component(y, se, specification):
    """Return log density and posterior mean; rank-one algebra stays positive."""
    y, se = validate(y, se)
    d = se * se
    kind, variance = specification["kind"], specification["tau"] ** 2
    if kind in ("zero", "iso"):
        total = d + variance
        logdet = np.log(total).sum(1)
        quad = np.sum(y*y/total, axis=1)
        mean = y * (variance/total)
    elif kind == "rank1":
        u = np.asarray(specification["u"], float)
        inv = 1/d
        a = np.sum(inv*u*u, axis=1)
        b = np.sum(inv*y*u, axis=1)
        denom = 1+variance*a
        # Residual plus parallel contribution avoids subtracting two large norms.
        residual = y-(b/a)[:, None]*u
        quad = np.sum(residual*residual*inv, axis=1) + (b*b/a)/denom
        logdet = np.log(d).sum(1)+np.log1p(variance*a)
        mean = (variance*b/denom)[:, None]*u
    else:
        raise ValueError(f"unknown component {kind}")
    logpdf = -.5 * (3*np.log(2*np.pi)+logdet+quad)
    if not np.isfinite(logpdf).all() or not np.isfinite(mean).all():
        raise ValueError("nonfinite component arithmetic")
    return logpdf, mean


def log_densities(y, se, specifications):
    return np.column_stack([component(y, se, c)[0] for c in specifications])


def predictive(logpdf, weights):
    weights = np.asarray(weights, float)
    if np.any(weights < 0) or not np.isclose(weights.sum(), 1., atol=1e-10):
        raise ValueError("weights outside probability simplex")
    with np.errstate(divide="ignore"):
        return logsumexp(logpdf + np.log(weights)[None, :], axis=1)


def fit(logpdf, max_iter=500, tolerance=1e-6):
    logpdf = np.asarray(logpdf, float)
    if logpdf.ndim != 2 or not len(logpdf) or not np.isfinite(logpdf).all():
        raise ValueError("invalid training densities")
    weights = np.full(logpdf.shape[1], 1/logpdf.shape[1])
    trace = [float(predictive(logpdf, weights).mean())]
    converged = False
    for iteration in range(1, max_iter+1):
        normalizer = predictive(logpdf, weights)
        with np.errstate(divide="ignore"):
            responsibility = np.exp(logpdf + np.log(weights)[None, :] - normalizer[:, None])
        weights = responsibility.mean(0)
        weights /= weights.sum()
        current = float(predictive(logpdf, weights).mean())
        change = current-trace[-1]
        if not np.isfinite(current) or change < -1e-9:
            raise ArithmeticError(f"EM likelihood decreased: {change}")
        trace.append(current)
        if iteration % 25 == 0:
            print(f"EM {iteration}: mean log density {current:.7f}, improvement {change:.3g}", flush=True)
        if iteration >= 10 and change < tolerance:
            converged = True
            break
    return {"weights": weights.tolist(), "iterations": iteration,
            "converged": converged, "mean_loglik_trace": trace}


def posterior(y, se, specifications, weights):
    values = [component(y, se, c) for c in specifications]
    lp = np.column_stack([p for p, _ in values])
    normalizer = predictive(lp, weights)
    with np.errstate(divide="ignore"):
        r = np.exp(lp+np.log(weights)[None, :]-normalizer[:, None])
    return sum(r[:, j, None]*m for j, (_, m) in enumerate(values))


def shuffled_targets(y, se, seed=20260909):
    """T by G by 3 arrays; preserve complete cases and paired marginals per gene."""
    y, se = np.asarray(y, float), np.asarray(se, float)
    if y.shape != se.shape or y.ndim != 3 or y.shape[2] != 3:
        raise ValueError("expected target by gene by condition")
    valid = np.all(np.isfinite(y) & np.isfinite(se) & (se > 0), axis=2)
    a, b = y.copy(), se.copy()
    rng = np.random.default_rng(seed)
    for g in range(y.shape[1]):
        idx = np.flatnonzero(valid[:, g])
        for c in (1, 2):
            perm = rng.permutation(idx)
            a[idx, g, c], b[idx, g, c] = y[perm, g, c], se[perm, g, c]
    return a, b
