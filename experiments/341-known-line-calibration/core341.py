"""Matched full-coordinate estimators for a different label regime."""

import numpy as np

SEED = 3410909
ARMS = ("source_mean", "crossline_diag", "crossline_full", "same_line_diag",
        "same_line_full", "labels_shuffled_full")


def diagonal(x, y):
    num, den = np.sum(x*y, axis=0), np.sum(x*x, axis=0)
    prior = np.median(den)
    return (num+prior)/np.maximum(den+prior, 1e-12)


def fit(x, y, gain=None):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 2 or x.shape != y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("finite matched matrices required")
    gain = diagonal(x, y) if gain is None else np.asarray(gain, dtype=float)
    scale = np.maximum(np.sqrt(np.mean(x*x, axis=0)), 1e-12)
    z = x/scale
    beta = np.linalg.solve(z.T@z+len(x)*np.eye(x.shape[1]), z.T@(y-x*gain))
    return {"gain": gain.copy(), "scale": scale, "beta": beta}


def predict(model, x):
    return x*model["gain"]+(x/model["scale"])@model["beta"]


def sign_endpoints(prediction, truth):
    weight = np.abs(truth)
    denominator = weight.sum(1)
    if np.any(denominator <= 0):
        raise ValueError("every target must have positive truth mass")
    correct = np.sign(prediction) == np.sign(truth)
    full = (weight*correct).sum(1)/denominator
    # Stable sorting preserves coordinate index for exact absolute-value ties.
    top = np.argsort(-weight, axis=1, kind="stable")[:, :min(64, truth.shape[1])]
    w = np.take_along_axis(weight, top, axis=1)
    c = np.take_along_axis(correct, top, axis=1)
    return {"weighted_sign": full, "top64_weighted_sign": (w*c).sum(1)/w.sum(1)}


def interval(delta, draws=10000):
    rng = np.random.default_rng(SEED)
    l, t = delta.shape
    target, crossed = [], []
    for start in range(0, draws, 100):
        n = min(100, draws-start)
        ti, li = rng.integers(t, size=(n, t)), rng.integers(l, size=(n, l))
        byline = delta[:, ti].mean(-1).T
        target.extend(byline.mean(1))
        crossed.extend(byline[np.arange(n)[:, None], li].mean(1))
    means = delta.mean(1)
    return {"mean": float(means.mean()),
            "target_ci95": np.quantile(target, [.025, .975]).tolist(),
            "crossed_ci95": np.quantile(crossed, [.025, .975]).tolist(),
            "line_means": means.tolist(), "positive_lines": int((means > 0).sum())}
