"""Sparse graph-local ridge on a frozen diagonal response baseline."""

import numpy as np

SEED, K = 3400909, 8
ARMS = ("source_mean", "diagonal", "own_refit", "graph", "random_graph",
        "correlation_graph", "labels_shuffled")


def neighbors_from_similarity(similarity, degrees=None, available=None, positive=False):
    similarity = np.asarray(similarity, dtype=float)
    n = len(similarity)
    if similarity.shape != (n, n) or not np.isfinite(similarity).all():
        raise ValueError("finite square similarity required")
    available = np.ones(n, bool) if available is None else np.asarray(available, dtype=bool)
    degrees = np.full(n, K) if degrees is None else np.asarray(degrees)
    out = np.full((n, K), -1, dtype=int)
    for gene in range(n):
        if not available[gene]:
            continue
        candidates = np.flatnonzero(available & (np.arange(n) != gene))
        if positive:
            candidates = candidates[similarity[gene, candidates] > 0]
        order = np.lexsort((candidates, -similarity[gene, candidates]))
        chosen = candidates[order[:int(degrees[gene])]]
        out[gene, :len(chosen)] = chosen
    return out


def random_neighbors(degrees, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(degrees)
    out = np.full((n, K), -1, dtype=int)
    for gene, degree in enumerate(degrees):
        choices = np.flatnonzero(np.arange(n) != gene)
        out[gene, :int(degree)] = rng.choice(choices, int(degree), replace=False)
    return out


def correlation_neighbors(source, degrees):
    centered = source-source.mean(0)
    lengths = np.linalg.norm(centered, axis=0)
    standardized = centered/np.maximum(lengths, 1e-12)
    return neighbors_from_similarity(np.abs(standardized.T@standardized), degrees=degrees)


def feature_columns(source, gene, neighbors):
    indices = np.r_[gene, neighbors[gene]]
    valid = indices >= 0
    x = np.asarray(source[:, np.maximum(indices, 0)], dtype=float).copy()
    x[:, ~valid] = 0.
    return x


def fit(source, response, gain, neighbors):
    source, response = np.asarray(source, dtype=float), np.asarray(response, dtype=float)
    if source.shape != response.shape or not np.isfinite(source).all() or not np.isfinite(response).all():
        raise ValueError("finite matched source/response matrices required")
    n, genes = source.shape
    if neighbors.shape != (genes, K) or np.any(neighbors < -1) or np.any(neighbors >= genes):
        raise ValueError("invalid neighbor shape/indices")
    beta = np.zeros((genes, K+1))
    scale = np.zeros_like(beta)
    residual = response-source*gain
    for gene in range(genes):
        x = feature_columns(source, gene, neighbors)
        scale[gene] = np.maximum(np.sqrt(np.mean(x*x, axis=0)), 1e-12)
        z = x/scale[gene]
        matrix = z.T@z+n*np.eye(K+1)
        beta[gene] = np.linalg.solve(matrix, z.T@residual[:, gene])
    return {"gain": np.asarray(gain).copy(), "neighbors": neighbors.copy(),
            "beta": beta, "scale": scale}


def predict(model, source):
    out = source*model["gain"]
    for gene in range(source.shape[1]):
        x = feature_columns(source, gene, model["neighbors"])/model["scale"][gene]
        out[:, gene] += x@model["beta"][gene]
    if not np.isfinite(out).all():
        raise ArithmeticError("nonfinite prediction")
    return out


def interval(delta, draws=10000):
    rng = np.random.default_rng(SEED)
    l, t = delta.shape
    target, crossed = [], []
    for start in range(0, draws, 100):
        n = min(100, draws-start)
        ti, li = rng.integers(t, size=(n, t)), rng.integers(l, size=(n, l))
        means = delta[:, ti].mean(-1).T
        target.extend(means.mean(1))
        crossed.extend(means[np.arange(n)[:, None], li].mean(1))
    byline = delta.mean(1)
    return {"mean": float(byline.mean()),
            "target_ci95": np.quantile(target, [.025, .975]).tolist(),
            "crossed_ci95": np.quantile(crossed, [.025, .975]).tolist(),
            "line_means": byline.tolist(), "positive_lines": int((byline > 0).sum()),
            "minimum_leave_one_line_out": float(((byline.sum()-byline)/(l-1)).min())}
