"""Source-only empirical marginal calibration with an explicit zero anchor."""

import numpy as np

SEED = 3390909
ARMS = ("source_mean", "diagonal", "marginal", "gaussian", "clipping",
        "curves_shuffled")


def fit_curves(training):
    training = np.asarray(training, dtype=float)
    if training.ndim != 3 or training.shape[1] < 2:
        raise ValueError("donor x training-target x gene array required")
    if not np.isfinite(training).all():
        raise ValueError("nonfinite training value")
    sorted_values = np.sort(training, axis=1)
    return {"sorted": sorted_values, "barycenter": sorted_values.mean(0),
            "sd": training.std(1)}


def curve_map(values, donor_quantiles, common_quantiles):
    x, first, counts = np.unique(donor_quantiles, return_index=True,
                                 return_counts=True)
    if len(x) < 2:
        return np.asarray(values, dtype=float).copy()
    # np.unique returns each run's start because donor_quantiles is sorted.
    y = np.add.reduceat(common_quantiles, first)/counts
    return np.interp(values, x, y)-np.interp(0., x, y)


def predict(curves, queries, mode="marginal", gene_permutation=None):
    queries = np.asarray(queries, dtype=float)
    q = curves["sorted"]
    bary = curves["barycenter"]
    if queries.ndim != 3 or queries.shape[0] != q.shape[0] or queries.shape[2] != q.shape[2]:
        raise ValueError("query donor/gene axes differ")
    if not np.isfinite(queries).all():
        raise ValueError("nonfinite query value")
    if mode not in ("marginal", "gaussian", "clipping"):
        raise ValueError(mode)
    ix = np.arange(q.shape[2]) if gene_permutation is None else np.asarray(gene_permutation)
    if not np.array_equal(np.sort(ix), np.arange(q.shape[2])):
        raise ValueError("gene permutation must be complete")
    result = np.zeros(queries.shape[1:])
    common_sd = curves["sd"].mean(0)
    for donor in range(q.shape[0]):
        for gene, curve_gene in enumerate(ix):
            values = queries[donor, :, gene]
            sorted_values = q[donor, :, curve_gene]
            if mode == "marginal":
                mapped = curve_map(values, sorted_values, bary[:, curve_gene])
            elif mode == "gaussian":
                sd = curves["sd"][donor, curve_gene]
                mapped = values*(common_sd[curve_gene]/sd if sd > 1e-12 else 1.)
            elif sorted_values[0] == sorted_values[-1]:
                mapped = values
            else:
                mapped = np.clip(values, sorted_values[0], sorted_values[-1])
                mapped = mapped-np.clip(0., sorted_values[0], sorted_values[-1])
            result[:, gene] += mapped/q.shape[0]
    if not np.isfinite(result).all():
        raise ArithmeticError("nonfinite output")
    return result


def interval(delta, draws=10000):
    rng = np.random.default_rng(SEED)
    n_lines, n_targets = delta.shape
    target, crossed = [], []
    for start in range(0, draws, 100):
        n = min(100, draws-start)
        ti = rng.integers(n_targets, size=(n, n_targets))
        li = rng.integers(n_lines, size=(n, n_lines))
        means = delta[:, ti].mean(-1).T
        target.extend(means.mean(1))
        crossed.extend(means[np.arange(n)[:, None], li].mean(1))
    byline = delta.mean(1)
    return {"mean": float(byline.mean()),
            "target_ci95": np.quantile(target, [.025, .975]).tolist(),
            "crossed_ci95": np.quantile(crossed, [.025, .975]).tolist(),
            "line_means": byline.tolist(),
            "positive_lines": int((byline > 0).sum()),
            "minimum_leave_one_line_out":
                float(((byline.sum()-byline)/(n_lines-1)).min())}
