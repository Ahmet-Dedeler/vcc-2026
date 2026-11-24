"""Per-(target, gene) reliability gate for cross-line response transfer.

The transferred response is the equal-weight mean of source-line log2 fold changes.
A gate in [0.15, 1] multiplies that mean gene by gene:

    w = 0.15 + 0.85 * sigmoid(b + <beta, z(features)>)

Eight truth-free features (sources only; the held-out line's response is never an input):

    log_z2, log_n_src, log_ctrl, log_abs_mu,
    logit_sign_agree, log_cross_target_var, log_n_cells, log_source_mad

z^2 is the squared z-statistic of the equal-weight mean, using within-source
sampling variances: z^2 = mu^2 / (mean_s var_s / n_src). Cross-source discordance
is not folded into z^2; it is logit_sign_agree and log_source_mad.

Fixed formula (Testh), with n := n_src (the small integer the +4 prior matches):

    w = 0.15 + 0.85 * z^2/(z^2+8) * n/(n+4)

Coefficients are fit by Nelder-Mead on leave-one-training-line-out energy-weighted
cosine. Theta = 0 is exactly the equal-weight mean on any cosine (constant scale).
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

FEATURE_NAMES = (
    "log_z2",
    "log_n_src",
    "log_ctrl",
    "log_abs_mu",
    "logit_sign_agree",
    "log_cross_target_var",
    "log_n_cells",
    "log_source_mad",
)
N_FEATURES = len(FEATURE_NAMES)
W_FLOOR = 0.15
W_SPAN = 0.85
LOG_EPS = 1e-8
Z2_EPS = 1e-12
VAR_FLOOR = 1e-8


def fixed_weight(z2: np.ndarray, n_src: np.ndarray) -> np.ndarray:
    """Testh fixed formula. z2 and n_src broadcast to the (target, gene) grid."""
    z2 = np.maximum(np.asarray(z2, dtype=np.float64), 0.0)
    n = np.maximum(np.asarray(n_src, dtype=np.float64), 0.0)
    return W_FLOOR + W_SPAN * (z2 / (z2 + 8.0)) * (n / (n + 4.0))


def fixed_weight_ncells(z2: np.ndarray, n_cells: np.ndarray) -> np.ndarray:
    """Sensitivity: the same formula with n := mean source cell count."""
    return fixed_weight(z2, n_cells)


def gate_weight(theta: np.ndarray, Z: np.ndarray) -> np.ndarray:
    """w = 0.15 + 0.85 * sigmoid(theta[0] + Z @ theta[1:]). Z is (..., 8)."""
    theta = np.asarray(theta, dtype=np.float64)
    eta = theta[0] + np.tensordot(np.asarray(Z, dtype=np.float64), theta[1:], axes=([-1], [0]))
    eta = np.clip(eta, -30.0, 30.0)
    return W_FLOOR + W_SPAN * expit(eta)


def _logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 1e-3, 1.0 - 1e-3)
    return np.log(p / (1.0 - p))


def aggregate_sources(
    delta: np.ndarray,
    var: np.ndarray,
    n_cells: np.ndarray,
    ctrl_log: np.ndarray,
    own_gene_col: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Combine source lines into one transferred response and the raw features.

    delta, var: (n_src, n_pert, n_gene), NaN where that source lacks the entry.
    n_cells: (n_src, n_pert), NaN where the source lacks the perturbation.
    ctrl_log: (n_src, n_gene) log(control CPM + 1).
    """
    delta = np.asarray(delta, dtype=np.float64)
    var = np.asarray(var, dtype=np.float64)
    n_cells = np.asarray(n_cells, dtype=np.float64)
    ctrl_log = np.asarray(ctrl_log, dtype=np.float64)
    valid = np.isfinite(delta)
    n_src = valid.sum(axis=0).astype(np.float64)  # (P, G)

    # Equal-weight mean over sources that measured the entry.
    mu = np.divide(
        np.nansum(np.where(valid, delta, 0.0), axis=0),
        n_src,
        out=np.full(delta.shape[1:], np.nan, dtype=np.float64),
        where=n_src > 0,
    )

    # Sampling variance of the mean: mean(var_s) / n_src.
    var_ok = valid & np.isfinite(var) & (var >= 0)
    var_sum = np.nansum(np.where(var_ok, var, 0.0), axis=0)
    n_var = var_ok.sum(axis=0).astype(np.float64)
    mean_var = np.divide(var_sum, n_var, out=np.full(mu.shape, np.nan), where=n_var > 0)
    v_mu = np.divide(mean_var, n_src, out=np.full(mu.shape, np.nan), where=(n_src > 0) & np.isfinite(mean_var))
    v_mu = np.where(np.isfinite(v_mu), np.maximum(v_mu, VAR_FLOOR), np.nan)
    z2 = np.divide(mu * mu, v_mu, out=np.zeros(mu.shape, dtype=np.float64), where=np.isfinite(v_mu) & np.isfinite(mu))
    z2 = np.where(np.isfinite(z2), z2, 0.0)

    # Sign agreement with the sign of the mean. Undefined (0.5) when the mean is ~0.
    mu_sign = np.sign(mu)
    src_sign = np.sign(np.where(valid, delta, np.nan))
    agree = valid & (src_sign == mu_sign) & (src_sign != 0) & (np.abs(mu) >= 1e-8)
    n_agree = agree.sum(axis=0).astype(np.float64)
    p_agree = np.divide(n_agree, n_src, out=np.full(mu.shape, 0.5), where=(n_src > 0) & (np.abs(mu) >= 1e-8))
    p_agree = np.where(np.abs(mu) < 1e-8, 0.5, p_agree)

    # MAD across sources, around the median. NaN sources ignored.
    # nanmedian / nanvar warn on all-NaN or single-point slices; those entries are unused.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(np.where(valid, delta, np.nan), axis=0)
        absdev = np.abs(np.where(valid, delta, np.nan) - med)
        mad = np.nanmedian(absdev, axis=0)
        # Cross-target variance of the downstream response: the cis knockdown
        # (the target's own gene) is not a property of that gene's responsiveness.
        mu_var = mu
        if own_gene_col is not None:
            mu_var = mu.copy()
            cols = np.asarray(own_gene_col)
            hit = cols >= 0
            mu_var[np.arange(mu_var.shape[0])[hit], cols[hit]] = np.nan
        cross_var = np.nanvar(mu_var, axis=0, ddof=1)  # (G,) gene-level
    mad = np.where(np.isfinite(mad), mad, 0.0)
    cross_var = np.where(np.isfinite(cross_var), cross_var, 0.0)

    # Cell counts: mean over sources that have the perturbation (broadcast across genes).
    pert_measured = np.isfinite(n_cells)  # (S, P)
    n_pert_src = pert_measured.sum(axis=0).astype(np.float64)  # (P,)
    n_cells_sum = np.nansum(np.where(pert_measured, n_cells, 0.0), axis=0)
    n_cells_mean = np.divide(
        n_cells_sum, n_pert_src, out=np.ones(delta.shape[1], dtype=np.float64), where=n_pert_src > 0
    )

    # Control expression: mean of source log(CPM+1), gene-level.
    ctrl_ok = np.isfinite(ctrl_log)
    ctrl_mean = np.divide(
        np.nansum(np.where(ctrl_ok, ctrl_log, 0.0), axis=0),
        ctrl_ok.sum(axis=0).astype(np.float64),
        out=np.zeros(delta.shape[2], dtype=np.float64),
        where=ctrl_ok.any(axis=0),
    )

    P, G = mu.shape
    raw = np.empty((P, G, N_FEATURES), dtype=np.float64)
    raw[:, :, 0] = np.log(z2 + LOG_EPS)
    raw[:, :, 1] = np.log(n_src + LOG_EPS)
    raw[:, :, 2] = ctrl_mean[None, :]
    raw[:, :, 3] = np.log(np.abs(mu) + LOG_EPS)
    raw[:, :, 4] = _logit(p_agree)
    raw[:, :, 5] = np.log(cross_var + LOG_EPS)[None, :]
    raw[:, :, 6] = np.log(n_cells_mean + LOG_EPS)[:, None]
    raw[:, :, 7] = np.log(mad + LOG_EPS)
    # Entries with no source are not features.
    raw = np.where(n_src[..., None] > 0, raw, np.nan)

    return {
        "mu": mu,
        "z2": z2,
        "n_src": n_src,
        "n_cells_mean": n_cells_mean,
        "raw": raw,
    }


def zscore_features(raw: np.ndarray, stat_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Z-score each feature on stat_mask (target, gene). Non-finite -> column median, then z.

    A feature with no finite values, or with ~0 std, becomes 0.
    Returns Z (same shape, 0 outside finite raw), mean, std.
    """
    raw = np.asarray(raw, dtype=np.float64)
    stat_mask = np.asarray(stat_mask, dtype=bool)
    P, G, F = raw.shape
    Z = np.zeros((P, G, F), dtype=np.float64)
    means = np.zeros(F, dtype=np.float64)
    stds = np.ones(F, dtype=np.float64)
    for j in range(F):
        col = raw[:, :, j]
        use = stat_mask & np.isfinite(col)
        if not np.any(use):
            continue
        vals = col[use]
        med = float(np.median(vals))
        filled = np.where(np.isfinite(col), col, med)
        mu = float(vals.mean())
        sd = float(vals.std(ddof=0))
        means[j] = mu
        stds[j] = sd
        if sd < 1e-8:
            Z[:, :, j] = 0.0
        else:
            Z[:, :, j] = (filled - mu) / sd
        Z[:, :, j] = np.where(np.isfinite(raw[:, :, j]) | stat_mask, Z[:, :, j], 0.0)
    return Z, means, stds


def cosine_rows(pred: np.ndarray, truth: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-target cosine on mask (True = include). Returns cos, energy=||truth||^2, ok."""
    pred = np.asarray(pred, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(pred) & np.isfinite(truth)
    p = np.where(mask, pred, 0.0)
    t = np.where(mask, truth, 0.0)
    pn = np.sqrt(np.sum(p * p, axis=1))
    tn = np.sqrt(np.sum(t * t, axis=1))
    dot = np.sum(p * t, axis=1)
    ok = (pn > 1e-12) & (tn > 1e-12) & mask.any(axis=1)
    cos = np.zeros(pred.shape[0], dtype=np.float64)
    cos[ok] = dot[ok] / (pn[ok] * tn[ok])
    energy = tn * tn
    return cos, energy, ok


def energy_weighted_cosine(cos: np.ndarray, energy: np.ndarray, ok: np.ndarray) -> float:
    if not np.any(ok):
        return float("nan")
    e = energy[ok]
    denom = float(e.sum())
    if denom <= 0:
        return float("nan")
    return float(np.sum(e * cos[ok]) / denom)


def unweighted_cosine(cos: np.ndarray, ok: np.ndarray) -> float:
    if not np.any(ok):
        return float("nan")
    return float(cos[ok].mean())


def normalised_rank(pred: np.ndarray, truth: np.ndarray, keep_genes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For each target, 1 - (n strictly closer)/(n-1), by cosine, on keep_genes.

    keep_genes is a 1-d bool over genes (panel genes already removed).
    Also returns top1 (nothing strictly closer).
    """
    pred = np.asarray(pred, dtype=np.float64)[:, keep_genes]
    truth = np.asarray(truth, dtype=np.float64)[:, keep_genes]
    pred = np.where(np.isfinite(pred), pred, 0.0)
    truth = np.where(np.isfinite(truth), truth, 0.0)
    pn = np.linalg.norm(pred, axis=1)
    tn = np.linalg.norm(truth, axis=1)
    # Guard empty rows.
    pn_safe = np.maximum(pn, 1e-12)
    tn_safe = np.maximum(tn, 1e-12)
    sim = (pred @ truth.T) / (pn_safe[:, None] * tn_safe[None, :])
    n = sim.shape[0]
    own = np.diag(sim).copy()
    n_closer = (sim > own[:, None]).sum(axis=1).astype(np.float64)
    # A zero-norm row is not a real match.
    bad = (pn <= 1e-12) | (tn <= 1e-12)
    if n <= 1:
        rank = np.ones(n, dtype=np.float64)
    else:
        rank = 1.0 - n_closer / (n - 1)
    rank = np.where(bad, np.nan, rank)
    top1 = (n_closer == 0) & ~bad
    return rank, top1


def prepare_task(
    delta: np.ndarray,
    var: np.ndarray,
    n_cells: np.ndarray,
    ctrl_log: np.ndarray,
    truth: np.ndarray,
    own_gene_col: np.ndarray,
    target_ok: np.ndarray,
) -> dict:
    """Build mu, z-scored features and the cosine mask for one prediction task.

    own_gene_col: (n_pert,) int, -1 if the target gene is not on the axis.
    target_ok: (n_pert,) bool, targets eligible for this task (present in the
    held-out line and in at least two sources — the caller enforces that).
    """
    agg = aggregate_sources(delta, var, n_cells, ctrl_log, own_gene_col=own_gene_col)
    P, G = agg["mu"].shape
    own = np.zeros((P, G), dtype=bool)
    cols = np.asarray(own_gene_col)
    rows = np.arange(P)
    hit = cols >= 0
    own[rows[hit], cols[hit]] = True
    finite_truth = np.isfinite(truth)
    # Cosine uses every gene the held-out line measured, except the target's own gene.
    cos_mask = finite_truth & ~own & np.isfinite(agg["mu"]) & target_ok[:, None]
    stat_mask = cos_mask & np.isfinite(agg["raw"]).all(axis=-1)
    # If a feature is non-finite on a gene we still score, zscore imputes it.
    Z, means, stds = zscore_features(agg["raw"], stat_mask)
    return {
        "mu": agg["mu"],
        "z2": agg["z2"],
        "n_src": agg["n_src"],
        "n_cells_mean": agg["n_cells_mean"],
        "Z": Z.astype(np.float32),
        "truth": np.asarray(truth, dtype=np.float64),
        "cos_mask": cos_mask,
        "target_ok": np.asarray(target_ok, dtype=bool) & cos_mask.any(axis=1),
        "z_mean": means,
        "z_std": stds,
        "own": own,
    }


def score_weight(task: dict, w: np.ndarray) -> dict[str, np.ndarray | float]:
    pred = np.asarray(w, dtype=np.float64) * task["mu"]
    # Only score eligible targets. Zero the others so they cannot enter a later mean.
    use = task["target_ok"]
    cos, energy, ok = cosine_rows(pred, task["truth"], task["cos_mask"])
    ok = ok & use
    return {
        "pred": pred,
        "cos": cos,
        "energy": energy,
        "ok": ok,
        "ew": energy_weighted_cosine(cos, energy, ok),
        "uw": unweighted_cosine(cos, ok),
    }


def _fast_objective(tasks: list[dict]):
    """Energy-weighted cosine of w(theta) * mu, averaged over tasks.

    Masked entries are zeroed once. The linear map stays in float32 so a
    Nelder-Mead step on the full gene axis is a single matmul.
    """
    packs = []
    for task in tasks:
        mask = task["cos_mask"] & task["target_ok"][:, None]
        # Z is not zeroed: masked genes are removed by zeroing mu, which keeps
        # the cosine identical and avoids a second copy of the feature tensor.
        Z = np.ascontiguousarray(task["Z"], dtype=np.float32)
        mu = np.ascontiguousarray(np.where(mask, task["mu"], 0), dtype=np.float32)
        truth = np.ascontiguousarray(np.where(mask, task["truth"], 0), dtype=np.float32)
        tn = np.sqrt(np.sum(np.square(truth), axis=1))
        packs.append(
            {
                "Z": Z.reshape(-1, N_FEATURES),
                "shape": Z.shape[:2],
                "mu": mu,
                "truth": truth,
                "tn": tn,
                "energy": tn * tn,
                "row_ok": (mask.any(axis=1) & task["target_ok"] & (tn > 1e-6)),
            }
        )

    def objective(theta: np.ndarray) -> float:
        beta = np.asarray(theta[1:], dtype=np.float32)
        b0 = np.float32(theta[0])
        scores = []
        for pack in packs:
            eta = np.clip(b0 + pack["Z"] @ beta, -30.0, 30.0).reshape(pack["shape"])
            w = (W_FLOOR + W_SPAN * expit(eta)).astype(np.float32, copy=False)
            pred = w * pack["mu"]
            pn = np.sqrt(np.sum(np.square(pred), axis=1))
            dot = np.sum(pred * pack["truth"], axis=1)
            ok = pack["row_ok"] & (pn > 1e-6)
            if not np.any(ok):
                continue
            cos = dot[ok] / (pn[ok] * pack["tn"][ok])
            scores.append(float(np.sum(pack["energy"][ok] * cos) / np.sum(pack["energy"][ok])))
        if not scores:
            return 0.0
        return -float(np.mean(scores))

    return objective


def fit_theta(tasks: list[dict], seed: int = 486, n_restarts: int = 4, maxiter: int = 800, maxfev: int | None = None) -> dict:
    """Nelder-Mead. Objective = mean over tasks of energy-weighted cosine. Maximise.

    Restarts: zeros (the equal-weight point) plus n_restarts N(0, 0.5) draws.
    """
    if not tasks:
        raise ValueError("no training tasks")
    objective = _fast_objective(tasks)
    if maxfev is None:
        maxfev = maxiter * 6

    rng = np.random.default_rng(seed)
    starts = [np.zeros(N_FEATURES + 1, dtype=np.float64)]
    for _ in range(n_restarts):
        starts.append(rng.normal(0.0, 0.5, size=N_FEATURES + 1))
    # SciPy's default simplex step is 0.00025 when x0 is 0, which never leaves
    # the flat region of the sigmoid. Use a step of 0.5 on every start.
    best = None
    traces = []
    n_dim = N_FEATURES + 1
    for x0 in starts:
        simplex = np.tile(x0, (n_dim + 1, 1))
        simplex[1:] += np.eye(n_dim) * 0.5
        options = {
            "maxiter": maxiter,
            "maxfev": int(maxfev),
            "xatol": 1e-3,
            "fatol": 1e-5,
            "adaptive": True,
            "initial_simplex": simplex,
        }
        res = minimize(objective, x0, method="Nelder-Mead", options=options)
        traces.append({"fun": float(res.fun), "nfev": int(res.nfev), "success": bool(res.success)})
        if best is None or res.fun < best.fun:
            best = res
    theta = np.asarray(best.x, dtype=np.float64)
    return {
        "theta": theta,
        "train_ew": float(-best.fun),
        "nfev": int(sum(t["nfev"] for t in traces)),
        "restarts": traces,
        "intercept": float(theta[0]),
        "coef": {FEATURE_NAMES[j]: float(theta[j + 1]) for j in range(N_FEATURES)},
    }


def paired_bootstrap(
    metrics: dict[str, dict[str, np.ndarray]],
    B: int = 4000,
    seed: int = 486,
) -> dict:
    """Target bootstrap of unweighted cosine, energy-weighted cosine, and mean normalised rank.

    metrics[name] = {cos, energy, rank} aligned on the same targets (ineligible already dropped).
    Ranks are the full-set normalised ranks; the bootstrap resamples those per-target scores
    (paired across methods). Energy weights are recomputed inside each draw.
    """
    names = list(metrics)
    n = len(metrics[names[0]]["cos"])
    rng = np.random.default_rng(seed)
    # Chunk so the index tensor stays small.
    chunk = 500
    acc = {name: {"uw": [], "ew": [], "rank": []} for name in names}
    done = 0
    while done < B:
        b = min(chunk, B - done)
        idx = rng.integers(0, n, size=(b, n))
        for name in names:
            cos = metrics[name]["cos"][idx]
            energy = np.maximum(metrics[name]["energy"][idx], 0.0)
            rank = metrics[name]["rank"][idx]
            acc[name]["uw"].append(cos.mean(axis=1))
            denom = energy.sum(axis=1)
            acc[name]["ew"].append(np.where(denom > 0, (energy * cos).sum(axis=1) / np.maximum(denom, 1e-12), np.nan))
            acc[name]["rank"].append(np.nanmean(rank, axis=1))
            if "top1" in metrics[name]:
                acc[name].setdefault("top1", [])
                acc[name]["top1"].append(metrics[name]["top1"][idx].mean(axis=1))
        done += b
    out = {}
    for name in names:
        out[name] = {k: np.concatenate(v) for k, v in acc[name].items()}
    return out


def bootstrap_summary(point: dict[str, dict[str, float]], draws: dict, contrasts: list[tuple[str, str]]) -> dict:
    """point[name][metric] is the full-sample estimate. CI is the percentile interval of the draws.

    contrasts are (a, b) pairs; the reported gain is a - b.
    """
    summary = {"methods": {}, "gains": {}}
    for name, series in draws.items():
        summary["methods"][name] = {}
        for metric, samples in series.items():
            samples = samples[np.isfinite(samples)]
            summary["methods"][name][metric] = {
                "point": float(point[name][metric]),
                "ci95": [float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5))],
                "boot_mean": float(samples.mean()),
                "boot_sd": float(samples.std(ddof=1)),
            }
    for a, b in contrasts:
        summary["gains"][f"{a}_minus_{b}"] = {}
        for metric in ("uw", "ew", "rank", "top1"):
            if metric not in draws[a] or metric not in draws[b]:
                continue
            diff = draws[a][metric] - draws[b][metric]
            diff = diff[np.isfinite(diff)]
            summary["gains"][f"{a}_minus_{b}"][metric] = {
                "point": float(point[a][metric] - point[b][metric]),
                "ci95": [float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))],
                "boot_sd": float(diff.std(ddof=1)),
                "frac_positive": float(np.mean(diff > 0)),
            }
    return summary
