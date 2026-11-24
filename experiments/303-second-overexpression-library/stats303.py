"""C7 readout v2: pure, strict array statistics; no files or reserve access.

Primary: centre both matrices on the same metadata-fixed TF cohort, then drop
the UNION of all cohort own-target columns from BOTH sides. Matrix entry [i,j]
is cos(KD_i, OE_j); a null draw indexes [i, permutation[i]], keeping KD and its
features fixed. Uniform permutations include fixed points. Pair-specific own
exclusion and derangements exist ONLY as explicitly labelled sensitivities.

NaN is permitted in *features* to mark frozen missingness, not in expression.
Zero-norm responses produce undefined cosines, never artificial zero signals.
Undefined statistics/null draws make a cell unestimable (p=1); the family
always contains four cells. Bootstrap intervals condition on frozen centred
rows: they are NOT between-study intervals and do not re-estimate centring.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata

VERSION = "C7-fixed-axis-uniform-v2"
SEED = 20260907
N_NULL = 2_000
N_BOOT = 20_000
CELLS = ("C7a", "C7b", "C7c", "C7d")
CI_SCOPE = ("95% percentile TF bootstrap, conditional on frozen centred rows "
            "and each cell's fixed feature-complete cohort; not between-study CI")


def numeric(a, name: str, ndim: int, *, missing: bool = False) -> np.ndarray:
    """Validate before conversion (no bool/object/string/complex coercion)."""
    a = np.asarray(a)
    if a.ndim != ndim or a.dtype.kind not in "fiu":
        raise ValueError(f"{name}: expected {ndim}-D real numeric array")
    if np.isinf(a).any() or (not missing and not np.isfinite(a).all()):
        raise ValueError(f"{name}: nonfinite values")
    return np.asarray(a, dtype=np.float64)


def index_vector(a, n: int, *, name: str = "indices",
                 upper: int | None = None, allow_absent: bool = False) -> np.ndarray:
    a = np.asarray(a)
    if a.shape != (n,) or a.dtype.kind not in "iu":
        raise ValueError(f"{name}: expected integer vector of length {n}")
    if (a < (-1 if allow_absent else 0)).any():
        raise ValueError(f"{name}: negative index")
    if upper is not None and (a >= upper).any():
        raise ValueError(f"{name}: index out of bounds")
    return a.astype(np.int64)


def residualise(M: np.ndarray, centre: np.ndarray | None) -> np.ndarray:
    M = numeric(M, "expression", 2)
    if centre is None:
        return M.copy()
    c = numeric(centre, "centre", 1)
    if c.shape != (M.shape[1],):
        raise ValueError("centre shape does not match gene axis")
    return M - c[None, :]


def cohort_residuals(OE, KD) -> tuple[np.ndarray, np.ndarray]:
    OE, KD = numeric(OE, "OE", 2), numeric(KD, "KD", 2)
    if OE.shape != KD.shape or min(OE.shape) < 1:
        raise ValueError("OE/KD must share a nonempty TF-by-gene shape")
    return OE - OE.mean(0), KD - KD.mean(0)


def fixed_axis_mask(own, n: int, n_gene: int) -> np.ndarray:
    own = index_vector(own, n, name="own", upper=n_gene, allow_absent=True)
    keep = np.ones(n_gene, dtype=bool)
    keep[own[own >= 0]] = False
    if not keep.any():
        raise ValueError("all measured coordinates were excluded")
    return keep


def _divide_cosine(dot, k2, o2):
    denom = np.sqrt(np.maximum(k2, 0)) * np.sqrt(np.maximum(o2, 0))
    out = np.full(np.broadcast_shapes(dot.shape, denom.shape), np.nan)
    np.divide(dot, denom, out=out, where=denom > 0)
    return np.clip(out, -1.0, 1.0)


def pair_cosine(a: np.ndarray, b: np.ndarray, drop: np.ndarray) -> float:
    """Literal reference implementation, also used in synthetic tests."""
    a, b = numeric(a, "a", 1), numeric(b, "b", 1)
    if a.shape != b.shape:
        raise ValueError("pair_cosine: unequal axes")
    drop = index_vector(drop, np.asarray(drop).size, upper=a.size)
    keep = np.ones(a.size, bool)
    keep[drop] = False
    a, b = a[keep], b[keep]
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.clip(a @ b / denom, -1, 1)) if denom > 0 else float("nan")


def cosine_matrix(OEr, KDr, own, *, exclusion: str = "cohort") -> np.ndarray:
    """One O(n^2 G) product, never repeated inside the permutation loop.

    exclusion='cohort' is the fixed-axis primary. 'pairwise' reproduces the
    own-two-coordinate draft as sensitivity (matched/shuffled axes differ).
    'none' is useful when the caller has already fixed a restricted axis.
    """
    O, K = numeric(OEr, "OE", 2), numeric(KDr, "KD", 2)
    if O.shape != K.shape or min(O.shape) < 1:
        raise ValueError("OE/KD shape mismatch or empty matrix")
    n, g = K.shape
    own = index_vector(own, n, name="own", upper=g, allow_absent=True)
    if exclusion == "cohort":
        keep = fixed_axis_mask(own, n, g)
        O, K = O[:, keep], K[:, keep]
    elif exclusion not in ("none", "pairwise"):
        raise ValueError("unknown exclusion convention")
    dot = K @ O.T
    k2, o2 = np.sum(K * K, axis=1)[:, None], np.sum(O * O, axis=1)[None, :]
    if exclusion == "pairwise":
        valid = own >= 0
        safe = np.maximum(own, 0)
        Ko, Oo = K[:, safe], O[:, safe]
        ki, oj = np.diag(Ko) * valid, np.diag(Oo) * valid
        second = valid[None, :] & (own[:, None] != own[None, :])
        dot = dot - ki[:, None] * Oo.T - Ko * oj[None, :] * second
        k2 = k2 - ki[:, None] ** 2 - Ko ** 2 * second
        o2 = o2 - Oo.T ** 2 * valid[:, None] - oj[None, :] ** 2 * second
    return _divide_cosine(dot, k2, o2)


def validate_permutations(perms, n: int) -> np.ndarray:
    p = np.asarray(perms)
    if p.ndim != 2 or p.shape[1] != n or p.shape[0] < 1 or p.dtype.kind not in "iu":
        raise ValueError("permutations: expected nonempty integer draws-by-TF array")
    if not np.all(np.sort(p, axis=1) == np.arange(n)):
        raise ValueError("pairing must be a bijection of the whole TF cohort")
    return p.astype(np.int64)


def uniform_permutations(n: int, draws: int, rng: np.random.Generator) -> np.ndarray:
    if not isinstance(n, (int, np.integer)) or n < 2 or draws < 1:
        raise ValueError("uniform permutation requires n>=2 and draws>=1")
    return np.stack([rng.permutation(n) for _ in range(draws)])


def derangements(n: int, draws: int, rng: np.random.Generator) -> np.ndarray:
    """Descriptive sensitivity only; not the full permutation group."""
    if n < 2 or draws < 1:
        raise ValueError("derangements require n>=2 and draws>=1")
    out = np.empty((draws, n), np.int64)
    for d in range(draws):
        while True:
            p = rng.permutation(n)
            if not np.any(p == np.arange(n)):
                out[d] = p
                break
    return out


def s_vector(OEr, KDr, own, pairing=None, *, exclusion="cohort") -> np.ndarray:
    """Matched/shuffled cosines on the SAME fixed primary coordinate axis."""
    C = cosine_matrix(OEr, KDr, own, exclusion=exclusion)
    p = np.arange(C.shape[0]) if pairing is None else pairing
    p = validate_permutations(np.asarray(p)[None, :], C.shape[0])[0]
    return C[np.arange(C.shape[0]), p]


def cell_stats(S: np.ndarray, thr: float = 0.05) -> dict:
    """Descriptive finite subset, explicitly counted; NOT an inferential gate."""
    S = numeric(S, "cosines", 1, missing=True)
    s = S[np.isfinite(S)]
    nan = float("nan")
    return dict(n=int(s.size), n_total=int(S.size), n_undefined=int(S.size - s.size),
                median=float(np.median(s)) if s.size else nan,
                mean=float(np.mean(s)) if s.size else nan,
                share_neg=float(np.mean(s < 0)) if s.size else nan,
                share_below=float(np.mean(s < -thr)) if s.size else nan,
                share_above=float(np.mean(s > thr)) if s.size else nan,
                n_below=int(np.sum(s < -thr)), n_above=int(np.sum(s > thr)),
                tail_asym=float(np.mean(s < -thr) - np.mean(s > thr)) if s.size else nan)


def _tail(S, thr=0.05):
    out = np.mean(S < -thr, axis=-1) - np.mean(S > thr, axis=-1)
    return np.where(np.isfinite(S).all(axis=-1), out, np.nan)


def null_distribution(OEr, KDr, own, draws=N_NULL, seed=SEED, thr=0.05,
                      *, exclusion="cohort"):
    C = cosine_matrix(OEr, KDr, own, exclusion=exclusion)
    perms = uniform_permutations(C.shape[0], draws, np.random.default_rng(seed))
    Sall = C[np.arange(C.shape[0]), perms]
    return dict(perms=perms, median=np.median(Sall, axis=1),
                tail_asym=_tail(Sall, thr), S=Sall)


def rank_correlation(x, y) -> np.ndarray:
    """Batched Spearman, average ranks for ties; undefined rows remain NaN."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    if x.ndim < 1:
        raise ValueError("correlation needs a target axis")
    if np.isinf(x).any() or np.isinf(y).any():
        raise ValueError("infinite feature or cosine")
    if x.shape[-1] < 3:
        return np.full(x.shape[:-1], np.nan)
    a, b = rankdata(x, axis=-1), rankdata(y, axis=-1)
    a, b = a - a.mean(-1, keepdims=True), b - b.mean(-1, keepdims=True)
    denom = np.sqrt(np.sum(a*a, axis=-1) * np.sum(b*b, axis=-1))
    out = np.full(x.shape[:-1], np.nan)
    np.divide(np.sum(a*b, axis=-1), denom, out=out, where=denom > 0)
    return np.clip(out, -1, 1)


def summary(values) -> dict:
    v = np.asarray(values, float).ravel()
    f = v[np.isfinite(v)]
    out = dict(n_total=int(v.size), n_finite=int(f.size))
    out.update({k: float("nan") for k in ("mean", "sd", "min", "q025", "median", "q975", "max")})
    if f.size:
        out.update(mean=float(f.mean()), sd=float(f.std(ddof=0)),
                   **dict(zip(("min", "q025", "median", "q975", "max"),
                              map(float, np.quantile(f, [0, .025, .5, .975, 1])))))
    return out


def permutation_p(observed: float, null, side: str) -> float:
    """Inclusive +1 Monte Carlo p; never discard undefined permutation draws."""
    v = numeric(null, "null", 1, missing=True)
    if side not in ("less", "greater"):
        raise ValueError("side must be less or greater")
    if not np.isfinite(observed) or not v.size or not np.isfinite(v).all():
        return 1.0
    count = np.sum(v <= observed) if side == "less" else np.sum(v >= observed)
    return float((count + 1) / (v.size + 1))


def holm_four(pvals: dict[str, float], alpha=0.05) -> dict:
    if set(pvals) != set(CELLS):
        raise ValueError("C7 Holm family must contain exactly C7a/b/c/d")
    vals = np.asarray([pvals[k] for k in CELLS], float)
    if not np.isfinite(vals).all() or np.any((vals < 0) | (vals > 1)):
        raise ValueError("p-values must be finite probabilities; use p=1 for unestimable")
    running, out = 0.0, {}
    for i in np.argsort(vals, kind="stable"):
        position = len(out)
        running = min(1.0, max(running, (4-position)*vals[i]))
        out[CELLS[i]] = dict(p_raw=float(vals[i]), p_holm=float(running),
                            reject=bool(running <= alpha))
    return {k: out[k] for k in CELLS}


def _bootstrap(values, feature=None, *, draws=N_BOOT, seed=SEED, batch=256):
    """Rerank each paired TF resample, including ties induced by resampling."""
    y = np.asarray(values, float)
    x = None if feature is None else np.asarray(feature, float)
    n = y.size
    names = ("median", "tail_asym") if x is None else ("rho",)
    result = {k: np.full(draws, np.nan) for k in names}
    if n < (1 if x is None else 3) or not np.isfinite(y).all():
        return result
    if x is not None and (not np.isfinite(x).all() or np.unique(x).size < 2):
        return result
    rng = np.random.default_rng(seed)
    for lo in range(0, draws, batch):
        hi = min(draws, lo+batch)
        idx = rng.integers(0, n, size=(hi-lo, n))
        if x is None:
            result["median"][lo:hi] = np.median(y[idx], axis=1)
            result["tail_asym"][lo:hi] = _tail(y[idx])
        else:
            result["rho"][lo:hi] = rank_correlation(x[idx], y[idx])
    return result


def _ci(values):
    v = summary(values)
    return dict(lo=v["q025"], hi=v["q975"], n_finite=v["n_finite"],
                n_draws=v["n_total"], scope=CI_SCOPE)


def evaluate_from_null(S, null_S, conservation, geometry, *, boot_draws=N_BOOT,
                       seed=SEED) -> dict:
    """Four confirmatory cells with absolute-direction AND relative-null gates."""
    S = numeric(S, "observed S", 1, missing=True)
    N = numeric(null_S, "null S", 2, missing=True)
    if N.shape[1] != S.size or N.shape[0] < 1 or boot_draws < 1:
        raise ValueError("null/bootstrap shape or count mismatch")
    for a in (S, N):
        if np.any(np.abs(a[np.isfinite(a)]) > 1 + 1e-12):
            raise ValueError("cosines outside [-1,1]")
    features = []
    for name, x in (("C_excl", conservation), ("gen_pc20", geometry)):
        x = numeric(x, name, 1, missing=True)
        if x.shape != S.shape:
            raise ValueError("feature/TF alignment mismatch")
        features.append(x)
    null = dict(median=np.median(N, axis=1), tail_asym=_tail(N))
    median, tail = float(np.median(S)), float(_tail(S))
    a_ok = bool(np.isfinite(S).all() and np.isfinite(N).all())
    p_med = permutation_p(median, null["median"], "less")
    p_tail = permutation_p(tail, null["tail_asym"], "greater")
    boots = _bootstrap(S, draws=boot_draws, seed=seed)
    cells = {"C7a": dict(
        n=int(S.size), estimable=a_ok,
        reason=None if a_ok else "undefined response cosine in observed or null cohort",
        observed=dict(median=median, tail_asym=tail),
        p_median_less=p_med, p_tail_greater=p_tail, p_raw=max(p_med, p_tail),
        direction_ok=bool(a_ok and median < 0 and tail > 0),
        ci={k: _ci(v) for k, v in boots.items()},
        null={k: summary(v) for k, v in null.items()})}
    boot_out = {f"C7a_{k}": v for k, v in boots.items()}
    for key, feat, y, yn in (("C7b", features[0], -S, -N),
                             ("C7c", features[1], -S, -N),
                             ("C7d", features[1], np.abs(S), np.abs(N))):
        keep = np.isfinite(feat)       # fixed by USED data, never by OE outcomes
        x, yy, nn = feat[keep], y[keep], yn[:, keep]
        rho = float(rank_correlation(x, yy))
        rnull = rank_correlation(x, nn)
        ok = bool(np.isfinite(rho) and np.isfinite(rnull).all())
        b = _bootstrap(yy, x, draws=boot_draws, seed=seed)["rho"]
        null[key] = rnull
        boot_out[key] = b
        cells[key] = dict(n=int(keep.sum()), feature_missing=int((~keep).sum()),
                          estimable=ok,
                          reason=None if ok else "missing/constant feature or undefined observed/null ranks",
                          observed=rho, p_raw=permutation_p(rho, rnull, "greater"),
                          direction_ok=bool(ok and rho > 0), ci=_ci(b), null=summary(rnull))
    adjusted = holm_four({k: cells[k]["p_raw"] for k in CELLS})
    for k in CELLS:
        cells[k].update(adjusted[k])
        cells[k]["verdict"] = ("UNESTIMABLE" if not cells[k]["estimable"] else
                               "PASS" if cells[k]["direction_ok"] and cells[k]["reject"] else "FAIL")
    return dict(cells=cells, null_stats=null, bootstrap=boot_out,
                observed_summary=cell_stats(S), S=S, null_S=N)


def evaluate_c7(C, conservation, geometry, *, draws=N_NULL, boot_draws=N_BOOT,
                seed=SEED, permutations=None) -> dict:
    C = numeric(C, "cosine matrix", 2, missing=True)
    if C.shape[0] != C.shape[1] or C.shape[0] < 2:
        raise ValueError("cosine matrix must be square with at least two TFs")
    n = C.shape[0]
    p = (uniform_permutations(n, draws, np.random.default_rng(seed))
         if permutations is None else validate_permutations(permutations, n))
    result = evaluate_from_null(np.diag(C), C[np.arange(n), p], conservation, geometry,
                                boot_draws=boot_draws, seed=seed)
    result.update(permutations=p, seed=seed, n_null=int(p.shape[0]), n_boot=boot_draws,
                  permutation_scope="uniform whole-cohort OE matching; fixed points allowed; KD/features fixed",
                  fixed_points=int(np.sum(p == np.arange(n))), ci_scope=CI_SCOPE)
    return result


def partial_spearman(x, y, controls) -> float:
    """Descriptive only: rank residuals with both row magnitudes as controls."""
    x, y = numeric(x, "x", 1, missing=True), numeric(y, "y", 1, missing=True)
    Z = numeric(controls, "controls", 2, missing=True)
    if x.shape != y.shape or Z.shape[0] != x.size:
        raise ValueError("partial correlation alignment mismatch")
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(1)
    if keep.sum() < Z.shape[1] + 3:
        return float("nan")
    X = np.column_stack((np.ones(keep.sum()), rankdata(Z[keep], axis=0)))
    a, b = rankdata(x[keep]), rankdata(y[keep])
    a -= X @ np.linalg.lstsq(X, a, rcond=None)[0]
    b -= X @ np.linalg.lstsq(X, b, rcond=None)[0]
    d = np.linalg.norm(a)*np.linalg.norm(b)
    return float(a @ b / d) if d > 1e-12 else float("nan")
