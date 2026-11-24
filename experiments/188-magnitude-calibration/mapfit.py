"""188 — the L1-optimal per-gene magnitude map, fitted leave-one-line-out.  Shared definition.

THE ESTIMAND, DERIVED RATHER THAN CHOSEN
----------------------------------------
169 §5's identity (gated against `debench` to 0.0e+00) says `nmae` is an L1 loss on log2 fold
changes over the reference's DE genes.  Write our emitted value on a gene as `x*s` with
`s = sign(p)` fixed by the direction model and `x >= 0` the magnitude we are choosing, and let
`r` be the truth.  Then

    |x*s - r| = |x - r*s|          because |s| = 1

so the loss is `E|x - y|` with **`y = r * sign(p)`**, the truth PROJECTED ONTO OUR OWN SIGN.
The minimiser of `E|x - y|` over `x >= 0` is

    m(p) = max(0, median(y | p))                                                       (M)

**not** the conditional median of `|r|`, which the mandate's phrasing suggests and which is a
different and strictly worse object: `median(|r|)` is positive for every gene, so it hands a
large magnitude to genes whose sign we get right 50% of the time, and those are exactly the
genes 02:15's landscape charges at slope +1.  (M) encodes sign accuracy automatically -- where
`P(sign correct) <= 1/2` the median of `y` is <= 0 and (M) says EMIT NOTHING.  Both are fitted
below and the difference between them is reported, because it is the whole reason a power law
is the wrong family: `|p|^q` is positive and expansive everywhere, including on the ~46% of the
DE pool whose sign we get wrong.

ESTIMATOR
---------
Low capacity, no parametric family: quantile-binned medians of `y` against `|p|`, then PAVA to
enforce monotone non-decreasing, then linear interpolation.  Monotonicity matters beyond
regularisation -- it preserves the `|effect|` ordering, so the shipped `keep_top 8000` block
and every `|effect|`-keyed policy are untouched, exactly as 169 §1 verified for the power law.

LEAKAGE
-------
The map is fitted on three PLOLO folds and read on the fourth, always.  A fold is one held
donor LIBRARY, so a line never contributes to the map that scores it.
"""
from __future__ import annotations

import numpy as np

# 60 quantile bins over ~2.5M (target, gene) pairs per fold: ~40k pairs per bin, so a bin
# median has a standard error of ~1.25*sd/sqrt(40000) ~ 0.006 truth units against a truth
# median |r| of ~0.36.  Fine enough to see shape, coarse enough not to fit noise.
NBINS = 60


def pava(y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Pool-adjacent-violators: the L2 projection onto monotone non-decreasing sequences."""
    y = np.asarray(y, dtype=np.float64).copy()
    w = np.asarray(w, dtype=np.float64).copy()
    n = y.size
    val = np.empty(n)
    wt = np.empty(n)
    idx = np.empty(n, dtype=np.int64)
    k = 0
    for i in range(n):
        val[k], wt[k], idx[k] = y[i], w[i], 1
        k += 1
        while k > 1 and val[k - 1] < val[k - 2]:
            tw = wt[k - 2] + wt[k - 1]
            val[k - 2] = (wt[k - 2] * val[k - 2] + wt[k - 1] * val[k - 1]) / tw
            wt[k - 2] = tw
            idx[k - 2] += idx[k - 1]
            k -= 1
    out = np.empty(n)
    j = 0
    for i in range(k):
        out[j:j + int(idx[i])] = val[i]
        j += int(idx[i])
    return out


class Map:
    """A fitted magnitude map `x -> m(x)`, monotone non-decreasing, `m(0) >= 0`."""

    def __init__(self, edges: np.ndarray, xs: np.ndarray, ms: np.ndarray, meta: dict):
        self.edges, self.xs, self.ms, self.meta = edges, xs, ms, meta

    def __call__(self, x):
        x = np.asarray(x, dtype=np.float64)
        return np.interp(np.abs(x), self.xs, self.ms, left=self.ms[0], right=self.ms[-1])

    def to_dict(self):
        return dict(xs=self.xs.tolist(), ms=self.ms.tolist(), meta=self.meta)

    @staticmethod
    def from_dict(d):
        xs = np.asarray(d["xs"], dtype=np.float64)
        return Map(xs, xs, np.asarray(d["ms"], dtype=np.float64), d.get("meta", {}))


def fit_map(x: np.ndarray, y: np.ndarray, nbins: int = NBINS, statistic: str = "median",
            clip0: bool = True, meta: dict | None = None) -> Map:
    """Bin `x` by its own quantiles, take the per-bin `statistic` of `y`, monotonise."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    qs = np.linspace(0.0, 1.0, nbins + 1)
    edges = np.quantile(x, qs)
    edges[0] = -np.inf
    edges[-1] = np.inf
    b = np.searchsorted(edges, x, side="right") - 1
    b = np.clip(b, 0, nbins - 1)
    cx = np.zeros(nbins)
    cy = np.zeros(nbins)
    cw = np.zeros(nbins)
    order = np.argsort(b, kind="stable")
    bs = b[order]
    starts = np.searchsorted(bs, np.arange(nbins), side="left")
    ends = np.searchsorted(bs, np.arange(nbins), side="right")
    for i in range(nbins):
        sel = order[starts[i]:ends[i]]
        if sel.size == 0:
            cx[i] = np.nan
            continue
        cx[i] = float(np.median(x[sel]))
        cy[i] = float(np.median(y[sel]) if statistic == "median" else np.mean(y[sel]))
        cw[i] = float(sel.size)
    ok = np.isfinite(cx)
    cx, cy, cw = cx[ok], cy[ok], cw[ok]
    ms = pava(cy, cw)
    if clip0:
        ms = np.maximum(ms, 0.0)
    return Map(edges, cx, ms, dict(meta or {}, nbins=int(nbins), statistic=statistic,
                                  n=int(x.size)))


def fit_map_budget(x: np.ndarray, y: np.ndarray, nbins: int = NBINS,
                   meta: dict | None = None, tol: float = 1e-4) -> tuple[Map, float]:
    """The L1-optimal map SUBJECT TO an L2 budget — the first-order condition, solved.

    The plain conditional median (M) is the unconstrained L1 optimum, and it is ~20-50x
    SMALLER than the amplitude we ship, because `mse` (096's rho knee, board-validated to
    +/-1%) sets our amplitude and `nmae` would rather we predicted almost nothing.  So the
    honest question for a SHAPE arm is not "what is the L1 optimum" but "given this L2
    budget, where should the L1 go".  Minimising `sum_g E|x_g - y_g|` s.t. `sum_g x_g^2 = C`
    gives, per gene,

        2 F_g(x_g) - 1 + 2 lam x_g = 0   ->   F_g(x_g) = 1/2 - lam x_g

    and since the budget forces us ABOVE the unconstrained optimum, `lam < 0`: writing
    `lam = -mu`, **the optimal magnitude is the conditional `1/2 + mu*x` QUANTILE of `y`**,
    a quantile level that rises with our own predicted magnitude.  `mu` is solved by
    bisection so the mapped L2 equals the input L2 -- no post-hoc rescale, the map hits the
    budget itself.

    Returns the map and the solved `mu`.  This is a first-order condition on a
    difference-of-convex problem, so it is a candidate, not a proof of optimality; it is
    checked against the plain rescaled median on HELD-OUT L1 in `budget188.py`.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    qs = np.linspace(0.0, 1.0, nbins + 1)
    edges = np.quantile(x, qs)
    edges[0] = -np.inf
    edges[-1] = np.inf
    b = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, nbins - 1)
    order = np.argsort(b, kind="stable")
    bs = b[order]
    starts = np.searchsorted(bs, np.arange(nbins), side="left")
    ends = np.searchsorted(bs, np.arange(nbins), side="right")
    ys, cx, cw = [], np.zeros(nbins), np.zeros(nbins)
    for i in range(nbins):
        sel = order[starts[i]:ends[i]]
        ys.append(np.sort(y[sel]) if sel.size else np.zeros(1))
        cx[i] = float(np.median(x[sel])) if sel.size else np.nan
        cw[i] = float(sel.size)
    ok = np.isfinite(cx)
    target = float((x ** 2).sum())

    def build(mu):
        tau = np.clip(0.5 + mu * cx, 0.5, 0.9995)
        m = np.array([float(np.quantile(ys[i], tau[i])) if cw[i] else 0.0
                      for i in range(nbins)])
        m = np.maximum(pava(m[ok], cw[ok]), 0.0)
        return m

    def energy(mu):
        m = build(mu)
        return float((cw[ok] * m ** 2).sum())

    lo, hi = 0.0, 1.0
    while energy(hi) < target and hi < 1e7:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if energy(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol * max(hi, 1.0):
            break
    mu = 0.5 * (lo + hi)
    m = build(mu)
    return Map(edges, cx[ok], m, dict(meta or {}, nbins=int(nbins), statistic="budget",
                                      mu=float(mu), n=int(x.size))), float(mu)


def apply_map(P: np.ndarray, m: Map, gate: np.ndarray, mode: str = "row") -> np.ndarray:
    """`p -> sign(p) * m(|p|)` on `gate`, then restore each row's own L2 norm (`mode='row'`).

    Row restoration is what makes this a SHAPE arm: `rA` stays at 1, `VCC_ALPHA` is untouched,
    096's whole `s`-parameterisation of `mse` is unchanged, and rule 13's threshold is the
    shape-only +0.0019 rather than the amplitude-changing +0.0089.
    """
    one = P.ndim == 1
    A = np.atleast_2d(np.asarray(P, dtype=np.float64))
    G = np.where(np.atleast_1d(gate)[None, :], A, 0.0)
    out = np.sign(G) * m(np.abs(G))
    out = np.where(np.atleast_1d(gate)[None, :], out, 0.0)
    if mode == "row":
        n0 = np.linalg.norm(G, axis=1)
        n1 = np.linalg.norm(out, axis=1)
        out = out * (n0 / np.maximum(n1, 1e-300))[:, None]
    elif mode == "none":
        pass
    else:
        raise ValueError(mode)
    return out[0] if one else out


def to_lfc(H: np.ndarray, b: np.ndarray, gate: np.ndarray) -> np.ndarray:
    """`h -> log2( expm1(b+h) / expm1(b) )`, zero off the gate.  Exact, invertible."""
    A = np.atleast_2d(np.asarray(H, dtype=np.float64))
    base = np.expm1(np.asarray(b, dtype=np.float64))
    live = np.asarray(gate, dtype=bool) & (base > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        u = np.where(live[None, :], np.expm1(b[None, :] + A) / np.where(base > 0, base, 1.0), 1.0)
    return np.where(live[None, :], np.log2(np.clip(u, 1e-12, None)), 0.0)


def from_lfc(L: np.ndarray, b: np.ndarray, gate: np.ndarray) -> np.ndarray:
    base = np.expm1(np.asarray(b, dtype=np.float64))
    live = np.asarray(gate, dtype=bool) & (base > 0)
    return np.where(live[None, :], np.log1p(base[None, :] * np.exp2(L)) - b[None, :], 0.0)


def apply_map_lfc(P: np.ndarray, m: Map, b: np.ndarray, gate: np.ndarray,
                  restore: str = "h") -> np.ndarray:
    """The same map, but applied to log2 fold change; the BUILDER-space row norm is restored.

    169 §6.4: the builder's space is `log1p(TS*p)`, linear under ~20 CPM and logarithmic
    above, so a magnitude transform of `h` concentrates on high-EXPRESSION genes while the
    four DE members are all paid in fold-change space.  `restore='h'` keeps the arm
    amplitude-neutral in the space `VCC_ALPHA` acts in (169 measured that restoring the lfc
    norm instead smuggles in a 9% amplitude cut, i.e. two variables in one flag).
    """
    one = P.ndim == 1
    A = np.atleast_2d(np.asarray(P, dtype=np.float64))
    base = np.expm1(np.asarray(b, dtype=np.float64))
    live = np.asarray(gate, dtype=bool) & (base > 0)
    G = np.where(live[None, :], A, 0.0)
    L = to_lfc(G, b, live)
    Lq = np.sign(L) * m(np.abs(L))
    Lq = np.where(live[None, :], Lq, 0.0)
    out = from_lfc(Lq, b, live)
    if restore == "h":
        n0 = np.linalg.norm(G, axis=1)
        n1 = np.linalg.norm(out, axis=1)
        out = out * (n0 / np.maximum(n1, 1e-300))[:, None]
    elif restore != "none":
        raise ValueError(restore)
    return out[0] if one else out
