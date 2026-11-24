"""193 — eigenvalue maps `f(ev)` that can replace production's `(C + lam I)^(-gamma/2)`.

Every whitening operator in this project has the form

    W = V f(ev) V^T ,     production:  f(ev) = (ev + lam)^(-gamma/2),  gamma=0.65, lam=0.0

so an alternative covariance regulariser is a drop-in: estimate a better set of eigenvalues
`d(ev)` and use `f = d^(-gamma/2)` with gamma held at production's value.

TWO INVARIANCES THAT MATTER FOR READING THE SCREEN
--------------------------------------------------
1. `plolo.Operator.apply` restores the input's median row norm after whitening, so a map and
   `c * map` are the SAME arm for every metric, not just for `pds`.  Any estimator that only
   differs from another by a global scale therefore cannot be screened apart.
2. Because of (1), **Ledoit-Wolf LINEAR shrinkage is algebraically inside the ridge family**:
       d = (1-rho) ev + rho mu   ==   (1-rho) * (ev + lam),   lam = rho*mu/(1-rho)
   so a lam-sweep already contains every linear-shrinkage intensity.  LW-linear, RBLW and OAS
   are therefore reported as *marked points on the lam axis*, not as separate arms.  The only
   estimators that can escape the ridge family are the NONLINEAR ones.

WHAT IS ACTUALLY IN THE SPECTRUM (measured, `work/floor-census.log`)
--------------------------------------------------------------------
`Whitener.fit` floors the eigenvalues at `ev.max() * 1e-12`, and every operator on disk hits
that floor exactly: condition number `1.000e+12` on all six.  The donor covariance has a
near-null tail (12.5 % of directions below `1e-6` relative on the production operator), and
because `f = ev^-0.325` is a *decreasing* map, **99.7 % of the operator's total squared gain
sits in that tail**.  Production's `lam = 0.0` means the only thing standing between the
prediction and a 4,450x amplification of the covariance's numerical null space is the relative
floor.  That is the object this experiment is testing.
"""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------------- linear family
def ridge(ev: np.ndarray, lam: float) -> np.ndarray:
    """d = ev + lam.  lam = 0 is production."""
    return ev + lam


def lw_linear_equivalent_lam(ev: np.ndarray, rho: float) -> float:
    """The ridge `lam` that reproduces linear shrinkage intensity `rho` up to a global scale."""
    mu = float(ev.mean())
    return rho * mu / (1.0 - rho)


def rblw_rho(ev: np.ndarray, n: int) -> float:
    """Rao-Blackwell Ledoit-Wolf shrinkage intensity (Chen et al. 2010, eq. 17).

    Depends only on `tr(C)` and `tr(C^2)`, so it is computable from a frozen spectrum.
    """
    p = len(ev)
    t1 = float(ev.sum())          # tr C
    t2 = float((ev * ev).sum())   # tr C^2
    num = ((n - 2.0) / n) * t2 + t1 * t1
    den = (n + 2.0) * (t2 - t1 * t1 / p)
    return float(min(1.0, num / den)) if den > 0 else 1.0


def oas_rho(ev: np.ndarray, n: int) -> float:
    """Oracle-approximating shrinkage intensity (Chen et al. 2010, eq. 23)."""
    p = len(ev)
    t1 = float(ev.sum())
    t2 = float((ev * ev).sum())
    num = (1.0 - 2.0 / p) * t2 + t1 * t1
    den = (n + 1.0 - 2.0 / p) * (t2 - t1 * t1 / p)
    return float(min(1.0, num / den)) if den > 0 else 1.0


# ------------------------------------------------------------------------- flooring / clipping
def floor_quantile(ev: np.ndarray, alpha: float) -> np.ndarray:
    """Hard floor at the alpha-quantile of the spectrum: d = max(ev, q_alpha).

    A ridge is a *soft* floor that also perturbs the head; this one leaves the head exactly
    alone.  If the near-null tail is what hurts, this separates the two effects.
    """
    q = float(np.quantile(ev, alpha))
    return np.maximum(ev, q)


def floor_relative(ev: np.ndarray, rel: float) -> np.ndarray:
    """Hard floor at `rel * ev.max()` -- production is this map at rel = 1e-12."""
    return np.maximum(ev, float(ev.max()) * rel)


def winsorize(ev: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.clip(ev, float(np.quantile(ev, lo)), float(np.quantile(ev, hi)))


def clip_rank(ev: np.ndarray, k: int) -> np.ndarray:
    """Keep the top-`k` eigenvalues, replace the rest by their mean (trace-preserving).

    This is eigenvalue clipping / the spiked-factor-model target: the bulk is declared noise
    and flattened rather than merely damped.  `ev` is assumed ASCENDING (numpy `eigh`).
    """
    d = np.array(ev, dtype=np.float64, copy=True)
    p = len(d)
    k = int(max(0, min(k, p)))
    if k < p:
        d[: p - k] = d[: p - k].mean()
    return d


def mp_sigma2(ev: np.ndarray, q: float, iters: int = 500) -> tuple[float, float, int]:
    """Iteratively fit the Marchenko-Pastur noise level.

    sigma^2 is the mean of the eigenvalues that fall below the current MP edge; the edge is
    `sigma^2 (1 + sqrt(q))^2`.  Returns (sigma2, edge, n_above).  On a spectrum with no MP
    bulk this iteration collapses toward the smallest eigenvalues, which is itself the
    diagnostic -- see `RESULTS.md`.
    """
    s2 = float(ev.mean())
    edge = s2 * (1.0 + np.sqrt(q)) ** 2
    for _ in range(iters):
        edge = s2 * (1.0 + np.sqrt(q)) ** 2
        bulk = ev <= edge
        if not bulk.any():
            break
        s2n = float(ev[bulk].mean())
        if abs(s2n - s2) <= 1e-15 * max(s2, 1e-300):
            s2 = s2n
            break
        s2 = s2n
    edge = s2 * (1.0 + np.sqrt(q)) ** 2
    return s2, edge, int((ev > edge).sum())


def clip_mp(ev: np.ndarray, q: float) -> np.ndarray:
    """MP denoising: everything at or below the fitted MP edge becomes the bulk mean."""
    _, edge, _ = mp_sigma2(ev, q)
    d = np.array(ev, dtype=np.float64, copy=True)
    m = d <= edge
    if m.any():
        d[m] = d[m].mean()
    return d


# ------------------------------------------------------------- analytical nonlinear shrinkage
def analytical_nls(ev: np.ndarray, n: int, block: int = 512,
                   h: float | None = None) -> np.ndarray:
    """Ledoit & Wolf (2020) analytical nonlinear shrinkage, `p <= n` branch.

    Direct port of the authors' `analytical_shrinkage.m`, chunked over rows so the p x p
    kernel matrix is never materialised (p = 17,942 would be 2.6 GiB per temporary).

        ftilde(l_i)  = (3/(4 sqrt5)) * mean_j  max(1 - x_ij^2/5, 0) / (h l_j)
        Hftilde(l_i) = mean_j [ (-3/(10 pi)) x_ij
                                + (3/(4 sqrt5 pi))(1 - x_ij^2/5) log|(sqrt5-x)/(sqrt5+x)| ] / (h l_j)
        x_ij = (l_i - l_j) / (h l_j),     h = n^(-1/3)

        d_i = l_i / [ (pi c l_i ftilde)^2 + (1 - c - pi c l_i Hftilde)^2 ],   c = p/n

    `h` may be overridden to probe the bandwidth's influence; the paper's choice is n^(-1/3).
    """
    lam = np.asarray(ev, dtype=np.float64)
    p = len(lam)
    c = p / float(n)
    if h is None:
        h = float(n) ** (-1.0 / 3.0)
    hl = h * lam                                   # (p,) bandwidth at each l_j
    inv_hl = 1.0 / hl
    ftil = np.empty(p, dtype=np.float64)
    hftil = np.empty(p, dtype=np.float64)
    s5 = np.sqrt(5.0)
    k1 = 3.0 / (4.0 * s5)
    k2 = -3.0 / (10.0 * np.pi)
    k3 = 3.0 / (4.0 * s5 * np.pi)
    for lo in range(0, p, block):
        hi = min(p, lo + block)
        x = (lam[lo:hi, None] - lam[None, :]) * inv_hl[None, :]
        ftil[lo:hi] = k1 * np.mean(np.maximum(1.0 - x * x / 5.0, 0.0) * inv_hl[None, :], axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            lg = np.log(np.abs((s5 - x) / (s5 + x)))
        tmp = k2 * x + k3 * (1.0 - x * x / 5.0) * lg
        on_edge = np.abs(np.abs(x) - s5) < 1e-300
        tmp[on_edge] = k2 * x[on_edge]
        tmp[~np.isfinite(tmp)] = 0.0
        hftil[lo:hi] = np.mean(tmp * inv_hl[None, :], axis=1)
        del x, lg, tmp
    denom = (np.pi * c * lam * ftil) ** 2 + (1.0 - c - np.pi * c * lam * hftil) ** 2
    d = lam / np.maximum(denom, 1e-300)
    return d


def rie_lp(ev: np.ndarray, n: int, eta_pow: float = 0.5, block: int = 512) -> np.ndarray:
    """Ledoit-Peche / Bun-Bouchaud-Potters rotationally-invariant estimator.

        xi_i = l_i / |1 - c + c * l_i * g(z_i)|^2,   z_i = l_i - i*eta,   eta = n^(-eta_pow)

    with `g` the empirical Stieltjes transform `(1/p) sum_j 1/(z - l_j)`.  Same family as
    `analytical_nls`, regularised by a Cauchy width instead of an Epanechnikov bandwidth; the
    two agreeing is evidence the estimate is real rather than a kernel artefact.
    """
    lam = np.asarray(ev, dtype=np.float64)
    p = len(lam)
    c = p / float(n)
    eta = float(n) ** (-eta_pow)
    out = np.empty(p, dtype=np.float64)
    for lo in range(0, p, block):
        hi = min(p, lo + block)
        # eta scales with the local eigenvalue so the tail is not over-smoothed
        e = eta * np.maximum(lam[lo:hi], 1e-300)
        z = lam[lo:hi][:, None] - 1j * e[:, None]
        g = np.mean(1.0 / (z - lam[None, :]), axis=1)
        denom = np.abs(1.0 - c + c * lam[lo:hi] * g) ** 2
        out[lo:hi] = lam[lo:hi] / np.maximum(denom, 1e-300)
        del z, g
    return out


# ------------------------------------------------------------------- two-sided shape family
def piecewise_power(ev: np.ndarray, theta_lo: float, theta_hi: float,
                    pivot_q: float = 0.5) -> np.ndarray:
    """Bend the spectrum about its own median without changing the level of flattening.

        d = m * (ev/m)^theta_lo   for ev <= m ;   d = m * (ev/m)^theta_hi   for ev > m

    `theta_lo = theta_hi = 1` IS production.  `theta_lo < 1` compresses the near-null tail --
    the same thing a ridge, a floor or any shrinkage estimator does, only smoothly.
    `theta_lo > 1` does the OPPOSITE: it stretches the tail so the whitener amplifies the
    covariance's weakest directions even harder than `ev^-gamma/2` already does.

    That second branch is the reason this family exists.  Every regulariser in the literature
    lives on the `theta < 1` side, so a sweep of them can only ever show that production is at
    a BOUNDARY.  Sweeping both sides is what turns "no shrinkage estimator helps" into
    "production sits at an interior optimum".
    """
    m = float(np.quantile(ev, pivot_q))
    r = np.maximum(ev, 1e-300) / m
    lo = ev <= m
    d = np.empty_like(ev)
    d[lo] = m * r[lo] ** theta_lo
    d[~lo] = m * r[~lo] ** theta_hi
    return d


# ------------------------------------------------------------------------------- registry
def build_maps(ev: np.ndarray, n: int, gamma: float = 0.65) -> dict[str, np.ndarray]:
    """`name -> f(ev)` for every screened map, at production's gamma."""
    p = len(ev)
    q = p / float(n)
    ev = np.asarray(ev, dtype=np.float64)
    out: dict[str, np.ndarray] = {}

    def add(name: str, d: np.ndarray) -> None:
        out[name] = np.asarray(d, dtype=np.float64) ** (-gamma / 2.0)

    add("prod(lam=0)", ridge(ev, 0.0))
    for lam in (1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 7e-3, 2e-2, 6e-2):
        add(f"ridge{lam:g}", ridge(ev, lam))
    for rel in (1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        add(f"floorrel{rel:g}", floor_relative(ev, rel))
    for a in (0.05, 0.10, 0.25, 0.50, 0.75):
        add(f"floorq{a:g}", floor_quantile(ev, a))
    for a in (0.90, 0.99):
        add(f"wins0.25-{a:g}", winsorize(ev, 0.25, a))
    for k in (128, 512, 2048, 8192, 15695):
        add(f"cliprank{k}", clip_rank(ev, k))
    add("clipMP", clip_mp(ev, q))
    add("nls-LW2020", analytical_nls(ev, n))
    add("nls-LW2020-h2x", analytical_nls(ev, n, h=2.0 * float(n) ** (-1.0 / 3.0)))
    add("rie-LP-eta0.5", rie_lp(ev, n, 0.5))
    add("rie-LP-eta0.33", rie_lp(ev, n, 1.0 / 3.0))
    return out


def build_maps_b(ev: np.ndarray, n: int, gamma: float = 0.65) -> dict[str, np.ndarray]:
    """Second registry: the two-sided shape scan and a finer rank/ridge grid.

    Run separately so the (already measured) first registry is not recomputed.
    """
    ev = np.asarray(ev, dtype=np.float64)
    out: dict[str, np.ndarray] = {}

    def add(name: str, d: np.ndarray) -> None:
        out[name] = np.asarray(d, dtype=np.float64) ** (-gamma / 2.0)

    add("prod(lam=0)", ridge(ev, 0.0))                        # identity row, re-measured
    for t in (0.80, 0.90, 0.95, 1.05, 1.10, 1.20):
        add(f"tail^{t:g}", piecewise_power(ev, t, 1.0))
    for t in (0.80, 0.90, 0.95, 1.05, 1.10, 1.20):
        add(f"head^{t:g}", piecewise_power(ev, 1.0, t))
    for lam in (3e-6, 3e-5):
        add(f"ridge{lam:g}", ridge(ev, lam))
    for k in (1024, 4096, 12000, 16000, 17500):
        add(f"cliprank{k}", clip_rank(ev, k))
    return out
