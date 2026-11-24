"""296 -- shared helpers.  Every piece of response arithmetic and every statistic comes from 285's
`lib285.py` verbatim (which imports 218's `_per_target` verbatim); nothing mathematical is new here.

New here, and only here:
  * `strata()`  -- the frozen knockdown stratification of PREREGISTRATION sec 2.3 (cuts C1/C2/C3),
                   re-implemented rather than imported from 293 so that KA-3 is a real check.
  * `two_sample_ratio_diff()` -- the STRATUM bootstrap of rho_guide(both) - rho_guide(weaker fails)
                   between two DISJOINT target sets (PREREGISTRATION sec 2.2).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/296-guide-efficacy"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "experiments/285-reliability-benchmark/work"))
import lib285 as L  # noqa: E402

NDRAW = L.NDRAW          # 8
SEED = L.SEED            # 20260906
NBOOT = L.NBOOT          # 20000


# ---------------------------------------------------------------------------- the frozen stratification
def rho_of(b: np.ndarray, a: np.ndarray) -> float:
    """rho_guide = median(b) / median(a) on the targets where both are finite (293's construction)."""
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() == 0:
        return np.nan
    ma = float(np.median(a[m]))
    return float(np.median(b[m]) / ma) if ma != 0 else np.nan


def rho_boot(b: np.ndarray, a: np.ndarray, n: int = NBOOT, seed: int = SEED) -> dict:
    """Bootstrap interval of the ratio of medians, resampling TARGETS (293 strat293.py / 285 k562_guide.py)."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size == 0:
        return dict(point=np.nan, lo=np.nan, hi=np.nan, n=0)
    rng = np.random.default_rng(seed)
    ii = rng.integers(0, a.size, size=(n, a.size))
    r = np.median(b[ii], 1) / np.median(a[ii], 1)
    return dict(point=float(np.median(b) / np.median(a)), lo=float(np.percentile(r, 2.5)),
                hi=float(np.percentile(r, 97.5)), n=int(a.size))


def two_sample_ratio_diff(b1, a1, b2, a2, n: int = NBOOT, seed: int = SEED) -> dict:
    """STRATUM bootstrap of rho_guide(stratum 1) - rho_guide(stratum 2) over two DISJOINT target sets.

    Targets are resampled independently inside each stratum; no pairing exists between disjoint sets."""
    m1 = np.isfinite(a1) & np.isfinite(b1)
    m2 = np.isfinite(a2) & np.isfinite(b2)
    a1, b1, a2, b2 = a1[m1], b1[m1], a2[m2], b2[m2]
    if a1.size == 0 or a2.size == 0:
        return dict(point=np.nan, lo=np.nan, hi=np.nan, n1=int(a1.size), n2=int(a2.size), p_le0=np.nan)
    rng = np.random.default_rng(seed)
    i1 = rng.integers(0, a1.size, size=(n, a1.size))
    i2 = rng.integers(0, a2.size, size=(n, a2.size))
    d = (np.median(b1[i1], 1) / np.median(a1[i1], 1)) - (np.median(b2[i2], 1) / np.median(a2[i2], 1))
    d = d[np.isfinite(d)]
    pt = float(np.median(b1) / np.median(a1) - np.median(b2) / np.median(a2))
    return dict(point=pt, lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5)),
                n1=int(a1.size), n2=int(a2.size), p_le0=float((d <= 0).mean()))


def two_sample_median_diff(x1, x2, n: int = NBOOT, seed: int = SEED) -> dict:
    """STRATUM bootstrap of median(x1) - median(x2) over two DISJOINT sets."""
    x1 = np.asarray(x1, float); x2 = np.asarray(x2, float)
    x1, x2 = x1[np.isfinite(x1)], x2[np.isfinite(x2)]
    if x1.size == 0 or x2.size == 0:
        return dict(point=np.nan, lo=np.nan, hi=np.nan, n1=int(x1.size), n2=int(x2.size), p_le0=np.nan)
    rng = np.random.default_rng(seed)
    i1 = rng.integers(0, x1.size, size=(n, x1.size))
    i2 = rng.integers(0, x2.size, size=(n, x2.size))
    d = np.median(x1[i1], 1) - np.median(x2[i2], 1)
    return dict(point=float(np.median(x1) - np.median(x2)), lo=float(np.percentile(d, 2.5)),
                hi=float(np.percentile(d, 97.5)), n1=int(x1.size), n2=int(x2.size), p_le0=float((d <= 0).mean()))


def strata_masks(kd_big: np.ndarray, kd_small: np.ndarray, cut: float,
                 kf_big: np.ndarray | None = None, kf_small: np.ndarray | None = None,
                 kf_cut: float = 0.5) -> dict:
    """PREREGISTRATION sec 2.3.  `kd_*` are the constructs' on-target log deltas (signed);
    the stratifier is |lfc| of the WEAKER construct."""
    w = np.minimum(np.abs(kd_big), np.abs(kd_small))
    out = {"both": w >= cut, "weak": w < cut}
    if kf_big is not None:
        wf = np.minimum(kf_big, kf_small)
        out["both_frac"] = wf >= kf_cut
        out["weak_frac"] = wf < kf_cut
    return out


def summarise(name: str, sel: np.ndarray, d: dict, n_for_sb: float | None = None) -> dict:
    """One stratum row: (a), (a2), (b), (x), rho_guide, rho_dis, pds (a)/(b)/(x)."""
    g = {k: np.asarray(d[k], float)[sel] for k in ("a", "a2", "b", "x", "pds_a", "pds_b", "pds_x")}
    med = lambda v: float(np.nanmedian(v)) if np.isfinite(v).any() else np.nan  # noqa: E731
    ma, ma2, mb = med(g["a"]), med(g["a2"]), med(g["b"])
    rho = rho_boot(g["b"], g["a"])
    rdis = (float(mb / np.sqrt(ma * ma2)) if np.isfinite(ma) and np.isfinite(ma2) and ma > 0 and ma2 > 0
            else None)
    row = dict(name=name, k=int(sel.sum()),
               a=L.boot_median(g["a"]), a2=L.boot_median(g["a2"]), b=L.boot_median(g["b"]),
               x=L.boot_median(g["x"]), rho_guide=rho, rho_dis=rdis,
               pds_a=med(g["pds_a"]), pds_b=med(g["pds_b"]), pds_x=med(g["pds_x"]),
               b_minus_x=L.boot_median_diff(g["b"], g["x"]))
    if n_for_sb:
        k = 400.0 / n_for_sb
        row["sb_a_400"] = float(np.nanmedian([L.spearman_brown(v, k) for v in g["a"] if np.isfinite(v)]))
        row["sb_b_400"] = float(np.nanmedian([L.spearman_brown(v, k) for v in g["b"] if np.isfinite(v)]))
    return row


def fmt_rho(r: dict) -> str:
    return f"{r['point']:+.3f} [{r['lo']:+.3f}, {r['hi']:+.3f}]"


def spearman(x, y) -> tuple[float, float]:
    from scipy.stats import spearmanr
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan
    r = spearmanr(x[m], y[m])
    return float(r.statistic), float(r.pvalue)
