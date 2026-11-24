"""813 KA-SG statistic (frozen by sha256 in PROPOSAL-2.md before any real TeloHAEC row exists).

Question: are estimator A's multi-guide rows biased relative to rows from single-guide cells? Two noise-unbiased
statistics over a fixed target set T (the 100 stage-0 targets), on the scored genes G (5-CPM rule, panel columns out):

  A_t   = estimator A on carriers with n_guides >= 2 vs stratum-matched non-carriers with n_guides >= 2
  S_t   = carriers with n_guides == 1 vs non-carriers with n_guides == 1 (same pool/rep mixture)
  The two cell sets are DISJOINT (n_guides == 1 vs >= 2), so A's and S's sampling noises are independent.
  Each of the four cell sets (A carriers, A background, S carriers, S background) is split once into halves h1/h2 by a
  seeded hash of the cell barcode (seed 813); X_t^h is the row from half h of both its carriers and its background.

  beta  = sum_t <A_t, S_t> / sum_t <S_t^h1, S_t^h2>
          the slope of A's SIGNAL on S's signal. Independent noises drop out of the numerator in expectation, and the
          split-half cross product estimates S's signal energy without its noise, so noise does not attenuate beta
          (unlike an OLS slope). Unbiased A => E[beta] ~ 1.
  dsh   = <m^h1, m^h2> / mean_t <A_t^h1, A_t^h2>,   m^h = mean_t (A_t^h - S_t^h)
          the energy of a SHARED offset between A and S, per unit of A's per-target signal energy, noise-free in
          expectation. Unbiased A => E[dsh] ~ 0.
  What dsh cannot see: mu_bar, the mean effect of one random perturbation. A (k-matched background) and S (the k = 1
  background carries one other random guide) both carry -mu_bar, so it cancels in A - S. dsh sees only k-dependent
  technical offsets (guide load, capture, size, stress). mu_bar is G0-BLOB's to catch.

Bounds [lo, hi] for each statistic are NOT set here. Stage 0a sets them from KA-SIM (before any real row is read):
  the central 95 % interval over R = 200 simulated replicates of the unbiased design (real assignment and library
  sizes, the pilot's b_ctrl as base, K562-GWPS rows as truth), and a power row: the naive zero-guide estimator on the
  simulation with an injected guide-load offset must fall outside [lo, hi] in >= 80 % of replicates. If it does not,
  KA-SG has no power at this depth; it is recorded as uninformative and is not a gate (KA-SIM stays the estimator KA).
"""
from __future__ import annotations

import numpy as np

TS = 5e4


def row(p_t: np.ndarray, p_bg: np.ndarray) -> np.ndarray:
    """TS-5e4 log-pseudobulk difference of two summed-count vectors (the donor-table and pds space)."""
    return np.log1p(TS * p_t / p_t.sum()) - np.log1p(TS * p_bg / p_bg.sum())


def beta(A: np.ndarray, S: np.ndarray, S1: np.ndarray, S2: np.ndarray) -> float:
    """A, S, S1, S2: (n_targets, n_genes) rows on the scored genes."""
    return float(np.einsum("tg,tg->", A, S) / np.einsum("tg,tg->", S1, S2))


def dsh(A1: np.ndarray, A2: np.ndarray, S1: np.ndarray, S2: np.ndarray) -> float:
    m1, m2 = (A1 - S1).mean(0), (A2 - S2).mean(0)
    return float(m1 @ m2 / np.einsum("tg,tg->t", A1, A2).mean())


def bounds(sim_values: np.ndarray) -> tuple[float, float]:
    """[lo, hi] = central 95 % of the statistic over the unbiased simulation's replicates."""
    lo, hi = np.quantile(np.asarray(sim_values, float), [0.025, 0.975])
    return float(lo), float(hi)


def _selftest() -> None:
    """Noise-free identity rows: beta == 1 and dsh == 0 exactly; a shared offset b gives dsh == |b|^2 / mean|A|^2."""
    rng = np.random.default_rng(813)
    D = rng.normal(size=(50, 200))
    assert abs(beta(D, D, D, D) - 1) < 1e-12 and abs(dsh(D, D, D, D)) < 1e-12
    b = rng.normal(size=200) * 0.3
    v = dsh(D + b, D + b, D, D)
    assert abs(v - (b @ b) / np.einsum("tg,tg->t", D + b, D + b).mean()) < 1e-9
    print("kasg813 selftest PASS")


if __name__ == "__main__":
    _selftest()
