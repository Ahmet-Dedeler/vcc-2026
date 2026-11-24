# 024 — Sizing the whitening operator: basis size matters, per-donor whitening does not

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations, equal donor weights.

## What the spectrum said first

The fitted covariance has top eigenvalue 3.90, median 0.0036, minimum 1.7e-4, condition number
23,000. The top 10 directions carry 16% of the variance, the top 100 carry 24%. A ridge of
0.01 — **above the median eigenvalue**, so it flattens more than half the spectrum — moved the
result by 0.0003.

That pins the mechanism from the operator side: **the gain comes from suppressing the dominant
directions, not from differentiating among the small ones.** Which is the hubness story, and
it made three specific refinements worth testing.

## 1. Basis size — the one that paid

| perturbations per donor | basis rows | best gamma | macro PDS |
|---|---|---|---|
| 2,000 | 4,000 | 0.55 | 0.7632 |
| 6,000 (previous default) | 12,000 | 0.65 | 0.7995 |
| **all** | **27,554–35,744** | **0.85** | **0.8078** |

**+0.008 over the 6,000 cap**, and the optimal gamma moves 0.65 → 0.85 with it. That coupling
is the expected one and is itself evidence the operator is doing what it claims: a
better-conditioned covariance can be inverted harder before noise in the small eigenvalues
starts to dominate. The 6,000 cap was arbitrary and was costing both the samples and the
gamma.

Full-basis gamma sweep, broad plateau across 0.75–1.05:

| gamma | 0.65 | 0.75 | **0.85** | 0.95 | 1.05 | 1.2 |
|---|---|---|---|---|---|---|
| macro PDS | 0.8043 | 0.8071 | **0.8078** | 0.8074 | 0.8068 | 0.8050 |

## 2. Per-donor whitening — a clear loss

Each donor line has its own dominant programme, reflecting its own growth rate and basal
state, so whitening each donor by *its own* covariance before averaging looked like the more
careful version. (Whitening every donor by the same operator is algebraically identical to
whitening the average, so only the per-donor-covariance form is a distinct thing to test.)

| | macro PDS |
|---|---|
| pooled operator | 0.7995 |
| per-donor operators | 0.7376 |

It costs 0.062 — nearly the entire whitening gain. Whitening each donor separately destroys
the *comparability* between them: after each is mapped by a different linear operator, the two
predictions no longer live in a common geometry, and averaging them mixes incompatible spaces.
The shared programme is a property of perturbation biology, not of any one line, and it should
be removed once, jointly.

## 3. Percentile flooring — inside noise

`w_i = max(ev_i, ev_q)^(-gamma/2)`, a scale-free version of the ridge: 0.8000 (q10), 0.8011
(q25), 0.8011 (q50) against 0.7995 unfloored. +0.0016 at best, which is inside the resolution
floor. Consistent with the ridge sweep, and with the mechanism — the small eigenvalues are not
where the action is. Not adopted; one fewer parameter.
