# 063 — PDS/MSE Pareto blend

Predeclared 2026-08-22, before any held-donor effect values were evaluated.

## Question

Production union whitening improves retrieval PDS but reduces raw alignment to
the true perturbation effect. Test whether a convex residual blend can recover
useful raw cosine, and therefore theoretical MSE headroom, without giving back
too much PDS:

```text
q_beta = (1 - beta) * q_whitened + beta * q_raw
```

`q_raw` is the exact equal-source, centered, raw donor consensus used before
the production union whitener. `q_whitened` applies the source-only union
operator at fixed `gamma = 0.65` with production median-norm preservation.
There is no row-wise alignment or post-hoc rescaling before the blend.

## Frozen grid

The raw-residual shares are fixed to:

```text
beta = [0.00, 0.05, 0.10, 0.20, 0.35, 0.50, 1.00]
```

No grid point will be added after inspecting held-donor outcomes. `beta=0` is
the production union-whitened geometry and `beta=1` is the unwhitened donor
consensus.

## Evaluation

- Outer holdouts: K562-GWPS, Orion-HCT116, and Orion-HEK293T, symmetrically.
- Sources: the other two donors, equal weight.
- Gallery: the VCC-panel targets shared by all three donors.
- Operator fit: all source-donor perturbation rows only.
- Operator/retrieval coordinates: official VCC gene order with **all 300 panel
  target coordinates removed before the dense transform**.
- Truth: untouched centered raw effect from the held donor, opened only after
  the raw and whitened predictions are frozen.
- No counts, `cell-eval2`, challenge outcomes, or submission artifacts.

For every beta and fold, report:

1. full-gallery rank PDS;
2. mean row-wise raw cosine;
3. mean per-perturbation theoretical minimum uncorrected normalized MSE,
   `mean(1 - cos_i^2)`;
4. the same MSE with the optimal scalar constrained non-negative, so
   anti-correlated rows cannot claim unattainable benefit.

The Pareto frontier maximizes PDS and cosine while minimizing both MSE
diagnostics. A candidate is considered capable of plausibly closing a `0.007`
overall leaderboard gap only if the optimistic two-metric translation

```text
(
  (PDS_beta - PDS_beta0) / 0.4525
  + (MSE_scaled_beta - MSE_scaled_beta0)
) / 6
```

is at least `0.007`, where `MSE_scaled = clip((0.989 - min_mse_nonnegative) /
0.9525, 0, 1)`. This deliberately assumes the other four metrics are
unchanged and therefore favors the blend.

## Resource rule

The faithful union eigendecomposition exceeds the 1 GiB local cap. Real-data
evaluation must run as one guarded process on a fresh isolated AWS Spot worker
with a 55 GiB process-tree limit. Local work is restricted to synthetic tests
under a 1 GiB guard.

