# 051 — CD4 breadth-adaptive blend results

**Status:** Complete. The adaptive mechanism failed its predeclared gate; the
fixed `0.50` CD4 blend remains the production candidate.

## Frozen protocol

- Development folds: held HCT116 and held HEK293T.
- One-time lineage holdout: held K562, opened only after selection.
- Prediction modes: raw and gene-normalized.
- CD4 sources: Stim8hr and equal Rest/Stim8hr/Stim48hr mean.
- Fixed fractions: `0.10, 0.20, 1/3, 0.50, 0.67`.
- Adaptive slopes: `0.5, 1.0, 2.0`.
- QC strengths: `0.0, 0.5, 1.0`.
- Weight bounds: `[0.02, 0.80]`, with exact mean matching to the fixed
  fraction.
- Shuffled-breadth null: 31 deterministic permutations for the selected
  adaptive rule.
- Optional whitening: production gamma `0.65`, only with audited fold-specific
  union operators fitted without the held donor.
- Metric: effect-level PDS only.
- No `cell-eval2`; no submission.

## Development selection

| family | mode | CD4 source | geometry | fraction | slope | QC strength | mean PDS | min PDS | HCT116 delta | HEK293T delta | mean delta |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| best fixed | gene_norm | cd4_mean3 | none | 0.50 | — | — | 0.72718 | 0.71616 | +0.03836 | +0.01641 | **+0.02738** |
| best adaptive | gene_norm | cd4_mean3 | none | 0.50 | 0.5 | 0.0 | 0.73004 | 0.72023 | +0.04002 | +0.02047 | **+0.03025** |

Selected overall before opening K562: **adaptive, gene-normalized,
three-condition mean, fraction 0.50, slope 0.5, no QC term**.

## Selected adaptive controls on development folds

| held donor | baseline PDS | matched constant PDS | adaptive PDS | adaptive − constant | shuffled mean | shuffled p95 |
|---|---:|---:|---:|---:|---:|---:|
| HCT116 | 0.69984 | 0.73819 | 0.73986 | +0.00166 | 0.73564 | 0.73964 |
| HEK293T | 0.69975 | 0.71616 | 0.72023 | +0.00406 | 0.71482 | 0.71993 |

## One-time K562 lineage holdout

| frozen role | PDS | baseline PDS | delta vs baseline |
|---|---:|---:|---:|
| selected fixed | **0.77273** | 0.75280 | **+0.01992** |
| selected adaptive | 0.76993 | 0.75280 | +0.01713 |

| adaptive mechanism check | value |
|---|---:|
| adaptive − matched constant | **−0.00279** |
| adaptive − independently selected fixed | **−0.00279** |
| shuffled-breadth mean PDS | 0.76956 |
| shuffled-breadth p95 PDS | 0.77483 |
| adaptive − shuffled mean | +0.00037 |
| one-sided empirical permutation p | 0.50 |

## Predeclared interpretation

Call breadth adaptation a meaningful win only if all are true:

1. adaptive beats its mean-matched constant by at least `+0.003` PDS on K562;
2. adaptive exceeds the 95th percentile of its shuffled-breadth controls;
3. adaptive's mean development delta versus baseline is positive;
4. neither development fold loses more than `0.002` versus baseline; and
5. adaptive is not worse than the independently selected fixed rule by more
   than `0.002` on K562.

Anything smaller is unresolved. A K562-only gain cannot change the frozen
hyperparameters. Passing this screen makes the rule eligible for a separate
production/count-level test; it does not authorize a submission.

## Result

Breadth adaptation does **not** pass:

1. it loses `0.00279` PDS to its mean-matched fixed control on K562;
2. it is below the shuffled-breadth 95th percentile;
3. its permutation result is indistinguishable from chance (`p = 0.50`).

The development-fold gain does not transfer. This is a clean negative result:
do not use CD4 DE breadth or the tested QC fields to choose per-perturbation
blend weights.

The fixed `0.50` blend remains positive on both development folds and the
K562 target-role holdout. Carry only that fixed rule into a production-union
whitening test.

## Validation and memory

- Shared perturbations: 242.
- Master non-panel union genes: 17,942.
- Synthetic tests: 4 passed.
- Real run internal peak RSS: 2.82 GiB.
- External process-tree guard: 4 GiB; run completed without intervention.

K562 is a target-role holdout rather than a fully independent donor holdout,
because its effects participate as source inputs during development LODO. The
experiment therefore supports rejecting the adaptive mechanism and retaining
the fixed rule, but should not be described as a completely independent
generalization estimate.
