# 052 — Source subspace versus union whitening results

Date: 2026-08-22

## Verdict

Do not replace or compose v9's union whitener with the rank-64 source
subspace transform.

The fixed diagnostic comparison showed a small apparent gain for applying the
subspace transform before union whitening, but the source-only nested selector
could not identify that candidate. Under both independent split seeds it
selected a mixture of `union_only` and `subspace_only` and lost on every outer
fold relative to union-only.

| fixed pipeline | macro PDS | delta vs union |
|---|---:|---:|
| identity | 0.70474 | -0.08073 |
| **union only (v9 geometry)** | **0.78547** | 0 |
| subspace only | 0.75234 | -0.03313 |
| union then subspace | 0.78788 | +0.00241 |
| subspace then union | 0.79050 | +0.00503 |

The last row is diagnostic only. It was not selected without access to the
outer VCC-panel outcomes.

## Nested selection

Both seeds (`52`, `152`) produced the same outer choices:

- held K562: `union_only`;
- held HCT116: `union_only`;
- held HEK293T: `subspace_only`.

For both seeds:

- selected macro PDS: `0.77148`;
- delta versus union-only: `-0.01399`;
- positive outer folds versus union-only: `0/3`;
- predeclared nested-selection gate: **failed**.

The matched random controls were worse, so experiment 048's source-subspace
signal was not merely a random-projection artifact. The failure here is more
specific: once compared against the full union operator, source-only inner
selection cannot choose a composition that generalizes.

## Leakage and validation

The corrected run:

- ranked each inner query against the full 600-candidate gallery;
- froze inner selections before loading the 259-target outer responses;
- excluded all 900 VCC plus inner-panel target coordinates from both
  transforms;
- did not center with held-donor background effects;
- fit every operator without the held donor;
- used no `cell-eval2`, count emission, or submission artifact.

Synthetic tests: 10 passed.

## Memory

The first cloud attempt deliberately stopped when cumulative peak RSS reached
3.96 GiB on a wide Orion fold. The complete cloud run used a 12 GiB external
guard and peaked at:

- in-process peak RSS: `3.97 GiB`;
- externally sampled process-tree peak: `2.98 GiB`.

Keep the 4 GiB local command fail-closed. Run this experiment only on a cloud
worker if it must be repeated.

## Decision

The v9 union whitener remains the selected geometry. The fixed
`subspace_then_union` number is too small and not source-selectable, so it is
not eligible for a challenge submission.
