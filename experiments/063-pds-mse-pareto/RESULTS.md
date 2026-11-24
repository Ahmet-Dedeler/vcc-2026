# 063 — PDS/MSE Pareto blend results

Date: 2026-08-22

## Verdict

Keep the production union-whitened direction unchanged (`beta = 0`).

Adding any amount of the raw donor-consensus residual monotonically recovers a
small amount of true-effect cosine, but it loses substantially more PDS. Even
the fully raw endpoint has a best-case nonnegative-alpha MSE floor of
`0.995052`, still worse than the official `0.989` baseline and therefore still
worth zero scaled MSE. No tested blend can plausibly close the current
approximately `0.007` overall leaderboard gap.

## Exact macro Pareto points

All seven points are formally non-dominated because PDS decreases while cosine
increases and theoretical MSE decreases. That geometric frontier is not a
score frontier: every MSE value remains beyond the metric's zero-score clamp.

| raw share beta | macro PDS | raw cosine | min MSE, unrestricted alpha | min MSE, alpha >= 0 | optimistic overall delta vs beta 0 |
|---:|---:|---:|---:|---:|---:|
| **0.00** | **0.794219469** | 0.025326209 | 0.998368269 | 0.998384798 | **0.000000000** |
| 0.05 | 0.789141301 | 0.026108973 | 0.998266475 | 0.998290012 | -0.001870412 |
| 0.10 | 0.783160237 | 0.026878705 | 0.998149741 | 0.998183198 | -0.004073382 |
| 0.20 | 0.770968643 | 0.028373786 | 0.997872969 | 0.997935644 | -0.008563840 |
| 0.35 | 0.753404567 | 0.030487954 | 0.997358594 | 0.997489638 | -0.015033113 |
| 0.50 | 0.738908344 | 0.032425023 | 0.996745278 | 0.996973991 | -0.020372422 |
| 1.00 | 0.705406403 | **0.037427361** | **0.994347422** | **0.995052101** | -0.032711995 |

The smallest residual, `beta = 0.05`, gains only `+0.000782765` cosine and
improves the nonnegative-alpha MSE floor by `0.000094786`, while losing
`0.005078168` raw PDS. That PDS loss alone translates to approximately
`-0.001870` overall.

The raw endpoint gains `+0.012101153` cosine and improves the MSE floor by
`0.003332697`, but loses `0.088813066` raw PDS. Its MSE floor is still
`0.006052101` above the official baseline, so none of that theoretical
improvement becomes a scored MSE gain.

## Fold consistency

The tradeoff has the same direction on all three held donors.

| held donor | PDS beta 0 | PDS beta 1 | cosine beta 0 | cosine beta 1 | nonnegative min MSE beta 0 | nonnegative min MSE beta 1 |
|---|---:|---:|---:|---:|---:|---:|
| K562-GWPS | 0.799751579 | 0.742509952 | 0.031266262 | 0.036373406 | 0.997629865 | 0.996447760 |
| Orion-HCT116 | 0.800574661 | 0.684969022 | 0.023589999 | 0.038437250 | 0.998646711 | 0.994288592 |
| Orion-HEK293T | 0.782332166 | 0.688740235 | 0.021122365 | 0.037471428 | 0.998877816 | 0.994419950 |

PDS falls monotonically and cosine rises monotonically at every grid step in
every fold. This is a clean tradeoff, not a fold-specific cancellation.

## Why it cannot close the gap

The optimistic translation assumes oracle per-perturbation nonnegative alpha
and leaves the other four metrics unchanged. It uses the project-standard
midpoint anchors:

- PDS span: `0.4525`
- MSE baseline: `0.989`
- MSE replicate: `0.0365`

Every blend still has theoretical MSE above `0.989`, so its scaled MSE delta is
exactly zero. The only translated change is the PDS loss. The least harmful
candidate is therefore `beta = 0.05` at `-0.001870` overall, already moving in
the wrong direction. The fully raw endpoint is `-0.032712` overall.

Conclusion: residual de-whitening is closed for this production geometry.
Improving MSE requires a genuinely better biological direction, not mixing
back the anisotropic donor consensus that whitening removed.

## Protocol and leakage audit

- Three symmetric leave-one-donor-out folds.
- 259 VCC-panel targets shared by K562-GWPS, Orion-HCT116, and
  Orion-HEK293T.
- Equal source weights and raw centered effects, matching the production
  pre-whitening model.
- Source-only union covariance; fixed `gamma = 0.65`.
- All 300 challenge target coordinates excluded from every operator and every
  retrieval axis.
- Held-donor values opened only after the complete blend grid was frozen.
- Full 259-candidate gallery used for every PDS query.
- No challenge outcomes, counts, `cell-eval2`, or submission artifacts.

Leakage checks in `out/results.json` are all zero/true as required:

- held rows in operator: `0`
- panel overlap with any operator: `0`
- panel overlap with any retrieval axis: `0`
- predictions frozen before held values opened: `true`

## Compute and validation

- Fresh isolated Spot worker: `i-0f382ff6e160980e9`
- Type/region: `r7i.2xlarge`, `us-east-1`
- One heavy experiment process
- Internal/process-tree guard: `55 GiB`
- In-process peak RSS: `16.797455 GiB`
- Local synthetic tests: `6 passed`, `0.11 GiB` peak
- Cloud exit code: `0`
- S3 output prefix:
  `s3://<your-s3-bucket>/output/exp063-pds-mse-pareto/20260822T215018Z/`

Operator fingerprints:

| held donor | operator genes | eigenvalues SHA-256 | eigenvectors SHA-256 |
|---|---:|---|---|
| K562-GWPS | 17,806 | `cf771a7d9ace77d55b53671aa445c48fe2651157636f5d08491e03394c438a70` | `c90ad6890627bee25be0ba0b1adcbcfdf30818d3a66ab72e04970ba151feb3e9` |
| Orion-HCT116 | 17,942 | `eb803b7bda0417dfdc3c10e92f5c04f68bf4e5e998d646dfc6ec0d7a1159d3a2` | `497676e8d7971e5550193649a7331785819318d058a01b8d95e158bbb9afd623` |
| Orion-HEK293T | 17,942 | `1c6dc2f58cf5054f9ef44b0925f3bbedc255b53c7b84342c0bff982adb7f6485` | `e33718b7a286c633511ebf1ebcd6adf871429d7794b589bfeb77ca90b55faf42` |

Artifact hashes are recorded in `out/hashes.sha256` and verified after
download. The primary machine-readable outputs are `out/results.json`,
`out/summary.csv`, `out/gap_translation.csv`, and
`out/per_target_metrics.csv`.

