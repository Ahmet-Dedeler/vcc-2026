# 100 — cross-fitted downstream sign confidence

Status: **preregistered before any experiment-100 evaluation output was read**.

Date: 2026-08-29.

## Question

Can a low-capacity model predict whether the frozen donor-transfer sign is
correct for a `(perturbation target, response gene)` pair better than the
existing donor-agreement signal alone, without using the held context's truth
at inference?

This is a confidence model, not a new expression predictor. It may only
reweight the magnitude of the existing effect. It may not flip signs.

## Hard leakage boundary

The three final evaluation surfaces are:

| surface | held truth | frozen inference donors |
| --- | --- | --- |
| `h1` | H1 test reference | K562-GWPS + Orion-HEK293T |
| `hct116` | clean Orion-HCT116 reference | K562-GWPS + Orion-HEK293T |
| `hek293t` | clean Orion-HEK293T reference | K562-GWPS + Orion-HCT116 |

For an outer surface `d`:

1. `d` is excluded from every training feature, label, donor statistic, and
   fitted parameter.
2. Training labels come only from the other two reference-DE caches.
3. A training context is itself predicted without its own donor. Its source
   pool is the available subset of K562-GWPS, H1, HCT116, HEK293T, and KOLF2
   after removing both the outer held donor and the training truth donor.
4. Final targets are assigned to four deterministic hash folds. For fold `k`,
   every occurrence of those target symbols is removed from training across
   all contexts. Predictions for fold `k` therefore come from a model that has
   never seen correctness labels for those perturbations.
5. Samples are whole `(target, response-gene)` observations from the
   reference-significant gate. Random pair-level train/test splitting is
   forbidden.
6. All 300 VCC panel genes are removed as response coordinates from features,
   labels, transformed vectors, NMAE-proxy evaluation, and PDS evaluation.
   The held harness's own target-gene coordinates are additionally removed
   from PDS, matching the gallery exclusion rule.
7. Held-surface truth is opened only for final scoring after all cross-fitted
   probabilities and transformed predictions for that surface are frozen.

Any failed assertion on these boundaries stops the run.

## Data

- Compact centered donor effects from `data/donors/` and
  `data/external/kolf2/`.
- H1, clean HCT116, and clean HEK293T reference-DE caches already produced by
  experiments 098/099.
- Target-context control pseudobulks already present in those harnesses.
- Pinned STRING 64-dimensional embedding:
  `data/external/slim/gene_string_embeddings.v0.3.h5`.

No challenge artifact, leaderboard result, hidden gallery statistic, or held
reference label is an inference feature.

## Features

For each target, source-donor effects are centered within donor, mapped to the
target gene axis, median-norm normalized, and uniformly averaged.

The **agreement-only** model receives:

1. `log1p(abs(source mean))`;
2. `abs(source mean) / mean(abs(source donor values))`;
3. source-sign unanimity;
4. minimum source magnitude divided by mean source magnitude;
5. source standard deviation divided by mean source magnitude;
6. number of donors carrying the perturbation;
7. mean pairwise source-donor cosine for the perturbation.

The **full** model adds:

1. sign of the source mean;
2. target-control `log1p(CPM)`;
3. STRING target/response cosine;
4. target STRING PCs;
5. response-gene STRING PCs;
6. coordinatewise target × response STRING-PC products;
7. target/response STRING missingness flags.

STRING PCA is unsupervised and fitted once to the pinned embedding table,
without response labels. The primary model uses rank 8. Rank 16 is a
preregistered sensitivity analysis, not an alternative selected after seeing
the final metrics.

## Model and fixed intervention

- Estimator: standardized ridge logistic regression, `C=1`, L2 penalty,
  deterministic `lbfgs`, at most 200 iterations.
- Training weight: each target has equal total weight; within a target,
  examples are weighted by the source prediction's absolute magnitude divided
  by the reference's total absolute LFC on its eligible DE gate.
- Output: probability that the frozen source sign is correct.
- Reweighting: within each target, standardize the logit over eligible response
  genes, clip to `[-2, 2]`, and multiply the frozen post-whitening effect by
  `exp(0.15 * z)`. Restore the row's pre-intervention L2 norm over eligible
  genes. Panel coordinates remain exactly zero. No sign can change.

The same intervention is applied to the agreement-only model. Identity uses a
multiplier of one.

## Metrics

All metrics are computed per held target and then aggregated:

1. **Weighted sign precision** on the top 20% of reference-DE pairs ranked by
   confidence. Weights are the experiment-098 linear-credit weights,
   `abs(identity_lfc_hat) / sum(abs(reference_lfc))`, with equal target weight.
   Identity's all-pair weighted sign precision is also reported.
2. **Experiment-098 linear NMAE proxy**:
   `1 - mean_p(dot(sign(reference_lfc), predicted_lfc) /
   sum(abs(reference_lfc)))`. Lower is better.
3. **Complete-gallery PDS**, from every target in the surface's gallery, with
   the VCC panel coordinates and held gallery target coordinates removed.

Target-block bootstrap uncertainty uses 2,000 deterministic resamples of whole
targets. Pair-level bootstrap is forbidden.

## Null

Twenty deterministic target-label permutations are run for the primary
rank-8 full model. Within each training context, the mapping from complete
target blocks of correctness labels to target-feature blocks is permuted.
Response-gene identity and class balance are retained; target correspondence
is destroyed. Every null repeats all four target-fold fits.

## Decision rule

The primary rank-8 model is a **GO** only if all are true:

1. weighted top-20% sign precision is no worse than agreement-only on H1,
   HCT116, and HEK293T;
2. linear NMAE proxy is no worse than identity on all three surfaces;
3. complete-gallery PDS is no worse than identity on all three surfaces;
4. the macro weighted-precision gain over agreement-only has a target-block
   bootstrap 95% lower bound above zero;
5. the observed macro weighted-precision gain exceeds the permutation-null
   95th percentile;
6. rank 16 does not reverse the sign on any surface.

Any negative surface is a NO-GO. Results smaller than `0.003` in PDS or linear
proxy are reported as unresolved even if numerically positive.

## Resource and execution boundary

- Local only; no AWS.
- One surface process at a time.
- Soft peak-RSS stop at 28 GiB.
- No `cell-eval2` CLI, count-level scorer, challenge build, submission, or git
  commit.
- All new code, logs, caches, and results stay under
  `experiments/100-sign-confidence/`.
