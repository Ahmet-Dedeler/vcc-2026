# 101 — cross-fitted selective wrong-sign correction

Status: **preregistered before any experiment-101 evaluation output is read**.

Date: 2026-08-29.

## Question

Can a low-capacity, cross-fitted model isolate a small downstream-gene tail
whose frozen donor-transfer signs are wrong, and can reversing only that tail
improve experiment-098 linear NMAE credit without harming complete-gallery
PDS?

For a predicted LFC `h[p,g]`, the signed training target is:

```text
y[p,g] = sign(reference_lfc[p,g]) * sign(base_lfc[p,g]) in {-1,+1}
```

The fitted ridge score `q(x)` estimates signed correctness. A candidate acts
only where `q < -tau`:

```text
flip:  h' = -lambda * h
zero:  h' = 0
other: h' = h
```

`zero` is a matched control. It is never eligible for model selection: a gain
from zeroing without a gain from flipping is ordinary shrinkage, not learned
directional correction.

## Frozen grid

- STRING/context feature rank: **16 only**.
- Weighted ridge penalties: `{0.1, 1, 10}`.
- Confidence thresholds `tau`: `{0.10, 0.25, 0.50}`.
- Flip strengths `lambda`: `{0.5, 1.0}`.
- Controls: identity and zero-selected-tail.
- Whole-target folds: four deterministic SHA-256 folds.
- Target-label permutation nulls: 20 deterministic complete-block mappings.
- Target-block bootstrap draws: 2,000.

No grid value may be added after an outer result is read.

## Hard leakage boundary

The outer surfaces are H1, HCT116, and HEK293T. For outer surface `d`:

1. `d` is absent from every training label, donor feature, inner prediction,
   fitted ridge parameter, and hyperparameter decision.
2. Training labels come only from the other two reference-DE caches.
3. A labeled context is predicted from donors excluding both `d` and that
   context's own truth donor.
4. Every fit uses four whole-target folds. For fold `k`, all occurrences of
   the held target symbols are removed from all fitting contexts.
5. All 300 challenge-panel genes are removed from labels, features used for
   intervention, linear-credit evaluation, and PDS coordinates.
6. The outer surface's perturbation values are not opened until the selected
   real model, selected candidate, null models, null-selected candidates, and
   all outer predictions or model parameters are frozen on disk.

Any failed assertion stops the run.

## Features and estimator

Experiment 100's exact rank-16 feature vector is reused:

- seven donor-magnitude/agreement features;
- source-mean sign;
- target-context control expression;
- STRING target/response cosine;
- 16 target STRING PCs;
- 16 response STRING PCs;
- 16 coordinatewise target-response PC products;
- two STRING missingness flags.

The feature table and donor loading code are imported from experiment 100
after redirecting all of its caches into this experiment directory.

Each target has equal reference scale. Pair weight is:

```text
abs(base_lfc[p,g]) / sum_g abs(reference_lfc[p,g])
```

Features are standardized with those weights. A deterministic weighted ridge
regression is fit directly to `y in {-1,+1}`.

## Nested selection

For each outer surface there are two inner held-context cells. In each cell:

1. fit on the other inner context only;
2. preserve the four whole-target folds;
3. predict the held inner context from a donor set excluding both the outer
   donor and held inner truth donor;
4. evaluate every frozen candidate on experiment-098 signed linear credit and
   complete-gallery PDS.

The available precomputed whitening operators include the outer held donor in
some nested cells. Using them would leak. Re-fitting a full dense whitener is
outside this experiment's 28 GiB boundary. Therefore all inner selection
cells use the same explicitly preregistered **unwhitened, median-norm-scaled**
base. The outer evaluation uses experiment 100's exact surface-specific
whitened base. This is a conservative proxy mismatch, not a donor leak.

A flip candidate is eligible only if its complete-gallery PDS delta versus
identity is nonnegative in both inner cells. Among eligible candidates, select
the largest macro signed-credit gain. If the best gain is not positive, select
identity. Exact ties prefer stronger ridge, larger `tau`, then smaller
`lambda`.

## Frozen outer evaluation

For the selected candidate and its matched zero-tail control, report:

1. experiment-098 linear NMAE proxy;
2. complete-gallery PDS;
3. weighted error rate inside the selected tail;
4. selected pair count and target coverage;
5. per-surface and macro target-block bootstrap intervals;
6. 20 complete target-label permutation-null outcomes.

The primary evidence is the selected flip's linear-credit gain over identity.
The directional correction is a **GO** only if:

- raw linear-credit gain is at least `0.003` on every surface;
- its macro 95% target-block-bootstrap lower bound is above zero;
- weighted wrong-sign rate in the flipped tail exceeds 60% on every surface;
- complete-gallery PDS delta is nonnegative on every surface;
- observed macro linear-credit gain exceeds permutation-null p95;
- flipping beats the matched zero-tail control in macro linear credit.

## Resource and execution boundary

- Local only, one surface at a time.
- Every substantive command runs under
  `scripts/run_with_memory_guard.py --max-rss-gib 28`.
- No `cell-eval2`, AWS, challenge build, submission, or commit.
- All writes stay below `experiments/101-wrong-sign-correction/`.

