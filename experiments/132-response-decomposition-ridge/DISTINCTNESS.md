# Experiment 132 — distinctness audit

Status: **scientifically distinct enough to run**

Date: **2026-08-30**

## Proposed mechanism

The candidate is the exact low-capacity model family used by
`xinyizhanglab/perturbation-decomposition` at commit
`a15214780619736d393f40240e56ba992fd416a3`:

1. represent each perturbation target by its DepMap CRISPR dependency profile;
2. fit PCA with 50 components on the source target features;
3. fit multi-output Ridge with `alpha=10` and an intercept from those components
   to raw, uncentered genome-wide perturbation effects;
4. apply the learned map to the same target in a strictly excluded outer
   context.

The authenticated Experiment 114 feature artifact is used as the target
representation. It is row-centered and L2-normalized before this experiment;
that upstream-input difference is explicit. The model call, capacity, and
hyperparameters remain the upstream PCA-50 plus Ridge-10 defaults.

## Difference from Experiment 097

Experiment 097's linear route used PCA scores of a **donor perturbation-response
matrix** as the target embedding. Its shippable form learned a map from donor
A's response PCs to donor B's response vector. That asks whether one response
geometry can be linearly translated into another.

Experiment 132 instead uses **external DepMap target features** as the input.
It pools source contexts and learns a target-feature-to-response map. The input
space, training examples, and biological hypothesis are different:

- Experiment 097 input: observed transcriptomic response of the target in a
  donor screen.
- Experiment 132 input: the target gene's dependency profile across 1,178
  DepMap cell lines.
- Experiment 097 hypothesis: donor response geometry transfers linearly.
- Experiment 132 hypothesis: DepMap target geometry predicts the conserved
  genome-wide response shared across contexts.

Experiment 097's negative result therefore does not test this model.

## Difference from Experiment 114

Experiment 114 used DepMap features only as a cosine-neighbor index. It
replaced 41 weak incumbent rows with a weighted average of neighboring source
responses, retained the incumbent norm, and fit no response map.

Experiment 132 is a supervised global low-rank regression. Every training
target contributes to a learned 50-dimensional linear map, and every eligible
row is predicted through that map. It does not use nearest neighbors,
weak-row routing, norm restoration, STRING features, or a residual blend.

Experiment 114's failed kNN replacement therefore does not test this model.

## Decision

Proceed with one source-only nested evaluation. No additional architecture,
rank, ridge, blend, mask, or dose candidate is allowed. Identity remains in
the candidate set and wins unless the Ridge candidate clears every frozen
gate.
