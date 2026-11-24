# 022 — Whitening does not resurrect the graveyard, and that is the point

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations, truth = raw held-out effects.

## Why run this

Equal donor weighting was the *worst* configuration in experiment 010 and the *best* in
experiment 019. Nothing about the data changed — whitening did. Several other rejected ideas
carry the **same stated mechanism**, so either they flip too, or whitening is doing something
specific. Those are very different worlds and the difference is measurable.

The prediction, written before running: ideas rejected for *shared-component / smoothing*
reasons should move; ideas rejected for *accuracy* or *information-loss* reasons should not.
A result where everything improves in rank would mean whitening merely adds slack.

## Result (macro PDS)

| combination | gamma = 0 | gamma = 0.65 | rank change |
|---|---|---|---|
| **gene_norm, equal** | 0.7210 | **0.7978** | worst → **best** |
| gene_norm, primary 4:1 | 0.7378 | 0.7870 | best → 2nd |
| gene_norm, inverse-variance | 0.7142 | 0.7910 | unchanged (−0.0068 both) |
| gene_norm + EB, equal | 0.7073 | 0.7826 | still loses |
| EB, equal | 0.6970 | 0.7668 | still loses |
| gene_norm, rank-200 truncation | 0.7046 | 0.7828 | still loses |
| gene_norm, rank-50 truncation | 0.6158 | 0.6677 | still loses badly |
| gene_norm, sign only | 0.6609 | 0.6823 | still loses badly |

Whitening lifts **every** configuration by +0.04 to +0.08 — but it re-ranks exactly one thing:
donor weighting. Everything else keeps its order, and inverse-variance pooling keeps its
deficit to the *third decimal* (−0.0068 at both gammas), which is as clean a statement of
orthogonality as this harness can produce.

## What that buys

1. **Whitening is not slack.** If it were, the weaker methods would close their gaps. They do
   not. It fixes one specific failure — the shared programme dominating a cosine — and leaves
   every other failure mode intact.
2. **The graveyard stays closed.** EB shrinkage, inverse-variance pooling, low-rank
   truncation and magnitude-discarding are all still rejected, now under the best available
   geometry. They do not need re-testing again.
3. **One refinement to the 019 conclusion.** The equal-weight flip is specific to
   `gene_norm`. Under EB, dominant-primary still wins at gamma = 0.65 (0.7745 vs 0.7668), and
   under `gene_norm + EB` it still wins narrowly (0.7842 vs 0.7826). EB introduces its own
   smoothing that whitening does not undo — consistent with EB shrinking toward zero
   preferentially on low-abundance genes, which *concentrates* the prediction on the abundant
   shared programme rather than dispersing it.

So the shipped configuration is `gene_norm`, no EB, equal weights, gamma = 0.65 — and the
three parts are not independently swappable.
