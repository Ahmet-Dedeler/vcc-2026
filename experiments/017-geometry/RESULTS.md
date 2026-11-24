# 017 — Spectral whitening of the prediction set: +0.060 macro PDS

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations, 3 donor lines, truth = raw
held-out effects. Incumbent to beat: **0.7378**.

## The idea, and why an old negative result predicted it

Perturbation effects are anisotropic. A few shared programmes — growth arrest, ribosome
biogenesis, stress — carry most of the variance, so two unrelated knockdowns that both slow
the cell down look alike under a cosine. `pds` is a **retrieval** question ("is this
prediction nearer its own true effect than to the other 299?"), and retrieval in an
anisotropic space suffers **hubness**: whichever true effects sit near the dominant programme
are near everything, and discrimination collapses.

Experiment 004 truncated the donor SVD and lost *monotonically* — rank 5 → cosine 0.039,
rank 800 → 0.052, full rank → 0.073. The standard reading was "the signal is
high-dimensional, stop compressing it". The stronger reading is that truncation **deletes**
the low-variance directions, so the opposite operation should **amplify** them:

```
P -> P (C + lam I)^(-gamma/2),      C = B^T B / n
```

`gamma = 0` is the incumbent, 1 is standard whitening, 2 a full Mahalanobis inverse. The
operator is invertible — nothing is discarded — and is fitted on donor screens only, so it
never sees the target line.

## Result

| gamma | 0 | 0.3 | 0.4 | **0.5** | **0.6** | 0.7 |
|---|---|---|---|---|---|---|
| macro PDS | 0.7378 | 0.7773 | 0.7834 | **0.7869** | **0.7873** | 0.7865 |

Consistent across all three folds (k562 0.7398→0.7870, hct116 0.7377→0.7978,
hek293t 0.7359→0.7770). Insensitive to the ridge `lam` over four orders of magnitude.

## The basis matters, and its behaviour is the evidence

| basis | shape | behaviour |
|---|---|---|
| the 300 panel predictions | 259 × 7,760 | rises monotonically to the edge of every grid tried (0.7848 at gamma=2, still climbing) |
| **all donor perturbations** | **12,000 × 7,760** | **clean interior optimum at gamma ≈ 0.5–0.65** |

A 7,760-dimensional covariance estimated from 259 samples is rank-deficient, and "keeps
helping as gamma grows" is the signature of an artefact rather than a mechanism — at the
limit it just orthogonalises 259 vectors, which improves a ranking by construction. Estimated
instead from every perturbation the donor lines measured, the covariance is full rank and the
sweep has a peak. **Use the donor basis.** That the two disagree in exactly this way is the
reason to trust the donor one.

## Harness correction made before believing any of this

`lodo.py` applied the same treatment (`gene_norm`) to the held-out line's effects as to the
predictors, making the *ground truth* a normalised object. The scorer does no such thing.
Whitening is precisely the kind of method that would be flattered by a normalised truth, so
the truth was switched to the raw held-out delta before the sweep. The effect survived it
(0.7378 → 0.7869 rather than 0.7412 → 0.7876); it is not an artefact of the harness.

## Two things that did not work, from the same run

| family | best | vs 0.7412 |
|---|---|---|
| rank-1 removal ("all-but-the-top") | 0.7420 at c=0.5 | +0.0008 |
| per-gene modulation by the **target's** basal profile, `p^beta` | 0.7443 at beta=0.25 | +0.0031 |
| hard CPM gate on the target's control profile | 0.7414 at >5 | +0.0002 |

Removing only the top direction captures almost none of what whitening does, which says the
hubness is spread across many directions rather than concentrated in one.

The basal-modulation result is the more interesting negative. Downweighting genes the target
does not express is biologically obvious — an unexpressed gene cannot move — and it is worth
+0.003 at `beta = 0.25` and actively harmful by `beta = 1` (0.7038). Combined with experiment
010's finding that basal similarity does not pick donors, **the control profile keeps failing
to be useful for anything except defining the gene universe**, which is worth stating plainly
because it is the only context-specific signal the Challenge gives us.

Gate > 5 edging gate > 20 is consistent with the scorer's own 5-CPM DE filter, but at +0.0002
it is inside noise; kept at 20 for v5 so that direction is the only change.
