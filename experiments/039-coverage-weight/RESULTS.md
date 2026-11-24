# 039 — The implicit shrinkage of under-covered genes is worth +0.019. Leave it alone.

**Date:** 2026-08-22 · production gene axis, union operator, 259 VCC-panel perturbations, folds
holding out an Orion line.

## The thing that looks like a bug

`MultiDonorModel.effect` divides by the total weight over donors carrying the *perturbation*,
then adds each donor's contribution only on the genes that donor *measures*. K562-gwps carries
7,679 of the official axis; each Orion line carries 18,111. So for any gene gwps does not carry:

```
(1/3)*0 + (1/3)*hct116 + (1/3)*hek293t  =  2/3 of the Orion average
```

**Roughly 59% of the scored gene axis is emitted at reduced magnitude**, because a donor that
never measured the gene is counted as a vote of zero. That was invisible while the whitening
operator only spanned the intersection; experiment 038 put those genes inside the operator, so
they now carry real weight in the direction.

## Result

| variant | hct116 | hek293t | macro |
|---|---|---|---|
| **as-is** — divide by n_donors, gamma 0.65 | **0.8053** | **0.7831** | **0.7942** |
| per-gene — divide by actual coverage, gamma 0.65 | 0.7888 | 0.7612 | 0.7750 |
| as-is, no whitening | 0.7030 | 0.7061 | 0.7045 |
| per-gene, no whitening | 0.6922 | 0.6957 | 0.6939 |

**"Fixing" it costs 0.019 macro PDS**, on both folds and at both gammas.

## Why, and it is the interesting part

Dividing by the donors that *actually cover* a gene gives every gene an unbiased mean of
whatever measured it — which is the statistically correct thing and is worse. Dividing by the
full donor count instead shrinks a gene in proportion to how few donors saw it, so a gene
measured by one donor is emitted at 1/3 magnitude and a gene measured by all three at full
magnitude.

That is **exactly a confidence weighting**, and it is the one the rest of this pipeline keeps
trying and failing to build deliberately. Experiment 022 found per-gene confidence schemes lose
to uniform treatment; experiment 034c found targeted emission variance lands inside the noise of
a uniform spread; experiment 035 found cross-donor agreement separates sign accuracy by only 3.7
points. Every explicit attempt to weight genes by how much we trust them has failed — and an
implicit one, arrived at by accident, is worth +0.019.

The difference is what the signal is. Our explicit schemes weighted by *agreement between
donors*, which is a noisy per-perturbation quantity. This weights by **how many donors measured
the gene at all**, which is a fixed, noiseless property of the assay panels. Coverage is a much
better-estimated confidence signal than concordance.

## Consequence

No change to the model. Recorded so nobody "corrects" it later — and as the one piece of
evidence that per-gene confidence weighting *can* work here, provided the confidence comes from
something we know exactly rather than something we estimate.
