# 019 — Whitening flips the donor-weighting result, and resolves 010 vs 016

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations, truth = raw held-out effects.

## Why this had to be checked

Experiment 017 measured whitening *on top of* the incumbent recipe — `gene_norm`, donors at
4:1:1 — and both of those were chosen with no whitening in the picture. Donor weighting in
particular had a mechanism attached to it: experiment 010 found equal weighting **worse than
a single donor**, because averaging smooths every prediction toward the shared response and
destroys the discrimination `pds` measures.

That is the *same* shared-component problem whitening attacks directly. So the prediction was
that whitening should shrink the penalty for averaging — and it does far more than that.

## Result: equal weighting goes from worst to best

| primary share | gamma = 0 | gamma = 0.55 |
|---|---|---|
| 1:0 (single donor) | 0.7190 | 0.7637 |
| 8:1 | 0.7351 | 0.7789 |
| **4:1 (incumbent)** | **0.7378** | 0.7873 |
| 2:1 | 0.7308 | 0.7950 |
| **1:1** | 0.7210 | **0.7963** |
| 1:2 | 0.7083 | 0.7894 |

Without whitening the curve peaks at 4:1 and equal weighting costs 0.017. With whitening the
curve is monotone toward equal weighting and 4:1 costs 0.009. **The optimum moved from
"dominant primary" to "equal".**

Refining gamma at equal weights gives a broad plateau, best at **gamma = 0.65 → 0.7978**:

| gamma | 0.4 | 0.5 | 0.55 | **0.65** | 0.75 | 0.9 | 1.1 |
|---|---|---|---|---|---|---|---|
| macro PDS | 0.7868 | 0.7940 | 0.7963 | **0.7978** | 0.7971 | 0.7932 | 0.7849 |

**Incumbent 0.7378 → 0.7978, +0.0600 macro PDS**, on the harness that retro-predicted the v3
regression.

`gene_norm` still edges `raw` under whitening (0.7963 vs 0.7946 at 1:1), so it is kept — but
the margin is now small, consistent with whitening doing the same job properly: `gene_norm`
divides each gene by its mean |effect|, which is the *diagonal* of the covariance whitening
inverts in full.

## This resolves the 010 / 016 contradiction

Experiment 010 (2 predictors, VCC panel) found dominant-primary beats equal by +0.015.
Experiment 016 (3 predictors, mostly non-panel genes) found the reverse by +0.016. That was
logged as open, with the gene population as the suspected confound.

The real variable is neither: it is **whether the shared component has been neutralised**.
Averaging has two effects that pull opposite ways — it cancels independent noise (good) and
it concentrates predictions on the shared response (bad for a ranking metric). Whichever
dominates decides the winner. Whitening removes the second, and averaging wins cleanly.

## What ships

`gene_norm`, no EB, **three donors at 1:1:1**, whitening on the donor basis at
**gamma = 0.65**, `lam = 0`. `alpha` re-derived to **1.0796** so the emitted median effect
norm is 4.801 — identical to v4b — leaving *direction* as the single change from the last
built submission.
