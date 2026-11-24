# 016 — MODZ per-perturbation donor weighting does not beat fixed weights

**Date:** 2026-08-21 · 4-line LODO (K562-gwps, Orion-HCT116, Orion-HEK293T, H1-2025),
226 perturbations shared across all four, 3 predictors per fold.

## Why this needed a fourth line

MODZ (CMap L1000, `cmapPy/math/agg_wt_avg.py`) weights each source by its mean correlation with
the *others*, per signature — the published prior art for experiment 011's agreement finding and
the top recommendation from the research sweep.

It was untestable with three lines: MODZ needs ≥3 sources to differentiate, because with two,
`corr(A,B)` is a single number and the weights are always equal. Holding one of three lines out
leaves two. Building H1-2025 into a donor table gives a fourth line and three predictors per fold.

## Result

| held out | equal | primary 4:1:1 | MODZ | MODZ × primary |
|---|---|---|---|---|
| k562_gwps | **0.8134** | 0.7633 | 0.7957 | 0.7546 |
| orion_hct116 | 0.7639 | 0.7412 | **0.7660** | 0.7438 |
| orion_hek293t | **0.7448** | 0.7423 | 0.7442 | 0.7376 |
| h1_2025 | 0.8397 | **0.8493** | 0.8232 | 0.8413 |
| **macro** | **0.7904** | 0.7740 | 0.7823 | 0.7693 |

**MODZ loses to plain equal weighting** (0.7823 vs 0.7904), and combining it with fixed weights is
worst of all (0.7693). The per-perturbation adaptivity does not pay.

Consistent with experiment 012: the agreement signal is real (Spearman 0.32–0.37) but what it
controls — how much to trust one source over another for one perturbation — is second-order next
to the directional error that dominates everything at our accuracy.

## The result that contradicts experiment 010, and why I do not trust the contradiction

Experiment 010 (3 lines, 2 predictors, **259 VCC-panel perturbations**) found a dominant primary
beats equal weighting by +0.015. This experiment (4 lines, 3 predictors, **226 mostly-non-panel
perturbations**) finds the opposite by +0.016.

Two differences, either of which could explain it:

1. **Number of predictors.** With two sources, equal weighting smooths away discrimination
   (exp 010's finding). With three, there is more independent noise to cancel and averaging may
   genuinely win. Plausible and would be a real effect.
2. **Gene population — the confound I would bet on.** H1-2025's panel overlaps the VCC panel by
   only 25 genes, so of the 226 perturbations shared across all four lines, **only 20 are on the
   VCC panel**. This test is largely measuring a *different class of genes* — H1's panel is
   essential-gene-enriched, exactly the class the VCC panel excludes.

**Experiment 010 is the one that matches our actual task**, so v4 keeps the dominant-primary
weighting. But this is an open question, not a settled one, and it is worth resolving with a
fourth line whose panel actually overlaps the VCC targets.

## What this closes

MODZ was the strongest single recommendation to come out of the research sweep — published,
battle-tested on ~1.3M CMap signatures, ~20 lines of code, and a direct formalisation of our own
experiment 011. It does not help here. Worth knowing that the prior art having a name does not
make it the right tool for this data.
