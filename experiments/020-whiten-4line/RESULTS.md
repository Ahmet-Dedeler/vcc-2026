# 020 — Whitening replicates on a second, independent LODO

**Date:** 2026-08-21 · 4-line LODO (K562-gwps, Orion-HCT116, Orion-HEK293T, H1-2025),
226 perturbations shared by all four, **three predictors per fold** — matching the submission,
which uses three donors, where experiments 017/019 had only two.

## Result

| config | macro PDS |
|---|---|
| equal, gamma = 0 | 0.7874 |
| **equal, gamma = 0.55** | **0.8307** |
| equal, gamma = 0.65 | 0.8291 |
| primary 4:1:1, gamma = 0 | 0.7732 |
| primary 4:1:1, gamma = 0.55 | 0.8092 |

**+0.0433 from whitening**, and equal weighting beats dominant-primary at *every* gamma
including 0. Consistent across all four folds.

## Why this is worth having as a separate run

This harness differs from 017/019 in two ways that could each have killed the effect:

- **Three predictors instead of two.** Averaging three sources cancels more independent noise,
  so the prediction handed to the whitener is already cleaner and the operator has less to fix.
  It still gains 0.043.
- **A different gene population.** Only 20 of the 226 shared perturbations are on the VCC
  panel — H1-2025's panel is essential-gene enriched, exactly the class the VCC panel
  excludes. So this is largely a *different biology* from the task.

The effect size is smaller here (+0.043 vs +0.060), which is what the first difference
predicts. The mechanism holds in both.

Experiment 019 remains the one that sets the shipped numbers, because its 259 perturbations
are actual VCC-panel targets. This one says the mechanism is not an artefact of that panel.

## Note on the 016 contradiction

Experiment 016 used this same 4-line harness and found equal weighting beating dominant
primary, contradicting experiment 010's VCC-panel result. Two candidate explanations were
logged, with the gene population favoured. Experiment 019 shows the real variable is
**whether the shared component has been removed** — and here at gamma = 0, with three
predictors, equal already wins by 0.014. So predictor count is the operative difference, not
gene population: the more sources you average, the more the noise-cancelling benefit of equal
weighting outruns its smoothing cost. Whitening then removes the cost entirely.
