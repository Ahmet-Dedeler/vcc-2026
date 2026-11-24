# 012 — Per-perturbation magnitude from agreement: mostly a dead end

**Date:** 2026-08-21 · LODO, 259 panel perturbations, agreement→α map fitted on the other two
folds only.

## Result

| held-out | uniform α | agreement α | oracle α*_p | fraction captured |
|---|---|---|---|---|
| k562_gwps | 2,995,662 | 2,995,669 | 2,988,223 | −0.1% |
| orion_hct116 | 3,954,224 | 3,952,706 | 3,942,967 | 13.5% |
| orion_hek293t | 3,521,905 | 3,520,422 | 3,513,925 | 18.6% |
| **mean** | | | | **10.7%** |

**The oracle itself is worth almost nothing: 0.25% squared-error reduction.** Even with perfect
per-perturbation magnitude — using the held-out truth — the gain is negligible, and agreement
captures a tenth of that.

## Why, and it generalises

At cos θ ≈ 0.08, `‖α·d − δ‖² ≈ ‖δ‖²` for any sensible α. The error is dominated by the
*directional* mismatch; magnitude is a second-order correction on top of a first-order failure.

**Per-perturbation magnitude calibration only has headroom once directional accuracy is decent.**
That is a general statement about where effort pays, and it argues against a whole family of
calibration refinements at our current cos.

Note the agreement signal itself is real and reproduces: Spearman(agreement, oracle α*) = **0.367**
(p = 4e-26), stronger than experiment 011's 0.324 against transfer cosine. The signal is sound;
the thing it controls just does not matter much yet.

## A bug in my own reasoning, caught by the data

I claimed per-perturbation scaling is `pds`-neutral because cosine distance is invariant to
scaling one perturbation's own effect. The k562 fold shows `pds` **dropping 0.7391 → 0.7186**
under agreement weighting.

The claim is true for *positive* scaling and false at zero. The fitted map produces negative α for
the lowest-agreement perturbations, and I clipped at 0 — which does not scale the effect, it
**deletes** it. A zero vector has no direction, ties against everything, and its `pds` collapses
to chance.

**Rule: clip per-perturbation α to a positive floor, never to zero.** Anything that can zero a
prediction is a `pds` hazard, not a magnitude adjustment.

## Curiosity worth noting

Oracle α*_p spread across perturbations: p10 **−1.46**, median 2.72, p90 11.10. The negative tail
is real — for those perturbations the donor consensus direction is *anti-correlated* with the
truth in the held-out line. Roughly 10–15% of perturbations would be better predicted by flipping
the donor's sign. We have no way to identify which in advance (agreement correlates at 0.367, not
enough to call individual signs), but it bounds how good uniform donor transfer can ever be.
