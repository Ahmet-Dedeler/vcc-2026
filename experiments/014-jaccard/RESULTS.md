# 014 — `jac` is not escapable at our precision (contradicts a research recommendation)

**Date:** 2026-08-21 · LODO, held out, 3 folds.

## The claim under test

A research sweep proposed the highest-value lever on the board: every team's raw `jac` sits at
0.029–0.033 against a baseline of 0.030, because minimum-variance emission calls the whole 5-CPM
gate (`R̂` ≈ 10,000 vs `R` ≈ 300, so `J → |R|/10,000 ≈ 0.03`). Call exactly `n_real` genes instead
and `J = p/(2−p)`, worth **+0.066 at p=0.3, +0.137 at p=0.5** — more than the entire gap from
rank 5 to rank 1.

The mechanism is correct. The arithmetic is correct. **The assumed precision is not ours.**

## Measured precision

| held-out line | precision at k=n_real |
|---|---|
| k562_gwps | 0.0522 |
| orion_hct116 | 0.0529 |
| orion_hek293t | 0.0507 |
| **macro** | **0.0519** |

Random selection of 300 from 7,581 genes gives 0.0396. **We are 1.31× chance at identifying
*which* genes respond.**

| precision | J | scaled | overall |
|---|---|---|---|
| **0.0519 (measured)** | 0.0266 | **−0.009** | **−0.0015** |
| 0.10 | 0.0526 | +0.061 | +0.0102 |
| 0.30 (assumed) | 0.1765 | +0.397 | +0.0662 |
| 0.50 (assumed) | 0.3333 | +0.822 | +0.1370 |

**At our real precision the proposed strategy scores worse than what we already do (+0.0036).**

## The sweep says the opposite of the recommendation

| k called | precision | recall | J | scaled |
|---|---|---|---|---|
| 300 | 0.0531 | 0.0531 | 0.0274 | −0.007 |
| 1,000 | 0.0441 | 0.1471 | 0.0352 | +0.014 |
| 3,000 | 0.0399 | 0.3993 | 0.0377 | **+0.021** |

Because precision is so low, **recall dominates the Jaccard** and calling *more* genes beats
calling fewer. Max realistic gain ≈ +0.021 scaled = **+0.0035 overall**, not +0.066.

## Why this matters beyond `jac`

`jac` asks *which genes respond*; `pds` asks *in what pattern*. They are different questions, and
we are near chance on the first one. That is the same directional-accuracy ceiling showing up in
a third guise — after `mse` (needs cos > 0.105, we have 0.08) and per-perturbation magnitude
(exp 012, oracle worth 0.25%).

**Three independent metrics now say the same thing: nothing downstream of the effect estimate can
be fixed by calibration. The binding constraint is the quality of the transferred direction.**

The board agrees — `score_jac` max across all 84 teams is +0.0076, and `score_fid` is negative
for **every single team**. These are not unexploited opportunities; they are where the biology
runs out.

## Method note

The proposal came from a research agent with correct mechanics and an assumed input. Checking the
assumption cost one script and 10 minutes; acting on it would have cost a submission slot and,
given the k-sweep, would have *lowered* our score. Assumed quantities in a recommendation are
where to look first.
