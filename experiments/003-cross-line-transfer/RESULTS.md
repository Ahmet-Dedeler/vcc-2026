# 003 — How well does a perturbation response transfer between cell lines?

**Date:** 2026-08-21 · **Data:** Replogle 2022 — K562 genome-wide vs RPE1, same lab, same
protocol, same year, different tissue of origin. **2,187 perturbations and 7,094 genes measured
in both lines.** This is the cleanest public read on the exact axis the 2026 Challenge tests.

## Why this is the number that matters

Experiment 002 established that every scored metric is a function of directional accuracy
`cos θ` between the predicted and true effect vector. The simplest possible model — look up the
same gene's response in a donor line and copy it — achieves exactly `cos(δ_donor, δ_target)`.
So this measurement *is* the score of naive donor transfer, and the headroom above it is what a
real model has to earn.

## Result

| quantity | what it is | mean cos | median | frac > 0 |
|---|---|---|---|---|
| **transfer** | same gene, K562 → RPE1 | **+0.1416** | +0.1247 | 0.90 |
| meanresp | vs RPE1's mean response | +0.0971 | +0.0881 | 0.77 |
| mismatch | different gene (null) | +0.0437 | +0.0314 | 0.69 |

**Naive donor transfer earns cos ≈ 0.14.** That is low in absolute terms, and worse than it
looks: most of it is not gene-specific.

## The decomposition that matters

Write a line's effect as `δ_g = m + s_g` — a shared mean response plus a gene-specific deviation.

**`m` earns exactly zero score.** It *is* the scored baseline: the context-mean response, the
0 end of every metric's scale. Any predicted component along `m` is wasted effort.

So the useful signal is `s`. Removing each line's mean response from both sides:

| | mean cos | median | p10 | p90 | frac > 0 |
|---|---|---|---|---|---|
| **centered transfer** | **+0.1275** | +0.1192 | −0.0182 | +0.2823 | 0.86 |
| centered null (permuted) | +0.0034 | +0.0032 | −0.1122 | +0.1257 | 0.52 |

**The gene-specific component transfers at cos ≈ 0.128 against a null of 0.003.** Forty times the
null, present for 86% of perturbations — unambiguously real, and unambiguously small.

Note what the raw number was hiding: of the +0.142 raw transfer, +0.097 is the shared mean
response that scores nothing. `transfer − meanresp` is only **+0.045**.

## Negative result: effect strength does not predict transferability

Experiment 002 implied a per-perturbation `α*_p` weighted by how well each direction is predicted.
The obvious proxy is donor effect magnitude — strong, clean perturbations ought to transfer better.

**They do not.** Spearman(donor `|δ_centered|`, centered transfer cos) = **0.063**.

| quintile of donor effect magnitude | n | mean transfer cos |
|---|---|---|
| Q1 (weakest) | 438 | 0.1280 |
| Q2 | 437 | 0.1089 |
| Q3 | 437 | 0.1017 |
| Q4 | 437 | 0.1373 |
| Q5 (strongest) | 438 | 0.1616 |

Non-monotone, and the Q5−Q1 gap (0.03) is small next to the spread within any quintile. **Donor
effect size is not a usable confidence signal.** Per-perturbation α needs a different predictor —
pathway membership, donor replicate agreement, or a learned one.

Best-transferring perturbations are dominated by core RNA-processing and mitochondrial machinery:
SMG5 (0.72), CMTR1 (0.62), GEMIN5 (0.61), UPF2 (0.61), CLOCK (0.61), TFAM (0.60), POLRMT (0.58),
CASC3 (0.57), MTOR (0.56), METTL3 (0.51). Housekeeping machinery behaves the same everywhere;
that is where transfer is easy, and it suggests **pathway identity is the confidence signal to
try next**.

## What this means for the plan

1. **cos ≈ 0.13 is the naive-transfer operating point**, not the 0.5 that experiment 002 swept.
   Re-run the calibration grid there — round 1's numbers are from an unrealistically easy regime.
2. **Predict `m` and `s` separately.** They have completely different roles: `m` sets the
   expression profile (drives `mse`) and earns no score of its own; `s` earns everything.
   A model that conflates them wastes capacity on the part that scores zero.
3. **RPE1 and K562 are far apart** (leukemia vs retinal epithelium). The 2026 lines may be closer
   to K562, in which case 0.128 is a floor. Worth checking against the H1 2025 data for the 25
   panel genes it shares.
4. The headroom is large in relative terms — the replicate anchor corresponds to near-perfect
   directional agreement — so almost all the achievable score is still unclaimed. But nobody
   should expect naive donor lookup to get there.

## Reproduce

```bash
uv run python scripts/build_donor_table.py --donor all
uv run python experiments/003-cross-line-transfer/measure_transfer.py
```
