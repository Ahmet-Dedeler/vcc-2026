# 015 — Inverse-variance pooling and magnitude dispersion both fail

**Date:** 2026-08-21 · LODO macro PDS, 259 perturbations, 3 folds.

## Per-gene inverse-variance pooling

Donors differ in how well they measure each *gene* — 200 counts in K562, 8 in HCT116 — and a
fixed donor weight ignores that. Inverse-variance pooling is the minimum-variance combination of
unbiased estimates, and Sidechain (public code, holder of the board's best `fid`) does exactly
this. Poisson variance from each donor's own control profile.

| scheme | macro PDS |
|---|---|
| **fixed 4:1:1** | **0.7413** |
| fixed 8:1:1 | 0.7396 |
| inverse-variance × fixed 8:1:1 | 0.7387 |
| inverse-variance | 0.7217 |

**Pure inverse-variance pooling costs 0.020 PDS** — the largest single loss of any variant tested.
Blending it with fixed weights recovers most of that but still does not beat fixed weights alone.

Why, most likely: inverse-variance weighting is optimal for combining estimates of *the same*
quantity. Donor effects are not estimates of one quantity — they are genuinely different
biological responses in different lines. Weighting by measurement precision then systematically
over-weights whichever donor happens to express a gene highly, which is a property of that
donor's basal state rather than of the transferable signal. It is a variance argument applied to
a bias problem.

## Nadig's magnitude-dispersion signal

Nadig 2025 found cell-type-preferential perturbations have ~2× larger cross-line differences in
effect *magnitude* (p = 4e-8), suggesting dispersion as a confidence feature orthogonal to the
directional agreement of experiment 011.

| feature | Spearman vs transfer |
|---|---|
| directional agreement (exp 011) | **+0.324** |
| −magnitude dispersion | **−0.050** |
| agreement vs −dispersion | +0.081 |

**Magnitude dispersion carries essentially no signal here** (−0.050, wrong sign and negligible),
and the 0.081 cross-correlation shows it is not merely redundant with directional agreement —
it is independently uninformative.

Possible reconciliation: Nadig measured this on **common-essential** genes across four lines,
which is exactly the gene class the VCC panel excludes. Their result may hold and simply not
extend to non-essential targets. Worth noting rather than treating as a contradiction.

## Small correction to the donor weights

With three donors, **4:1:1 (0.7413) edges 8:1:1 (0.7396)**. The 2-donor sweep peaked at 4:1
(primary 80%); the 3-donor equivalent by primary share is 8:1:1, but the direct 3-donor test
prefers 4:1:1 (primary 67%). The gap is 0.0017 raw PDS ≈ +0.0006 overall — inside the plateau
and probably inside the noise — but 4:1:1 is the measured optimum, so v4 was rebuilt with it.
