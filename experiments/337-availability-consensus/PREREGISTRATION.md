# 337 — available-weight normalization of a zero-filled donor consensus

September 9, 2026. Freeze before construction or scoring. No response-model
fit, new source, reserve, production default, package or submission.

## Question and prior distinction

`MultiDonorModel.effect` divides by the total weight of donors carrying a
target, even where a carrying donor does not measure a particular gene.
Zero-filling therefore lowers the local aggregate weight on partial-support
coordinates. Does availability normalization improve the current family's
actual-count output, or is that attenuation useful regularization?

This is NOT a proven implementation defect: missing low-expression genes
are not missing at random. It adds no biological information. It changes
an aggregation assumption, not a neural architecture or scalar confidence.

097's pairwise agreement modifies both missingness and dissent at fully
measured coordinates; it is not this denominator policy. Moreover 319's
actual config has `shrink_key=none`, `keep_key=effect`: changing the old
agreement function would not fix its active path. 327 changed a coordinate
estimand toward a median only on full-support coordinates. 228 fitted
gene-specific reliability gains; 054 added a CD4 union source. This test
does none of those. Their negative outcomes remain contrary evidence for
expecting a free gene-magnitude gain, not measurements of this construction.

## Fixed construction and controls

Use the exact 323 cap15 H1 composition, H1-excluding operator and unchanged
actual-count emitter. H1 union OFF, alpha .982, CD4 f .34/k 5, no mass fix.
The frozen 327 composer supplies its three alpha/global/renorm-scaled raw
source rows, target-and-gene observation masks and exact cap15 pre-PM row.

Recompute the cap weights from those scaled source rows in float64:
per donor, divide a target's full-axis L2 norm by the median among carried
targets; weight = max(1, ratio)^(-1.5). The per-donor scalar cancels in
this ratio mathematically. This is an explicit float64 reconstruction,
not a claim of bit-exact original float32 intermediate weights. Base
donor weights are uniform, and a known-answer fixture checks the uniform
weight limit. Missing target rows do not enter the norm medians.

Let c(t,g) be the sum of normalized cap weights from donors measuring
gene g, and n(t,g) their unweighted count. Let k(t) be carrying donors.
Only the operator's non-target coordinates with n>=2 may receive an
availability multiplier. One-donor coordinates remain unchanged before
norm restoration. Full-support coordinates have multiplier exactly one.

Four arms:

1. `cap15`: unchanged baseline.
2. `available`: multiplier min(4,1/c) on partial-support coordinates.
3. `count_only`: multiplier k/n, no cap-weight information.
4. `shuffled`: the available arm's multiplier vector cyclically permuted
   within each target's operator coordinates with n>=2, seed 3370909.
   This preserves the multiset of multipliers, not the resulting energy.

For every nonbaseline arm restore the complete operator-subvector L2 norm
to cap15, leaving all coordinates outside it bit-identical. Full-support
coordinates can change through that shared norm restoration; do not call
them bit-identical after it. If either subvector norm is zero, retain base.
Then apply production's existing PM correction, whitening and CD4
composition unchanged. Pre-PM norm equality is NOT delivered invariance.
No dose grid, outcome-based coordinate selection or alternate policy.

## Verification

Pure fixtures: all-observed identity, true-zero versus absent coordinate,
one-source handling, uniform-weight count-control identity, deterministic
derangement, cap bound, outside-axis identity and subvector norm identity.
Before scoring any new arm, baseline pre/PM/final hashes, emitted-count and
panel-total digests and all six metrics must reproduce 323 exactly at its
frozen tolerances. Freeze original composer/evaluator dependencies as well
as this adapter. All edits and outputs stay in 337.

The adapter reuses 327's construction and scoring functions in-process,
with only its arm constructor, output directory, names and input verifier
replaced. Remove the inherited median-dose metadata field from 337's
receipt and record this policy explicitly. No shared file is edited.

## Analysis and gates

Actual counts on all 85 reused H1 targets; all six official metrics.
10,000 paired target resamples with seed 20260909, fixed gallery/anchors/
models, panel MSE and global claim cap recomputed in each draw.
Keep raw MSE visible even when its local scaled score clamps to zero.

`available` earns further development only if:

- overall H1-local 95% lower bound >0 versus cap15;
- overall lower bound >0 versus shuffled multipliers;
- raw NMAE point improves versus cap15 and its upper bound <=+.005;
- count-only does not beat available's overall point estimate;
- all construction/scoring/identity gates pass.

Otherwise record no advancement and do not retune. All control comparisons
are conditional development diagnostics, not multiplicity-adjusted
discoveries. One reused context does not establish A/B/C transfer or
depth robustness; no slot is cleared by this result.

90-minute resource bound, one heavy process, maximum 22 GiB RSS with the
repository launch gate, no paid cloud or external calls. Stop on any
known-answer, source, operator, finite-value or scoring failure.
