# 327 — Bounded gene-coordinate median correction on the current family

September 9, 2026. Freeze before constructing or scoring these arms.
Coordinator: Codex. No external model call or new source data.

## Hypothesis and novelty

Whole-row cap weighting can suppress an oversized donor response but cannot
remove a donor's isolated gene-coordinate outlier without changing its other
genes. A coordinate median can change which downstream genes move, including
signs where the source mean is driven across zero by one donor. Test a bounded
25% interpolation toward that statistic, not a new scalar confidence model.

This adds **no new biological information** and assumes neither that the three
donors are independent nor that their majority is correct. It changes the
estimand from a cap-weighted mean toward an unweighted coordinate median after
the existing per-donor scale normalization. Robust pooling is not claimed as
a new statistical invention. It is a distinct construction not found in the
reviewed experiments.

History checked: 097 (whole-row weights), 132 (DepMap replacement), 161
(context residual and power-law magnitude control), 195/080 (learned retrieval
operator), 224/228 (per-gene diagonal and delivered NMAE cost), 234/101
(learned/selective sign flipping), 268/269 (context residuals and router
leakage), 271/286 (replicate-Wiener/cross-product estimators), plus 160/174
(control geometry). The median uses neither a noise estimate, a fitted
diagonal, a target lookup, nor a learned flip key. Three-source gene support
is essential; missing coordinates are **not zero-valued observations**.

## Fixed construction

Baseline is exactly 323's cap15 H1 composition. Reuse its frozen source
preparation, H1-excluding operator, CD4 complete-mean3 bank, and emitter:
alpha .982, gamma .65, CD4 f .34/k 5, H1 union OFF, mass fix unset.

For each target, use the current three main donors only. A coordinate is
eligible only if all three carry that target AND measure that gene AND the
gene belongs to the H1 operator's non-target axis. Else leave it unchanged.

Let `b` be cap15's **pre-PerturbMatch** row. Let `m` be the coordinate median
of the three source rows after the same `_renorm`, `_global`, and alpha
scalars used by `MultiDonorModel.effect`. Beta must be exactly zero.

On eligible coordinates, set `q = .75*b + .25*m`. Rescale this eligible
subvector to the eligible subvector's baseline L2 norm. All ineligible
coordinates remain bit-identical. If either subvector norm is zero, retain
the baseline subvector. This identity is only **pre-PM**; make no claim of
delivered-amplitude invariance.

Then use production's existing PM correction, whitening, CD4 residual
composition and actual-count emission unchanged.

## Frozen arms and controls

1. `cap15`: baseline reproduction, must pass before new scoring.
2. `median25`: construction above, the only candidate.
3. `magnitude25`: median25's exact pre-PM magnitudes, but baseline signs on
   every nonzero baseline coordinate. At a baseline zero retain median25.
   This isolates pre-PM sign differences, not delivered sign differences.
4. `power25`: the 161-style magnitude-only comparator on the same eligible
   support. Norm-match `sign(b)*b^2` to the baseline subvector, interpolate
   25%, then restore that subvector norm. No donor median.
5. `shuffled25`: independently permute median25's correction `q-b` within each
   row's eligible gene coordinates, seed 20260909; add to baseline and restore
   eligible norm. This preserves pre-restoration correction energy and support,
   not final delivered perturbation size. Measure the latter explicitly.

No dose grid, adaptive sign region, target filtering or model selection on H1.
No additional candidate is authorized by this registration.

## Known answers and frozen inputs

Unit tests: exact zero-dose identity; missing-source/coordinate fallback;
coordinate-median fixture including a sign disagreement; eligible-subvector
norm identity; outside-support identity; sign-control magnitude identity;
power-control original signs; deterministic shuffled correction.

The baseline must reproduce 323's pre/PM/final array hashes, emitted-count and
panel-total digests and all six raw metrics at the unchanged 323 tolerances.
Freeze protocol, implementation, tests and the transitive 323 input hashes
in git before construction. No writes into 323 or production modules.

## Readouts, uncertainty and decision

Actual counts and all six official metrics on the 85-target H1 panel.
Primary: median25 minus cap15 on H1-local scaled overall and raw NMAE.
Controls: median25 minus magnitude25, power25 and shuffled25.
Use 10,000 paired target draws, seed 20260909, conditioning on the frozen
gallery/anchors/models. Recompute panel MSE and global claim cap per draw.
Record raw MSE separately wherever the local score clamps it.

Report pre-PM, post-composition and emitted pseudobulk sign changes and norm
ratios, and fraction of coordinates/targets eligible. Do not infer a delivered
identity from a pre-PM identity.

This one reused context cannot clear a submission. It earns additional
multi-context **development only** if the primary overall 95% lower bound
is positive, raw NMAE upper bound is negative, and raw NMAE versus BOTH
magnitude25 and shuffled25 has an upper bound below zero. Power25 must not
beat median25's overall point estimate. Otherwise record the scoped negative
or unresolved result and do not tune until it passes. Control comparisons
are descriptive, not multiplicity-controlled biological discoveries.

## Resource and authority bound

First construction/measurement within 90 minutes; ≤22 GiB process-tree RSS,
≤5 GiB scratch; one heavy process at a time with repository memory gates.
Stop on construction/hash/scorer failure. No reserve, external downloads,
production changes, package build, submission, commits outside own work,
push, public disclosure or paid resources.
