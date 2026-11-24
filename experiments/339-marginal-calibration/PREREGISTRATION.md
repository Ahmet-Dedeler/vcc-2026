# 339 — zero-anchored donor marginal calibration

September 9, 2026. Source-direction development pilot, not a VCC score.
Freeze identities, code and raw-source hashes before new comparisons.

## Question and novelty

Can donor-specific, gene-wise nonlinear response-range calibration improve
transfer of the same target to a wholly held line?

006 ranked GENES WITHIN each perturbation and lost magnitude information.
327 replaced coordinate-wise donor means with medians. 330 fitted shared
per-gene scales; 337 normalized missing-source support. None estimates a
separate empirical response marginal for each donor and gene across
non-evaluation training targets. This is a different statistical
construction on existing data, not a new biological information source.
080/195's learned retrieval adapters and 286's cross-product operators
are not reopened. No parameter selection or post-result retuning.

## Fixed data and splits

Use precisely 334's six CRISPRi lines, 1,024 non-panel targets and 512
response coordinates, and its two complementary 512-target blocks.
Both H1 evaluation and A/B/C panel coordinates/targets remain excluded.
Load raw delta arrays and divide each row by its coordinate RMS exactly
as in 334, without across-target centering. Use no text vectors or labels.
This is reused development data; study and line remain confounded.

For each outer held line and target block, fit marginal curves on the
other five lines' TRAINING targets only. Query the other five lines'
responses for TEST targets. The held line never contributes a fit or
query source. Its held targets are read only for endpoints.

The diagonal comparator is exactly 334's inner-held-line estimator on
training targets (four donors per inner training recipient), applied to
the five-donor query mean. Reproduce its saved predictions and the source
mean on all twelve splits.
The successful 334 artifacts are under `recovery-work/`, not its original
coverage-failed `work/`. This path was corrected after an initial freeze
attempt stopped on a missing file, before an input manifest or any fit
existed. The metadata selection, hypothesis and gates were unchanged.

## Frozen calibration

For one gene, sort each donor's 512 training values into Q_j(r). Their
equal-weight mean at each empirical rank is the common quantile curve
Q_bar(r). Define T_j by piecewise linear interpolation of
Q_j(r) -> Q_bar(r). Equal donor values share the mean Q_bar across
their tied ranks. Outside the fitted value range, clamp to the endpoint.
If a donor curve has fewer than two distinct values, use identity.

Apply the ZERO-ANCHORED map T_j(x)-T_j(0), then average over donors.
The zero anchor prevents creating an effect from an all-zero input.
No query-target values fit or extend the marginal curves.

Arms:

1. `source_mean`: untouched mean; known-answer.
2. `diagonal`: 334's fixed learned diagonal; known-answer.
3. `marginal`: primary zero-anchored empirical marginal transport.
4. `gaussian`: zero-anchored Gaussian-marginal control, i.e.
   x_j * mean_j(sd_j)/sd_j; population SD on training targets only.
   For sd <=1e-12 use identity. No clipping.
5. `clipping`: x_j clipped to its OWN training range, minus clipped zero.
   Constant training curves use identity. This isolates tail clipping.
6. `curves_shuffled`: identical marginal construction with all donor
   curve banks cyclically deranged across gene identities, one fixed
   permutation, while the query coordinates stay fixed. This is a matched
   correspondence control, not a sampled null distribution.

All six arms predict every target and gene; no favorable subset. Curves
need no cross-donor training-target pairing: independent permutation of
each donor's training rows must leave them identical, as a known answer.

## Checks, endpoints and gates

Pure fixtures: zero anchor; identical marginal identity inside support;
analytic linear-scale transport; tie averaging; clipping-only equality
when all source marginal curves are identical; constant-curve fallback;
training-target permutation invariance; held/query noninterference;
shuffled curve identity mapping; finite outputs.

Reproduce the five historical raw-loader arrays used by 334 and every
source-mean/diagonal prediction. Poison the held line and require identical
calibration and comparator predictions. Verify each held target has one
out-of-fold prediction. Never fit a calibration on held target rows.

Report 334's per-target own cosine, full 1,024-target gallery retrieval,
decoy cosine, unit-RMS absolute and squared error. They are NOT VCC raw
MSE/NMAE. Use 10,000 paired target and crossed line/target bootstraps,
fixed seed 3390909, conditional on fitted curves. No training bootstrap.

Advancement requires ALL:

- marginal retrieval improvement >=.005 over source_mean;
- both retrieval 95% lower bounds >0 versus source_mean, diagonal,
  gaussian, clipping and curves_shuffled;
- own-cosine point loss <=.005 versus source_mean AND diagonal;
- squared-error point no worse than source_mean;
- retrieval versus source_mean positive on >=5/6 lines, in every
  leave-one-line-out average, and separately in 334's frozen lower/higher
  median-cell-number line halves;
- every integrity and known-answer check passes.

All six competition members remain unpriced and potentially movable.
A pass earns only a separate current-family count-stage registration.
A failure closes this fixed empirical-marginal construction, not all
possible nonlinear calibrations. No sweep or outcome-selected blend.

## Budget and authority

One guarded CPU job, six threads, 8 GiB ceiling, 45-minute total budget,
repository memory launch gate. No H1 response read, reserve, new data
download, paid API/cloud, worker, production edit, package or submission.
