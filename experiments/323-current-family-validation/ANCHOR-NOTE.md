# H1 anchor reconstruction — September 9, 2026

The missing H1 anchor artifacts can be rebuilt locally after the six-metric
known-answer gate passes. This fulfills the original normalization requirement;
it is not a change to a candidate, target cohort, or qualification threshold.

Use cell-eval2's `generic_response_profile` and `build_baseline_prediction`
with dispersed emission and seed 0 for the baseline. Its fractional counts are
an official baseline-only allowance, never a submission-format allowance.

For the replicate anchor, reuse the package's `_derive_seeds(0, 5)` and
`_disjoint_halves`, with disjoint controls and each half's own DE reference.
Use the arithmetic mean of the five metric aggregates. NMAE uses the package's
`_nmae_ref_from_tables` on the **full-reference** gate and denominator, not the
gate selected by a half. Use the exact-normalization fast rank kernel validated
in amendment 3 and the official downstream metric functions.

Before H1 anchor measurement, a synthetic three-target/40-gene fixture must
reproduce both the ordinary full baseline scoring pipeline and
`compute_replicate_anchor` on all six metrics under the existing known-answer
tolerances. Hash-pin the code, upstream constructors, fixture result and
historical gate result, then commit.

These are **H1-local reconstructed anchors**, not the competition's A/B/C
anchors or a server-attested real bundle. MSE is clamped to [0,1] and scaled
NMAE floored at -6 in the local composite. Overall comparisons concern one
recipient and remain retrospective.
