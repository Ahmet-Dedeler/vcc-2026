# Implementation amendment 1 — September 9, 2026

The first reference pass stopped while reading `var/_index`, before reading
`X` or computing reference statistics. The H1 file uses a nullable-string
`values` group for its gene index; the initial helper handled categorical
groups and plain strings only.

Add nullable-string handling with a fail-closed missing-label check. A synthetic
HDF5 test covers all three encodings and missing labels. No cohort, candidate,
seed, metric definition, tolerance or scientific gate changes. The original
evaluator hash remains in commit `4493f327`; the evaluation input manifest is
updated in the amendment commit to bind the repaired implementation.

Failure evidence: `work/reference.log` (first launch, retained).
