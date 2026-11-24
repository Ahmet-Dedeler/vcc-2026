# 129 — broad Pan response-neighborhood transfer

Status: **source-only protocol and implementation in preparation**

Date: **2026-08-30**

## Correction to the original brief

GSE295214 is the Pan2025 H1 source already used by Experiments 117 and 118.
Those experiments close direct same-target Pan residuals, their fractions,
their target masks, and their Pan-derived confidence rankings. Experiment 129
does not rerun any of them.

This experiment instead asks whether the 2,354 broad non-panel Pan rows define
a useful **target-similarity space** for transferring a correction learned on
other public source contexts.

For a query target, Pan is used only to retrieve biologically similar training
targets. The query target's Pan effect is never added to its prediction.
Training labels are cross-context tangent residuals from RPE1 and Nadig
Jurkat. RPE1 is predicted only from Jurkat labels and Jurkat only from RPE1
labels. K562-GWPS and both Orion lines are predictors and whitening sources,
never outcome surfaces.

## Fresh target boundary

The eligible target pool is the exact-symbol intersection of Pan, RPE1,
Nadig Jurkat, K562-GWPS, Orion HCT116, Orion HEK293T, and the pinned STRING
embedding, excluding:

- all 300 VCC targets;
- every Nadig target whose outcome was consumed by Experiment 107; and
- every source target consumed by Experiment 123.

This leaves exactly 155 targets. A deterministic hash order assigns the first
100 to source-only selection and the final 55 to an untouched source test.
No held HCT116, HEK293T, KOLF, Pan, DLD1, HCASMC, or challenge outcome is
opened.

## Learner

The source-only incumbent is an a14/arm-D-style geometry surrogate:

- K562-GWPS, Orion HCT116, and Orion HEK293T have equal source slots;
- PerturbMatch replaces 50% of the K562 slot where both support the target;
- source effects are donor-centered and equalized by source median row norm;
- a source-only signal-covariance whitener uses `gamma=0.65`;
- the operator excludes all VCC response coordinates and all 155 experiment
  target coordinates;
- CD4 is an exact no-op because these are non-panel perturbation targets.

For each training target, the teacher correction is the tangent component of
the opposite context's unit truth direction relative to the incumbent unit
direction. Pan unit-effect cosine retrieves the `k` nearest training targets.
Their teacher tangent residuals are averaged with equal weight, projected
orthogonally to the query incumbent, and applied as a fixed five-degree
rotation. The only selected hyperparameter is:

```text
k in {5, 15, 30}
```

Identity is always present and wins exact ties. Unsupported or degenerate
queries are bit-exact incumbent fallback.

## Validation

- Five deterministic whole-target folds produce out-of-fold predictions on
  all 100 selection targets.
- Candidate selection requires nonnegative complete-gallery raw uncentered
  PDS and signed-credit deltas in both RPE1 and Jurkat.
- The selected candidate is sealed before any of the 55 source-test rows are
  opened.
- Test PDS is computed over each complete 55-target gallery with exact
  midranks.
- Directional signed credit is
  `sum(sign(truth) * prediction) / sum(abs(prediction))` per target.
- A 20,000-draw reciprocal-STRING-cluster bootstrap and 199 complete-procedure
  target-label permutations are frozen before source effects are loaded.

A source-only pass authorizes at most a separately reviewed held evaluation.
It does not authorize production, counts, a challenge build, or submission.

