# Experiment 141 — cross-context Pan residual calibration

Status: **preregistered development experiment**

Date: **2026-08-30**

## Question

Experiment 117 showed that the same target's Pan2025 H1 residual contains a
small but real cross-context signal (`+0.00119484` raw PDS), while a direct
fixed-fraction blend was too weak to promote. Can a low-capacity
cross-context calibration retain the transferable coordinates of that
residual and suppress H1-specific coordinates?

## Candidate family

For destination biological family `h`, permitted teacher family `c`, and
target `t`, let:

```text
m[h,t]   = the frozen Experiment 140 production-shaped non-PerturbMatch surrogate
q[h,t]   = norm-matched Pan row in h's fixed operator geometry - m[h,t]
r[c→h,t] = teacher-context truth - m[h,t]
```

For each destination, fit only inside that one operator geometry:

```text
b[g] = sum(q * r) / sum(q^2)
b0   = sum(q * r) / sum(q^2), pooled over all coordinates
a[g; rho] = clip((1-rho) * b0 + rho * b[g], 0, 1)
candidate = m + a * q
```

The frozen grid is:

- `identity`
- `scalar_control` (`rho=0`; diagnostic only, never promotable)
- `diag_rho050` (**the sole primary candidate**)
- `diag_rho025` and `diag_rho075` (descriptive diagnostics only)

There is no second dose, target router, target subset, confidence threshold,
or context-specific production choice.

## Isolation

- Reuse only Experiment 132's already-consumed 613-target development split.
- Keep its 410-target source test sealed.
- Treat K562-GWPS and K562-essential as one biological family.
- When a biological family is scored, exclude it from every calibration fit.
- K562 is also excluded as a calibration teacher for RPE1 and Jurkat because
  the non-K562 surrogate contains a K562 donor slot. Thus K562-held learns
  from RPE1+Jurkat, RPE1-held learns from Jurkat, and Jurkat-held learns from
  RPE1. Each fitted coordinate map stays inside one fixed operator geometry.
- When a target fold is scored, exclude those targets from every calibration
  fit in every context.
- Pan is an input only, never a teacher or outcome.
- Pan never enters a whitening fit. Experiment 140's frozen leave-family-out
  operators and baselines are reused byte-for-byte.
- Pan centering excludes all 1,023 Experiment 132 gallery targets, including
  the sealed 410. The sealed rows are used as names only and their effect
  values are never read.
- Remove all 300 VCC target coordinates and all 1,023 gallery-target
  coordinates by reusing Experiment 132's 4,992-gene axis.
- Experiment 140's baseline intentionally omits a15's PerturbMatch K562-slot
  replacement. This is a development surrogate, not an exact-a15 comparator.
  A positive result must first reproduce on an exact-a15 613-target comparator
  before any independent confirmation can be considered.

## Selection

Five frozen target folds are reused from Experiment 140.

For each target fold, fit each fixed candidate separately for each destination
family using the frozen teacher policy and all other target folds, then
evaluate that fold against the complete 613-target gallery. No candidate is
selected from development outcomes: `diag_rho050` is the sole primary and the
other arms are controls or diagnostics.

## Gates

Development may progress to expensive nulls and bootstrap only if the fixed
cross-fitted `diag_rho050` primary satisfies:

- PDS delta nonnegative in every biological family;
- signed-credit delta nonnegative in every biological family;
- macro PDS delta at least `+0.003`;
- macro signed-credit delta at least `+0.003`;
- `diag_rho050` beats `scalar_control` by at least
  `+0.001` macro PDS.

If the raw gate passes, the already-implemented frozen follow-on requires:

- 49 cross-fitted `diag_rho050` Pan target-label nulls;
- 49 cross-fitted `diag_rho050` Pan coordinate-correspondence nulls;
- 20,000 target-block bootstrap draws;
- both observed endpoints above both null p95 values;
- both bootstrap lower 95% bounds positive;
- positive-gain concentration at most 0.5 by family and target block.

Both null mappings are shared across families and folds and are applied to
the transformed Pan rows before destination-row norm matching. Incumbents and
outcomes remain fixed.

Only a separate review and an exact-a15 surrogate follow-up could authorize
the sealed 410-target confirmation.
No result from this experiment authorizes counts, `cell-eval2`, a challenge
build, or a submission.

## Distinctness

- Experiment 117 directly interpolated toward Pan with fixed fractions.
  Experiment 141 learns a gene-coordinate reliability map; its scalar
  equivalent is only a non-promotable control.
- Experiment 129 used Pan only as cross-target neighborhood geometry.
  Experiment 141 uses the query target's own Pan discrepancy and no neighbors.
- Experiment 132 replaced the response with a DepMap-feature Ridge prediction.
  Experiment 141 preserves the incumbent and uses no DepMap response model.
- Experiment 136 scaled or routed the rejected Ridge residual. Experiment 141
  contains no Ridge output or margin router.
- Experiment 096 found that an unrestricted per-gene model-to-truth slope was
  largely cell-line-specific. That is a strong negative prior. Experiment 141
  differs by fitting a sign-constrained reliability map for an external
  same-target residual around an unchanged incumbent, and it must beat its
  scalar ancestry control by `+0.001` PDS to demonstrate that distinction.
