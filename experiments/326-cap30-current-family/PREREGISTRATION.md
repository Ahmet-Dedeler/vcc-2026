# 326 — Current-family H1 cap-dose diagnostic

Registered September 9, 2026, before cap30 outcomes on this reference.
Coordinator: Codex. This is a local diagnostic of the existing, built
319d3 cap-at-exponent-minus-three construction, **not a novel model**.
Source/operator fitting and emitter follow experiment 323 exactly.

## Question

On the repaired H1 current-family composition, how does changing only the
norm-weighting exponent from -1.5 to -3.0 in cap mode change all six metrics?
Does this reference show a distinguishable cap dose response, and how does it
compare descriptively with 323's already-measured plain dose response?

Do not change 319's frozen board bins or interpret a local result as its
A/B/C outcome. 323 did not qualify this one-context panel as a standalone
small-arm selector. The unknown A/B/C cap30 score stays unknown.

## Construction and known answers

- Reuse 323's hash-pinned, H1-target-excluding operator and re-extracted
  complete-mean3 CD4 bank. H1 donors/union remain absent.
- Alpha .982, gamma .65, lambda 0, CD4 f .34/k 5, mass fix unset, production
  PM placement, and the exact rc_callset emission settings of 323.
- Prepare the three non-H1 production donors with the same functions as 323.
- Reproduce cap15's staged and final effect hashes before scoring cap30.
- Reproduce cap15's actual emitted count digest and all six metric values
  at 323's frozen tolerances. Never relax a tolerance to proceed.
- On A/B/C, compare cap30 pre-whitening/PM-corrected effect hashes with the
  build metadata of existing `submissions/319d3-cap30.vcc`. No A/B/C reference
  perturbation outcomes or new submission are accessed.
- Freeze the implementation and explicit input hashes before measurement.
  Nothing under experiment 323 may be edited or overwritten.

## Frozen arms and comparisons

1. Cap15 reproduction (known-answer arm).
2. Cap30 (-3.0, cap, PM default).
3. Cap30 target derangement, seed 20260909 (large-effect negative control).

Primary descriptive contrast: cap30 minus cap15.
Secondary descriptive interaction:
`(cap30 - cap15) - (plain30 - plain15)`, using 323's frozen plain arms.
This contrasts complete constructions; it does not isolate noise suppression
or a pure delivered-amplitude mechanism.

Compute actual emitted counts, all six raw metrics, H1-anchor-scaled overall
values, per-target rows, and 10,000 paired target bootstraps with seed
20260909. Recompute panel MSE and its global noise-credit cap on each draw
using 323's tested bootstrap implementation. Condition intervals on the fixed
gallery, local H1 anchors and models. Report nonfinite draws as failures.

Report the point sign and whether the paired 95% interval excludes zero,
regardless of direction. No outcome-dependent dose, seed, target filtering,
additional candidate or threshold. Both registered comparisons are descriptive,
not multiplicity-controlled discoveries. The shuffled arm must be worse in PDS
with the paired lower bound above zero for unshuffled minus shuffled.

## Bound and authority

One guarded heavy process, peak ≤22 GiB, scratch ≤5 GiB, implementation and
first measurement ≤90 minutes unless an explicit coordinator amendment records
the reason. Stop on hash/axis/known-answer failure before interpreting cap30.
Only own experiment outputs; no production changes, packaging, uploads,
reserve access, external model calls or paid resources by this script.

Expected deliverables: `RESULTS.md`, `work/results.json`, exact metrics/digests,
paired intervals, and a clear statement of what this reused H1 diagnostic does
and does not say about 319d3.
