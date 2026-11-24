# 109 — Gauge-invariant global guide nuisance

## Status

**SOURCE-ONLY NO-GO.** The repaired guarded run completed on 2026-08-30 and
selected exact identity (`eta=0`). Untouched source-test cosine and normalized
MSE deltas were both `0.0`; the bootstrap lower bound and both permutation-null
p95 values were also `0.0`. Every positive-source gate therefore failed and
held-context evaluation remained disabled. See `RESULTS.md`.

The scientific design was preregistered before any experiment-109 expression
extraction, fitting, held-context scoring, count generation, or challenge
submission. An independent adversarial review found deployment/objective and
authorization defects in the first scaffold; the protocol, validators, and
source orchestrator were hardened before execution.

The first guarded cache-generation attempt on 2026-08-30 streamed the complete
K562 source matrix, then failed the compact-reconstruction gate before cache
publication or fitting. The extractor had incorrectly allowed the per-guide
>=20-cell threshold to redefine the target-level compact estimand. This
explicit pre-fit repair preserves guide filtering for guide effects but
reconstructs compact target pseudobulks from all cells assigned to targets with
>=20 total cells. The `2e-5` gate and every scientific parameter are unchanged;
no held outcome was opened.

The repaired rerun reconstructed all 9,787 compact targets with maximum
absolute error `4.76837158203125e-07`, then completed the frozen source test,
20,000-target bootstrap, and 100+100 source nulls. It rejected the family
before any held-context outcome access.

This directory contains:

- the frozen scientific protocol and all candidate grids/gates;
- a pure NumPy stable-subspace model;
- fit/freeze helpers and a guarded source-only orchestrator that operate only
  on experiment-083 compact K562 arrays and frozen source artifacts;
- strict schemas for source inputs, source-test outputs, singleton
  corrections, and held-label nulls;
- adversarial synthetic tests for the structural and authorization claims that
  make singleton deployment identifiable.

The source-only orchestrator can generate missing compact experiment-083
caches from the pinned local K562 source or verify existing caches. It does
not open held-context truth, launch cloud work, generate counts, or submit
anything. Held evaluation remains deliberately unreachable in `evaluate.py`.

## Question

Can guide-sequence features explain a reproducible component of the measured
K562 guide response, allowing that component to be removed from singleton
source rows without changing the source-wide mean?

Experiment 083 tested a multi-guide target estimator. Its within-target
contrasts cannot identify the absolute origin of a nuisance function:

```text
guide_effect(t, g) = target_effect(t) + nuisance(features(t, g)) + error
```

For any response vector `a`, replacing `nuisance(x)` by
`nuisance(x) + a` and the target effect by `target_effect - a` leaves every
within-target contrast unchanged. Directly subtracting experiment 083's
absolute nuisance prediction from singleton targets would therefore be
invalid.

Experiment 109 tests only the identifiable, gauge-invariant estimand:

```text
correction(t) =
    -eta * (
        nuisance(P_R x_t)
        - mean over every sequence-supported singleton source target
    )
```

`P_R` is the stable row space of within-target feature contrasts. The
corrections are forced to sum to zero, to floating-point audit tolerance,
over the complete singleton centering population. Unsupported and
multi-guide targets receive exact zero.

## Frozen source audit

The metadata/sequence-only design audit found:

- 506 identifiable non-panel multi-guide development targets;
- deterministic target split: 303 fit / 101 tune / 102 untouched test;
- 229 guide features, 225 varying on fit;
- 188 stable contrast-supported dimensions at
  `singular_value / largest_singular_value > 1e-6`;
- 223 VCC-panel targets with sequence support;
- 212 with both sequence and PerturbMatch support:
  206 singleton and 6 multi-guide;
- 6,689 sequence-supported targets with exactly one retained source construct
  in the complete K562 centering population. The earlier source-only count
  6,692 included `CDK11B`, `DST`, and `TAL1`, which each have two retained
  constructs but only one sequence-annotated construct; they are now excluded;
- median panel-singleton feature-norm fraction inside the stable contrast
  subspace: 0.995.

These are preregistered expectations, not results. A future real-data run must
reproduce the pinned source hashes and audit counts or stop.

## Leakage-safe fit

1. Reuse experiment 083's compact guide-effect and sequence-feature formats.
   The extractor and feature builder are not copied here. Guide effects retain
   only >=20-cell guide categories, while compact-target integrity
   reconstruction independently includes all cells for targets with >=20 total
   cells.
2. Exclude all 300 VCC target coordinates and all 506 development-target
   coordinates from the learned response.
3. Hash-split identifiable non-panel multi-guide targets 60/20/20 by target.
   No guide from a target may cross partitions.
4. Standardize features using fit-target rows only.
5. Center features and effects within target.
6. SVD the centered fit-feature matrix and retain only directions satisfying
   `sigma / sigma_max > 1e-6`.
7. Learn response bases of rank `{8,16,32}` from centered fit effects.
8. Ridge-fit the stable feature scores to response-basis scores with
   `lambda in {0.1,1,10,100}`.
9. Tune `eta in {0,0.25,0.5,0.75,1}` on tune targets only.
   Guide-cell weight power is fixed at **0**. It is part of the validation
   definition, not a deployed parameter, and therefore may not be selected.
10. Validate by holding out every guide in each target as a pseudo-singleton.
    Form the raw held-versus-sibling residual, clean the held and sibling rows
    through the deployment formula, and score the nuisance component removed
    by cleaning against the raw residual. This is translation-invariant and
    gives exactly zero credit to a target-constant nuisance. Metrics are first
    averaged within target so targets with more constructs cannot dominate.
    Raw-guide reconstruction is retained only as a diagnostic and cannot
    select a model.
11. Refit the selected configuration on fit+tune and freeze that exact model as
    `source-test-model.npz`.
12. Test that model once on the untouched test targets and freeze every
    per-guide raw residual, recovered-nuisance vector, label, and guide row in
    `source-test-holdouts.npz`. Metrics and bootstrap samples are recomputed
    from those arrays; JSON metric claims are never trusted.
13. Only a passing source test may be refit on all 506 development targets as
    the distinct deployment `selected-model.npz`, before held truth is opened.

The implementation imports only pure numerical helpers from experiment 083:
deterministic splitting, robust centers, randomized response bases, cosine/NMSE
metrics, and stable seeds. Experiment 083's absolute singleton adjustment path
is intentionally not reused.

## Deployment contract

- Center the nuisance prediction over **all** 6,689 sequence-supported
  singleton K562 targets, never only the VCC panel. The sorted target axis,
  count, and SHA-256 are mandatory inputs; subsets, duplicates, and reordered
  axes fail.
- Apply a correction only to a sequence-supported singleton K562 row.
- Keep unsupported and multi-guide rows exactly unchanged.
- Force the raw K562 correction sum over the complete singleton support back to
  zero after floating-point arithmetic. This invariant is defined **before**
  target-by-gene PerturbMatch retention. No post-retention zero-mean claim is
  made.
- Keep all excluded response coordinates exactly zero.
- Retain K562's existing 50% share only on the target-by-gene intersection
  where PerturbMatch supports both the target and the gene. PerturbMatch-
  unsupported targets and genes retain the full K562 contribution.
- Leave PerturbMatch unchanged.
- Require correction rows to arrive in the exact canonical 300-target order.
  A caller with a different complete order must explicitly request and test a
  remap to canonical order; set equality alone is not accepted.
- Then use the unchanged donor average, frozen whitening, norm preservation,
  and CD4 residual from the production arm-D path.

## Preregistered nulls

- exact identity (`eta=0`);
- within-target guide-feature derangement, with every eligible guide row moved
  in every replicate, repeating complete selection;
- split-local cross-target guide-row derangement. Strata use only identifiable
  guide count and sequence-feature bundle norm; they never use expression,
  effect estimates, held-test metrics, or cell-count outcomes. Strata are
  deterministically coarsened from guide-count plus split-local feature-norm
  decile, to guide count, to the complete split until every eligible target
  and guide row can move. Fit, tune, and untouched-test rows never mix;
- direction-preserving shrink control;
- frozen singleton-correction target-label permutation in held-context
  evaluation. The 100 permutations use seed `1092027`, are exact bijections,
  keep unsupported/multiguide rows fixed, and permute only within the
  PerturbMatch-supported and PerturbMatch-unsupported singleton strata;
- structural gauge shift: adding a finite constant response vector bounded by
  `1e6` in absolute value to every nuisance prediction must leave every
  deployed correction unchanged within the frozen scale-aware float64
  roundoff tolerance;
- every source-null row assignment is frozen in `source-nulls.npz`.
  Authorization regenerates each assignment from its seed and split, verifies
  every eligible test target moved, checks the permuted-feature hash, requires
  the complete candidate grid for every replicate, and independently reruns
  deterministic candidate selection. No `complete=true` Boolean is trusted;
- complete candidate-grid selection, target-cluster bootstrap, source-null
  derangements, and direction-preserving shrink gates are implemented as pure
  helpers;
- `source-input-manifest.json` binds every experiment-083 compact output, the
  extractor and feature-builder sources, a shared guide-axis hash, the
  canonical gene axis, and feature eligibility defined by normalized valid
  protospacer and seed sequences. The experiment-109 feature cache must carry
  a byte-identical copy of the extractor's `guide_ids.txt`; a feature cache
  with only a matching row count is rejected;
- the source runner copies every authenticated compact-cache member into
  `source-cache/` inside the staging bundle and fits only from that snapshot.
  Original and copied members are hashed before array loading and rehashed
  after all source computation. Atomic publication is refused if either side
  changes. `source-cache-snapshot.json` makes the frozen arrays and provenance
  inseparable;
- `memory-guard.json` is a validated receipt from the pinned repository memory
  guard. Monitoring starts before source preflight/cache generation and ends
  only after post-compute cache rehashing and receipt freeze;
- the panel-singleton feature matrix used by the projection diagnostic is the
  exact ordered subset selected from the frozen full singleton matrix; its
  shape, dtype, content hash, and source indices are rechecked;
- `source-test-model.npz` has a strict role/schema and must authenticate the
  exact fit+tune target axis;
- `singleton-corrections.npz` carries canonical singleton features,
  development features, panel mappings, response corrections, final
  source-slot corrections, and deterministic stable-score extrapolation
  diagnostics. Every correction is recomputed during authorization;
- a held evaluation cannot be authorized unless the source test and every
  individual source gate pass, the exact artifact set exists, all source and
  artifact hashes match, and the full centering/response axes authenticate.

## Promotion gates

All gates in `preregistered_protocol.json` are mandatory. In summary:

- untouched source-test cosine delta positive;
- untouched source-test normalized-MSE delta nonpositive;
- target-cluster bootstrap lower 95% bound above zero;
- source result above both fitted-null p95 values;
- exact identity at `eta=0`;
- zero correction outside supported singleton rows and excluded coordinates;
- zero-mean correction over the full singleton centering population;
- feature-nullspace perturbations change nuisance predictions by at most
  `1e-10`;
- exact reconstruction of the current production baseline;
- nonnegative H1-test and HEK293T PDS deltas, with HCT116 retained as a
  PDS-only safety diagnostic;
- nonnegative H1-test and HEK293T directional-credit deltas;
- equal-donor macro raw-PDS delta at least `+0.003`;
- singleton-only held-panel raw-PDS delta at least `+0.003`;
- observed held gain above the frozen correction-label-null p95;
- unchanged target coverage;
- KOLF nonnegative if reported, labelled as a previously inspected diagnostic.

No count-level run is allowed until every effect-level gate passes. No
submission is authorized by this experiment.

## Source-only CLI

Read-only preflight; this does not generate caches or fit:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy,h5py python \
  experiments/109-global-guide-nuisance/source_only.py \
  --preflight --cache-mode generate --allow-seed-download
```

The guarded execution command is intentionally separate and heavy. It
generates only absent caches, repeats the complete frozen candidate selection
inside both null families, monitors the complete process tree under the pinned
repository memory guard, and publishes the final source bundle by atomic
directory rename:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy,h5py python \
  experiments/109-global-guide-nuisance/source_only.py \
  --run-source --cache-mode generate --allow-seed-download --guarded
```

Partial caches, changed source/dependency hashes, existing output directories,
changed cache members during computation, invalid memory receipts, and failed
provenance checks stop publication. A failed source test can
publish a `source_rejected` evidence bundle, but never a deployment model,
freeze manifest, or held authorization.

The 2026-08-30 pre-fit extraction failure published neither a guide cache nor a
source bundle. The repaired extractor is pinned by its new SHA-256 in the
amended protocol; rerunning the guarded command still requires explicit human
review and preserves the original `2e-5` reconstruction threshold.

Experiment 107 is not imported, read, or awaited by this command. A failure or
absence of experiment 107 does not affect experiment 109's source-only path;
the guarded command depends only on the pinned experiment-083 source/cache
contracts and the files listed by its own preflight.

## Lightweight checks

```bash
PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy python \
  experiments/109-global-guide-nuisance/test_model.py

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy python \
  experiments/109-global-guide-nuisance/fit_freeze.py --check-contract

PYTHONDONTWRITEBYTECODE=1 uv run --no-project --python 3.12 \
  --with numpy python \
  experiments/109-global-guide-nuisance/evaluate.py --check-contract
```

HCT116 is PDS-only because its directional/DE labels are depth-confounded.
H1-test and HEK293T are the directional surfaces; both are
hypothesis-specific rather than pristine holdouts.
