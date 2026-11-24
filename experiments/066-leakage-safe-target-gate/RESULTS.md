# 066 — Leakage-safe per-perturbation selector over the 064 joint arm

**Date:** 2026-08-23  
**Decision:** **REJECT for shipping. Keep experiment-064 arm D fixed.**

## Question

Can a low-capacity gate use only information available at challenge inference
time to select, per perturbation, among experiment 064's frozen arms?

- `A`: incumbent production union-whitened baseline
- `B`: `A` plus the fixed 50:50 incumbent/PerturbMatch K562 slot
- `C`: `A` plus the fixed post-whitening CD4 residual
- `D`: both changes, the experiment-064 winner

The primary estimand is effect-level PDS improvement over always using `D`.
There was no count emission, `cell-eval2`, challenge prediction, or submission.

## Bottom line

The nested gate found a small, repeatable-looking **PM exclusion** signal:

| outer train → test | selected gate | switched arm | switched targets | fixed D PDS | gated PDS | gain |
|---|---|---|---:|---:|---:|---:|
| HEK293T → HCT116 | `pm_shift_size <= 0.0472149253` | `C` | 85 / 221 | 0.825092555 | 0.826532291 | **+0.001439737** |
| HCT116 → HEK293T | `pm_slot_companion_advantage >= 0.00351278763` | `C` | 108 / 221 | 0.800658165 | 0.801851090 | **+0.001192925** |
| **macro** | | | | **0.812875360** | **0.814191691** | **+0.001316331** |

Both selected models were one-feature decision stumps, and both only chose
between `C` and `D`: retain CD4, but omit PerturbMatch for a subset of targets.
No CD4-removal or four-arm rule survived nested selection.

This is not enough to ship:

1. The gain is below the predeclared `+0.003` practical floor.
2. The 20,000-replicate target-cluster bootstrap interval crosses zero:
   `[-0.001912793, +0.004781983]`.
3. The selected feature and direction are not stable across donor directions.
   The HEK-trained rule uses PM shift size; the HCT-trained rule drops PM when
   PM improves agreement with the companion donor, which is not a stable
   mechanistic rule.
4. There are only two valid outer donors. A final gate trained on both would
   have no untouched donor on which to validate its model-family choice.

The signal is worth retaining as a hypothesis for a genuinely independent
third somatic holdout, but it should not replace fixed arm `D`.

## Oracle headroom

Oracle routing has substantial in-fold headroom:

| oracle | HCT116 gain vs D | HEK293T gain vs D | macro gain |
|---|---:|---:|---:|
| choose `D` or `B` | +0.005861785 | +0.006787330 | +0.006324558 |
| choose `D` or `C` | +0.015240642 | +0.014376800 | **+0.014808721** |
| choose `D` or `A` | +0.015528589 | +0.017174002 | +0.016351296 |
| choose best of `A/B/C/D` | +0.021287536 | +0.021143562 | **+0.021215549** |

The nested gate captures only:

- **8.89%** of the `D`-versus-`C` oracle headroom;
- **6.20%** of the full four-arm oracle headroom.

The gap is not lack of oracle opportunity. It is lack of transferable target
ranking. Even a 221-parameter same-target lookup—choose each target's best arm
on the other donor—loses:

| train → test | gain vs D |
|---|---:|
| HEK293T → HCT116 | −0.006993007 |
| HCT116 → HEK293T | −0.006478815 |
| **macro** | **−0.006735911** |

Cross-donor target-level stability is weak:

| quantity | Pearson | Spearman | same-sign fraction |
|---|---:|---:|---:|
| `A − D` | +0.0025 | +0.0132 | 0.5520 |
| `B − D` | −0.2855 | −0.0795 | 0.5158 |
| `C − D` | +0.0921 | +0.0484 | 0.5566 |
| `D − A` | +0.0025 | +0.0132 | 0.5520 |

## Validation protocol

- Outer folds:
  - train only on HCT116 outcomes, test only on HEK293T;
  - train only on HEK293T outcomes, test only on HCT116.
- Inner selection: deterministic five-fold target CV inside the outer training
  donor.
- Candidate class: 126 predeclared models:
  - always use one fixed arm;
  - one-feature decision stumps;
  - ridge gain models over fixed PM, CD4, context, and joint feature families.
- Maximum ridge capacity: 14 inputs and three predicted alternative-arm gains.
- Feature scaling, thresholds, coefficients, margins, and candidate selection
  were fit from inner-training rows only.
- Experiment-064 outer-test outcomes were not used for tuning.

## Inference-time features

The 14 features use a 6,102-gene common public-source axis with all 300 panel
gene coordinates removed:

- K562/companion-donor agreement;
- incumbent-K562/PerturbMatch agreement and relative geometry;
- PM-slot agreement improvement against the non-held companion donor;
- CD4/source and cross-condition agreement;
- CD4 response norm, breadth, and QC;
- target expression in the held donor's unperturbed control and its relation
  to source controls.

For HCT116 evaluation, HCT116 perturbation effects were excluded and only its
control target expression was used. For HEK293T evaluation, HEK293T
perturbation effects were excluded analogously. The held control is a valid
challenge input.

Target identity was a join key only and was never encoded. PDS columns, held
truth effects, experiment-058 per-target PDS, and experiment-062 per-target PDS
were not features.

## Permutation and null controls

Despite the small effect, the selected assignments were better than random
routing:

| control | replicates | null mean | null p95 | one-sided p |
|---|---:|---:|---:|---:|
| rerun complete nested selection after jointly permuting training-target `A/B/C` gains | 500 | −0.002002365 | +0.000380502 | **0.011976** |
| permute selected arm assignments within each outer test fold | 20,000 | −0.002037096 | +0.000462773 | **0.011849** |

These controls say the source-side PM geometry contains some target-routing
signal. They do not overcome the bootstrap uncertainty, small effect size,
unstable selected rule, or two-donor model-selection limit.

## Compact-output audit

The earlier compact outputs reproduce their corresponding 064 main effects,
which validates target alignment but does not make them admissible features:

| source | held | comparison | Pearson | Spearman |
|---|---|---|---:|---:|
| exp058 | HCT116 | PM f0.50 gain vs 064 `B − A` | 0.84599 | 0.76974 |
| exp058 | HEK293T | PM f0.50 gain vs 064 `B − A` | 0.89171 | 0.79001 |
| exp062 | HCT116 | CD4 f0.05 gain vs 064 `C − A` | 0.99754 | 0.99190 |
| exp062 | HEK293T | CD4 f0.05 gain vs 064 `C − A` | 0.99724 | 0.98459 |

Experiment 012 motivated source agreement as an inference-safe feature.
Experiment 050 motivated CD4 breadth/QC. Experiment 051's failed adaptive
breadth gate was respected: breadth was only one candidate feature, not a
privileged routing rule. Experiments 058 and 062 supplied diagnostic compact
per-target outputs; experiment 064 supplied the four-arm labels.

## Resource and safety audit

The authoritative experiment run:

- used one Python process and one BLAS thread;
- peaked at **0.212966919 GiB RSS**, below the 0.75 GiB cap;
- used compact CSV/JSON artifacts and tiny direct `.npy` streams;
- ran no full eigendecomposition, `cell-eval2`, cloud job, download, or
  submission.

Safety incident: before the experiment implementation, two read-only mmap
probes of an Orion matrix retained file-backed pages and reached about
**2.602 GiB peak RSS**, violating the requested cap. They exited immediately
and wrote nothing. The experiment reader was then changed to direct
eight-row/one-column streams; all authoritative and verification runs stayed
under the cap. This is disclosed so the execution record is accurate.

## Decision

**REJECT. Do not ship a per-perturbation gate from this evidence. Keep fixed
experiment-064 arm D.**

Reopen only with a third independent valid somatic holdout. The narrow
hypothesis to test is a predeclared PM-inclusion stump based on source-side PM
geometry; do not reopen broad CD4 breadth/QC routing without a new mechanism.

