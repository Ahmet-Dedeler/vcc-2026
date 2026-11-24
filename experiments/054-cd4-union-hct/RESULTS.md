# 054 — CD4 under the clean HCT union whitener

Date: 2026-08-22

## Verdict

Do not add the primary-CD4 donor to v9.

The frozen `0.50` three-condition CD4 blend that won before whitening loses
badly after the exact clean-HCT union operator:

| candidate | all 132 HCT targets | delta vs v9-shaped baseline |
|---|---:|---:|
| **K562 + HEK293T union-whitened baseline** | **0.80095** | — |
| CD4 raw, preserve missing genes | 0.73791 | **-0.06303** |
| CD4 raw, zero-vote missing genes | 0.73710 | -0.06384 |
| CD4 gene-normalized, preserve missing genes | 0.76544 | **-0.03551** |
| CD4 gene-normalized, zero-vote missing genes | 0.76498 | -0.03597 |

The loss remains on the 130 CD4-covered targets and the 125 targets present
in all three CD4 conditions:

| primary candidate subset | baseline PDS | candidate PDS | delta |
|---|---:|---:|---:|
| all HCT targets | 0.80095 | 0.76544 | -0.03551 |
| CD4-covered | 0.80493 | 0.76888 | -0.03605 |
| all three CD4 conditions | 0.80763 | 0.77814 | -0.02950 |

Production-style self-renormalization and a fixed baseline renormalization
gave the same PDS values. Missing-gene handling is therefore not the cause.

## Interpretation

Experiments 047 and 051 showed a repeatable positive CD4 contribution in the
unwhitened space. Experiment 054 shows that this benefit does not compose with
the model's dominant production lever. The union whitener was fitted on
K562/HEK source covariance; adding a lineage-distant CD4 direction before that
operator moves predictions into low-variance source directions that are
strongly amplified but do not retrieve the HCT truth.

This is an interaction failure, not evidence that the CD4 data are poor. It
also reinforces the project rule that every donor or normalization change
must be tested after the exact production geometry.

## Protocol and safety

- Held target context: Orion HCT116, 132 VCC-panel perturbations.
- Source donors: K562-GWPS and Orion-HEK293T.
- Operator: experiment 040's clean held-HCT union whitener, gamma `0.65`,
  17,942 columns.
- CD4 coverage: 130/132 HCT targets in at least one condition, 125/132 in all
  three.
- All 300 challenge target coordinates excluded from scoring and whitening.
- No `cell-eval2`, count emission, or submission.
- Synthetic tests: 6 passed.
- Cloud peak RSS: 2.86 GiB.

## Decision

Close fixed CD4 blending under the current union-whitening pipeline. Keep the
compact CD4 donor for future operator retraining or a separately gated
ensemble, but do not modify v9 from this experiment.
