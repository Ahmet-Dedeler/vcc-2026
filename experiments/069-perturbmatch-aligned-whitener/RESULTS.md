# 069 — PerturbMatch-aligned whitener

Date: 2026-08-23

## Decision

**Reject. Keep the incumbent source-only whitener used by experiment 064 arm
D. Do not run a count-level confirmation.**

Fitting covariance to the same frozen 50:50 incumbent/PerturbMatch K562
ensemble used by the prediction made retrieval worse on both valid held
somatic lines:

| held line | incumbent arm-D PDS | aligned-whitener PDS | delta |
|---|---:|---:|---:|
| Orion HCT116 | 0.8023716840 | 0.7943939831 | -0.0079777009 |
| Orion HEK293T | 0.7802287582 | 0.7708093041 | -0.0094194541 |
| **macro** | **0.7913002211** | **0.7826016436** | **-0.0086985775** |

The paired target-cluster bootstrap interval was
`[-0.0133902826, -0.0041270185]`, with positive probability `0.0001`.
All three predeclared gates failed:

- both fold deltas nonnegative: false;
- macro gain at least `+0.003`: false;
- bootstrap lower bound positive: false.

The incumbent K562 covariance geometry transfers better than the apparently
more internally consistent PerturbMatch-aligned geometry. PerturbMatch remains,
at most, a small prediction-side ensemble component; it should not enter
covariance fitting.

## Protocol and leakage audit

- Valid held folds: Orion HCT116 and Orion HEK293T.
- Full common gallery: 289 targets.
- Candidate count: one.
- Gamma: fixed `0.65`.
- PerturbMatch fraction: fixed `0.50`.
- CD4 residual: complete mean3, `gene_norm`, fixed fraction `0.05`.
- Candidate covariance sources: one equal-mass K562 processing ensemble plus
  the non-held Orion donor.
- All 300 challenge perturbation rows were excluded from covariance fitting.
- All 300 panel coordinates were excluded from both operators and retrieval.
- Held truth was never used in candidate covariance construction.
- No count emission, `cell-eval2`, challenge build, or submission.

## Cloud execution

- Successful run id: `20260823T041514Z`.
- Worker: `i-08d90ae88c7b210ed`, `r7i.2xlarge`.
- Peak guarded process-tree RSS: `8.862 GiB`.
- Exit code: `0`.
- Result SHA-256:
  `afce25dee17713aa686e4f9f3f903c8f537d0a4adafcb57b4a459303164c0bf8`.
- Durable S3 prefix:
  `s3://<your-s3-bucket>/output/exp069-perturbmatch-aligned-whitener/20260823T041514Z/`.
- Compact artifacts and whitener metadata are preserved under
  `cloud-output/20260823T041514Z/`; the two 1.2 GiB operators remain in S3.

The first dispatch failed before tests or scientific computation because the
source archive omitted the repository-level `README.md` required by
`pyproject.toml`. The launcher was corrected; only the successful run above is
scientific evidence.

## Local verification

- Synthetic harness: `11 synthetic tests passed`.
- Harness self-reported Python peak RSS: `0.051 GiB`.
- External process-tree peak RSS: `0.078537 GiB`.
- External process-tree guard: `1.0 GiB`; return code `0`, no kill.
- AWS shell syntax: `bash -n` passed for `aws/run.sh` and `aws/launch.sh`.
- Ruff: all checks passed under the `1.0 GiB` guard; process-tree peak
  `0.006058 GiB`.
- Python compile check: all four Python files compiled under the `1.0 GiB`
  guard; process-tree peak `0.006012 GiB`.

The synthetic suite covers the frozen constants, exact 300-row/coordinate
exclusion in a tiny fitted whitener, equal source mass, K562 target/gene
fallback, complete-three-condition CD4 support, paired target bootstrap,
conjunctive stopping gate, experiment-061 artifact/hash/exit-code gate, and
AWS launch/run refusal invariants.
