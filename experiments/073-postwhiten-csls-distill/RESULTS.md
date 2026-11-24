# 073 — Post-whitening CSLS-density distillation

**Date:** 2026-08-23

## Status

**Completed successfully and rejected.**

The guarded biological run completed on AWS with exit code `0`. The nested
selection gate failed decisively, so this correction must not be added to a
submission.

## What is ready

- Strict three-fold outer donor LODO over K562-GWPS, Orion HCT116, and Orion
  HEK293T.
- Strict two-direction inner teacher LODO over RPE1 and Nadig-Jurkat.
- Disjoint 600-perturbation non-VCC training panel and 259-target VCC outer
  panel.
- Both panels' target-gene coordinates excluded from the effective operator,
  learned correction, and retrieval axis.
- Dual-axis geometry matching experiment 044: 6,821 outer retrieval genes and
  a 5,582-gene teacher-common correction basis. The 1,239 valid score-only
  genes are retained in cosine normalization/ranking.
- Exact reuse and hash validation of experiment-062's three gamma-`0.65`
  operators.
- Low-memory, chunked query-cache staging with no whitening/full covariance
  refit. The downstream rank-8/16/32 distiller still uses randomized SVD.
- Experiment-044 affine, spectral-band, and combined distillation families.
- Exact-CSLS oracle diagnostics, per-target PDS output, nested macro result,
  clustered target bootstrap, and a predeclared promotion gate.
- Guarded AWS runner and explicit-approval launcher.

## Local validation

- Synthetic protocol tests: **4 passed**
- Peak guarded process-tree RSS: **0.0818 GiB** under a 1 GiB ceiling
- Ruff lint: **passed**
- Ruff format check: **passed**
- Python AST compilation: **passed**
- `bash -n` for both AWS scripts: **passed**
- Real artifact SHA-256 validation: **passed**
- Real metadata/axis audit: **passed**
  - held K562 operator: 17,806 coordinates
  - held HCT116/HEK293T operators: 17,942 coordinates each
  - outer score axis: 6,821 genes
  - correction axis: 5,582 genes
  - score-only coordinates retained: 1,239 genes
  - training rows: 600
  - outer VCC rows: 259
- Every referenced S3 donor, operator, and compact-cache object: **content
  SHA-256 verified against the frozen manifest**
- Intended worker `i-0e05556f09944111f`: **running, SSM Online,
  `r7i.2xlarge`, 65,536 MiB**

## Promotion gate

Retain the method only if:

1. every outer donor delta is nonnegative;
2. nested macro PDS gain is at least `+0.003`;
3. the 95% target-cluster bootstrap lower bound is above zero.

The fixed-grid outer table is emitted as an explicitly optimistic diagnostic
and cannot override the nested verdict.

## Biological result

| held donor | baseline PDS | selected PDS | delta |
|---|---:|---:|---:|
| K562-GWPS | 0.795412 | 0.795696 | +0.000284 |
| Orion HCT116 | 0.792493 | 0.776286 | -0.016207 |
| Orion HEK293T | 0.783155 | 0.766140 | -0.017015 |
| **macro** | **0.790353** | **0.779374** | **-0.010979** |

The target-cluster bootstrap interval was
`[-0.019355, -0.002459]`, with probability of a positive delta `0.00595`.
All three predeclared promotion conditions failed.

The exact-CSLS oracle retained `+0.008216` macro headroom, but the
leakage-safe teacher distillation could not learn a transferable correction.
Post-whitening density information is therefore not a usable production
signal under this protocol.

Artifacts:

- run id: `20260823T060938Z`
- local mirror: `cloud-output/20260823T060938Z/`
- peak process-tree RSS: `2.708 GiB`
- exit code: `0`

## Safety state

- Cloud launched: **yes; completed**
- Whitening/full covariance eigensystems refitted: **0**
- Randomized low-rank distillation SVD: **run on real data**
- `cell-eval2` invoked: **no**
- Counts generated: **no**
- Submission built or sent: **no**
