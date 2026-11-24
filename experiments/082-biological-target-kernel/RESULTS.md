# 082 — Leakage-safe biological target-kernel residual

Status: **NO-GO at the non-panel pre-gate**.

## Valid run — 2026-08-30

Run ID: `20260830T042050Z`

The guarded AWS run completed with exit code 0. All 22 synthetic/contract
tests passed, every leakage audit flag passed, and all 27 downloaded output
hashes verify. Peak process-tree RSS was 4.44 GiB for staging and 0.52 GiB
for evaluation. The `r7i.8xlarge` worker was stopped after artifact recovery.

The frozen selector chose **identity**. No nonidentity candidate passed the
requirement of nonnegative PDS delta in every outer-donor × inner-teacher
cell:

| Candidate | Macro raw-PDS delta | Worst cell |
| --- | ---: | ---: |
| rank 32, ridge 1, beta 0.02 | **-0.005120** | -0.005718 |
| rank 16, ridge 1, beta 0.02 | -0.005158 | -0.005746 |
| rank 64, ridge 1, beta 0.02 | -0.005214 | -0.005913 |

Larger beta values degraded monotonically. There were zero eligible
nonidentity configurations among the 60 tested.

Because the source-only pre-gate failed:

- HCT116, HEK293T, and KOLF panel truth was **not opened**;
- no held-panel result was produced;
- no count generation or submission was allowed.

Conclusion: a linear ESM2 target-kernel residual is not a transferable
correction to arm D under this protocol. Keep the family closed unless a
materially different biological representation or mechanism is introduced.

Artifacts are under `cloud-output/20260830T042050Z/`.

## Earlier invalid attempts

The implementation follows hypothesis 2 and the common leakage contract in
experiment 078. It includes:

- direct reconstruction of the frozen experiment-064 arm-D formula with
  strict HCT116/HEK293T outer holdout;
- exact `1e-6` arm-D reproduction on experiment 064's pinned 221-target
  intersection, with fallback extensions audited separately;
- separate development and final score axes: all available historical target
  coordinates are excluded from source normalization and a principal
  restriction of the development whitener, while all 109 targets absent from
  the compact master axis remain valid perturbation rows (49 also lack an
  official coordinate) and the final panel uses the pinned 17,806-gene axis
  excluding exactly the 300 VCC target coordinates;
- a strict 5,582-gene teacher-common development axis, so genes absent from
  either RPE1 or Nadig-Jurkat are excluded rather than zero-filled as truth;
- deterministic non-panel fit/tune/test partitions;
- fit-only ESM2 standardization;
- residual PCA ranks 16/32/64;
- ridge strengths 1e-3 through 10;
- beta 0 through 0.20, with exact identity fallback;
- complete 600-target development galleries and complete 290-target HCT116,
  294-target HEK293T, and 273-target KOLF final galleries;
- embedding-permutation and equal-dimensional Gaussian nulls that each repeat
  the full selection procedure;
- independently source-derived clusters for each HCT116/HEK293T outer-fold
  bootstrap, with KOLF excluded from the primary macro/bootstrap;
- fold-specific source-only final null strata for HCT116, HEK293T, and KOLF,
  with KOLF null results reported separately;
- duplicate-preserving Orion axis handling: all 38,606 matrix columns remain
  aligned while unique projections use the deterministic last occurrence of
  each symbol;
- complete-source, partial-source, direct-row, and missing-source reporting;
- exact final model weights frozen and hashed before role-specific held-panel
  rows are selected and loaded for evaluation;
- independent KOLF confirmation using a three-donor source-only arm-D query,
  KOLF-specific weights fitted only on the 600 non-panel rows, and the complete
  273-target KOLF gallery;
- all real and null final models serialized before any role-specific held
  panel is numerically evaluated, with no post-panel refit;
- hashes, provenance, memory guards, fail-closed held-truth access,
  source/output S3 prefix-reuse rejection, strict timestamp run IDs, canonical
  path containment, and refusal to delete or reuse an existing local run path.

Local validation: **22 tests passed** under the 1 GiB process-tree guard;
recorded peak process-tree RSS was **0.161 GiB** and the guard did not kill the process.
Ruff format/check, Python compilation, shell syntax, JSON parsing, and pinned
dependency hashes also pass.

Two cloud attempts on 2026-08-29 are invalid and carry no scientific result:
`20260829T215910Z` failed while staging a missing cache file, and
`20260829T221225Z` was stopped before held-truth evaluation after the
across-perturbation truth-centering defect was found. Their S3 logs and the
archived local launch metadata are retained. They do not alter the valid
2026-08-30 NO-GO result above.
