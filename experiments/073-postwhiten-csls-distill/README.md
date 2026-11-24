# 073 — Post-whitening CSLS-density distillation

## Question

Experiment 044 rejected donor-trained CSLS-density distillation **before**
spectral whitening. It explicitly left one regime untested: learn and evaluate
the query-only correction after the source prediction has passed through the
same source-only union whitener used by the incumbent.

This experiment tests that regime. It reuses the three immutable
experiment-062 gamma-`0.65` eigensystems and does not refit whitening or a full
gene-covariance eigensystem. The low-rank distillation model still uses
experiment 044's bounded randomized SVD.

## Strict nested protocol

The outer folds are K562-GWPS, Orion HCT116, and Orion HEK293T.

For each outer donor:

1. The query is the equal consensus of the other two main donors, using
   experiment 062's full-catalogue centering, source norm calibration, and
   fold-specific frozen union operator.
2. The 600 experiment-044 non-VCC perturbations form the distillation panel.
   The 259 shared 2026 VCC targets form the outer panel. Their overlap is zero.
3. RPE1 and Nadig-Jurkat are the only density teachers. In one inner direction,
   fit on RPE1 and validate on Jurkat; in the other, fit on Jurkat and validate
   on RPE1.
4. Select one experiment-044 candidate from those two held-teacher scores.
   Then refit that model family on both teachers and evaluate it once against
   the held outer donor's VCC-panel truth.
5. The outer donor is absent from the source query, the frozen operator, both
   density fits, and hyperparameter selection. Outer scores are diagnostic
   only and never choose a candidate.

This is stricter than training the transform on the other two experiment-062
outer folds: those folds' operators contain the current outer donor as a
source, which would violate outer LODO.

## Coordinate exclusion

The axes follow experiment 044's split geometry:

- **Outer retrieval score axis: 6,821 genes.** Start from its common
  7,312-gene outer axis and remove the 491 training-target symbols represented
  on that axis.
- **Learned correction axis: 5,582 genes.** Apply 044's stricter saved
  `train_mask`, which additionally removes 1,239 genes without complete
  five-donor teacher coverage.

The 1,239 score-only genes remain in baseline cosine normalization and
retrieval. The learned affine/spectral correction is embedded only on the
5,582-gene teacher-common subset.

- All 300 VCC target coordinates were already absent from each experiment-062
  operator axis.
- All 600 distillation-target coordinates are zeroed before operator
  application and discarded afterward.
- Therefore the effective action is the active principal map `W[A,A]`; neither
  panel's own-target coordinate can contribute to a query, correction basis,
  or retrieval score.
- Held truth remains raw-centered and untransformed.

The experiment-062 eigensystem was fitted before the 600-target panel was
defined, but no excluded coordinate is supplied to or retained from its frozen
linear action. Refitting after deleting those coordinates would test a
different upstream model and is intentionally forbidden here.

## Model family and selection

The candidate family and default grid are intentionally inherited from 044:

```text
affine:   q' = q + eta * b
spectral: q' = q + eta * V diag(a) V^T q
combined: both terms
```

- CSLS neighborhood `k`: `5, 10, 20`
- rank: `8, 16, 32`
- ridge: `0.01, 0.1, 1.0`
- eta: `0.5, 1.0`
- methods: affine/bias, spectral/bands, combined

Selection uses mean inner-donor PDS. Candidates within `5e-4` of the best use
044's deterministic lower-complexity tie-break. A candidate must beat the
inner baseline by at least `0.001`; otherwise the fold selects identity.

Promotion requires all three outer deltas nonnegative, macro PDS gain at least
`+0.003`, and a target-cluster bootstrap lower bound above zero.

## Resource design

The heavy step is three dense applications of saved 17.8k-gene operators. It
is cloud-staged because each eigenvector array is about 1.2 GiB and the matrix
multiplication is compute-heavy. The implementation:

- memory-maps every donor `delta.npy`;
- prepares one donor at a time in bounded row chunks;
- loads only one frozen operator at a time;
- applies it in 16-row chunks;
- retains all 6,821 leakage-safe outer score coordinates;
- uses a 48 GiB process-tree kill guard;
- never invokes `cell-eval2`, emits counts, builds a submission, or calls VCC.
- does not refit whitening/full covariance eigensystems; the small randomized
  SVD used to fit the rank-8/16/32 distillation basis still runs.

The S3 object keys and SHA-256 values for the reused donor banks, experiment-044
truth cache, challenge CSVs, and experiment-062 operators are frozen in
`artifact_manifest.json`.

## Local verification

Synthetic tests do not open real donor data or the 1.2 GiB operators:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/073-postwhiten-csls-distill/test_run.py

ruff check \
  experiments/073-postwhiten-csls-distill/run.py \
  experiments/073-postwhiten-csls-distill/test_run.py \
  experiments/073-postwhiten-csls-distill/memory_guard.py

bash -n experiments/073-postwhiten-csls-distill/aws/run.sh
bash -n experiments/073-postwhiten-csls-distill/aws/launch.sh
```

## Cloud run

The launcher requires an already-running isolated worker with at least 64 GiB
RAM and checks the EC2 type before dispatch. Each run gets a clean
`/work/exp073/<RUN_ID>` directory; the runner independently checks visible
memory and required tools. Packaging sets `COPYFILE_DISABLE=1`, rejects
AppleDouble metadata, and uses `uv --no-project`, so the root README is included
for provenance but is not needed for dependency resolution. The launcher
deliberately cannot create or start an instance.

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
bash experiments/073-postwhiten-csls-distill/aws/launch.sh
```

The launcher currently defaults to the verified worker
`i-0e05556f09944111f`, but `INSTANCE_ID` may override it. Do not run that
command without explicit approval. The current task only implements and
validates the experiment; it does not launch cloud work.
