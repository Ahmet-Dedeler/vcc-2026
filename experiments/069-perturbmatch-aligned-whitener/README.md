# 069 — PerturbMatch-aligned whitener

Status: **implemented and locally tested; not launched**.

## Question

Experiment 064's arm D uses a 50:50 incumbent/PerturbMatch K562 processing
ensemble, but its whitening operator was fitted from incumbent K562 alone.
Does fitting the covariance to the actual frozen K562 ensemble improve
effect-level PDS?

## Frozen comparison

Valid held folds are Orion HCT116 and Orion HEK293T. The common full gallery is
all 289 challenge targets measured in both held donors.

Both arms use:

- the frozen 50:50 incumbent/PerturbMatch K562 slot;
- incumbent fallback for PerturbMatch-missing targets and genes;
- the non-held Orion donor;
- the frozen post-whitening CD4 three-condition mean;
- CD4 mode `gene_norm`, fraction `0.05`, and complete-three-condition target
  support only;
- gamma `0.65`;
- raw centered, untransformed held truth.

The only scientific difference is the operator:

- **baseline:** existing experiment-064 arm-D whitener, fitted from incumbent
  K562 plus the non-held Orion donor;
- **candidate:** covariance fitted from the 50:50 K562 processing ensemble plus
  the non-held Orion donor.

The candidate covariance gives the K562 ensemble and Orion donor exactly equal
total source mass. It uses the incumbent K562 non-panel row universe, with
incumbent fallback where PerturbMatch lacks a row or gene. All 300 challenge
perturbation rows and all 300 target-gene coordinates are excluded. The two
K562 views are independently centered and median-norm scaled before blending;
the resulting K562 ensemble is then centered and median-norm scaled as one
covariance source, matching the treatment of the Orion source.

There is one candidate and one gamma. No target, fold, or hyperparameter
selection occurs.

## Stopping gate

The candidate advances to at most one count-level confirmation only if:

1. both held-fold PDS deltas are nonnegative;
2. macro PDS gain is at least `+0.003`; and
3. the paired target-cluster bootstrap 95% lower bound is positive.

Failure of any condition means stop.

## Prerequisite

The AWS launcher refuses to dispatch unless experiment 061 has:

- exit code `0`;
- status `scored`;
- nonnegative count-level six-metric average delta;
- baseline and candidate aggregate files matching the hashes in
  `results.json`;
- a candidate H5AD checksum matching `build-report.json`; and
- an S3 candidate object size matching the build report.

The gate is checked locally before dispatch and again on the worker.

## Local tests

```bash
PYTHONDONTWRITEBYTECODE=1 \
uv run python experiments/069-perturbmatch-aligned-whitener/memory_guard.py \
  --max-rss-gib 1.0 \
  --poll-seconds 0.25 \
  --report /tmp/exp069-test-memory.json \
  -- \
  uv run python experiments/069-perturbmatch-aligned-whitener/test_run.py
```

## Authorized cloud launch

The launcher never creates or starts an instance and performs no submission.
It requires an already-running worker plus explicit approval:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
INSTANCE_ID=<existing-running-worker> \
bash experiments/069-perturbmatch-aligned-whitener/aws/launch.sh
```

The worker runs synthetic tests under a 1 GiB guard, then exactly one
effect-level experiment under a 55 GiB guard. It does not invoke
`cell-eval2`, build a challenge submission, stop an instance, or terminate an
instance.
