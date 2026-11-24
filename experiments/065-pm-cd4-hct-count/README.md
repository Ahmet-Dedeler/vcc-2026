# 065 — frozen PerturbMatch + CD4 HCT116 count-level test

Status: **completed successfully; positive HCT116 count-level result**.

Completed run:

- run id: `20260823T040540Z`
- worker: `i-0343d04439564b709`
- systemd unit: `vcc-exp065-20260823T040540Z.service`
- S3 output:
  `s3://<your-s3-bucket>/output/exp065-pm-cd4-hct-count/20260823T040540Z/`

The candidate build and single `cell-eval2` scorer completed with exit code
`0`. The candidate improved raw PDS by `+0.0046264168` and the reference-scaled
six-metric average by `+0.0020082200` over the exact clean gamma-0.65 HCT116
baseline. See `RESULTS.md` and
`cloud-output/20260823T040540Z/results.json`.

This is the production-shaped HCT116 count-level confirmation of experiment
064's frozen joint candidate:

- reconstruct experiment 061's exact 50:50 incumbent/PerturbMatch K562-slot
  candidate by calling its implementation directly;
- fail closed unless all four reconstructed pre/post effect hashes match
  experiment 061's successful build report;
- retain experiment 061's exact behavior for PerturbMatch-missing targets and
  genes;
- use experiment 040's existing candidate-independent union whitener at gamma
  `0.65`;
- add the experiment-062 CD4 Rest/Stim8hr/Stim48hr mean residual after
  whitening, using `gene_norm`, fraction `0.05`, and the same global
  norm-preserving scalar;
- apply the CD4 residual only where a target is present in all three CD4
  conditions, matching experiment 064's frozen mean3 gallery; targets with
  incomplete CD4 support retain the PerturbMatch candidate before the one
  global preserve scalar;
- keep alpha `1.08`, CPM gate `5`, quasi exact-total spread `0.25`, and seed
  `0`;
- compare against the exact clean gamma-0.65 HCT116 baseline and report all six
  `cell-eval2==0.16.0` metrics.

All 300 VCC target coordinates are excluded from the whitener and from CD4
residual support. Existing values on those coordinates may only receive the
single global preserve scalar, matching experiment 062.

## Safety and serialization

The local launch wrapper and cloud runner together:

- require experiment 061's successful `results.json`, candidate aggregate,
  and zero exit code in S3 before dispatch; the worker then rechecks the
  already-local exit code and artifacts because its IAM role cannot read
  `output/*`;
- reuses only experiment 061's existing worker and refuses any other EC2 state
  instead of creating a new instance;
- acquires both `/work/exp061.lock` and `/work/vcc-heavy.lock`;
- refuses to start if a scorer or experiment-061/065 builder process exists;
- completes the candidate build before starting exactly one scorer;
- guards the build at 32 GiB process-tree RSS and scoring at 55 GiB;
- uploads the candidate before scoring and uploads all outputs on exit;
- does not submit, terminate, stop, or shut down the worker.

The current frozen exp061 worker default is `i-0343d04439564b709` in
`us-east-1`. If that instance is stopped after experiment 061, the launch
wrapper starts the same instance and waits for it; if it is terminated or in
any unexpected state, the wrapper fails.

## Launch after explicit clearance

Do not run this until Ahmet explicitly clears the heavy execution:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
  bash experiments/065-pm-cd4-hct-count/aws/launch_reuse_exp061.sh
```

No challenge submission is created or sent.
