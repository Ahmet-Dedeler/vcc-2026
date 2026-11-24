# 071 — exploratory PerturbMatch + CD4 packaging

Status: **exploratory package build active; not authorized for challenge
submission**.

Experiment 065 completed successfully and passed this experiment's
exploratory packaging gate. The exact terminal evidence is bound to run
`20260823T040540Z`.

Active packaging run:

- run ID: `20260823T080948Z`;
- worker: `i-0343d04439564b709`;
- systemd unit: `vcc-exp071-20260823T080948Z.service`;
- output:
  `s3://<your-s3-bucket>/output/exp071-production-pm-cd4/20260823T080948Z/`.

The first authorization attempt failed closed before cloud dispatch because
the validator omitted the literal `exp` in experiment 065's S3 prefix. The
prefix contract and its regression fixture were corrected, after which all 15
local tests passed under a 1 GiB process-tree guard at 0.16 GiB peak RSS.

## Scope and interpretation

This path may construct a `.vcc` candidate for inspection after the exact
experiment-065 run produces a complete, positive, internally consistent
terminal result. It never invokes `vcc submit`.

Experiment 065 is evidence for the PerturbMatch-plus-CD4 direction, but it is
not mathematically identical to the exact v9 production path:

- experiment 065 builds effects at alpha `1.08`, then uses the frozen
  baseline-derived whitening/global scaling from experiment 061;
- exact v9 production uses alpha `1.0796` inside `MultiDonorModel.effect` and
  candidate-specific on-axis preservation;
- experiment 065 uses one RNG stream for control-library picks and count
  emission;
- exact v9 production uses separate deterministic RNG streams for those two
  operations; and
- experiment 064 validates candidate-specific normalization only at effect
  level on the 221-target complete intersection.

Therefore, a positive experiment-065 result authorizes exploratory packaging
only. It does not prove that the exact 300-target A/B/C package improves v9, and
it does not authorize challenge submission. Submission always requires a
separate explicit decision after inspecting the finished package and evidence.

## Fail-closed authorization and provenance

The launcher requires an explicit `EXP065_RUN_ID` and downloads these artifacts
only from:

`output/exp065-pm-cd4-hct-count/${EXP065_RUN_ID}/`

- `exit-code.txt`
- `artifact-sha256.txt`
- `results.json`
- `build-report.json`
- `baseline-agg_results.csv`
- `candidate-run/agg_results.csv`
- `candidate.h5ad.sha256`
- `score-memory-guard.json`
- `runtime.json`

`validate_authorization.py` rejects a wrong run ID or prefix, missing or
duplicate manifest entries, unsafe paths, hash disagreements, non-finite
metrics, incorrect experiment/status/protocol fields, inconsistent score or
effect hashes, and failed exploratory promotion gates. Its output records
`challenge_submission_authorized: false`.

The launcher stages the verified terminal evidence and creates a source
manifest binding the source archive, worker runner, validator, builder, tests,
memory guard, and authorization artifacts by SHA-256. The worker verifies the
downloaded runner before execution, verifies the source archive and manifest
before extraction, checks selected archive members, and reruns the
authorization validator before reading model inputs.

All frozen challenge controls, donor banks, PerturbMatch files, CD4 banks,
production whitener files, and imported helper implementations are
checksum-pinned and verified before use.

## Packaging safety

The worker requires at least 128 GiB nominal RAM. The H5AD build remains
guarded at 55 GiB process-tree RSS; both `vcc prep` passes use a separate
110 GiB guard because the CLI's validation/rewriting phase transiently crossed
55 GiB on a 64 GiB worker.

The worker:

1. holds both `/work/exp061.lock` and `/work/vcc-heavy.lock`;
2. refuses overlap with a scorer, candidate builder, submission builder, or
   `vcc prep`;
3. runs the local synthetic/adversarial tests;
4. builds the H5AD under a 55 GiB process-tree memory guard;
5. runs `vcc prep --dry-run`;
6. writes the real package to a hidden `.partial-*.vcc` path;
7. atomically renames it only after successful prep;
8. creates the artifact manifest and success marker only after checksum
   verification; and
9. uploads a `.vcc` only when the process exits successfully and the success
   marker exists.

Failed runs may upload logs and diagnostic metadata, but exclude every H5AD,
`.vcc`, `.vcc.sha256`, and partial package.

The S3 bucket has a short one-day lifecycle. Any staged evidence or completed
package is temporary, not durable; preserve an explicitly approved artifact
elsewhere before expiry.

No launch has been performed. If experiment 065 later becomes terminal and a
packaging launch is explicitly approved, the invocation requires the exact run
ID:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
EXP065_RUN_ID=20260823T040540Z \
  bash experiments/071-production-pm-cd4/aws/launch.sh
```
