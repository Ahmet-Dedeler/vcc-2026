# 086 — exact-margin sparse-responder emission

Status: implementation and bounded synthetic validation complete. The real HCT116
comparison is intentionally **NO-GO** until experiment 077 has completed and its
terminal files have been bound to an AWS-authoritative SSM attestation.

## Question

Experiment 077 emits every low-alpha perturbation with the exact-total quasi
emitter at spread `0.25`. That preserves each gene's aggregate count but makes
most cells sit close to the aggregate mean. Against zero-inflated real controls,
this can make weak or null genes rank-shift in nearly every predicted cell and
become false-positive Wilcoxon calls.

This experiment asks whether the same selected experiment-077 aggregate profile
can instead be represented by:

- a deterministic responder subset carrying nearly all composition change; and
- nonresponders that remain exact control-like template rows whenever the
  negative effects make that mathematically possible.

The frozen production candidates are responder fractions `0.35`, `0.50`, and
`0.65`. At most the two smallest fractions passing prediction-only construction
gates may be scored. There is no submission path.

## Construction

For one perturbation, let:

- `B` be a `400 × G` non-negative integer control-like baseline;
- `r_i = sum_g B[i,g]` be its fixed row/library margins;
- `c_g = sum_i B[i,g]` be its baseline columns;
- `t_g` be experiment 077's exact quasi-emitter columns; and
- `sum_g c_g = sum_g t_g`.

Responder rows are selected by a SHA-256 ordering of `(seed, target, row)`.
For every gene with `t_g < c_g`, counts are removed from responder rows first,
using deterministic integer Hamilton apportionment without creating negative
counts. If responders do not contain enough of that gene, only the unavoidable
remainder is removed from nonresponders and recorded as leakage. The removed
counts create row-specific deficits. For genes with `t_g > c_g`, the exact
positive mass is allocated into those deficits, again preferring responders.

Because total positive and negative deltas are equal, all deficits close exactly.
The result therefore has:

```text
sum_i X[i,g] = t_g       for every gene
sum_g X[i,g] = r_i       for every cell
X[i,g] >= 0 and integer
```

Rows are rejected if their fixed library is outside `[1, 1,000,000]`.

### When exact control-like nonresponders are possible

Keeping every nonresponder entry unchanged requires, for every down-regulated
gene,

```text
sum_{i in responders} B[i,g] >= c_g - t_g.
```

This condition is both necessary and sufficient for this construction. If it
fails, no allocation can reach `t_g` while leaving every nonresponder unchanged:
the missing negative mass must be removed somewhere. The emitter fails closed
only on impossible global margins; otherwise it minimally clips nonresponders
and reports exact all-gene and scored-gene leakage.

`representative_subset_screening_fraction_* = max_g(1 - t_g/c_g)` is only a
screening statistic under a representative-subset assumption. It is not a
feasibility verdict. The authoritative diagnostics are
`selected_subset_exact_control_feasible_*`, computed from the actual selected
rows and emitted matrix: they are true exactly when no nonresponder entry
changed in the relevant gene set.

## Metric semantics

- **PDS is exactly fixed.** `cell-eval2` bulk-log pseudobulk is a function only
  of each perturbation's column totals and grand total. Both are bit-identical
  to experiment 077.
- **The plug-in MSE distance is exactly fixed** for the same reason.
- **Corrected MSE is not automatically invariant.** `cell-eval2` 0.16.0 also
  subtracts a capped prediction-side delete-one-cell jackknife correction,
  bounded in total by the fixed across-perturbation pseudobulk budget.
  Experiment 086 computes the official jackknife quantity for every target and
  permits scoring only if each sparse correction is no smaller than its quasi
  counterpart. With pseudobulks fixed, this makes the aggregate allowed
  correction non-decreasing, so corrected MSE cannot worsen. The scorer result
  is checked again after evaluation.
- **NMAE is not fixed.** It uses arithmetic means of per-cell CPM. Exact raw
  columns do not fix that quantity when row libraries differ.
- **Fidelity, reach, and Jaccard are the intended moving metrics.** They depend
  on per-cell CPM ranks and BH significance.

This is a distributional emission experiment, not a new biological predictor.

## Gates before any real score

`gate.py` rejects experiment 077 unless all of the following are true:

1. `exit-code.txt` exists and is exactly `0`;
2. the byte-identical local snapshot of exp077's launch metadata matches the
   frozen SHA-256 and fields in `exp077-trusted-provenance.json`;
3. AWS SSM's launch-command record proves that the expected instance received
   the exact frozen run id, bucket, input/output prefixes, and runner command;
4. a separate successful terminal SSM command runs the content-addressed
   `aws/attest_exp077.py` and returns canonical JSON binding the exact source
   archive, exp077 runner, grid, summarizer, exit marker, result files, and
   checksum manifest;
5. local copies of the terminal outputs match those SSM-attested hashes;
6. result and summary tags, raw metrics, scaled metrics, deltas, and gate
   booleans agree when independently recomputed;
7. the only candidates are alpha `0.35` and `0.40`;
8. both candidates improve experiment 077's six-metric average;
9. the selected candidate passes experiment 077's individual promotion gate;
10. `submission_authorized` is exactly `false`.

The co-located `artifact-sha256.txt` remains a consistency check only. It is not
treated as an authentication root because the exp077 bucket has neither
versioning nor Object Lock.

### Required exp077 terminal provenance addition

The exp077 run launched at `20260823T060938Z` predates this requirement, so
`exp077-trusted-provenance.json` intentionally has
`terminal_attestation: null`. Experiment 086 must remain fail-closed until all
of these steps happen after the exp077 systemd unit exits successfully:

1. Upload the exact `aws/attest_exp077.py` bytes to the content-addressed S3 key
   recorded in `required_terminal_attestation`, using `If-None-Match: *`.
2. Issue one dedicated `AWS-RunShellScript` command to the original exp077
   instance. Its command list must exactly match the four commands constructed
   by `gate.py`: download the attestor, verify its pinned SHA-256, then run it
   against `/work/exp077-20260823T060938Z`.
3. Confirm AWS reports both the command and invocation as `Success`, with
   response code `0`. The attestor's canonical one-line JSON must be the only
   stdout and must report a successful inactive systemd unit.
4. Record only that terminal SSM command id under `terminal_attestation` in
   `exp077-trusted-provenance.json`. At launch, 086 fetches both SSM records
   directly from AWS and verifies the stdout hashes against the downloaded
   terminal files.

For the pinned run, the SSM command list must be byte-for-byte:

```text
set -Eeuo pipefail
aws s3 cp s3://<your-s3-bucket>/input/exp086-provenance-attestor/c28cc6948d78fc830c22a891415d706f3f6ce58af3071aae9608eb52a3310475/attest_exp077.py /root/exp086-attest-exp077.py --only-show-errors
printf '%s  %s\n' 'c28cc6948d78fc830c22a891415d706f3f6ce58af3071aae9608eb52a3310475' /root/exp086-attest-exp077.py | sha256sum --check - >/dev/null
python3 /root/exp086-attest-exp077.py --run-id 20260823T060938Z --instance-id i-0e4070204e83188fb --region us-east-1 --source-prefix input/exp077-low-alpha-union/20260823T060938Z --output-prefix output/exp077-low-alpha-union/20260823T060938Z --work-root /work/exp077-20260823T060938Z
```

Until this exists, neither a rewritten S3 manifest nor internally consistent
rewritten result files can open the gate.

After that gate, `run.py` verifies the frozen experiment-077 source archive,
source files, HCT116 reference, whitener, challenge panel, experiment-043
template emitter, and every donor file by SHA-256.

For each responder fraction, scoring is allowed only if:

- every target column is exact and all row margins are preserved;
- every row remains within the legal library range;
- scored-gene nonresponder clipping is at most `1%`;
- at least `95%` of absolute composition change lies in responders; and
- every target's jackknife correction is non-decreasing versus exact-total
  quasi.

Only the two smallest passing fractions are scored, sequentially. The runner
contains no VCC preparation or submission command.

## Local validation

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  experiments/086-sparse-responder-emission/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.05 \
  --report experiments/086-sparse-responder-emission/test-memory.json \
  -- uv run python experiments/086-sparse-responder-emission/test_emitter.py

PYTHONDONTWRITEBYTECODE=1 python3 \
  experiments/086-sparse-responder-emission/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.05 \
  --report experiments/086-sparse-responder-emission/synthetic-memory.json \
  -- uv run python experiments/086-sparse-responder-emission/synthetic_study.py
```

See `RESULTS.md` and `synthetic-results.json` for measured evidence.

## Guarded AWS runner

`aws/launch.sh` can target only an already-running worker with at least
120,000 MiB RAM. It refuses to run unless
`I_HAVE_EXPLICIT_CLEARANCE=YES`, performs the experiment-077 terminal gate
before packaging or dispatch, and never creates or starts an instance.

The remote `aws/run.sh` takes a non-blocking global heavy-process lock, rejects
pre-existing scorer/training jobs, requires at least 100 GiB currently
available RAM and 45 GiB disk, re-runs the terminal gate before downloading
large inputs, verifies all frozen dependencies, runs unit tests under the local
1 GiB guard, then runs experiment 086 under the repository's 110 GiB process
tree guard. Predictions (`*.h5ad`) are not uploaded.

Both launch and worker scripts validate every externally supplied identifier
before interpolation. The launcher rejects existing local metadata and
non-empty S3 input/output prefixes, uses a unique temporary source archive, and
claims both prefixes with conditional `If-None-Match: *` reservation objects.
The worker accepts only canonical paths, rejects symlinked or pre-existing work
directories, and never recursively deletes a work directory.

No AWS workload has been launched.
