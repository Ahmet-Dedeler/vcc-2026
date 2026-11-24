# Experiment 144 — cross-donor consensus covariance

Status: **corrected source-only v2 scaffold; unauthorized and unevaluated**

Date: **2026-08-31**

## Question

The incumbent whitener estimates covariance from stacked individual donor
effects, while the submitted predictor averages donors. Donor-private variance
therefore enters the operator even though averaging suppresses it. Experiment
144 asks whether the operator should instead be fitted to matched, same-target
source consensus rows:

```text
m[t] = mean_d normalized(centered(effect[d,t]))
C_consensus = mean_t m[t] m[t]^T
```

This is a source-only geometry experiment. It does not yet estimate a score or
authorize access to any development or held outcome.

## Why v1 is quarantined

`work/preoutcome-v1/` and
`SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json` are immutable historical artifacts.
Their SHA-256 values are frozen in `preregistered_protocol.json`.

The v1 implementation is not reusable because:

- each donor was centered on its complete perturbation table before forbidden
  rows were removed, so excluded rows affected every retained row indirectly;
- donor renormalization came from upstream full-table statistics;
- “full rank” was inferred from rows being greater than columns;
- it materialized a ten-candidate grid and 613-target prediction arrays;
- its later loader bypassed RSS enforcement with a zero limit.

The corrected path writes only to `work/source-only-v2/`. Existing v1 files
are verified but never opened as NumPy arrays, copied, renamed, or overwritten.

## Frozen v2 candidate

Exactly two choices exist:

| name | gamma | eta |
| --- | ---: | ---: |
| `identity` | 0.00 | 0.00 |
| `consensus:g=0.65:eta=0.25` | 0.65 | 0.25 |

There is no grid and no tuning. The future candidate application is:

```text
q_consensus = W_consensus(gamma=0.65) q_pre
u = normalize(0.75 normalize(q_identity) + 0.25 normalize(q_consensus))
candidate = ||q_identity|| u
```

The v2 prepare stage fits and records the source operators and this candidate
specification only. It does not load or materialize development predictions.

## Leakage-safe source fit

For each role:

1. Intersect perturbation names across every source donor in that role.
2. Remove the complete union of:
   - all 1,023 frozen Experiment-132 613/410 target names;
   - all 300 VCC panel target names;
   - all targets in `data/donors/k562_essential/perts.txt`.
3. Only then map and load the retained rows from each donor `delta.npy`.
4. Compute each donor center from those retained rows only.
5. Compute each donor median centered-row norm from those retained rows only.
6. Renormalize donors to the median of the role-local donor medians.
7. Average matched same-target rows equally across donors.
8. Center the consensus basis using its retained rows only.
9. Fit covariance and its eigensystem from that basis only.

Excluded rows cannot affect centering, consensus means, normalization,
covariance, rank checks, or another fitted statistic.

## Numerical-rank gate

Rows greater than columns is not a rank test. V2 diagonalizes the float64
covariance and sets:

```text
tol = lambda_max * max(n_rows, n_columns) * eps(float64) * 10
```

Every one of the 4,992 raw eigenvalues must be finite and strictly greater than
`tol`. A materially negative eigenvalue or numerical rank below 4,992 fails the
run before the final directory is written.

The operator then uses:

```text
floor = max(lambda_max * 1e-12, tol)
```

and verifies that the regularized spectrum is finite and strictly positive.
Both the raw rank decision and regularized-spectrum diagnostics are recorded.

## Memory boundary

The frozen ceiling is **18 GiB process-tree RSS** and the free-memory floor is
**25%**.

- `memory_guard.py` launches the command in a new process group and samples
  the complete group every 0.25 seconds across covariance construction, the
  LAPACK eigensolver, compression, hashing, and finalization.
- A process-table error, unavailable free-memory measurement, orphaned
  descendant, interrupt, or limit violation terminates the complete guarded
  process group and fails closed.
- `run.py` refuses preparation unless the exact experiment-local wrapper
  completes a live two-way pipe handshake and enforces limits at least as
  strict as the authorization.
- Internal checkpoints use the same `ps` and `memory_pressure` measurements
  around every source-array loader, fitted-statistic pass, heavy numerical
  boundary, hash pass, and serializer.
- A zero/disabled RSS limit or a lower free-memory floor is rejected.
- The historical shared wrapper remains unchanged because completed
  experiments bind its previous hash.

## Authorization boundary

No `SOURCE_ONLY_PREPARE_AUTHORIZATION.v2.json` exists in this scaffold.
Source preparation must not run until an independent no-truth audit passes and
a hash-bound source-only v2 authorization is created.

There is no development/evaluation command. A later evaluation path must be a
separately reviewed artifact and independently authorized. This scaffold does
not access development or held outcomes, 613/410 values, KOLF values, reserved
expression, challenge controls, `cell-eval2`, counts, builds, or submissions.

## Commands allowed now

Synthetic/static tests:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/144-consensus-covariance/test_run.py
```

Read-only scaffold audit:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/144-consensus-covariance/run.py audit
```

Future source-only preparation, only after the independent audit and v2
source-only authorization:

```bash
uv run python experiments/144-consensus-covariance/memory_guard.py \
  --max-rss-gib 18 --min-free-percent 25 --poll-seconds 0.25 \
  --report-file experiments/144-consensus-covariance/work/source-only-v2-memory.json \
  -- \
  env OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 \
  PYTHONDONTWRITEBYTECODE=1 python \
  experiments/144-consensus-covariance/run.py prepare \
  --max-rss-gib 18 --min-free-percent 25
```

Expected v2 output:

```text
work/source-only-v2/
  candidate-spec.json
  consensus-operator-k562.npz
  consensus-operator-k562.basis.json
  consensus-operator-non_k562.npz
  consensus-operator-non_k562.basis.json
  result.json
  artifact-manifest.json
```

All writes are staged atomically and fail if the final v2 path already exists.
