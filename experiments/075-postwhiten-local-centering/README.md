# 075 — Post-whitening reciprocal local query centering

## Question

Can a query-only local hubness correction improve cosine retrieval after the
frozen union whitener?

For each post-whitened prediction `q_i`, this experiment finds reciprocal
cosine-nearest neighbors among the other predicted perturbations, subtracts a
fraction of their unit-vector centroid, renormalizes, and restores the exact
input row norm:

```text
u_i  = q_i / ||q_i||
R_i  = reciprocal-kNN(i) among the submitted queries
m_i  = mean(u_j for j in R_i)
q'_i = ||q_i|| * normalize(u_i - beta * m_i)
```

The transform uses no hidden reference gallery, no perturbation label, no
challenge context control, and no truth vector. It is directly deployable on a
complete post-whitened panel before count generation. Exact row-norm
preservation isolates direction from the already calibrated effect magnitude.

## Why this is not a repeated experiment

The audit covered all retrieval/geometry results from 017 through 073.

- 017 tested pooled covariance whitening and a single global top-component
  removal, not local query centering.
- 030 tested entrywise row/column power sharpening, not neighborhood
  residualization.
- 044 learned a pre-whitening surrogate for the hidden gallery-side CSLS
  density penalty.
- 048/052 tested source-covariance subspaces before or around whitening.
- 067 tested control-variance diagonal transport.
- 069 changed the covariance fit to match PerturbMatch.
- 073 stages post-whitening CSLS-density distillation, which still learns a
  hidden-gallery penalty from auxiliary truths.

Experiment 075 instead changes each query vector by a locally estimated common
direction obtained only from the query set. Reciprocal neighborhoods are used
to reduce the influence of query hubs. Rows with fewer than two reciprocal
neighbors fall back deterministically to their ordinary top-`k` neighbors.

## Frozen protocol

This experiment consumes experiment 073's authenticated post-whitening query
cache. It does not fit or modify an eigensystem.

- Outer roles: K562-GWPS, Orion HCT116, Orion HEK293T.
- Query for each outer role: equal consensus of the other two main donors,
  transformed by that role's frozen experiment-062 union whitener at
  `gamma=0.65`.
- Selection panel: 600 public perturbations with zero overlap with the 2026
  VCC panel.
- Selection truths: RPE1 and Nadig-Jurkat.
- Outer panel: 259 shared 2026 VCC targets, opened only after every selector is
  frozen.
- Inner axis: 5,582 genes covered by both auxiliary truth donors.
- Outer axis: the full leakage-safe experiment-073 score axis.
- All 300 VCC target coordinates and all 600 selection-target coordinates are
  absent from the effective query cache.
- Held outer truth remains raw-centered and untransformed.

The fixed candidate grid is:

```text
k:    8, 16, 32, 64
beta: 0.10, 0.25, 0.50, 0.75
```

For each outer source role, a candidate is eligible only if its PDS delta is
nonnegative against both auxiliary truths and its mean gain is at least
`+0.001`. Selection maximizes the worst teacher delta, then mean delta, with a
`5e-4` tie window favoring smaller `beta` and broader `k`.

## Predeclared nested promotion gate

Promotion requires every item below:

1. The all-inner-data deployment selector chooses a non-identity candidate.
2. Leave-one-source-role-out inner validation is nonnegative in all six held
   role × teacher cells.
3. Its pooled held-role inner gain is at least `+0.001` PDS.
4. The frozen deployment candidate is nonnegative on all three outer donor
   folds.
5. Its outer macro PDS gain is at least `+0.003`.
6. The 95% target-cluster bootstrap lower bound is above zero.

The same deployment candidate is used for all three outer folds. Outer
candidate grids are not evaluated, so the VCC-panel outcomes cannot select a
configuration.

## Memory and execution

Synthetic tests stay below a 1 GiB process-tree ceiling and open no real donor
matrix or dense whitener:

```bash
PYTHONDONTWRITEBYTECODE=1 \
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 1 --min-free-percent 30 --poll-seconds 1 -- \
  uv run --no-project --python 3.12 --with numpy python \
  experiments/075-postwhiten-local-centering/test_run.py
```

The real evaluator itself is bounded at 2 GiB. The only heavy operation is
upstream creation of experiment 073's post-whitening query cache, which belongs
on an isolated cloud worker. No `cell-eval2`, counts, submission artifact, or
VCC API call is involved.

The guarded launcher reuses experiment 073's audited staging code and frozen
artifact manifest. Their current hashes are pinned in
`dependency-sha256.txt`; launch fails closed if either dependency changes. It
requires an already-running, idle worker with at least 64 GiB RAM and cannot
create or start an instance:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES INSTANCE_ID=i-... \
bash experiments/075-postwhiten-local-centering/aws/launch.sh
```

When an authenticated query cache exists:

```bash
PYTHONDONTWRITEBYTECODE=1 \
uv run --no-project --python 3.12 --with numpy python \
  experiments/075-postwhiten-local-centering/run.py \
  --truth-cache /path/to/exp052-cache \
  --query-cache /path/to/exp073-query-cache \
  --upstream-manifest \
    experiments/073-postwhiten-csls-distill/artifact_manifest.json \
  --genes-csv data/challenge/gene_names.csv \
  --panel-csv data/challenge/pert_counts.csv \
  --guarded
```

That command evaluates only; it does not build or submit a challenge file.
