# 052 — Rank-64 source subspace versus v9 union whitening

This is the decisive follow-up to experiment 048. It compares the deployable
rank-64, gamma-2 source-only subspace transform against a leak-free refit of the
full-rank union whitener used by v9.

Five predeclared pipelines are compared:

1. identity
2. union only
3. subspace only
4. union then subspace
5. subspace then union

Matched random-subspace variants are controls and are never eligible for model
selection.

## Leakage boundary

For each outer donor among K562-gwps, Orion-HCT116, and Orion-HEK293T:

- the outer donor is absent from every fitted operator and every inner score;
- the operator-fit background excludes all 300 VCC targets and all 600
  non-VCC inner-panel perturbations;
- all 900 target-gene coordinates are removed before either dense or low-rank
  transform;
- composition selection uses only five-fold PDS on the 600 non-VCC panel;
- each fold aggregates held query rows, but every query is ranked against the
  full 600-candidate catalog;
- all inner selections are frozen before any 259-target outer panel is loaded;
- held-donor truth is not centered using held-donor perturbation backgrounds;
- the five fixed pipelines and matched nulls are evaluated once on the outer
  panel after selection. Their scores are diagnostics, never a selector.

The outer PDS itself is scored on every held-donor gene except the 300 VCC
target genes. The extra 600 coordinates are merely left untouched, rather than
being removed from the outer score.

## Memory design

- Original donor matrices are opened with `mmap_mode="r"`.
- The union basis and covariance are temporary experiment-local memory maps.
- Covariance is accumulated in float32.
- The symmetric eigendecomposition uses LAPACK `evr` in float32.
- Only one outer fold and one full-rank operator exist at a time.
- Scratch files are removed after every fold.
- The run is serial and has an internal 3.5 GiB soft stop plus the external
  4 GiB process-tree guard.
- Primary and null transforms are materialized separately to avoid retaining
  unused full-panel copies.

No `cell-eval2`, submission command, or production artifact is used.

The first real-data attempt showed a transient process peak of about 3.96 GiB
on the wide Orion fold, even though the external one-second sampler observed
2.76 GiB. Therefore the 4 GiB command below remains the local fail-closed
command, but a complete run should use a larger cloud instance and a 10–12 GiB
guard. Do not weaken the local Mac guard.

```bash
uv run --frozen python scripts/run_with_memory_guard.py \
  --max-rss-gib 4 \
  --min-free-percent 25 \
  --poll-seconds 0.25 \
  -- \
  env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 \
  uv run --frozen python experiments/052-subspace-vs-union/run.py \
  --soft-rss-gib 3.5
```

```bash
uv run pytest -q \
  experiments/052-subspace-vs-union/test_subspace_vs_union.py
```

Cloud command:

```bash
uv run --frozen python scripts/run_with_memory_guard.py \
  --max-rss-gib 12 \
  --min-free-percent 10 \
  --poll-seconds 1 \
  -- \
  env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 \
  uv run --frozen python experiments/052-subspace-vs-union/run.py \
  --soft-rss-gib 10
```

For a synthetic-only audit under a strict 1 GiB process-tree limit:

```bash
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 1 \
  --min-free-percent 20 \
  --poll-seconds 0.1 \
  -- \
  env PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
  python3 -m pytest -q -p no:cacheprovider \
  --basetemp=experiments/052-subspace-vs-union/.pytest_tmp \
  experiments/052-subspace-vs-union/test_subspace_vs_union.py
```
