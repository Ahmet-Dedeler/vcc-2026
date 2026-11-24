# 048 — Source-only low-rank geometry

Experiment 046 found a large retrieval gain from a low-rank basis learned from
the held-out target donor, but permutation controls showed that perturbation
pairing was irrelevant. The only defensible follow-up is to ask whether the
same kind of geometry can be learned entirely from source donors.

This experiment uses strict nested leave-one-donor-out evaluation over:

- `k562_gwps`
- `orion_hct116`
- `orion_hek293t`

For each outer donor:

1. The outer donor is unavailable during fitting and hyperparameter selection.
2. The two source donors are evaluated in both inner directions.
3. Each inner direction is split into five perturbation folds. The PCA basis is
   fit on one source donor's non-VCC training perturbations and scored against
   the other source donor on held-out non-VCC perturbations.
4. A source-only candidate is selected from those inner scores.
5. The selected transform is refit on both source donors' 600 non-VCC
   perturbations and evaluated once against the outer donor on the 259 shared
   VCC targets.

The source preprocessing is also leakage-safe: donor centering and per-gene
scale are recomputed by streaming each original `delta.npy` over background
perturbations that exclude both the VCC panel and the complete 600-perturbation
selection panel. The original arrays are memory-mapped.

Comparators and controls:

- Identity (untransformed equal-source consensus).
- A fixed `gamma=0.65`, rank-256 source-covariance transform: a bounded,
  low-rank analogue of the current donor whitener.
- Random-subspace controls with matched rank and attenuation.
- A perturbation-label permutation audit. Since the source transform is
  covariance-only, row permutation must leave it unchanged up to numerical
  precision.

The bounded source grid uses ranks 8–384, spectral powers 0.25–8.0, and simple
uniform suppression of the same principal subspaces. Rank 384 stays below the
480-row inner-fit limit. Source PCs are exact via the small sample-space Gram
matrix, so no randomized source-basis instability or full gene covariance
decomposition is involved.

No `cell-eval2`, count generation, submission build, or production artifact is
used.

Run:

```bash
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 4 \
  --min-free-percent 20 \
  -- \
  uv run python experiments/048-source-only-subspace/run.py
```

Tests:

```bash
uv run pytest -q \
  experiments/048-source-only-subspace/test_source_only_subspace.py
```
