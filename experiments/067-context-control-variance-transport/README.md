# 067 — Context-control variance transport

Test the public leaderboard disclosure

```text
effect[g] *= sd_target_control[g] / sd_source_control[g]
```

as a distinct effect-level transfer method. This is not experiment 010's
basal-mean donor weighting: donor weights remain fixed and only the per-gene
shape of an already-built effect vector changes.

## Protocol

- Outer donors: K562-GWPS, Orion HCT116, and Orion HEK293T.
- Outer panel: the 259 VCC targets measured by all three donors.
- Every outer fold uses the fixed gamma-0.65 union whitener from experiment
  062 that was fit from the other two donors only.
- All 300 VCC target coordinates remain excluded from the whitening and
  retrieval axes.
- Held-donor perturbation effects are unavailable during selection.
- Held-donor control cells are allowed because the challenge exposes target
  controls at inference.
- Inner selection predicts each of the two source donors from the other on a
  deterministic 96-target subset.
- Control-scale estimators:
  - `sd`: 1%/99% winsorized standard deviation of per-cell `log1p(CPM)`;
  - `mad`: `1.4826 * median(abs(x - median(x)))`.
- Ratio:

  ```text
  r[g] = clip(scale_target[g] / geometric_mean(scale_source[g]), 0.5, 2.0)
  effect[g] *= r[g] ** eta
  ```

- Powers: `eta = 0, 0.25, 0.5, 1.0`.
- Placement: immediately before or immediately after the fixed whitener.
- `eta=0` is represented once as the exact identity baseline.

Two disjoint, deterministic control samples are extracted per donor. Replicate
A drives source-only selection; replicate B is a frozen stability evaluation.
The selected ratio is also gene-permuted 100 times on the outer fold.

## Resource contract

All commands must run through the process-tree guard:

```bash
python scripts/run_with_memory_guard.py \
  --max-rss-gib 0.70 --min-free-percent 20 --poll-seconds 0.05 -- \
  uv run python experiments/067-context-control-variance-transport/run.py \
  --stage orion-workers --max-rss-gib 0.68

python scripts/run_with_memory_guard.py \
  --max-rss-gib 0.70 --min-free-percent 20 --poll-seconds 0.05 -- \
  uv run python experiments/067-context-control-variance-transport/run.py \
  --stage controls --max-rss-gib 0.68

python scripts/run_with_memory_guard.py \
  --max-rss-gib 0.70 --min-free-percent 20 --poll-seconds 0.05 -- \
  uv run python experiments/067-context-control-variance-transport/run.py \
  --stage evaluate --max-rss-gib 0.68

python scripts/run_with_memory_guard.py \
  --max-rss-gib 0.70 --min-free-percent 20 --poll-seconds 0.05 -- \
  pytest -q \
  experiments/067-context-control-variance-transport/test_variance_transport.py
```

Orion's Parquet files each contain one 370–490 MB row group. The first stage
therefore extracts one file per fresh child process, never keeping two decoded
nested columns resident together. The second stage consumes those tiny local
row caches and removes them after writing `control_stats.npz`.

The large fold-specific `.npz` eigenvector arrays are ZIP-stored rather than
compressed. `run.py` reads their payload directly in row blocks; it never maps
or loads a complete 1.2 GB operator into RAM. Donor `delta.npy` tables use the
same bounded direct-row reader because macOS retained pages from ordinary
memory maps.
