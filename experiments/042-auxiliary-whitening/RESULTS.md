# 042 — Auxiliary screens do not improve the whitening covariance

**Date:** 2026-08-22. Two Orion leave-one-donor-out folds, 259 shared VCC-panel
perturbations, 17,806 scored genes. The whitening operator was screened on a
fixed random 3,000-gene subset of the donor intersection to keep peak RSS below
8 GiB. Prediction effects still came only from the production donors; auxiliary
screens contributed covariance rows only.

## Result

Best gamma shown for each basis:

| covariance basis | best gamma | macro PDS | delta vs base |
|---|---:|---:|---:|
| **production donors, row pooled** | **0.85** | **0.70538** | — |
| + Nadig Jurkat | 0.85 | 0.70527 | -0.00011 |
| + H1-2025 | 0.85 | 0.70499 | -0.00039 |
| + K562-essential + RPE1 | 0.85 | 0.70484 | -0.00055 |
| + all auxiliary except H1 | 0.85 | 0.70466 | -0.00073 |
| production donors, equal source mass | 0.85 | 0.70466 | -0.00073 |
| + all auxiliary, row pooled | 0.85 | 0.70444 | -0.00094 |
| + all auxiliary, equal source mass | 0.85 | 0.70040 | -0.00498 |

Every auxiliary variant was worse than the base on the macro average. Jurkat
was effectively tied, but its change is far below the project's resolution
threshold and had opposite fold-level signs too small to justify a full-axis
run.

## Decision

Do not add H1, Replogle essential, RPE1, or Nadig Jurkat as whitening-only rows
to v9. The useful covariance appears to be panel/task-distribution specific;
more rows from mismatched perturbation panels smooth the operator rather than
improve retrieval.

This does **not** rule out using a genuinely matched screen as a prediction
donor. It rules out the cheaper idea that unrelated perturbations from extra
screens improve the global covariance by sample size alone.

## Memory correction

The first implementation loaded and duplicated complete Orion matrices before
selecting columns. It reached 15.2 GiB RSS and was killed by the new memory
guard. The final implementation memory-maps `delta.npy`, centers selected
columns in chunks, and peaked at **6.75 GiB** with system free memory at 85–86%.

Reproduce:

```bash
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 VECLIB_MAXIMUM_THREADS=4 \
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 8 --min-free-percent 35 --poll-seconds 2 -- \
  uv run python experiments/042-auxiliary-whitening/aux_whiten.py \
  --axis inter --max-genes 3000 \
  --gammas 0.45 0.55 0.65 0.75 0.85 \
  --variants base_rows base_equal_source plus_h1_rows \
  plus_essential_rows plus_jurkat_rows plus_all_no_h1_rows \
  plus_all_rows all_equal_source
```
