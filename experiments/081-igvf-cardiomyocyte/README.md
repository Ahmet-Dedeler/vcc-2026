# 081 — IGVF WTC11 cardiomyocyte sparse-direction transfer

This experiment asks one narrow question: can the released compact IGVF
WTC11-derived cardiomyocyte CRISPRi differential-expression directions improve
the frozen production-shaped, source-only union-whitened predictor when K562,
HCT116, and HEK293T are held out in turn?

It is deliberately bounded:

- the only downloaded payloads are `IGVFFI0830FXFI.csv.gz` (filtered global
  differential expression) and `IGVFFI8270UPKB.csv.gz` (guide manifest);
- it never reads or downloads a large IGVF H5AD/H5MU;
- it never fits a whitener or runs `cell-eval2`;
- it uses the already frozen experiment-062 gamma-0.65 held-donor operators;
- all writes stay in this directory or
  `data/external/igvf-cardiomyocyte/`;
- its process-tree RSS is externally capped below 1 GiB.

## Leakage-safe protocol

1. Verify the official IGVF release status, byte sizes, and MD5 checksums
   recorded in `PROVENANCE.json`.
2. Parse the DE CSV and tab-delimited guide manifest, aggregate duplicate
   target-response pairs by cell-count-weighted mean log2 fold change, and
   reject malformed or non-finite rows.
3. Require at least 1,000 IGVF targets with an official-axis response shared
   by all three held-truth banks outside the VCC panel. Deterministically
   hash-split those targets 60/20/20 into fit/tune/test sets.
4. Use fit targets only to freeze robust clipping and two source-side scale
   rules. Use tune targets only to select:
   - representation: clipped log2FC, signed square-root log2FC, or sign;
   - scale: global L2 matching or capped per-target L2 matching;
   - residual fraction: 0, 0.01, 0.02, 0.05, 0.10, 0.20, or 0.35.
5. Evaluate the frozen choice once on the untouched non-panel test split.
6. Only after selection and the non-panel test are complete, load VCC-panel
   effects once for the final evaluation.

Every held fold excludes that donor from prediction construction. The
corresponding frozen union operator was fitted from the other two donors only.
The full 300-gene VCC panel is absent from every operator axis and retrieval
axis. For the non-panel proxy task, all 1,525 proxy target coordinates are
also zeroed before transformation and excluded from retrieval, preventing
on-target knockdown from selecting the adapter.

The IGVF contribution is added after incumbent whitening and only at
target-response pairs actually present in the released filtered DE table.
Missing targets and missing genes therefore retain the incumbent exactly.
Every split uses its complete, predeclared gallery. The final panel gallery is
259 targets; 81 retain at least one usable IGVF edge after excluding the 300
VCC target coordinates, and the other 178 remain exact incumbent fallbacks.

Uncertainty uses a target-cluster bootstrap. Clusters are derived from
candidate-independent incumbent prediction similarity, never from held-truth
outcomes. Two frozen nulls are reported:

- an exact-sparsity sign null, which randomizes the sign of every released
  target-response edge while preserving targets, genes, magnitudes, and nnz;
- a matched label null, which permutes whole IGVF rows within target sparsity
  and norm strata.

Promotion requires nonnegative PDS delta in every held fold, macro delta at
least `+0.003`, a positive target-cluster bootstrap lower 95% bound, and an
observed delta above both null 95th percentiles on both the untouched
non-panel test and final complete VCC-panel evaluation.

## Run

```bash
PYTHONDONTWRITEBYTECODE=1 \
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
VECLIB_MAXIMUM_THREADS=2 \
.venv/bin/python scripts/run_with_memory_guard.py \
  --max-rss-gib 0.98 --min-free-percent 25 --poll-seconds 1 \
  -- .venv/bin/python experiments/081-igvf-cardiomyocyte/run.py \
     --guarded --max-rss-gib 0.90
```

Tests:

```bash
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest \
  -p no:cacheprovider experiments/081-igvf-cardiomyocyte/test_run.py
```
