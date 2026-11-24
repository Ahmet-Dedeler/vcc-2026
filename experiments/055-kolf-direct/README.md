# 055 — Direct KOLF2.1J donor transfer

## Goal

Test the official KOLF2.1J genome-scale CRISPRi atlas as one additional
perturbation-response donor. This is a bounded effect-level experiment:

- fetch only the 273 VCC-panel perturbation rows exposed by the official web
  atlas, not the 392,507,202-byte full matrix;
- map those rows onto the official VCC gene axis;
- choose a fixed KOLF blend fraction from `0.1`, `0.25`, and `0.5` using only
  held HCT116 and HEK293T folds;
- after freezing that choice, evaluate exactly one K562 target-role holdout;
- compare raw and mean-absolute gene-normalized effects;
- exclude all 300 VCC target coordinates before centering, normalization, or
  retrieval; and
- rank every query against the complete shared-target gallery.

It does not run `cell-eval2`, emit counts, build a submission, or alter the
production model.

## Verified official binary contract

The source is pinned to official repository commit:

```text
8f84531a784f24b14d7e2c9b80c375d7e78fd46c
```

Before decoding rows, `extract_kolf.py` downloads and SHA-256 verifies the
official metadata, matrix builder, browser loader, and license. The builder
states and implements:

```text
source: zscore_filt_log1p.parquet
matrix: perturbations × genes
dtype: int8
layout: C-order row-major
scale: 0.04
quantization: clip(round(value / 0.04), -127, 127)
```

The browser loader independently computes:

```text
chunk = floor(perturbation_index / 2500)
offset = (perturbation_index % 2500) * 38606
Range = bytes offset .. offset + 38606 - 1
decoded value = int8 * 0.04
```

The metadata has 10,167 unique perturbations, 38,606 unique genes, five
chunks, and exactly one byte per gene. The extractor rejects any response
other than an exact HTTP `206` with the expected `Content-Range`; it does not
read a full-chunk `200` fallback.

## Run

Both commands have an external process-tree cap below the requested 4 GiB
ceiling, as well as an internal 3.5 GiB peak-RSS check.

```bash
uv run --frozen python scripts/run_with_memory_guard.py \
  --max-rss-gib 3.9 \
  --min-free-percent 25 \
  --poll-seconds 0.25 \
  -- \
  env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 \
  uv run --frozen python experiments/055-kolf-direct/extract_kolf.py
```

```bash
uv run --frozen python scripts/run_with_memory_guard.py \
  --max-rss-gib 3.9 \
  --min-free-percent 25 \
  --poll-seconds 0.25 \
  -- \
  env OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 \
  uv run --frozen python experiments/055-kolf-direct/evaluate_lodo.py
```

Tests:

```bash
uv run --frozen pytest -q experiments/055-kolf-direct/test_kolf_direct.py
uv run --frozen ruff format --check experiments/055-kolf-direct
uv run --frozen ruff check experiments/055-kolf-direct
```

## Outputs

`data/external/kolf2/`:

- `source_rows_int8.npy`: the exact 273 fetched source rows on all 38,606
  source genes;
- `delta.npy`: decoded float32 effects on the exact-overlap VCC gene axis;
- `perts.txt`, `genes.txt`, and compatibility-only `b_ctrl.npy`;
- `target_row_map.csv` and `gene_map.csv`;
- `meta.json` and `provenance.json`; and
- `official/`: byte-for-byte pinned metadata, builder, loader, and license.

This directory contains no full source chunk.

Experiment outputs:

- `results.json`: protocol, all development candidates, frozen selection, and
  the one K562 holdout result;
- `results.csv`: HCT116/HEK293T development grid only;
- `RESULTS.md`: interpreted result and decision.
