# 119 — Altair A549 source acquisition and compact source

Status: **official source authenticated; source semantics verified; compact
source built**

Date: **2026-08-30**

## Scope

This experiment is intentionally limited to source acquisition, audit, and
lossless compacting of the official perturbation-level Altair CROP-seq H5AD.
It does not open VCC held perturbation truth or existing outcome directories,
fit a model, evaluate predictions, run `cell-eval2`, emit counts, build a
submission, or submit.

The already-downloaded official source is:

```text
data/external/altair-a549/altair-scrna.h5ad
```

It is authenticated against:

- size: `281,474,102` bytes;
- SHA-256:
  `4632972e6c2871603f4ca9108489787d7761aae05dc6a4bf3488f1ed78a63c9f`;
- S3 ETag: `cb152c5b6fb206a010d327174e224493-34`;
- S3 version:
  `faMp5Esq6WYP5QbbKEwXEbneXEMTyMyk`;
- last modified: `2026-08-25T17:01:21Z`.

All source, API, axis, and expected-audit bindings are frozen in
`SOURCE_MANIFEST.json`.

## Source semantics

The H5AD is a perturbation-level differential-expression object, not a
single-cell count matrix:

- rows: 1,000 gene-targeting aggregates plus one pooled non-targeting control;
- columns: 9,999 Ensembl genes with unique gene symbols;
- `.X`: `float32` log2 fold change;
- layers: `baseMean`, `lfcSE`, `stat`, `p_value`, `adj_p_value`, and
  `neg_log10_fdr`.

The `.X` interpretation is verified independently:

1. the official deployed OPS API returns the same row under the field
   `log2_fold_change`;
2. the API `ALG12` vector and `.X['ALG12']` agree within the API's five-decimal
   serialization precision;
3. source values satisfy `.X / lfcSE == stat` within floating-point error.

Finite `.X` is the measured-effect definition. Adjusted-p-value availability
is separate: DESeq2-style filtering leaves many `adj_p_value` and
`neg_log10_fdr` entries as NaN even when the LFC is measured.

Four rows are entirely unsupported in the released differential-expression
matrix: `COQ10B`, `HSD17B12`, `NEK1`, and the pooled `Nontargeting` row.
Therefore 997 of 1,000 targeting rows have source effects. The source contains
21 VCC-panel target labels, but only 20 have measured effects because
`COQ10B` is all-NaN.

## Exact-symbol projection

The compact artifact uses the canonical 18,533-gene VCC axis from
`data/challenge/gene_names.csv`.

- 7,843 source response symbols map by exact string equality;
- 2,156 source symbols are not on the canonical axis and remain documented as
  unmapped;
- 10,690 canonical symbols are absent from the source;
- aliases are never resolved or merged;
- unsupported output coordinates remain NaN;
- `measured_mask.npy` records finite source LFC support explicitly;
- `fdr_mask.npy` separately records finite source FDR support.

The artifact retains all 1,000 targeting rows in official source order. It
does not convert log2 fold changes into project-space effects and must not be
treated as a ready-made donor `delta.npy`.

## Reproduce

Verify the source, official API receipts, schema, and semantic identities:

```bash
uv run python experiments/119-altair-a549/source_audit.py audit --write-report
```

Run tests:

```bash
uv run python -m unittest \
  experiments/119-altair-a549/test_source_audit.py \
  experiments/119-altair-a549/test_extract_compact_source.py
```

Build the source-only artifact:

```bash
uv run python experiments/119-altair-a549/extract_compact_source.py build
```

Verify an existing artifact without reopening the public H5AD:

```bash
uv run python experiments/119-altair-a549/extract_compact_source.py verify
```

Default output:

```text
experiments/119-altair-a549/work/compact-source/
```

The output is immutable-by-default: the builder refuses to overwrite an
existing directory. It contains the canonical LFC and FDR matrices, explicit
support masks, exact source-to-canonical mapping, target metadata, a complete
artifact manifest, and a receipt binding the artifact to the source,
experiment manifest, audit, and implementation.

