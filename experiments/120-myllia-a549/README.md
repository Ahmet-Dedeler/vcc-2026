# 120 — Myllia A549 low-MOI compact donor

Status: **source acquisition authenticated; source-only compact donor built**

Date: **2026-08-30**

## Scope

This experiment authenticates and converts the two downloaded A549 channels
from `GSE337804` into a compact, canonical-axis source donor. It is strictly
source-only:

- no VCC held perturbation truth or existing outcome directory is opened;
- no model is evaluated and `cell-eval2` is never invoked;
- no predicted cells, challenge package, or submission is produced;
- no commit or push is performed.

The four source files, canonical 18,533-gene axis, and 300-target panel are
content-addressed in `SOURCE_MANIFEST.json`. `SOURCE_PROVENANCE.json` records
the assay and file semantics from the supplied description and library table.

## Source

The source is a low-MOI CRISPR knockout Perturb-seq screen in monoclonal A549
cells. Two Chromium X channels were processed independently. Retained cells
have exactly one called sgRNA in `obs['sgrna_id']`, with the target or control
class in `obs['sgrna_symbol']`.

The two H5ADs contain:

| Channel | Cells | Genes | NTC cells | Gene-desert cells |
| --- | ---: | ---: | ---: | ---: |
| `ch_1` | 25,702 | 18,711 | 1,033 | 988 |
| `ch_2` | 27,809 | 18,808 | 1,060 | 1,075 |
| Total | **53,511** | — | **2,093** | **2,063** |

Both `X` and `raw.X` are raw, non-negative integer-valued CSR matrices.
Extraction uses `raw.X` and verifies that `X` is logically identical.

The sgRNA library has 694 rows: 654 targeting guides for 218 genes, 20
non-targeting guides, and 20 gene-desert guides. Ten targeting guides are not
observed in either retained-cell object.

## Exact VCC overlap

The library and observed target labels have 15 exact-symbol overlaps with the
VCC-2026 panel:

```text
BRD1, BRPF1, CDK8, E2F1, HDAC1, KDM5A, KMT2A, MTM1, PCGF1, PCGF6,
PSMB9, SIN3B, SMARCA5, SYVN1, UPF3B
```

No aliases, prefixes, Ensembl remapping, case folding, or punctuation changes
are allowed.

## Frozen extraction and QC

For each channel, raw counts are summed independently for every observed
targeting guide and for the pooled non-targeting controls. Each pseudobulk is
normalized over its complete source feature axis to 50,000 counts and then
mapped to the canonical axis by exact symbol.

The guide-channel effect is:

```text
log1p(TS50000 guide pseudobulk) - log1p(TS50000 channel NTC pseudobulk)
```

The frozen structural QC is deliberately outcome-independent:

1. a guide-channel row needs at least 10 retained cells;
2. a guide is eligible only when the same exact guide passes in both channels;
3. the two channel effects are averaged with equal weight;
4. a target is eligible only with at least two eligible guide sequences;
5. the target effect is the coordinate-wise median of its eligible guides.

This retains 491 replicated guides and 175 of 218 targets. Fourteen of the 15
VCC-overlap targets pass. `KMT2A` fails because only `KMT2A_2` has at least 10
cells in both channels.

The build records channel-level and guide-level replication diagnostics but
does not use those biological outcomes to select targets. Cell-count and
replication structure alone determine the frozen measured masks.

## Panel-coordinate exclusion

All 300 VCC panel coordinates are excluded from every transferable response:

- they are exact zero in `guide_channel_delta.npy`, `guide_delta.npy`, and
  target-level `delta.npy`;
- they are false in every measured mask;
- `transferable_gene_mask.npy` is the intersection of exact canonical coverage
  in both channels, minus the full 300-gene panel.

The source covers 16,777 canonical genes in both channels, leaving 16,477
transferable non-panel coordinates.

## Compact artifacts

`work/compact-donor/` contains:

- `delta.npy`, `b_ctrl.npy`, `genes.txt`, `perts.txt`: target-level donor;
- `measured_mask.npy`, `target_qc_mask.npy`: frozen target support/QC;
- `source_targets.txt` and `vcc_design_overlap.txt`: axes for all source
  targets and all 15 exact panel overlaps, including QC failures;
- `guide_delta.npy`, `guide_rows.tsv`, `guide_measured_mask.npy`: replicated
  guide effects after equal channel pooling;
- `guide_channel_delta.npy`, `guide_channel_rows.tsv`, and
  `guide_channel_measured_mask.npy`: qualifying channel-specific guide effects;
- `channel_b_ctrl.npy`, channel/source coverage masks, and
  `transferable_gene_mask.npy`;
- complete `guide_qc.tsv` and `target_qc.tsv`, including rejected rows;
- `meta.json`, `artifact-manifest.json`, and a separately hashed receipt.

All arrays use the canonical gene order from
`data/challenge/gene_names.csv`. Unsupported coordinates are exact zero.

## Reproduce

Authenticate and audit the downloaded source without writing donor arrays:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/120-myllia-a549/build_compact_donor.py audit --write-report
```

Run the synthetic/fail-closed tests:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python -m unittest \
  experiments/120-myllia-a549/test_build_compact_donor.py
```

Build once into the default immutable output:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/120-myllia-a549/build_compact_donor.py build
```

Verify without reopening either large H5AD:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run python \
  experiments/120-myllia-a549/build_compact_donor.py verify
```
