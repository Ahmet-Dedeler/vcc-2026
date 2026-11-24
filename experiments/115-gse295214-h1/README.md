# 115 — GSE295214 H1 metadata-first audit

Status: **metadata audit complete; author pseudobulk acquisition authorized;
raw GEO archive remains locked**

Date: **2026-08-30**

## Decision

Acquire the author-produced pseudobulk H5AD. Keep the raw GEO archive locked.

This is the strongest newly identified direct donor:

- same H1 hESC context as the existing VCC-2025 donor;
- CRISPRi rather than knockout or chemical perturbation;
- 85 released 10x H5 libraries and approximately 1.3 million profiled cells;
- 64 of the 300 VCC-2026 targets survive the authors' final pseudobulk QC;
- 57 of those 64 are absent from the existing 255-target VCC-2025 H1 donor.

The 4,425,260,319-byte pseudobulk is the minimum sufficient source because it
contains the post-QC perturbation effect, matched control profile, and measured
coordinate mask for every retained row. The 8,340,520,960-byte GEO archive is
not needed for this first donor experiment and remains unauthorized.

## Exact overlap result

Supplementary Table 6 contains 43,451 oligo-pool rows representing 6,638
unique sequences, including 370 unique non-targeting controls. Literal parsing
of the target labels yields 2,995 symbols; the paper reports 2,978 mapped genes.
That naming difference does not affect the VCC overlap.

The library contains 66 VCC-panel targets. The authors' public derived
pseudobulk database contains 2,354 post-QC H1 perturbations. Exactly 64 of the
66 VCC targets survive:

```text
ACLY, ADNP, ANKRD52, ANP32B, BRD1, BRPF1, C5orf22, COIL, CSNK2A1,
DCAF16, DHX35, DLAT, DNTTIP1, EAPP, EPC2, EZH2, GIGYF2, GMEB2,
GPATCH8, GTF2IRD1, HMGN3, HMGN4, HSP90AB1, KHDRBS1, KMT2A, MED15,
MED25, NABP2, PAN2, PAPOLA, PCGF1, PHF23, PI4KB, PIN1, PKN3, POR,
PPIF, PQBP1, PRMT6, PRNP, RSRC1, SCAF8, SFXN4, SIN3B, SLIRP,
SMARCA5, SNRK, SSR4, STT3A, TAF4, TEAD4, TIMM17B, TIMM50, TMEM41B,
TRPC4AP, UPF3B, UROS, USP28, ZNF219, ZNF22, ZNF397, ZNF462, ZNF670,
ZNF714
```

The two library-overlap targets absent after final QC are:

```text
PARP3, SYVN1
```

The new screen's retained axis overlaps 96 of the existing VCC-2025 H1 donor's
255 targets. Within the VCC-2026 panel, seven retained targets are already in
that donor (`ACLY`, `MED15`, `MED25`, `SIN3B`, `SLIRP`, `SMARCA5`,
`ZNF714`), leaving **57 genuinely new direct H1 rows**.

## Sources and provenance

- GEO: `GSE295214`, public May 5, 2026; last updated August 6, 2026.
- Paper: `10.64898/2026.05.16.725005`, posted May 18, 2026, CC BY 4.0.
- Analysis code: `xingjiepan/SCMG_scripts` at
  `aa848b9d517071d27b141e058e9bdae1966a869f` (MIT).
- SCMG code: `xingjiepan/SCMG` at
  `29f44c98c5621d575a058b549e0fdde0ba31a730` (MIT).
- Public derived data: `xingjiepan/SCMG_data` at
  `1d3f28d966772dcf75870beac052923efdc078c5`.

The exact assets, hashes, repository revisions, local axes, and expected audit
values are frozen in `SOURCE_MANIFEST.json`.

The post-QC target list was obtained from only
`/obs/dataset` and `/obs/perturbed_gene_name` in the authors' public
4,425,260,319-byte pseudobulk H5AD. A bounded HTTP range reader fetched 5 MiB;
it never indexed `/X` or `/layers`, and the full file was not downloaded.
`extract_post_qc_metadata.py` reproduces this operation and writes a receipt.

## Paper QC represented by the retained axis

The authors:

1. retained cells with at least 500 detected genes, exactly one sgRNA, and at
   least 20 sgRNA-barcode counts;
2. required at least 50 cells per perturbed gene;
3. required target expression below 90% of non-targeting controls and a
   one-sided Mann-Whitney U p-value below 0.01;
4. removed 259 perturbations with shift norm above 4 and cross-guide
   correlation below 0.2.

The 2,354-row `Pan2025_hESC` axis is the authors' released post-QC pseudobulk
set, not an inference from the 259 count.

## Safety boundary

This experiment may:

- inspect and hash small public metadata;
- read the two allowed H5AD `obs` columns through bounded HTTP ranges;
- compare target names with the public VCC-2026 panel and the local
  VCC-2025 H1 donor target list;
- run synthetic metadata tests.

It may not:

- download either multi-gigabyte source by default;
- read expression arrays from the public pseudobulk H5AD;
- open any VCC held perturbation truth;
- fit a model, run `cell-eval2`, generate counts, build a submission, or submit.

The original metadata manifest remains unchanged and denies both large assets.
`PSEUDOBULK_ACQUISITION_AMENDMENT.json` separately authorizes only the pinned
author-produced pseudobulk. The raw GEO archive remains denied. The downloader
requires the asset-specific confirmation token
`I_CONFIRM_GSE295214_PSEUDOBULK_DOWNLOAD`, validates the Hugging Face revision,
size, and LFS SHA-256 before transfer, resumes through a stable partial file,
uses an exclusive lock, and publishes the verified file plus receipt
atomically.

## Reproduce

```bash
uv run python experiments/115-gse295214-h1/extract_post_qc_metadata.py
python experiments/115-gse295214-h1/audit.py
python -m unittest experiments/115-gse295214-h1/test_audit.py
```

The first command requires `h5py` and `numpy` and performs metadata-only range
reads. The audit and tests use only the standard library.

Authorized pseudobulk acquisition:

```bash
uv run python experiments/115-gse295214-h1/acquire.py \
  --asset public_pseudobulk_database \
  --confirm-large-download I_CONFIRM_GSE295214_PSEUDOBULK_DOWNLOAD
```

After acquisition has atomically published both the exact H5AD and its
download receipt, build the compact source donor:

```bash
uv run python experiments/115-gse295214-h1/extract_compact_donor.py
```

The extractor fails closed unless the source manifest, acquisition amendment,
download receipt, 4,425,260,319-byte H5AD, H5AD SHA-256, frozen 2,354-target
axis, 300-target panel, and canonical 18,533-gene axis all authenticate. It
then validates dense float32 `X`, `layers/control`, and
`layers/measure_mask`; reconstructs perturbed and control CPM10k in linear
space; collapses duplicate symbols there; and converts effects to project
TS=50,000 log1p space without renormalizing away unmapped source genes.

The default output is `work/compact-donor/`. It contains donor-compatible
`delta.npy`, `b_ctrl.npy`, `genes.txt`, and `perts.txt`, plus:

- `measured_mask.npy`: row-by-canonical-coordinate support;
- `source_covered_gene_mask.npy`: canonical genes present in the source;
- `vcc_target_mask.npy` and `vcc_targets.txt`: the exact 64 VCC rows;
- row/source mapping, duplicate-symbol mapping, conversion metadata;
- a deterministic artifact manifest and a separately SHA-256-pinned receipt.

Verify an existing compact donor without reopening the 4.4 GB source:

```bash
uv run python experiments/115-gse295214-h1/extract_compact_donor.py \
  --verify-only
```

Synthetic tests do not require the downloaded H5AD:

```bash
uv run python -m unittest \
  experiments/115-gse295214-h1/test_extract_compact_donor.py
```
