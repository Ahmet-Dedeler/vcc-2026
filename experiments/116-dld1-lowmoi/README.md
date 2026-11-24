# 116 — GSE337988 DLD-1 low-MOI source audit

Status: **source-only compact artifact built and verified**

Date: **2026-08-30**

## Decision

The narrow acquisition amendment authorized only the two author-produced
low-MOI LFC matrices:

| Role | Official file | Exact bytes |
| --- | --- | ---: |
| Low1 reference split | `GSE337988_sublib2_de_matrices_lfc_matrix_Low1.csv.gz` | 155,201,822 |
| Low2 holdout split | `GSE337988_sublib2_de_matrices_lfc_matrix_Low2.csv.gz` | 155,251,513 |

NCBI exposes stable HTTPS URLs, exact byte counts, and a `Last-Modified` value
of `Tue, 30 Jun 2026 21:16:28 GMT` for both files, but no checksum or ETag.
Both objects were therefore reproduced by two distinct complete transfers,
gzip-validated, SHA-256 matched, and compared byte-for-byte before atomic
acceptance. Accepted hashes:

- Low1: `aaa70afe2dfefd22b72fe6294b87887bc6ee928590e419f459bcdf9616fba517`
- Low2: `4e9bd64d8a3223dcb50d9c2a5b844698099a81dc2f86c38e4c0e270b4c87c390`

The 26,495,651,840-byte `GSE337988_RAW.tar` archive is explicitly denied.
Every other processed asset is denied by omission from the two-item allowlist.

The source-only artifact is complete at `work/compact-source-v1/`. Its
67-target Low1/Low2 retrieval macro PDS is 0.9471, but the exact preregistered
held-context promotion gate is a NO-GO because broad canonical response
coverage and median split cosine fail. See `RESULTS.md`.

## Why this source matters

The scaled experiment is a direct CRISPRi Perturb-seq screen in DLD-1, a
colorectal adenocarcinoma line. It uses an inducible Zim3-dCas9 system and a
6,455-guide library. The low-MOI arm is approximately MOI 0.3 and spans 24
10x Chromium channels. The paper reports 1,869,846 cells after QC across the
scaled screen and 1,432,753 cells with at least one assigned guide.

Low1 and Low2 are not different cell lines or treatments. The authors split
the low-MOI cells in half: Low1 is the reference split and Low2 is a
same-condition holdout/noise-floor split. That pairing makes the two matrices
especially useful for estimating perturbation reliability before attempting
cross-line transfer.

## Assay and matrix semantics

The official GEO record describes:

- DLD-1 cells with `TET_DD_inducible_Zim3dCas9`;
- 250 ng/mL doxycycline plus Shield-1 at 1:500 for five days;
- 10x Chromium Single Cell 3' GEM-X v4 GEX, sgRNA capture, and MULTI-seq CMO
  hashing libraries;
- Cell Ranger 9.0.1 against GRCh38-2020-A with intron counting;
- guide assignment, hashing demultiplexing, and low-quality-cell removal;
- limma-voom differential expression on pseudobulk counts with
  perturbation-matched controls.

The candidate CSVs are genes by promoter-level perturbation labels. Their
values are author-produced limma-voom log-fold-change estimates, not raw
counts, single-cell observations, or a ready-made project-space delta. Any
future conversion into the VCC donor representation needs a separate,
preregistered source-only contract.

The Low1 and Low2 headers are byte-identical:

- 6,451 unique target labels;
- 62,617 bytes in the complete first CSV row;
- header SHA-256
  `de441b374d2f020f67993fa36922e2195fb34e8de82e885c5091b601a0a878a4`;
- newline-delimited target-label SHA-256
  `d49fb70eb8a6c9c39f5ff72a86ca71fcf47cb66a507b8d1b1d5ed7804b4d3a0f`.

The header includes 496 `NTC_*` controls and 5,955 non-NTC promoter labels
covering 4,954 gene prefixes. The paper's downstream analyses further remove
controls, `TRPC1_P1`, and columns without complete signal; this source audit
does not inspect matrix values and therefore does not claim a post-effect-QC
row count.

## Exact VCC overlap

The overlap is computed directly from the official Low1/Low2 CSV labels, not
from prose counts or a guessed gene list. A VCC gene matches only an identical
label or a label beginning with `<GENE>_`; this avoids substring collisions.

Exactly **67 of 300 VCC target genes** have at least one official LFC column,
represented by **71 promoter labels**:

```text
AGO1, AKT2, ALG12, API5, BCL6, BRPF1, C5orf22, CAPRIN2, CARM1, CHD3,
COIL, CYTH2, DCP2, DEK, DTNBP1, ELK1, EPC2, EZH2, FCGRT, GTF2IRD1,
HMGXB3, KMT2A, LDLR, LENG1, LPCAT3, MAP2K2, MAPK7, MSANTD4, MTF1,
MYPOP, NAIF1, PAN2, PARP3, PBLD, PHF19, PHF23, PIN1, PPIF, PRMT7,
PRNP, PSKH1, RABEP2, RBM26, RFK, RNF141, SEMA4F, SH2D3A, SMARCA5,
STT3A, TARBP2, THAP7, TRIM35, TTBK2, UPF3B, UROS, ZBTB34, ZC3HAV1,
ZHX3, ZKSCAN2, ZKSCAN3, ZNF219, ZNF25, ZNF251, ZNF324, ZNF410,
ZNF821, ZRANB3
```

Four VCC genes have two promoter-level columns:

```text
DEK: DEK_P1, DEK_P2
LDLR: LDLR_P1, LDLR_P2
PBLD: PBLD_P1, PBLD_P2
SMARCA5: SMARCA5_P1, SMARCA5_P2
```

This is label coverage only. No LFC values were opened, so this experiment
does not assert that all 71 columns contain finite or sufficiently reliable
effects.

## Licensing

- The bioRxiv preprint is CC BY 4.0.
- The pinned `Genentech/perturbmatch_analysis` repository licenses code under
  MIT and its `data/` directory under CC BY 4.0.
- The pinned `Genentech/perturbmatch` R package declares MIT licensing.
- GEO does not expose a separate per-file license field for these assets.
  Reuse should retain attribution to the paper, authors, GEO accession, and
  the author repository's data-license statement.

## Historical metadata-audit safety boundary

Before the acquisition amendment, this experiment could:

- download and verify the three small metadata files in `metadata/`;
- issue `HEAD` requests for the two candidate LFC matrices and denied raw
  archive;
- read at most 131,072 compressed bytes from each candidate LFC matrix to
  recover its first CSV row;
- compare those public labels with `data/challenge/pert_counts.csv`;
- write audit artifacts only inside this experiment;
- run synthetic and offline metadata tests.

At that stage it could not:

- download either LFC matrix in full;
- download the raw archive or any other processed object;
- read any matrix values beyond the first-row metadata probe;
- open held perturbation truth;
- fit or evaluate a model, run `cell-eval2`, generate counts, submit, commit,
  or push.

`source_audit.py` remains metadata-only. The later acquisition and conversion
are governed by `ACQUISITION_AMENDMENT.v1.json`,
`PREVALUE_FREEZE.v1.json`, and `CONVERSION_PROTOCOL.v1.json`. They do not
authorize the raw archive, held truth, outcomes, scoring, counts, submission,
commit, or push.

## Reproduce

Acquire only the three small, hash-pinned metadata files:

```bash
uv run python experiments/116-dld1-lowmoi/source_audit.py acquire-metadata
```

Re-run the live `HEAD` and bounded-header audit:

```bash
uv run python experiments/116-dld1-lowmoi/source_audit.py audit --write-report
```

Reproduce the audit without network access:

```bash
uv run python experiments/116-dld1-lowmoi/source_audit.py audit --offline
```

Run synthetic and frozen-metadata tests:

```bash
uv run python -m unittest \
  experiments/116-dld1-lowmoi/test_source_audit.py
```

The complete acquisition, conversion, verification, diagnostics, and hashes
are recorded in `RESULTS.md`.
