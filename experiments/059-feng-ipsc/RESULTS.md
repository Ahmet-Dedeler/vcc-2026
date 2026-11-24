# 059 — Feng 34-line iPSC atlas audit and direct-donor screen

**Date:** 2026-08-22  
**Scope:** official-source audit, bounded-memory extraction, and raw-space
effect-level PDS screening only. No cell-level expression download, AWS use,
`cell-eval2`, count emission, or submission.

## Decision

**The dataset is real and directly usable, but the tested Feng blend is a
practical no-go.**

The earlier repository note missed the most useful release:
`GenomeWideScreen_LFC_byGene.tsv.gz`, an official pooled target-level effect
table covering both the fitness and non-fitness genome-wide screens. It yields
a clean compact artifact for **182/300 VCC targets** and **6,151/18,533 VCC
genes**.

The predeclared raw-space LODO gate technically passed. The selected candidate
was a **2% Feng blend limited to targets with significant on-target
downregulation**:

- HCT116 delta: **+0.000014965 raw PDS**
- HEK293T delta: **+0.000074826**
- development mean: **+0.000044895**
- frozen K562 delta: **+0.000119721**
- three-fold macro delta: **+0.000069837**

That gain is only a few discrete retrieval-rank steps over the entire
evaluation and is far below a practically resolvable effect. Larger blend
fractions were consistently harmful. Therefore:

- **GO:** retain the verified compact matrix as a research artifact.
- **NO-GO:** do not add Feng to v9/v10, do not refit whitening with it, and do
  not spend AWS/scorer budget on it now.
- **Conditional next step:** only if a future method supplies a stronger,
  independently justified confidence gate, run exactly one fixed-clean-
  whitener HCT check of `self_sig_0_1`, fraction `0.02`, preserving the
  baseline for missing targets and genes. The present result does not warrant
  that escalation.

## Authoritative provenance

### Publication

- Claudia Feng et al., *A genome-scale single-cell CRISPRi map of trans gene
  regulation across human pluripotent stem cell lines*.
- Cell Genomics 6(2):101076; online publication **2025-12-01**.
- DOI: `10.1016/j.xgen.2025.101076`
- PMID: `41330380`; PMCID: `PMC12903452`.
- The article is CC BY 4.0.

The paper reports:

- **7,226 genes** targeted with three Dolcetto guides per gene plus 40
  non-targeting guides;
- **34 iPSC lines from 26 healthy donors**;
- **219,206** cells assigned to a targeting guide and **499,998** controls;
- median **8 cells/guide** and **25 cells/target**;
- **6,673** targets passing the minimum-10-cell effect threshold;
- effects computed across **6,471 biological expressed genes**.

### Official processed data: Figshare

Transcriptional Changes:

- article `26819743`, version 1
- DOI: `10.6084/m9.figshare.26819743.v1`
- license: **MIT**
- published 2025-10-08

| File | Figshare ID | Bytes | MD5 |
|---|---:|---:|---|
| `GenomeWideScreen_LFC_byGene.tsv.gz` | **49291654** | **794,610,300** | `03cb8a3b4328de50d45746a9d784bbc5` |
| `GenomeWideScreen_LFC_byGuide.tsv.gz` | 49291672 | 1,812,579,303 | `adb0418b44f95866331b5c72259a26a1` |
| `GenomeWideScreen_Coregulation-Target-Correlation.tsv.gz` | 49291633 | 155,213,386 | `1b0183ed08a415f019aba59d338651d9` |
| `GenomeWideScreen_Coperturbation-Expressed-Gene-Correlation.tsv.gz` | 49291645 | 529,575,408 | `236092cfd03c78bcf13718cd58c1b044` |
| `TargetedScreen_LFC_byGene.tsv.gz` | 49291630 | 77,563,731 | `f75a54e972e14dad2465070e2108f58a` |
| `TargetedScreen_LFC_byGuide.tsv.gz` | 49291639 | 231,621,570 | `925fb7a41c5152993007b625300e13c2` |
| `TargetedScreen_LFC_byGene-perLine.tsv.gz` | 49291675 | 1,919,661,029 | `27ae8d0109b9c42f1f972c45987db77e` |
| `TargetedScreen_LFC_byGuide-perLine.tsv.gz` | 49291687 | 3,592,732,872 | `30ebe001680aeab38ae2483c0ed0e95a` |
| `00_raw_data_id_mapping.csv` | 50865516 | 122,631 | `f698f9da456bbe480992fb2b7f509947` |
| `scripts.zip` | 58551544 | 3,303,491 | `51e9b796e1fba046a7305f17fc7ba592` |
| `Peer Review.pdf` | 58551748 | 2,453,785 | `418b84a003550805b7bbc596f03aa202` |

Stable download form:
`https://ndownloader.figshare.com/files/<FIGSHARE_ID>`.

Count Data:

- article `27989294`, version 2
- DOI: `10.6084/m9.figshare.27989294.v2`
- license: **MIT**
- nine files, **7,925,882,368 bytes** total
- includes RNA UMI counts, guide UMI counts, and cell metadata for fitness,
  non-fitness, and targeted screens

Only the three metadata files were retained locally. The RNA and guide matrices
were unnecessary once the pooled LFC artifact was found.

### Raw reads, code, mirrors, and negative searches

- Raw sequencing: ENA/SRA study `ERP165335` / BioProject `PRJEB81502`.
  ENA exposes **564 runs**, totaling **15,852,618,838,414 FASTQ bytes** and
  **11,585,012,548,046 submitted CRAM bytes**. These are raw reads, not a
  practical effect matrix.
- Official code:
  `https://github.com/claudiafeng123/crispri_scrnaseq_hipsci`, MIT, pinned
  commit `16ee5c09b7e5fa6410305160eec73860b884bc73`.
- The exact source code bundle was downloaded from Figshare and hashes to
  SHA-256
  `e9670eb5927394134682d58c1c56cb8e14fc5f5fd4a7b3609597f1cb9f40c3e8`.
- Hugging Face mirror:
  `altoslabs/primeflow-vcc-datasets/feng24_preprocessed.h5ad.gz`,
  **10,252,894,481 bytes**, Xet SHA-256
  `40a6c426bad1e8353c5ee8e72d2622bb4c7728d52162866481f0669644048bac`.
  It is a large reprocessed cell-level mirror and adds no advantage for this
  audit. The compilation is CC BY 4.0; its retained Feng notice is MIT.
- A third-party GitHub importer,
  `Katsevich-Lab/import-feng-2025`, only reconstructs a guide-count matrix
  (6,824 guides × 322,746 cells); it has no gene-expression effect matrix.
- No GEO series was indexed for the title/accession, and an exact-title Zenodo
  search found no official copy. The authoritative archival routes are
  Figshare plus ENA/SRA.
- Official paper supplements include Table S1 (target selection, 877,590
  bytes) and Table S2 (assigned-cell counts, 309,560 bytes). The Figshare
  metadata files expose the same cell-level coverage needed here more directly.

## The actual combined effect artifact

The selected source is:

```text
https://ndownloader.figshare.com/files/49291654
GenomeWideScreen_LFC_byGene.tsv.gz
794,610,300 bytes
MD5 03cb8a3b4328de50d45746a9d784bbc5
```

Schema:

```text
Target
Expressed_Gene_Symbol
Expressed_Gene_Ens_ID
lfc
pval_adj
```

The downloaded byte stream is a rectangular table of:

- **6,673 targets**
- **6,472 stored features**
- **43,187,656 data rows**

The paper says 6,471 genes and 43,180,983 effects. The apparent discrepancy is
fully explained by one extra non-biological feature, `PURO`, present once for
every target:

```text
6,673 × (6,471 + PURO) = 43,187,656
```

The extractor verified the exact source byte size and MD5, the complete row
count, target count, feature count, rectangular feature order, finite selected
values, p-value range, and absence of duplicate selected rows.

The published `lfc` is the authors' covariate-adjusted mean perturbation effect.
Their method first fits control-cell expression against cell line, inlet,
mitochondrial percentage, cell-cycle scores, and total UMI count, then averages
the residual perturbation effect over cells assigned to each target. Values are
preserved as published; no guessed log-base conversion or p-value threshold was
applied.

## Exact VCC coverage

### Effect table

| Quantity | Value |
|---|---:|
| VCC panel targets | 300 |
| Targets with published pooled LFCs | **182 (60.7%)** |
| VCC genes | 18,533 |
| VCC genes measured by Feng | **6,151 (33.2%)** |
| Retained finite target-gene effects | **1,119,482** |
| Effects per covered target | **6,151** |
| Covered targets with their own gene measured | 140 |
| Self-LFC negative | **140/140 (100%)** |
| Self-LFC significant at adjusted p < 0.1 | **80/140 (57.1%)** |
| Median self-LFC | **-0.2773** |

### Screen provenance and depth

| Screen | Cell lines | All targeting cells | Targets in metadata | VCC targets in metadata | VCC targets in released LFC table | Median VCC cells/target |
|---|---:|---:|---:|---:|---:|---:|
| Genome-wide fitness | 24 | 104,094 | 2,254 | 31 | **31** | 50 |
| Genome-wide non-fitness | 34 | 115,112 | 4,955 | 166 | **151** | 22 |
| Combined genome-wide | 34 | **219,206** | 7,207 | 197 | **182** | **25** |
| Targeted | 19 | 626,781 | 444 | **5** | separate table | **1,544** |

The five targeted-screen overlaps are `EZH2`, `GIGYF2`, `NT5DC1`, `STT3A`,
and `ZNF526`; all five are already present in the pooled genome-wide effect
table. The targeted per-line artifact is therefore potentially useful only as
a high-depth five-gene calibration check, not as a broad donor.

The 182 usable genome-wide panel targets consist of **31 fitness-screen** and
**151 non-fitness-screen** targets. The earlier 166-target count described
metadata membership only and omitted the fitness screen; 15 of those 166
non-fitness targets also fail the authors' minimum-10-cell effect-release
threshold.

Among the 182 released panel targets:

- depth range: **10–102 cells**
- median: **25**
- at least 25 cells: **94**
- at least 50 cells: **22**
- at least 100 cells: **2**
- at least 150 cells: **0**
- median distinct cell lines represented per target: **13**

This is shallow by our donor-bank standard, but the released values are already
pooled, covariate-adjusted model estimates rather than naïve pseudobulks.

## Compact local artifact

Location:
`data/external/feng-ipsc/genomewide_vcc_panel/`

| File | Meaning | Bytes |
|---|---|---:|
| `lfc.npy` | published LFCs, 182 × 18,533; NaN where unavailable | 13,492,152 |
| `centered_lfc.npy` | LFC minus each gene's mean over all 6,673 source targets | 13,492,152 |
| `pval_adj.npy` | published adjusted p-values | 13,492,152 |
| `observed.npy` | target-gene availability mask | 3,373,134 |
| `source_gene_mean.npy` | global source centering vector on VCC axis | 74,260 |
| `measured_gene_mask.npy` | 6,151 measured VCC genes | 18,661 |
| `source_gene_axis.tsv` | all 6,472 source features and Ensembl mappings | 208,966 |
| `targets.txt` | 182 targets in VCC panel order | 1,131 |
| `genes.txt` | canonical 18,533-gene VCC order | 119,285 |
| `meta.json` | provenance, checks, coverage, and artifact hashes | 8,124 |

`extract_feng_effects.py` streamed the full official gzip twice during
development but retained neither copy. The final extraction used **0.138 GiB
peak RSS**. Metadata auditing used **0.052 GiB peak RSS**.

## Raw-space LODO screen

### Protocol

- Main donors: K562-GWPS, Orion HCT116, Orion HEK293T.
- Gallery: **259** VCC targets shared by all three main donors.
- Features: **7,312** shared non-panel genes.
- Feng availability within that gallery: **165 targets × 5,137 features**.
- Each main donor was centered over all of its perturbations.
- Feng was centered over all 6,673 source targets.
- Missing Feng targets and genes preserved the baseline exactly.
- Feng was norm-matched on its measured coordinates before blending.
- Fractions: `0.02, 0.05, 0.10, 0.20, 0.35`.
- Confidence policies:
  - all 165 covered targets;
  - source depth ≥25 cells (84 targets);
  - source depth ≥50 cells (19 targets);
  - significant on-target knockdown, adjusted p <0.1 (76 targets).
- Candidate selection used HCT116 and HEK293T only. K562 was frozen until one
  candidate had been selected.
- No whitening, count emission, `cell-eval2`, or submission.

### Complete development grid

Values are full-gallery raw PDS deltas.

| Policy | Fraction | HCT116 | HEK293T | Development mean | Worst fold |
|---|---:|---:|---:|---:|---:|
| all | 0.02 | -0.000149651 | +0.000149651 | 0.000000000 | -0.000149651 |
| all | 0.05 | +0.000044895 | -0.000119721 | -0.000037413 | -0.000119721 |
| all | 0.10 | -0.000419024 | -0.000329233 | -0.000374128 | -0.000419024 |
| all | 0.20 | -0.001152315 | -0.000613570 | -0.000882943 | -0.001152315 |
| all | 0.35 | -0.003307294 | -0.004160306 | -0.003733800 | -0.004160306 |
| depth ≥25 | 0.02 | -0.000074826 | +0.000104756 | +0.000014965 | -0.000074826 |
| depth ≥25 | 0.05 | -0.000059861 | +0.000044895 | -0.000007483 | -0.000059861 |
| depth ≥25 | 0.10 | -0.000209512 | +0.000134686 | -0.000037413 | -0.000209512 |
| depth ≥25 | 0.20 | -0.000778187 | -0.000553710 | -0.000665948 | -0.000778187 |
| depth ≥25 | 0.35 | -0.001601269 | -0.002095118 | -0.001848194 | -0.002095118 |
| depth ≥50 | 0.02 | -0.000059861 | -0.000029930 | -0.000044895 | -0.000059861 |
| depth ≥50 | 0.05 | 0.000000000 | -0.000134686 | -0.000067343 | -0.000134686 |
| depth ≥50 | 0.10 | +0.000014965 | -0.000239442 | -0.000112238 | -0.000239442 |
| depth ≥50 | 0.20 | 0.000000000 | -0.000374128 | -0.000187064 | -0.000374128 |
| depth ≥50 | 0.35 | -0.000179582 | -0.000882943 | -0.000531262 | -0.000882943 |
| self p <0.1 | **0.02** | **+0.000014965** | **+0.000074826** | **+0.000044895** | **+0.000014965** |
| self p <0.1 | 0.05 | +0.000029930 | -0.000029930 | 0.000000000 | -0.000029930 |
| self p <0.1 | 0.10 | -0.000179582 | -0.000014965 | -0.000097273 | -0.000179582 |
| self p <0.1 | 0.20 | -0.000733291 | -0.000583640 | -0.000658466 | -0.000733291 |
| self p <0.1 | 0.35 | -0.001945467 | -0.003112747 | -0.002529107 | -0.003112747 |

### Predeclared gate

The evaluator declared a pass when:

```text
development worst-fold delta >= -0.002
and frozen K562 delta >= 0
and three-fold macro delta > 0
```

The selected `self_sig_0_1`, fraction `0.02` candidate met all three:

```text
development worst fold  +0.000014965
frozen K562              +0.000119721
three-fold macro         +0.000069837
```

So the machine-readable field is correctly `passes_screen: true`.

This was deliberately a permissive anti-catastrophe gate, not a practical
promotion threshold. Here one full-gallery rank step is approximately
`1 / (259 × 258) = 0.000014965`; the development worst-fold gain is exactly
one such step. The selected macro result is too small to distinguish from
rank granularity and validation noise. It must not be interpreted as evidence
that Feng improves the production model.

## Interpretation

The most informative pattern is monotonic harm as Feng weight increases:
fractions of 0.10–0.35 are negative in almost every fold and policy. Source
depth alone does not recover a useful signal. Restricting to knockdowns with
significant self-repression changes the sign at 2%, but not the effect size.

This supports two conclusions:

1. The pooled LFC table is substantially better than trying to reconstruct
   effects from 22–25 cells/target ourselves; the authors already performed
   the appropriate covariate adjustment and pooling.
2. Even that best available representation contributes essentially no
   transferable somatic-cell retrieval signal beyond the current donors.
   Pluripotent-context mismatch and shallow per-target depth remain binding.

The prior blanket statement that “Feng is not usable because no effect matrix
exists” is wrong. The corrected statement is:

> A verified pooled effect matrix exists and covers 182 VCC targets, but its
> incremental LODO signal is practically zero under a conservative direct
> blend. Retain it; do not promote it.

## Reproduction

```bash
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 0.95 --min-free-percent 25 --poll-seconds 5 -- \
  uv run python experiments/059-feng-ipsc/extract_feng_effects.py

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 0.95 --min-free-percent 25 --poll-seconds 1 -- \
  python3 experiments/059-feng-ipsc/audit_metadata.py

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 0.95 --min-free-percent 25 --poll-seconds 2 -- \
  uv run python experiments/059-feng-ipsc/evaluate_lodo.py

python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 0.95 --min-free-percent 25 --poll-seconds 0.5 -- \
  uv run pytest -q \
    experiments/059-feng-ipsc/test_extract_feng_effects.py \
    experiments/059-feng-ipsc/test_evaluate_lodo.py
```

Observed peaks:

- extractor: **0.138 GiB**
- metadata audit: **0.052 GiB**
- LODO evaluator: **0.23 GiB**
- tests: **0.01 GiB**, **4 passed**

Primary machine-readable outputs:

- `data/external/feng-ipsc/genomewide_vcc_panel/meta.json`
- `data/external/feng-ipsc/metadata_coverage.json`
- `experiments/059-feng-ipsc/lodo_results.json`
- `experiments/059-feng-ipsc/lodo_results.csv`
