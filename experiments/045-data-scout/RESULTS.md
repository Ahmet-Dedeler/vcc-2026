# 045 — Data scout: CD4 T-cell, Altos VCC bundle, and lineage-matched screens

Checked **2026-08-22** using repository APIs, GEO/ArrayExpress metadata, public
S3 listings, and small files only. No artifact larger than 250 MB was
downloaded. No scorer, eigendecomposition, or donor-table build was run.

## Bottom line

The **primary CD4+ T-cell CRISPRi atlas is the clear acquisition priority**.
Its released differential-expression object covers **291/300 VCC targets**,
including **27 targets absent from K562-GWPS**. The union of CD4, K562-GWPS,
and the Feng iPSC screen covers **296/300** targets; only `ABCD1`, `EPHB2`,
`KIF21B`, and `NICN1` remain uncovered in the metadata checked here.

For this project, do **not** begin by downloading the 12 cell-level CD4 files
(1.843 TB total as stored in the source bucket, 197.52 GB in Altos' gzip
repackaging). Instead, copy the single **16.786 GB**
`GWCD4i.DE_stats.h5ad` to an ephemeral cloud disk and extract only:

- the 291 panel targets across `Rest`, `Stim8hr`, and `Stim48hr`;
- `log_fc`, `zscore`, `lfcSE`, and adjusted p-values;
- the 10,282-gene axis and relevant `.obs` quality columns.

That panel slice is only about 36 MB per dense float32 layer before
compression. It supports both direct T-lineage donor transfer and a
T-cell-specific whitening covariance without ingesting 22 million cells.

The currently available broad Jurkat screen (Nadig) has **0/300 target
overlap** because it targets common-essential genes. The two older Jurkat
CROP/scifi screens also have **0/300 overlap** and only 20–32 target genes.
They remain useful only as small, lineage-matched covariance or modality
sanity checks.

No broad, transcriptome-readout CRISPRi atlas for **HeLa, CAL27, or CAL33**
was found in the authoritative sources searched. The available HeLa and
head-and-neck records are only a few bulk-RNA perturbations or pooled
fitness screens; they are not serious whitening datasets.

## 1. Primary CD4+ T-cell CRISPRi atlas

### Authoritative sources

- Dataset landing page:
  https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Public S3 prefix:
  `s3://genome-scale-tcell-perturb-seq/marson2025_data/`
- HTTPS prefix:
  https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/
- GEO raw-data record:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342
- Analysis repository and detailed schema:
  https://github.com/emdann/GWT_perturbseq_analysis_2025
- Manuscript:
  https://doi.org/10.64898/2025.12.23.696273

The S3 metadata identifies the assay as 10x Gene Expression Flex, human
primary blood CD4+ T cells, four donors, and three states (`Rest`,
`Stim8hr`, `Stim48hr`). The dataset is MIT-licensed. GEO reports 22 million
cells and became public on 2026-02-06.

### Best artifacts

| Artifact | Exact bytes | Purpose |
|---|---:|---|
| `GWCD4i.DE_stats.h5ad` | 16,786,240,107 | Best first acquisition: pooled target-level effects |
| `GWCD4i.DE_stats.by_donors.h5mu` | 16,866,278,447 | Donor-pair effects for confidence and nested validation |
| `GWCD4i.DE_stats.by_guide.h5mu` | 29,424,424,894 | Guide-level robustness/off-target analysis |
| `GWCD4i.pseudobulk_merged.h5ad` | 44,566,657,140 | Raw summed UMI pseudobulks; needed only for custom contrasts |

Direct URL pattern:

`https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/<ARTIFACT>`

The full cell-level release contains 12 files named
`D{1-4}_{Rest,Stim8hr,Stim48hr}.assigned_guide.h5ad`, ranging from
118,609,904,492 to 172,796,432,870 bytes. Their total is
**1,843,493,571,254 bytes (1.843 TB / 1.717 TiB)**. They are unnecessary for
the first whitening/donor-transfer experiment.

### Schema

`GWCD4i.DE_stats.h5ad`:

- shape: **33,983 perturbation-condition rows × 10,282 measured genes**;
- `.obs`: target symbol/Ensembl ID, condition, cell count, up/down/total DE
  counts, on-target effect/significance, off-target flags, guide count,
  downstream count, guide correlations, and donor-pair correlations;
- `.var`: `gene_ids`, `gene_name`;
- `.layers`: `log_fc`, `p_value`, `adj_p_value`, `baseMean`, `lfcSE`,
  `zscore`;
- `.varm`: condition-specific measured-gene statistics.

`GWCD4i.pseudobulk_merged.h5ad`:

- one row per guide × donor × condition pseudobulk;
- 18,129 genes;
- `.X` is summed UMI counts;
- `.obs` includes donor, condition, guide, target, cell count, total counts,
  and DE eligibility flags.

The cell-level files contain sparse UMI counts in `.X`, with donor/condition,
guide assignment, target symbol/ID, QC, and guide quality fields in `.obs`.

### Measured VCC overlap

Using the current `data/challenge/pert_counts.csv`:

| Condition | CD4 targets | VCC targets present | On-target significant among panel | Median cells / panel target |
|---|---:|---:|---:|---:|
| Rest | 11,287 | 283 | 246 | 512 |
| Stim8hr | 11,415 | 287 | 255 | 490 |
| Stim48hr | 11,281 | 284 | 250 | 495.5 |
| Any condition | 11,526 unique | **291** | **267 in at least one condition** | — |

Of the 291 overlapping panel genes, 276 have at least one reported downstream
effect. The nine absent targets are `ABCD1`, `C5orf22`, `EEF1A2`, `EPHB2`,
`FZD2`, `KIF21B`, `LENG1`, `NICN1`, and `REV1`.

K562-GWPS overlaps 269 panel genes. CD4 and K562 share 264 panel targets;
CD4 supplies 27 K562-missing targets, while K562 supplies five CD4-missing
targets. Adding Feng raises union coverage to 296/300.

### Small files worth acquiring immediately

All are below the requested 250 MB ceiling and live in:
https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/master/metadata/suppl_tables

| File | Bytes | Utility |
|---|---:|---|
| `DE_stats.suppl_table.csv` | 7,276,533 | Target/condition coverage and QC; no gene-level effect vectors |
| `guide_kd_efficiency.suppl_table.csv` | 13,319,852 | Guide effectiveness and high-confidence-no-effect flags |
| `sgrna_library_metadata.suppl_table.csv` | 9,935,249 | Target/guide mapping and off-target annotations |
| `K562_comparison.suppl_table.csv` | 969,086 | Cross-line transferability by target |
| `clustering_results_and_annotations.csv` | 121,625 | Maps 1,692 genes into 112 response-program clusters |
| `clustering_downstream_genes.csv.gz` | 19,694,461 | Cluster/condition downstream-gene ranks and sign coherence |

The cluster summaries alone cover **56/300 VCC target genes** as cluster
members and may support a cheap low-rank T-cell prior before the full DE
object is extracted.

## 2. Altos Labs PRiMeFlow VCC-preprocessed bundle

### Authoritative sources and license

- Dataset:
  https://huggingface.co/datasets/altoslabs/primeflow-vcc-datasets
- Repository:
  https://github.com/altoslabs/primeflow
- Paper:
  https://arxiv.org/abs/2604.13986

Dataset revision checked: `d773a3476d15e5d32d385b13cebf87e62580097e`,
last modified 2026-08-11.

The **dataset compilation** is CC BY 4.0. The card explicitly says this does
not replace underlying licenses. CD4 and Feng retain MIT notices. The
PRiMeFlow **software** has a separate academic/non-commercial license; do not
confuse that with the dataset license.

The current non-legacy `.h5ad.gz` payload is **244,368,801,634 bytes
(244.37 GB / 227.59 GiB)**. CD4 alone is **197,523,199,492 bytes
(197.52 GB / 183.96 GiB)**.

### Relevant exact files

Each file is downloadable as:

`https://huggingface.co/datasets/altoslabs/primeflow-vcc-datasets/resolve/main/<FILENAME>?download=true`

| File | Exact bytes | SHA-256 / Xet hash |
|---|---:|---|
| `nadig24_jurkat_preprocessed.h5ad.gz` | 1,935,671,130 | `b154ce033e21e3accd1a6de6be165dc8840191bc10f84ae14172cec0149f2166` |
| `jiang24_preprocessed.h5ad.gz` | 10,043,729,912 | `0ae2b8e3e971448b0b656a3164a9d5d110659e4a35831d5cb6697c038710e855` |
| `feng24_preprocessed.h5ad.gz` | 10,252,894,481 | `c1a9c8da415727043a12a4f2ccee763e6523e851d8ffc1b0ed912638e3a0d652` |
| `replogle22_k562_preprocessed.h5ad.gz` | 13,678,926,390 | `573b646b816d0501a2e91fc1a28aa62095f31de06073cc0b187d866f331c4df8` |
| `replogle22_rpe1_preprocessed.h5ad.gz` | 1,888,128,718 | `14a12a58f579fcd8e48795f863c169330bf6994fcb99df0f485d34c974a480b6` |
| `replogle20_preprocessed.h5ad.gz` | 252,963,483 | `45e5802f304e417e83924e1949b68db85a25a40b3dc800d20a0bdc086a1cd123` |
| `vcc_train_preprocessed.h5ad.gz` | 4,291,253,394 | `936e9a9a9d54556368e34ab4043155b6432263b97bd16f300cb9680ddd04e2d2` |
| `vcc_filtered_genes_18001.csv` | 115,478 | — |

The 12 CD4 gzip files are, in donor/condition order:

| File | Exact bytes |
|---|---:|
| `cd4t_perturbseq_D1_Rest_preprocessed.h5ad.gz` | 13,947,693,582 |
| `cd4t_perturbseq_D1_Stim8hr_preprocessed.h5ad.gz` | 15,249,649,141 |
| `cd4t_perturbseq_D1_Stim48hr_preprocessed.h5ad.gz` | 16,698,313,911 |
| `cd4t_perturbseq_D2_Rest_preprocessed.h5ad.gz` | 16,122,694,401 |
| `cd4t_perturbseq_D2_Stim8hr_preprocessed.h5ad.gz` | 20,424,701,893 |
| `cd4t_perturbseq_D2_Stim48hr_preprocessed.h5ad.gz` | 17,152,298,762 |
| `cd4t_perturbseq_D3_Rest_preprocessed.h5ad.gz` | 16,162,768,261 |
| `cd4t_perturbseq_D3_Stim8hr_preprocessed.h5ad.gz` | 16,991,236,780 |
| `cd4t_perturbseq_D3_Stim48hr_preprocessed.h5ad.gz` | 17,711,508,406 |
| `cd4t_perturbseq_D4_Rest_preprocessed.h5ad.gz` | 13,198,211,441 |
| `cd4t_perturbseq_D4_Stim8hr_preprocessed.h5ad.gz` | 15,161,706,400 |
| `cd4t_perturbseq_D4_Stim48hr_preprocessed.h5ad.gz` | 18,702,416,514 |

### Standardized schema

PRiMeFlow consumes these as AnnData with:

- `.obs['condition']`: perturbed gene, with `control` as the control label;
- covariates: `cell_type`, `treatment`, `dataset`, `vector`,
  `scrna_method`;
- `.obs_names`: cell IDs used by the split CSVs;
- gene-symbol feature axis, optionally filtered by
  `vcc_filtered_genes_18001.csv`;
- normalized/log1p expression in `.X`; source counts may remain in
  `.layers['counts']` where retained by the source-specific curation.

The exact preprocessing transform is not fully documented in the dataset card,
so this should be verified from one file before mixing it numerically with our
raw-count-derived donor tables.

### Utility notes

- **CD4:** highest value; direct coverage and likely lineage match for context A.
- **Jiang:** 1,628,476 cells × 15,476 genes, 218 perturbations plus control,
  six cell lines and five signaling treatments. No HeLa/CAL27/CAL33, but strong
  multi-line covariance and context-transfer training value.
- **Feng:** 34 iPSC lines. Local metadata shows 166 panel overlaps in the
  genome-wide non-fitness screen and five in the targeted screen; useful for
  cross-genotype stability and missing-gene coverage, less lineage-matched.
- **Nadig Jurkat:** lineage matched, but 0/300 direct target overlap. The
  already-built local donor table has 2,270 retained perturbations and 8,881
  genes, so downloading Altos' duplicate is unnecessary unless its
  normalization is specifically required.

## 3. Other Jurkat resources

### Nadig et al. — strongest existing Jurkat atlas

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667
- Direct file:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad
- Size: approximately 8.7 GB on GEO.
- Schema: raw single-cell AnnData; 10x 3' v3; dual-sgRNA vectors; target labels
  in `obs['gene']`, symbols in `var['gene_name']`.
- Design: Jurkat and HepG2 screens run in parallel, primarily DepMap
  common-essential genes, assayed seven days after transduction.
- License: GEO does not attach a machine-readable per-dataset license. Check the
  associated publication terms for redistribution.
- Current local state: already downloaded and reduced to
  `data/donors/nadig_jurkat/`; 2,270 retained perturbations, 8,881 genes,
  **0/300 VCC target overlap**.

Use it as a Jurkat covariance/basal-state control, not a direct target donor.

### Datlinger 2017 Jurkat CROP-seq

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92872
- Raw count matrix:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92872/suppl/GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz
  (about 17 MB)
- Guide count matrix:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92872/suppl/GSE92872_CROP-seq_Jurkat_TCR.count_matrix.csv.gz
  (about 1.5 MB)
- Ready curated h5ad:
  https://zenodo.org/api/records/13350497/files/DatlingerBock2017.h5ad/content
  (**39,085,166 bytes**, MD5
  `c797d8e0320afcc0cc2a94e94805b142`)
- Curated license: CC BY 4.0.
- Verified schema: 5,905 cells × 36,722 genes, sparse `.X`;
  `.obs` includes `target`, guide-level `perturbation`, stimulated/unstimulated
  state, replicate, and Jurkat annotations.
- Coverage: 32 targets, **0/300 VCC overlap**.

### Datlinger 2021 Jurkat scifi-RNA-seq CRISPR screen

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168620
- Relevant source file:
  `GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells.h5ad.gz`,
  **292,102,388 bytes** (over the present 250 MB download limit).
- Ready curated h5ad:
  https://zenodo.org/api/records/13350497/files/DatlingerBock2021.h5ad/content
  (**33,600,780 bytes**, MD5
  `5cfbcb4770e9c859a07b2cc44fd12066`)
- Curated license: CC BY 4.0.
- Verified schema: 39,194 cells × 25,904 genes, sparse `.X`;
  `.obs` includes guide-level `perturbation`, stimulated/unstimulated state,
  and Jurkat annotations.
- Coverage: 20 genes represented by 40 targeting guides, **0/300 VCC overlap**.

These two small h5ads are safe diagnostic datasets, but adding 20–32
TCR-pathway directions to a whitening covariance is unlikely to move a
12,000-plus-vector basis materially.

## 4. HeLa resources

No broad HeLa single-cell CRISPRi/Perturb-seq expression atlas was found in
GEO, ArrayExpress/BioStudies, scPerturb, Hugging Face, GitHub, or the literature
indexes searched.

The best exact match is a small bulk-RNA study:

- Record:
  https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-5308
- Processed counts:
  https://www.ebi.ac.uk/biostudies/files/E-MTAB-5308/lncRNA_counts.tsv
  (**103,418,822 bytes**)
- Additional counts:
  https://www.ebi.ac.uk/biostudies/files/E-MTAB-5308/additional_counts.tsv
  (**8,103,171 bytes**)
- SDRF:
  https://www.ebi.ac.uk/biostudies/files/E-MTAB-5308/E-MTAB-5308.sdrf.txt
  (**1,320,372 bytes**)
- Study: 126 HeLa samples comparing CRISPRi, RNAi, and LNA.
- CRISPRi targets: `LOC100289019`, `MALAT1`, `H19`, and `CKAP5`, plus
  negative-guide and dCas9-KRAB clone controls.
- VCC panel overlap: **0/300**.
- Schema: tab-delimited bulk gene-count matrices plus SDRF sample metadata;
  hg38/Ensembl 90 featureCounts.
- License: the paper is CC BY 4.0; the BioStudies record does not expose a
  separate machine-readable dataset license.

This is too small and clone-confounded for a whitening basis. Its only plausible
use is estimating HeLa-specific background/off-target directions to regress
out, and even that should be tested very cautiously.

`GSE207330` adds bulk RNA-seq after CRISPRi of `EIF3D` in HeLa and Jurkat, but
provides processed matrices only for K562; the HeLa/Jurkat data are raw reads.
It has no VCC target overlap and is not worth acquiring.

## 5. Head-and-neck / squamous resources

No broad transcriptomic CRISPRi atlas was found for CAL27, CAL33, or another
close HNSCC line.

Available low-value records:

| Record | Readout | Artifact | Size | VCC target overlap |
|---|---|---|---:|---:|
| `GSE316531` | CAL27 kinome CRISPR-Cas9 survival screen under radiation | `GSE316531_Cal27_CRISPR_Screening_Results_JAK1_Paper.xlsx` | 414 KB | Not expression |
| `GSE153354` | CAL27 TCF19 knockout bulk RNA-seq, 3 WT + 3 KO | `GSE153354_RAW.tar` | 1,792,000 B | 0 |
| `GSE226355` | SCC15/SCC25 YBX1 knockout bulk RNA-seq | two normalized count tables | 1.0 MB + 762 KB | 0 |
| `GSE316529` | Detroit562 JAK1 knockout bulk RNA-seq | one XLSX | 3.7 MB | 0 |

Direct GEO URL pattern:

`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ACCESSION>`

These can inform lineage-specific pathway interpretation or viability gating,
but they do not supply enough perturbation directions for covariance
estimation and should not delay the CD4 extraction.

## 6. Prioritized acquisition and experiment plan

1. **Cloud-extract `GWCD4i.DE_stats.h5ad` first.**
   Use the public S3 object directly on an ephemeral AWS instance. Export a
   compact artifact containing panel rows, gene symbols, `log_fc`, `zscore`,
   `lfcSE`, `adj_p_value`, and `.obs` QC. Do not copy the full object back to
   the Mac.

2. **Build two distinct CD4 candidates.**
   - `Rest` direct donor for context A only.
   - All-condition CD4 covariance as an auxiliary whitener, keeping the
     production K562/Orion direction source fixed.
   Tune any blend only with donor/condition holdouts and existing H1/K562/RPE1
   validation—not leaderboard probing.

3. **If step 1 helps, extract the panel slice from
   `GWCD4i.DE_stats.by_donors.h5mu`.**
   Use donor-pair agreement to weight targets and to test whether an apparent
   gain survives donor holdout.

4. **Use the small CD4 cluster files immediately as a cheap prototype.**
   They can test a cluster-program residual or DE-ranking prior without large
   storage. A positive result justifies the full DE extraction; a negative
   result does not rule it out because the cluster summaries discard
   target-specific magnitude.

5. **Do not reacquire Nadig Jurkat.**
   The raw source and compact donor table already exist locally. It is
   covariance-only and has zero panel overlap.

6. **Do not bulk-download the Altos bundle.**
   Pull individual artifacts only. For current goals, the Altos CD4 cell-level
   files are inferior to the compact source DE object; Jiang and Feng are the
   next worthwhile multi-context acquisitions after CD4.

7. **Defer HeLa/HNSCC bulk and fitness datasets.**
   They contain too few transcriptomic perturbations and no panel-target
   matches. Basal CCLE similarity is more useful than these response datasets
   unless a new broad Perturb-seq release appears.

## Safety / process note

During this audit, an orphaned `experiments/042-auxiliary-whitening` process
was discovered at roughly 5 GB RSS and terminated. A subsequently appearing
`experiments/044-csls-distill` process was also terminated. No local VCC
Python/scorer process remained after the cleanup; the only repository-related
background operation observed was an existing AWS upload.
