# 482c — Literature and data hunt (2026-09-26)

**Decision.** No new public genome-scale CRISPRi transcriptome is available for a cancer line outside the banks already on disk. The one new ≥1,000-target cancer-line screen that could be opened, CZ Biohub's A549 Cas9 knockout CROP-seq (GSE345058), knocks the target's own mRNA down only part way (median residual 0.58 of control CPM). Xaira’s promised genome-scale CRISPRi for Jurkat and HepG2 is still an empty Hugging Face repo. HeLa and CAL-33 still have no public single-cell gene-knockdown atlas. The 2026 challenge sources do not name the final three lines.

Same-line mode stays limited to the lines already in the registry. A549 gains a deep knockout atlas beside the existing CRISPRi banks; it is a weak CRISPRi donor because the cis mRNA effect is partial.

Facts below are **VERIFIED** (this session opened the file or the primary page) or **REPORTED** (secondary, or a primary page that returned an error).

---

## Ranked shortlist

### 1. A549, 1,000-gene Cas9 knockout CROP-seq — GSE345058

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345058 · preprint https://doi.org/10.64898/2026.06.01.728087 · FTP `geo/series/GSE345nnn/GSE345058/suppl/` |
| Line | A549 (lung). CZ Biohub, Leonetti lab. Public 26 Aug 2026. |
| Modality | CRISPRko (dox-inducible Cas9), CROP-seq vector. Guides: 4,000 targeting (Brunello 3,478, TKOv3 609, GeCKO 124) + 211 nontargeting. **VERIFIED** from `sgRNA.csv.gz` and the GEO overall design. |
| Targets / cells | 1,000 genes + Nontargeting. 606,075 cells (560,755 perturbed, 45,320 NTC). Cells/target: median 565, p10 349, min 73 (RPL9), max among targets 1,196 (MXRA8). Mean library 28,723 (NTC). **VERIFIED** from `cell_annotations.csv.gz`. |
| Chemistry | 10x GEM-X, four batches × 8 lanes. **VERIFIED** from GEO overall design. Gene axis 29,892 symbols, GRCh38. |
| Licence | Not stated on the GEO record. bioRxiv page returned HTTP 429, so the preprint licence was not read. |
| Size | FTP listing: h5ad 23G, normalized mtx 12G, raw mtx 10G, annotations 38M + 370K + 67K. Opened raw mtx is 10,956,707,618 bytes. |
| Why it could help | A549 is already a same-line line (Jiang + Myllia CRISPRi). This is a different modality with ~565 cells/target, deep enough to estimate a knockout response. Cis knockdown is too weak to treat it as a CRISPRi donor. |

**Knockdown (VERIFIED).** Streamed `GSE345058_raw_counts.mtx.gz` with `kd_scan.c`. Header `%%MatrixMarket matrix coordinate integer general`, dimensions 606,075 × 29,892, 3,598,670,133 entries; the scanner counted the same 3,598,670,133 entries. Residual = pseudobulk CPM in cells assigned to the target / pseudobulk CPM in Nontargeting cells. CPM = gene-count sum / sum of per-cell `total_counts` × 10^6. Denominator: exact symbol match to the gene axis, and NTC CPM ≥ 5.

| Floor | n | fraction residual < 1 | fraction ≤ 0.5 | median residual |
|---|---:|---:|---:|---:|
| NTC CPM ≥ 5 | 832 | 0.974 | 0.343 | 0.575 |
| NTC CPM ≥ 1 | 873 | 0.964 | 0.329 | 0.582 |
| NTC CPM ≥ 10 | 780 | 0.980 | 0.351 | 0.567 |

951 of 1,000 target symbols sit on the gene axis; 49 do not (HGNC renames such as QARS→QARS1, or symbols absent from this 29,892-gene list). Two on-axis genes have NTC CPM = 0 (TRIM43, KRTAP4-2). Examples at CPM ≥ 5: RPL9 residual 0.44 (73 cells), KIF11 0.32 (86 cells), RPL26 0.43, KRAS 0.97. UBC, HSPA5 and DYNLL1 go up (residuals 4.5, 4.2, 3.0). Numbers are in `work/a549_kd.json`.

### 2. X-Atlas / Pisces — Jurkat and HepG2 genome-scale CRISPRi, not released

| | |
|---|---|
| URL | https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces · paper https://doi.org/10.64898/2026.03.18.712807 |
| What is public today | **VERIFIED** 2026-09-26 via the Hugging Face API: `usedStorage` 0, `lastModified` 2026-03-25T22:58:18Z, files `.gitattributes`, `LICENSE.md`, `README.md`, card licence `cc-by-nc-sa-4.0`. |
| What the paper describes | **REPORTED** (paper PDF was read earlier in this hunt; not re-downloaded this hour). Seven screens, 25.6M cells: HCT116 18,924 perts / 3.4M cells; HEK293T 18,312 / 4.5M; HepG2 9,735 / 2.6M; iPSC 10,095 / 4.2M; Jurkat resting 10,872 / 2.8M; Jurkat CD3/CD28 10,878 / 2.8M; iPSC multi-lineage Flex 12,175 / 5.1M. Paper median knockdowns are high (Jurkat resting ~79%, active ~71%, multi-diff ~96%, HCT116 ~70%). |
| Why it could help | Genome-scale CRISPRi in Jurkat would turn on same-line mode for validation context A, and HepG2 is already a registry line. Nothing is downloadable. |

### 3. A375 CRISPRi, PerturbFate — GSE291147

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291147 · code https://github.com/JunyueCaoLab/PerturbFate |
| Line | A375 melanoma, dCas9-BFP-KRAB, BRAF V600E. Public 6 Feb 2026. **VERIFIED** from the GEO series record and GSM8829735. |
| Modality | CRISPRi, combinatorial-indexing PerturbFate (not 10x). Overall design: 143 gRNA pairs targeting 143 genes plus no-target, DMSO or 1 µM vemurafenib. **VERIFIED.** |
| What was opened | `GSE291147_Dual_omics_RNA_gene_count_matrix.RDS` is a `dgCMatrix`, 60,706 × 100,590, 135,443,057 nonzeros, **no dimnames**. `GSE291147_Dual_omics_sgRNA_count_matrix.RDS` is a `dgTMatrix`, 155 × 78,510, **no dimnames**. RNA cell metadata has `conditions` (PLX vs the other arm) and no guide identity. The two matrices have different cell counts and cannot be joined. |
| Knockdown | Not computed. The deposited count objects do not carry gene or cell names. |
| Size opened | RNA RDS 302 MB, sgRNA RDS 1.9 MB, metadata 1.1 MB. Other supplements (ATAC h5ads, labeled/unlabeled RNA) were not downloaded. |
| Licence | Not read. |
| Why it could help | A375 is a cancer line that is not in the same-line registry. The 2025 Arc “Behind the Data” post names A375, with K562, as a line that already dominates public Perturb-seq. 143 genes is far below a genome-scale bank, and the screen is drug-split (vemurafenib vs DMSO). |

### 4. DLD-1 CRISPRi, ~5,000 genes, three MOIs — GSE337988

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337988 |
| Line | DLD-1. Public 10 Jul 2026. Sublibrary 2 at MOI ~0.3, ~3 and ~10, plus a 35-gene pilot. This is the public release of the low-MOI bank already held, plus the medium and high MOI of the same sublibrary. |
| What was opened | `GSE337988_sublib2_de_matrices_lfc_matrix_Low1.csv.gz` (149 MB). 6,451 perturbation columns, 4,955 gene symbols after stripping the last suffix (4,934 `P1`, 489 `P2`, …). **2,827 gene rows** — a DE subset, not the transcriptome. Guide file `NGS6475` feature barcodes, column `crispr`. |
| Knockdown | Among the 747 symbols that are both a row and a column, taking the `_P1` column when present: fraction LFC < 0 = 0.987, fraction LFC < −0.5 = 0.487, median LFC −0.488. **This is not a library-wide rate.** 4,208 target symbols are absent from the DE table. **VERIFIED.** |
| Why it could help | DLD-1 is already in the registry. Medium and high MOI are multiply infected cells of the same genes, not a new target set. |

### 5. NCI-H358 Cas12a Perturb-seq — GSE283574

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283574 |
| Line | NCI-H358 (a few samples are K562). Public 7 Apr 2026. GEO summary: optimized Cas12a vectors and a degron-enhanced Cas12a, Perturb-seq across “multiple cell types and target genes.” |
| What was opened | `GSE283574_RAW.tar` (299 MB) and the feature TSVs. GEX feature lists are 36,751 genes. The full-length feature files named both `EpiGenFactor` and `CellSurfaceReceptors` contain the **same 54 gene symbols** (epigenetic regulators: ACTL6A, ARID1A, DNMT1, …, plus non-targeting barcodes; 61 rows). Combinatorial dTag feature files have 79 symbols, 19 of them pairs. A distinct receptor gene list was not in the files that were extracted. |
| Knockdown | Not computed. The count matrices inside the tar were left packed. |
| Why it could help | New cancer line, but ~50 genes. Useful only if a final line is NCI-H358. |

### 6. hTERT RPE-1 viral microprotein ORF screen, 10x Flex — GSE344535

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344535 |
| What | Public 21 Aug 2026. GEO summary: dox-inducible ORF library of viral microproteins in hTERT RPE-1, two sublibraries by insert length, MOI < 30%, 10x Flex. 16 samples. **VERIFIED** from the GEO summary. RAW.tar is listed at 27 GB and was not downloaded. |
| Knockdown | Not applicable (overexpression of exogenous ORFs). |
| Why it could help | Flex chemistry in a line we already have (Replogle RPE1), but the perturbation is not gene CRISPRi. |

### 7. 769P kidney-cancer CROP-seq — GSE279287

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279287 |
| What | Public 10 Feb 2026. Four samples. GEO summary: CROP-seq in 769P cells to link renal-cell-carcinoma susceptibility variants to target genes, with a stated association of rs28684409 and RPL4. **VERIFIED** from the GEO summary. Counts not opened. |
| Why it could help | A cancer line (769P) that is not in the registry. The design is a variant-to-gene screen, not a gene-knockdown panel of useful size. |

### 8. HAP1 regulatory-element scMulti CRISPRi — GSE286226

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286226 |
| What | Public 10 Jan 2025. SuperSeries “Endogenous fine-mapping of functional regulatory elements in complex genetic loci in HAP1.” Child samples include single-cell multiome (ATAC + gene expression) labeled CRISPRi (GSM8720909–911). **VERIFIED** from GEO titles in the esummary catalog. Matrices not opened. |
| Why it could help | HAP1 is already a Jiang line. This is enhancer fine-mapping, not a new gene-target bank. |

### 9. K562 direct-capture TAP-seq, 97 genes — GSE303901

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303901 |
| What | Public 15 Aug 2025. “Unbiased survey of distal element–gene regulatory interactions with direct-capture targeted Perturb-seq.” Pilot sample GSM9137723: 97 genes in K562, CROP-HyPR, 10 guides. **VERIFIED** from the catalogued GEO title. Counts not opened. |
| Why it could help | K562 is already covered at genome scale. 97 enhancer-focused targets do not add a line. |

### 10. KOLF2.2J, 225 genes, embryoid bodies — GSE343369

| | |
|---|---|
| URL | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343369 |
| What | Public 12 Aug 2026. GEO summary: “phenotypes the knockout of 225 genes by single cell RNAseq during embryoid body formation,” guides against genes expressed in the pluripotent state or during EB formation. Four samples. **VERIFIED** from the GEO summary. h5ads were not downloaded (several are ~3 GB). |
| Why it could help | Same iPSC family as the KOLF atlas already held. 225 genes, differentiation time course, not a cancer line. |

---

## How Arc chose the 2026 lines and targets

**VERIFIED, Arc blog 20 Aug 2026** (https://arcinstitute.org/news/virtual-cell-challenge-2026), full text extracted:

- Zero-shot. No challenge training set. Models predict CRISPRi knockdown responses “in six cell lines they have never seen perturbed,” using the unperturbed state and a list of genes.
- “Perturb-seq across six cell lines from different tissues of origin,” CRISPRi, 10x Flex, Ultima UG100, “perturbations chosen to provide a strong set of perturbations and responses.”
- Three lines for validation and the live leaderboard; three held back for the final test.
- Each line ships non-targeting profiles plus the gene identifiers of the knockdown targets. Predictions are post-perturbation profiles “after CRISPRi knockdown followed by 10x Flex profiling.”
- Validation live 20 Aug 2026. Final test set 22 Oct 2026. Deadline 5 Nov 2026, 23:59 UTC. Prizes $100k / $50k / $25k. Final ranking is an aggregate of six metrics.
- The blog points teams at the Altos Labs flow-matching preprint before building.
- The page does not name a line, a tissue, or any of the 300 genes. String counts for Jurkat, HeLa, A549, K562, HCT116, HepG2, A375, MCF7 are zero.

**VERIFIED, Europe PMC** for DOI 10.1016/j.cell.2026.08.004, “Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts.” Abstract only (full text is not open on Europe PMC; `cell.com` HTML and PDF returned HTTP 403 from this VM): zero-shot prediction across multiple independent cellular contexts, gene knockdown responses in unseen cell lines, no line names. The commentary body was not re-read this session, so anything beyond that abstract is **REPORTED**.

**VERIFIED, Arc “Behind the Data,” 22 Aug 2025** (https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge). This post is about the **2025 H1** challenge, and it is the only public selection recipe:

- Modality: dual-guide CRISPRi (dCas9-KRAB). The post says Cas9 knockout “isn’t always readily reflected in its mRNA levels,” which matches the A549 measurement above.
- H1 was chosen because K562 and A375 “dominate existing Perturb-seq datasets,” and H1 is a shift relative to that pretraining mass.
- 300 genes came from an H1 pilot, binned by number of DE genes (negligible / subtle / strong), sampled for diversity with ContrastiveVI+, plus clean guide assignment, high basal target expression, and strong repression in Flex.
- Pilot comparison of 3′, 5′ and Flex on matched H1 perturbations; Flex won on UMI depth, genes detected, guide assignment, and perturbed-vs-control separation.
- The 2025 wrap-up (below) restates >50,000 UMIs per cell and ~1,000 cells per perturbation for that H1 dataset.

The 2026 blog does not say the 2026 panel was built the same way. The H1-vs-K562/A375 sentence is a statement about 2025 pretraining, not a list of 2026 test lines.

**Negative searches (VERIFIED):**

- GEO `HeLa AND (Perturb-seq OR CROP-seq OR CRISPRi)`, 13 series, all titles read. None is a single-cell gene-knockdown atlas. The only genome-scale hit is GSE85011 (public 1 Dec 2016): bulk CRISPRi growth phenotypes for 16,401 lncRNA loci in 7 lines, including transformed lines. GSE207330 is bulk RNA-seq of a few CRISPRi knockdowns in K562, HeLa and Jurkat. GSE228798 is a dCas9 roadblock study, not a Perturb-seq bank.
- GEO `(CAL-33 OR CAL33 OR Cal33) AND (Perturb-seq OR CROP-seq OR single-cell)`: count 0.

---

## 2025 winners and cross-line methods

**VERIFIED, Arc wrap-up 6 Dec 2025** (https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up):

- First place, team BM_xTVC (Yucheng Guo, Qirong Yang, BioMap), model `xTrimoSCPerturb`. Improved scFoundation, protein embeddings, a curated subset of the public Perturb-seq listed on the challenge site, pseudo-bulk, explicit training-set DEG frequency and mean expression, fused loss of PDS + DES + a small MAE term. Their quoted conclusion: purely AI approaches did not consistently beat statistical baselines on DES and MAE. Full methods “forthcoming.” No public GitHub was found.
- Second place, XLearning Lab (Xi Peng, Sichuan University), model `X`. Pseudo-bulk, fully connected net on control expression + ESM-2 + UMI indicators, residual/delta learning, MSE plus auxiliary losses, weighted toward PDS. Trained on public Perturb-seq including a small H1 set in PerturbAtlas. No official code repository was found.
- Third place, team Outlier (Chicago, Dartmouth, HKUST), model `TransPert`. Summary-level only: pseudo-bulk and Wilcoxon results, similarity-aware aggregation across reference lines, a call of which genes are significant in the target line, then global linear scaling fit by cross-validation on the challenge training data. arXiv 2511.16954 is cited by Arc for the PDS scale/metric point. Full TransPert methods “forthcoming.” No code repository was found.
- Generalist prize, Altos Labs, model string `go-with-the-flow`. Flow matching in gene-expression space, U-Net velocity. Pretrained on ~7M cells including Altos internal screens (none in H1), then fine-tuned on h1ESC. Highest average rank across seven metrics among the top 50. The 2026 blog tells teams to read this preprint.

**VERIFIED, PRiMeFlow.** https://arxiv.org/abs/2604.13986, “PRiMeFlow: Capturing Complex Expression Heterogeneity in Perturbation Response Modelling,” authors Zichao Yan, Yan Wu, and the same Altos team as the Generalist prize. Code https://github.com/altoslabs/PRiMeFlow (GitHub API 2026-09-26: 24 stars, pushed 2026-08-11, licence SPDX `NOASSERTION`).

**VERIFIED, TxPert.** https://arxiv.org/html/2505.14919 (v1, 20 May 2025), licence CC BY-NC-SA 4.0. Code https://github.com/valence-labs/TxPert (52 stars, pushed 2026-03-25). Data: Replogle 2022 essential CRISPRi in K562 and RPE1, Nadig 2025 CRISPRi in Jurkat and HepG2, Norman 2019 CRISPRa doubles in K562. On unseen single perturbations inside a known line, TxPert beats GEARS, scLAMBDA, and the non-learned general baseline (strongest of additive and mean); GEARS falls below that baseline. On leave-one-out cell lines (all perturbation examples of the target line held out, controls of that line still used in training), TxPert “exceeded the general baseline in all four held out cell lines” and beat scLAMBDA. The best reported graph mix is STRING + GO + Recursion PxMap and TxMap. A nearest-cell-line transfer comparison was not in the sections read.

**VERIFIED as existing code, claims not re-scored here.**

- Arc State, https://github.com/ArcInstitute/state (704 stars, pushed 2026-07-24). This repo already closed STATE on A/B/C (pds ≈ 0.50). Not re-run.
- CellFlow, https://github.com/theislab/cellflow (MIT, 160 stars, pushed 2026-09-21). Flow matching for drugs, genetic perturbations, cytokines, organoids. Whether it beats a linear model on a held-out cancer line was not re-measured.
- GEARS, https://github.com/snap-stanford/GEARS (MIT). TxPert’s own comparison puts it under the non-learned baseline on their single-perturbation task.
- CPA, https://github.com/facebookresearch/CPA (MIT, archived, last push 2023-09-07).
- Broad 16-line Cas9 screen code https://github.com/broadinstitute/perturb-seq-depmap-public (pushed 2026-09-19). The screen itself is already on disk (100 genes × 16 lines). Not a new download.

PerturBench (arXiv 2408.10609) is the standing caveat, **REPORTED** from the known preprint: on average-effect metrics many deep models do not beat linear or additive baselines. TxPert’s “beats the general baseline” claim is the one that was re-read in the paper HTML.

---

## Inventory

### Opened (counts or a DE table)

| Accession | Line | What was opened | Result |
|---|---|---|---|
| GSE345058 | A549 | raw mtx 3.599e9 entries, annotations, guides | 1,000-gene Cas9 KO, median residual 0.58. See shortlist 1. |
| GSE291147 | A375 | RNA and sgRNA RDS | 143 genes, CRISPRi, no dimnames, KD not computed. |
| GSE337988 | DLD-1 | Low-MOI LFC table | DE-subset only; 747 genes, median LFC −0.49. Already held. |
| GSE283574 | NCI-H358 | feature TSVs inside the RAW tar | ~54-gene Cas12a panel. KD not computed. |

### Primary text read, matrices not a new bank

| Source | Date read | Fact |
|---|---|---|
| HF `Xaira-Therapeutics/X-Atlas-Pisces` | 2026-09-26 | `usedStorage` 0 since 2026-03-25. |
| HF `arcinstitute/Perturb-Sapiens` | 2026-09-26 | 1.335e12 bytes, cytokine h5ads (`4-1BBL.h5ad`, …), CC-BY-4.0. Not CRISPR. |
| HF `Xaira-Therapeutics/X-Cell` | 2026-09-26 | API returned HTTP 401. Earlier in this hunt the dataset card was README + a PNG and no weights (**REPORTED** for the file list). |
| GSE344535 | GEO summary | RPE-1 ORF overexpression, 10x Flex. Not gene CRISPRi. |
| GSE343369 | GEO summary | KOLF2.2J, 225-gene knockout, embryoid bodies. |
| GSE279287 | GEO summary | 769P CROP-seq, 4 samples, variant screen. |
| GSE286226 | GEO titles | HAP1 enhancer CRISPRi multiome. |
| GSE303901 | GEO title | K562 TAP-seq, 97 genes. |
| GSE85011 | GEO summary | 2016 bulk CRISPRi of 16,401 lncRNA loci in 7 lines. Not single-cell. |
| GSE337804 | catalog title, public 2026-07-10 | “Single cell perturb seq of A549 cells.” This is the Myllia A549 screen already held. Not re-downloaded. |
| Arc blog 2026, Behind the Data 2025, wrap-up 2025 | full text | Line choice and winners, above. |
| Europe PMC abstract, Cell 2026 commentary | abstract | No line names. |
| arXiv 2604.13986, arXiv 2505.14919 HTML | pages | PRiMeFlow and TxPert, above. |
| GitHub API | 2026-09-26 | PRiMeFlow, TxPert, State, CellFlow, GEARS, CPA, Broad depmap repo. |

### Catalogued, too small or the wrong modality (titles from the GEO esummary pull; counts not opened)

The pull is `work/geo_catalog.json` (NCBI esearch/esummary, 2025–2026 Perturb-seq / CROP-seq / CRISPRi / TAP-seq, capped at 200 records per query). Human cancer-keyword hits that are not already in the shortlist:

- GSE295649 / GSE295657 / GSE295658 / GSE295659 — in vivo lncRNA CRISPRi in cutaneous squamous cell carcinoma, plus one scRNA sample. lncRNAs, not a gene panel.
- GSE308682 — K562 CRISPRi of 7 transcription factors.
- GSE281160 — HepG2 CROP-seq of 20 MASLD variants.
- GSE339461 / GSE339463 — RNA fingerprinting; K562 CRISPRi and a K562 drug FlexPlex. Known line.
- GSE299590 — “Clonal dynamics… [CROP-Seq],” 2 samples. Earlier in this hunt the RAW tar was 10 KB (**REPORTED**).
- GSE306853 — pancreatic CAF Perturb-seq, patient-derived fibroblasts. Not a standard cancer line; earlier read of the design was 5 sgRNAs (**REPORTED** for the guide count).
- Mouse and spatial series (GSE348514, Spatial Perturb-Seq in mouse brain, public 24 Sep 2026) are the wrong species.

### Negatives

- No public single-cell gene CRISPRi or CROP-seq atlas for HeLa or CAL-33 (searches above).
- Pisces is empty. Perturb-Sapiens is cytokines.
- Illumina’s Billion Cell Atlas was not found as a new public dump. **REPORTED**: no new accession turned up in the GEO pull or the Hugging Face org listing used here.
- HT29/SW480 enhancer CRISPRi (Parse, dCas9) was identified earlier in this hunt as EGA-controlled (EGAD50000000375). **REPORTED.** Not downloaded; an application would be required, and the targets are enhancers. HT29 gene CRISPRi is already held via Jiang.

### Already held — skipped except where a new release was checked

Replogle K562 genome-wide and essential, RPE1; Nadig HepG2 and Jurkat; X-Atlas/Orion HCT116 and HEK293T; CD4+ T-cell genome-scale CRISPRi; KOLF2.1J; Feng iPSC; Jiang A549, HT29, MCF7, BxPC3, HAP1; Boyle 2026 Broad 16-line; VIPerturb K562; Pan H1 (GSE295214); VCC 2025 H1; Tahoe-100M; DLD-1 low-MOI (now confirmed as the public GSE337988 low-MOI arm); Myllia A549 (GSE337804); Sunshine Calu-3; Wessels THP-1; CMap/LINCS; IGVF cardiac and endothelial screens; scPerturb; Arc State and Stack. No newer genome-scale release of these was found beyond the public GEO landing of DLD-1 and Myllia A549.

---

## What changed

One new cancer-line screen of ≥1,000 targets is actually on the public internet and was opened: A549 Cas9 knockout, 606k cells, median on-target residual 0.58 of control CPM. It does not replace a CRISPRi donor. The genome-scale CRISPRi release that would matter for Jurkat (Pisces) is still empty, and nothing names the October lines.
