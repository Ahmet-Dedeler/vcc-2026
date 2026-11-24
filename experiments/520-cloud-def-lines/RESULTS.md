# 520 — Which lines are D, E, F, and which same-line banks do we lack?

Date: 2026-09-28. Read-only against the rest of the repo. Nothing outside `experiments/520-cloud-def-lines/` was edited. `vcc` was not run. This VM has no gitignored matrices, so coverage numbers below are either recomputed here from tracked or small public files, or cited from an earlier experiment that measured them on the effect tables.

## Decision

**No public Arc document names the three held-out lines.** The Cell commentary and the 20 August 2026 Arc news post say only that Arc ran CRISPRi Perturb-seq (10x Flex, Ultima UG100) on six human cancer lines from different tissues, that three are the validation set, and that three are reserved for the test. The CLI states the labels are opaque on purpose. The practical October fact is stronger than the ranking: **every public genome-scale non-essential cancer CRISPRi bank that covers an Arc-style 300-gene panel is already in `scripts/def_sameline.py`'s registry** (K562 GWPS 272/300, Orion HCT116 290/300, DLD-1 67/300). The lines that are most plausible as D/E/F (A549, MCF7, HepG2, HT-29, BxPC3) are also already registered, and their public banks cover about 0–14 genes of the A/B/C panel. Nothing released 2026-09-20 through 2026-09-28, among the records actually opened, is a new human cancer genome-scale CRISPRi or CRISPR-KO single-cell screen.

Subjective probability that any specific trio is correct is low. Confidence that the high-coverage banks are already registered is high.

## 1. What the sources actually say

Opened 2026-09-28 unless noted.

**Cell commentary**, doi 10.1016/j.cell.2026.08.004, “Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts.” Europe PMC (PMID 42648290) abstract, opened this session: “Participants will build models to predict gene knockdown responses in a new Arc-generated dataset comprising unseen cell lines.” Full text was fetched earlier the same day from the Cell HTML (Cloudflare blocked a repeat fetch; the saved page is the source of the sentences below). Published online 26 August 2026, CC BY-NC-ND 4.0. Relied-on sentences:

- “At Arc Institute, we have performed Perturb-seq experiments on six cell lines originating from different tissues using CRISPR interference (CRISPRi) and 10x Flex single-cell RNA-seq (scRNA-seq).”
- Figure 1C: “The 2026 Virtual Cell Challenge is built on six distinct CRISPRi cancer cell line datasets. Three of these will be used during the validation phase, while the remaining three will be reserved for the testing phase.”
- Supplement listing on that page: “Document S1. Consortium authors” and “Document S2. Article plus supplemental information.” No cell-line table in the HTML full text.
- Declaration: “L.A.G. is a cofounder of nChroma Bio and a scientific advisory board member of Myllia Biotechnology.” Luke Gilbert is a commentary author. That is a conflict note, not evidence that a Myllia line is in the challenge.

**Arc news**, https://arcinstitute.org/news/virtual-cell-challenge-2026, 20 August 2026, fetched this session:

- “At Arc we ran Perturb-seq across six cell lines from different tissues of origin, using CRISPR-interference, 10x Flex single-cell profiling, and high-throughput sequencing on Ultima’s UG100, with perturbations chosen to provide a strong set of perturbations and responses for a robust Challenge. Three of these cell lines will be used for the validation phase and the live leaderboard. The remaining three are held back for final testing and determining our winners.”
- “Models must predict CRISPRi knockdown responses in six cell lines they have never seen perturbed, using only the unperturbed state of those cells and a list of genes to knock down.”
- No line names.

**CLI wiki** (saved fetch of the contexts page): “Each context is a different cell line. We don't publish which one. `A`, `B` and `C` are deliberately opaque labels.”

**2025 Behind the Data**, https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge, fetched this session. This is the H1 design, not the 2026 lines, and it is the reason a genome-scale public line is a worse prior than a line with a thin bank:

- “Unlike previously profiled immortalized cell lines (like K562 or A375 cells), which dominate existing Perturb-seq datasets, the pluripotent H1 ESCs regulatory architecture, chromatin state, and response dynamics represent a true distributional shift relative to most public pretraining data. This choice was strategic: we wanted a dataset that requires models to generalize.”
- Panel construction (the sentence FINDINGS H.1 compresses): targets came from a large-scale CRISPRi pilot in H1, binned by number of DE genes (negligible / subtle / strong), sampled with ContrastiveVI+ for response-mode diversity, plus “clean guide assignment, sufficiently high basal expression of the target gene … and strong target gene repression.” The 2025 post’s knockdown table says “83% of cells with >80% knockdown.”
- “Adequate representation in existing public datasets” in FINDINGS H.1 is about **target genes** (why K562-gwps covers 269/300 and Orion 290–294). It is not a statement that each challenge line has its own Perturb-seq bank.

**Validation identities** (FINDINGS, exp 041 / 356; not re-derived here): A = Jurkat (T-ALL, CCLE haematopoietic and lymphoid), B = HeLa (cervix), C = CAL-33 (head and neck squamous carcinoma). A/B/C are absent from Tahoe-100M (FINDINGS). Jurkat’s public CRISPRi bank is essential-only (0/300 of this panel). HeLa and CAL-33 have no gene-level bank in the registry.

**What that implies, as inference, not as a quoted rule.** “Different tissues of origin” plus the three known CCLE tissues makes a second haematopoietic line (K562, THP-1, HAP1), a second cervix line, or a second upper-aerodigestive line less likely. The 2025 choice to avoid K562 and A375 because they dominate public Perturb-seq makes HCT116 (Orion ~290/300) and K562 (GWPS ~272/300) worse candidates if the 2026 benchmark is built the same way: copying a public genome-scale non-essential screen would not be a zero-shot test. CAL-33 being an uncommon DSMZ line shows Arc will also pick lines that are not the usual Perturb-seq staples, so the prior is flat below the top few.

Jiang et al. use almost the same sentence and a **different** set. From the supplementary/full text fetched this session (bioRxiv 10.1101/2024.01.29.576933; Nature Cell Biology version is the same study): “we aimed to perform Perturb-seq experiments in six different cancer cell lines from different tissues of origin: A549 (lung), MCF7 (breast), HT29 (colon), HAP1 (bone marrow), BxPC3 (pancreas), and K562 (bone marrow).” VCC validation lines are not those six. Treat Jiang as the closest published analogue of the task structure, not as the VCC line list.

## 2. Ranked candidates

Probabilities are a subjective chance that the line is **one of the three**, conditional on the public text above. Lines in the same tissue are mutually exclusive under the tissue reading, so the numbers are not an independent draw and do not sum to 3. Registry status is `scripts/def_sameline.py` `REGISTRY` as of this checkout.

| Rank | Line | Tissue | P(one of D/E/F) | In registry | Best public gene-level bank vs A/B/C panel |
|---|---|---|---:|---|---|
| 1 | A549 | lung | 0.25 | yes (`jiang_A549`, `myllia_a549`) | Myllia CRISPRi 14/300 (exp 225); Jiang ≤9 |
| 2 | MCF7 | breast | 0.20 | yes (`jiang_MCF7`) | Jiang ≤9 |
| 3 | HepG2 | liver | 0.18 | yes (`nadig_hepg2`) | essential CRISPRi, 0/300 on this panel |
| 4 | HT-29 | colon | 0.15 | yes (`jiang_HT29`) | Jiang ≤9; enhancer CRISPRi is not a gene bank |
| 5 | HCT116 | colon | 0.12 | yes (`orion_hct116`) | Orion genome-wide CRISPRi 290/300 (exp 225) |
| 6 | MDA-MB-231 | breast | 0.12 | no | enhancer CRISPRi only; 0 genes by construction |
| 7 | PANC-1 | pancreas | 0.10 | no | none found (Tahoe is chemical) |
| 8 | BxPC3 | pancreas | 0.08 | yes (`jiang_BXPC3`) | Jiang ≤9 |
| 9 | PC-3 | prostate | 0.08 | no | none found this sweep |
| 10 | A375 | melanoma | 0.07 | no | CRISPRa, 2/300 recomputed here |
| 11 | SW480 | colon | 0.06 | no | none found (cultured in the HT29 enhancer paper, not screened) |
| 12 | NCI-H460 | lung | 0.05 | no | none found; Tahoe alternative to A549 |
| 13 | DLD-1 | colon | 0.05 | yes (`dld1_lowmoi`) | 5,000-gene CRISPRi, **67/300 recomputed here** |
| 14 | U2OS | bone | 0.05 | no | Myllia 218-gene Cas9 screen, licensed, not public |
| 15 | BT-474 | breast | 0.04 | no | none found; Tahoe |

Modal trios, if the tissue reading holds and Arc repeats the 2025 “don’t pick the line everyone has memorised” choice: **A549 + MCF7 + HT-29**, or **A549 + HepG2 + HT-29**. If they instead pick lines with no public gene bank, MDA-MB-231, PANC-1, and PC-3 or NCI-H460 fill the same tissue slots. A trio built from HCT116 + K562 + A375 would contradict both the tissue reading and the 2025 design note; that trio is the one our registry would reward most, and it is the one the design is least likely to hand us.

**Down-weighted, not in the top 15.** K562 (same CCLE tissue as Jurkat, and the line the 2025 post named as dominating public Perturb-seq; registry holds GWPS). THP-1 and HAP1 (haematopoietic / near-haploid myeloid; registry holds thin banks). HEK293T (Orion genome-wide, ~294/300 in exp 225, donor exists in-repo as `orion_hek293t`, **not** in `REGISTRY`; not a classic cancer line and not a D/E/F candidate). RPE1 (Replogle essential CRISPRi; not cancer). Calu-3 (Sunshine, registry, 6/300 in exp 225; lung slot already taken by the A549 prior).

**Tahoe does not narrow the list.** Tahoe-100M is a chemical perturbation atlas (50 lines, 47 analysed, 13 organs). A/B/C are absent (FINDINGS). Membership of a line in Tahoe is not evidence it was CRISPRi-profiled by Arc.

## 3. Per-line public single-cell perturbation data

Panel proxy: `data/challenge/pert_counts.csv`, 300 symbols, column `target_gene`. “≤9” for a Jiang line is because exp 225 measured the **union** of the six Jiang lines at 9/300 (`jiang_6line`), so each line is at most 9. Per-line Jiang overlaps were not recomputed here (the Seurat objects are multi-GB). Essential-library 0/300 is exp 225’s measurement on the effect tables, consistent with FINDINGS H.1 (this panel excludes common essentials). Exp 268 measured the Broad 90-target list at **0/300** against the same CSV.

### A549 (registered)

| Dataset | Lab | Modality | Chemistry | Genes | Cells | Licence | Size | URL | Panel |
|---|---|---|---|---|---|---|---|---|---|
| Jiang pathway screens, GEO GSE281048 (public 2024-11-04, 4 samples) | Satija / NYGC | CRISPRi (dCas9-KRAB-MeCP2), stimulated | Parse Evercode WT Mega; Ultima, subset Illumina | 44–61 regulators per pathway × 5 pathways, 3 sgRNAs each | ~2.6 million cells across six lines and two replicates (paper total, not A549-only) | paper CC; GEO public | Zenodo Seurat objects 2.6–5.6 GB each — not downloaded | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE281048 · https://doi.org/10.5281/zenodo.14518762 | ≤9 (exp 225 union) |
| Myllia A549, GEO GSE337804 (public 2026-07-10, 4 samples) | Myllia | CRISPRi (sample titles say “CRISPIRi”) | not stated on the GEO summary | **unverified** (a Zenodo search snippet said 220 genes / 53,511 cells; the GEO page does not) | unverified | GEO public | not downloaded | GEO accession above | 14/300 (exp 225, `myllia_a549`) |
| Comparative CROP-seq, 8 conditions | Myllia | **Cas9 KO**, not CRISPRi | WTA | 218 | “almost 300,000” across 8 conditions, not A549-only | **licensed, not a public download** | n/a | brochure PDF, see manifest | not measurable; not public |

Relied-on Jiang sentences (full text): “for each pathway, we selected 44 to 61 genes”; “We adapted the Perturb-seq workflow to be compatible with the Parse Biosciences Evercode™ Whole Transcriptome Mega kit”; “we applied our combined workflow to sequence approximately 2.6 million cells across two experimental replicates.” GEO esummary title: “Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens.” GEO esummary for GSE337804: “Single cell perturb seq of A549 cells from Myllia Biotechnology.” Sample title: “CRISPIRi guide library from single cell mRNA sequenicing of A549 perturbed cells, channel 1.”

Altair optical pooled screen on A549 is 20/300 in exp 225 and is not CRISPR. Not a same-line Perturb-seq bank.

### MCF7 (registered)

Jiang only, same study as A549. No other gene-level single-cell CRISPR screen was opened for MCF7. Panel ≤9. Enhancer or other MCF7 screens were not re-verified this session.

### HepG2 (registered)

Nadig / TRADE, Nature Genetics (preprint full text PMC11244993, opened). Relied-on sentences: “Jurkat-Essential (2,393 essential gene perturbations in a T-cell leukemia cell line) and HepG2-Essential (2,393 essential gene perturbations in a hepatocellular carcinoma cell line).” Table 1: HepG2-Essential, 2,393 perturbations, **45 cells/perturbation**; Jurkat-Essential, 2,393, **83 cells/perturbation**. Methods: 10x, Cell Ranger 4.0.0, GRCh38 2020-A; HepG2 retained 145,473 cells after QC. GEO GSE264667, public 2024-05-05, 224 samples, esummary title “Transcriptome-wide characterization of genetic perturbations.” Figshare of the Replogle companion (which this paper reanalyses) is 10.25452/figshare.plus.20029387, CC BY 4.0, total ~171 GB; not downloaded. Panel overlap of essential libraries: **0/300** (exp 225). This is the Jurkat analogue: a real same-line CRISPRi bank that does not cover this panel’s gene class.

### HT-29 (registered)

- Jiang, same as A549. Panel ≤9. Registered as `jiang_HT29`.
- Enhancer CRISPRi, not genes. Wellcome Open Research page opened: “we performed genome-wide CRISPR/dCas9-ZIM3-KRAB perturbation of 12,117 candidate enhancer regions in the HT29 CRC cell line, followed by single-cell transcriptomic profiling.” “Ten days post-transduction, we profiled 252,485 single-cell transcriptomes.” Library: 35,139 gRNAs. High MOI. **0/300 by construction** (elements, not genes). GEO/ArrayExpress accession was not on the page extract. **Accession unverified.**

### HCT116 (registered)

X-Atlas/Orion, Xaira, bioRxiv 10.1101/2025.06.11.659105 (posted 16 June 2025), full text opened. Relied-on sentences: “We engineered HCT116 (human colorectal carcinoma) and HEK293T … with dCas9-KRAB-Zim3”; “41,780 sgRNAs targeting 18,903 genes”; “We thawed 4-8 million cells per batch for scRNAseq library generation on the 10x Genomics GEM-X 5’ Platform”; “5,387 and 5,871 median genes per cell on HCT116 and HEK293T respectively.” Not Flex. Exp 225 panel overlap **290/300**. Already `orion_hct116`. The same Myllia licensed 218-gene Cas9 screen includes HCT116 (brochure; not public).

### MDA-MB-231 (not registered)

PMC page “Enhancer regulatory networks globally connect non-coding breast cancer loci to cancer genes” opened. Relied-on sentences: “We performed single-cell CRISPRi screens targeting 3,512 regulatory elements in two diverse breast cancer cell lines. This dataset spans 31,878 sgRNAs (8–10 sgRNAs per regulatory element) in 558,137 cells”; “We sequenced 558,137 cells from (ER+) MDA-MB-361 and (ER−) MDA-MB-231”; “After all quality control filters, 383,513 cells remained”; “on average 1,063 cells for each perturbed region.” dCas9-KRAB. **0 gene-panel overlap by construction.** GEO accession not captured on the extract. **Accession unverified.**

### PANC-1, PC-3, SW480, NCI-H460, BT-474 (not registered)

No gene-targeting single-cell CRISPR screen (CRISPRi, KO, or CRISPRa) was found on a landing page opened this session. Tahoe-100M contains PANC-1, SW480, NCI-H460, and BT-474 as **drug** perturbations. SW480 is named as a culture line in the HT29 enhancer paper’s methods, not as a screened line. **These negatives are absence-of-evidence from this sweep, not a proof that no dataset exists.**

### BxPC3 (registered)

Jiang only. Panel ≤9. Registered as `jiang_BXPC3`.

### A375 (not registered)

Jerby lab, Nature Genetics, published 7 April 2026, “High-content CRISPR activation screens…” (page opened). Relied-on sentences: “a CRISPRa Perturb-seq screen (276 sgRNAs targeting 61 genes or hits…)”; “5′-droplet-based single-cell RNA sequencing”; “20,207 high-quality scRNA profiles from monoculture (12,640 cells) and coculture (7,567 cells)”; “A375 cells to stably express CRISPR–dCas9-VPR.” Figshare 10.6084/m9.figshare.31119196.v1, CC BY 4.0, description: “CRISPR activation Perturb-seq screen in melanoma (A375) cells in monoculture and coculture.” Metadata xlsx downloaded. Guide IDs of the form `GENE_<pos>_23_P1` yield **59 gene symbols** (the paper’s 61 includes guides that may not all appear as that token). Intersection with `pert_counts.csv`: **E2F1, VIM — 2/300.** The bulk CRISPRa screen (32,670 sgRNAs, 2,919 genes) is not single-cell.

A separate “Perturb-ME” A375 MHC-I screen was seen only as a search snippet. **Not opened. Unverified.**

### DLD-1 (registered)

GEO GSE337988, public 2026-07-10, 48 samples. Es summary title: “Joint analysis of multiply perturbed cells improves statistical power and cost efficiency in Perturb-seq screens.” Relied-on summary sentence: “In a 5000-gene Perturb-seq screen across MOI levels, high MOI reduced per-cell costs by 81% while keeping information loss within 1.5-fold of a low-MOI baseline.” Sample titles: “DLD1 CRISPRi sublib2 … MOI Low (~0.3) / Medium (~3) / High (~10).” This session recomputed coverage from the tracked label file `experiments/116-dld1-lowmoi/metadata/lfc_target_labels.txt` (6,451 guide labels, 4,955 unique symbols before the `_P1` suffix): **67/300**. Already `dld1_lowmoi`. RAW.tar is 26.5 GB. Not downloaded.

### U2OS (not registered)

Only the Myllia licensed comparative screen (below). No public download.

### Myllia brochure (covers several candidates; not a public bank)

PDF “Large-Scale CROP-Seq Datasets Available at Myllia Biotechnology, Q1 2026,” downloaded. Page 2: “Datasets available for licensing at Myllia.” Comparative screen (pages 52–54), relied-on lines:

- “Comparative CROP-Seq screen in 6 cell lines: A549, HCT116, Jurkat, K562, THP-1, U2OS” plus “One cell line in 3 conditions: THP-1 monocytes, THP-1 M0 macrophages, THP-1 M1 macrophages.”
- “218 target genes perturbed using Cas9.”
- “Unbiased whole transcriptome amplification.”
- “Almost 300.000 single-cell transcriptomes.”

Also listed, targeted readout not WTA: genome-wide CRISPRi THP-1, panel of 357 mRNAs; genome-wide CRISPRi Jurkat, panel of 374 mRNAs. The Jurkat targeted genome-scale set is the class of GSE249595 (esummary, public 2024-08-22, 16 samples: “Genome-scale CRISPRi Perturb-seq on both unstimulated and stimulated Jurkat cells…”). Exp 166 rejected that readout as a donor because the transcriptome is a targeted panel. **U2OS is the only line in the comparative list that is not already in `REGISTRY`. There is no public URL.**

### Broad DepMap Perturb-seq pilot (16 lines; not D/E/F priors, included because it is the largest new public KO set)

Figshare 10.6084/m9.figshare.33273600.v1, CC BY 4.0, “Dissecting context-dependent cancer vulnerabilities using Perturb-seq - Data.” bioRxiv 10.64898/2026.08.24.746802 posted 2026-08-25 (before the 20 September cutoff). Cas9 KO, not CRISPRi. Metadata zip downloaded. `cell_line_metadata.csv` lines and Oncotree lineages:

HCT15 bowel; KYSE140 esophagus/stomach; KYSE450 esophagus/stomach; RPMI7951 skin; SKGII cervix; SLR23 kidney; U343 CNS/brain; C4I cervix; JVE127 bowel; MG63 bone; ONE58 pleura; T3M5 head and neck; KMRC20 kidney; SNU8 ovary; UMRC3 kidney; HKGZCC biliary.

Exp 268, relied-on sentence: “Overlap of the 90 Broad targets with the 300-target VCC panel: 0 (`data/challenge/pert_counts.csv`, verified).” Exp 283: each reserved line has 100 KOs in the matrix, median cells/KO from 27 (SLR23) to 117 (JVE127). Deep rescreens of some lines used 10x GEM-X Flex (stated in the preprint; not re-quoted from a methods page this session — **Flex-chemistry sentence for the rescreens is unverified here**). `data.zip` (6.9 MB, downloaded) contains gene locations and genesets, **not** the 100-gene target list. Overlap 0 is exp 268’s measurement, not recomputed from a gene list opened today. Code: https://github.com/broadinstitute/perturb-seq-depmap-public. These 16 lines are unusual DepMap lines. SKGII and C4I are cervix (HeLa’s tissue). T3M5 is head and neck (CAL-33’s tissue). They are poor D/E/F priors, and a 0/300 essential-skewed library would not feed same-line mode on an A/B/C-style panel.

## 4. Releases since 2026-09-20

Window checked: 2026-09-20 through 2026-09-28. **No new human cancer genome-scale or large CRISPRi/KO single-cell screen was opened.**

| Source | Query | Result | Relied-on text |
|---|---|---|---|
| GEO esearch, title Perturb-seq OR CROP-seq OR “CRISPR interference”, PDAT that window | `geo_tight.json` | **1 hit: GSE348514** | esearch `count: 1`, id 200348514 |
| GEO esummary of the broader same-day set | `geo_new.json`, 4 series | none is a human cancer screen | GSE348514 (2026-09-24): “Xenium In Situ Spatial Perturb-Seq… mouse brain,” 2 samples, *Mus musculus*. GSE348133 (2026-09-24): lncRNA CYTOR/MORRBID bulk-style KD in MDA-MB-231, not a screen. GSE335154 (2026-09-23): PDAC CAF trogocytosis, not a screen. GSE312628 (2026-09-23): TTC5 iPSC neurons, not a screen |
| Europe PMC | TITLE Perturb-seq, FIRST_PDATE that window | hitCount **0** | `epmc.json` |
| Europe PMC | TITLE CRISPRi OR Perturb-seq OR abstract “genome-scale Perturb-seq” | **1 hit** | “Genome-wide CRISPRi maps pneumococcal attachment and intracellular survival…” doi 10.64898/2026.09.18.752567, first publication 2026-09-21. Bacterial, not human cancer |
| bioRxiv API | interval 2026-09-20:2026-09-28 | 1,435 new papers; **only the first page of 30 was read** | no human cancer Perturb-seq on that page. **Sweep incomplete** |
| Zenodo | most-recent “Perturb-seq OR CRISPRi single-cell” | first records were unrelated reviews | query was too loose. **Not a complete Zenodo sweep** |
| ArrayExpress, Hugging Face | not re-queried this session | — | exp 225 had swept them earlier. **Unverified for this window** |

The Broad pilot (25 August 2026) and GSE337988 / GSE337804 (10 July 2026) are outside this window. They are already in section 3.

## 5. Download manifest — banks the registry does not already use

**Nothing in this list adds a public same-line gene bank with useful coverage of an A/B/C-style panel.** Do not download the multi-GB matrices for same-line mode.

| What | URL | Size | Hash | Why it does not fill a hole |
|---|---|---|---|---|
| A375 metadata (downloaded) | https://ndownloader.figshare.com/files/61206028 | 1,872,823 B | sha256 `972a111575f5f033433ec4744ab997ccbee9394b2a8f5e4fb26bf04e461dbd08`; figshare md5 `040029c3db2baa929ef902169fd9ce01` | 2/300. Counts mtx is 1,496,011,467 B (over the 1 GB cap); Seurat rds 5,379,051,121 B. Not worth it |
| Broad metadata.zip (downloaded) | https://ndownloader.figshare.com/files/67643358 | 1,838,227 B | sha256 `c7af8e03fb719efe350918d61ceeeaa85944023d33a92de854603ef98c3c360e`; figshare md5 `11ad26cc9c5ab272b80b20e0b97f86e0` | line list only |
| Broad data.zip (downloaded) | https://ndownloader.figshare.com/files/67643355 | 6,891,220 B | sha256 `860b07aa03c76bf65fff0d2a9c8fb861de4817af04f32f7149b9bacd8a27fc94`; figshare md5 `6bbe267ba864ef7ede1e15847bb6af90` | genesets, not the KO list. Exp 268 overlap is 0/300 |
| Broad `single_cell_data.zip` (not downloaded) | https://ndownloader.figshare.com/files/67643568 | 3,171,208,158 B; figshare md5 `78b6506f8cb40dd111ccc8480641ecb3` | — | over the 1 GB cap; 0/300 |
| Broad `h2h_data_formatted.zip` | https://ndownloader.figshare.com/files/67643532 | 1,835,216,779 B; md5 `8aa4600fbd27478841340fb9f56d9d70` | — | not downloaded |
| Broad `processed.zip` | https://ndownloader.figshare.com/files/67643565 | 3,233,927,690 B; md5 `8ad0b34b854add4b28c29f7425470894` | — | not downloaded |
| Broad `downloaded.zip` | figshare file id 67643361 area of the same record | 461,671,192 B; md5 `64d39184ee3edc6c815f84828f807dc8` | — | not downloaded; still no panel overlap to gain |
| Myllia Q1 2026 brochure (downloaded) | https://myllia.com/resources/large-scale-crop-seq-perturbation-datasets-available-at-myllia/ (HTML says the sets are available for licensing; the six line names are in the PDF, not on that HTML) | 8,559,015 B | sha256 `900c341294f65de1ce51012e5126b9646f5ec9406b6b45af1a017b65b31be726` | U2OS would be a new line. **No public count matrix** |
| HEK293T Orion | already a donor path in-repo (`orion_hek293t`), not in `REGISTRY` | not re-downloaded | — | not a cancer-line candidate. Exp 225: 294/300 |

Registered banks that already exist and should not be re-fetched for this question: Orion HCT116, K562 GWPS / essential / VIPerturb, DLD-1 low-MOI, Nadig Jurkat, Nadig HepG2, Jiang A549/MCF7/HT29/BxPC3/HAP1, Myllia A549 (GSE337804), Sunshine Calu-3, Wessels THP-1. Sunshine and Wessels primary pages were **not re-opened this session**; their panel counts (6/300 and 1/300) are exp 225’s, and both keys are already in `REGISTRY`.

## 6. What to do with this on 22 October

Same-line rows turn on only when `def_identity.py` is unanimous and the CCLE name is in `REGISTRY` (DEF-RUNBOOK). If D, E, or F is A549, MCF7, HT-29, HepG2, BxPC3, HAP1, Calu-3, or THP-1, the bank is already wired and covers a handful of an essential-free panel (HepG2/HAP1/THP-1: about zero on this panel’s gene class). If it is HCT116, K562, or DLD-1, the bank is also already wired and is the only case with real panel coverage — and those three are the less likely identities under the tissue-plus-zero-shot reading. If it is MDA-MB-231, PANC-1, PC-3, A375, U2OS, or a Broad-16 line, there is no public gene-level bank to add, and downloading A375 or the Broad matrices would not change that.

The one registry gap that is real but irrelevant to D/E/F is HEK293T Orion: the donor exists and is not a `REGISTRY` key. Do not add it for October unless the identity call returns a HEK293T label, which a cancer-line panel should not.
