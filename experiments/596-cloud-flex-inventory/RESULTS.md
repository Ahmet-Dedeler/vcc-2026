# 596 — Public 10x Flex profiles of named human lines

2026-10-01. Inventory only. No matrices were downloaded, no `vcc`, no edits outside this directory.
Facts are tagged **VERIFIED** (the page, API record, sample table, or file listing was opened on 2026-10-01)
or **UNVERIFIED** (a field the opened record did not state). 10x Genomics dataset HTML is behind a
Cloudflare challenge from this host (HTTP 429); those rows use the indexed page text returned by search,
and anything that text omitted stays UNVERIFIED. Direct `cf.10xgenomics.com` HEAD requests returned 403,
so 10x file sizes are UNVERIFIED.

## 1. Verdict

**Eight datasets can support a control-pseudobulk identity test, covering eight distinct named human
lines.** Five of those datasets carry an explicit CC BY 4.0 licence. Two more are public GEO / ArrayExpress
records with no CC licence named. One more (EndoC-βH1) is Flex but not a basal culture.

| Line (as written) | CCLE / DepMap | In a clean unperturbed Flex profile? |
|---|---|---|
| Jurkat, ATCC TIB-152 | `JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE` (590's call for context A) | Yes. 10x's own 145,628-cell culture. Same line as validation context A, different lab. |
| K562-r, ATCC CRL-3344 | Parental symbol K562 (`K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE`). The resistant clone is not a separate CCLE entry that this inventory re-checked. | Yes, but it is the imatinib-resistant clone of K562, not the parental line. |
| A375 | A375 (melanoma; full CCLE sample string not re-queried against the local cache) | Yes. Flex v2, one population split across probe barcodes. The deposited object is described as normalized counts. |
| MCF7 | MCF7 (breast; full CCLE sample string not re-queried) | Yes, 7,127 QC cells, whole-transcriptome Flex plus 24 custom fusion probes. Multiplexed in one capture with paediatric B-ALL samples. |
| PC-9 | PC-9 / PC9 (lung adenocarcinoma; full CCLE sample string not re-queried) | Yes, one human-probe-only culture. Also a vehicle-treated xenograft. |
| REH, DSMZ ACC 22 | REH (B-ALL; full CCLE sample string not re-queried) | Vehicle only: 0.9% DMSO for 72 h, not an untreated culture. |
| K562-CRISPRi (KRAB-dCas9) | Same K562 symbol, engineered | Non-targeting cells only, in three public screens (10x GEM-X Flex, 10x Ultima, Satija VIPerturb-seq). |
| hTERT RPE-1 | Not a cancer line. 590's 3′ RPE1 profile was called SW982 at 0.466. | Only the ORF-unassigned cells inside a doxycycline-induced screen. |
| HEK293 and HAP1 (CRISPRi) | Not CCLE cancer lines. 590 called HEK293 → HAP1_ENGINEERED at 0.554–0.606. HAP1 is in the Jiang registry the matcher already uses. | Yes, as hashing-separated non-targeting cells in one Flex v1 pilot. This is the paired "line absent from CCLE" case. |

**No public Flex profile identifies a final-test line.** Jurkat is already validation context A, so the 10x
Jurkat set is a same-chemistry check of a known line, not a candidate for D, E, or F. HeLa and CAL33
(the other two validation lines) have **no** public Flex profile in this search. MCF7, K562, HAP1, and
RPE1 are lines that already have public perturbation data of other chemistries; a Flex profile of them
tests the 0.65 floor, it does not say the final panel is that line.

The only Flex profiles already on disk remain the challenge's three validation contexts and the H1
embryonic-stem-cell screen. H1 is not a cancer line and is not re-listed below as a new find.

## 2. Datasets, most useful first

### 2.1 Usable for a basal (or non-targeting) pseudobulk

**U1. 10x Genomics — 128k Jurkat cells, 16 probe barcodes (Next GEM Flex).**
- Landing page: https://www.10xgenomics.com/datasets/128k-jurkat-cells-multiplexed-samples-16-probe-barcodes-1-standard
  (indexed page text VERIFIED; live HTML 429). No GEO accession on that text (UNVERIFIED whether one exists).
- Line: Jurkat, ATCC TIB-152. CCLE: `JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE`.
- Chemistry: Chromium Fixed RNA Profiling, multiplexed-samples user guide CG000527; fixation CG000478, 1 h at room temperature. Product labelled Flex, Next GEM. Human whole-transcriptome probe set (part number not stated). 16 probe barcodes are **sub-pools of one culture**, not 16 lines.
- Cells: **145,628 detected** (VERIFIED, indexed page). All of them are the unperturbed line.
- Controls: the whole dataset is unperturbed.
- Files: Cell Ranger `multi` outputs (filtered feature-barcode matrix in H5 and MEX, per-sample outs, web summary, Loupe). Sizes UNVERIFIED. Sequenced ~35k read pairs/cell, NovaSeq 6000, 28/10/10/90.
- Licence: **CC BY 4.0** (VERIFIED, indexed page).
- Release date: UNVERIFIED (not in the indexed text that was returned).
- Caveat: none for an identity test. This is the closest public stand-in for context A (same line, same probe chemistry, different lab and a larger cell count).

**U2. 10x Genomics — 10k human K562-r cells, singleplex (Next GEM Flex).**
- Landing page: https://www.10xgenomics.com/datasets/10k-human-k562-r-cells-singleplex-sample-1-standard
- Line: K562-r, ATCC CRL-3344. ATCC describes it as an imatinib-resistant clone of parental K562, paired with K562-s (VERIFIED, https://www.atcc.org/products/crl-3344).
- Chemistry: Fixed RNA Profiling singleplex **with Feature Barcode technology for cell-surface protein**, user guide CG000477; fixation CG000478, 1 h at room temperature; Chromium X. Whether an antibody panel was actually stained is not stated (UNVERIFIED). The gene-expression matrix is the whole-transcriptome probe set.
- Cells: **11,645 detected** (VERIFIED). One library, one line.
- Controls: unperturbed (no drug treatment stated).
- Files: Cell Ranger outputs as above. Sizes UNVERIFIED. ~33k read pairs/cell, NovaSeq 6000.
- Licence: **CC BY 4.0** (VERIFIED).
- Release date: UNVERIFIED.
- Caveat: a resistant subclone can move the transcriptome relative to parental K562. A call of K562 below 0.65 would be about the clone, not about the floor.

**U3. Weir, Krebs, Chen (Broad) — A375, Flex v2. Zenodo 19363777.**
- Landing: https://zenodo.org/records/19363777 — README.txt (2.7 KB) and the API record opened (VERIFIED). DOI 10.5281/zenodo.19363777.
- Line: A375 melanoma. CCLE symbol A375.
- Chemistry: **10x Flex v2**, as named in the README. One population split across probe barcodes (the paper's point is barcode-driven false DE).
- Cells: not stated in the README (UNVERIFIED count).
- Controls: unperturbed; barcodes are technical splits. For an identity test, sum every barcode.
- Files (VERIFIED, API sizes):
  - `seurat_A375_flexv2_filtered.qs` — 5,830,989,043 bytes. Seurat v5. README: normalized counts in the default RNA assay; metadata column `probe_barcode`.
  - `flex_v2_probe_expr_sparse_cache_10k_feat_filt.rds` — 386,881,722 bytes. Per-probe (not per-gene) matrix, top 10,000 features only.
  - `vcc_probe_expr_sparse_cache_10k_feat_filt.rds` — 1,231,755,232 bytes. **This is Virtual Cell Challenge cells, Flex v1, top 10,000 probes.** It is not a new line. Do not use it as an external test.
- Licence: **CC BY 4.0** (VERIFIED). Published 2026-04-01.
- Caveat: the gene-level object is described as normalized. 590's matcher wants a sum of raw counts. Whether a raw `counts` layer is also stored was not opened (UNVERIFIED). The per-probe cache is the wrong axis (probe IDs, top 10k).

**U4. Oshlack lab — MCF7 Flex validation, plus two B-ALL patient cohorts. Zenodo 22092839.**
- Landing: https://zenodo.org/records/22092839 (VERIFIED). Preprint https://www.biorxiv.org/content/10.64898/2026.08.26.747171v1 (VERIFIED). Code https://github.com/Oshlack/flex-fusion-probes (VERIFIED).
- Line: MCF7. The other samples are paediatric B-ALL bone marrow, not cell lines.
- Chemistry: standard 10x Flex whole-transcriptome probe set (">18,000 curated genes") **plus custom junction probes** (24 fusions in MCF7). Kit generation (Next GEM vs GEM-X vs Flex v2) is not named in the preprint text that was read (UNVERIFIED).
- Cells: **7,127 high-quality MCF7 cells** after QC (VERIFIED, preprint). MCF7 was co-multiplexed with three B-ALL samples in one capture (Pool-240917-FLEX is the pool name in the Zenodo description; which probe barcode is MCF7 is not in that description — UNVERIFIED until the tar's sample names are listed). Cohort 2 captures do not include MCF7.
- Controls: the MCF7 culture is unperturbed aside from carrying its endogenous fusions. Ambient MCF7 fusion UMIs were detected in B-ALL cells; the preprint treats that as ambient, not as mixed identity.
- Files: `flex_fusion_probes_data.tar` — 2,894,115,840 bytes (VERIFIED). Contains, per pool, the raw Cell Ranger matrix and filtered per-sample matrices. FASTQ is not in the deposit.
- Licence: **CC BY 4.0** on the Zenodo record (VERIFIED). Published 2026-08-25. Patient FASTQ is restricted; the count matrices in the tar are part of the public record.
- Caveat: demultiplex MCF7 before summing. Custom fusion features must be dropped before the sum so they are not treated as genes.

**U5. Single Cell Discoveries / AstraZeneca — PC-9 in vitro, human probe set. E-MTAB-13962.**
- Landing: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13962 (API and the 9.1 KB study TSV opened, VERIFIED). ENA project ERP159157. Paper: Llorà-Batlle et al., Genome Biology, https://doi.org/10.1186/s13059-026-04092-0 (preprint 10.1101/2024.01.25.577066).
- Line: PC-9 (lung adenocarcinoma). NIH/3T3 is the mouse partner in the other arms.
- Chemistry: 10x Gene Expression Flex, human WTA probes PN-2000510, fixation 16–24 h at 4 °C (PN-2000517), Chromium X, cellranger 7.1.0, GRCh38-2020-A. Next GEM Flex era, not GEM-X and not Apex.
- The pure-human arm, from `sample_metadata_pilt_g012.txt` (401 bytes, VERIFIED): source name **PC9**, library PS-PILT-g012, probe barcode **BC001**, matrix file **PS-SCDI-009n**, description "Human cells profiled with human WTA probe set". The other libraries are mixtures (25/50% NIH/3T3, or human+mouse probe mixes) and are the wrong input for a line pseudobulk.
- Cells: not stated in the sample table (UNVERIFIED). Sequencing depth stated as 15,000 reads for the spiked/PC9 samples (protocol text).
- Files (VERIFIED, directory listing): `PS-SCDI-009n.tar` 83 MB, plus the other mixture tars (40–777 MB). ENA runs for g012: ERR12826688 and ERR12826690, first public 2024-04-16, FASTQ ~2.5 GB + ~4.8 GB each. Release date on the study: 2024-04-15.
- Licence: **no licence field** on the BioStudies record (VERIFIED by absence). Public ArrayExpress release. Not a CC grant.
- Caveat: use only g012 / BC001 / PS-SCDI-009n. A file named `sample_metadata_pilt_g023.txt` on the same FTP actually describes an LXFS xenograft (library PS-PILT-g027), which disagrees with the SDRF row for g023. Do not infer g023's contents from the SDRF alone.
- The matched xenograft study **E-MTAB-13963** (VERIFIED) is PC-9 grown in female SCID mice, vehicle vs osimertinib 25 mg/kg for 10 days, human and mouse probes on different barcodes (vehicle: BC002 human, BC006 mouse; matrix PS-ASZE-022n.tar, 35 MB). In-vivo stromal admixture and drug treatment. Useful only as a stress test, and only the human barcode of the vehicle sample.

**U6. Nordlund lab — REH, Fixed RNA Profiling, drug barcode panel. GSE229617.**
- Landing: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229617 (SOFT brief opened, VERIFIED). PMID 38288374, PMC10823582 (paper HTML opened). Public 2024-01-04.
- Line: REH, DSMZ ACC 22, glucocorticoid-resistant E/R+ pre-B ALL.
- Chemistry: Chromium Fixed RNA Profiling, probe set named in the sample record as `Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv`, cellranger 7.0.1, four-barcode human transcriptome kit, Chromium iX. The same series also has Parse Evercode and MULTI-seq (10x 3′ v3.1) arms; those are not Flex.
- Cells: target 2,500 per condition, 7,500 per library. After QC, **9,586 Flex cells** across six conditions; library F-2 (imatinib, fludarabine, DMSO) had 5,329 cells called (Table 1). The DMSO share of those 5,329 is not in Table 1 (UNVERIFIED).
- Controls: **0.9% DMSO, 72 h**, probe barcode BC002, file `GSM7167646_fixed2_Control.tar.gz`. Not untreated. The other barcodes are drugs (5-azacytidine, dexamethasone, XRP44X in pool 1; imatinib and fludarabine in pool 2).
- Files (VERIFIED, `filelist.txt`): `GSM7167646_fixed2_Control.tar.gz` 64,472,876 bytes inside `GSE229617_RAW.tar` (757,155,840 bytes). Contents are Cell Ranger `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`.
- Licence: GEO public release. No CC licence on the series record.
- Caveat: the paper states this v1 probe set has **no ribosomal-protein probes** (98.7% of mapped reads exonic; zero ribosomal-protein counts). Those genes will be NaN on the challenge axis. DMSO for 72 h is a vehicle, and the paper's own DE shows fludarabine moves the cells while dexamethasone barely does; the DMSO arm is the one to use.

**U7. 10x Genomics — ~1.2 million CRISPRi K562 cells, GEM-X Flex, 16 probe barcodes.**
- Landing: https://www.10xgenomics.com/datasets/16-plex_GEM-X_Flex_1M_human_K562_CRISPR_aggregate (indexed text VERIFIED).
- Line: K562 stably expressing KRAB-dCas9.
- Chemistry: GEM-X Flex multiplex, user guide CG000787; fixation CG000782, **24 h at 4 °C**; custom sgRNA probes spiked in (CG000621). Whole transcriptome plus a CRISPR feature library. 16 probe barcodes are sub-pools of one screen, run across four GEM lanes. Sequenced NovaSeq X, 10,000 GEX reads/cell and 5,000 CRISPR reads/cell.
- Library: 6,903 sgRNAs, of which **177 are non-targeting**. 4,045 sgRNAs have named targets (280 lncRNAs, 567 protein-coding genes); 2,681 are "ignore".
- Cells: title says approximately 1.2 million (the called-cell total was not in the indexed paragraph; UNVERIFIED exact count). Non-targeting cell count was not stated. If guides were even, 177/6,903 of cells would be non-targeting; that fraction is a design estimate, not a count.
- Controls: cells whose CRISPR assignment is one of the 177 non-targeting guides. Perturbed cells are not a basal profile.
- Files: aggregated Cell Ranger outs, including the CRISPR feature matrix needed to select non-targeting cells. Sizes UNVERIFIED.
- Licence: **CC BY 4.0** (VERIFIED, indexed page). Release date UNVERIFIED.
- Sibling, same screen, different sequencer: **320k Human CRISPRi K562, Ultima**, https://www.10xgenomics.com/datasets/320k_K562_Flex_CRISPR_Ultima — GEM-X Flex, Cell Ranger 9.0.1, Ultima UG 100, 8,789 GEX reads/cell. **Date published 2025-05-28** (VERIFIED). CC BY sentence was not in the indexed text retrieved for this page (licence UNVERIFIED). Same 177 non-targeting guides.
- Not Flex, do not use: "Aggregate of ~84k Human CRISPRi K562" (https://www.10xgenomics.com/datasets/84k_K562_5pv3_GEMX_CRISPR_aggregate) is **GEM-X 5′ v3** with CRISPR feature barcodes (CG000735), CC BY 4.0, same biological screen. Indexed page text VERIFIED as 5′, not Flex.

**U8. Satija lab — VIPerturb-seq. Zenodo 18460279.**
- Landing: https://zenodo.org/records/18460279 (API record VERIFIED). DOI 10.5281/zenodo.18460279. Preprint https://www.biorxiv.org/content/10.64898/2026.02.12.705613v1 (full text VERIFIED). Site https://satijalab.org/viperturb-seq/.
- Lines: monoclonal **K562, HAP1, and HEK293**, each expressing KRAB-dCas9-MeCP2. HEK293FT (Thermo R70007) was used only to make virus, not as a profiled line. K562 was kept in IMDM + 10% FBS; HAP1 in the same; HEK293 in high-glucose DMEM + 10% FBS.
- Two experiments:
  1. **Pilot, Flex v1.** A pool of all three lines, 100-sgRNA mTOR library plus non-targeting guides, custom barcode probes, intracellular pRPS6 ADT, and **hashing (HTO)** so the lines can be split. Fixation in the methods is 1% formaldehyde for 10 min, not the standard 10x fixation protocol. File `multimodal_cell_line_mixing_pilot.rds`, 1,579,169,440 bytes. Per-line cell counts are not in the Zenodo description (UNVERIFIED). The preprint reports >95% of HEK293 cells assigned a guide (87% singlets) and 6,700 RNA UMI/cell at 10,300 reads/cell.
  2. **Genome-wide, K562 only, GEM-X Flex v2** (user guide CG000834, fixation CG000782, custom probes CG000839). Non-targeting cells are in every deposited bin. Files: `genome_wide_binA.RDS` 3.56 GB, `binB` 3.78 GB, `binC` 2.89 GB, `genome_wide_filtered.rds` 3.61 GB, `genome_wide_manifest.txt` 398,905 bytes. The description says each bin contains **all non-targeting cells** plus a subset of perturbations, so the NTC cells must be taken from one bin, not summed across bins. `vimentin_screen.rds` (1.88 GB) is a FACS-enriched VIM-high population (12,000 gRNA singlets) and is the wrong basal.
- Licence: **CC BY 4.0** (VERIFIED). Published 2026-02-02.
- Caveat: CRISPRi machinery is present in the non-targeting cells. The pilot's fixation is non-standard. HEK293 and HAP1 are the absent-from-CCLE pair; K562 NTC cells are the in-CCLE control in the same object.

**U9. Weissman lab — hTERT RPE-1, viral microprotein ORF screen. GSE344535.**
- Landing: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE344535 (SOFT brief and GSM9981026 opened, VERIFIED). Public 2026-08-21. BioProject PRJNA1516028.
- Line: hTERT RPE-1, immortalized retinal pigment epithelium, not a cancer line.
- Chemistry: Chromium Single Cell Gene Expression Flex, 16-plex, GEM-X modifications in technical note CG000814. Human whole-transcriptome probes (PN-2001249 through PN-2001274) plus a custom ORF probe panel (9,874 targets). Fixation 16–24 h at 4 °C. Eight GEM wells, each a GEX library and an ORF feature library.
- Cells: design target ~320,000 per GEM well (VERIFIED, sample protocol). Called-cell total not in the SOFT record (UNVERIFIED). Infection was 30% (300mer) and 17% (400mer), so most cells never received an ORF, but the whole culture saw doxycycline (24 h or 48 h; the GEO record still contains a `<FILL>` note on the dox concentration).
- Controls: cells with no ORF feature assignment, if the feature table separates them. There is no uninduced, untransduced library in the deposit (an untransduced well is mentioned only as a culture control).
- Files (VERIFIED, `filelist.txt`): `GSE344535_RAW.tar` 28,851,425,280 bytes. Per GEM well, a `sample_filtered_feature_bc_matrix.tar.gz` of 3.17–3.97 GB and a `feature_df.csv.gz` of 27–53 MB. Eight GEX wells: GSM9981026, 028, 030, 032, 034, 036, 038, 040.
- Licence: GEO public. No CC licence on the series record.
- Caveat: doxycycline on every profiled cell, and only the ORF-negative subset is even approximately basal. This is still the only public Flex RPE-1 profile found.

### 2.2 Flex, but not a usable basal cancer-line profile

**N1. EndoC-βH1 xenograft. GSE277402 / GSM8521718.** Public 2025-05-26 (VERIFIED, SOFT brief). NLuc-EndoC-βH1-mCherry, Chromium Fixed RNA Profiling multiplex (4-barcode human transcriptome kit), but the Flex arm is cells treated with an LIP cocktail, transplanted subcutaneously into SCID-beige mice, recovered after a week, and sorted on mCherry. Not basal. The other two samples in the series are 3′ CellPlex, not Flex. Engineered beta-cell line, not in CCLE, and not a cancer line.

**N2. Mouse cell-line mix, 10x.** https://www.10xgenomics.com/datasets/40k-mixture-of-mouse-cell-lines-multiplexed-samples-4-probe-barcodes-1-standard — A20, NIH/3T3, Neuro-2a, EL4; 41,082 cells; Flex; published 2022-11-30 (VERIFIED, indexed page). Mouse probe set. Useless for a human CCLE match.

**N3. 10x Flex primary / FFPE sets** (indexed pages or CELLxGENE API, VERIFIED as not cell lines): NSCLC and adjacent lung (~900k, 16 probe barcodes, flash-frozen tissue); PBMCs; kidney nuclei; 320k FFPE from eight tissues; ovarian cancer FFPE dissociated tumor; breast, DLBCL, and lung FFPE tumors on CELLxGENE (collection `bd552f76-…`, 10x Next GEM Flex v1). No named line.

**N4. GSE318840** osteosarcoma Fixed RNA Profiling is patient tumors across treatment stages (public 2026-06-30), not a cell line. **GSE295527** / CELLxGENE collection `398e34a9-…` is a PBMC platform bake-off that includes Flex (7,750 cells in each of `10X_FRP-rep1` and `rep2`) plus Parse, Scale, and PIP-seq. Not a cell line.

### 2.3 Explicitly not Flex

- **X-Atlas/Orion** (https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion; preprint 10.1101/2025.06.11.659105, VERIFIED). HCT116 and HEK293T, genome-wide Perturb-seq. Chemistry is FiCS: fix, cryopreserve, then **GEM-X 5′** with direct-capture sgRNA, not a probe-based Flex library. Non-targeting cells exist and would test the floor on a different chemistry. Out of scope here.
- **Parse Evercode, Scale Biosciences, Fluent PIP-seq, BD Rhapsody, HIVE.** Probe-based or partitioned fixed-cell alternatives appear in GSE295527 and in the CELLxGENE assay list. They are not Flex. The REH series' Parse arm is the same cells on a different chemistry and is the wrong comparison for a Flex floor.
- **Replogle K562/RPE1, Nadig, the challenge H1 screen.** Already in the 590 read. H1 is Flex but embryonic, and it is already on disk. Not a new cancer-line find.
- **Zenodo 11283842** (A375, HAP1, THP1, HEK293FT, K562, A549 mixtures) is poly-A bulk plus single-cell, and the record does not say Flex. Not claimed as Flex.

## 3. Proposals (not run)

These are for the CC BY datasets only. GEO and ArrayExpress rows in §2 are public but do not name a CC licence; the same procedure would apply to REH DMSO (`GSM7167646_fixed2_Control.tar.gz`), PC-9 g012 (`PS-SCDI-009n.tar`), and RPE-1 ORF-negative cells, after someone confirms reuse.

The scoring path is the one `experiments/590-identity-validation/match590.py` already uses. Do not refit anything.

1. Sum raw integer counts for the control cells of one line. Map features to gene symbols. Drop the target gene's own row only later, inside the matcher (it already does). Drop CRISPR, ORF, fusion, and antibody features before the sum.
2. Place that sum on `data/challenge/gene_names.csv` (18,533 symbols, column `gene_name`). A symbol the probe set does not contain stays NaN. Do not impute ribosomal-protein genes the v1 probe set omitted.
3. Turn the sum into the query the matcher expects, which is what `counts_query` in `match590.py` does: log2(1 + 1e6 × count / sum of finite counts) on the genes that are present, then restrict to the official ∩ CCLE axis (`def_identity.load_cache`, 17,900 genes in the 590 run).
4. Score with the eight methods, the three anchors, and the Borda consensus by importing `scripts/def_identity.py`, the same way `match590.py` does. Read `min_top1_score` (the minimum of the eight top-1 correlations). The preregistered gate is: registered, unanimous, and `min_top1_score` ≥ 0.65.
5. Write the row next to 590's A/B/C Flex values (0.772 / 0.719 / 0.733) and next to the two false HEK293 calls (0.554 and 0.606 plain; 0.506 and 0.561 with the simulated Flex bias).

Per dataset:

| Dataset | Which cells | Expected reading |
|---|---|---|
| U1 Jurkat | All 145,628 cells, all 16 probe barcodes summed (they are one culture). | Call `JURKAT_…`, unanimous. A `min_top1_score` near context A's 0.772, and at least ≥ 0.65, is the same-chemistry positive control. A value under 0.65 would mean the floor does not survive a second Flex lab. |
| U2 K562-r | All 11,645 cells. | Call K562 if the resistant clone is close to the parental CCLE profile. Record the score either way; do not treat a miss as a floor failure. |
| U3 A375 | All cells, all `probe_barcode` values, **only if a raw counts layer exists**. If the object is normalized only, stop; do not feed normalized values to `counts_query`. | Call A375. In-CCLE positive control on Flex v2, which the challenge contexts are not. |
| U4 MCF7 | The filtered per-sample matrix named MCF7 inside `flex_fusion_probes_data.tar` (Pool-240917-FLEX). Drop the 24 fusion features. | Call MCF7. Confirms a carcinoma line that is not Jurkat. |
| U7 K562 CRISPRi, and U8 genome-wide bins | Cells assigned to a non-targeting guide. One bin only for U8, because every bin repeats the NTC cells. | Call K562. The cells carry KRAB-dCas9, so a weaker score than U2 is informative rather than a contradiction. |
| U8 pilot, split by HTO | Non-targeting cells of each line. | **K562:** in-CCLE positive control inside the same object. **HAP1:** should call HAP1_ENGINEERED if that sample is in the matcher's reference, at a score that can be compared with the gate. **HEK293:** the absent-from-CCLE case. 590's 3′ HEK293 profiles were unanimous for HAP1_ENGINEERED at 0.554 and 0.606. The floor holds on real Flex if this HEK293 `min_top1_score` stays under 0.65. If it clears 0.65, the simulated-bias transfer in 590 was not conservative. |
| U9 RPE-1, only if used | ORF-unassigned cells, and say in the row that they were doxycycline-exposed. | Second absent-from-CCLE case. 590's 3′ RPE1 called SW982 at 0.466. A Flex value ≥ 0.65 would be a floor failure on a non-cancer line; a low false call would agree with 590. |

Do not score the VCC probe cache in U3, the VIM-high VIPerturb arm, the PC-9/NIH-3T3 mixtures, the osimertinib xenograft, or the REH drug arms other than DMSO.

## 4. What was searched, and what was empty

Searched 2026-10-01. Empty means no named human cancer or immortalised line on Flex, not that the portal has no Flex data at all.

| Source | Query | Result |
|---|---|---|
| GEO DataSets, `esearch`/`esummary` | `"10x Flex"[All Fields]` OR `"10X Flex"` OR `"10x Genomics Flex"`: 763 records, every title containing Flex/fixed-RNA scanned. `"Fixed RNA"[All Fields]` OR `"fixed-RNA"` OR `"fixed RNA profiling"`: 220 records, same scan. | Phrase `"Fixed RNA Profiling"[All Fields]` returned **0** (the words are not indexed as a phrase). Cell-line series that survived the title/summary scan are GSE344535, GSE229617, GSE277402. The rest are tissue, PBMC, or mouse. |
| GEO sample titles | The Flex-titled GSM list from those 983 records. | REH (GSM7167645/46), EndoC-βH1 (GSM8521718), RPE-1 (the GSE344535 samples). No HeLa, CAL33, A549, HCT116, MCF7, or HepG2 GSM title in that scan. MCF7 and A375 are on Zenodo, not GEO. |
| 10x datasets | Site search for Flex plus Jurkat, K562, HEK293, A549, HeLa, Raji, MCF7, HCT116, 293T, "cell line". | Human cell-line Flex hits are U1, U2, U7, and the Ultima sibling. The 84k K562 CRISPR aggregate is GEM-X 5′. Mouse four-line mix is Flex but mouse. No HeLa, A549, MCF7, HCT116, HEK293, or Raji Flex demo page. Live HTML could not be opened (Cloudflare 429), so a facet browse of "Biomaterial = Cell Line" was not completed. A later sweep that can pass that challenge should export the facet rather than repeat these name queries. |
| CELLxGENE curation API `GET /curation/v1/datasets` | All 2,237 datasets. Assay labels kept if they were `10x gene expression flex`, `10x Next GEM Flex v1`, `10x GEM-X Flex v1`, or `10x Flex Apex` (17 datasets). | None is a named cell line. They are PBMC (aging cohort; the commercial-platform bake-off), endometrium, olfactory epithelium, CTCL skin (Flex mixed with 3′/5′), and FFPE breast, lung, and DLBCL. |
| Human Cell Atlas Azul | `organPart is "cell line"` | **0 projects.** A `projectTitle contains "Flex"` filter was rejected by the schema (`is` is required). Not re-queried with a valid schema. |
| Single Cell Portal API | `type=study&terms="Chromium Flex"` | **0 studies.** An unquoted "Fixed RNA Profiling" query tokenised into three words and returned 492 unrelated studies; that result was discarded. |
| Zenodo API `type=dataset` | `"10x Genomics Flex" OR "Fixed RNA Profiling" OR "Chromium Flex"` | 9 hits. Cell-line hits: 19363777 (A375), 22092839 (MCF7). The other seven are intestine, liver slices, colon, rhabdomyosarcoma, glioblastoma, dermatomyositis skin, or an unrelated 2020 barcode paper. A looser `"cell line"` query returned 595 unrelated records and was not useful. |
| figshare API | `search_for=Chromium Flex Gene Expression cell line`, `item_type=3` | **0 articles.** |
| BioStudies / ArrayExpress | `query=Fixed RNA Profiling cell line` | The endpoint tokenised the query and returned **953,300** hits; page 1's only Flex-titled study was E-MTAB-14782 (primary HPAEpiC air-liquid interface, not a cancer line). E-MTAB-13962 and E-MTAB-13963 were found from the paper, not from this ranking. A phrase-constrained repeat is worth doing; this query is not. |
| Arc Virtual Cell Atlas / 2025 H1 | Not re-downloaded. | Already held. H1 hESC, Flex, not a cancer line. |
| X-Atlas / Orion | Preprint 10.1101/2025.06.11.659105 and the Hugging Face dataset card. | HCT116 and HEK293T are FiCS on GEM-X 5′, not Flex. Excluded on purpose. |
| Names with no Flex hit in any of the above | HeLa, CAL33, CAL-33, A549, HCT116, HepG2, Raji, HT-29, DLD-1, U2OS, PC-3, NCI-H1299, HEK293 as a 10x demo (HEK293 appears only inside VIPerturb-seq). | Empty as of this date. |
