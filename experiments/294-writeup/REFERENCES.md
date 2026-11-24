# REFERENCES.md — sources cited in `PAPER.md`

Base list from `experiments/287-novelty-audit/work/refs.md` (all web items accessed **2026-09-06** in a single session, 14:0x–15:5x Z), plus the project artefacts this write-up cites directly. Nothing new was fetched for experiment 294; the access dates below are those of the audit session that retrieved them.

**Round 6 (2026-09-07).** Nothing was fetched from the web for this revision either. The one new bibliographic entry, **[HORLBECK16]**, is assembled entirely from records already on disk: doi, locus, licence, file sizes and sha256 digests from `experiments/301-metadata-acquisitions/work/sources301.md`, and the author list and title from a cached PMC page in the same directory. Its row says so.

Quotations in `PAPER.md` are ≤ 15 words and attributed. One quotation the audit reached through an automated page fetch (a set of median cross-guide correlations attributed to [GWCD4i-P]) **failed to reproduce against the authors' own table** under five plausible filters and is **not used** anywhere in this manuscript.

---

## A. Source papers of the datasets the findings are measured on

| key | citation | route | accessed |
|---|---|---|---|
| **ORION** | *X-Atlas/Orion: Genome-wide Perturb-seq Datasets via a Scalable Fix-Cryopreserve Platform for Training Dose-Dependent Biological Foundation Models.* Xaira Therapeutics. bioRxiv 2025.06.11.659105 (posted 2025-06-16). | https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1.full | 2026-09-06 |
| **REPLOGLE22** | Replogle JM *et al.* *Mapping information-rich genotype–phenotype landscapes with genome-scale Perturb-seq.* Cell 185:2559–2575 (2022). doi:10.1016/j.cell.2022.05.013. | https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/ ; figshare+ 20029387 | 2026-09-06 |
| **HORLBECK16** | Horlbeck MA, Gilbert LA, Villalta JE, Adamson B, Pak RA, Chen Y, Fields AP, Park CY, Corn JE, Kampmann M, Weissman JS. *Compact and highly active next-generation libraries for CRISPR-mediated gene repression and activation.* eLife 5:e19760 (2016). doi:10.7554/eLife.19760. CC BY 4.0. The parent library of the challenge's 2025 constructs. **Supplementary file 3** ("Library composition of hCRISPRi-v2 and hCRISPRi-v2.1", 209,078 sgRNAs) supplies the protospacer pool rescanned for §3(g)'s off-target propensity and fills the 14 non-targeting pairs Replogle's table does not carry; **Supplementary file 2** ("TSS annotations for hg19 and mm10", 20,883 human rows, FANTOM-CAGE P1/P2/P1P2) is the annotation the library was designed against and the basis of the promoter windows. **Defect reported in §3(g):** the gene labels in these tables are Excel-escaped as `="SYMBOL"`, which silently stops an own-target exclusion from firing (per-spacer 12-nt hit rate 0.107 → 0.597 if uncorrected). | supp3 https://cdn.elifesciences.org/articles/19760/elife-19760-supp3-v2.xlsx (sha256 `afcc9b4c…`, 31,182,240 B); supp2 https://cdn.elifesciences.org/articles/19760/elife-19760-supp2-v2.xlsx (sha256 `2808fde8…`, 2,790,673 B) — both recorded with licence and digest in `experiments/301-metadata-acquisitions/work/sources301.md` | 2026-09-07 (files); the author list and title are read from the cached PMC page `E301/work/pmc9380471.html`, **not re-fetched from the web for this revision** |
| **GWCD4i-P** | *Genome-scale perturb-seq in primary human CD4⁺ T cells maps context-specific regulators of T cell programs and human immune traits.* bioRxiv 10.64898/2025.12.23.696273 (2025-12-23). | https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1.full | 2026-09-06 |
| **GWCD4i-C** | Same work, journal version: Cell, article S0092-8674(26)00929-3 (2026). Full text returned 403; the preprint and the authors' supplementary tables were used instead. | https://www.cell.com/cell/fulltext/S0092-8674(26)00929-3 (403) | 2026-09-06 |
| **GWCD4i-TAB** | The authors' data-sharing readme and supplementary tables: `DE_stats.suppl_table.csv` (33,983 rows; `donor_correlation_all_mean`, `donor_correlation_hits_mean`, `guide_correlation_all`, `guide_correlation_signif`), `K562_comparison.suppl_table.csv` (5,640 rows; `logfc_pearson_r`, `random_r1..3`, `donor_correlation_mean`, `n_degs_MASH_*`). | on disk: `data/external/cd4_cell2026/data_sharing_readme.md`; `data/external/cd4_cell2026/suppl_tables_2026-08-24/` | 2026-09-06 |
| **BROAD16** | Maffa S, Boyle IA, Ward L, Colgan WN, Borck P *et al.* *Dissecting context-dependent cancer vulnerabilities using Perturb-seq.* bioRxiv 10.64898/2026.08.24.746802 (2026-08-24). 16 cancer lines × 100 Cas9 KO targets, 6-day screens; Figshare article 33273600, CC BY 4.0. | https://www.biorxiv.org/content/10.64898/2026.08.24.746802v1.full | 2026-09-06 |
| **PAN26** | Pan X, Saunders R, Replogle J, Weissman J, Zhuang X. *Global cell-state and gene-program representations reveal conserved and context-specific perturbation responses of cells.* GEO GSE295214 (public 2026-05-05); bioRxiv 10.64898/2026.05.16.725005. Genome-scale CRISPRi Perturb-seq in hESCs. | on disk: `data/external/gse295214/GSE295214_family.soft.gz`; https://www.biorxiv.org/content/10.64898/2026.05.16.725005v1 | 2026-09-06 |
| **NADIG25** | Nadig A, Replogle JM, Pogson AN, Murthy M, McCarroll SA, Weissman JS, Robinson EB, O'Connor LJ. *Transcriptome-wide analysis of differential expression in perturbation atlases.* Nat Genet 57:1228–1237 (2025). doi:10.1038/s41588-025-02169-3. | https://pmc.ncbi.nlm.nih.gov/articles/PMC11244993/ | 2026-09-06 |
| **PEIDLI-NV** | Peidli S. *Quantifying reproducible multivariate differential expression in genetic perturbation screens.* Nat Genet 57:1066–1067 (2025) — News & Views on NADIG25. doi:10.1038/s41588-025-02177-3. | nature.com landing (paywalled body) | 2026-09-06 |
| **IPSC34** | *A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines.* Cell Genomics (2025), S2666-979X(25)00332-5. 7,226 genes × 34 iPSC lines from 26 donors. | cell.com landing (403 to automated fetch); summary via search index | 2026-09-06 |
| **SCPERTURB** | Peidli S, Green TD, Shen C *et al.* *scPerturb: harmonized single-cell perturbation data.* Nat Methods 21:531–540 (2024). doi:10.1038/s41592-023-02144-y. | nature.com landing; bioRxiv 2022.08.20.504663v3 | 2026-09-06 |
| **PISCES/X-CELL** | *X-Cell: Scaling Causal Perturbation Prediction Across Diverse Cellular Contexts via Diffusion Language Models* (X-Atlas/Pisces, 25.6 M cells, 16 contexts). bioRxiv 10.64898/2026.03.18.712807. Cited only as unreleased data. | on disk: `data/external/xcell/` | 2026-09-06 |

## B. Reliability, reproducibility and evaluation of perturbation data

| key | citation | route | accessed |
|---|---|---|---|
| **WANG26** | Wang X, Kuipers J, Hugi F, Platt RJ, Beerenwinkel N. *Reliable single-cell perturbations explain and improve model performance.* bioRxiv 10.64898/2026.08.11.744177 (2026-08-11). | https://www.biorxiv.org/content/10.64898/2026.08.11.744177v1.full | 2026-09-06 |
| **TXPERT** | *TxPert: using multiple knowledge graphs for prediction of transcriptomic perturbation effects.* Nat Biotechnol (2026), s41587-026-03113-4; preprint arXiv:2505.14919. Uses split-half experimental reproducibility as an evaluation ceiling. | arXiv:2505.14919 | 2026-09-06 |
| **GEOM26** | *Geometric coherence of single-cell CRISPR perturbations reveals regulatory architecture and predicts cellular stress.* arXiv:2604.16642. Split-half cosine of half-shift vectors; explicitly a sampling-precision, not a reproducibility, statistic; uses independent guides as the non-circular check. | https://arxiv.org/html/2604.16642v3 | 2026-09-06 |
| **PDSGEOM** | Liu Q, Zhang Q, Du J, Zhao S, Wang J. *Effects of Distance Metrics and Scaling on the Perturbation Discrimination Score.* arXiv:2511.16954 (2025-11-21). | https://arxiv.org/html/2511.16954v1 | 2026-09-06 |
| **CDS26** | *Score Distributions, Not Cells: Evaluating Single-Cell Perturbations Under Class Overlap.* arXiv:2607.04595. | https://arxiv.org/abs/2607.04595 | 2026-09-06 |
| **SYSTEMA** | Viñas Torné R, Wiatrak M, Piran Z, Fan S, Jiang L, Teichmann SA, Nitzan M, Brbić M. *Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation.* Nat Biotechnol (2025), s41587-025-02777-8; PMID 40854979. | nature.com landing | 2026-09-06 |
| **AHLMANN25** | Ahlmann-Eltze C, Huber W, Anders S. *Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines.* Nat Methods 22:1657–1661 (2025). doi:10.1038/s41592-025-02772-6. | nature.com; PMC12328236 | 2026-09-06 |
| **PERTURBENCH** | Wu Y, Wershof E, Schmon SM, Nassar M, Osiński B, Eksi R, Yan Z, Stark R, Zhang K, Graepel T. *PerturBench: Benchmarking Machine Learning Models for Cellular Perturbation Analysis.* arXiv:2408.10609 (2024). | https://arxiv.org/abs/2408.10609 | 2026-09-06 |
| **BENDIDI24** | Bendidi I *et al.* *Benchmarking Transcriptomics Foundation Models for Perturbation Analysis: one PCA still rules them all.* arXiv:2410.13956 (2024). | https://arxiv.org/abs/2410.13956 | 2026-09-06 |
| **CSENDES25** | Csendes G, Sanz G, Szalay KZ, Szalai B. *Benchmarking foundation cell models for post-perturbation RNA-seq prediction.* BMC Genomics 26:393 (2025). | link.springer.com/article/10.1186/s12864-025-11600-2 | 2026-09-06 |
| **CALIB25** | *Deep Learning-Based Genetic Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics.* bioRxiv 2025.10.20.683304. | https://www.biorxiv.org/content/10.1101/2025.10.20.683304v1.full | 2026-09-06 |
| **TRANSPORT26** | *Mechanisms Matter: Transportability of Cellular Perturbation Effects.* bioRxiv 10.64898/2026.05.08.723625 (2026-05-08). Frames the virtual-cell task as causal transportability (Bareinboim & Pearl 2013); when mechanisms differ across contexts the bias is irrecoverable at any sample size. | https://www.biorxiv.org/content/10.64898/2026.05.08.723625v1.full | 2026-09-06 |
| **PERTURBMAP** | Cui P, Liu Y, Sun W. *PerturbMap: Cross-Context Transfer of Single-Cell Perturbation Responses.* arXiv:2607.28090 (2026-07-30). Weights transported source responses by route reliability on validation anchors. | https://arxiv.org/abs/2607.28090 | 2026-09-06 |

## C. The challenge

| key | citation | route | accessed |
|---|---|---|---|
| **VCC-CELL** | *Virtual Cell Challenge: Toward a Turing test for the virtual cell.* Cell 188(13) (2025). | https://www.cell.com/cell/fulltext/S0092-8674(25)00675-0 | 2026-09-06 |
| **VCC25-WRAP** | Arc Institute. *Virtual Cell Challenge 2025 Wrap-Up: Winners and Reflections.* | https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up | 2026-09-06 |
| **VCC26-NEWS** | Arc Institute. *The 2026 Virtual Cell Challenge: predicting perturbation responses in cell contexts a model has never seen.* | https://arcinstitute.org/news/virtual-cell-challenge-2026 | 2026-09-06 |
| **VCC-BTD** | Arc Institute. *Behind the Data of the Virtual Cell Challenge* (panel construction). | https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge | 2026-09-06 |
| **VCC-SPEC** | `cell-eval2` `vcc2026` metric profile, vendored. Source of the split-half score-1 anchor quoted in §1.2. | on disk: `docs/reference/vcc2026-metrics-FULL.md`, `-BRIEF.md` | 2026-09-06 |
| **VCC-MANIFEST** | The 2026 challenge manifest and perturbation counts: 300 constructs, 400 cells per perturbation, six CRISPRi lines on 10x Flex, zero-shot. Basis for the 25/300 library-overlap fact in §1.3 and §3(g). | on disk: `data/challenge/manifest.json`, `data/challenge/pert_counts.csv` | 2026-09-06 |

## D. Estimator and metric-geometry precedents

| key | citation | route | accessed |
|---|---|---|---|
| **LEDOIT04** | Ledoit O, Wolf M. *A well-conditioned estimator for large-dimensional covariance matrices.* J Multivar Anal 88:365–411 (2004); and *Analytical nonlinear shrinkage of large-dimensional covariance matrices*, Ann Statist 48:3043–3065 (2020). | scikit-learn docs + search index | 2026-09-06 |
| **GAVISH14** | Gavish M, Donoho DL. *The optimal hard threshold for singular values is 4/√3.* IEEE Trans Inf Theory 60(8):5040–5053 (2014); arXiv:1305.5870. | https://arxiv.org/abs/1305.5870 | 2026-09-06 |
| **RMT-SC** | Aparicio L, Bordyuh M, Blumberg AJ, Rabadan R. *A Random Matrix Theory Approach to Denoise Single-Cell Data.* Patterns 1(3):100035 (2020). | https://www.cell.com/patterns/pdfExtended/S2666-3899(20)30040-4 | 2026-09-06 |
| **BIWHITEN** | Biwhitening (Sinkhorn–Knopp row/column scaling) so that a count matrix's spectrum follows Marchenko–Pastur (Landa, Kluger *et al.*). | via search index; see RMT-SC thread | 2026-09-06 |
| **CVPCA** | Stringer C, Pachitariu M, Steinmetz N, Carandini M, Harris KD. *High-dimensional geometry of population responses in visual cortex.* Nature 571:361–365 (2019). Cross-validated PCA: independent trial noise contributes zero in expectation. | via search index | 2026-09-06 |
| **CVPCA-BIAS** | *Revisiting the high-dimensional geometry of population responses in visual cortex.* bioRxiv 2024.02.16.580726; PNAS 2025, doi:10.1073/pnas.2506535122. cvPCA is a biased estimator of the signal eigenspectrum; proposes a moment-based alternative — the precedent for the cross-product moment identity used in §3(f). | https://www.pnas.org/doi/10.1073/pnas.2506535122 | 2026-09-06 |
| **HUBNESS** | Radovanović M, Nanopoulos A, Ivanović M. *Hubs in space: popular nearest neighbors in high-dimensional data.* JMLR 11:2487–2531 (2010). | via search index | 2026-09-06 |
| **ALLBUTTOP** | Mu J, Viswanath P. *All-but-the-Top: Simple and Effective Postprocessing for Word Representations.* ICLR 2018. | https://www.researchgate.net/publication/313398040 | 2026-09-06 |
| **WHITEN-CAVEAT** | Nielsen *et al.* (2023) on ZCA/marginal whitening with Mutual Proximity for hubness; Forooghi *et al.* (2024) on PCA/ZCA whitening degrading supervised classification, worse at higher dimension. The prescription tension T1 is stated against these. | via search index (OFAI TR-2014-01 and topic pages) | 2026-09-06 |

## E. Annotation and reference resources used as features or filters

| resource | use | on disk |
|---|---|---|
| DepMap (24Q2 segments; 26Q1) | gene essentiality and gene effect; arm-level copy number for 289's state features | `data/external/depmap26q1`, `data/external/ccle_broad_2025` |
| CCLE | expression breadth and level | `data/external/ccle_broad_2025` |
| HGNC gene families | the family axis and every family-held split | as used by `E281/work/feat281.csv`, `E284/work/splits284.npz` |
| MSigDB hallmark gene sets | the nine programme projections in 290; P53_PATHWAY / E2F_TARGETS / G2M_CHECKPOINT in 268 and 289 | as used by `E290/work/feat290.csv` |
| hu.MAP 2.0, Reactome complexes | complex membership features | `E281/work/humap2_complexes.txt`, `reactome_complexes_human.txt` |
| Ensembl BioMart | paralogue counts and orthologue features | `E281/work/bm_paralog.tsv`, `bm_orth_*.tsv` |
| LINCS / CMap L1000 | the negative control in §3(b): cross-line signature agreement does not predict transfer | via `E270`, read by `E276/work/analyse276.json` → `cmap` |
| hCRISPRi-v2 library composition + its CAGE TSS annotation ([HORLBECK16] supp. 3 and 2) | §3(g): the 209,070-protospacer rescan behind the off-target propensity, and the 18,387 promoter windows over 16,936 genes the seed rule scans | `E301/work/horlbeck2016_supp{2,3}_*.xlsx`, `promoters_cage_tss_hg19.fa.gz`; digests in `E301/work/sources301.md` |
| Replogle 2022 Table S1, tab `TabA_K562_day8_library` | §3(g): the 2025 challenge library's own protospacers (328 of 348 guide ids; all 34 constructs of the 17 two-construct targets resolve) | `E301/work/vcc2025_spacers.csv` |
| Ensembl GRCh37 REST (`sequence/region/human`) | §3(g): promoter sequence for the accepted hg19 window set; a first GRCh38 gene-level-TSS attempt was rejected by its own known-answer row | `E301/work/fetch_promoters_hg19.py` |

## F. Project artefacts cited as evidence

These are the objects `CLAIMS.md` points at; they are internal to this repository and are listed so every number in `PAPER.md` can be traced without re-deriving it.

| artefact | what it holds |
|---|---|
| `docs/FINDINGS-2026-09.md` (v3.0 at this revision; v1.9 when the base list was written) | F1–F16; the v1.x content had 447 numbers re-verified against their source JSON/CSV by experiment 275. **Every number this manuscript takes from F15/F16-era work was read from the experiment's own `RESULTS.md`, not from this summary** |
| `docs/CONFIRMATION-LEDGER.md` | claims C1–C9 with results, the reserve declarations (R1, R3, R4, R5 — all spent) and the once-per-claim access rule. C8 and C9 declared no reserve |
| `experiments/287-novelty-audit/NOVELTY-AUDIT.md` | the 24-row provenance table, the Tier 1 / Tier 2 split, tensions T1–T5, and the coordinator's 14:45Z addendum reconciling C-1 |
| `docs/RESEARCH-PROGRAMME-2026-09.md` | the organising question, rules 1–10, the workstreams and the Round 4 results |
| `E276/work/DEFINITION-276.md` | the conservation statistic, frozen before any cosine was computed |
| `E276/work/conservation276.csv`, `E279/work/znf279.csv` | the 300-target and 17,106-gene conservation maps |
| `E285/work/bench285.csv`, `analyse285.json`, `equiv285.json` | the matched reliability benchmark, the depth law and the unit decomposition |
| `E284/work/{fam284,rem284,fit284}.json`, `pred284.csv` | the strict-fold family axis and the remainder test |
| `E290/work/fit290.json`, `pred290.csv`; `E291/work/{surfA291,surfB291}.json` | response geometry inside and outside the map |
| `E268/work/decomp268_cpm5.json`, `lolo268_K2.csv`; `E283/work/c2_results.json`; `E289/work/{lolo289.csv,p1p2_289.json,p1_scores289.csv}` | the Broad-16 decomposition, its reserve confirmation and the state-feature rerun |
| `E272/work/analysis272.json`; `E286/work/analysisB286.json` | the operator spectrum, the damping arms and the cross-product estimator |
| `E262/work/{analysis262.json,tables262.md}`; `E273/work/tables273.md` | the deficit decomposition and the knee on both references |
| `E288/work/check288.json` | the magnitude decomposition of F8's association |
| `E292/work/reserve292.json`, `E292/RESERVE.md` | reserve R3's manifest and declaration (metadata only) |
| `E293/RESULTS.md`, `work/construct293.csv`, `strat293.json`, `ceiling293.json` | the construct ceiling in the challenge's own regime (2025 H1 library, Flex, 400 cells), its knockdown stratification and the zero-shot ceiling |
| `E295/RESULTS.md`, `work/reserve295.csv`, `stats295.json`, `paired295.json` | claim C5 on reserve R3: F8's score and F13's sign on eight untouched datasets |
| `E301/RESULTS.md`, `work/{vcc2025_spacers.csv,sources301.md}` | the 2025 library's protospacers and promoter windows, and the failing-versus-working seed scan (4/13 vs 4/21, Fisher p = 0.679) that §4.5's T7 says still stands |
| `E305/RESULTS.md` and `COORDINATOR-REVIEW-2026-09-07.md` | the per-pair expression test T7 superseded: zero of eleven pass Holm (p = 0.10476), and with 8 testable pairs the attainable minimum Holm-adjusted p was 0.0506 > 0.05 |
| `E307/RESULTS.md`, `work/{enrichment307.csv,results307.json,sens307.json,extra307.json}` | claim C8: the aggregate candidate-repression test with the construct as the unit, and the two known-answer rows it aborted on |
| `E309/RESULTS.md`, `work/{enrichment309.csv,results309.json,sens309.json,ka309.json}` | claim C9: the calibration repair, its nine known-answer rows, the partner-construct specificity test and the real-spacer resampling null |
| `E308/RESULTS.md`, `work/kd308.csv` | §3(h)'s closure of donor quality weighting and the measured reason it fails |
| `E306/RESULTS.md`, `work/{sources306.md,lb_live_306.json}` | §3(h)'s reconnaissance: the post-freeze data sweep (S1–S32, each with the sentence relied on) and the 769-entry leaderboard read of 2026-09-07 15:48:55Z |

---

## Added from the audit extension (297, 2026-09-06; routes and access dates in `experiments/297-audit-extension/work/refs297.md`)

# References consulted — 297 novelty-audit extension

Session **2026-09-06, 16:47–17:05Z**. All web items accessed **2026-09-06** unless stated. On-disk items are cited by
repo path. Every quotation in `AUDIT-EXTENSION.md` is ≤ 15 words and attributed; everything else is paraphrase.

## 0. Provenance discipline (inherited from 287, which caught one bad automated quote)

Three routes, marked per row below:

- **[F]** = the paper's **own text** was fetched and the quote read out of it this session (WebFetch on the full text /
  HTML / PMC page).
- **[S]** = **search-index summary only** — the paper's own text was NOT retrieved (PDF binary, HTTP 429, or paywall).
  Used for **direction and existence of a prior report**; any number taken from an [S] row is labelled "unverified
  against the source text" at its point of use in `AUDIT-EXTENSION.md`.
- **[287]** = inherited verbatim from `experiments/287-novelty-audit/work/refs.md`; not re-fetched today.

**Known-answer row for this audit (KA-297).** No new computation was run (reading-and-citation mandate). Instead the
two numbers this audit re-uses from 287's read-only checks were re-read out of
`experiments/287-novelty-audit/work/checks287.log` and had to match `NOVELTY-AUDIT.md` and
`docs/FINDINGS-2026-09.md` exactly:

| row | expected (as recorded by 287) | read today from `checks287.log` | verdict |
|---|---|---|---|
| KA-1 (GWCD4i K562 table, real vs random cross-cell-type r) | +0.3030 / +0.3147 / +0.3117 vs +0.0040 / +0.0037 / +0.0038 | identical, lines 1–5, "KA-1 PASS" | **PASS** |
| ρ(conservation, `n_degs_MASH`) in the same table | +0.51 … +0.54 (FINDINGS F8) | +0.507 (Rest) / +0.530 (Stim8hr) / +0.543 (Stim48hr), lines 35–37 | **PASS** |
| `guide_correlation_all` medians (GWCD4i DE_stats) | +0.062 / +0.072 / +0.061 (287 row 3) | +0.0617 / +0.0715 / +0.0607, lines 16–20 | **PASS** |

No `uv run python` was run; no file under `data/` or outside `experiments/297-audit-extension/` was written.

## A. New sources read this session

| key | citation | route | accessed |
|---|---|---|---|
| **WANG26** | Wang X, Kuipers J, Hugi F, Platt RJ, Beerenwinkel N. *Reliable single-cell perturbations explain and improve model performance.* bioRxiv 10.64898/2026.08.11.744177 (2026-08-11). Defines reliability ρ and shared-axis alignment φ; the 65 / 11 / 24 % split; Shared > Specific > Unreliable model performance. | **[F]** https://www.biorxiv.org/content/10.64898/2026.08.11.744177v1.full | 2026-09-06 |
| **GEOM26** | *Geometric coherence of single-cell CRISPR perturbations reveals regulatory architecture and predicts cellular stress.* arXiv:2604.16642 (v3). Coherence S_p = mean cosine of per-cell shift vectors to the mean direction; split-half reproducibility across six datasets. **Provenance note:** the arXiv landing/HTML gives this title; the PDF's embedded metadata reads *"Directional coherence and effect magnitude in single-cell CRISPR perturbation responses"* (Prashant C. Raju) — same arXiv id, title apparently changed between versions. | **[F]** https://arxiv.org/html/2604.16642v3 (HTML); https://arxiv.org/abs/2604.16642 (abstract) | 2026-09-06 |
| **REPLOGLE22** | Replogle JM *et al.* *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq.* Cell 185:2559–2575 (2022). doi:10.1016/j.cell.2022.05.013. Re-fetched today for the guide/knockdown numbers (median knockdown, dual-sgRNA element design, knockdown↔phenotype). | **[F]** https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/ | 2026-09-06 |
| **HARTMAN26** | Hartman A, Blair JD, Nguyen TP, … Satija R, Satpathy AT, Roth TL. *Systematic identification of seed-driven off-target effects in Perturb-seq experiments.* bioRxiv 10.64898/2026.03.27.714658v2, posted 2026-03-30. | **[F]** https://www.biorxiv.org/content/10.64898/2026.03.27.714658v2.full | 2026-09-06 |
| **NADIG25** | Nadig A, Replogle JM, … O'Connor LJ. *Transcriptome-wide analysis of differential expression in perturbation atlases.* Nat Genet 57:1228–1237 (2025). doi:10.1038/s41588-025-02169-3. Re-fetched for the cross-cell-type consistency numbers (TI correlations 0.40–0.75; 56 % / 44 % split) and π_DEG. | **[F]** https://pmc.ncbi.nlm.nih.gov/articles/PMC11244993/ (preprint version on PMC) | 2026-09-06 |
| **ZHANG25-CRISPRi** | Zhang H, Zhang P, Bindels E, Mulugeta E. *Insights from pooled CRISPRi single-cell screens in K562 cells reveal gene functions, regulatory networks, and highlight opportunities and limitations.* bioRxiv 2025.10.23.684112 (2025-10-23). Reports that a large fraction of sgRNAs in published CRISPRi screens fail an efficacy criterion (≈ 97/240 iPSC, ≈ 56/377 K562). **PDF fetch returned binary; `.full` returned HTTP 429 twice — numbers are search-index-derived and UNVERIFIED against the source text.** | **[S]** https://www.biorxiv.org/content/10.1101/2025.10.23.684112v1.full | 2026-09-06 |
| **MIXSCALE** | *Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens.* Nat Cell Biol (2025), s41556-025-01622-z; preprint bioRxiv 2024.01.29.576933v2. Continuous per-cell perturbation score replacing mixscape's binary call; weights cells by knockdown strength. | **[S]** https://www.nature.com/articles/s41556-025-01622-z | 2026-09-06 |
| **HETERO25** | *Decoding heterogeneous single-cell perturbation responses.* Nat Cell Biol (2025), s41556-025-01626-9. Review; states that mixscape-family methods were designed for Cas9 KO and do not fully account for partial CRISPRi perturbation. | **[S]** https://www.nature.com/articles/s41556-025-01626-9 | 2026-09-06 |
| **SEED24** | *Seed sequences mediate off-target activity in the CRISPR-interference system.* Cell Genomics (2024), S2666-979X(24)00322-7; preprint bioRxiv 2024.04.10.588881. Off-target CRISPRi activity is pervasive and mostly explained by PAM-proximal seed complementarity. | **[S]** https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00322-7 | 2026-09-06 |
| **DUALSG** | *Maximizing CRISPRi efficacy and accessibility with dual-sgRNA libraries and optimal effectors.* eLife 81856. Dual-sgRNA library median knockdown 82 % vs 65 % for single-sgRNA Dolcetto in K562. | **[S]** https://elifesciences.org/articles/81856 | 2026-09-06 |
| **DEPTH-TSE** | Zhang MJ, Ntranos V, Tse D. *Quantifying the tradeoff between sequencing depth and cell number in single-cell RNA-seq* (bioRxiv 762773) and *One read per cell per gene is optimal for single-cell RNA-Seq* (bioRxiv 389296). Optimal budget allocation maximises cells at ≈ 1 read per cell per gene. | **[S]** https://www.biorxiv.org/content/10.1101/762773v1 ; https://www.biorxiv.org/content/10.1101/389296v2.full | 2026-09-06 |
| **SVENSSON-DEPTH** | Svensson V, da Veiga Beltrame E, Pachter L. Depth-vs-cells analysis in a VAE framework; benefit of added depth minor above ≈ 15,000 reads per cell. | **[S]** via search index | 2026-09-06 |
| **SCPOWER** | Schmid KT *et al.* *scPower accelerates and optimizes the design of multi-sample single cell transcriptomic studies.* Nat Commun 12:6625 (2021). doi:10.1038/s41467-021-26779-7. Pseudobulk DE power as a function of samples × cells × depth. | **[S]** https://www.nature.com/articles/s41467-021-26779-7 | 2026-09-06 |
| **JOLLIFFE82** | Jolliffe IT. *A note on the use of principal components in regression.* J R Stat Soc Ser C 31(3):300–303 (1982). Low-variance principal components can be the predictive ones; four published examples. | **[S]** https://academic.oup.com/jrsssc/article/31/3/300/6985100 | 2026-09-06 |
| **WHITEN-TASK** | Task-dependence of whitening: gains on STS/retrieval, losses on classification and at higher embedding dimension (Su et al. 2021 whitening-for-retrieval; Forooghi et al. 2024 degradation); the stated failure mechanism is that uniform rescaling promotes low-variance directions where covariance estimation error is worst; truncated / shrinkage whitening is the standard hedge. Extends 287's **WHITEN-CAVEAT**. | **[S]** https://www.emergentmind.com/topics/pca-whitening ; https://www.emergentmind.com/topics/whitening-sentence-representations | 2026-09-06 |
| **TXPERT-CEIL** | TxPert (Nat Biotechnol 2026, s41587-026-03113-4; arXiv:2505.14919): split-half experimental reproducibility as a benchmark ceiling, explicitly **not** a theoretical upper bound; a sample-based extrapolation to full-size reproducibility is also given. | **[S]** https://www.nature.com/articles/s41587-026-03113-4 | 2026-09-06 |
| **STRAIGHT26** | *Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward.* bioRxiv 10.64898/2026.02.14.705879. Reference/evaluation split as an "ideal model" proxy; bootstrapped positive controls as an upper bound. | **[S]** https://www.biorxiv.org/content/10.64898/2026.02.14.705879v1.full | 2026-09-06 |
| **URBUT19** | Urbut SM, Wang G, Carbonetto P, Stephens M. *Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions.* Nat Genet 51:187–195 (2019). doi:10.1038/s41588-018-0268-8. mash; cis-eQTL sharing across 44 tissues, 16,069 genes; sign largely shared, effect size varies. | **[S]** https://www.nature.com/articles/s41588-018-0268-8 | 2026-09-06 |
| **DEPPRIOR26** | Zhou F, Wang C, Mo L, Sun X, Zhang J. *DepPrior: integrating CRISPR dependency predictability with multi-omics reproducibility…* Front Pharmacol (2026). doi:10.3389/fphar.2026.1859926. Used only as the route to the DepMap portal's **Predictability** framing (selective dependencies with no interpretable predictive marker). | **[S]** https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2026.1859926/full | 2026-09-06 |

## B. Sources inherited from 287 and relied on here without re-fetching

`ORION`, `REPLOGLE22` (also re-fetched, above), `GWCD4i-P`, `GWCD4i-TAB`, `BROAD16`, `PAN26`, `PEIDLI-NV`, `IPSC34`,
`SCPERTURB`, `TRANSPORT26`, `PERTURBMAP`, `PDSGEOM`, `CDS26`, `SYSTEMA`, `PERTURBENCH`, `VCC-SPEC`, `VCC25-WRAP`,
`VCC-BTD`, `LEDOIT04`, `CVPCA` (Stringer 2019), `CVPCA-BIAS` (PNAS 2025, doi:10.1073/pnas.2506535122), `GAVISH14`,
`RMT-SC`, `HUBNESS`, `ALLBUTTOP`, `WHITEN-CAVEAT`. Full citations and routes: `experiments/287-novelty-audit/work/refs.md`.
**[287]** in the audit table means the quote is theirs, checked against their file, not re-fetched today.

## C. On-disk programme documents read (no writes)

`docs/READ-FIRST.md` (v5.4), `docs/FINDINGS-2026-09.md` (v2.0, F8–F14 and §3–§5),
`experiments/287-novelty-audit/NOVELTY-AUDIT.md`, `experiments/287-novelty-audit/work/refs.md`,
`experiments/287-novelty-audit/work/checks287.log`, `experiments/297-audit-extension/MANDATE.md`.
Experiment RESULTS.md files for 284 / 286 / 290 / 291 / 293 / 295 were read through their consolidated,
already-amended statements in `docs/FINDINGS-2026-09.md` v2.0 (F8, F12/284, F13/290+291+C5b, F14/285+293) and
`docs/READ-FIRST.md` §5–§6, which is the wording the mandate says to audit. Sibling 296's directory was not read.

## D. Searches run (queries kept for reproducibility)

Perturb-seq shared/generic response conserved across cell types vs magnitude · Replogle 2022 K562/RPE1 cophenetic ·
Nadig 2025 TRADE polygenicity · Mixscale CRISPRi knockdown heterogeneity · sgRNA seed off-target transcriptional
signature · DepMap selective dependency predictability · PerturbMap cross-context transfer · sequencing depth vs cell
number pseudobulk power · low-variance principal components predictive (Jolliffe) · whitening amplifies noise subspace
retrieval · reproducibility ceiling perturbation benchmark · Wang/Kuipers/Beerenwinkel shared-specific-unreliable ·
"Insights from pooled CRISPRi single-cell screens in K562" · scPerturb E-distance · Urbut/Stephens mash.


---

## Added from the F16 novelty audit (320, 2026-09-08)

Routes, provenance codes [F]/[P]/[S]/[X], access dates and the exact sentence relied on for every row below are in
`experiments/320-f16-novelty/work/refs320.md`. Nothing was fetched from the web for this revision of
`REFERENCES.md`; each row is transcribed from that file.

| key | citation | route | accessed |
|---|---|---|---|
| **ROHATGI24** | Rohatgi N, Fortin J-P, Lau T, Ying Y, Zhang Y, Lee BL, Costa MR, Reja R. *Seed sequences mediate off-target activity in the CRISPR-interference system.* Cell Genomics 4(11):100693 (2024). doi:10.1016/j.xgen.2024.100693. PMID 39510079. **The full citation for what earlier revisions cited as [SEED24].** Screens + RNA-seq + **dCas9 ChIP-seq occupancy** (four guides, 98–1,141 reproducible peaks each, 10–33 % promoter-proximal); the 10–12 nt PAM-proximal seed rule and a ±2,000 bp promoter window; 21 of ≈ 562 seed-matching sgRNAs at the GSDMD promoter behaved like on-target guides; run in the **Horlbeck/Weissman genome-wide CRISPRi library**, the parent of the challenge's 2025 constructs. Establishes the phenomenon and its mechanism at the level of binding, two years before §3(g). | **[F]** on the preprint https://www.biorxiv.org/content/10.1101/2024.04.10.588881v1.full (all quotes are from the preprint); journal version **[X]** HTTP 403 at sciencedirect.com and cell.com; bibliographic record [F] via EuropePMC | 2026-09-08 |
| **JACKSON03** | Jackson AL, Bartz SR, Schelter J, Kobayashi SV, Burchard J, Mao M, Li B, Cavet G, Linsley PS. *Expression profiling reveals off-target gene regulation by RNAi.* Nature Biotechnology 21(6):635–637 (2003). doi:10.1038/nbt831. PMID 12754523. The founding observation of the seed off-target field: different siRNAs against one transcript give siRNA-specific rather than target-specific signatures, and multiple duplexes per target are prescribed as the specificity control. **The source of the partner-construct control used in §3(g), which is therefore not ours.** | **[F]** read in full | 2026-09-08 |
| **SEEDMATCHR24** | Cazares T, Higgs RE, Wang J, Ozer HG. *SeedMatchR: identify off-target effects mediated by siRNA seed regions in RNA-seq experiments.* Bioinformatics 40(1):btae011 (2024). doi:10.1093/bioinformatics/btae011. The published aggregate seed-shift test — ECDF of seed-matching transcripts against a background set, one-sided KS — with the **gene** as the unit and no covariate matching. One of the two components §3(g)'s design recombines. | **[F]** https://pmc.ncbi.nlm.nih.gov/articles/PMC10799297/ | 2026-09-08 |
| **MATCHRANGES23** | Davis ES, Mu W, Lee S, Dozmorov MG, Love MI, Phanstiel DH. *matchRanges: generating null hypothesis genomic ranges via covariate-matched sampling.* Bioinformatics 39(5):btad197 (2023). doi:10.1093/bioinformatics/btad197. PMID 37084270. Covariate-matched null sets as a published method (propensity score; nearest-neighbour, rejection and stratified sampling), generating **one** matched null set rather than a repeated-draw empirical p. The second component §3(g)'s design recombines. | **[F]** preprint https://www.biorxiv.org/content/10.1101/2022.08.05.502985v1.full; bibliographic record [F] via EuropePMC | 2026-09-08 |
| **GOEMAN07** | Goeman JJ, Bühlmann P. *Analyzing gene expression data in terms of gene sets: methodological issues.* Bioinformatics 23(8):980–987 (2007). Competitive versus self-contained nulls; **gene sampling versus subject sampling**; gene-sampling nulls are anti-conservative under correlation between the resampled items (simulated rejection rates up to 50,000× the nominal level). The source for both the aggregate-over-per-item logic §3(g) uses and the correction that its headline p is a gene-sampling p. | **[F]** https://academic.oup.com/bioinformatics/article/23/8/980/198511 | 2026-09-08 |
| **HARTMAN26** *(row already present in §A above; the audit adds these readings)* | Systematic seed off-target workflow across three genome-wide Perturb-seq datasets; the 10–12 PAM-proximal seed rule with a **≥12 bp** high-confidence threshold and a ±2,000 bp TSS window; **47.2–59.8 %** of guides carry at least one ≥12 bp promoter match against **0.83 %** (0.2–1.8 %) with a nominated event; the partner-guide specificity control used on CRISPRi Perturb-seq for exactly §3(g)'s purpose; unit = the guide–gene pair, no covariate-matched resampling, no global aggregate statistic. **Preprint; the author list is [S], the title, venue and all quoted sentences [F].** | **[F]** https://www.biorxiv.org/content/10.64898/2026.03.27.714658v2.full | 2026-09-06 (re-read 2026-09-08) |

**One withdrawn row.** **[SEED24]** carried provenance **[S]** (search-index only), no author list and no DOI. It is
the same paper as **[ROHATGI24]**, which 320 read in full on the preprint; `PAPER.md` and `CLAIMS.md` now cite
`ROHATGI24` and `SEED24` is retained here only so that earlier revisions remain resolvable.

**Added 2026-09-25 for §3(j) (F18); provenance [S] until read in full:**

| key | reference |
|---|---|
| **SYSTEMA25** | Viñas Torné R, et al. *Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation.* Nat Biotechnol (2025), doi:10.1038/s41587-025-02777-8; PMID 40854979. Names "systematic variation" (consistent perturbed-vs-control differences) across ten datasets and five lines, and treats it as an evaluation confounder. §3(j) measures its line-specificity and its transfer across labs within a line. |
| **MIXSEQ20** | McFarland JM, et al. *Multiplexed single-cell transcriptional response profiling to define cancer vulnerabilities and therapeutic mechanism of action.* Nat Commun 11:4296 (2020). Line-specific drug responses with a shared component across pooled cancer lines. |
| **PAN25** | Pan X, Saunders R, Replogle JM, Zhuang X, Weissman JS, et al. GEO GSE295214 (2025): a genome-scale CRISPRi Perturb-seq screen in H1 hESCs (10x 5′), the same-line bank of §3(j). |
| **VIPERTURB26** | Bradu A, …, Satija R. *Genome-wide single-cell perturbation screens with VIPerturb-seq.* bioRxiv 2026.02.12.705613; Zenodo 18460279 (the K562 truth of §3(j)'s replication). |

## Added for Draft v15 (533, 2026-09-30)

Nothing was fetched from the web for this revision. Each row is transcribed from the experiment that opened the source; the
route, access date and the sentences relied on are in that experiment's `RESULTS.md`.

| key | reference | read by |
|---|---|---|
| **VCC26-CELL** | *Virtual Cell Challenge 2026: Benchmarking zero-shot generalization across cellular contexts.* Cell (2026), doi:10.1016/j.cell.2026.08.004, PMID 42648290, published online 2026-08-26, CC BY-NC-ND 4.0. Six CRISPRi Perturb-seq lines from different tissues on 10x Flex; three validation, three test; no line names. | 520 (full text via the Cell HTML and the Europe PMC abstract, 2026-09-28) |
| **JIANG25** | Jiang *et al.* (Satija laboratory, NYGC; the author list was not transcribed by 520). *Systematic reconstruction of molecular pathway signatures using scalable single-cell perturbation screens* (the GEO series title). bioRxiv 10.1101/2024.01.29.576933 (Nature Cell Biology version the same study); GEO GSE281048, Zenodo 10.5281/zenodo.14518762. Six cancer lines (A549, MCF7, HT29, HAP1, BxPC3, K562), CRISPRi (dCas9-KRAB-MeCP2), 44–61 regulators per pathway. The "Jiang 2025" banks of §3(j) and F-10; the closest published analogue of the challenge's design. | 520 (full text and supplement, 2026-09-28) |
| **ATLASSHIFT** | *AtlasShift* emitter by Kaipeng Zheng (board author, as 525 records), github.com/kaipengm2/Virtual-Cell-Challenge-2026 at commit `d24ce4fdae0cd7cba3ba29546cd8737094eae98a` (2026-09-23); `predict.py`, `model.py` (dual-moment count fitting). The most reproduced public emission design on the 2026 board (§3(h)). | 525 (code read at the pinned commit, blob sizes matched, 2026-09-29) |
| **DEPMAP-CHRONOS** | DepMap 24Q4 CRISPR gene effect (Chronos), `CRISPRGeneEffect` subset for K562 (ACH-000551) and HCT116 (ACH-000971), CC BY 4.0. §3(h)'s dependency test. | on disk: `experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv`; used by 514 |

## Added for Draft v16 (592, 2026-10-01)

| key | reference |
|---|---|
| **CELLOSAURUS** | Bairoch A. *The Cellosaurus, a cell-line knowledge resource.* J Biomol Tech 29:25–38 (2018). Records read via api.cellosaurus.org on 2026-10-01 [F]: CVCL_1640 PC-14 ("Problematic cell line: Contaminated. Shown to be a PC-9 derivative"), CVCL_2090 KCI-MOH1 (an HPAC derivative), CVCL_2094 KPL-1 (an MCF-7 derivative), CVCL_2450 GP2d (same individual as CVCL_1235 GP5d). §5.9. |
| **CCLE-2025** | Cancer Cell Line Encyclopedia, Broad 2025 release as mirrored by cBioPortal (`ccle_broad_2025`, DepMap Public 25Q3): `data_mrna_seq_tpm.txt` (1,600 lines with expression) and `data_clinical_sample.txt`; on disk under `data/external/ccle_broad_2025` (cache sha in `ccle_cache.meta.json`). §5.9's reference panel. |
