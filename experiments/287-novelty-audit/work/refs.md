# References consulted — 287 novelty audit

All web items accessed **2026-09-06** (single session, 14:0x–15:5x Z). On-disk items are cited by repo path.
Quotations elsewhere in this experiment are ≤ 15 words and attributed; everything else is paraphrase.

## A. The datasets our findings are measured on (source papers)

| key | citation | route | accessed |
|---|---|---|---|
| **ORION** | *X-Atlas/Orion: Genome-wide Perturb-seq Datasets via a Scalable Fix-Cryopreserve Platform for Training Dose-Dependent Biological Foundation Models.* Xaira Therapeutics. bioRxiv 2025.06.11.659105 (posted 2025-06-16). | https://www.biorxiv.org/content/10.1101/2025.06.11.659105v1.full | 2026-09-06 |
| **PISCES/X-CELL** | *X-Cell: Scaling Causal Perturbation Prediction Across Diverse Cellular Contexts via Diffusion Language Models* (X-Atlas/Pisces, 25.6 M cells, 16 contexts). bioRxiv 10.64898/2026.03.18.712807. | on disk: `data/external/xcell/xcell_biorxiv_full.txt`, `xcell_preprint.txt`, `xcell_README.md`, `pisces_README.md` | 2026-09-06 |
| **REPLOGLE22** | Replogle JM *et al.* *Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq.* Cell 185:2559–2575 (2022). doi:10.1016/j.cell.2022.05.013. | https://pmc.ncbi.nlm.nih.gov/articles/PMC9380471/ ; figshare+ 20029387 | 2026-09-06 |
| **GWCD4i-P** | *Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators of T cell programs and human immune traits.* bioRxiv 10.64898/2025.12.23.696273 (2025-12-23). | https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1.full | 2026-09-06 |
| **GWCD4i-C** | same work, journal version: Cell, article S0092-8674(26)00929-3 (2026). **Full text 403 to WebFetch**; used the preprint + the authors' own supplementary tables instead. | https://www.cell.com/cell/fulltext/S0092-8674(26)00929-3 (403) | 2026-09-06 |
| **GWCD4i-TAB** | The authors' data-sharing readme and supplementary tables: `DE_stats.suppl_table.csv` (33,983 rows; columns `donor_correlation_all_mean`, `donor_correlation_hits_mean`, `guide_correlation_all`, `guide_correlation_signif`), `K562_comparison.suppl_table.csv` (5,640 rows; `logfc_pearson_r`, `random_r1..3`, `donor_correlation_mean`, `n_degs_MASH_*`). | on disk: `data/external/cd4_cell2026/data_sharing_readme.md`, `data/external/cd4_cell2026/suppl_tables_2026-08-24/` | 2026-09-06 |
| **BROAD16** | Maffa S, Boyle IA, Ward L, Colgan WN, Borck P *et al.* *Dissecting context-dependent cancer vulnerabilities using Perturb-seq.* bioRxiv 10.64898/2026.08.24.746802 (2026-08-24). 16 cancer lines × 100 Cas9 KO targets, 6-day screens; Figshare article 33273600, CC BY 4.0. | https://www.biorxiv.org/content/10.64898/2026.08.24.746802v1.full | 2026-09-06 |
| **PAN26** | Pan X, Saunders R, Replogle J, Weissman J, Zhuang X. *Global cell-state and gene-program representations reveal conserved and context-specific perturbation responses of cells.* GEO **GSE295214** (public 2026-05-05), bioRxiv 10.64898/2026.05.16.725005. Genome-scale CRISPRi Perturb-seq in hESCs. | on disk: `data/external/gse295214/GSE295214_family.soft.gz`; https://www.biorxiv.org/content/10.64898/2026.05.16.725005v1 | 2026-09-06 |
| **NADIG25** | Nadig A, Replogle JM, Pogson AN, Murthy M, McCarroll SA, Weissman JS, Robinson EB, O'Connor LJ. *Transcriptome-wide analysis of differential expression in perturbation atlases.* Nat Genet 57:1228–1237 (2025). doi:10.1038/s41588-025-02169-3; preprint bioRxiv 2024.07.03.601903. | https://pmc.ncbi.nlm.nih.gov/articles/PMC11244993/ ; nature.com landing | 2026-09-06 |
| **PEIDLI-NV** | Peidli S. *Quantifying reproducible multivariate differential expression in genetic perturbation screens.* Nat Genet 57:1066–1067 (2025) — News & Views on NADIG25. doi:10.1038/s41588-025-02177-3. | nature.com landing (paywalled body) | 2026-09-06 |
| **IPSC34** | *A genome-scale single-cell CRISPRi map of trans gene regulation across human pluripotent stem cell lines.* Cell Genomics (2025), S2666-979X(25)00332-5. 7,226 genes × 34 iPSC lines from 26 donors; targeted arm 444 genes × 19–20 lines from 10 donors, 1,161,864 cells, median 74 cells/target. | cell.com landing (403 to WebFetch); summary via search index | 2026-09-06 |
| **SCPERTURB** | Peidli S, Green TD, Shen C *et al.* *scPerturb: harmonized single-cell perturbation data.* Nat Methods 21:531–540 (2024). doi:10.1038/s41592-023-02144-y. E-distance / E-test. | nature.com landing; bioRxiv 2022.08.20.504663v3 | 2026-09-06 |

## B. Reliability, reproducibility and evaluation of perturbation data

| key | citation | route | accessed |
|---|---|---|---|
| **WANG26** | Wang X, Kuipers J, Hugi F, Platt RJ, Beerenwinkel N. *Reliable single-cell perturbations explain and improve model performance.* bioRxiv 10.64898/2026.08.11.744177 (2026-08-11). | https://www.biorxiv.org/content/10.64898/2026.08.11.744177v1.full | 2026-09-06 |
| **TXPERT** | *TxPert: using multiple knowledge graphs for prediction of transcriptomic perturbation effects.* Nat Biotechnol (2026), s41587-026-03113-4; preprint arXiv:2505.14919. Uses split-half experimental reproducibility as an evaluation ceiling. | nature.com (redirect wall); arXiv:2505.14919 | 2026-09-06 |
| **GEOM26** | *Geometric coherence of single-cell CRISPR perturbations reveals regulatory architecture and predicts cellular stress.* arXiv:2604.16642. Split-half cosine of half-shift vectors, 50 splits; Replogle n=1,832 and Norman n=236; explicitly a sampling-precision, not a reproducibility, statistic. | https://arxiv.org/html/2604.16642v3 | 2026-09-06 |
| **PDSGEOM** | Liu Q, Zhang Q, Du J, Zhao S, Wang J. *Effects of Distance Metrics and Scaling on the Perturbation Discrimination Score.* arXiv:2511.16954 (2025-11-21). | https://arxiv.org/html/2511.16954v1 | 2026-09-06 |
| **CDS26** | *Score Distributions, Not Cells: Evaluating Single-Cell Perturbations Under Class Overlap.* arXiv:2607.04595. Reports PDS variants degrading as cells become scarce. | https://arxiv.org/abs/2607.04595 | 2026-09-06 |
| **SYSTEMA** | Viñas Torné R, Wiatrak M, Piran Z, Fan S, Jiang L, Teichmann SA, Nitzan M, Brbić M. *Systema: a framework for evaluating genetic perturbation response prediction beyond systematic variation.* Nat Biotechnol (2025), s41587-025-02777-8; PMID 40854979. | nature.com landing | 2026-09-06 |
| **AHLMANN25** | Ahlmann-Eltze C, Huber W, Anders S. *Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines.* Nat Methods 22:1657–1661 (2025). doi:10.1038/s41592-025-02772-6; preprint bioRxiv 2024.09.16.613342. | nature.com; PMC12328236 | 2026-09-06 |
| **PERTURBENCH** | Wu Y, Wershof E, Schmon SM, Nassar M, Osiński B, Eksi R, Yan Z, Stark R, Zhang K, Graepel T. *PerturBench: Benchmarking Machine Learning Models for Cellular Perturbation Analysis.* arXiv:2408.10609 (2024). Rank metrics; mode collapse. | https://arxiv.org/abs/2408.10609 | 2026-09-06 |
| **BENDIDI24** | Bendidi I, Whitfield S, Kenyon-Dean K, Ben Yedder H, El Mesbahi Y, Noutahi E, Denton AK. *Benchmarking Transcriptomics Foundation Models for Perturbation Analysis: one PCA still rules them all.* arXiv:2410.13956 (2024). | https://arxiv.org/abs/2410.13956 | 2026-09-06 |
| **CSENDES25** | Csendes G, Sanz G, Szalay KZ, Szalai B. *Benchmarking foundation cell models for post-perturbation RNA-seq prediction.* BMC Genomics 26:393 (2025). | link.springer.com/article/10.1186/s12864-025-11600-2 | 2026-09-06 |
| **CALIB25** | *Deep Learning-Based Genetic Perturbation Models Do Outperform Uninformative Baselines on Well-Calibrated Metrics.* bioRxiv 2025.10.20.683304. Rebuttal; metric calibration with positive/negative controls. | https://www.biorxiv.org/content/10.1101/2025.10.20.683304v1.full | 2026-09-06 |
| **TRANSPORT26** | *Mechanisms Matter: Transportability of Cellular Perturbation Effects.* bioRxiv 10.64898/2026.05.08.723625 (2026-05-08). Frames the virtual-cell task as causal transportability (Bareinboim & Pearl 2013); when mechanisms differ across contexts the bias is irrecoverable at any sample size. | https://www.biorxiv.org/content/10.64898/2026.05.08.723625v1.full | 2026-09-06 |
| **PERTURBMAP** | Cui P, Liu Y, Sun W. *PerturbMap: Cross-Context Transfer of Single-Cell Perturbation Responses.* arXiv:2607.28090 (2026-07-30). | https://arxiv.org/abs/2607.28090 | 2026-09-06 |

## C. The challenge itself

| key | citation | route | accessed |
|---|---|---|---|
| **VCC-CELL** | *Virtual Cell Challenge: Toward a Turing test for the virtual cell.* Cell 188(13) (2025). | https://www.cell.com/cell/fulltext/S0092-8674(25)00675-0 | 2026-09-06 |
| **VCC25-WRAP** | Arc Institute. *Virtual Cell Challenge 2025 Wrap-Up: Winners and Reflections.* | https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up | 2026-09-06 |
| **VCC26-NEWS** | Arc Institute. *The 2026 Virtual Cell Challenge: predicting perturbation responses in cell contexts a model has never seen.* | https://arcinstitute.org/news/virtual-cell-challenge-2026 | 2026-09-06 |
| **VCC-BTD** | Arc Institute. *Behind the Data of the Virtual Cell Challenge.* (panel construction; quoted in `AGENTS.md` §2.1) | https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge | 2026-09-06 |
| **VCC-SPEC** | `cell-eval2` `vcc2026` profile, vendored: `docs/reference/vcc2026-metrics-FULL.md`, `-BRIEF.md`. | on disk | 2026-09-06 |

## D. Estimator precedents (for 286's cross-product estimator and F3's spectrum)

| key | citation | route | accessed |
|---|---|---|---|
| **CVPCA** | Stringer C, Pachitariu M, Steinmetz N, Carandini M, Harris KD. *High-dimensional geometry of population responses in visual cortex.* Nature 571:361–365 (2019). Introduces cross-validated PCA (cvPCA): project one repeat onto PCs of the other; independent trial noise contributes zero in expectation. | via search index; follow-up below | 2026-09-06 |
| **CVPCA-BIAS** | *Revisiting the high-dimensional geometry of population responses in visual cortex.* bioRxiv 2024.02.16.580726; published PNAS 2025 (doi:10.1073/pnas.2506535122). **cvPCA is a biased estimator of the signal eigenspectrum**; proposes a moment-based (eigenmoment) alternative; the V1 spectrum is a broken power law, not a single power law. | https://www.biorxiv.org/content/10.1101/2024.02.16.580726v1.full ; https://www.pnas.org/doi/10.1073/pnas.2506535122 | 2026-09-06 |
| **GAVISH14** | Gavish M, Donoho DL. *The optimal hard threshold for singular values is 4/√3.* IEEE Trans Inf Theory 60(8):5040–5053 (2014); arXiv:1305.5870. | https://arxiv.org/abs/1305.5870 | 2026-09-06 |
| **RMT-SC** | Aparicio L, Bordyuh M, Blumberg AJ, Rabadan R. *A Random Matrix Theory Approach to Denoise Single-Cell Data.* Patterns 1(3):100035 (2020). MP bulk vs Tracy–Widom edge; sparsity breaks universality for ≈5 % of eigenvalues. | https://www.cell.com/patterns/pdfExtended/S2666-3899(20)30040-4 | 2026-09-06 |
| **BIWHITEN** | Biwhitening (Sinkhorn–Knopp row/column scaling) so that the count-matrix spectrum follows Marchenko–Pastur and outliers are those above λ₊ = (1+√q)². (Landa, Kluger *et al.*) | via search index; see RMT-SC thread | 2026-09-06 |
| **LEDOIT04** | Ledoit O, Wolf M. *A well-conditioned estimator for large-dimensional covariance matrices.* J Multivar Anal 88:365–411 (2004); and *Analytical nonlinear shrinkage…*, Ann Statist 48:3043–3065 (2020). Sample eigenvalues are biased away from 1; the bias shrinks in n/p. | scikit-learn docs + search index | 2026-09-06 |
| **HUBNESS** | Radovanović M, Nanopoulos A, Ivanović M. *Hubs in space: popular nearest neighbors in high-dimensional data.* JMLR 11:2487–2531 (2010). | via search index | 2026-09-06 |
| **ALLBUTTOP** | Mu J, Viswanath P. *All-but-the-Top: Simple and Effective Postprocessing for Word Representations.* ICLR 2018. Mean removal + top-PC removal ≈ isotropy; improves similarity/retrieval. | https://www.researchgate.net/publication/313398040 | 2026-09-06 |
| **WHITEN-CAVEAT** | Nielsen *et al.* (2023) on ZCA/marginal whitening + Mutual Proximity for hubness (skewness −69…−83 %); Forooghi *et al.* (2024) on PCA/ZCA whitening degrading supervised classification, worse at higher dimension. | via search index (emergentmind topic pages, OFAI TR-2014-01) | 2026-09-06 |

## E. Small local checks run for this audit (`uv run python`, read-only)

Script `work/checks287.py`, output `work/checks287.log`, seed 287. Known-answer row KA-1 (random controls at chance,
real cross-cell-type correlation ≈ 0.30) **PASS** before any other number is read.

1. `data/external/cd4_cell2026/suppl_tables_2026-08-24/DE_stats.suppl_table.csv` — distributions of the authors' own
   `donor_correlation_all_mean`, `donor_correlation_hits_mean`, `guide_correlation_all`, `guide_correlation_signif`,
   whole screen and restricted to the 300 challenge panel targets, by culture condition.
2. Same directory, `K562_comparison.suppl_table.csv` — the authors' own K562-vs-CD4 cross-cell-type
   `logfc_pearson_r` against its three random controls; the ZNF-vs-rest contrast inside that table
   (Mann–Whitney + 4,000-draw bootstrap, seed 287); Spearman of cross-line conservation against
   `donor_correlation_mean` and against `n_degs_MASH_*`.
3. `data/external/gse295214/GSE295214_family.soft.gz` — series title, contributors and design for PAN26.
4. `experiments/268-broad16-context-residuals/work/inventory.json` — the 16 line names, for BROAD16 identification.

No dataset under `data/` was modified; nothing was written outside `experiments/287-novelty-audit/`.
