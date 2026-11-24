# 532 — Literature and data refresh, 15–29 September 2026

**Decision.** Nothing that appeared in this window changes a modelling choice. No public document names the final lines. No new same-lab Flex Perturb-seq matrix is downloadable. No new genome-scale CRISPRi transcriptome of a cancer line was deposited. The new method papers predict within a dataset we already hold, or they predict composition, drugs, or patient hold-outs. The teams ahead of us on the public board still publish no method, including nothing on fidelity, reach, Jaccard, or log-fold-change error.

Facts are **VERIFIED** (this session opened the page, the API record, or the PDF text) or **UNVERIFIED** (the page was a script shell, the PDF extract missed a section, or the claim is a competitor’s free-text description). Earlier scans (306, 482c, 520) are not re-listed except where a watch item was re-checked and is still empty.

Window: preprints and deposits dated 2026-09-15 through 2026-09-29. bioRxiv API `total` 3,443 and medRxiv `total` 1,493 were title-and-abstract filtered. Public leaderboard API fetched the same day (1,212 entries).

---

## Ranked by expected value

### 1. Final lines and new Arc Flex data — still unnamed. No action.

Re-checked 2026-09-29. This is the item that would have changed October if it had moved.

| Source | What was opened | Result |
|---|---|---|
| Arc news index, https://arcinstitute.org/news | Page text | Latest Virtual Cell Challenge post is still **20 August 2026**. The 2 September and 3 September posts are a 3D-genome profile and an anaerobic-metabolism story. No line name. **VERIFIED.** |
| Cell commentary, doi [10.1016/j.cell.2026.08.004](https://www.cell.com/cell/fulltext/S0092-8674(26)00931-1) | HTML full text (same sentences 520 relied on) | Six cancer lines, CRISPRi, 10x Flex. Three validation, three held for the test. No line names in the fetched text. **VERIFIED** that the fetched page still does not name D, E, or F. |
| Hugging Face `Xaira-Therapeutics/X-Atlas-Pisces` | API | `usedStorage` **0**, `lastModified` **2026-03-25T22:58:18Z**, siblings `.gitattributes`, `LICENSE.md`, `README.md`, card licence CC BY-NC-SA 4.0. Jurkat and HepG2 genome-scale CRISPRi are still “Coming Soon”. **VERIFIED.** |
| Hugging Face `arcinstitute` datasets, sorted by `lastModified` | API, 15 rows | Newest perturbation-relevant set is unchanged (Perturb-Sapiens 2026-01-09; State-Replogle 2025-12-06). `opengenome2` moved 2026-09-02 and is not Perturb-seq. **VERIFIED.** |
| `ArcInstitute/cell-eval2`, `ArcInstitute/state`, `Xaira-Therapeutics/X-Cell` | GitHub API | Last pushes **2026-08-22**, **2026-07-24**, **2026-03-16**. X-Cell preprint is still v1, dated 2026-03-20 (bioRxiv API). **VERIFIED.** |
| GEO, PDAT 2026-09-15 to 2026-09-30, Perturb-seq / CROP-seq / CRISPRi | esearch + esummary, 9 series | None is a human cancer gene-knockdown screen. The only Perturb-seq title is **GSE348514** (2026-09-24): “Xenium In Situ Spatial Perturb-Seq… mouse brain”, *Mus musculus*, 2 samples. Others: GSE348133 lncRNA knockdown in MDA-MB-231 (8 samples, not a screen), GSE335154 PDAC CAF trogocytosis, GSE312628 TTC5 iPSC neurons, GSE347678 *E. coli*, GSE291754 mouse microproteins, GSE347564 and GSE306808 iPSC morphogenesis, GSE347395 coronary variant-to-gene. **VERIFIED** from esummary titles. |
| BioStudies query `Perturb-seq` with release date in September 2026 | API | `totalHits` **0**. **VERIFIED.** |
| Figshare search `CRISPRi Perturb-seq`, newest first | API, first page | Newest records are **26 August 2026** (HSPC CRISPRi, figshare 31276981 / 31277494 / 32531301), before this window. **VERIFIED.** |
| Zenodo query `"Perturb-seq"`, most recent | API | In-window hits are not a new human screen. The closest biological record is 20 September, “An arrayed genome-wide CRISPR imaging screen of cGAS…”, imaging of cGAS foci, not a transcriptome. **VERIFIED.** |

**How this could have changed a decision.** A named D/E/F line with a public gene-level bank, or a new Arc Flex matrix, would turn on same-line rows (520, 451). Neither exists in the records opened here. October still rests on general transfer, above all signs.

### 2. Participant write-ups of the DE metrics — the leaders are still silent. No action.

Public leaderboard API, 2026-09-29, 1,212 current entries. Descriptions are free text. Treat every sentence as a self-report. **VERIFIED** as text present or absent on that API response. The biology in them was not re-measured.

The five entries ahead of the field’s mid-pack have **empty descriptions**:

| Team | Submitted | Overall | Model string | Description |
|---|---|---:|---|---|
| Illumina AI | 2026-09-29 | 0.4185 | `PerturbationAI` | empty |
| cqawesome | 2026-09-29 | 0.3282 | `m0949r` | empty |
| Testh | 2026-09-28 | 0.3061 | `test` | empty |
| CYH | 2026-09-29 | 0.2721 | `Cyan` | empty |
| chill | 2026-09-26 | 0.2708 | `Rain All Night` | empty |

Illumina’s raw members on that row: fidelity 0.561, reach 0.441, Jaccard 0.053, log-fold-change error 0.804. No text says how.

The `cellpred` cluster 525 already called private is still undescribed. Ava Brown, `cellpred-0.9.27`, submitted 2026-09-29, overall 0.2600, raw fidelity 0.528, Jaccard 0.034, log-fold-change error 0.890, description empty. Same fingerprint 525 reported (fidelity ~0.528, log-fold-change error ~0.88, Jaccard ~0.033). **No action:** still no public emitter to copy.

The highest score that does describe a DE-relevant change is **William, 0.2343, 2026-09-29, `transfer-protocol-v16a`**. The new sources named are KOLF2.1J (Nourreddine/Yeo 2026) and HipSci multi-donor iPSC CRISPRi (Feng et al.), each at weight 0.25, excluded from magnitude. Both are already on disk. KOLF is a standing no-go as a donor (FINDINGS). Feng is the iPSC bank already held. His raw fidelity is 0.509 (scaled −0.013). **No action.**

Below that, the descriptions that actually talk about fidelity, reach, Jaccard, or log-fold-change error restate mechanisms already closed, and they score 0.11–0.20:

- **Ming-Chang Chiu, 0.1618, `idea1-s6-callset-decoder`.** Call-set ranked by cross-source sign agreement, budget scaled to 0.25 of the previous budget. Raw fidelity falls to 0.447 (scaled −0.222) while log-fold-change error improves to 0.873. That is the yield penalty for under-calling. **No action** (192, 521).
- **Plan A, 0.1594, 2026-09-25, `v6.0.0-vcc82-c4b16-log2fc`.** Cross-source DE-consensus on log2 fold-change only. Their own note says reach +0.059 and log-fold-change error −0.017 on their recipe. Overall 0.159. Consensus signs are already at the public-data limit (516). **No action.**
- **OmicsMan, 0.1536, `R76_FlatCalls`.** Extra mid-expression genes flattened across the 400 cells so Wilcoxon calls them, with the pseudobulk unchanged. **No action.** 525 already found that the public emitters which can be opened lose the DE members.
- **Daniel Murray, 0.1905, `24-EchoCis50`.** Claims that on the Arc Flex platform, CRISPRi silencing extends past 2 kb, and caps predicted log2 change using 2025 H1 training bins: −0.22 at 2–10 kb and −0.12 at 10–50 kb, 505 neighbour pairs. **UNVERIFIED** as a measurement (leaderboard text only). Experiment 425 already measured this distance range on the Flex H1 truth (<1 kb mean log2 fold-change of significant calls −2.13; 5–20 kb −0.24) and closed an explicit cis term: an oracle on the <1 kb pairs is −0.0022 raw log-fold-change error, and the donor rows already carry the neighbours. **No action.**

GitHub search `virtual cell challenge 2026` (15 repositories). Repos pushed 27–29 September (`ferry-ai/virtual-cell-challenge-2026`, `mmyen/virtual-cell-challenge-2026`, `MMatinGerami/zero-shot-perturbation`) are the public emitters 525 already opened. Matin’s README, fetched this session, reports its own best overall at **0.079** (weighted transfer, 26 September). **No action.**

### 3. New papers on sign, programmes, or magnitude — read, then no action

Ordered by how close the abstract came to a cross-line CRISPRi sign or a per-target magnitude.

**SHERLOCK.** Mingxuan Zhang, Joshua D. Myers, Lingting Shi, Ross M. Giglio, Sharanya Chatterjee, José L. McFaline-Figueroa, Elham Azizi. Columbia. Posted **2026-09-26**. https://doi.org/10.64898/2026.09.25.754573. **VERIFIED** from the PDF (60 pages).

Shared transcriptional programmes here are clusters of perturbations *inside one screen*, plus condition-dependent drug responses in glioma, plus combinatorial CRISPRa magnitudes in K562 (Norman: 105 singles and 131 doubles). Held-out combination average-treatment-effect correlation is reported as r = 0.836. That is combination prediction in a line we already model, not the sign of a weak downstream gene in an unseen cancer line.

Data, from the PDF’s availability section: Norman via figshare; Replogle via the existing public release; spatial data **GSE261618**; “GBM kinase genetic and drug screening datasets will be deposited in GEO” — **not deposited in this paper**. Code https://github.com/azizilab/SHERLOCK (GitHub API: pushed **2026-07-26**, before this window, 337 MB, licence field empty) and https://github.com/azizilab/SHERLOCK_Reproducibility. Preprint licence **CC BY-NC-ND 4.0**. **No action.** Do not fetch the code or wait on the GBM deposit for October: the assay is not a new cancer-line CRISPRi bank of this panel.

**SPECTRA.** Michele Calabrò, Patrick Sheehan, Francesco Cambuli, Andrea Sottoriva. Human Technopole. Posted **2026-09-22**. https://doi.org/10.64898/2026.09.21.752624. **VERIFIED** from the HTML full text.

A directed-GRN autoencoder, trained on the **2025** Virtual Cell Challenge H1 set (about 300,000 cells, 150 CRISPRi targets), split by perturbation *inside that one line*. They report higher DEG F1 and AUPRC than GEARS, scLAMBDA, and a pseudobulk, and they report **worse unweighted MSE** than that pseudobulk (about 0.0044 versus under 0.003). On the 2025 `vcc eval` suite their highlighted numbers are PR-AUC 0.210 and ROC-AUC 0.472. They explicitly de-emphasize overlap-at-N. This does not measure fidelity, reach, Jaccard, or log-fold-change error on an unseen cancer line, and it does not beat a same-target transfer.

Code https://github.com/MCalabroCode/SPECTRA. GitHub API: pushed **2026-09-23**, **MIT**, about **100 MB**. Training data is the 2025 challenge set, already held. **No action.** Not worth the download.

**DELPHAI v2.** Xian Zhang, Hui Wu, Haikun Liu. aiPTO / DKFZ. Version posted **2026-09-17**. https://doi.org/10.64898/2026.07.01.735965. **VERIFIED** from the PDF text that extracted.

Ranks first on DEG overlap inside scPerturBench and scArchon: 51.3% of the top 100 DEGs in the in-distribution split, **33.2%** when the cellular context is held out, against scPRAM as runner-up. The out-of-distribution examples they discuss are **IFN-β on Kang (cross-patient)** and **CEP18770 on Kaggle (cross-patient)** — cytokine and drug, patient hold-out, not CRISPRi in an unseen cancer line. Licence on the PDF footer: **CC BY-NC 4.0**. The only code URL in the extract is the benchmark harness https://github.com/bm2-lab/scPerturBench. A DELPHAI repository was **not found** in the extracted text (**UNVERIFIED** whether one exists). **No action.**

**P2P.** Haojie Yang, Ran Su. Tianjin University. arXiv [2609.34391](https://arxiv.org/abs/2609.34391), posted **2026-09-28**. **VERIFIED** from the PDF. Licence in the PDF metadata: **CC BY 4.0**.

Population-level denoising on Adamson, Norman, Replogle K562, and Replogle RPE1. On Replogle K562, effect Pearson 0.643 → 0.702 and DEG F1 0.067 → 0.178 versus Scouter. Those are within-dataset genetic screens we already use as donors. No code URL was present in the extracted text (**UNVERIFIED**). **No action.**

**IGNITE / CD4 paper.** Angela L. Han and colleagues, including Joshua Welch. University of Michigan. Posted **2026-09-16**. https://doi.org/10.64898/2026.09.12.748195. **VERIFIED** from the HTML full text.

This is the write-up of the **22 million primary CD4+ T-cell** Perturb-seq already held (GSE314342 / the CZI genome-scale T-cell set). Data-availability sentence: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq. The CZI landing page itself returned only a script shell, so the **file list on that URL is UNVERIFIED** this session; the paper’s sentence and the match to the known 22-million-cell study are verified. Code https://github.com/alhanster/IGNITE_public (pushed 2026-09-14, MIT, about 10 MB) ranks immune drug targets. Scores also at https://ignite.eecs.umich.edu/. **No action.** It does not add counts, and it does not estimate knockdown signs in a cancer line.

**PerturbBridge.** Zhihao Liu, Changzhi Jiang, Can Yang, Xiangrong Liu. Xiamen University. Posted **2026-09-25**. https://doi.org/10.64898/2026.09.24.754259. **VERIFIED** from the PDF (9 pages). Norman CRISPR and Sci-Plex3 chemical benchmarks only. No code URL in the extracted text. Footer: “All rights reserved. No reuse allowed without permission.” **No action.**

**Pop-Corn.** Jiapeng Chen, Yan Cui, Yanjun Shao, Na Sun, María Rodríguez Martínez. Yale / Whitehead. Posted **2026-09-26**. https://doi.org/10.64898/2026.09.25.754560. **VERIFIED** from the PDF. Licence **CC BY 4.0**.

Predicts **cell-state composition**, and the abstract says expression-prediction models that get the average transcriptome right still fail at composition. Our score is the transcriptome. Data: CRISPRa T cells **GSE190604** (already public for years) and spatial CRISPR **GSE275483** / Hugging Face `xingjiepan/PerturbMulti`. Stated code https://github.com/GabbyKoki/pop returned **HTTP 404** this session. **No action.**

**One-line no-actions** (abstract verified from the bioRxiv API; full text not opened):

- Speciesformer, Wang et al., 2026-09-23, doi 10.64898/2026.09.22.752128. Cross-species foundation model on 131 million observational cells. Not a CRISPRi bank.
- Targeted Pol III CRISPRi, Simeonov et al., 2026-09-27, doi 10.64898/2026.09.24.754042. Fitness screen in fibroblasts, HEK293T, and glioblastoma. Not a transcriptome panel.
- T-promoter CRISPRi of CTNNB1, Jin and Palecek, 2026-09-28, doi 10.64898/2026.09.25.754435. One circuit in human pluripotent stem cells.
- RAS+ digenic knockout, 2026-09-25, doi 10.64898/2026.09.24.754062. Growth phenotypes in 12 lines, not single-cell expression magnitudes.
- EBV host-factor CRISPR, Wang, Yeo, et al., 2026-09-28, doi 10.64898/2026.09.28.754942. B-cell infection fitness. The “Yeo” here is not the KOLF atlas.
- ISP / Geneformer platform, Sanaki et al., 2026-09-24, doi 10.64898/2026.09.17.752497. In-silico overexpression screens of Yamanaka factors. No new counts.
- PRAXIS-VirtualCell, arXiv 2609.29221, 2026-09-24. An agent framework. PDF opened; it does not release a perturbation matrix. **No action.**

---

## Download manifest

**Fetch nothing for a decision.** The rows below are the only objects that were candidates, and each fails the bar.

| Item | URL | Size | Licence | Why it stays unfetched |
|---|---|---|---|---|
| X-Atlas/Pisces | https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces | 0 bytes (`usedStorage` 0) | CC BY-NC-SA 4.0 | Still three metadata files. Re-check `usedStorage` before October; do not download today. |
| SPECTRA code | https://github.com/MCalabroCode/SPECTRA | ~100 MB | MIT | 2025 H1, within-line split, metrics that are not the 2026 six. |
| SHERLOCK code | https://github.com/azizilab/SHERLOCK | ~337 MB | unset | Code predates the paper. Training data are Norman, Replogle, and an undeposited GBM drug screen. |
| IGNITE code | https://github.com/alhanster/IGNITE_public | ~10 MB | MIT | Drug-target ranker on the CD4 atlas already held. |
| Pop-Corn code | https://github.com/GabbyKoki/pop | — | CC BY 4.0 paper | Repository 404. Composition task. |
| CZI CD4 counts | https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq | already held | paper points here | Same 22-million-cell CD4 study. File list on the landing page unverified this session; do not re-download. |
| DELPHAI, P2P, PerturbBridge weights | no URL in the extracted text | — | CC BY-NC 4.0 / CC BY 4.0 / all rights reserved | No artifact to fetch. Wrong benchmark. |

Watch, do not download: Pisces `usedStorage`, a GEO accession for the SHERLOCK GBM screens, and any Arc news post after the 20 August challenge announcement.
