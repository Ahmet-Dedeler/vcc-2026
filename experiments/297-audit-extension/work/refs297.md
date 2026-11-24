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
