# 320 — sources for the F16 novelty audit

All access dates **2026-09-08**. Provenance codes, following 287/297's convention:

- **[F]** full text fetched and read this session (WebFetch on the URL given; quotes are extractions from that fetch).
- **[P]** PDF fetched and read as pages this session (quotes read off the page images).
- **[S]** search-index summary only — the source itself was **not** opened. Treated as weaker evidence and
  never used as the sole support for a verdict.
- **[X]** unreachable (paywall / 403 / IdP redirect). Contents are **not** inferred from the title.

No claim in `RESULTS.md` rests on a source I did not open, except where explicitly labelled [S] in the
justification column.

---

## A. The two papers that decide items 1, 2 and 4

### ROHATGI24 — the primary CRISPRi seed off-target report

Rohatgi N, Fortin J-P, Lau T, Ying Y, Zhang Y, Lee BL, Costa MR, Reja R.
**"Seed sequences mediate off-target activity in the CRISPR-interference system."**
*Cell Genomics* **4**(11):100693 (2024). doi:10.1016/j.xgen.2024.100693. PMID 39510079.
First published 2024-11-06.

- Journal version: `https://www.sciencedirect.com/science/article/pii/S2666979X24003227` — **[X] HTTP 403.**
  `https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00322-7` — [X] (287 recorded cell.com 403 as well).
- Preprint read in full: `https://www.biorxiv.org/content/10.1101/2024.04.10.588881v1.full` — **[F]**.
- Bibliographic record: `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:39510079&resultType=core&format=json` — **[F]**.
- **All quotes below are from the bioRxiv preprint, not the journal version.** This is `REFERENCES.md`'s
  `SEED24`, which carried provenance **[S]** and no author list; 320 upgrades it to [F] on the preprint and
  supplies the full citation.

Sentences relied on:

1. *The rule.* "10-12 nucleotides located at the 3'-end of sgRNA sequence, often referred as the 'seed'
   sequence and presence of adjacent Protospacer Adjacent Motif (PAM) sequence in the genome are crucial for
   target recognition."
2. *The promoter window.* They search the "9 bp PAM proximal seed sequence" against promoters "defined as
   +/- 2000 bp centered around the TSS", keeping "sgRNA matches that were within 2000 bp of a TSS".
3. *The rule evaluated as a predictor, at one locus.* "There were approximately 562 sgRNAs in the library with
   seed sequences that mapped to the GSDMD promoter (defined as -2000 to +2000 bps relative to the GSDMD TSS)
   in the pyroptosis screen, however, only 21 of these displayed fold-change values equivalent to that
   observed for the on-target sgRNAs."  → **21/562 ≈ 3.7 %.**
4. *Why the rule under-predicts.* "Many of the sgRNAs that did not display a high fold change had seed
   sequence matches close to the TSS and in regions of open chromatin, suggesting that there must be factors
   in addition to 9-nt seed sequence complementarity that determine the binding stability of a sgRNA to
   genomic DNA."
5. *Occupancy evidence.* "All four guides yielded multiple reproducible peaks ranging from 98 for sgGRK4 to
   1141 for sgADAT1. Of these peaks, the proportion of promoter proximal peaks ranged from 10% for sgLRRC34
   to 33% for sgGRK4."
6. *The 12-nt threshold, with occupancy.* "Using ChIP-seq we showed that sequence complementarity that
   extends as far as 12 bp of seed sequence can sometimes be crucial for a strong occupancy of dCas9 on the
   DNA and for off-target transcriptional repression."
7. *Same parent library as ours.* They "infected an EA.hy926 cell line…with the Horlbeck and Weissman *et al.*
   whole-genome CRISPRi sgRNA library."
8. *The partner-guide logic as their discovery device.* "Upon examining other enriched guides in each screen,
   we noticed certain guides displayed a high fold change that deviated from the other guides designed to
   target the same gene."
9. *Mismatch tolerance.* "depending on the sequence context, mismatches in the seed sequence of the guide do
   not necessarily impact its binding to the DNA."
10. *Derivation was post hoc.* Off-targets were identified after screening; no prospective predictive model
    was built (fetch's explicit finding).

### HARTMAN26 — the systematic Perturb-seq seed off-target workflow

Hartman A, Blair JD, Nguyen TP, Dyson K, Bradu A, Takacsi-Nagy O, Santostefano K, Boade T, Bolaños M, Zhu R,
Dann E, Marson A, Gitler AD, Satija R, Satpathy AT, Roth TL.
**"Systematic identification of seed-driven off-target effects in Perturb-seq experiments."**
bioRxiv 10.64898/2026.03.27.714658, v2 posted 2026-03-30.
`https://www.biorxiv.org/content/10.64898/2026.03.27.714658v2.full` — **[F]** (fetched twice, different prompts).
Code: `https://github.com/AustinHartman/perturb_seed` (not opened).
Author list from the search result for the preprint title — **[S] for the author list only**; the title,
venue and all quoted sentences are [F].

Sentences relied on:

1. *Their rule.* Seed matches are "the 10-12 PAM-proximal bases", searched as "≥5 contiguous matching bp
   adjacent to a PAM" within "a ±2,000 bp window (totaling 4,000 bp) surrounding the annotated TSS";
   high-confidence calls filter to "seed matches of 12 bp or more (and excluded matches of 19 and 20 bp)".
2. *Background rate of the sequence rule.* "47.2–59.8% of guides had at least one match" (≥12 bp, within
   1,000 bp of a protein-coding TSS).
3. *Rate of actual events.* "0.83% (range: 0.2%-1.8%) of guides were affected when we filtered to only include
   seed matches of 12-18 bp" — across the K562, HCT116 and CD4⁺ T-cell GWPS datasets.
   → **implied precision of the sequence rule alone ≈ 1.4–1.8 %.**
4. *Sequence alone is not enough.* "transcriptional readouts from existing Perturb-seq datasets would provide
   a means of filtering the search space"; the workflow is applied **post hoc** to expression-nominated events
   ("we asked whether it repressed any of the genes targeted by the neighboring guides" and only then
   "searched for candidate off-target binding sites near TSSs").
5. *The partner-guide specificity control, in CRISPRi Perturb-seq.* "the TCR inactivation phenotype was
   observed exclusively for the APPL2 and LRBA guides…and not for the remaining four guides (WDR53-1,
   WDR53-2, LRBA-1, and APPL2-1), despite repression of their intended targets"; and "only the shared APPL2
   guide reduced LAT levels, while both guides effectively repressed APPL2".
6. *Their unit and test.* The unit is the **guide–gene pair**; "we ran single cell two-sided Wilcoxon rank-sum
   tests"; the fetch found **no** covariate-matched resampling of candidate gene sets and **no** global
   aggregate statistic across all guides.
7. *Positional enrichment.* "TSS enrichment was markedly stronger among high-confidence candidates (seed match
   ≥12 bp and log2 fold-change < -0.5)"; "even among candidates with shorter seed matches…we still observed
   significant enrichment within approximately 1,000 bp of the annotated TSS".
8. *Their own caveat.* "the number of off-target events identified in a dataset should not be interpreted as a
   measure of data quality."

---

## B. The RNAi prior art for items 3, 4 and the "protospacer not target" headline

### JACKSON03 — the founding sequence-not-target observation

Jackson AL, Bartz SR, Schelter J, Kobayashi SV, Burchard J, Mao M, Li B, Cavet G, Linsley PS.
**"Expression profiling reveals off-target gene regulation by RNAi."**
*Nature Biotechnology* **21**(6):635–637 (2003). doi:10.1038/nbt831. PMID 12754523.
Received 10 Feb 2003; accepted 18 March 2003; published online 18 May 2003.
`https://arep.med.harvard.edu/pdf/Jackson03.pdf` — **[P]** (WebFetch returned the raw PDF; read as 3 pages).

Sentences quoted off the page:

1. Abstract: "Transcript profiles revealed siRNA-specific rather than target-specific signatures, including
   direct silencing of nontargeted genes containing as few as eleven contiguous nucleotides of identity to
   the siRNA."
2. "different siRNAs to the same target transcript elicit a small number of gene regulations in common (data
   not shown), but the vast majority of the transcript expression patterns were siRNA-specific rather than
   target-specific."
3. "Thus, we have observed patterns of gene regulation that are specific for the siRNA sequence used for
   silencing, rather than for the intended target."
4. "Surprisingly, our analysis revealed few genes regulated in common by different siRNAs to the same target
   gene. Each of the 8 siRNA duplexes targeted to *MAPK14* produced a distinct expression pattern (Fig. 1a).
   Likewise, each of the 16 siRNAs targeted to *IGF1R* produced a unique expression pattern (Fig. 1b)."
5. Their own prescription, which is the partner-construct control as standard practice: "incorporating
   multiple siRNA duplexes to silence a target gene of interest will increase the confidence with which an
   observed phenotype and expression pattern can be linked to target gene silencing."

### BIRMINGHAM06 — seed matches, not identity, drive RNAi off-targets

Birmingham A, Anderson EM, Reynolds A, Ilsley-Tyree D, Leake D, Fedorov Y, Baskerville S, Maksimova E,
Robinson K, Karpilow J, Marshall WS, Khvorova A.
**"3′ UTR seed matches, but not overall identity, are associated with RNAi off-targets."**
*Nature Methods* **3**(3):199–204 (2006). doi:10.1038/nmeth854. PMID 16489337.

- `https://www.nature.com/articles/nmeth854` — **[X]** 303 redirect to `idp.nature.com` (not followed).
- Citation + abstract fragment from EuropePMC REST — **[F]**:
  "off-targeting is associated with the presence of one or more perfect 3′ untranslated region (UTR) matches
  with the hexamer or heptamer seed region."
- **I could not read its Methods**, so I do **not** quote its statistical design; the "seed-match genes vs a
  background set" idiom is attributed below to a source I did open (SEEDMATCHR24).

### SEEDMATCHR24 — the standard aggregate seed-shift test, and its unit

Cazares T, Higgs RE, Wang J, Ozer HG.
**"SeedMatchR: identify off-target effects mediated by siRNA seed regions in RNA-seq experiments."**
*Bioinformatics* **40**(1):btae011 (2024). doi:10.1093/bioinformatics/btae011.
`https://pmc.ncbi.nlm.nih.gov/articles/PMC10799297/` — **[F]**.

1. "Comparison of the ECDFs for genes with a seed match to a background set has been applied to assess the
   effects of seed-mediated binding on gene expression profiles in miRNA and siRNA experiments."
2. "we expect the distribution of log2(fold changes) for transcripts with seed matches to be shifted to the
   left (downregulated) compared to a background set. We use a one-sided KS test to determine if the ECDF of
   transcripts with a seed match is less than the background distribution."
3. Fetch's explicit findings: **no covariate matching** is used or discussed; the **unit of analysis is the
   gene**; the design is credited to Grimson 2007, Janas 2018, McGeary 2019, Schlegel 2021/2022.

---

## C. Statistical-design prior art for items 3 and 5

### GOEMAN07 — competitive vs self-contained nulls; gene- vs subject-sampling

Goeman JJ, Bühlmann P. **"Analyzing gene expression data in terms of gene sets: methodological issues."**
*Bioinformatics* **23**(8):980–987 (2007).
`https://academic.oup.com/bioinformatics/article/23/8/980/198511` — **[F]**.

1. Competitive null: "The genes in G are at most as often differentially expressed as the genes in Gc."
   Self-contained null: "No genes in G are differentially expressed."
2. The two sampling models: gene sampling "uses genes as independent sampling units"; subject sampling "uses
   subjects as independent sampling units". "A replication of the experiment under a subject-sampling model
   would involve a new sample of subjects…the experiment is repeated for new subjects, but with the same
   genes."
3. Why gene sampling is anti-conservative: "Correlations between gene expression measurements of genes tends
   to result in positive correlations between their…P-values, which in turn causes their indicators of
   differential expression to be correlated…the use of the hypergeometric test is anti-conservative; it may
   greatly understate the true P-values if the genes in the gene set are not independent." Their simulation
   reached rejection rates "up to 50,000 times the nominal level" under gene correlation.

### MATCHRANGES23 — covariate-matched null loci as a published method

Davis ES, Mu W, Lee S, Dozmorov MG, Love MI, Phanstiel DH.
**"matchRanges: generating null hypothesis genomic ranges via covariate-matched sampling."**
*Bioinformatics* **39**(5):btad197 (2023). doi:10.1093/bioinformatics/btad197. PMID 37084270.
Preprint read: `https://www.biorxiv.org/content/10.1101/2022.08.05.502985v1.full` — **[F]**.
Bibliographic record via EuropePMC REST — **[F]**.

1. "Accurate inference requires that null sets exhibit similar distributions of covariates observed in the
   focal set, to mitigate interpretability issues due to confounding."
2. "matchRanges computes for each range a propensity score, the probability of assigning a range to focal or
   background groups, given a chosen set of covariates. It provides three methods including nearest-neighbor
   matching, rejection sampling, and stratified sampling for null set selection."
3. **Fetch's explicit finding:** the preprint "does not discuss computing empirical p-values through
   resampling multiple matched null sets. It focuses on generating a single matched null set." So the
   *matching* is published as a method; the *repeated-draw empirical p on matched sets* is not part of it.

### AGGREGATION — why a set-level test can reject where per-item tests cannot

**[S]** — search-index reading only; used as background, not as a load-bearing citation. The canonical
references surfaced were the rare-variant aggregation literature (Lee S, Abecasis GR, Boehnke M, Lin X,
*Rare-Variant Association Analysis: Study Designs and Statistical Tests*, Am J Hum Genet 95:5–23, 2014;
`https://www.sciencedirect.com/science/article/pii/S0002929714002717`) and a 2025 AJHG power analysis
(`https://pmc.ncbi.nlm.nih.gov/articles/PMC12414686/`). The standard points recovered: single-item tests are
underpowered and carry the multiplicity burden, aggregation both accumulates signal and reduces the number of
tests; a set-level rejection is evidence for the set and not for any member; results are contingent on the
pooling definition ("mask sensitivity"); burden-style statistics assume same-direction effects.
Also recovered, on permutation resolution: "permutations give no basis to distinguish between two gene sets
that both hit the floor p-value of 1/(M+1)".

---

## D. Library, mechanism and challenge-data sources

### HORLBECK16 — the parent library, and what it did about off-targets

Horlbeck MA, Gilbert LA, Villalta JE, Adamson B, Pak RA, Chen Y, Fields AP, Park CY, Corn JE, Kampmann M,
Weissman JS. **"Compact and highly active next-generation libraries for CRISPR-mediated gene repression and
activation."** *eLife* **5**:e19760 (2016). doi:10.7554/eLife.19760.
`https://elifesciences.org/articles/19760` — **[F]**.

1. Off-target filter at design time: "sgRNAs were required to have no more than 1 alignment (the sgRNA target
   site itself) with 'mismatch score' less than 31 proximal to the TSS and under 21 in the genome."
2. "For hCRISPRi-v2, 96.6% of sgRNAs incorporated passed at this threshold."
3. Fetch's explicit findings: **no empirical off-target rate is reported**, and **no direct measurement of
   off-target transcriptional effects was performed**. The predictive scoring used nucleosome positioning,
   TSS distance, sequence features, secondary structure, protospacer length and chromatin accessibility —
   an *on-target* activity model.

### JOST20 — why the seed is the seed

Jost M, Santos DA, Saunders RA, Horlbeck MA, Hawkins JS, Scaria SM, Norman TM, Hussmann JA, Liem CR,
Gross CA, Weissman JS. **"Titrating gene expression using libraries of systematically attenuated CRISPR guide
RNAs."** *Nat Biotechnol* **38**(3):355–364 (2020). doi:10.1038/s41587-019-0387-5. PMID 31932729.
**[S]** — search-index summary of the paper and its bioRxiv version (`10.1101/717389`); not opened.
Point used: sgRNAs with mismatches in the PAM-proximal seed region had strongly compromised activity — the
mechanistic basis for treating the seed as the binding determinant. Corroborating only.

### OFFTARGET-PREDICTORS — CFD, Elevation and why they are the wrong family here

**[S]** — search-index reading; not opened. Doench JG *et al.* 2016 (CFD, Nat Biotechnol 34:184–191) and
Listgarten J *et al.* 2018 (Elevation, Nat Biomed Eng) are **cleavage-outcome** models: CFD is calibrated on
off-target *cutting* and validated with GUIDE-seq. The recovered field statements are that "Cas9-DNA binding
requirements are different from Cas9-mediated DNA double strand breaks" and that "occupancy is not a good
predictor of cleavage… mismatches in the PAM-distal region do not affect binding, but impair cleavage".
Affinity-based scores (CRISPRoff/CRISPRspec) and chromatin-aware binding predictors (CROP-IT) are the
closer family for dCas9. Used in `RESULTS.md` only to say what the *off-target-score* literature is **not** —
it is not a published version of the promoter-seed repression rule.

### ARC-BTD — what the challenge's own data notes say about their library

Arc Institute, **"Behind the Data of the Virtual Cell Challenge."**
`https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge` — **[F]**.

1. "Using a dual-guideRNA expression construct, where two guides targeting a gene of interest are
   simultaneously expressed from the same vector."
2. Protospacers come from elsewhere: "protospacer sequences and cloning strategy from (Replogle, Bonnar,
   et al. 2022)"; "Having selected 300 target genes and their protospacer sequences" the library "was cloned
   in a similar manner as our pilot screen and validated for uniformity and coverage."
3. QC is efficacy, not specificity: "83% of cells with >80% knockdown".
4. **Fetch's explicit finding: "Off-target effects, specificity filtering, and off-target impact assessments
   are not mentioned in this document."**

Also read at the search-index level **[S]**: the 2025 and 2026 challenge overviews
(`https://www.cell.com/cell/fulltext/S0092-8674(25)00675-0`,
`https://www.cell.com/cell/fulltext/S0092-8674(26)00931-1`,
`https://arcinstitute.org/news/virtual-cell-challenge-2026`) — no dedicated off-target/specificity analysis of
the challenge guide library was surfaced by any query.

### ESCAPEEXCEL17 — the origin and the documented class of the `="SYMBOL"` defect

Welsh EA, Stewart PA, Kuenzi BM, Eschrich JA. **"Escape Excel: A tool for preventing gene symbol and
accession conversion errors."** *PLoS One* **12**(9):e0185207 (2017). doi:10.1371/journal.pone.0185207.
PMID 28953918. Citation via EuropePMC REST — **[F]**; the paper itself **[S]**.
Point used, from the search-index reading: the `="MARCH1"` formula-quoting style *is* Escape Excel's escaping
strategy, and "some supplementary tables ship with literal `="GENE"` strings: the file was escaped for safe
Excel import, but downstream parsers reading the raw text (rather than through Excel) see the formula wrapper
as part of the value — a defect in itself if not stripped." Searches for this specific instance in
Horlbeck 2016's supplements returned nothing.

---

## E. On-disk documents read (no writes outside `experiments/320-f16-novelty/`)

`docs/READ-FIRST.md` (v5.4 + the 09-08 validity limit), `docs/FINDINGS-2026-09.md` §F16 and §3,
`experiments/307-offtarget-enrichment/RESULTS.md`, `experiments/307-offtarget-enrichment/PREREGISTRATION.md`
(the covariate table and the window/filter definitions),
`experiments/309-c9-calibration-repair/RESULTS.md`, `experiments/301-metadata-acquisitions/RESULTS.md`
(the `[−50, +300]` CRISPRi active window and the 10 kb distal filter),
`experiments/287-novelty-audit/RESULTS.md` and `work/refs.md`,
`experiments/297-audit-extension/RESULTS.md`, `experiments/294-writeup/PAPER.md` §3(g) and §4,
`experiments/294-writeup/REFERENCES.md`.

## F. Compute

**None.** No `uv run python` was run; every number in `RESULTS.md` is either quoted from a source above or
copied from 301/307/309's own committed RESULTS/PREREGISTRATION files. No metadata field needed reading, so
no known-answer row of my own was required — the audit's known-answer rows are §G below, which are
reproductions of numbers the programme already published, checked against the files that hold them.

## G. Known-answer rows for this audit (reading checks, no compute)

| # | row | expected (as published) | read today from the cited file | verdict |
|---|---|---|---|---|
| KA-A | 307's primary `S` and its interval | +0.78567 [+0.38669, +1.16604], n = 34 | `E307/RESULTS.md` §4 table, identical | **PASS** |
| KA-B | 309's repaired `S`, p and KA count | +0.71293 [+0.32419, +1.07351], p 1.99998e-05, 9/9 | `E309/RESULTS.md` "THE LINE", identical | **PASS** |
| KA-C | the partner-construct difference | +0.74500 [+0.22896, +1.25766], Wilcoxon p 0.00459 | `E309/RESULTS.md` abstract; `E307` §6 P3 gives +0.74461 [+0.22730, +1.24514] on C8's pools — the two differ because C9 uses repaired pools, as both files state | **PASS** |
| KA-D | 301's selection reproduced by 307's KA-2 | failing 4/13, working 4/21, NTC 4/31 | `E307/RESULTS.md` §2 KA-2 and `E309` KA-2(a), identical | **PASS** |
| KA-E | the seed rule's window, as implemented | 12/11/10/9-nt seed + NGG, offset in **[−50, +300]**, 10 kb distal filter | `E301/RESULTS.md:149-152`; `E307/PREREGISTRATION.md:67-69`, identical | **PASS** |
| KA-F | the measured per-spacer L = 12 hit rate | real 0.132353 / pool 0.160879; bug-equivalent 0.597 | `E309/RESULTS.md` KA-2(b), identical | **PASS** |

KA-E and KA-F are the two rows the literature comparison depends on: KA-E is why our per-spacer rate is not
comparable to Hartman's, and KA-F is the number the Horlbeck escape defect moved.
