# 287 — Novelty audit of F1–F12, C1, C2 against the original dataset papers and the transfer literature

**Verdict: NEW 4 · EXTENDED 11 · KNOWN 7 · CONTRADICTED 2.** Reading and citation only, no compute beyond four
read-only checks on published supplementary tables. Deliverable: `NOVELTY-AUDIT.md` (the 24-row table, the
"what can be claimed as new" section, the contradictions section). References with URLs and access dates:
`work/refs.md`. Checks: `work/checks287.py`, log `work/checks287.log`. Nothing written outside this directory;
no `vcc`, no `cell-eval2`, no `data/` write, no heavy job.

## 0. The four things a coordinator needs

1. **Two of our four surfaces publish the reliability numbers we report.** F6's 85/272 = **31 %** of panel targets at
   or above the 10-gene DE gate is, to the digit, Replogle 2022's own published **30.5 %** of 9,608 K562 perturbations
   causing > 10 DEGs. X-Atlas/Orion's own preprint reports **48.36 / 40.37 / 54.30 %** of on-target perturbations
   (HCT116 / HEK293T / K562) with a detectable single-cell transcriptomic change — the complement is F1's "a third of
   targets ≤ 0", in the same size regime. And a 2026 preprint (Wang, Kuipers, Hugi, Platt, Beerenwinkel, bioRxiv
   2026-08-11) already reports **65 % of 7,170 perturbations across 29 datasets unreliable**, by split-half Pearson with
   a Spearman–Brown correction — our construction, on a superset of our data. **F1a and F6 are KNOWN.**
2. **The Broad-16 source paper (identified: Maffa et al., *Dissecting context-dependent cancer vulnerabilities using
   Perturb-seq*, bioRxiv 10.64898/2026.08.24.746802) contradicts F5's headline and pre-explains C2a's failure.** It
   predicts line-specific responses from cell state — p53 activation on UBE3A/MDM2 only in the HPV-E6 lines SKGII and
   C4I; IER3IP1 dependence from a hemizygous 18q/11q co-deletion at p = 3.28 × 10⁻⁶ — while our null is for a linear
   ridge on two basal *expression* programmes. It also publishes the depth–power law that 283 rediscovered post hoc
   ("strongest predictor of statistical power was the total number of UMIs recovered", r = 0.85; elbow ≈ 5 M
   UMIs/perturbation): **C2a's threshold should have been depth-conditioned before it was frozen.** Its two dominant
   programmes are the same E2F/G2M-vs-p53/stress pair F5 reports as its two largest components.
3. **One sentence in `docs/FINDINGS-2026-09.md` must be deleted.** F10's "*What it changes*" still says "the depth axis
   is closed on this disk: cells are not the binding constraint at any genome-scale budget." That survived the 13:40Z
   retraction of 280's depth curve two paragraphs above it, and it is contradicted by Wang et al.'s measured
   84 % → 10 % unreliable law over cell count and by Broad-16's UMI–power relation. **CONTRADICTED; delete now,
   285 supplies the replacement.**
4. **Four items survive as new** (`NOVELTY-AUDIT.md` §2): F3b (rank retrieval *needs* the noise-dominated directions of
   its own whitening operator — the one result that runs against the standard shrinkage-whitening prescription),
   F11c (cross-line conservation does not predict within-individual state dependence, preregistered), F10b (the donor
   is a marginal replicate unit; three donors beat one at matched UMI budget), and F9 as a *quantification* of a
   qualitatively known axis. F7 is new but is private-scorer calibration, not science.

## 1. Verdicts

| verdict | items |
|---|---|
| **NEW (4)** | F3b (damping the noise-dominated directions loses `pds`); F7 (board class factors — not a scientific claim); F10b (donor as a marginal unit; replicates-vs-cells at matched budget); F11c (conservation does not predict across-state agreement) |
| **EXTENDED (11)** | F1b (Flex-400 contrast), F2 (retained fraction of a cross-validated oracle), F3a (operator spectrum), F4a (the `pds` knee), F4b (the deficit decomposition), F5a (three-way energy decomposition), F8 (conservation predicts transfer with the line held out), F9 (family-resolved conservation axis), F11b (crossed donor × guide × state 2 × 2), F12a (family-level predictability), F12b (reliability dominates the gene-level residual) |
| **KNOWN (7)** | F1a (150–200-cell reliability), F1c (guide construct), F6 (K562 DE sparsity), F10a (low across-donor agreement in CD4), F11a (state dependence exists), C1 (a split-half robustness check), C2's depth explanation |
| **CONTRADICTED (2)** | F5b's headline "not basal-state-predictable" vs Broad-16's own genotype/viral-status predictions; F10's surviving "cells are not the binding constraint" sentence vs Wang 2026 and Broad-16 |

## 2. Checks run (read-only, `work/checks287.log`)

**KA-1 PASS** — the GWCD4i K562-comparison table's three random controls read +0.0040 / +0.0037 / +0.0038 against a
real cross-cell-type median of +0.3030 / +0.3147 / +0.3117 (Rest / Stim8hr / Stim48hr, n = 1,880 each).

| check | result | which finding it touches |
|---|---|---|
| The authors' own CD4 donor and guide reproducibility columns | `donor_correlation_all_mean` median **+0.107 / +0.112 / +0.114** all-genes and **+0.103 / +0.091 / +0.126** on the 300 panel targets; `donor_correlation_hits_mean` **+0.705 / +0.725 / +0.758**; `guide_correlation_all` **+0.062 / +0.072 / +0.061** (share ≤ 0 = 19 / 11 / 13 %) | F10a → **KNOWN**; F1c, F11b (cross-guide is worse than cross-donor in their table too, the same ordering as 282's crossed cells) |
| ZNF vs the rest inside the authors' K562-vs-CD4 table | Δmedian **−0.0304 / −0.0317 / −0.0320**, AUC(other > ZNF) **0.581 / 0.583 / 0.587**, one-sided p 0.019 / 0.017 / 0.013; bootstrap interval below 0 in Stim48hr only | F9 — an **independent-estimator corroboration** (DESeq2 z-scores, one cross-line pair, 56 ZNFs) of 279's AUC 0.564. Latent in a published table; not reported by its authors |
| Conservation vs cross-donor reproducibility, their two columns | ρ **+0.167 [+0.084, +0.247] / +0.112 [+0.033, +0.189] / +0.110 [+0.033, +0.186]** | F11's "conservation predicts the across-donor ceiling (ρ +0.301)" — same sign, weaker, on their estimator |
| The magnitude confound | ρ(conservation, MASH DE count) **+0.507 … +0.543**; ρ(donor reproducibility, DE count) **+0.242 … +0.316** | F8 — conservation tracks response magnitude ~2× more strongly than it tracks reproducibility. **F8 has no magnitude-matched arm; F9 does. Add one.** |
| Reproducibility of a fetched quote | The preprint figure "median cross-guide correlation 0.49 / 0.47 / 0.52" does **not** reproduce from the authors' table under five plausible filters (all read +0.60 … +0.76 on significant sets, +0.06 … +0.07 all-genes) | Provenance: quotes obtained through automated page fetches are flagged in `NOVELTY-AUDIT.md`; this one is not used |

## 3. Actions for the coordinator, in priority order

1. **Delete** F10's "cells are not the binding constraint at any genome-scale budget" —
   `docs/FINDINGS-2026-09.md:302`. It is the *only* surviving copy: `docs/READ-FIRST.md:69` already carries the
   13:40Z retraction verbatim ("the '4,500 cells' ceiling and 'cells are not the binding constraint' are NOT
   established … The depth axis is OPEN pending 285"), so FINDINGS and READ-FIRST currently disagree with each
   other as well as with Wang 2026 and Broad-16.
2. **Narrow F5's headline** to "not predictable from basal *expression* by a linear model at 11 training lines", and
   record the Broad-16 counter-examples. Optional and cheap on disk: re-run 268's LOLO with lineage, HPV status and
   arm-level copy number (`data/external/ccle_broad_2025`, `data/external/depmap26q1`) as the state features. That is
   the measurement that would settle it.
3. **Add a magnitude-matched arm to F8** (recompute the LOLO conservation → per-target `pds` relation within DE-count
   strata using `E276/work/conservation276.csv`). The source's own table says the confound is large.
4. **Re-frame F1a, F6, F10a and F11a as re-measurements** with citations, not as findings. In particular F1 must state
   explicitly that Orion's headline median Spearman of 0.99 is NTC *profile* stability, not response reproducibility —
   otherwise a reviewer will read F1 as contradicting the source paper.
5. **Separate F2 from Wang 2026** in wording: F2 rejects per-target selection of a donor row by proximity to the truth;
   it does not reject reliability-based filtering of training data, which that paper reports as a positive.
6. **Hold F12a's +0.811** until 284 lands (the first-HGNC-group blocking caveat).
7. **Bind the next preregistration** on any multi-line set to a depth-stratified or depth-adjusted pass rule
   (the C2a lesson).
8. **Cite as textbook, not as discovery:** F3a's eigenvalue-bias mechanism (Ledoit–Wolf; Marchenko–Pastur), F4a's
   saturation of a rank-among-decoys statistic, and F5a's two named programmes (they are Broad-16's own).

## 4. Files

- `NOVELTY-AUDIT.md` — the deliverable: 24-row table, "what can be claimed as new", "contradictions to resolve", tally.
- `work/refs.md` — 40 references in five groups with URLs/DOIs and 2026-09-06 access dates; §E lists the local checks.
- `work/checks287.py`, `work/checks287.log` — the four read-only checks with KA-1.

**Limits of this audit.** Two hours of reading; the literature on perturbation-response prediction moves weekly and I
searched it, not exhausted it. Three sources were unreachable to the fetcher (`cell.com` 403 for the CD4 journal
version and the Cell Genomics iPSC map; `nature.com` behind an IdP redirect for TxPert and Systema) and were read
through preprints, PMC or search-index summaries — those rows carry the weaker provenance. bioRxiv rate-limited three
fetches, each retried once. I did not read the Broad-16, Pan and Wang preprints in full; I extracted targeted
quantitative claims. Every number I computed myself is behind KA-1 and reproducible from `work/checks287.py`.
