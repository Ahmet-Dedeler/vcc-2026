# 297 — RESULTS

**Mandate:** extend 287's novelty audit to the five results that landed after 14:07Z — F13 (290 / 291 / C5b),
F14 (285), 293's amendment, 284's remainder test, 286's mechanism, and the map's off-map transfer (291 / C5a).
Reading and citation only. Started 2026-09-06 16:47Z, finished 17:05Z — 18 min of the 1.5 h box.

**Deliverables (this directory only; nothing outside it was written):**
`AUDIT-EXTENSION.md` (15-row audit table in 287's format, per-item "what the write-up can now say", tensions, tally)
and `work/refs297.md` (every source with URL/DOI/route/access date, the provenance marks, and KA-297).

## Headline

**NEW: 3, EXTENDED: 7, KNOWN: 5, CONTRADICTED: 0 (none).**

| verdict | rows |
|---|---|
| **NEW (3)** | 1c F13's identity-vs-geometry ablation (hallmark sets → random matched sets costs +0.023 of a +0.294 model) · 2e a different STUDY on the same line (ρ 0.79) agrees better than a second CONSTRUCT (0.09–0.41) · 4 the 284 remainder test (ρ +0.328 [+0.241, +0.373] vs ceiling +0.382, null +0.044 ± 0.018, 0/200, 20/20 partitions) |
| **EXTENDED (7)** | 1b 290's source-disjoint out-of-fold design · 2a cells and UMI depth as one axis for a pseudobulk statistic · 2d batch free (ρ 0.87–1.12; Orion HEK293T 0.42–0.47) · 3b a failing construct still yields a reproducible transcriptome · 3c the construct-level zero-shot ceiling (cos 0.213 … ≈ 0.57) · 5 286's unbiased-estimator control on F3b · 6 the map's off-map / reserve-R3 transfer (ρ +0.502, file-block [+0.153, +0.566]) |
| **KNOWN (5)** | 1a F13's SIGN (generic/broad/proliferation-loaded → more conserved) · 2b Spearman–Brown as the reliability law in cells · 2c the guide construct as the dominant irreducible term · 2f the line 0.35–0.45 / state 0.49–0.53 / donor 0.02–0.07 ladder · 3a knockdown failure as the mechanism of guide disagreement |
| **CONTRADICTED (0)** | none found; three new tensions logged (T6 cells-vs-depth objective, T7 the "one truth in five" extrapolation, T8 the 0.79-vs-0.50 statistic mismatch) |

## The three findings that changed the reading

1. **F13's surviving sign is published.** Wang et al. 2026 label every perturbation by φ, *"alignment with the
   transcriptional axis common"* across perturbations, and report *"in 39 of 41 method–benchmark combinations, mean
   model performance followed Shared > Specific > Unreliable"* — 7,170 perturbations, 29 datasets. F13's direction is
   theirs; 290's design and target are ours.
2. **The breadth half of F13 is in a table 287 already read.** In the GWCD4i K562-comparison table the published
   conservation column tracks the MASH DE-gene count at ρ **+0.507 / +0.530 / +0.543** — a stronger association than
   290's fitted model reaches on its own target.
3. **F14's CD4 guide term is the source paper's own column.** `guide_correlation_all` medians +0.0617 / +0.0715 /
   +0.0607 against our ρ_guide 0.02–0.06; and the mechanism 293 measures is the field's standing explanation
   (Mixscale weights cells by knockdown; a 2025 re-analysis finds ≈ 56/377 K562 sgRNAs efficient).

## Known-answer row (KA-297)

No computation was run. The guard is that the three numbers this audit re-uses from 287's read-only checks still read
out of `experiments/287-novelty-audit/work/checks287.log` exactly as recorded: KA-1 real +0.3030 / +0.3147 / +0.3117
vs random +0.0040 / +0.0037 / +0.0038 → **PASS**; ρ(conservation, `n_degs_MASH`) +0.507 / +0.530 / +0.543 → **PASS**;
`guide_correlation_all` +0.0617 / +0.0715 / +0.0607 → **PASS**. Sources are marked **[F]** (paper text fetched),
**[S]** (search index only, unverified) or **[287]** (inherited); no verdict rests on an [S] number alone.

## Limits of this audit

Reading only, 18 minutes, one session, ~20 searches and 8 full-text fetches. Three sources refused retrieval (bioRxiv 429 / binary PDF) and are used at
search-index strength: `ZHANG25-CRISPRi`'s sgRNA-efficacy counts, `MIXSCALE`, `SEED24`, `DUALSG`, `TXPERT-CEIL`,
`STRAIGHT26`, `JOLLIFFE82`, `DEPTH-TSE`, `SVENSSON-DEPTH`, `SCPOWER`, `URBUT19`, `DEPPRIOR26`. Absence of a prior
report is evidence from ~20 searches, not proof. The experiments' own RESULTS.md files were audited through the
consolidated, already-amended wording in `docs/FINDINGS-2026-09.md` v2.0 and `docs/READ-FIRST.md`, which is the
wording the mandate names. Sibling 296 was not read. The coordinator fills the write-up's
`[[MISSING: provenance verdict]]` markers; this experiment wrote nothing outside its own directory.
