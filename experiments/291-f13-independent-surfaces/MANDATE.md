# 291 — F13 on surfaces OUTSIDE the 28-bank map: does a response's geometry, measured on the banks, predict its cross-line conservation (i) across the 12 Broad-16 working lines (disjoint lines, Cas9 knockouts) and (ii) in the GWCD4i authors' own published CD4-vs-K562 conservation (1,773 genes × 3 states)? Two preregistered claims, ≤ 3 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/291-f13-independent-surfaces/`. **Edit nothing outside it** except
appending the two claims and their results to `docs/CONFIRMATION-LEDGER.md` as instructed below. No build, no `vcc`,
no `cell-eval2`, nothing under `submissions/` or `data/`; no download (everything is on disk). `uv run python`. Peak RSS
under ~8 GB. **Time box: 3 hours.** Two sibling experiments (285, 286) run in parallel; do not read or write their
directories. `docs/RESEARCH-PROGRAMME-2026-09.md` binds you (rules 1, 2, 5, 9, 10): write the claims and their pass
rules BEFORE computing any feature–outcome number on either surface. The four reserved Broad-16 lines (JVE127, MG63,
SLR23, UMRC3) are SPENT and must not be read (`load268.py` asserts this; keep the assertion).

## Why (F13, `docs/FINDINGS-2026-09.md` v1.5; R290 §3–§9)
290 found, inside 279's 28-bank conservation map, that a gene's cross-line conservation remainder is predicted at the
ceiling by the GEOMETRY of its response measured on a disjoint half of the banks — breadth and spectral centrality
(energy in the bank's leading PCs) — with programme identity adding little, and with the sign opposite to the
preregistered mechanism: broad, spectrally central, proliferation-loaded responses are MORE conserved. Both halves were
the same 28 banks, the same pipeline and the same 2,739-gene axis (R290 §10.1). **The claim needs a surface the map
never touched.** Two exist on disk: (A) Broad-16 — 16 cancer lines (12 working), Cas9 knockouts, 100 targets (92 in
the held map), lines and dataset disjoint from all 28 banks, a different perturbation modality; (B) the GWCD4i
authors' own supplementary table `K562_comparison.suppl_table.csv` — per gene and activation state, the Pearson
correlation of DESeq2 log-fold-changes between primary CD4 T cells and K562 (`logfc_pearson_r`, with three random
controls `random_r1..3`, `donor_correlation_mean`, and MASH DE-gene counts per condition), 1,880 genes of which 1,773
are in the held map; CD4 is in none of the 28 banks. Neither surface is untouched (268 / 283 / 289 analysed Broad-16;
287's checks read the GWCD4i table against family labels) — say so; neither has ever been compared with these
features, and the claims below are frozen before that comparison.

## READ FIRST
290's `RESULTS.md` §2 (the feature formulas), §8.2–8.3, §9, and `lib290.py` / `build290.py`; **the per-source feature
store `experiments/290-…/work/persrc290.npz`** (keys `genes` 17,106, `keys` / `lines` 28 sources, `half1` / `half2`,
`on_axis`, and per-feature arrays over sources × genes: `mag`, `breadth`, `gen_cos`, `pc5`, `pc20`, `pe_*` …) — every
feature for any gene over ANY source subset is an aggregation of that store; 284's `work/pred284.csv` (the annotation
predictions, `component`, `gene_group_all`); 279's `work/znf279.csv` (`C_held`, `r_bar`, `class`); 268's `RESULTS.md`
§2 and `work/{load268.py, decomp268.py, decomp268.json, halves/<line>.npz, inventory.json}` (the 12 working lines' full
responses on 12,199 genes, cell halves, guide halves; the per-target shared / line-specific / noise decomposition); 276's
`DEFINITION-276.md` (the conservation statistic: mean pairwise cosine of residual rows across lines); 287's
`work/checks287.py` (how it read the GWCD4i table; its KA: the random controls read at chance) and
`data/external/cd4_cell2026/data_sharing_readme.md`; `data/external/broad16line/metadata/cell_line_metadata.csv`,
`data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv` (Broad-16 depth per line, rule 9).

## The claims (write `PREREGISTRATION.md`, hash it into `work/PREREGISTRATION.sha256`, append both claims to
## `docs/CONFIRMATION-LEDGER.md` as C3 and C4 in C1's format — BEFORE any feature–outcome number)
**Features (frozen; 290's 18 columns, aggregated over the stated source subset; standardised on the training genes):**
for surface A, over ALL 28 banks (no Broad-16 line is among them); for surface B, over the banks that contain NO K562
source (exclude every source whose line is K562: `k562_gwps`, `k562_essential`, `gasperini_k562`, `jiang_K562`,
`perturbmatch` if its line is K562 — read `lines` and state the excluded list), so the K562 side of the outcome shares
no data with the features; a labelled sensitivity arm includes them.

- **C3 (Broad-16).** Outcome: per target, the cross-line conservation across the 12 working lines — 276's statistic on
  268's full responses (mean pairwise cosine of the line-centred residual rows over the 66 line pairs; the same gene
  filter as 268; ALSO 268's per-target shared-energy fraction as a second outcome). Unit: targets (≈ 92 with features).
  Prediction: the 290 ridge (refit on the 17,106-gene held map with family-component folds, predicting `C_held`, then
  APPLIED to the Broad-16 targets — no Broad-16 outcome enters any fit) correlates with the Broad-16 conservation at
  Spearman > 0 with the paired target-bootstrap interval above 0; and the single features `pc20` and `breadth`
  correlate POSITIVELY (the inverted-sign finding) with intervals above 0. Baselines: 284's annotation ridge applied
  the same way; `C_held` itself (the map's own conservation — the reference: does the map predict Broad-16 at all?);
  magnitude alone. Null: 200 target shuffles. Depth (rule 9): report the per-target conservation's Spearman with the
  targets' cells per KO and the lines' UMIs per cell; nothing enters the model.
- **C4 (GWCD4i, the authors' own statistic).** Outcome: `logfc_pearson_r` per gene, in each of Rest / Stim8hr /
  Stim48hr (three cells; the frozen primary is Stim48hr, 280's run-clean condition). Unit: genes with features
  (≈ 1,773). Prediction: the 290 ridge applied to these genes (features from the K562-free banks) correlates with
  `logfc_pearson_r` at Spearman ≥ +0.15 with the component-block interval above 0 in the primary condition, beats
  the annotation ridge (paired Δ above 0), and survives partialling out the authors' own magnitude (`n_degs_MASH_K562`
  and the condition's `n_degs_MASH_*`) — partial Spearman > 0 with the interval above 0; and `pc20` / `breadth` are
  positively associated with the outcome (intervals above 0). KA (287's): the three `random_r*` columns read at chance
  (|median| < 0.01) before anything else. Baselines: annotation; `C_held` itself; magnitude (`n_degs_MASH_K562`);
  the authors' `donor_correlation_mean` (their cross-donor reproducibility, an upper-bound-like reference); the
  sensitivity arm with K562 sources included (labelled).
- **Pass rules, frozen:** C3 PASS = the 290-ridge Spearman interval above 0 AND `pc20` positive with interval above 0
  (on the mean-pairwise-cosine outcome); C4 PASS = all four conditions above in Stim48hr. Report the other cells as
  secondary. A FAIL is a FAIL; write it.

## Protocol
1. Known answers: rebuild 290's primary cell from the store (+0.317385 at 0.0e+00 on the 2,848-gene unit, the
   frozen half-1 aggregation); 279's `c_held` from `held279.npz` (a read, not a recompute); 287's GWCD4i random-control
   KA; 268's per-target decomposition read off `decomp268.json` (state the row you reproduce).
2. Fit the 290 ridge on the held map (family-component grouped CV for α; then a single fit on all 17,106 genes) with
   the surface's source subset; apply to the surface's genes; compute the outcomes; the rules; the baselines; the
   null; the depth / magnitude covariates.
3. Append the results to the ledger under C3 and C4.
4. Descriptive (post-hoc, labelled): on surface B, is the DNA-binder exception present (positive residual for
   `class` = ZNF / TF genes, 279's labels)? On surface A, which of the 92 targets are most over- and under-predicted?

## Deliverable
`PREREGISTRATION.md` (hashed), the ledger entries (claims before, results after), `RESULTS.md`: KA rows; C3 table
(outcome, prediction, baselines, single features, null, depth); C4 table per condition (with the partials and the
K562-included sensitivity); the descriptive tail; `work/pred291_broad16.csv`, `work/pred291_gwcd4i.csv`; one line:
**C3 (Broad-16, disjoint lines, Cas9): PASS / FAIL (ρ … [lo, hi]; pc20 …; annotation …; C_held …); C4 (GWCD4i
CD4-vs-K562, Stim48hr): PASS / FAIL (ρ … [lo, hi]; partial given magnitude …; annotation …; pc20 …)**. No model use for
the competition is proposed. Then stop.
