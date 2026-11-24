# 360 — DepMap/CCLE covariates for the three context lines — PLAN

## Hypothesis
The three challenge contexts are DepMap/CCLE cancer lines (356: A=Jurkat, B=HeLa, C=CAL-33).
If DepMap 26Q1 (CC BY 4.0) carries basal expression + CRISPR dependency + mutations for exactly
these lines, we get the only same-line covariates available for B and C, plus an independent
check of the 041 identities (DepMap basal vs challenge control pseudobulk correlation) and a
per-target dependency table the PI can use to reason about response strength per context.

## Inputs (exact paths)
- Held: `data/external/ccle_broad_2025/ccle_cache.npz` (1600 samples × 19209 genes, TPM) +
  `ccle_cache.meta.json`; `data/external/depmap26q1/extracted/{metadata,processed}/` (16-line
  Broad-16 Cas9 extract — NOT the portal release; checked 2026-09-10).
- Challenge: `data/challenge/context_{A,B,C}.h5ad` (backed/chunked reads only),
  `data/challenge/gene_names.csv` (18,533), `data/challenge/pert_counts.csv` (300 targets).
- Donor line lists: `experiments/332-conditional-median/inputs.json` (12 Broad-16 working lines);
  CRISPRi donor names K562/HCT116/HEK293T/HepG2/RPE1/H1 from selection.json/341 core.
- Downloads (portal https://depmap.org/portal/data_page/?tab=allData, release 26Q1, CC BY 4.0):
  `Model.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (or 26Q1 equivalent),
  `CRISPRGeneEffect.csv`, `OmicsSomaticMutations.csv` if < 1 GB. Subset to 3 context + 12
  Broad-16 donor rows on the fly, keep only the subset under `work/`, delete full files.
  Record URL, bytes, sha256 of everything kept.

## Estimator / analysis
1. Inventory: exact presence/absence of the 3 lines (+ donors) in CCLE cache and depmap26q1
   extract, with sample IDs and ACH IDs.
2. `work/three_lines.npz`: per line, basal log-TPM and CRISPR gene effect on the challenge
   18,533-gene axis (symbol mapping, losses reported); same for every identifiable donor line.
3. (a) Pearson (and Spearman) of DepMap log-TPM vs challenge control mean-CPM (log) per
   context, matched-identity vs the 2 shuffled identities — identity validation.
   (b) Nearest donors per context by expression correlation and by dependency-profile
   correlation (common-gene axes, reported).
   (c) For the 300 panel targets: context-line Chronos score vs donor scores; flag targets
   with strong context dependency (score < −0.5 / prob-dependent) and donor-divergent ones.

## Controls / falsifiers
- KNOWN ANSWER first: recompute 343's mean NMAE 0.950836 over the 58 eligible H1 targets from
  323's saved `stats-cap15.npz` + `reference.npz` with my own code (read-only). If mismatch,
  stop and report.
- Identity check must beat shuffled pairings; a matched correlation that does not exceed both
  shuffles kills the "independent validation" claim.
- Negative result (lines absent / correlations flat / dependencies uninformative) is a result.

## Kill criterion
If the three lines are not identifiable in any licensed DepMap/CCLE matrix (no ACH match and
no expression match above shuffles), OR the licensed files cannot be fetched under the
constraints (< 1 GB kept, CC BY 4.0, no paid API), stop with a NEGATIVE + inventory.

## Budget
2 h. Memory gate before heavy reads; ≤12 GiB RSS; backed/chunked h5ad; delete scratch.
No builders, no submit, no git ops, no reserve paths.
