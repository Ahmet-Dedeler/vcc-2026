# 297 — Extend the novelty audit (287) to the findings that landed after it: F13 (response geometry and its sign; 290 / 291 / C5), F14 (the construct term; 285) and 293's amendment (guide efficacy; a failing construct still yields a reproducible transcriptome; the zero-shot ceiling), 284's remainder test, 286's mechanism, and the map's off-map transfer (291 / C5). Reading and citation only, ≤ 1.5 hours. No compute, no build, no submission.

You are a research subagent. Work in `experiments/297-audit-extension/`. **Edit nothing outside it.** You MAY use web
search / web fetch for papers and preprints; record every source (URL, DOI, date accessed) in `work/refs297.md`. No
`vcc`, no `cell-eval2`, nothing under `data/` or `submissions/`; `uv run python` only for reading a supplementary
table if a check needs it (state any such check and its known-answer row as 287 did). **Time box: 1.5 hours.** One
sibling (296) runs in parallel; do not read or write its directory.

## Why
287 (`experiments/287-novelty-audit/NOVELTY-AUDIT.md`, its `work/refs.md`, and the coordinator's addendum) audited
F1–F12, C1 and C2 at 14:07Z. Five results landed after it and the write-up (`experiments/294-writeup/PAPER.md`,
`CLAIMS.md`) carries three `[[MISSING: provenance verdict]]` markers and 13 `n/a (post-audit)` rows because of that.
Programme rule 10: no claim enters a write-up without its NEW / EXTENDED / KNOWN / CONTRADICTED verdict and the prior
report it extends or repeats.

## The items to audit (read the amended wording in `docs/FINDINGS-2026-09.md` v2.0 and the experiments' RESULTS.md)
1. **F13 and its narrowing** (290 → 291 → C5 in 295): a response's cross-line conservation is predicted by its own
   geometry measured in other lines; broad, spectrally central, proliferation-loaded responses are MORE conserved
   (the preregistered opposite sign refuted); programme identity adds little; the sign confirmed on Broad-16, on
   GWCD4i's CD4-vs-K562 table and on untouched reserve data; the "beats annotation" headline withdrawn off-map. Prior
   work to check: generic / shared stress components in perturbation atlases (Replogle 2022's shared clusters;
   Nadig 2025; Peidli 2024 E-distance vs "genericness"), hub / "housekeeping response" literature, the relation of
   effect size and breadth to cross-context conservation.
2. **F14** (285): cells and UMI depth as one axis obeying Spearman–Brown on real cell ladders; the guide construct as
   the dominant irreducible term (ρ_guide by dataset), batch free, a different study on the same line ρ 0.79, a
   different line 0.35–0.45, donor ≈ guide in CD4. Prior work: Replogle 2022's guide-level consistency, GWCD4i's
   `guide_correlation_*`, scPerturb / Wang 2026's reliability laws, Datlinger's CROP-seq guide replicates, the
   "cells vs depth" trade-off literature (Zhang / Svensson power analyses for scRNA-seq; Perturb-seq design papers).
3. **293's amendment**: in the challenge's own library the construct term is knockdown failure (ρ_guide +0.687 where
   both constructs work, −0.003 where the weaker fails; two-construct targets the library's weakest); a construct
   that misses its gene still produces a REPRODUCIBLE transcriptome; the zero-shot ceiling cos 0.21–0.57 vs the
   field's 0.04. Prior work: guide efficacy and off-target / integration-site / vector responses in CRISPRi screens
   (Replogle 2020 / 2022 guide validation; Horlbeck 2016 library design; the "non-targeting vs safe-targeting"
   control literature; reports that ineffective guides still produce transcriptional signatures), and any published
   estimate of a reproducibility ceiling on a perturbation benchmark.
4. **284's remainder test**: the per-gene conservation remainder after annotation reproduces across disjoint source
   sets — prior work on gene-level reproducibility of context dependence beyond family / essentiality.
5. **286's mechanism**: `pds` is paid in directions where an unbiased signal covariance finds no signal, surviving
   re-estimation of the basis from independent views (cross-validated PCA / cvPCA and its bias — Stringer 2019,
   the 2025 PNAS revision; Gavish–Donoho; Ledoit–Wolf) — is the "retrieval rewards noise directions" statement
   anticipated anywhere (e.g. whitening in embedding retrieval, hubness)?
6. **The map's off-map transfer and C5** (291 / 295): a cross-line conservation score predicting transfer into
   unseen studies and a new cell type — prior work: PerturbMap / transportability papers, GWCD4i's K562 comparison,
   any published "conservation score" used prospectively.

For each item: our claim as amended (number, unit, data); prior reports (paper, section / figure, ≤ 15-word quotes
with attribution); verdict NEW / EXTENDED / KNOWN / CONTRADICTED; what would settle a contradiction. Be adversarial:
prefer KNOWN when a paper reports it.

## Deliverable
`AUDIT-EXTENSION.md` (the table for items 1–6 in 287's format; a tally; "what the write-up can now say" for each),
`work/refs297.md`; one line: **NEW: n, EXTENDED: n, KNOWN: n, CONTRADICTED: n (list)**. No edits outside this
directory (the coordinator fills the write-up's markers). Then stop.
