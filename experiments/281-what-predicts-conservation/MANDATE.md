# 281 — What predicts a gene's cross-line conservation? Gene-level features on disk against the 17,106-gene map, held-out and preregistered. Analysis only, ≤ 4 hours. No build, no submission.

You are a research subagent. Work in `experiments/281-what-predicts-conservation/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/` (fetch only small public gene-annotation tables into
`work/`, ≤ 500 MB, cite sources and terms). `uv run python`. Peak RSS under ~6 GB. **Time box: 4 hours.**
`docs/RESEARCH-PROGRAMME-2026-09.md` binds you: preregister the model and the held split before fitting; strong
baselines; report the held-set numbers only as the result.

## Why (F8, F9, F10 in `docs/FINDINGS-2026-09.md`)
279 built a 17,106-gene cross-line conservation map (`experiments/279-znf-state-dependence/work/znf279.csv`, 276's frozen
statistic over 28 sources and 21 lines) and found a family axis from core machinery (spliceosome, proteasome, ribosome
biogenesis, most conserved) to cell-identity families (receptors, homeoboxes, zinc fingers, least). 280 found that
reproducibility across donors tracks response magnitude, not cells. **The organising question's next step is
predictive: from what is known about a gene BEFORE any perturbation is measured — essentiality, expression breadth,
regulon size, complex membership, evolutionary age, TF class — how much of its conservation can be predicted, and
which features carry it?** A predictable conservation is a prior a model could carry; an unpredictable one says the
map itself is the necessary object.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; `docs/FINDINGS-2026-09.md` F8–F10; 279's `RESULTS.md` §1.1 (the held map's
construction), §4 (family nulls), §6 (the reliability terciles — reliability is a confounder to model, not a feature to
predict from), `work/znf279.csv` (columns: the conservation statistic per gene, its pair count, reliability where
available, family labels); 276's `work/conservation276.csv` (`ccle_breadth`, `ccle_median_tpm`, `essential`,
`shared_programme_share` — the feature conventions); `data/external/depmap26q1/extracted/` (gene effect / essentiality
across lines; expression); `data/external/omnipath-dorothea/`, `data/external/regulon/` (TF–target regulons: regulon
size as a feature); `data/external/ccle_broad_2025/` (expression breadth); HGNC family tables (as 279 used them);
public gene-age / paralogue tables only if small and citable.

## Protocol (write `PREREGISTRATION.md`, hashed, before any fit)
1. **Target and units:** y = the map's conservation statistic per gene (raw; and the reliability-corrected version on
   the genes that have it) on the 17,106 held genes of 279 (the panel's 300 are kept apart as a second held set, with
   their reliability from 267/271). Unit: genes. **Held split: by HGNC family** (whole families held out, ~20 % of
   genes; declare the seed) so a feature cannot memorise family identity; also a random-gene split for comparison.
2. **Features, frozen (all build-time-observable, none derived from perturbation responses):** DepMap common-
   essentiality and mean gene effect; expression breadth and median TPM across CCLE lines; regulon size (DoRothEA /
   the regulon table) and whether the gene is a TF; HGNC family class (one-hot of the ~90 families with ≥ 30 genes,
   used only in the random split, excluded in the family-held split); protein-complex membership count if a complex
   table is on disk or small to fetch; number of paralogues if available; number of lines/sources measured (a
   nuisance covariate: state whether you include it).
3. **Models, small:** ridge; a shallow gradient-boosted regressor (depth ≤ 3); a Spearman per feature. Inner CV for
   hyperparameters inside the training split only. **Baselines:** mean predictor; the single best feature; the
   reliability alone (the confounder — report how much of the variance it explains and whether features add to it).
4. **Metrics on the held split:** Spearman and R² between predicted and observed conservation; AUC for the bottom
   decile (the state-dependent tail); the same for the panel's 300 targets as a second held set; feature attribution
   (permutation importance) on the held split only.
5. **Descriptive (post-hoc, labelled):** the residual — genes whose conservation is far from predicted; are they
   concentrated in particular families or pathways? Name the top 20 each way.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the feature table with sources; the held-split results with intervals and
baselines; the panel as a second held set; the attribution table; the residual tail; one line: **CONSERVATION IS
PREDICTABLE FROM GENE FEATURES: YES / PARTLY / NO (held-family Spearman …, R² …, over reliability alone by …)**. Write
`work/pred281.csv` (per gene: observed, predicted, split). No model use for the competition is proposed here. Then stop.
