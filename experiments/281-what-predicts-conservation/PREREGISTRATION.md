# 281 — PREREGISTRATION: what predicts a gene's cross-line conservation?

Frozen before a single model was fitted and before any feature–target association was computed.
The feature table `work/feat281.csv` already exists (built by `feat281.py`, which copies the target
out of 279's deliverable and computes no association with it); its four known-answer rows KA-1,
KA-2, KA-2b, KA-3 are recorded in `work/feat281.json` and all PASS. Nothing beyond feature
*coverage* has been read.

Seed **20260906** everywhere. All numbers reported as *the result* come from held splits.

---

## 1. Question and target

> From what is known about a gene **before any perturbation of it is measured** — essentiality,
> expression breadth, regulatory role, complex membership, evolutionary depth — how much of its
> cross-line response conservation can be predicted, and which features carry it?

**Target `y_raw`** = `C_held`, the cross-line conservation statistic of 279's held map
(`experiments/279-znf-state-dependence/work/znf279.csv`), 276's frozen statistic over 28 sources and
21 cell lines, on the **17,106 held (non-panel) perturbed genes**. Unit of analysis and of every
resampling: **the gene**.

**Secondary target `y_rel`** = `C_rel`, the reliability-corrected conservation on the 10,695 held
genes that have it (279 §6ii; `C_raw` divided by √(r·r′) with an r ≤ 0.02 floor).

**Second held set**: the challenge panel's **300 targets**, `C_raw` from
`experiments/276-conservation-map/work/conservation276.csv`, with per-target reference reliability
from 271's `pertarget271-<fold>.npz` `ref_rel`. The two sets are disjoint by 279's construction
(verified: overlap 0). The panel map is built on a **different axis and source set** (26 sources,
2,994-gene axis) from the held map (28 sources, 2,739-gene axis), so the panel is a *transfer* held
set, not an identically-distributed one; this is stated with its result.

## 2. Splits (declared before fitting; seed 20260906)

* **PRIMARY — family-held.** Each held gene is assigned to a block: its **first HGNC `gene_group`**
  if it has one (14,234 of 17,406 rows do), else `NOFAM:<HGNC cytoband>`. Blocks are permuted with
  `np.random.default_rng(20260906)` and taken in that order until ≥ 20 % of the 17,106 held genes
  are covered; those whole blocks are the test set. A feature therefore cannot memorise family
  identity, and tandem paralogue arrays (279 §7.1) cannot straddle the split.
* **COMPARISON — random-gene.** 20 % of genes, `np.random.default_rng(20260907)`.
* **SECOND HELD SET** — the 300 panel targets, scored by the model trained on the family-held
  split's training genes.
* Every preprocessing step (median imputation, standardisation, hyperparameter choice) is fitted
  **inside the training split only** (programme rule 3).

## 3. Features — frozen, 22 columns, all build-time observable

None is derived from a measured perturbation response. Sources and terms in `work/feat281.json`.

| block | columns | source |
|---|---|---|
| essentiality / fitness | `dm_mean_effect`, `dm_sd_effect`, `dm_frac_dependent` (effect < −0.5), `dm_frac_strong` (< −1.0), `dm_min_effect`, `dm_common_essential` (dependent in ≥ 90 % of models), `replogle_essential` | DepMap **24Q2 Public** Chronos `CRISPRGeneEffect.csv`, 1,150 models (sha256 `d15514918130…`, on disk from 269); Replogle 2022 essential library (`data/donors/k562_essential/perts.txt`, 276's convention) |
| expression | `ccle_breadth` (share of lines with log2 TPM > 1), `ccle_median_tpm`, `ccle_sd` | CCLE / DepMap 2025 mRNA TPM, 1,600 lines (`data/external/ccle_broad_2025/ccle_cache.npz`; 276's convention, KA-2) |
| regulatory role | `log_collectri_regulon_size`, `collectri_indegree`, `is_tf_collectri`, `log_dorothea_regulon_size`, `is_tf_dorothea` | CollecTRI via OmniPath (`data/external/regulon/`), DoRothEA A–C (`data/external/omnipath-dorothea/`) |
| protein complexes | `humap_n_complexes`, `humap_max_complex_size`, `reactome_n_complexes`, `reactome_max_complex_size` | hu.MAP 2.0 (Drew et al. 2021, `humap2_complexes_20200809.txt`); Reactome `ComplexParticipantsPubMedIdentifiers_human.txt` (CC0) — both fetched into `work/`, 4.2 MB total |
| evolution | `log_n_paralogues`, `ortholog_depth` (0–5: mouse, zebrafish, fly, worm, yeast) | Ensembl BioMart, current release, homologs page (fetched into `work/`, 70 MB) |

**Declared additions, reported separately, never in the primary model:**

* `fam::<family>` — one-hot of the HGNC families with ≥ 30 held genes. **Random split only**;
  meaningless in the family-held split by construction.
* `log_n_pairs`, `n_lines` — the nuisance covariate (how many source pairs measured this gene).
  **Decision, stated in advance: EXCLUDED from the primary model.** It is a property of the map's
  measurement design, not of the gene, and the question is what is knowable about the gene. It is
  reported as a declared variant (`family_raw_nuisance`) and as a per-feature Spearman.
* `r_bar` — the gene's measured split-half reliability (279 §6). **The confounder, not a feature.**
  Reported three ways: alone (the baseline to beat), added to the features, and per-feature.

## 4. Models (small, as specified) and baselines

* **ridge** — `sklearn.linear_model.Ridge` on median-imputed, standardised features; `alpha` chosen
  by 5-fold inner CV **inside the training split** over `np.logspace(-3, 4, 15)`, criterion =
  mean inner-fold Spearman.
* **gbm** — `HistGradientBoostingRegressor(max_depth=3, min_samples_leaf=40, l2_regularization=1.0,
  early_stopping=False)`; inner 5-fold CV over `learning_rate ∈ {0.03, 0.10} × max_iter ∈ {200, 500}`.
* **per-feature Spearman** on the held split, with 1,000-draw gene bootstrap intervals.

**Baselines, all on the same held genes:** (1) the **mean predictor** (training mean); (2) the
**single best feature**, chosen by |Spearman| *inside the training split* and applied sign-corrected;
(3) **reliability alone** (`r_bar`); and (4) **features + reliability**, which is how "do the
features add to the confounder?" is answered.

**Control (KA-4):** the identical pipeline with the *training* target shuffled. Its held Spearman
must have a 95 % interval containing 0 for both models, or the screen is aborted.

## 5. Metrics — held split only

Spearman ρ, R² (against the held-set mean), and AUC for the **bottom decile** of observed
conservation (the state-dependent tail; decile threshold taken on the held set). All three with
2,000-draw gene bootstrap 95 % intervals. Permutation importance (10 repeats per feature, on the
held split only, GBM) reported as ΔSpearman and ΔR².

## 6. Reporting rule (the verdict sentence, fixed now)

**CONSERVATION IS PREDICTABLE FROM GENE FEATURES:**

* **YES** iff, on the **family-held** split with `y_raw`, the better of ridge/gbm has Spearman
  ρ ≥ 0.20 with its 95 % interval entirely above 0 **and** held R² > 0 with its interval entirely
  above 0.
* **PARTLY** iff the Spearman interval is entirely above 0 but ρ < 0.20 or the R² interval touches 0.
* **NO** iff the family-held Spearman interval contains 0.

The margin over the confounder is reported as Δρ (features vs `r_bar` alone) and as
(features + `r_bar`) vs (`r_bar` alone), on the genes where `r_bar` exists.

## 7. Post-hoc, and labelled as such

The residual tail (top 20 genes each way, family/pathway concentration of the residual) is
descriptive, computed after the held numbers are read, and cannot change the verdict.

## 8. Multiplicity

Cells that decide the verdict: **2** (family-held ridge, family-held gbm, on `y_raw`). Everything
else — the random split, the family one-hot, the nuisance variant, `y_rel`, the panel, per-feature
Spearman, permutation importance, the residual tail — is secondary or descriptive and is labelled.

## 9. Scope

Analysis only. No build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`, no
production file touched. No competition arm is proposed from this experiment regardless of outcome.
