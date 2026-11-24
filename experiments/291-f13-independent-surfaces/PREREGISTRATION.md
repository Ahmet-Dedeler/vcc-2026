# 291 — PREREGISTRATION: F13 (response geometry predicts cross-line conservation) on two surfaces
# outside the 28-bank map

Frozen before any feature–outcome number on either surface existed. Written 2026-09-06 by the 291
subagent. `work/PREREGISTRATION.sha256` carries this file's digest; every analysis script aborts on
a mismatch. Anything not in this file is POST-HOC and is labelled as such in `RESULTS.md`.

Programme rules 1, 2, 5, 9, 10 (`docs/RESEARCH-PROGRAMME-2026-09.md`) bind this experiment.

## 0. What was read before freezing (and what was NOT)

Read: `docs/READ-FIRST.md`, `AGENTS.md`, the programme, `docs/CONFIRMATION-LEDGER.md`, R290 (all),
`lib290.py`, `build290.py`, `DEFINITION-276.md`, R268 §0–§2, `load268.py`, `decomp268.py`,
`checks287.py`, the GWCD4i `data_sharing_readme.md`, and column headers / row counts of every input.

Computed before freezing (counts and availability only — programme rule: no feature–outcome
number): the 28 source keys and their lines; the K562 source list; the number of Broad-16 targets
available on all 12 working lines (97; 87 non-OR, 82 of those in the 17,106-gene held map; 10 OR
controls, 7 in the held map; `KIF11`, `PSMA1`, `RRM1` unavailable); the number of GWCD4i genes in
the table (1,880) and in the held map (1,773); the number of 284 family components spanned by each
unit (21 and 408); the per-line depth covariates. No correlation, no fit, no outcome distribution
was computed.

**Neither surface is untouched data and neither is a declared reserve.** 268 / 283 / 289 analysed
Broad-16; 287 read this exact GWCD4i supplementary table (against ZNF labels, donor reproducibility
and MASH DE counts). What is new is that **neither surface has ever been compared with 290's
response-geometry features or with any model fitted on the 28-bank conservation map**, and the two
claims below are frozen before that comparison. This is an out-of-map generalisation test, not an
untouched confirmation.

The four reserved Broad-16 lines (`JVE127`, `MG63`, `SLR23`, `UMRC3`) are SPENT and are not read;
`load268.py`'s assertion is kept.

## 1. The frozen features

290's 18-column model feature set (`lib290.FEATS290`), computed per source in
`experiments/290-remainder-response-features/work/persrc290.npz` and aggregated by the unweighted
mean over a stated source subset exactly as `build290.agg` does:

`logmag`, `breadth`, `gen_cos`, `gen_abscos`, `gen_pc5`, `gen_pc20`, the nine `pe_<HALLMARK>`
columns, `pe_prolif_abs`, `ontarget`, `on_axis`.

`gen_abscos` = mean over sources of `|gen_cos_s|`; `pe_prolif_abs` = mean over sources of
`|proj_prolif_s|`; `logmag` = `log(mean_s mag_s)`. Missing values are median-imputed by the
standardiser fitted on the training genes (290's `Standardiser`).

**Source subsets, frozen:**
* **Surface A (Broad-16): all 28 banks.** No Broad-16 cell line is among the 28, and Broad-16 is not
  one of the 28 datasets.
* **Surface B (GWCD4i): the 23 banks whose line is NOT K562**, i.e. excluding `gasperini_k562`,
  `jiang_K562`, `k562_essential`, `k562_gwps`, `perturbmatch` (its line is K562). The K562 side of
  the outcome therefore shares no data with the features. A **labelled sensitivity arm** repeats
  everything with all 28 banks.

## 2. The frozen model

`M_resp` = 290's ridge. Target: `C_held`, 279's cross-line conservation over the 28 banks on the
17,106 held genes (`held279.npz` key `c_held`; identical to `znf279.csv` column `C_held`).
Procedure: α chosen from `lib290.ALPHAS` (`np.logspace(-3, 4, 15)`) by 5-fold **grouped** CV whose
groups are 284's family-component folds (`rem284.npz` key `fold`), scored by Spearman; then a
single Ridge fit on all 17,106 genes with that α, features standardised as in `lib290.Standardiser`.
The fitted model is APPLIED to the surface's genes. **No surface outcome enters any fit.**

`M_ann` = the same procedure with 284's 21 annotation features (`lib290.FEATS`, read from
`experiments/281-what-predicts-conservation/work/feat281.csv`, `set == "held"`).

Declared secondary for both: the leak-free out-of-fold prediction (`oof_ridge` over the same
component folds), reported next to the single-fit prediction.

## 3. Claim C3 — Broad-16 (12 working lines, Cas9 knockouts, disjoint lines and dataset)

**Unit.** The Broad-16 analysis targets (the 10 olfactory-receptor controls EXCLUDED) that have a
valid `full` response on all 12 working lines and are in the 17,106-gene held map: **n = 82**.
Secondary units, declared now: the same with the 7 held-map OR controls added (n = 89), and the OR
controls alone as a negative-control row.

**Outcome 1 (primary), `C_broad16`.** DEFINITION-276 §3 applied to 268's full responses. For each
working line `l`, take its `full` responses (`load268.build_responses`, `min_cpm = 5.0`, 12,199
genes, the same gene filter as 268) for the 82 unit targets; subtract the mean over exactly those
82 rows (line-centring, DEFINITION-276 §2); unit-normalise. `C_broad16(t)` = mean over the 66
line pairs of `cos(R_l[t], R_l'[t])`.

**Outcome 2 (second outcome).** 268's per-target shared-energy fraction, column `shared_full` of
`work/per_target_decomp_cpm5.csv` (the 12-line, 12,199-gene run), on the unit targets present there.

**Prediction.** `M_resp` (features over all 28 banks) applied to the 82 unit genes.

**Statistic.** Spearman ρ; paired target bootstrap, 2,000 draws, seed 20260906, percentile
[2.5, 97.5]. Paired deltas from the same resample.

**Baselines, in the table:** `M_ann` applied the same way; `C_held` itself (does the 28-bank map
predict Broad-16 at all?); `logmag` (magnitude) alone; the mean predictor.

**Single features:** `gen_pc20` and `breadth`, aggregated over the 28 banks — the inverted-sign
finding of 290 predicts these are POSITIVELY correlated with `C_broad16`.

**Null:** 200 independent shuffles of the outcome across targets; report the null mean ± sd and the
number of draws reaching the observed ρ.

**Depth (programme rule 9), reported, nothing enters the model:** Spearman of `C_broad16` with the
target's median and minimum cells per KO over the 12 lines; and, over the 66 line pairs, Spearman of
the pair's mean cosine (averaged over unit targets) with the pair's geometric-mean
`avg_umis_per_ko` (`data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv`).
Declared secondary: the C3 ρ after partialling log median cells/KO out of both sides.

**PASS RULE (frozen).** C3 PASSES iff BOTH:
1. ρ(`M_resp` prediction, `C_broad16`) has its 95 % paired-target-bootstrap interval entirely
   above 0; AND
2. ρ(`gen_pc20`, `C_broad16`) > 0 with its 95 % interval entirely above 0.

Anything else is a FAIL and is written as a FAIL. Outcome 2 and every baseline are reported but do
not enter the rule.

## 4. Claim C4 — GWCD4i, the authors' own published CD4-vs-K562 conservation

**Unit.** The genes in `data/external/cd4_cell2026/suppl_tables_2026-08-24/K562_comparison.suppl_table.csv`
with a non-null `logfc_pearson_r` that are in the 17,106-gene held map: **n = 1,773** per condition.

**Outcome.** `logfc_pearson_r` — the authors' Pearson correlation of DESeq2 log-fold-changes between
primary CD4 T cells and K562 for that perturbation — in each of `Rest`, `Stim8hr`, `Stim48hr`.
**The frozen primary condition is `Stim48hr`** (280's run-clean condition); the other two are
secondary.

**Prediction.** `M_resp` with features from the 23 K562-free banks, applied to the 1,773 unit genes.

**Statistic.** Spearman ρ; **component-block bootstrap** over 284's `component` labels
(`rem284.npz` key `comp`), 2,000 draws, seed 20260906. Paired deltas from the same resample.

**Partial.** Partial Spearman of prediction with `logfc_pearson_r` given the authors' own magnitude
columns `n_degs_MASH_K562` and `n_degs_MASH_<condition>`: rank-transform all four vectors, residualise
prediction and outcome on the two rank covariates by OLS, and take the Pearson correlation of the
residuals; same component-block bootstrap.

**KNOWN-ANSWER, before any other number on this surface (287's KA-1):** the three `random_r*`
columns read at chance — `|median| < 0.01` in every condition (287's own bar was `< 0.02`; the
stricter bar is used here), and the real `logfc_pearson_r` median in (0.25, 0.36). Abort on failure.

**Baselines, in the table:** `M_ann`; `C_held` itself; magnitude (`n_degs_MASH_K562`, and the
condition's own `n_degs_MASH_*`); the authors' `donor_correlation_mean` (their cross-donor
reproducibility, an upper-bound-like reference, available on a subset — n reported); the sensitivity
arm with the 5 K562 banks included (labelled).

**Single features:** `gen_pc20` and `breadth` over the 23 K562-free banks.

**PASS RULE (frozen).** C4 PASSES iff ALL FOUR hold in `Stim48hr`:
1. ρ(`M_resp` prediction, `logfc_pearson_r`) ≥ **+0.15** with the component-block 95 % interval
   entirely above 0;
2. it beats `M_ann`: paired Δρ with its 95 % interval entirely above 0;
3. the partial Spearman given `n_degs_MASH_K562` + `n_degs_MASH_Stim48hr` is > 0 with its 95 %
   interval entirely above 0;
4. `gen_pc20` AND `breadth` are each positively associated with the outcome, each with its 95 %
   interval entirely above 0.

`Rest` and `Stim8hr` are reported as secondary and do not enter the rule. Anything else is a FAIL
and is written as a FAIL.

## 5. Known-answer rows (all before any claim number)

* **KA-1** — 290's primary cell rebuilt from the per-source store: aggregating the 18 columns over
  the 14 half-1 sources and running the frozen out-of-fold ridge against `remainder_half2` on 284's
  2,848-gene unit reproduces **ρ = +0.317385** (290 §8.5's own internal known-answer row) at
  **0.0e+00** to six decimals. Abort otherwise.
* **KA-2** — 279's `c_held` read (not recomputed) from `held279.npz`: n = 17,106, median
  **+0.018605**, and it equals `znf279.csv`'s `C_held` column.
* **KA-3** — 287's GWCD4i random-control row (§4 above), the gate on surface B.
* **KA-4** — 268's per-target decomposition recomputed from `work/halves/<line>.npz` on the 12
  working lines: the `shared_full`, `line_full` and `noise_full` columns of
  `work/per_target_decomp_cpm5.csv` (82 targets) reproduce at **0.0e+00**, and
  `decomp268_cpm5.json`'s `cell_halves.pooled_full_depth` reads
  shared **0.09068076312541962** / line **0.11537937074899673** / noise **0.7939398288726807** —
  the row reproduced.
* **KA-5** (internal) — the half-1 aggregation from the per-source store reproduces
  `work/feat290.csv`'s `*_h1` columns to `< 1e-6` with an identical NaN pattern (290's KA-5 bar and
  its stated reason: the frozen artifact is a float32 CSV).

Abort on any failure (`AGENTS.md`:1375-1382).

## 6. Multiplicity, and what is post-hoc

Verdict cells: **2** (C3, C4), decided by the rules in §3 and §4 only. Everything else — outcome 2,
the OR-control row, the secondary conditions, the secondary units, the sensitivity arm, the
out-of-fold variant, the depth partials, the descriptive tail (§7) — is secondary or post-hoc and
cannot change a verdict. The number of reported cells is counted in `RESULTS.md`.

## 7. Descriptive tail (declared now, POST-HOC by construction, no verdict)

* On surface B: is 290's DNA-binder exception present — do `class` = ZNF (279's `znf279.csv` label)
  and TF genes sit above the model's prediction (positive residual)?
* On surface A: which of the unit targets are most over- and under-predicted?

## 8. Operational

Analysis only. No build, no `vcc`, no `cell-eval2`, no download, nothing under `submissions/` or
`data/`, no file written outside `experiments/291-f13-independent-surfaces/` except the two claim
entries appended to `docs/CONFIRMATION-LEDGER.md`. `uv run python`. Peak RSS under ~8 GiB.
Deterministic, seed **20260906**.
