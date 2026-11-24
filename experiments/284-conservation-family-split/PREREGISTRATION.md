# 284 — PREREGISTRATION: 281 rerun with EVERY membership of a held family excluded, component bootstraps, and the remainder test

Frozen **before any model was fitted on a new split** and before any number from S1, S2 or the
remainder test was computed. The only things read beforehand were structural: the column list of
`experiments/281-what-predicts-conservation/work/feat281.csv`, the multi-membership count of
`gene_group_all` on the 17,106 held genes (0 groups 3,135 · 1 group 9,237 · ≥ 2 groups 4,734), the
key/shape list of `experiments/279-znf-state-dependence/work/held279.npz`, its 28 source keys and
21 line labels, and 281's own already-published numbers (`work/fit281.json`, `work/boot281.json`).
No association between a feature and a target on any new split has been computed.

Seed **20260906** everywhere unless stated. Every number reported as *the result* comes from a held
split. Unit of analysis and of every resampling: **the gene** (bootstrapped in blocks — see §5).

---

## 1. Question

281 concluded "conservation is a FAMILY-level property, not a gene-level property". Two defects,
both named by the third external review and verified in the code:

1. Its family-held split blocks on each gene's **first** HGNC group only
   (`feat281.py:122`, `fit281.py:197`), so a gene with two or more HGNC groups can share **another**
   family across train and test. "Whole families held out" is weaker than stated.
2. "Family property, not gene property" conflates (i) *our features cannot predict the gene-specific
   remainder* with (ii) *that remainder is measurement noise*. (ii) was never tested.

This experiment (a) fixes the split, (b) re-runs 281's preregistered analysis under it with
cluster-aware bootstraps, and (c) tests (ii) directly: **is the gene-level remainder reproducible
across disjoint sets of source datasets?**

## 2. Known answers (run first; ABORT on failure)

| row | requirement |
|---|---|
| **KA-1** | 281's primary cell re-derived on 281's OWN split from the copied feature table and copied functions: held set identical to `pred281.csv`'s `split_family` (3,422 genes, 481 blocks), ridge Spearman **+0.320420**, R² **+0.111332**, gbm **+0.315974**, reliability alone (ridge on `r_bar`) **+0.351664**, all at **0.0e+00** |
| **KA-2** | `c_held` of `held279.npz` re-derived by my own half-code over ALL 364 cross-line pairs (`np.nanmean(Ccross, axis=0)`) — max abs diff **0.0e+00**, n 17,106, median +0.018606 |
| **KA-3** | 279's balanced core read off `znf279.csv`: n 8,541, ZNF median +0.011813, non-ZNF +0.029663, Δ −0.017850 |
| **KA-4** | shuffled-training-target control on S1, 200 independent shuffles (281's corrected KA-4′ form): the observed S1 ρ must lie outside the null's spread (0/200 draws reaching it) |

## 3. The leak, quantified before it is fixed (item 2 of the mandate)

Reported as counts, not as a fit: (a) held genes with ≥ 2 HGNC groups; (b) under 281's split, the
number of TEST genes that share any HGNC family with a TRAIN gene; (c) the number of TRAIN genes that
carry one of the 481 held blocks as a **non-first** group; (d) the same two counts under S1 and S2,
which must be 0 by construction for the family part.

## 4. Splits — frozen

Let `G(g)` = the set of HGNC groups of gene `g` from `gene_group_all` (pipe-separated, HGNC
`data/external/xcell/hgnc_protein_coding.txt`), and `B0(g)` = 281's block = first group if any, else
`NOFAM:<cytoband>`.

* **S1 — ALL-MEMBERSHIP (PRIMARY).** The 2,425 blocks `B0` are permuted with
  `np.random.default_rng(20260906)` exactly as in `fit281.py:201-210` (identical permutation).
  Blocks are taken in that order, with 281's check-then-append loop, until **≥ 20 % of the 17,106
  held genes are in TEST by the rule below**:
  * `HELD_FAM` = the taken blocks that are real HGNC families; `HELD_NOFAM` = the taken
    `NOFAM:<cytoband>` blocks.
  * **TEST** = every gene with `G(g) ∩ HELD_FAM ≠ ∅`, plus every family-less gene whose
    `NOFAM:<cytoband>` block is in `HELD_NOFAM` (281's cytoband blocks and their draw kept).
  * **TRAIN** = every other gene: no membership in any held family, and (if family-less) a cytoband
    block that was not taken.
  Because each held family now pulls in more genes, S1 takes **fewer** blocks than 281's 481 to
  reach 20 %. S1's TEST ⊇ nothing is asserted in advance; the counts are the result.
* **S2 — COMPONENT (stricter comparison, reported next to S1).** Graph on the 17,106 held genes:
  an edge between two genes when they share any HGNC group; plus, for the 3,135 family-less genes
  only, an edge when they share a cytoband (this preserves 281's `NOFAM` blocking so tandem
  paralogue arrays cannot straddle). Connected components of that graph are the units. Components
  are permuted with `np.random.default_rng(20260906)` and taken with the same check-then-append loop
  until ≥ 20 % of genes are held. The component-size distribution and the largest component's share
  of genes are reported. Whole components are held.
* **REFERENCE — random-gene**, 281's seed 20260907 (`np.random.default_rng(SEED+1)` permutation,
  20 % of genes), reported unchanged next to S1/S2.
* **Every preprocessing step** (median imputation, standardisation, ridge α, GBM hyperparameters by
  5-fold inner CV) is fitted **inside TRAIN only**.

**Declared in advance: S1 is PRIMARY.** S2 is a stricter comparison; the random split is the
reference row. The verdict sentence is read off S1.

## 5. Models, baselines, metrics, uncertainty — exactly 281's

* **Models:** ridge (`Ridge`, α by 5-fold inner CV over `np.logspace(-3,4,15)`, criterion = mean
  inner-fold Spearman); gbm (`HistGradientBoostingRegressor(max_depth=3, min_samples_leaf=40,
  l2_regularization=1.0, early_stopping=False)`, inner 5-fold over `learning_rate ∈ {0.03,0.10} ×
  max_iter ∈ {200,500}`).
* **Baselines:** mean predictor; calibrated best single feature (chosen by |Spearman| inside TRAIN);
  **reliability alone** (`r_bar`); features + reliability. The same 21 frozen features as 281.
* **Targets:** `C_held` (primary), `C_rel`, and the balanced core `c3_balanced_core` (8,541 genes).
* **Metrics:** Spearman ρ, R² against the held-set mean, AUC for the bottom decile of observed
  conservation; the paired margins vs reliability alone.
* **Uncertainty — frozen:** a **block bootstrap over the CONNECTED COMPONENTS of the test genes**
  (§4's S2 graph; the cluster is the component, not the first family), 2,000 draws, resampling whole
  components with replacement. This is the target-cluster bootstrap the review asked for. It is
  reported next to 281's first-family block bootstrap SD on the same cell.

## 6. The family axis under S1

Grouped 5-fold out-of-fold predictions in which a gene in fold k's test shares **no** family with
fold k's training: folds are whole **components** (§4's graph), assigned to 5 folds greedily by
descending size to balance gene counts. On the balanced core, observed vs out-of-fold predicted
family medians over the families with ≥ 20 genes (**primary: all memberships** — a gene counts in
every family it belongs to; 281's first-group version reported next to it), family bootstrap 2,000
draws. The zinc-finger deficit predicted vs observed (281 §6: −0.01506 vs −0.01785).

## 7. The remainder test — frozen (the distinction the review asked for)

* **Partition.** The 28 sources of `held279.npz` carry 21 distinct cell lines (A549 ×3, K562 ×5,
  H1 ×2, the rest ×1); how lines repeat across sources is reported. A partition assigns whole
  **line groups** to half 1 or half 2, so the two halves share **no source and no cell line**.
  **Admissible** = each half has ≥ 8 distinct lines and ≥ 20 cross-line pairs. Draws are made with
  `np.random.default_rng(20260906)` (each line group to a half with p = ½, redrawn until admissible).
  **Partition #1 is the primary partition**; 20 admissible partitions in total are drawn and their
  distribution reported.
* **C⁽ⁱ⁾** = `np.nanmean` of `Ccross` over the cross-line pairs whose two sources are both in half i
  — 279's statistic restricted to that half (KA-2 fixes the code on all 364 pairs).
* **Remainder.** For each half: the S1 grouped out-of-fold **ridge** feature prediction fitted on
  **that half's map** (target `C⁽ⁱ⁾`, component folds of §6, α re-chosen by inner 5-fold CV inside
  each fold's training rows); `remainder⁽ⁱ⁾ = C⁽ⁱ⁾ − prediction`. **Second remainder:** minus the
  family-median-only prediction, computed out-of-fold as the mean over the gene's families of the
  median `C⁽ⁱ⁾` of that family's training genes (training global median for family-less genes).
* **Statistic:** Spearman(remainder⁽¹⁾, remainder⁽²⁾) over genes with ≥ 1 pair in **both** halves,
  with the component block bootstrap (2,000 draws).
* **Ceiling:** the raw Spearman(C⁽¹⁾, C⁽²⁾) on the same genes.
* **Attenuation correction:** split each half's pairs into two disjoint sub-halves, recompute the
  remainder in each, correlate → `ρ_ii`, Spearman–Brown it to full length (`2ρ/(1+ρ)`); report
  `ρ_obs / sqrt(ρ_11 ρ_22)`. Because sub-halves of a half share sources, `ρ_ii` is biased **up** and
  the correction is therefore conservative — stated with the number. Reported next to it, descriptively:
  the remainder Spearman by `r_bar` tercile on the genes where `r_bar` exists.
* **Null:** permute `remainder⁽²⁾` **within components** (destroys gene identity, keeps family
  structure), 200 draws; report the null mean ± SD and the share reaching the observed.
* **PASS RULE (frozen).** On the primary partition: the remainder's cross-set Spearman 95 % interval
  **excludes 0 and the point estimate is ≥ 0.05** → **"the remainder is a GENE PROPERTY the features
  do not predict"**. Interval **covering 0** → **"not distinguishable from measurement noise at this
  depth"**. If the median over the 20 partitions disagrees with partition #1 in sign, or falls on the
  other side of the 0.05 threshold, the verdict is reported as **EQUIVOCAL** with both numbers.
  The verdict states which of 281's two readings survives.

## 8. Post-hoc, labelled as such

The 20 genes with the largest reproducible remainder each way (both halves agreeing in sign), with
their HGNC families, complex membership and essentiality. Descriptive; cannot change any verdict.

## 9. Multiplicity

Cells that decide a verdict: **3** — S1 ridge and S1 gbm on `C_held`, and the remainder test's
primary partition. Everything else (S2, the random split, `C_rel`, the balanced core, the family
axis, the 20 partitions, the second remainder, the strata, the tail) is secondary, comparative or
descriptive and is labelled.

## 10. Scope

Analysis only. No build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`, no
download, no production file touched. `uv run python`. Peak RSS under ~6 GB. Nothing is written
outside `experiments/284-conservation-family-split/`. **No model use for the competition is proposed
regardless of outcome.** The reserves of `docs/CONFIRMATION-LEDGER.md` are SPENT and are not read:
neither `experiments/267-*/work/reserve/` nor the Broad-16 lines JVE127 / MG63 / SLR23 / UMRC3.
