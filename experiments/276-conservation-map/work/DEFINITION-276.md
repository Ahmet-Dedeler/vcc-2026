# 276 — the conservation statistic, FROZEN before any cosine was computed

Written 2026-09-06, experiment 276, after the inventory (`work/inventory276_*.csv`) and before
`build276.py` computed a single pairwise cosine. Nothing below may be changed after this file's
sha256 is recorded; anything added later is labelled POST-HOC in `RESULTS.md`.

## 1. The unit of a source

A **source** is one (dataset, cell line[, cell state]) pair that measured genetic perturbations of
panel targets and has a compact-donor table on disk (`delta.npy` / `b_ctrl.npy` / `genes.txt` /
`perts.txt`, the object `scripts/build_donor_table.py` writes and production uses everywhere), or
one (target, line) LINCS/CMap consensus signature from 270's extraction. Aggregates over several
lines (`jiang_6line`, `jiang_5line_noK562`, the line-collapsed CMap tables), mouse tables, CRISPRa
tables (the opposite perturbation direction) and sources carrying < 5 panel targets are EXCLUDED.

## 2. The common residual space

* **Common gene axis G**: the intersection, over the included single-cell/DE sources, of the gene
  symbols each measures, intersected with the challenge's 18,533 (`data/challenge/gene_names.csv`).
  Frozen once; |G| reported. The panel's own target-gene columns are NOT removed (the scorer
  removes them; here the target gene's own knockdown is part of the response direction — the
  variant with them removed is a declared secondary).
* **Residual**: for source s, take its rows for the panel targets it carries, restricted to G, and
  subtract the mean over exactly those rows. (Each source centred on its own across-target mean.)
* No whitening, no scaling: the statistic is a cosine, so per-row scale is irrelevant, and the
  across-target mean is the only thing removed.

## 3. The statistic

For target t:

* **`C_raw(t)`** = mean over all pairs (s, s') of sources that both carry t and come from
  **DIFFERENT cell lines** of `cos(R_s[t], R_s'[t])`. `n_pairs(t)` is recorded.
  Pairs from the same cell line (K562-GWPS vs PerturbMatch vs Jiang-K562; the three CD4 states;
  the two A549 screens) are EXCLUDED from the primary — they measure reproducibility, not
  conservation across cell state. Declared secondary `C_within(t)`: the same over same-line pairs.
* **`C_rel(t)`** = mean over the pairs where BOTH sources have a measured per-target split-half
  reliability, of `cos(R_s[t], R_s'[t]) / sqrt(r_s(t) · r_s'(t))`, with r the Spearman–Brown
  full-depth reliability `2r_half/(1 + r_half)` and pairs dropped when either `r ≤ 0.02`
  (the correction is unstable there); the number of dropped pairs is reported. Sources with
  reliability: K562-GWPS, Orion HCT116, Orion HEK293T (271's batch halves), DLD-1 (its two
  low-MOI replicates). Every other source's reliability is UNKNOWN and only `C_raw` is reported
  for its pairs.
* **`C_ext(t)`**: `C_raw` over the SC sources plus the CMap per-(target, line) signatures, whose
  cosines are computed on the pairwise gene intersection restricted to CMap's 978 measured
  landmark genes. Declared SECONDARY (a different modality; 270 measured CMap cross-line sign
  stability at chance against per-target `pds`).

## 4. Intervals and nulls

* **Interval on `C(t)`**: paired bootstrap over that target's pairs, B = 2,000, seed 20260906,
  percentile [2.5, 97.5]. Targets with `n_pairs` < 3 get no interval.
* **Interval on a distribution statistic or a Spearman**: paired bootstrap over targets, B = 2,000.
* **Null**: target-shuffled pairings — for each pair (s, s') of target t, replace source s''s row
  by that of a uniformly random OTHER panel target carried by s' (seed 20260906), keeping the pair
  structure and the gene sets exactly; 200 draws; the null mean and its [2.5, 97.5] band are
  reported per target and pooled.

## 5. Replicate validation of the SCORE (item 3)

Strict version: rebuild the map from **half A of every source that has halves and from nothing
else** — {K562-GWPS batch-A, Orion HCT116 batch-A, Orion HEK293T batch-A, DLD-1 replicate 1} — and
again from half B / replicate 2, and report the Spearman between the two per-target scores over the
targets where both are defined, with a paired target bootstrap. The two maps share no cells.
A score whose halves do not agree is labelled NOISE.
Declared secondary: the same with the full-row (half-less) sources included in both halves — an
UPPER BOUND, because those sources are shared.

## 6. Does conservation predict transfer? (item 4)

Outcomes, all pre-existing and untouched:
* `pds_t` of `cd4svd5` per (target, fold) on the three PLOLO cancer folds
  (`experiments/262-where-pds-is-lost/work/diag262-{k562_gwps,orion_hct116,orion_hek293t}.npz`,
  key `cd4svd5|pds`) and `cov|base_cos_truth` from the same files;
* the anti-aligned label (`base_cos_truth ≤ 0`) per (target, fold);
* `INCUMBENT|pds_eff` on h1-test (`experiments/273-flex-reference-instrument/work/instrument273.npz`).

Predictors: `C_raw`, `C_rel`, and the honest **LOLO** version `C_lolo(t | fold)` = the same
statistic with EVERY source from the held line removed (K562 fold: k562_gwps, perturbmatch,
jiang_K562 and CMap's K562 columns; HCT116 fold: orion_hct116; HEK293T fold: orion_hek293t and
xucao_hek293; h1-test: gse295214_h1, h1_2025).

Baselines, in the same table: 262's donor–donor agreement (`agree_donor_donor`, the two-line
version of the same idea, ρ +0.38 pooled in R262), the CD4 DE count (`cov|cd4_count`… the
number of CD4 conditions carrying the row, R262 §4.3), and the target gene's own expression in
the challenge control profile.

Tests: Spearman with a paired target bootstrap, per fold and pooled; AUC for rest-vs-anti-aligned.
No fitting of any kind; no free parameter anywhere in the statistic.

## 7. Known answers (must reproduce before anything is read)

* **KA-1**: 262's per-fold `cd4svd5` `pds` means = 0.815384 / 0.801659 / 0.780978 (R262 KA-1).
* **KA-2**: 262's per-fold ρ(base·truth cosine, `pds`) = +0.913 / +0.885 / +0.855 (R262 §4.1).
* **KA-3**: the pair cosine of the two production donor rows reproduces 262's
  `cross_line_truth_cosine` scale for K562–HCT116 / K562–HEK293T / HCT116–HEK293T
  (that table is TRUTH-level, ours is DONOR-row level; the KA is that our
  `agree_donor_donor` re-derivation reproduces R262's pooled ρ +0.381 with `pds` to ±0.02).
* **KA-4**: 271's per-target split-half reliability medians +0.047 / +0.027 / +0.028 for
  K562 / HCT116 / HEK293T (F1, R267 §3), read from the half banks used here.
