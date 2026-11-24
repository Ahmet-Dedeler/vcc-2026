# 270 — Frozen target-level feature: cross-study SIGN STABILITY of a target's knockdown response (GPSAdb 2.0)

Written 2026-09-05 ~16:25Z, BEFORE any GPSAdb signature (DE) table was opened. The zip listing of the figshare
archive (`gpsadb_zip_listing.tsv`, member names and sizes only) and the site/paper metadata were the only GPSAdb
inputs read at freezing time. This file's sha256 is recorded in `RESULTS.md` and the file is not edited afterwards.

## Inputs
- GPSAdb 2.0 per-group differential-expression tables (figshare 29435822, `GPSAdiffResults.zip`, CC BY 4.0):
  one table per *perturbation group* g = (study, cell line, perturbed gene, method) with, per gene i, a signed
  effect `lfc_g,i` (log2 fold change treat vs control) and, if present, a test statistic / adjusted p-value.
- GPSAdb per-group metadata: perturbed gene symbol, cell line, study accession (GSE), perturbation method.
- The 300-target panel `data/challenge/pert_counts.csv` (symbols as given; no aliasing beyond exact match and
  the HGNC previous-symbol table only if a symbol is absent — any alias used is listed in RESULTS).

## Definition (primary feature, `SS_cross`)
For a panel target t with GPSAdb groups G_t:
1. Per group g: `S_g` = the top **N = 100** genes ranked by |stat_g,i| if a Wald/t statistic is present, else by
   |lfc_g,i| among genes with adjusted p available (else plain |lfc|), **excluding the target gene t itself** and
   excluding genes with `lfc = 0` or NA.
2. For every pair (g, h) of groups of t from **different studies (different GSE accessions)**:
   `U_gh` = (S_g ∪ S_h) ∩ genes measured (non-NA lfc) in both;
   `c_gh` = mean over i ∈ U_gh of 1[sign(lfc_g,i) = sign(lfc_h,i)]  (a zero on either side counts as disagreement).
   The pair is used only if |U_gh| ≥ **20**.
3. `SS_cross(t)` = the unweighted mean of `c_gh` over all usable cross-study pairs. Undefined (NA) if no usable pair.
   Chance level is 0.5 (independent random signs); 1.0 = perfect sign agreement on the union of the two top sets.
Also recorded per target: the number of groups, of distinct studies, of distinct cell lines, of usable pairs,
and the perturbation-method mix.

## Declared secondary variants (reported, NOT the primary claim)
- `SS_line`: the same statistic over pairs of groups of t from different **cell lines** (any study).
- `RHO_cross`: Spearman ρ of `lfc` over `U_gh` averaged over cross-study pairs.
- N ∈ {50, 200} sensitivity of `SS_cross`.

## Frozen tests (targets are the unit; both tests use only targets with a defined `SS_cross`)
Reference arrays (read-only): `experiments/262-where-pds-is-lost/work/diag262-<fold>.npz` for the three cancer
folds `k562_gwps`, `orion_hct116`, `orion_hek293t` (kolf2 reported separately, never pooled), with target names
from the same experiment's fold caches / arrays.
- (a) `pds_t` = mean over the cancer folds in which t is present of `cd4svd5|pds`. Test: Spearman
  ρ(`SS_cross`, `pds_t`); a **paired target bootstrap** (resample targets with replacement, B = 2,000, percentile
  95 % interval, the same resampled targets index both arrays) and a permutation null (shuffle `SS_cross` across
  targets, 2,000 draws, two-sided p).
- (b) Anti-aligned label: t is **group (e)** if `cov|base_cos_truth` ≤ 0 on at least half of the cancer folds in which
  it is present (262 §5's per-fold definition, majority-pooled). Test: AUC of `SS_cross` for (rest = 1, anti-aligned
  = 0), i.e. P(SS_cross(rest) > SS_cross(anti)) with ties at 0.5, bootstrap interval as above. A feature that
  separates conserved from state-dependent responses should read AUC > 0.5 and ρ > 0.
- Known-answer rows (must reproduce before the feature rows are read): (i) 262's pooled ρ(CD4 DE count, `pds`) and
  ρ(held-vs-donor agreement, `pds`) recomputed from the same arrays by the same within-fold-rank pooling (expected
  +0.381 and +0.630, R262 §4.1); (ii) the per-fold `cd4svd5` `pds` means 0.815384 / 0.801659 / 0.780978 (R262 KA-1).
- Multiplicity: one primary feature, two primary tests; the secondary variants are exploratory and labelled so.
- Decision rule written now: the feature "separates" only if BOTH the bootstrap interval of ρ excludes 0 AND the
  AUC interval excludes 0.5, on the targets with a defined feature. Otherwise the result is a well-supported
  negative for THIS construction (N = 100 top-set sign agreement across bulk knockdown studies), not for the idea.

## What is deliberately NOT allowed after freezing
Changing N, the pair rule, the exclusion rule, the label rule, the pooling of folds, or the set of tests after
seeing any GPSAdb DE table or any test result. Sensitivity rows are reported as such.
