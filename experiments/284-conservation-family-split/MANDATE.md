# 284 — 281 rerun with EVERY membership of a held family excluded from training, target-cluster (component) bootstraps, and the test that separates "our features cannot predict the gene-level remainder" from "the remainder is noise". Analysis only, ≤ 3 hours. No build, no submission.

You are a research subagent. Work in `experiments/284-conservation-family-split/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`; no download (every input is on disk). `uv run python`.
Peak RSS under ~6 GB. **Time box: 3 hours.** Two sibling experiments (285, 286) run in parallel; do not read or write
their directories. `docs/RESEARCH-PROGRAMME-2026-09.md` binds you: preregister (hash) before any held number; name the
unit; strong baselines; report held numbers only as the result. Both reserves in `docs/CONFIRMATION-LEDGER.md` are SPENT:
never read `experiments/267-*/work/reserve/` or the four Broad-16 lines JVE127 / MG63 / SLR23 / UMRC3.

## Why (third external review, 2026-09-06; verified in the code)
281 (`experiments/281-what-predicts-conservation/RESULTS.md`) reported "conservation is predictable at the FAMILY level
(ρ +0.81 over 51 families) and only weakly at the gene level (held-family Spearman +0.320 [+0.256, +0.379], R² 0.111),
not better than measured reliability alone". Its family-held split blocks on each gene's FIRST HGNC group
(`feat281.py:122` `gene_group1 = gene_group_all.split("|")[0]`; `fit281.py:197`), so a gene with two or more HGNC groups
can share ANOTHER family across train and test. That weakens "whole families held out", and 281's conclusion
"conservation is a family property, not a gene property" conflates two different statements: (i) the features cannot
predict the gene-specific remainder; (ii) that remainder is measurement noise. (ii) was never tested. This experiment
fixes the split, re-runs 281's preregistered analysis under it with cluster-aware bootstraps, and tests (ii) directly:
is the gene-level remainder REPRODUCIBLE across disjoint sets of source datasets?

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; `docs/FINDINGS-2026-09.md` F8, F9, F12 (with the 13:40Z caveats); 281's
`PREREGISTRATION.md`, `RESULTS.md` §0–§6, §8, and its scripts (`feat281.py`, `fit281.py`, `boot281.py`, `core281.py`,
`fam281.py`, `rel281.py`) — reuse its functions by copying what you need into this directory (do not import from a
path that checks its own preregistration digest); its feature table `experiments/281-…/work/feat281.csv` (17,406 × 43;
`gene_group_all` carries EVERY HGNC group, pipe-separated; `gene_group1` the first) and `work/pred281.csv` (its
splits, `split_family`); 279's `experiments/279-znf-state-dependence/work/held279.npz` (keys: `genes` 17,106, `keys` /
`lines` 28 sources, `pairs_cross` 364 × 2 source-index pairs, `Ccross` 364 × 17,106 per-pair conservation cosines,
`c_held` = the held map, `n_pairs`, `null_mean`, `null_sd`) and `work/znf279.csv` (`C_held`, `c3_balanced_core`, `r_bar`,
`C_rel`); HGNC `data/external/xcell/hgnc_protein_coding.txt` (`gene_group`).

## Protocol (write `PREREGISTRATION.md`, hashed into `work/PREREGISTRATION.sha256`, BEFORE any fit on a new split)
1. **Known answers first.** Re-derive 281's primary cell on 281's own split from the copied feature table: the held
   set must equal `pred281.csv`'s `split_family` (3,422 genes, 481 blocks), ridge Spearman **+0.3204** and R² **+0.1113**,
   gbm **+0.3160**, reliability alone **+0.3517** — at 0.0e+00 (same seed 20260906, same inner CV). Abort if not.
2. **Quantify the leak before fixing it:** genes with ≥ 2 HGNC groups; under 281's split, the number of TEST genes that
   share any family with a TRAIN gene, and the number of TRAIN genes that carry a held family as a non-first group.
3. **Splits, frozen:** (S1, PRIMARY) ALL-MEMBERSHIP split — draw held families as 281 did (same seed, same permutation
   of family blocks, taken until ≥ 20 % of genes are in TEST by the rule below); TEST = every gene with ANY membership in
   a held family; TRAIN = genes with NO membership in any held family; family-less genes keep 281's cytoband blocks and
   their draw. (S2, stricter comparison) COMPONENT split — connected components of the gene–family graph (two genes
   linked when they share any HGNC group); report the component-size distribution and the largest component's share;
   hold whole components. Declare in advance which is primary (S1) and that S2 is reported next to it. 281's random
   split (seed 20260907) stays as the reference row. Every preprocessing step and hyperparameter is fitted inside TRAIN.
4. **Models, baselines, metrics exactly 281's** (ridge; gbm depth ≤ 3; mean predictor; calibrated best single feature;
   reliability alone `r_bar`; features + reliability) on `C_held`, on `C_rel`, and on the balanced core
   (`c3_balanced_core`, 8,541 genes); Spearman, R², bottom-decile AUC; the paired margins vs reliability alone.
   **Uncertainty: block bootstrap over the CONNECTED COMPONENTS of the test genes** (the cluster is the component, not
   the first family), 2,000 draws — the target-cluster bootstrap the review asked for; state how it compares with
   281's first-family block bootstrap.
5. **The family axis under S1:** grouped 5-fold out-of-fold predictions where a gene in fold k's test shares NO family
   with fold k's training; observed vs predicted family medians over the 51 families with ≥ 20 genes (family
   bootstrap); the zinc-finger deficit predicted vs observed (281 §6: −0.01506 vs −0.01785).
6. **The remainder test (the distinction the review asked for), frozen:** partition the 28 sources into two halves that
   share NO source and NO cell line (`held279.npz` `lines`; check how lines repeat across sources), keep the cross-line
   pairs inside each half (`pairs_cross`, `Ccross`), and compute per-gene C⁽¹⁾ and C⁽²⁾ as 279's statistic over the
   pairs of each half (KA: over ALL pairs the same code reproduces `c_held` at 0.0e+00). For each half, the gene-level
   remainder = C⁽ⁱ⁾ − the S1 grouped out-of-fold feature prediction fitted on that half's map (and, as a second
   remainder, minus the family-median-only prediction). Statistic: Spearman(remainder⁽¹⁾, remainder⁽²⁾) over genes
   with ≥ 1 pair in both halves, with the component block bootstrap; the raw Spearman(C⁽¹⁾, C⁽²⁾) as the ceiling; an
   attenuation-corrected value where `r_bar` exists. Null: permute remainder⁽²⁾ within components (destroys gene
   identity, keeps family structure). Repeat over 20 random admissible partitions and report the distribution.
   **Pass rule (frozen):** the remainder's cross-set Spearman interval excludes 0 and the point estimate is ≥ 0.05 →
   "the remainder is a GENE PROPERTY the features do not predict"; interval covering 0 → "not distinguishable from
   measurement noise at this depth". Report which of 281's two readings survives.
7. **Descriptive (post-hoc, labelled):** the 20 genes with the largest reproducible remainder each way (both halves
   agree), their families, complexes, essentiality.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the KA table; the leak quantification; the S1 / S2 / random tables next
to 281's; the family axis under S1; the remainder test with its rule; the descriptive tail; `work/pred284.csv` (per
gene: observed, S1/S2 predictions, split, component id, remainder⁽¹⁾/⁽²⁾ where defined); one line: **FAMILY-HELD (ALL
MEMBERSHIPS): ρ … [lo, hi], R² … (281: +0.320 / 0.111); vs reliability alone …; REMAINDER REPRODUCIBLE ACROSS DISJOINT
SOURCE SETS: YES / NO (ρ … [lo, hi]; ceiling …)**. No model use for the competition is proposed. Then stop.
