# 328 — Joint, standard-error-aware CD4 effect estimation: source-only pilot

September 9, 2026. Codex. Register/select from metadata and commit code, selected
identities and input hashes **before reading the selected response layers**.

## Decision and distinctness

Can a small multivariate empirical-Bayes mixture learn cross-condition
effect sharing that generalizes to held-out CD4 perturbation targets?
The immediate decision is whether this estimator earns a separately frozen
current-family H1 count test. **This pilot cannot clear a package or submission.**

The source is the already-used local `GWCD4i.DE_stats.h5ad` (experiment 199):
three conditions, reported DESeq2 log2 fold changes and per-pair standard errors.
Metadata-only inspection established layer shapes 33,983 × 10,282; no selected
response values were read to design this pilot.

125/203 fitted a single per-response-gene Gaussian shared-effect model with
noise inferred from cross-family moments of row-normalized responses.
271/286 used replicate differences and gene-space Wiener/cross-product filters.
Here the mixture acts in **three-condition space**, uses each observation's
reported SE, and permits shared, condition-specific, and opposite-sign effects.
This is a MASH-style canonical covariance mixture, **not an invocation or
reproduction of the full mashr package**, and not a new statistical invention.
It adds measurement-error information not used by the current CD4 mean3
preparation. It does not add an independent biological study or recipient data.

## Fixed metadata selection

- Source: `experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad`.
- Keep target IDs with exactly one row in each of Rest, Stim8hr, Stim48hr and
  one unambiguous target name. Exclude H1's 85 and A/B/C's 300 target symbols,
  case-insensitively, **and their mapped Ensembl IDs**.
- Order remaining target IDs by SHA-256 of `328|<id>`.
  First 1,000 are fit; next 250 are held-out source-development targets.
  This is not a new reserve or independent-study confirmation.
- Use the existing 273 H1 compact bank's frozen gene map, removing H1/A/B/C
  target coordinates. Select **256 genes** by SHA-256 of `328-gene|<symbol>`.
  The same genes across all targets allow a gene-preserving target shuffle.
- Bind the exact rows, genes, columns and split in `selection.json`.
  Metadata ambiguities are excluded before any effect is read.

Read only `log_fc` and `lfcSE` for this selection. A target-gene observation
requires all three finite fold changes and finite SEs strictly greater than
zero. Missing values are not zero observations; no imputation or SE calibration
from `baseMean` (which is not a control-only expression measure).
Require at least 800 fit targets and 200 development targets with ≥64 complete
gene triples each; otherwise stop for insufficient support. Do not replace
excluded targets or genes based on response values.

## Model

For each target/gene triple `y`, sampling covariance
`D = diag(SE_Rest^2, SE_8h^2, SE_48h^2)`.
The prior is a zero-mean mixture with fixed covariance components:

1. Point mass at zero.
2. `tau^2 I`, tau in {.05, .1, .2, .4, .8, 1.6, 3.2} log2FC units.
3. `tau^2 u u'` for the 13 distinct nonzero vectors in {-1,0,1}^3,
   identifying u and -u, normalized to unit L2 length, and the same tau grid.

Thus **99 components**, with only mixture weights learned. Gaussian
likelihood `N(y; 0, D+U)`; posterior component mean `U(D+U)^-1 y`.
Fit by EM from uniform weights, up to 500 iterations, stopping after at least
10 iterations when mean log-likelihood improvement is <1e-6.
Likelihood must be nondecreasing within 1e-9 per observation.
Nonconvergence or nonfinite arithmetic is an operational NO-GO, not a claim
that joint estimation is scientifically impossible. No hyperparameter sweep.

Three fixed fits:

- **joint**: all 99 components, true aligned triples.
- **diagonal**: point mass plus seven isotropic components only. This shares
  strength but contains no cross-condition covariance or sign borrowing.
- **shuffled_joint**: same 99 components and capacity; independently permute
  the target pairing for Stim8hr/Stim48hr within each gene among fully observed
  fit targets, moving each (effect, SE) together. Rest fixed; seed 20260909.
  Marginals, gene identities, complete-case support and sample size preserved.

There is no source target identity, gene identity or held-context outcome
feature in the fitted mixture. Source-development values select nothing.

## Known answers, readouts, gates

Pure tests before data: zero component posterior; rank-one Gaussian density and
posterior versus direct 3×3 linear algebra; isotropic closed form; finite/stable
large signal arithmetic; component count and sign symmetry; EM simplex and
monotonicity; synthetic held-data signal against the zero prior; deterministic
shuffle preserving paired SE/effect marginals and support.

Before fitting, reproduce a frozen compact-bank source row's `log_fc`,
`lfcSE`, and its float32 `ln(2)*log_fc` delta from the parent H5AD on each
condition. Exact comparisons, fixed first row by compact `perts.txt` order;
no cherry-picking. Fail before interpretation if this does not reproduce.

Read source-development log predictive density on the untouched 250 targets.
Average within a target, then across eligible targets (equal target weight).
Primary differences: **joint minus diagonal**, **joint minus shuffled_joint**.
10,000 paired target bootstrap draws, seed 20260909. Keep the fixed genes,
conditions, source study and fitted models fixed. Report coverage and every
fit's convergence/weights; do not call gene triples independent replicates.

Earn a H1 count-stage registration only if all fits converge, all known
answers pass, both paired 95% lower bounds exceed zero and at least 60% of
eligible held targets favor joint over diagonal. Otherwise stop this
construction; no retuning to pass. Passing predicts observed source data
better, not unobserved biological truth or A/B/C competition performance.

## Bound and authority

One guarded heavy process, ≤4 GiB RSS and ≤2 GiB scratch, ≤45 minutes for
extraction plus the three fits. Work only in this directory. No external
response download, new reserve access, held H1 outcomes, counts, production
edits, builds, submissions, public disclosures, cloud/model API cost or push.
If earned, all-six-metric emitted-count testing requires a separate frozen
count-stage protocol and implementation before measurement.
