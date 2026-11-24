# C7 — observed cross-study overexpression/knockdown alignment

**Coordinator final protocol, September 7, 2026.** This supersedes
PREREGISTRATION-DRAFT.md. It is written without observing R5 expression.
Freeze means: this file and complete tested analysis/access code are committed,
their exact hashes are entered in the C7 ledger record, and a separately
issued coordinator release binds those hashes. This file alone is not a release.

## Question and scope

For the 100 targets shared by the Nakatake overexpression row subset and the
Pan2025 H1 knockdown bank, is the **observed downstream response alignment**
more antisymmetric than chance target matching? Do previously measured
conservation and KD geometry predict its sign or absolute size beyond generic
alignment to unrelated targets?

This is a cross-study comparison, treated as cross-line (SEES3 versus H1).
Exact instrument provenance of the producer's Nakatake input table remains
unresolved. Different line/state, duration, intervention, normalization and
measurement error are confounded. The processed OE table has one summary
per target and a shared Emerald control, with no replicate structure.
No latent-cosine bound, disattenuation, equivalence/orthogonality claim,
unlimited bulk reliability, or causal modality comparison is available.

## Inputs and untouched subset

- R5: `data/external/gse295214/pseudo_bulk_perturbation_database.h5ad`,
  source SHA256
  `0754d6206fe045744dfd9b21eb37c8e3435152c76b158b9810a84c64a029c0ac`.
  Only `obs.dataset == Nakatake2020_hESC_TF_screen`: 510 ordered rows,
  184–693. Other unused datasets are not authorized.
- KD: `experiments/115-gse295214-h1/work/compact-donor/`.
- Frozen cohort: `work/COHORT.tsv`, 100 symbols in lexical order, with
  exact source/KD row mappings. Symbol-list hash fca12254… is distinct
  from the actual TSV file hash; the release binds both content and mapping.
- Shared measured axis: `work/AXIS.txt`, 15,309 canonical genes, SHA256
  `ec1a5719fc5a0bb624e15c13d1a3787a4d89d9c59af31c3dc5e447ef7ecdda20`.
  Every coordinate is measured on both sides for every cohort target,
  including the duplicate-collapse AND-mask rule. Padded zeros are not values.
- Conservation: `held279.npz`, average finite Ccross over the 338
  cross-line pairs **excluding gse295214_h1 on both sides**; 97/100
  cohort targets have a finite value. The source-sharing C_held is only a
  descriptive comparator.
- Geometry: gse295214_h1's `pc20` column in `persrc290.npz`.
- Spent 302 saved arrays/table are descriptive context only. All 58 overlapping
  scored targets and the separate 45 primary-subset targets are reported.
  No raw Joung atlas re-read.

The final manifest binds all files actually used, including bank axes/masks,
feature inputs, metadata/known answers, complete code and tests. Metadata and
measurement-mask observations preceded this protocol; no outcome did.

## Fixed statistic

Reconstruct source profiles by the independently known-answer-checked
115 transform:

```
perturbed_log = control_log + X
source_cpm   = expm1(max(log_profile, 0))
canonical    = sum duplicate source-symbol counts in linear space
delta       = log1p(5 * canonical_perturbed) - log1p(5 * canonical_control)
```

Stored log profiles are normalized to CPM10k; 5 maps to CP50k. Each source
row's reconstructed CPM sum must be 10,000 ± 1.

Primary OE and KD residuals subtract their respective means over the **same
100 cohort targets**. Then remove **all 98 cohort target-transcript columns**
present on the axis from both matrices. The primary scoring axis is therefore
15,211 columns for **every matched or shuffled pair**. S_i is the cosine of
the matched residual vectors on that fixed axis.

This is intentionally different from the draft's variable one/two-column
exclusion. The latter is a sensitivity. Primary inference concerns
target-specific downstream residuals; a shared response removed by centering
is not tested by it.

## Reference distribution and four-cell family

Seed 20260907. Draw **2,000 uniform permutations** of the OE target assignment
from the full permutation group, fixed points allowed. KD rows, target
features, centering and the scoring axis stay fixed. Recompute every statistic
from the permuted matching. Report plus-one Monte Carlo tail probabilities:
`(1 + number at least as extreme) / 2001`.

| cell | statistic and predeclared direction | additional sign condition |
|---|---|---|
| C7a | median S below null **and** tail excess `share(S < -.05) - share(S > +.05)` above null; cell p = max of both one-sided p values | observed median < 0 and observed tail excess > 0 |
| C7b | Spearman(C_excl, -S) above its matching-null distribution | observed rho > 0 |
| C7c | Spearman(gen_pc20, -S) above its matching-null distribution | observed rho > 0 |
| C7d | Spearman(gen_pc20, abs(S)) above its matching-null distribution | observed rho > 0 |

Holm step-down controls the **four-cell** family at .05. A PASS requires
both its adjusted p ≤ .05 and its sign condition. Unestimable cells enter
the full family at p=1 and cannot pass; report their reason and finite n.
No alternative cutoff, selected tail, dose, cohort or axis substitutes for
the primary result.

These are conditional target-matching tests under an exchangeability
assumption, not tests that randomize biological interventions. Correlated
gene families and the single study per arm limit generalization. No
family-held, study-replicated or causal conclusion is claimed.

## Intervals and descriptive secondaries

20,000 paired target-bootstrap draws give unadjusted 95% intervals for the
cell statistics. These are **conditional on frozen centered rows**; they do
not re-estimate preprocessing or provide between-study uncertainty and do
not decide PASS. Both null and observed distributions must be shown.

Always report the following as descriptive, never as replacement verdicts:

1. Uncentered control-referenced effects on the same fixed scoring axis.
2. Bank-wide centering (OE over 510 rows; KD over 2,354), same scoring axis.
3. Draft own-coordinate exclusion: remove both pairing members' own genes,
   one column for matched pairs, usually two for unmatched pairs.
4. Own-transcript change and share positive on OE / negative on KD on
   measured own coordinates. These are intervention checks, not score axes.
5. S association with 302's S_g on all 58 overlapping scored targets and
   separately its 45 primary-subset targets, with sample counts and limitations.
6. Source-sharing conservation and magnitudes as descriptive comparators.
   Any additional plots/tails are explicitly exploratory.

The draft's HGNC block bootstrap and direct OE–OE reliability interpretation
are **not adopted**. The design cannot identify biological replicate
reliability across these different studies. No unimplemented optional
secondary is grounds for suppressing or rerunning the primary outcome.

## Regression checks, validity and stop rules

Before release:

- Reproduce the 20 Pan known-answer delta/control/mask rows and 2,354-row
  identity; reproduce 279/302 reference arithmetic; no R5 expression needed.
- Synthetic exact mirror, alignment, null, direction-sign refusal,
  four-cell multiplicity/missingness and vectorized-vs-literal checks pass.
- Access refusal tests cover stale hash sidecars, live-source/metadata/cohort
  changes, wrong types/missing release, direct function invocation, and
  interrupted/second access.
- Final source/row/cohort/axis/code/input hashes and manifest-bound test
  receipts reproduce before the atomic access marker is created.

At access and analysis, shape/order, masks, finite values, product hash and
CPM10k premise are checked. Invalid data or a failed premise is a validity
failure: **no biological PASS/FAIL is inferred from it**. Record VOID or
ABORT and the reason, retaining whatever sufficient arrays were safely saved.
An unestimable statistical cell is explicitly unestimable, never positive.

One guarded source access; the marker is atomically created and fsynced
before the first expression slice. Its presence spends R5 even if interrupted.
No second source read to repair a failure. All subsequent computations use
the hashed sufficient product `work/nak303.npz`.

## Deliverables and decision

Commit the protocol and full implementation before issuing the release.
Record access time, source/product hashes, validity, all four decisions,
nulls/intervals and deviations in `RESULTS.md` and the confirmation ledger,
regardless of outcome. This is a scientific measurement, not a production
change, build approval or leaderboard submission.
