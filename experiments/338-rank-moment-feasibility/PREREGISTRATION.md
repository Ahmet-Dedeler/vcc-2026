# 338 — rank-sum-constrained arithmetic-CPM feasibility

September 9, 2026. Engineering/source-prediction feasibility only. No
perturbed reference values, accuracy measurement, new reserve or submission.

## Why this is distinct

147 demonstrated two-summary emission freedom; 150 found a negative HCT116
score due mainly to reach loss. 151/152's target/gene consensus gates did
not fix it. Do not rerun those transforms or their dose selectors.
185 subsequently demonstrated rank-controlled emission at fixed library
sizes. The question here is whether useful arithmetic-CPM freedom remains
when a VARIABLE-library emitter also constrains the original rank sum.

This is not yet an emitter. Continuous per-gene histogram solutions do not
establish integer realizability, joint row-capacity feasibility, identical
tie-corrected z-scores, BH call sets, improved signs, or a competition gain.
Those are separate barriers if this relaxed bound passes.

## Inputs and fixed selection

Use only 323's cap15 predicted effects/count summaries and real CONTROL
fields. Never load real perturbed sums, jackknifes or DE values.
Four targets are selected by SHA-256 order from the 85 metadata names.
Select 32 non-target response genes: eight hash-ordered genes in each of
four control-CPM quartiles. Exclude both H1 evaluation targets and the
A/B/C panel. Exact identities are frozen in `selection.json`.

Re-emit the entire baseline to maintain its exact RNG stream. Require all
85 per-gene sums AND its whole count digest to reproduce before solving
any LP. Retain only the four selected count blocks for this experiment.
All cells remain ordinary baseline predictions; no new prediction is emitted.

## Relaxed feasibility problem, per target/gene

Choose the original 147 two-depth row-margin constructor with contrast .25,
400 cells and the baseline's exact grand count total. Its integer rounding
can produce up to four distinct depth groups. For each group j and allowed
integer count k, x(j,k) is a NONNEGATIVE CONTINUOUS histogram multiplicity.

Allowed counts: 0..K, K=min(gene total, minimum row depth,
max(64, four times the baseline maximum single-cell gene count)).
This finite-support conditional bound is not a universal mathematical ceiling.

Constraints:

- sum_k x(j,k) = number of cells in depth group j;
- sum_jk k*x(j,k) = original exact gene total;
- total multiplicity at k=0 = baseline zero-cell count;
- sum_jk Phi_control(k, depth_j)*x(j,k) = baseline sum of Phi,
  where Phi is the normalized midrank control CDF.

Phi uses the scorer's EXACT arithmetic order:
log1p(count / (library / 1e6)), not count*(1e6/library).
Minimize/maximize arithmetic CPM, sum_jk (1e6*k/depth_j)*x(j,k)/400.
Solve with HiGHS LP, two-second time limit per optimization; record all
status codes and equality/bound residuals. No outcome-based retries.
Before any outcomes, solver primal/dual tolerances are fixed at 1e-9.
The count equation is scaled by maximum support count, not total mass.
Successful solutions must have absolute gene-total residual <=1e-4,
rank-sum-Phi residual <=1e-6, all scaled equality residuals <=1e-6,
nonnegativity/reduced-cost residuals <=1e-7, and scaled primal-dual gap
<=1e-6. Save equality duals for an independent optimality audit.

Comparators: the original near-constant row depths under the same
constraints (known-answer feasible, nearly fixed arithmetic mean); and
the variable-depth problem WITHOUT its rank-sum constraint.
Preserving zero counts holds the dominant zero tie fixed, but other
pooled ties can change, so rank-sum equality is NOT z-score equality.

## Checks and declared engineering gate

Pure known answers: constant-depth fixed mean, analytically feasible
variable-depth example with fixed zeros/rank sum, measured-zero behavior,
midrank CDF versus brute pooled ranks, and equality residual enforcement.
On every selected real-control gene, baseline Phi U must match independent
pooled midranks; its tie-corrected z must match the scorer kernel.

Pass only if all baseline/identity checks pass, >=95% of the 128 pairs
solve BOTH constrained extrema, and at least 25% of ALL 32 selected genes
in EACH target have both upward and downward room >=5% of
max(baseline arithmetic CPM,1). Failed/unknown LPs count against that
fraction; never drop them to improve the denominator.

PASS means `EARNS_INTEGER_ROW_COUPLING_FEASIBILITY_ONLY`.
It does not earn response scoring, a candidate, or a submission.
Failure closes this support/depth/constraint construction, not all possible
dual-moment methods. Report every LP and comparison, not only feasible ones.

No inference of biological benefit; all six metrics remain unpriced and
potentially movable (PDS would be pinned only after an EXACT column-total
construction). Work locally, one heavy process, 8 GiB RSS ceiling, the
repository memory launch gate and a 45-minute overall budget. No paid API,
worker, cloud, production edit, model switch or reserve.
