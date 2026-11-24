# Experiment 147 — dual-moment directional emitter

Status: **preregistered; source-only implementation**

User-reported start time: **2026-08-31 21:45 Europe/Istanbul**

## Hypothesis

The VCC scorer summarizes the same 400 predicted cells in two materially
different ways:

- PDS and the MSE plug-in term use pooled per-gene column totals.
- NMAE and the Wilcoxon metrics use per-cell CPM values, including their
  arithmetic mean and ranks.

Therefore a count matrix can preserve every incumbent pooled column total
exactly while changing the arithmetic per-cell CPM profile. The proposed
emitter uses two cell-depth groups and allocates each gene's fixed total
between them so that the arithmetic CPM profile moves toward the less-whitened
gamma-0.45 prediction, whose direction is the fixed source-only target.

This is not the rejected control-template emitter from Experiments 043/049/053.
That method pursued recipient-control realism and worsened NMAE/fidelity. This
experiment explicitly solves for a model-derived directional second moment.

## Frozen candidate family

- incumbent pooled-count source:
  `experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`
- directional CPM target:
  `experiments/040-union-gamma/out/runs/g0.45_q1.0_a1.08_spread0.25/pred.h5ad`
- development reference:
  `data/harness/h1-fast.h5ad`
- lambdas: `0.25`, `0.50`, `1.00`
- depth contrast: `0.50`
- 200 low-depth and 200 high-depth cells per perturbation
- within-depth quasi spread: `0.25`
- seed: `14700`, target-stable

For baseline column total `C_g`, low-group total `x_g`, average low/high
depths `r_L, r_H`, and desired arithmetic CPM `M_g`:

```text
M_g = 1e6 / n * (x_g / r_L + (C_g - x_g) / r_H)
```

The unconstrained solution for `x_g` is projected onto
`0 <= x_g <= C_g` with the exact low-group grand total. Counts are then
allocated inside each group with the frozen quasi spread and repaired only by
within-gene transfers until every requested row total is exact.

## Access boundary

The source-only stage may read:

- both prediction artifacts, including their count matrices;
- `h1-fast` observation labels and gene names only.

It must not read H1 perturbation truth. Candidate construction may later copy
only H1 `non-targeting` rows, because controls are model inputs. H1
perturbation rows may be opened only by the guarded `cell-eval2` command after
the source-only report and independent audit both pass.

## Source-only gate

Every target and lambda must satisfy:

- raw finite non-negative integer counts;
- exact row margins;
- exact incumbent column totals;
- bit-exact scorer pseudobulk;
- no cell above 1,000,000 counts;
- deterministic reproduction under the same seed.

Median relative L1 error between achieved and desired arithmetic CPM must be:

- at most `0.002` for lambda `0.25`;
- at most `0.005` for lambda `0.50`;
- at most `0.020` for lambda `1.00`.

Failure closes the experiment before truth access.

## Development decision

If the source-only gate passes, build one exact incumbent mini-artifact and
the three frozen candidate artifacts on the 40-target H1-fast panel. Run one
`cell-eval2` process at a time under a 24 GiB process-tree guard.

A candidate can advance only if:

- raw PDS is exactly equal to the incumbent;
- MSE movement does not erase the aggregate gain;
- the reference-scaled six-member average improves by at least `0.002`;
- NMAE improves;
- neither fidelity nor reach suffers a material regression.

No challenge build or submission is authorized by this experiment.
