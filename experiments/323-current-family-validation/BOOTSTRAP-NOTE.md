# Paired uncertainty implementation — September 9, 2026

Use the preregistered 10,000 target draws and seed 20260909, paired across
all six arms. PDS conditions on the fixed reference gallery. The constructed
models and local H1 anchors are fixed. These intervals do not include
between-context variation, donor-training uncertainty, or anchor uncertainty.

DE and PDS aggregate the valid per-target readings, with each metric's own
reference-defined cohort. For MSE, recompute the **ratio of panel sums** and
the global across-perturbation sampling-credit cap in every resample.
No per-target MSE scores are created.

The weighted across-perturbation variance can be computed from the small
prediction row Gram matrix plus corrections for each target's excluded gene.
The synthetic regression test compares this calculation with the official
`_across_pert_budget` on explicitly duplicated rows, including a one-target-only
degenerate resample. Every arm must also reproduce its existing full-panel MSE
and the official variance budget before its intervals are accepted.

Convert aggregate raw values to H1-local scaled values, applying MSE's [0,1]
clamp and NMAE's -6 floor before taking the six-member average. Use the
same paired resample on each side of a contrast. Report wrong point signs and
intervals crossing zero separately; neither earns standalone qualification.
No threshold or candidate changes follow from these computations.
