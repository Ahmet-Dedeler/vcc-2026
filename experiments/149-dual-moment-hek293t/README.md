# Experiment 149 — held-HEK293T dual-moment confirmation

Status: **preregistered; source-only build pending**

Date: **2026-08-31**

## Question

Does Experiment 148's fixed lambda-1 dual-moment emitter improve the current
a15-shaped prediction in a genuinely held biological context?

The primary surface is the clean, depth-balanced Orion HEK293T harness:

- 132 perturbations;
- 150 real cells per perturbation;
- 18,400 control cells;
- control/perturbation median-depth ratio approximately 1.0.

HEK293T perturbation expression is forbidden during candidate construction.
The build may read only the harness's 18,400 contiguous non-targeting rows
(rows 19,800 through 38,199), plus frozen target/gene metadata, K562-GWPS and
Orion-HCT116 source donors, PerturbMatch, CD4-T-cell residuals, and the
held-HEK source-only covariance operator.

## Frozen construction

- incumbent effect: current a15-shaped held-context analogue;
- source donors: K562-GWPS and Orion-HCT116 only;
- PerturbMatch: fixed 50:50 replacement of the K562 slot where supported;
- whitener: essential-pruned held-HEK operator;
- incumbent gamma: `0.65`;
- arithmetic-CPM direction gamma: `0.45`;
- alpha: `0.982`;
- CD4 residual: fixed solved fraction `0.09472706196346517`;
- profile gate: control CPM `>10`;
- incumbent count emission: current a15 block-spread policy;
- dual-moment lambda: `1.00`;
- depth contrast: `0.50`;
- within-group quasi spread: `0.25`;
- dual-moment seed: `14700`;
- 400 predicted cells per perturbation.

The candidate preserves every incumbent per-target gene total exactly, so its
PDS and pooled profile must be identical. Only the within-target allocation of
those counts changes.

## Source-only bridge checks

Before any HEK perturbation truth access, the builder must:

1. reproduce the frozen Experiment 104 HEK baseline using the unpruned
   gamma-0.65 operator and solved CD4 fraction;
2. reproduce the frozen Experiment 106 HEK baseline using the pruned
   gamma-0.65 operator and CD4 fraction `0.05`;
3. preserve every incumbent column total and scorer pseudobulk exactly in the
   final dual-moment candidate;
4. read expression only from rows whose preflighted metadata label is
   `non-targeting`; no perturbation row may be materialized.

## Frozen confirmation gate

On all 132 targets:

- every per-target raw PDS value must be exactly equal;
- official-anchor projected six-member delta must be at least `+0.004`;
- raw NMAE delta must be at most `-0.004`;
- raw reach delta must be at least `+0.010`;
- raw fidelity delta must be at least `-0.003`;
- raw Jaccard delta must be at least `-0.010`;
- raw normalized MSE delta must be at most `+0.010`;
- the target-bootstrap 95% lower bound for projected average delta must exceed
  zero.

Passing authorizes a production challenge build and local validation only. It
does not authorize a VCC submission.
