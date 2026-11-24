# 085 — Frozen PM inclusion gate on untouched KOLF

Date: **2026-08-23**

Status: **completed 2026-08-30; source-only training gate failed, so KOLF was not opened.**

## Question

Does one deliberately low-capacity, source-only rule identify VCC targets for
which PerturbMatch should be omitted from experiment-064's combined arm?

- comparator `D`: PerturbMatch K562-slot replacement plus post-whitening CD4;
- alternative `C`: the same query and CD4 residual, without PerturbMatch.

This is the one path experiment 066 left open: train on the HCT116 and HEK293T
experiment-064 outcomes, freeze one PM-inclusion/exclusion gate, then evaluate
that one gate exactly once against independent KOLF2.1J truth.

## Predeclared gate

The only feature is computed before any held truth is opened:

```text
raw_pm_displacement(p) = 1 - cosine(C(p), D(p))
feature(p) = within-fold percentile rank(raw_pm_displacement(p))
```

The cosine uses the complete frozen whitening-operator axis. The feature uses
only public source perturbation effects and the already-frozen arm C/D
transformations. It does not use target identity, any PDS value, any held
perturbation response, or any held control.

Percentile rank is the deterministic midrank `(average_rank - 0.5) / n`,
with exact ties receiving the same value. HCT116 and HEK293T are ranked
separately on their complete 221-target experiment-064 galleries. Production
is ranked once on the complete 248-target joint K562/PerturbMatch support set
inside the already-built 300-target VCC panel. The other 52 targets receive
an explicit feature sentinel and are routed directly to exact arm D. All of
this occurs before the KOLF artifact is authenticated.

The model family is one monotone decision stump:

```text
choose C if feature <= q, otherwise D
choose C if feature >= 1-q, otherwise D
q in {0.20, 0.30, 0.40, 0.50, 0.60}
```

`always D` is the mandatory identity fallback. Candidate selection uses only
the 221-target HCT116 and HEK293T experiment-064 outcomes. A nonidentity rule
must have nonnegative realized raw-PDS gain in both donors. Among eligible
rules, choose maximum two-donor macro gain; ties within `1e-12` prefer fewer
switches, then smaller `q`, then `low`.

The complete selection is repeated under 5,000 joint target-label
permutations, preserving the HCT116/HEK293T outcome pair for each permuted
target. KOLF may be opened only if the frozen training rule is nonidentity,
its macro gain is at least `+0.0015`, both donor gains are nonnegative, and its
one-sided full-selection permutation p-value is at most `0.05`.

No feature, model family, direction, quantile, threshold, or fallback may be
selected or changed after KOLF is authenticated.

## Exact arm contract

Before training:

1. Authenticate the experiment-064 implementation, results, per-target
   outcomes, donor banks, PerturbMatch bank, CD4 banks, and both held-donor
   whitening operators.
2. Rebuild experiment-064 arms C and D on its exact ordered 221-target
   intersection.
3. Require reconstructed per-target C and D PDS to match the pinned
   experiment-064 CSV within absolute tolerance `1e-12`.
4. Require all 300 VCC target coordinates to be absent from each whitening
   operator, CD4 residual, and scored axis.

For KOLF, build all 300 production queries before reading KOLF:

- equal K562-GWPS, HCT116, and HEK293T source weights;
- alpha `1.0796`, beta `0`;
- frozen three-donor union whitener, gamma `0.65`, lambda `0`;
- candidate-specific on-axis norm preservation;
- PerturbMatch replaces 50% of the K562 slot only when both target rows exist;
- arm C omits PerturbMatch;
- arm D includes PerturbMatch;
- both retain only complete three-condition CD4 mean residuals in `gene_norm`
  mode at fraction `0.05`;
- all 300 VCC target coordinates are excluded from whitening and CD4 support.

The complete arm C matrix, arm D matrix, gate assignment vector, and gated
matrix are serialized and hashed before KOLF truth is authenticated. Every
sealed file is then uploaded with no-overwrite semantics to a fresh run-specific
S3 prefix and its SHA-256 metadata is read back. Only after that durable freeze
passes may the runner atomically create the single global experiment-085 KOLF
access marker. The KOLF scorer uses all 273 shared targets against the complete
273-target gallery and all 17,806 shared non-panel genes.

The source-derived bootstrap clusters are also frozen before KOLF is opened.
They are reciprocal-5-nearest-neighbor connected components of the complete
300-row production arm-D query on the production whitening-operator axis,
requiring pairwise cosine at least `0.80`. KOLF rows and values never enter
the graph. A valid interval requires at least eight represented clusters
after the frozen labels are restricted to the 273 KOLF targets.

KOLF's `b_ctrl.npy` is not a biological control profile: its pinned metadata
describes it as an all-ones compatibility placeholder because the source
release has no basal profile. It is therefore forbidden as a feature.

## One-time KOLF promotion gate

Promotion requires every condition:

- complete 273-target gallery and exact 17,806-gene score axis;
- raw PDS delta over fixed arm D at least `+0.0015`;
- 20,000-replicate source-query-cluster bootstrap lower 95% bound above zero;
- one-sided 20,000-replicate matched-routing permutation p-value at most
  `0.05`;
- at least one C and one D assignment among KOLF targets with PerturbMatch
  support;
- exact identity for every target without joint K562/PerturbMatch support;
- all provenance, freeze-order, target-coordinate, and arm-reproduction checks
  pass.

The matched-routing null uniformly reallocates the frozen number of C routes
among the PM-supported KOLF targets. It preserves the C/D route count, leaves
every PM-unsupported target on exact D, and recomputes the complete-gallery
mean from the fixed per-target C and D PDS values. The add-one one-sided
p-value is `(1 + # null >= observed) / (R + 1)`.

This is effect-level evidence only. It does not authorize count generation,
`cell-eval2`, challenge packaging, or submission.

## Execution safety

Synthetic tests must run under a 1 GiB process-tree RSS guard and cannot open
real donor matrices, whitening matrices, or KOLF.

The real path is cloud-only. Its launcher must require explicit clearance,
an already-running SSM-online worker with at least 64 GiB RAM and 48 GiB
available, an uncontended global VCC lock, a new S3 run prefix, pinned source
and dependency hashes, and fixed package versions. It may not create/start an
instance or submit to VCC.

The launcher also refuses if any prior experiment-085 KOLF access marker
exists. Immediately before the first KOLF file is authenticated, the runner
durably persists the full frozen bundle, writes a local access marker containing
that freeze record, atomically creates a global no-overwrite S3 marker, and
reads the marker back to verify its SHA-256. Competing workers therefore cannot
both cross the KOLF boundary. KOLF remains unopened if either durable step
fails. A failed or completed KOLF-opening attempt consumes the one-time
evaluation.

KOLF is hypothesis-specific held truth here, not globally unseen project data:
earlier disjoint experiments used it for direct-donor and local-centering
questions. Experiment 066 predeclared this narrow PM-exclusion follow-up before
this implementation, and experiment 085 never reads KOLF outcomes while
choosing its feature, rule family, threshold, or fallback.
