# 104 — context-held-out direct-credit adapter

Status: **complete — NO-GO under the preregistered held-surface gates**.

Date: 2026-08-30.

## Question

Can a low-capacity model trained directly on signed downstream credit use
source-donor agreement, STRING relationships, and the supplied target-context
control profile to redistribute an existing
a15-shaped effect toward more reliable genes without sacrificing
complete-gallery PDS?

The competition hides the target **context**, not the perturbation identities.
The primary design is therefore context-held-out but target-visible: a target
may have labels in another public source context, but the held context is
absent from every label, feature, statistic, operator fit, model fit, and
selection decision.

## Frozen primary intervention

The intervention is a bounded, sign-preserving reweighting of the held
context's frozen baseline effect. It cannot introduce a nonzero coordinate or
flip a sign.

For one target row `h` and fitted score `s`:

```text
z_g   = clip((s_g - median(s_active)) / max(MAD(s_active), 1e-6), -2, 2)
m_g   = exp(beta * z_g)
h'_g  = h_g * m_g
h'    = h' * ||h_active||_2 / ||h'_active||_2
```

`active` excludes all 300 VCC panel-gene coordinates. Coordinates unsupported
by the feature builder receive multiplier one. Exact zero coordinates remain
zero. The primary `beta` is **0.10**; `beta=0` is identity.

The estimator is weighted ridge regression with fixed penalty **1.0**. Its
label is signed correctness, `sign(reference_lfc) * sign(base_lfc)`, and its
sample weight is the exact marginal signed-credit weight
`abs(base_lfc) / sum(abs(reference_lfc))`, with equal total weight per target.
This targets the linear NMAE functional directly rather than binary accuracy
or a rare wrong-sign tail.

Two matched models are frozen:

1. **agreement-only ablation** — the seven donor delta/agreement features used
   by experiment 100;
2. **full primary model** — the same seven features plus target-control
   expression and fixed rank-8 STRING features.

The full model must beat the agreement-only arm; an identity-relative gain
alone is insufficient.

## A15-shaped held baselines

The challenge inference baseline is the exact staged a15 direction:

- `alpha=0.982`;
- equal real-donor weights;
- 50:50 PerturbMatch replacement of the K562 slot where supported;
- gamma-0.65 union whitening;
- post-whitening three-condition CD4 `gene_norm` residual using the exact
  a15 solved fraction `0.09472706196346517`;
- all 300 VCC target coordinates zero.

For held contexts, use the closest leakage-safe analogue:

| held context | real source donors | operator |
| --- | --- | --- |
| H1 | K562-GWPS, Orion-HCT116, Orion-HEK293T | exact production pruned three-donor operator, projected to the H1 gene axis |
| HCT116 | K562-GWPS, Orion-HEK293T | frozen held-HCT116 operator |
| HEK293T | K562-GWPS, Orion-HCT116 | frozen held-HEK293T operator |

PerturbMatch and CD4 are external fixed sources and remain available in every
fold. The held real donor is never used in its own baseline. For targets or
genes unsupported by PerturbMatch/CD4, the incumbent real-donor contribution
is retained exactly.

## Features available at submission time

Agreement-only features, per `(target, response gene)`:

1. `log1p(abs(source mean))`;
2. `abs(source mean) / mean(abs(source values))`;
3. source-sign unanimity;
4. minimum source magnitude divided by mean source magnitude;
5. source standard deviation divided by mean source magnitude;
6. number of carrying real donors;
7. mean pairwise donor cosine for the target.

The full model additionally receives:

1. sign of the baseline effect;
2. target-context control `log1p(CPM)`;
3. STRING target/response cosine;
4. eight target STRING PCs;
5. eight response-gene STRING PCs;
6. eight coordinatewise target × response STRING-PC products;
7. STRING missingness flags.

PerturbMatch and CD4 affect the a15-shaped baseline but are not model
features. This keeps the adapter low-capacity and tests the narrow positive
residue from experiment 100 without adding an unvalidated feature family. No
leaderboard value, challenge perturbation truth, hidden gallery statistic, or
KOLF artifact is an input.

## Leakage boundary and whole-target folds

The outer surfaces are H1, clean HCT116, and clean HEK293T.

For outer surface `d`:

1. `d` is absent from all training labels and all source features.
2. A labeled training context is predicted without its own donor and without
   `d`.
3. The outer baseline uses no held donor.
4. All pair rows for one perturbation target remain in one deterministic
   SHA-256 fold. Pair-level splitting is forbidden.
5. Four target folds are used to produce cross-fitted training diagnostics.
   The final outer model is then refit on both non-held contexts. Retaining
   the same target in a different context is intentional and
   production-valid.
6. The 300 VCC panel coordinates are excluded from model features,
   transformations, and PDS.
7. Directional credit excludes only the row's own target coordinate, matching
   the metric. It does **not** remove all 300 coordinates. Historical genes
   absent from the 18,533-gene challenge axis are unavailable rather than
   treated as predicted zeros and are excluded from surrogate scoring.
8. Outer perturbation truth is read only by `score` after the baseline,
   fitted models, candidate matrices, null models, hashes, and manifest have
   been frozen.

Any failed invariant aborts.

## Nulls and uncertainty

- **20 response-label permutation nulls.** Within every complete target block,
  signed labels are independently permuted within sample-weight quintiles.
  This preserves the target, gate size, class balance, and credit-weight
  distribution while destroying response-gene correspondence. Every null
  refits both agreement-only and full models.
- **2,000 target-block bootstrap draws** per surface and 20,000 macro draws.
  Pair-level bootstrap is forbidden.
- Report target-cluster sensitivity using reciprocal-5-nearest-neighbour
  components of the frozen baseline, cosine threshold 0.80.

## Metrics

Effect-level primary metrics:

1. experiment-098 linear NMAE credit on the held reference-DE gate, excluding
   only each row's own target coordinate;
2. weighted sign precision in the top 20% ranked by fitted confidence;
3. complete-gallery PDS, excluding all 300 VCC panel coordinates;
4. predicted-call count proxy and signed precision among predicted calls.

If and only if the effect-level gate passes, a later separately reviewed
stage may realize counts on all three clean harnesses. This experiment does
not authorize count generation.

## Promotion gate

Every condition is required:

1. full-model linear-credit gain over identity is nonnegative on every outer
   surface and at least `+0.003` macro;
2. its macro target-block-bootstrap 95% lower bound is above zero;
3. full beats agreement-only in macro linear credit with lower bound above
   zero;
4. complete-gallery PDS point delta is nonnegative on every surface;
5. macro PDS non-inferiority lower bound is above `-0.001`;
6. top-20% weighted sign precision is nonnegative versus agreement-only on
   every surface and improves by at least `+0.01` macro;
7. predicted-call coverage changes by no more than 5%;
8. the observed full-minus-agreement credit gain exceeds all 20 target-label
   permutation nulls;
9. unsupported targets and coordinates are bit-identical to identity;
10. row norms are restored to `rtol=1e-12, atol=1e-12`, signs never flip, and
    all 300 panel coordinates remain exactly zero;
11. no source-derived target cluster supplies more than half the positive
    macro gain, and every leave-one-cluster-out macro delta is nonnegative.

Any failed condition is a **NO-GO**. Diagnostics cannot rescue a failed
primary gate.

## Scope and execution

- All writes must remain under
  `experiments/104-context-heldout-credit-adapter/`.
- Local execution only, one substantive process at a time.
- Every substantive process runs under a 24 GiB process-tree memory guard.
- No KOLF access, AWS, `cell-eval2`, count generation, challenge package,
  submission, commit, or modification outside this directory.
