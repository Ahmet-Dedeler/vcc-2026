# 103 — residual-only a15-direction confirmation on KOLF2.1J

Status: **completed; preregistered primary gate is a NO-GO**.

Initially preregistered 2026-08-29; adversarially audited and amended
pre-truth on 2026-08-30.

## Question

Does removing each real source donor's projection onto its own donor-mean
response axis improve the staged a15 direction after the otherwise unchanged
PerturbMatch, whitening, and CD4 path?

This is a production-shaped follow-up to experiment 102's post-hoc
residual-only diagnostic. It intentionally changes the axis-fit population and
the downstream transform to match the staged challenge path, so it is not a
bit-for-bit replication of experiment 102. Positive or negative evidence
applies only to this **donor-axis surrogate**. It does not validate or refute
official COMPASS, whose target axis requires target-context perturbation
responses that VCC does not provide.

## Held surface

The only held surface is KOLF2.1J.

- No KOLF file is accepted by or visible to the freeze command.
- All 300 production target rows, masks, source axes, projection loadings,
  source clusters, baseline, candidate, null variants, hashes, and config are
  sealed before KOLF authentication.
- KOLF contributes no feature, fit, threshold, target selection, axis
  selection, cluster, whitening row, or fallback.
- The score command creates an exclusive local access marker before hashing or
  opening the KOLF directory.
- After truth access, only fixed axis validation, row restriction, PDS,
  bootstrap, null comparison, and reporting are permitted.

The numerical KOLF effect values remain untouched during implementation and
freeze.

## Exact staged a15 direction baseline

The baseline is the direction matrix recorded in
`submissions/a15-cd4f010.vcc.meta.json`:

1. K562-GWPS, Orion-HCT116, and Orion-HEK293T with equal weights among donors
   carrying each target;
2. centered donor effects with `gene_norm=False`, `eb=False`;
3. `alpha=0.982`, `beta=0`;
4. current 50:50 PerturbMatch replacement of the K562 slot where both source
   rows exist;
5. the pinned common-essential-pruned three-donor union whitener,
   `gamma=0.65`, `lambda=0`, with candidate-specific on-axis norm
   preservation;
6. the frozen complete-three-condition CD4 residual in `gene_norm` mode at
   experiment-099 target ratio `f=0.10` (solved fraction `0.094727...`), with
   its existing post-add global preservation.

The runner constructs this baseline twice:

- once through `build_production_effects`;
- once by explicitly invoking the production pre-whitening, PerturbMatch,
  whitening, and CD4 stages.

The exact a15 effect matrix was re-exported from its pinned Linux build tree
with the original eight-thread BLAS settings and preserved under
`work/staged-a15-linux/`. Its array hash reproduces a15 exactly
(`2980d706...`). The local Apple-silicon reconstruction must match a15's
pre-whitening direction bit-for-bit (`0601336f...`) and the preserved final
matrix elementwise at `rtol=1e-12, atol=1e-12`; the observed maximum absolute
difference before freezing was `1.34e-14`. The exact preserved Linux matrix,
not the local reconstruction, is the frozen baseline scored against KOLF.
Current implementation, protocol, experiment-102 evidence, exporter, a15
metadata, and all model inputs are pinned. Failure is a hard blocker; no nearby
baseline is allowed.

## Exact score axis without KOLF access

The final eligible score axis is reconstructed from non-KOLF sources:

```text
official genes
intersect (Orion-HCT116 genes union Orion-HEK293T genes)
minus all 300 VCC panel genes
```

It must contain exactly 17,806 genes and match the preregistered ordered-axis
SHA-256:

```text
280aab34a0e86f27db8ade213e000f9066c8d5f9b89edd3ea966aba061ba9ae9
```

This is the previously established KOLF/VCC non-panel score-axis hash, but the
axis used during freeze is derived only from the two source-donor gene lists.
The 17,942-column production whitening mask remains unchanged.

## One-change candidate

For each real donor `c`, fit on raw donor responses:

```text
S_c = donor perturbation rows whose target is not in the 300-gene VCC panel
m_c = mean(q in S_c) delta_cq
u_c = normalize(m_c)
```

All vectors are first mapped to the official axis and restricted to the exact
17,806 score coordinates. Forbidden coordinates therefore cannot affect
`m_c`, `u_c`, a target projection, a residual, a norm, or a null.

For a carried VCC target `p`, the candidate replaces only the real donor's
centered row:

```text
s'_cp = s_cp - <s_cp, u_c> u_c
```

The production donor scaling, source weights, global effect scale, target
coverage, and fallback are unchanged. On targets receiving PerturbMatch, the
same 50:50 replacement is retained: the surviving real K562 half is
residualized while the external PerturbMatch half is untouched.

On a target receiving PerturbMatch, only genes represented in PerturbMatch
retain the 50% real-K562 share; PM-missing genes retain the full real K562
slot. The projection subtraction follows that exact per-gene retention mask.

The candidate then runs through the exact same pruned production whitener and
CD4 functions. Their masks, support, parameters, and renormalization rules are
not refitted or modified.

## Null and clusters

Twenty deterministic null variants independently permute each donor's scalar
projection loadings across the VCC target labels carried by that donor.
Axes, loading marginals, source coverage, PM support, residual-free baseline
rows, whitening, and CD4 remain fixed.

Target clusters are frozen from the complete 300-row production baseline on
the exact 17,806-gene score axis. They are reciprocal-5-nearest-neighbor
connected components requiring cosine at least `0.80`. KOLF rows and values
never enter the graph.

## Evaluation

After the sealed freeze:

1. authenticate the pinned KOLF artifact hashes;
2. verify the complete ordered 273-target VCC-panel subset;
3. verify the exact ordered 17,806-gene score axis;
4. read only KOLF's stored effect matrix on that axis;
5. compute raw, uncentered complete-gallery PDS.

No per-gene or per-target centering is applied to KOLF truth or predictions.
KOLF's placeholder `b_ctrl.npy` is authenticated but never loaded as biology.

The target-cluster bootstrap uses 20,000 deterministic draws of whole frozen
source-derived clusters.

## Promotion gate

The candidate is a **GO** only if all are true:

1. mean raw PDS candidate minus baseline is at least `0.003`;
2. the target-cluster bootstrap 95% lower bound is above zero;
3. the observed mean delta exceeds the 95th percentile of the 20 source-wise
   target-label permutation nulls, using the conservative empirical upper
   order statistic, and its add-one one-sided empirical p-value is at most
   `0.05`;
4. source coverage is exactly unchanged;
5. every implementation hash, input hash, score-axis hash, target-axis hash,
   OFF reconstruction, PM/CD4 support, whitening mask, freeze-order, and
   no-post-truth-fit check passes.

No diagnostic can rescue a failed primary gate.

A GO confirms this residual-only mechanism on the frozen staged-a15 direction.
It remains effect-level evidence because KOLF has no biological control
profile, so it does not by itself authorize count generation or challenge
submission.

## Commands for later review

Freeze without any KOLF path:

```bash
uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/103-residual-only-production/run.py freeze \
  --max-rss-gib 24
```

Only after reviewing the sealed manifest, score exactly once:

```bash
uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/103-residual-only-production/run.py score \
  --kolf data/external/kolf2 --max-rss-gib 24
```

Synthetic tests:

```bash
uv run python experiments/103-residual-only-production/test_run.py
```

## Hard scope

- Local only; no AWS.
- No KOLF access during implementation or freeze.
- No `cell-eval2`.
- No count build, challenge package, submission, or commit.
- Every generated artifact stays under
  `experiments/103-residual-only-production/`.

KOLF is hypothesis-specific held truth here, not globally unseen project data:
earlier disjoint experiments used it for other questions. Experiment 103's
candidate, score axis, null, bootstrap, and promotion gate were fixed without
reading KOLF outcomes for this residual-only question.
