# 082 — Leakage-safe biological target-kernel residual

Status: **NO-GO at the frozen non-panel pre-gate; held-panel truth was never
opened**.

## Question

Can ESM2 protein embeddings predict a small, transferable part of the
cross-donor residual that frozen experiment-064 arm D systematically misses?

The candidate is deliberately residual and near-identity:

```text
r(d,p)       = truth(d,p) - arm_D(excluding d,p)
r_hat(p)     = ridge(standardized ESM2(p)) in a source-only residual PCA basis
prediction   = arm_D + beta * r_hat(p)
```

Arm D remains the comparator and supplies every prediction. A target with no
usable ESM2 row receives exactly the arm-D prediction.

## Frozen upstream geometry

`stage_cache.py` reconstructs experiment-064 arm D directly from the pinned
experiment-062/064 implementations:

- K562 source slot: 50:50 incumbent/PerturbMatch where both rows exist;
- source-only union whitening at gamma `0.65`;
- post-whitening mean3 CD4 `gene_norm` residual at fraction `0.05`;
- exact experiment-062 held-HCT116 and held-HEK293T operators;
- the pinned three-donor production operator for independent KOLF2.1J
  confirmation.

All available historical development-target coordinates are excluded from
source means, scale statistics, effects, CD4 statistics, and the principal
block of the dense operator for development fit/tune/test, so they cannot mix
into safe output genes. The 109/600 targets absent from the compact master
axis remain valid perturbation rows; 49 of those are also absent from the
official axis and therefore have no coordinate to remove. Both sets are
counted and hashed. The final panel path removes exactly the 300 VCC target
coordinates and uses the pinned complete 17,806-gene held-donor score axis,
retaining every available historical target coordinate. On the pinned
221-target intersection, the saved query must match the direct arm-D algebra
within `1e-6` on that complete score axis.

Residual fitting and all development PDS calculations use the frozen
5,582-gene `train_mask`: the genes measured by both RPE1 and Nadig-Jurkat
after excluding all development and VCC target coordinates. A gene absent
from either teacher is excluded, never represented as zero-valued truth.

The Orion `genes.txt` files have 38,606 ordered columns, 38,584 unique
symbols, 21 duplicated symbols, and 22 duplicate rows. Their full ordered axes
are retained for matrix-shape and centering semantics. Unique official-gene
projection deterministically uses the last occurrence, matching the frozen
upstream arm-D mapping, and the duplicate-axis contract is recorded in cache
metadata and held-truth audits.

## Leakage-safe protocol

1. Remove all 300 VCC target coordinates and all 600 non-panel target
   coordinates before arm-D transformation or retrieval.
2. Hash-split the 600 shared non-panel perturbations into 60% fit, 20% tune,
   and 20% test partitions.
3. Standardize the 5,120-dimensional ESM2 features using only the current fit
   partition. Missing embeddings are retained as identity predictions.
4. In each outer role, arm D and its operator exclude the held outer donor.
5. Inner donor holdout alternates RPE1 and Nadig-Jurkat. Residual PCA and
   ridge fitting use the other teacher only.
6. Search identity plus exactly:

```text
residual PCA rank: 16, 32, 64
ridge strength:    1e-3, 1e-2, 1e-1, 1, 10
beta:              0.02, 0.05, 0.10, 0.20
```

   Ridge uses the standardized linear target kernel `X X^T / 5120`.
7. Every fit/tune/test query is ranked against the complete 600-target truth
   gallery. Partition membership controls query rows only. Truth is the raw,
   control-referenced donor delta; subtracting an across-perturbation mean is
   forbidden because official PDS uses uncentered effects.
8. Select one configuration by macro tune PDS, requiring a nonnegative delta
   in every outer-role × inner-held-teacher cell; otherwise select identity.
9. Freeze the selection, refit on fit+tune, and evaluate the untouched
   non-panel test split.
10. Repeat the complete selection for two null families:
   - ESM2 rows permuted within split, arm-D response-norm decile, and
     source-coverage strata;
   - an equal-dimensional deterministic Gaussian embedding with the same
     missing-row mask.
11. Only a real candidate that beats both non-panel null p95 values may
    continue. Exact fold-specific final real and null fitted weights are then
    serialized and hashed on all 600 non-panel targets, still before held
    panel rows are selected and loaded for evaluation. HCT116 and HEK293T
    files may already have been authenticated and used as source donors in
    the opposite outer fold; they are never used in their own held fold.
12. Source-query clusters are also frozen before held-panel evaluation. The exact
    saved weights are reloaded and evaluated without refitting against the
    complete 290-target HCT116, 294-target HEK293T, and 273-target KOLF
    galleries, always against uncentered control-referenced truth. HCT116 and
    HEK293T alone define the primary macro, bootstrap, and null thresholds.
    KOLF is reported and gated separately.

No panel target selects an embedding feature, rank, ridge, beta, stratum,
cluster, or model family.

## Decision gate

Promotion requires every common experiment-078 condition:

- nonnegative raw-PDS delta in both valid outer folds;
- macro raw-PDS delta at least `+0.003`;
- source-derived target-cluster bootstrap 95% lower bound above zero;
- observed panel gain above both embedding-permutation and Gaussian null p95;
- unchanged target coverage;
- no obvious dependence on one source-derived target cluster.

It additionally requires, separately in every outer fold, median angular
displacement below 10 degrees and nonnegative complete-source and
partial-source stratum deltas. Per-target PDS is always computed against the
complete gallery before any stratum mean is taken. Bootstrap clusters are
derived independently from each outer fold's source-only arm-D queries.

KOLF is required as an independent confirmation fold. It is never averaged
into the primary HCT116/HEK293T macro, cluster bootstrap, aggregate source
strata, or null p95. Its query uses all
three main source donors and the pinned production arm-D operator, while its
ESM2 residual weights are fitted only from the same 600 non-panel rows.
KOLF truth is neither authenticated nor opened until every real/null model
and source-query cluster file has been frozen. Promotion requires nonnegative
KOLF delta plus the same per-fold angular and complete/partial-source gates.

## ESM2 artifact

Only this file is allowed:

```text
dataset:  altoslabs/primeflow-vcc-datasets
revision: d773a3476d15e5d32d385b13cebf87e62580097e
file:     ESM2_pert_features.parquet
bytes:    595319311
sha256:   f05fe94e3a592bfa0973b3816b956af8c0a99387484d4c68ee89e0eac0b21db3
```

The associated expression H5AD bundle is forbidden by both the launcher and
worker. The Parquet has 19,790 rows, 5,120 float features named `0` through
`5119`, and one string index column.

## Local synthetic validation

The tests do not open real donor matrices or download the ESM2 artifact:

```bash
uv run --no-project --python 3.12 --with numpy --with pyarrow python \
  experiments/082-biological-target-kernel/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report experiments/082-biological-target-kernel/local-test-memory.json \
  -- \
  uv run --no-project --python 3.12 --with numpy --with pyarrow python \
    experiments/082-biological-target-kernel/test_run.py
```

## Cloud execution

The launcher accepts only an already-running, SSM-online instance with at
least 64 GiB RAM. It cannot create or start an instance and requires explicit
clearance:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES INSTANCE_ID=i-... \
  bash experiments/082-biological-target-kernel/aws/launch.sh
```

The launcher rejects a run if either its source or output S3 prefix already
contains a key, and the worker independently rejects a reused output prefix.
Worker run IDs must be exact UTC timestamp strings (`YYYYMMDDTHHMMSSZ`);
canonical work-root containment is verified and any existing local run path,
including a symlink, is rejected without deletion.
The worker enforces a global VCC-heavy lock, 48 GiB for arm-D staging, 24 GiB
for evaluation, and 1 GiB for synthetic tests. It downloads exactly one
Hugging Face file and rejects any H5AD artifact. The launcher computes the
source-archive and runner SHA-256 values at dispatch; the worker verifies both,
checks available memory, and rejects broad process/service collisions. There
is no count generation, `cell-eval2`, `vcc prep`, or submission command.

## Result

The first valid cloud run, `20260830T042050Z`, selected exact identity before
held-panel access. All 60 nonidentity ESM2 residual candidates were negative
in all four outer-donor × inner-teacher cells. The least harmful candidate,
rank 32 / ridge 1 / beta 0.02, had macro raw-PDS delta `-0.005120` and minimum
cell delta `-0.005718`.

All 22 tests and every leakage audit passed. All 27 recovered artifact hashes
verify. Peak process-tree RSS was 4.44 GiB for staging and 0.52 GiB for
evaluation. The AWS worker was stopped after recovery. No held-panel truth,
count generation, `cell-eval2`, or challenge submission was used.
