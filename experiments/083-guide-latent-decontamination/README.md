# 083 — Guide-latent decontamination

## Status

Implemented and locally testable, but **hard-stopped before fitting or held
evaluation for the current Replogle K562 GWPS source**.

On 2026-08-30, experiment 109 invoked this extractor in a guarded source-only
cache-generation attempt. The source matrix was streamed once, but the compact
reconstruction gate failed before cache publication or fitting because target
reconstruction incorrectly inherited the per-guide cell threshold. The
pre-fit repair keeps the >=20-cell rule for guide effects while reconstructing
each compact target from all cells assigned to every target having >=20 total
cells. The `2e-5` integrity tolerance is unchanged.

The read-only source audit found:

- 1,989,578 cells, 8,248 measured genes;
- 11,187 guide-pair categories assigned unambiguously to 9,867 target labels;
- 775 perturbation targets with at least two retained (>=20-cell) guide-pair
  observations;
- 272/300 VCC-panel targets in the existing pooled compact bank, of which 269
  retain at least one >=20-cell guide-pair row;
- only **8/300** panel targets with at least two guide-pair observations;
- the same 8 have at least two guide-pair rows identifiable in the pinned
  sequence/seed table.
- only **6/300** of those enter the exact arm-D K562/PerturbMatch source slot;
  the other two retain the unchanged fallback.

Hypothesis 3 in experiment 078 predeclared a minimum of 100 evaluable
multi-guide panel targets. The current source therefore cannot pass its own
promotion gate. `aws/launch.sh`, `aws/run.sh`, and `fit.py` all fail closed
before expensive extraction/modeling when this condition is false.

## Question

Can a source-wide learned decomposition preserve target-shared K562 biology
while removing the feature-predictable guide/seed component that hard
whole-row deletion discarded in experiment 068?

```text
guide_effect(target, guide_pair)
    = target_shared_on_target
    + low_rank(sequence_and_seed_features)
    + residual
```

This experiment never deletes a source row or manually zeros a named
off-target coordinate. Unsupported, single-guide, or otherwise
unidentifiable targets retain the original K562 compact row byte-for-byte.

## Data interpretation

The scPerturb K562 GWPS H5AD is guide-level, but each `guide_id` is usually a
*paired-guide construct*, for example:

```text
A1BG_+_58858964.23-P1|A1BG_-_58858788.23-P1
```

Most targets have one such construct. A small minority have independent P1
and P2 constructs. Those independent construct rows are the replication unit
needed to separate a target-shared latent effect from guide-specific signal.
The two protospacers inside one construct are not treated as independent
expression observations.

The pinned Hartman et al. source table from experiment 068 supplies
protospacer and seed strings for 7,738 local guide IDs. Features use sequence,
seed-match geometry, direction, and candidate burden only. The expression
columns (`target_log2fc`, `offtarget_log2fc`) and named off-target genes are
explicitly forbidden as predictors.

## Leakage-safe protocol

1. Inspect only H5AD metadata and the pinned seed table. Abort before matrix
   extraction unless at least 100 panel targets have two identifiable
   guide-pair constructs.
2. Stream the dense HDF5 matrix in row blocks. Persist >=20-cell guide-pair
   pseudobulks, but independently pool every cell assigned to targets having
   >=20 total cells for compact reconstruction; never materialize the full
   cell-by-gene matrix.
3. Reproduce the existing compact K562 pooled effect for every compact target
   within `2e-5` maximum absolute error. Per-guide filtering must not alter
   this target-level integrity estimand.
4. Remove all 300 VCC target-gene coordinates from every fit, transform,
   held-guide score, and retrieval axis. Also remove the non-panel development
   target coordinates from the learned correction so they cannot dominate
   held-guide selection or mix through the dense whitener. The final held-panel
   PDS axis still excludes **exactly** the 300 VCC target coordinates, not the
   development targets.
5. Hash-split identifiable **non-panel targets** 60/20/20 into fit, tune, and
   test. No target crosses partitions.
6. Fit output ranks `{8,16,32}`, ridge strengths `{0.1,1,10,100}`, blend
   fractions `{0.25,0.5,0.75,1}`, and guide-cell weight powers `{0,0.5}` on
   fit targets only. The output basis is learned from within-target guide
   residuals. The nuisance map is ridge-fitted from within-target-centered
   sequence/seed features to those residual scores.
7. Select on deterministic leave-one-guide-out prediction in tune targets.
   A nonidentity candidate must improve mean cosine, not worsen normalized
   MSE, and be nonnegative in every observed guide-count stratum; otherwise
   identity wins.
8. Refit the frozen configuration on fit+tune, test once on untouched
   non-panel targets, and require positive cosine delta, nonworse normalized
   MSE, and a source-cluster bootstrap lower bound above zero.
9. Only then refit once on all non-panel targets and freeze:
   model parameters, panel K562 adjustments, shrink controls, two null
   families, and hashes. This all occurs before held panel truth is opened.
10. A separate evaluator authenticates the freeze and performs no fitting.
    It reconstructs experiment-064 arm D independently through the
    experiment-068 baseline path, requiring exact per-target PDS fidelity.
    It uses the full available VCC-panel truth gallery for held HCT116,
    HEK293T, and independent KOLF. KOLF is never a source.

Experiment 080 caches, coordinate masks, galleries, null-selection results,
and adapters are prohibited dependencies because independent review found
leakage defects in that path.

## Learned decomposition

For fit-target guide effects `Y` and guide features `X`, both are centered
within target. A randomized SVD of `Y` supplies a rank-constrained response
basis `V`; ridge regression learns `W`:

```text
Y_centered ≈ X_centered W V
```

For a new multi-guide target:

```text
adjusted_guide = guide_effect - X W V
theta_target = robust_weighted_center(adjusted_guides)
cleaned = pooled_target + eta * (theta_target - pooled_target)
```

The correction is continuous. It is inserted only into K562's portion of the
frozen 50:50 K562/PerturbMatch source slot, then passed through the frozen
source-only gamma-0.65 operator and unchanged post-whitening CD4 residual.

## Predeclared nulls and controls

- **Ordinary averaging / identity:** exact arm-D baseline.
- **Guide-to-target assignment null:** complete feature bundles are permuted
  among targets matched on guide count and source-response norm decile.
- **Seed-feature null:** feature rows are permuted within the same strata.
- **Shrink control:** each original K562 row is direction-preservingly scaled
  to the cleaned row norm. The learned candidate must beat this control.

Each permutation replicate repeats the complete candidate-grid selection and
untouched non-panel test gate. For computational fairness, both null families
receive the same response-derived low-rank bases as the real search; only the
feature/assignment link is broken. One hundred frozen replicates per family
are produced before panel truth is opened.

## Source-derived bootstrap

Panel targets are clustered separately in each outer fold, without held
truth, by reciprocal 8-nearest-neighbor connected components of that fold's
frozen arm-D source query. Each of 20,000 bootstrap replicates independently
resamples the HCT-fold and HEK-fold clusters, then averages the two fold
means. This prevents one fold's held donor from influencing the other fold's
cluster definition.

## Promotion gates

All are mandatory:

- at least 100 identifiable multi-guide panel targets;
- nonidentity candidate;
- nonnegative PDS delta for held HCT116, HEK293T, and KOLF;
- equal-donor macro raw-PDS delta at least `+0.003`;
- source-cluster bootstrap 95% lower bound above zero;
- observed gain above both permutation-null p95 values;
- exact experiment-064 arm-D baseline PDS reconstruction;
- unchanged target coverage and exact zero adjustment for every unsupported
  target;
- nonnegative mean delta in both the seed-evidence and unflagged
  multi-guide strata;
- gain greater than the direction-preserving shrink control;
- every leave-one-source-cluster-out macro delta nonnegative.

Only a passing effect-level result may authorize a separate HCT116
count-level confirmation. This package contains no count generation,
`cell-eval2`, VCC submission, or model-name path.

## Local checks

```bash
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report experiments/083-guide-latent-decontamination/local-test-memory.json \
  -- \
  uv run --no-project --python 3.12 --with numpy,h5py python \
    experiments/083-guide-latent-decontamination/test_model.py
```

The current structural audit is deliberately cheap and read-only:

```bash
uv run --no-project --python 3.12 --with numpy,h5py python \
  experiments/083-guide-latent-decontamination/inspect_source.py \
  --allow-network \
  --hash-h5ad \
  --output experiments/083-guide-latent-decontamination/source-audit.json
```

Adding `--require-promotion-coverage` exits 78 for the present source.

## Cloud launcher

The launcher can use only an already-running, SSM-online worker with at least
64 GiB RAM. It requires explicit clearance, reruns the local structural audit,
and refuses to upload or dispatch anything when the 100-target coverage gate
fails:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES INSTANCE_ID=i-... \
  bash experiments/083-guide-latent-decontamination/aws/launch.sh
```

No cloud run, fitting, held evaluation, or submission has been launched. The
single local source-streaming attempt described above produced no cache.
