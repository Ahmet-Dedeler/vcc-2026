# 080 — Leakage-safe listwise metric adapter

## Question

Can same-target cross-donor supervision learn a small residual correction to
the frozen post-whitening geometry that improves PDS retrieval beyond arm D?

The deployable map is:

```text
z' = normalize(z + beta * B U V^T B^T z)
```

`B` is a source-only randomized right-singular basis and `UV^T` has rank
4, 8, or 16. The learned residual is embedded only on genes observed by both
public teacher donors. Its listwise loss assigns the same perturbation in a
different donor as the positive and every other perturbation in that donor as
a negative.

## Frozen upstream geometry

The cache stage keeps three deliberately separate representations.

The saved final-panel baseline exactly reconstructs experiment-064 arm D:

- K562 source slot: 50:50 incumbent/PerturbMatch where both target rows exist;
- source-only union whitening at gamma `0.65`;
- post-whitening mean3 CD4 `gene_norm` residual at fraction `0.05`;
- exact experiment-062 held-HCT116 and held-HEK293T operators.

For that fixed baseline only, every main donor is centered and norm-calibrated
over all mapped source coordinates, exactly as in experiment 064, and the 300
VCC coordinates are blocked afterward. The stage then performs a second,
independent reconstruction directly from the authenticated raw compact banks
through experiment 064's original functions. The saved 221-row panel query
must match this independently reconstructed arm D over the complete 17,806-gene
retrieval axis with maximum absolute error at most `1e-6`.

This baseline reproduction is not a learned feature path. Training, source
clustering, and final-panel application use separate contracts:

- The fit/tune/test branch removes the 300 VCC coordinates and all 600
  historical target coordinates before estimating source means, source norms,
  whitening restrictions, adapter bases, weights, hyperparameters, nulls, or
  source-query clusters. It applies the principal restriction
`W_RR = E_R diag(lambda^-gamma/2) E_R.T` directly. For forbidden coordinates
`x_F = 0`, `(xW)_R = x_R W_RR`; no historical target coordinate can therefore
mix back into a retained coordinate through the dense whitener.
- A separate `cluster_query.npy` uses the final panel rows but the same
  retained-coordinate-only centering, norm calibration, and principal
  whitener restriction as the training branch. It excludes all 300 VCC and
  all 600 historical target coordinates and is the only input to source-query
  clustering, bootstrap grouping, and leave-one-cluster-out dependence checks.
- The final-panel branch starts from the fixed, independently verified exp064
  baseline, removes exactly the 300 VCC coordinates from the saved score
  vector, and applies only the already frozen adapter. Historical non-panel
  target coordinates are restored as ordinary score coordinates.

Thus excluded target coordinates can affect exp064's already-frozen global
source scale, but they cannot affect any learned statistic or model-selection
decision in experiment 080. The final HCT116/HEK293T/KOLF score axis is pinned
to exactly 17,806 genes with order and position hashes. The 49 historical
targets absent from the official gene axis are valid labels and simply have no
coordinate to remove. Extensions outside exp064's 221 targets are audited
explicitly: incumbent K562 is retained when PerturbMatch is absent, only
available main donors are averaged, and CD4 is added only where measured.

K562 is not a valid outer fold because PerturbMatch reanalyzes the same K562
assay. The strict outer folds are Orion HCT116 and Orion HEK293T.

## Leakage-safe nested protocol

1. Remove all 300 VCC target coordinates throughout. Remove the 600 non-panel
   target coordinates only from non-panel fitting and validation; retain them
   for final VCC-panel scoring.
2. Hash-split the 600 non-panel perturbations into 360 fit, 120 tune, and 120
   test targets.
3. For each outer donor, its prediction and frozen operator use only K562 plus
   the other Orion donor.
4. Inner teacher holdout alternates RPE1 and Nadig-Jurkat. The held teacher is
   absent from the basis and adapter fit.
5. Search identity plus the exact grid:

```text
rank:            4, 8, 16
beta:            0.02, 0.05, 0.10, 0.20
temperature:     0.03, 0.07, 0.15
regularization:  1e-4, 1e-3, 1e-2
```

The sole regularizer is `0.5 * lambda * ||UV^T||_F^2`; there is no separate
penalty on `U` or `V`. Every tune subset is ranked against the complete
600-target truth gallery. A candidate must improve or tie every outer-role ×
held-teacher cell and remain below 10 degrees median angular displacement.
Among eligible candidates, selection exactly maximizes macro tune PDS. Only an
exact floating-point tie uses the deterministic ordering: lower rank, lower
beta, stronger regularization, temperature nearest `0.07`, then candidate key.
There is no near-tie window. Otherwise identity is selected.

6. Refit the frozen configuration on fit+tune and evaluate the untouched
   non-panel test queries against the complete 600-target gallery.
7. Freeze one coherent per-replicate permutation plan and reuse it across
   inner selection, pretest, and final refit while repeating the complete
   candidate grid.
8. If the pre-gate passes, fit and serialize the exact real and null-selected
   final weights on all 600 non-panel targets. Also build and serialize
   fold-specific clusters only from each fold's staged leakage-safe
   `cluster_query.npy`.
9. Authenticate panel truth only after those model and cluster artifacts are
   hashed. Final evaluation reloads the frozen weights; no post-panel refit is
   permitted.

The primary macro, bootstrap, and permutation threshold use only HCT116 and
HEK293T. KOLF is a separate required nonnegative confirmation and cannot rescue
an otherwise failing primary result. The bootstrap independently resamples
each primary fold's source-query clusters.

Nine candidate-invariant fit contexts cache normalization, bases, projections,
and base scores. The 36 configurations sharing a rank use batched linear
algebra with separate parameters and Adam state.

## Independent KOLF confirmation

The current artifacts support a required untouched KOLF2.1J confirmation:

- all 273 KOLF/VCC shared targets form the complete KOLF gallery;
- its arm-D query uses K562/PerturbMatch plus both Orion donors and CD4;
- its operator is the pinned three-donor production gamma-0.65 whitener;
- KOLF contributes zero rows to prediction, fitting, selection, whitening, or
  clustering.

The chosen configuration and KOLF-specific weights are fitted only from the
600 non-panel rows and frozen before KOLF files are authenticated. Promotion
fails closed if KOLF support is absent or the KOLF delta is negative.

## Promotion gate

Every item is required:

- selected adapter is non-identity and rank at most 16;
- nonnegative PDS delta in both valid outer folds;
- nonnegative PDS delta in the independent KOLF fold;
- macro raw-PDS delta at least `+0.003`;
- target-cluster bootstrap 95% lower bound above zero;
- observed gain above the 95th percentile of the response-norm-stratified
  label-permutation null;
- median angular displacement below 10 degrees separately in every fold;
- unchanged target coverage;
- no single source-derived target cluster carries most of the positive gain,
  and every leave-one-cluster-out macro delta remains nonnegative.

No count generation, `cell-eval2`, challenge prediction, or submission is
performed. A passing result would authorize at most one separate guarded
HCT116 count-level confirmation.

## Local synthetic validation

The tests do not open held truth or run heavy donor processing. Permanent
regressions include exact cached/uncached semantics and batched-versus-serial
optimizer equivalence for ranks 4, 8, and 16 with two teachers, frozen
permutations, and forced gradient clipping:

```bash
uv run --no-project --python 3.12.12 --with 'numpy==2.5.2' python \
  experiments/080-listwise-metric-adapter/memory_guard.py \
  --max-rss-gib 1 \
  --poll-seconds 0.1 \
  --report experiments/080-listwise-metric-adapter/local-test-memory.json \
  -- \
  uv run --no-project --python 3.12.12 --with 'numpy==2.5.2' python \
    experiments/080-listwise-metric-adapter/test_run.py
```

## Cloud execution

The launcher accepts only an already-running, SSM-online instance with at
least 64 GiB nominal RAM and 48 GiB currently available. It performs a
read-only lock/process/service preflight and requires an explicit clearance
environment variable:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES INSTANCE_ID=i-... \
  bash experiments/080-listwise-metric-adapter/aws/launch.sh
```

Before dispatch, the launcher rejects reused S3 RUN_ID prefixes and hashes the
exact source archive and runner. The worker verifies both, refuses reused local
work directories, pins Python 3.12.12 and NumPy 2.5.2, checks collisions and
available memory, and enforces one global VCC-heavy lock. It contains no
submission command.
