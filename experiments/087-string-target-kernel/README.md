# 087 — preregistered STRING target-kernel residual

Status: **completed; non-panel pre-gate failed before held-panel evaluation**.

## Question

Does a compact STRING-network embedding predict a transferable component of
the residual left by frozen experiment-064 arm D?

This is an independent biological-prior ablation parallel to experiment 082's
ESM2 test. It is preregistered before either embedding has produced a held-panel
result. No embedding family may be selected, blended, or tuned after opening
HCT116, HEK293T, or KOLF panel truth.

## Pinned external artifact

- repository: `RasmussenLab/SLIM`
- revision: `5a7e9ade5d0a`
- file: `data/gene_string_embeddings.v0.3.h5`
- bytes: `19,502,608`
- SHA-256:
  `789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d`
- shape: 21,688 gene symbols × 64 float64 features
- VCC-panel coverage: 298/300; missing `EEF1A2`, `KIF21B`
- license: MIT for the repository; retain the upstream STRING-embedding
  citation recorded in `data/external/slim/README.md`

No SLIM expression dataset, fitted model, or prediction is an input. Only the
pinned target embeddings are eligible.

## Frozen protocol

Use experiment 082's corrected infrastructure and exact leakage boundary, but
replace its ESM2 feature block with the 64-dimensional STRING block:

1. Reconstruct frozen experiment-064 arm D independently in each outer donor
   fold, excluding the held donor from prediction, normalization, whitening,
   and every fitted statistic.
2. Use the same deterministic 600 shared non-panel targets and the same
   60/20/20 fit/tune/test hash split as experiment 082. The exact target lists
   and hashes must be serialized before any fit.
3. Remove all 600 development-target coordinates before development
   preprocessing. Use the pinned 17,806-gene final score axis with exactly the
   300 VCC target coordinates removed for held-panel evaluation. Teacher
   fitting and development PDS use only the frozen 5,582-gene `train_mask.npy`
   axis, which is exactly the ordered gene support common to RPE1 and Jurkat;
   unmeasured teacher genes are excluded, never encoded as zero-valued truth.
   Learned residuals are mapped from this teacher axis onto the final axis, and
   every other final-axis coordinate remains the arm-D identity prediction.
4. Standardize STRING features using fit rows only. Missing embeddings are
   zero rows with an explicit missingness flag; they may never be silently
   dropped.
5. Fit only low-rank residual ridge models. Candidate output ranks are
   `{8, 16, 32, 64}`, ridge strengths `{1e-3, 1e-2, 1e-1, 1, 10}`, and residual
   fractions `{0, 0.02, 0.05, 0.10, 0.20}`. Identity is mandatory.
6. Select on non-panel tune PDS, run the untouched non-panel test pre-gate,
   then serialize and hash the exact real and null final weights before
   authenticating panel truth. Every PDS comparison uses raw,
   control-referenced donor deltas; across-perturbation truth centering is
   forbidden.
7. Evaluate the frozen candidate against complete HCT116, HEK293T, and
   273-target KOLF galleries. KOLF remains an independent confirmation and is
   excluded from the primary two-fold macro/bootstrap.
8. Nulls must repeat complete selection using:
   - STRING rows permuted within response-norm and source-coverage strata; and
   - equal-dimensional Gaussian features.

## Promotion gate

Every condition is required:

- nonnegative raw-PDS delta in HCT116 and HEK293T;
- primary HCT116/HEK293T macro delta at least `+0.003`;
- primary source-cluster bootstrap 95% lower bound above zero;
- real gain above both full-selection null p95 values;
- nonnegative KOLF delta;
- median angular displacement below 10 degrees in every held fold;
- nonnegative complete-source and partial-source stratum deltas;
- no single source-derived target cluster supplies more than half of the
  positive gain, and leave-one-cluster-out delta remains nonnegative;
- unchanged target/gallery coverage.

A pass authorizes only one separately reviewed HCT116 count-level
confirmation. It does not authorize a challenge build or submission.

## Execution boundary

Local work is limited to synthetic tests under a 1 GiB process-tree guard.
Real staging and evaluation belong on one isolated AWS worker under the global
VCC-heavy lock. No implementation may launch until experiment 082's shared
staging/axis/leakage defects are repaired and independently reviewed.
