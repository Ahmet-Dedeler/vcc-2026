# 078 — Remaining-hypothesis audit after experiments 001–077

Date: **2026-08-23**

## Conclusion

Only three hypotheses survive this audit:

1. **Cross-donor listwise metric adapter** — highest priority and executable from
   existing compact effect tables.
2. **Biological-prior target-kernel residual** — tests whether target identity
   contains transferable information absent from expression-only donor pooling.
3. **Guide-level latent on-target decontamination** — replaces the failed
   all-or-nothing off-target deletion with a learned source-wide decomposition.

They act on different objects: the first learns the **output metric**, the
second learns across **perturbation identities**, and the third improves the
**source measurement** before donor transfer. None is another alpha, emission,
whitening-strength, donor-weighting, target-routing, or donor-addition sweep.

Separately, experiment 075 is now a **near-threshold existing candidate for
independent confirmation**, not a rejected family and not one of the three new
hypotheses below.

The default comparison for all three is the frozen experiment-064 arm D
(PerturbMatch replacement plus the post-whitening CD4 residual), using its
existing source-only union whitener and preprocessing unless a protocol below
explicitly says otherwise.

## What is already occupied or closed

| Family | Evidence through 077 | Audit decision |
| --- | --- | --- |
| Global/per-target magnitude, confidence, per-gene calibration | 001, 028–037, 063, 066, 077 protocol | Closed as a new search family. Experiment 077 is an unfinished frozen confirmation, not a new hypothesis. |
| SVD truncation, EB, MODZ, inverse-variance pooling, ordinary per-gene weighting | 004–016, 022, 039 | Closed. |
| Basal similarity and direct context adaptation | 010, 041, 046, 067 | Closed. Full covariance/TTA is not selected without a new pretrained model. |
| Whitening variants and source-subspace composition | 017–025, 038, 040, 042, 048, 052, 057, 060, 062, 069 | The production geometry is strong; auxiliary, noise, per-donor, aligned and composed variants have been tested or rejected. |
| Gallery-density/CSLS corrections | 044, 073 | Pre- and post-whitening CSLS distillation are rejected. |
| Reciprocal local query centering | 075, 079 | The original result selected `k=64, beta=0.5` and improved all three outer folds by `+0.002943142`, but independent KOLF confirmation was only `+0.0002962724` with CI `[-0.0022893773,+0.0028549881]`. Closed; do not promote. |
| Fixed donor additions or target routing | 050–060, 064, 066, 071–074, 081 | KOLF, Jiang, Feng, SIGrecon, and sparse IGVF cardiomyocyte transfer are not production donors; fixed arm D remains selected. Experiment 071 is a packaging path, not a scientific hypothesis. |
| Count-template and emission variants | 043, 049, 053, 061, 065 | Template emission is rejected. Experiment 065 was staged/scoring but has no recorded aggregate, so it remains pending rather than reusable evidence. |
| Hard seed/off-target deletion | 068 | Rejected. The result specifically says whole-row removal is too blunt and leaves a learned guide-level decomposition open. |
| VIPerturb | 070 | Corrected evaluation completed and was negative in both held folds: macro PDS delta `-0.0069624664`, bootstrap 95% CI `[-0.0108373282,-0.0034542003]`. Closed. |

Pending 065 and 077 must be resolved under their already-declared protocols.
Experiments 070, 075/079, and 081 are now closed. These occupied paths do not
justify a fourth new experiment here.

## Live-resource addendum: useful evidence, not a new ready donor

Metadata was checked on 2026-08-23 without downloading payloads:

- [`cfy2yue/perturbseq_normalized`](https://huggingface.co/datasets/cfy2yue/perturbseq_normalized),
  revision `3f96df93440d35d6e2511c1406550e6863b85039`, contains 254
  normalized H5AD files. Its gene-perturbation additions beyond the existing
  scPerturb catalogue are the two already-used X-Atlas/Orion contexts,
  HCT116 and HEK293T. It contains no Pisces release.
- [`altoslabs/primeflow-vcc-datasets`](https://huggingface.co/datasets/altoslabs/primeflow-vcc-datasets),
  revision `d773a3476d15e5d32d385b13cebf87e62580097e`, repackages the
  already-audited CD4, Feng, Jiang, Nadig, Replogle and VCC-family data. Its
  genuinely convenient new compact asset is
  `ESM2_pert_features.parquet` (595,319,311 bytes), which enables experiment 2
  without acquiring another perturbation-expression atlas.
- [`Xaira-Therapeutics/X-Atlas-Pisces`](https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces),
  revision `8ace335dfe4d23d8fb36cf0df55f7a643991475e`, currently has only
  README, LICENSE and `.gitattributes` (24,308 bytes total). Its card says 16
  biological contexts and that data are “coming soon.”
- The [Pisces preprint](https://doi.org/10.64898/2026.03.18.712807) and
  [`X-Cell`](https://huggingface.co/Xaira-Therapeutics/X-Cell) card describe
  seven CRISPRi screens; the larger context count includes differentiated
  populations within the iPSC Multi-Diff screen rather than 16 independent
  broad donor screens. X-Cell revision
  `07737b494a3ea06c6a579d389fda0216ba61bf02` still has no weights or
  inference implementation.

Pisces/X-Cell therefore motivates target priors and control-conditioned
adaptation, but neither is executable today. Monitor releases; do not duplicate
experiment 076's data scout or plan around unreleased artifacts.

### 2026-08-23 addendum — SLIM provides a compact STRING prior

The newly released MIT-licensed
`RasmussenLab/SLIM` repository (revision `5a7e9ade5d0a`) provides a
19,502,608-byte STRING-network embedding file and a closed-form bilinear ridge
model for unseen genetic perturbations. The pinned artifact has 21,688 genes,
64 dimensions, and covers **298/300** VCC panel targets; only `EEF1A2` and
`KIF21B` are absent. Its local SHA-256 is
`789416877b8701ef6f800106d26bf7bb97ea8e72744e6ab93e24933a717f247d`.

This is a stronger and much smaller independent target-prior ablation than
ESM2 alone. It does not supersede experiment 082 or authorize combining
embeddings after seeing held-panel outcomes. The leakage-safe next use is a
separately preregistered STRING-only residual model under the same non-panel
fit/tune/test and held-donor contract, followed by a fixed ESM2-versus-STRING
comparison only if both implementations pass their own null and KOLF gates.
The artifact and provenance are recorded under `data/external/slim/`.

## Common leakage and decision contract

These requirements apply to every experiment below:

1. Remove all 300 VCC target-gene **coordinates** from effect vectors before
   fitting, transformation, similarity calculation or retrieval.
2. Create a deterministic hash split of shared **non-panel perturbations** into
   fit/tune/test partitions. No VCC-panel perturbation may select a feature,
   rank, regularizer, blend fraction, stopping epoch or model family.
3. In every outer held-donor fold, exclude that donor from prediction
   construction and from every fitted statistic. Freeze the candidate using
   only the remaining donors and the non-panel fit/tune targets.
4. Open held-donor VCC-panel effects once, only for the final frozen
   effect-level test, and rank every query against the complete shared
   VCC-panel truth gallery.
5. Report per-fold and per-target rank deltas. Bootstrap by target clusters
   (gene family/pathway clusters where available), not by treating correlated
   targets as independent.
6. Promotion requires all of:
   - nonnegative PDS delta in every valid outer fold;
   - macro raw-PDS delta at least **+0.003**;
   - cluster-bootstrap 95% lower confidence bound above zero;
   - observed gain above the 95th percentile of the experiment-specific null;
   - no material collapse in target coverage or an obvious dependence on one
     target family.
7. Only a passing, fully frozen effect-level candidate may receive one
   HCT116 count-level `cell-eval2` confirmation. Use one isolated guarded
   scorer; never run parallel local scorers.

## 1. Cross-donor listwise metric adapter

### Mechanism

Whitening is unsupervised and improves retrieval despite lowering ordinary
cosine to truth. The remaining mismatch is therefore plausibly **metric
supervision**: donor data reveal which low-amplitude directions consistently
distinguish the same perturbation from the other perturbations.

Learn a near-identity, low-rank residual map on the post-whitened prediction:

```text
z' = normalize(z + beta * U(V^T z)),   rank(UV^T) in {4, 8, 16}
```

Optimize a temperature-scaled listwise softmax in which the positive is the
same target in another donor and every other target in that donor is a
negative. Penalize `||UV^T||`, constrain `beta`, and include exact identity as
a candidate. This directly optimizes the ordering underlying PDS rather than
mean cosine or reconstruction error.

This is not another whitener, PCA truncation, CSLS correction or diagonal
gene-weight sweep: it uses same-target cross-donor labels and a listwise loss,
and it is constrained to make only a small residual correction after the
frozen production geometry.

### Exact leakage-safe protocol

1. Load the compact donor effect tables and reconstruct arm D in each outer
   fold using no held-donor data.
2. On non-panel fit targets, form cross-donor pairs
   `(prediction_from_other_donors, held_donor_truth)`. Train ranks
   `{4,8,16}`, `beta in {0.02,0.05,0.10,0.20}`, temperatures
   `{0.03,0.07,0.15}`, plus identity.
3. Select one configuration by macro listwise PDS on the non-panel tune split
   with an inner donor holdout. Require the same sign of improvement in every
   inner fold; otherwise select identity.
4. Run the frozen configuration on the non-panel test split as a pre-gate.
5. If the pre-gate passes, evaluate once on held HCT116, HEK293T and every
   other donor with adequate shared-panel coverage, using the full panel
   gallery and the common contract above.
6. Null: repeat the complete fit with same-target labels permuted within
   response-norm deciles. The real delta must exceed null p95.

### Compute and data

- **Data:** existing compact donor effect tables and frozen arm-D artifacts;
  no new download and no single-cell matrix.
- **Compute:** CPU or Apple silicon; expected peak RAM below 16 GB if features
  are projected/chunked. A single GPU is optional, not required.
- **Time:** hours, not a cloud-scale training run. This is the preferred next
  experiment.

### Promotion gate

Apply the common `+0.003`, every-fold, confidence-bound and null gates. Add two
guards: adapter rank must be at most 16, and median angular displacement from
arm D must remain below 10 degrees. A larger unconstrained transform is a new
model and is not supported by the available number of independent donors.

### Why prior evidence does not refute it

Experiments 017–025 learned an unlabeled covariance geometry. Experiments 044
and 073 use gallery-density distillation, while experiment 075 applies a
query-only reciprocal-neighborhood centering rule. None uses same-target
cross-donor supervision. Experiment 066 routes among fixed arms per target
rather than learning one globally regularized retrieval metric. No prior
experiment directly trained a held-donor-safe listwise objective.

## 2. Biological-prior target-kernel residual

### Mechanism

The incumbent predicts each target mostly from that target's measured donor
rows. It does not learn that perturbations of biologically related proteins
can share a transferable residual. Use a target embedding only to predict the
part systematically missing from the expression-only incumbent:

```text
r(d,p) = truth(d,p) - arm_D(excluding d, p)
r_hat(p) = low-rank ridge/kernel regression from target_embedding(p)
prediction = arm_D + beta * r_hat(p)
```

Start with the released ESM2 perturbation features. STRING or GenePT may be
predeclared later as independent ablations, not combined in the first search.
Fit the residual output in a rank-16/32/64 basis learned only from non-panel
training residuals. This can denoise directly observed targets and provide a
principled fallback for targets missing from one source donor.

### Exact leakage-safe protocol

1. Acquire only `ESM2_pert_features.parquet` and verify its revision and hash;
   do not download the associated 244+ GB expression bundle.
2. Intersect embeddings with donor targets. Standardize embeddings using the
   non-panel fit partition only.
3. In each outer donor fold, construct arm-D residuals on non-panel targets
   with the held donor excluded everywhere.
4. Fit residual PCA ranks `{16,32,64}` and ridge strengths
   `{1e-3,1e-2,1e-1,1,10}` on the fit split. Select rank, ridge and
   `beta in {0,0.02,0.05,0.10,0.20}` on the tune split via inner-donor macro
   PDS. Identity (`beta=0`) is mandatory.
5. Freeze, run the non-panel test pre-gate, then perform the one-time
   full-gallery panel evaluation under the common contract.
6. Report direct-row and missing/partial-source target strata separately.
7. Nulls: permute embeddings among targets within baseline response-norm and
   source-coverage strata; also compare an equal-dimensional random Gaussian
   embedding. The biological embedding must beat both null p95 values.

### Compute and data

- **Data:** existing compact donor tables plus the 595 MB ESM2 parquet from
  the Altos collection. No large expression download.
- **Compute:** CPU is sufficient for kernel/ridge variants after PCA; one
  modest GPU may accelerate linear models. Target peak RAM below 24 GB.

### Promotion gate

Use the common gate, and additionally require a nonnegative delta in both the
complete-source and partial-source strata. Do not promote a result carried
only by targets absent from the incumbent or by one protein family.

### Why prior evidence does not refute it

The failed confidence, donor-weighting and target-routing experiments use
response magnitude, agreement, basal state or arm diagnostics. None maps
target identity to a cross-donor residual. X-Cell's reported prior ablations
make the mechanism credible, but its unreleased weights cannot answer whether
a small leakage-safe target kernel helps this incumbent.

## 3. Guide-level latent on-target decontamination

### Mechanism

Experiment 068 showed that deleting an entire K562 source row when one guide
has a suspected seed off-target removes useful on-target biology and usually
hurts. The unresolved alternative is to decompose guide-level effects:

```text
guide_effect(target, guide)
    = target_shared_on_target
    + learned_seed/guide_component
    + residual_noise
```

Estimate the target-shared latent component across all multi-guide K562
targets, using guide sequence/seed features and robust shrinkage. Replace only
the K562 contribution to arm D with that component; keep single-guide or
unidentifiable targets unchanged. This is source denoising, not hand-curated
editing of four panel rows.

### Exact leakage-safe protocol

1. Stream the already-local Replogle K562 genome-wide H5AD and build
   guide-level pseudobulk effects. Do not materialize a dense cell-by-gene
   matrix.
2. On non-panel targets only, fit a hierarchical model with a target-shared
   latent effect and a low-rank guide/seed component. Hyperparameters are
   selected by holding out one guide per target and predicting its effect;
   target labels never cross the fit/tune/test boundary.
3. Freeze the model. For panel targets, infer the target-shared component from
   their public K562 guide rows using fixed global parameters. This uses no
   held-context truth and no literature lookup.
4. Insert cleaned K562 rows into the otherwise frozen arm-D predictor. Keep
   the existing whitener fixed in the primary test so the only changed object
   is the K562 measurement.
5. Evaluate on held HCT116, HEK293T and independent KOLF truth where coverage
   permits. Those held donors are excluded from every fit and are opened only
   for final full-gallery ranking.
6. Nulls:
   - permute guide-to-target assignments within guide-count and expression
     strata;
   - use guide-count-matched ordinary averaging;
   - apply the model with seed features permuted.
   The learned decomposition must beat all null p95 values.

### Compute and data

- **Data:** the existing approximately 8.2 GB local K562 genome-wide H5AD plus
  guide metadata/sequences; no new atlas download.
- **Compute:** stream/chunk locally only if measured peak RSS stays below
  24 GB. Preferred execution is one isolated 64 GB cloud worker. Persist only
  compact guide pseudobulks and fitted factors.
- **No `cell-eval2`** until the effect-level gate passes.

### Promotion gate

Use the common gate across all applicable held truths. Also require:

- at least 100 evaluable multi-guide panel targets;
- nonnegative delta on both targets flagged by seed evidence and unflagged
  multi-guide targets;
- no gain attributable merely to shrinking every K562 row toward zero.

### Why prior evidence does not refute it

Experiment 068 tested four hard row deletions and explicitly concluded that
whole-row removal was too blunt. PerturbMatch's success shows that processing
multiple source estimates can improve transfer, but it does not separate
guide-specific from target-shared signal inside K562. No experiment through
077 fits a guide-level latent decontamination model.

## Deliberately not selected

- **Finish 065 and 077:** required bookkeeping/scientific closure, but not new
  hypotheses. Experiment 070 is now closed.
- **Experiment 075:** independent experiment-079 confirmation failed its
  confidence and practical-effect gates. Do not retune or promote it.
- **More alpha, spread, CPM, whitening-gamma, donor-weight or arm-routing
  sweeps:** repetitions of mature or closed families.
- **Full control-covariance transport/TTA:** technically distinct from
  experiment 067's diagonal transport, but basal/control conditioning has
  repeatedly failed and the sample size is inadequate for a stable
  high-dimensional map. Reopen only if X-Cell weights/code are released and
  its exact control-only TTA can be reproduced.
- **Pisces as a donor or X-Cell as a benchmark:** high priority immediately
  after a real data/model release; not executable from the current placeholder
  repositories.
- **Another data scout or bulk download:** experiment 076 and the live metadata
  check already establish the current public inventory.

## Recommended order

Run experiment 1 first. It is the cheapest, uses only existing compact data,
and attacks the actual dominant metric with supervised information not used by
the incumbent. Run experiment 2 next if the ESM2 artifact is acceptable to
fetch. Run experiment 3 only on an isolated worker after a streaming
guide-pseudobulk dry run establishes the memory ceiling.
