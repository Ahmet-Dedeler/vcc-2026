# 329 — Joint-CD4 estimation through current-family H1 counts

September 9, 2026. Freeze this protocol, code, tests, fitted models and input
hashes in git before constructing candidate coefficients or scoring them.

## Question and evidence boundary

328's canonical three-condition covariance mixture passed its source-density
gates after the explicitly post-result numerical continuation. Does applying
its posterior means improve the current cap15 family's emitted H1 predictions?
The original source pilot remains a numerical NO-GO; the separate recovery
result earns this experiment, not a submission.

Only 328's recovered, source-only fitted models are used. No refitting, component
selection, new donor data, SE calibration, H1 response fitting or dose sweep.
The 328 fit excluded H1 and A/B/C targets by symbols and mapped IDs; it also
excluded their response coordinates. Its development targets are not refitted.
Posterior inference on CD4 observations of H1 targets is donor inference, not
access to the recipient H1 perturbation outcomes.

This is a reused H1 development panel, not independent confirmation. Its
reference, anchors and previous results have already been used extensively.
No inference here establishes cancer-context transfer or submission clearance.

## Fixed construction

Use the H1 compact CD4 banks from experiment 273, in fixed order Rest, Stim8hr,
Stim48hr. Align by exact target and gene symbols. For every gene/target with
all three source rows, three finite log2FCs and three finite strictly positive
reported SEs, replace the three coefficients by their fitted-mixture posterior
means. Process batches of 8,192 triples to bound memory.

- Apply the model to all mapped **non-H1-panel response genes** in these banks,
  not only the 256 genes used to fit the gene-agnostic prior.
- Convert posterior means to float32, then use the extractor's existing
  float32 `ln(2) * log2FC` conversion. This is not log-pseudobulk reconstruction.
- Missing condition, invalid triple, or excluded response gene: retain the
  original three source deltas. Do not treat absence as an observed zero.
  If a source logFC is nonfinite, the original extractor's stored zero delta
  remains the fallback. No support expansion or imputation.
- Recompute the **same production preprocessing** on each modified bank:
  per-bank centering and gene normalization, then mean3 and complete-row freeze.
  The preprocessing rule stays fixed; its estimated means/scales change with
  the coefficients. Thus raw fallback is not a final prediction identity.
- Recompute the same k=5 target-row-space deflation from the resulting source
  residual, then the same f=.34 sizing and global preserve step. Do not freeze
  a baseline projector or change to a gene-space projector in this experiment.
- Keep the main-donor/PerturbMatch cap15 stage fixed, alpha .982, corrected
  H1-excluding whitener, gamma .65, H1 union OFF and mass fix unset.
  Use 323's unchanged actual-count emitter and scorer.

Neither compact `b_ctrl` nor DESeq2 baseMean is a control-only expression
estimate; they are not used for basal matching or SE/noise calibration.

## Arms, in fixed order

1. **cap15**: identity coefficient conversion and original preparation.
2. **joint**: the only candidate, full posterior means under 328's recovered
   99-component fit.
3. **diagonal**: identical path, 328's eight-component isotropic/zero fit.
4. **shuffled_fit**: identical path, 328's same-capacity 99-component fit
   trained on target-misaligned conditions. Inference uses the real aligned
   triples; it is the fitted prior, not inference labels, that differs.
5. **shuffled_targets**: fixed complete derangement of joint's final target
   rows, seed 20260909. This deliberately breaks target identity and is a
   harness negative control, not a candidate.

328 component grids and weights are immutable. No alternate interpolation,
SE scaling, selected gene set, amount of deflation or CD4 fraction may be tried
under this registration.

## Known answers and construction gates

Pure tests cover exact extractor conversion including nonfinite values;
batch versus unbatched posterior; isotropic closed form and zero prior;
missing-condition/invalid-SE fallback; excluded-gene identity; support and
target alignment errors; probability/component validation; deterministic
target derangement and geometry identities.

Before any candidate score:

- Every compact bank's identity coefficient conversion must equal its stored
  delta exactly; the newly materialized identity banks must reproduce 323's
  prepared CD4 arrays (mean, votes and counts), baseline U hash and final hash.
- All arms have the same CD4 votes and counts: 80 carried targets, 77 complete;
  incomplete rows have zero residual/support before and after deflation.
- H1 target coordinates never enter the CD4 residual or dense whitener.
- Corrected operator SHA is
  `a46fd3239a74752d5c72f5c6ae7169d354793f04f7c22b9b3fd39fe6f74b9741`.
- Identity final counts and panel-total digests match 323 cap15 exactly;
  all six metrics reproduce at unchanged 323 tolerances.
- Posterior estimates and all final effects are finite. Refuse nonconverged
  fits or changed component specifications. Input manifest changes fail closed.

The fitted-model and source-bank bytes, the production imports, 323 instrument
and prior-stage manifests are hash-bound. No frozen shared experiment is edited.

## Readouts and frozen advancement rule

Score all six official members from actual 400-cell emitted blocks for all
85 targets. Report raw values and H1-anchor scaled values separately.
Record raw MSE and the local clamp; a zero scaled MSE difference cannot
hide a raw profile cost in the report.

Primary contrasts: joint minus cap15, diagonal, and shuffled_fit. Use 10,000
paired target bootstrap draws (seed 20260909), fixed gallery/anchors/models.
Recompute panel-ratio MSE and the global jackknife claim cap in each draw.
Record source, composed and final emitted pseudobulk sign changes and norm
ratios; sign changes are not automatically sign corrections.

Further multi-context development is earned only if:

1. all known-answer/construction gates pass;
2. joint's H1-local overall lower 95% bound is positive versus **each** of
   cap15, diagonal and shuffled_fit;
3. no clearly adverse primary profile-member interval: joint-cap15 PDS upper
   bound is nonnegative, raw MSE lower bound is nonpositive, and raw NMAE lower
   bound is nonpositive;
4. deranging joint's targets reduces raw PDS by more than .10 at the point
   estimate.

Failure yields a scoped negative/unresolved result, not retuning until it passes.
These development intervals are conditional and not multiplicity-adjusted
biological discoveries. Even all gates passing does not authorize packaging,
production promotion or submission.

## Resources and authority

First construction and all five scores within 90 minutes of launch.
One heavy process at a time, ≤22 GiB tree RSS, ≤5 GiB scratch, repository launch
memory gates. Save per-stage execution receipts; guard command-mode logs are
not JSON reports. Stop on drift, nonfinite output, construction or scorer failure.
No new reserve, external responses, paid resource, production change, package,
upload, submission, public disclosure or push.
