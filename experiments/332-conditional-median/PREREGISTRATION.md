# 332 — conditional median of a cross-context response

September 9, 2026. Development pilot, not independent confirmation or a
competition candidate. Freeze this protocol, implementation and all inputs
before building new response features or running a biological fit.

## Question and novelty boundary

Can a **shared, identity-free conditional-median predictor** reduce
target-balanced absolute response error in an unseen cell line and on unseen
targets, without losing retrieval, beyond source-only and simple-model
controls?

This extends, rather than erases, earlier negative results. 161 already
regressed a pooled scalar mean residual on donor/recipient features; its
orthogonalized, norm-restored addition was mostly context-free reshaping.
184's basal gains likewise failed their recipient control on CRISPRi. 104
learned sign-preserving signed-credit reweighting. 224 fit gene-space
operators/diagonals and a target-embedding/basal model. 268/289 fit separate
per-target basal/state regressions from eleven lines; 269 routed among
existing actions. None establishes that merely swapping a model library is
new information.

Here the estimand changes to the **conditional median of the full scalar
effect**, normalized by source response RMS, minimizing absolute error
rather than a squared-error residual direction followed by orthogonalization
and norm matching. A pooled model shares one response law across both target
and response genes, with no identities. The new substrate has twelve working
Cas9 lines, not 161's five CRISPRi lines; it is reused development data with
only 82 targets. TabICL supplies a synthetic learned prior, not new biological
labels. Its conditional **mean from the same fit** is a mandatory ablation.
The difference from the mean is measured, not assumed beneficial.

## Data and isolation

Only 268's twelve working half-pseudobulk files:
C4I, HCT15, HKGZCC, KMRC20, KYSE140, KYSE450, ONE58, RPMI7951,
SKGII, SNU8, T3M5, U343. JVE127, MG63, SLR23, UMRC3 remain spent and unread.
No new raw data, reserve, donor bank, challenge response or H1 response.

Use the **82 already fixed targets** in 268's `lolo268_K2.csv`; read that
historical table only for target IDs and known-answer baseline verification.
Two target blocks: sort by SHA256("332-target|" + symbol), alternate positions.
For every outer held line and target block, fit only the other target block
in the eleven other lines. A training recipient's source set excludes both
itself and the outer held line. A test target's source responses in the
training lines are legitimate build-time inputs, not training labels.
There are 24 fits per learned arm.

Full log-profile frame: log1p(50,000 * counts / total allowed counts).
Controls are the ten OR knockouts, as in 268. Gene measurement/finite mask
uses controls only; expression >=5 mean CPM uses **outer training controls**
only. Remove all Broad targets, all 300 A/B/C targets, and all 2025 H1
panel target coordinates before normalization/features/scoring. H1 target
IDs only come from 323's registered panel metadata; no truth is read.
Keep **512 genes**, ordered by SHA256("332-gene|" + symbol), from each
outer-training eligible axis. This is a limited-coordinate pilot; no
full-transcriptome or A/B/C score claim.

Labels and source-row RMS use the full eligible axis before slicing.
For a source set S and target t, scale is
max(sqrt(mean over sources and full-axis genes of response squared), 1e-3).
Each training label is recipient response / its source-only scale.
No held response enters a filter, transform, normalization, feature,
training label, stopping rule, model choice or training sample selection.
Training standardization is fitted on training features only.

Per split, sample **2,048** distinct (recipient, target, gene) rows, balanced
over the 11 recipients by deterministic seeded sampling without replacement.
Predictions cover every held target and all 512 selected genes, not just
sampled pairs. Donor features may use all source responses for that query
target. Cross-target independence is therefore out-of-label, not a claim
that the model has never seen the query target's donor effect.

## Features and arms

Fourteen scalar columns, all source-only except the last four:

1. source mean / scale;
2. source median / scale;
3. source SD / scale;
4. fraction of source coefficients >0;
5. source minimum / scale;
6. source maximum / scale;
7. log1p(source scale);
8. mean source control log expression of the response gene;
9. SD of that source control expression;
10. mean source covariance of scaled effect and control expression;
11. recipient control log expression;
12. recipient minus mean source control expression;
13. column 1 times column 12;
14. column 1 divided by (1 + column 11).

No target/gene/line ID, assay batch, genotype, depth, measured perturbation
cell count or recipient perturbation outcome feature.

Fixed arms, no dose/ensemble/feature sweep:

- Source mean baseline; zero-effect baseline.
- **Primary:** TabICL posterior median, converted back with source-only scale.
- Same-fit TabICL posterior mean (estimand ablation).
- TabICL median with only first ten source columns (recipient ablation).
- Same full fitted model, median with a **fixed wrong recipient**: next
  working line in cyclic order. Source features remain unchanged. This is
  one deterministic negative control, not an estimated permutation null.
- TabICL median with training target-label blocks deranged cyclically
  within each training recipient; label includes the alternate target's own
  scale. This preserves gene/line marginals, breaks target correspondence.
- Feature-matched standardized Ridge(alpha=10), including intercept.
- Feature-matched HistGradientBoostingRegressor(loss="absolute_error",
  max_iter=100, learning_rate=.05, max_leaf_nodes=15, min_samples_leaf=20,
  l2_regularization=1, early_stopping=False, random_state=332).

TabICL version/checkpoint/environment are exactly 331's. CPU, four threads,
one estimator, normalization none, AMP/cache/FA3 off, auto-download false,
seed 332; inference batches of 256 rows. Offline and telemetry-disabled.
No clipping, output blending, orthogonalization, row-norm restoration or
post-result tuning. Predictive intervals are not claimed calibrated.

## Tests, outcomes and decision

Before fitting: pure synthetic tests for split/source isolation, source-only
features unchanged by recipient substitution, exact baseline/zero metrics,
tie-handled retrieval, deterministic sample/target splits. A poisoned held
response must not change any training/query feature or training label.
Reproduce all 268 historical per-target baseline cosines/retrieval at <=1e-6
in its original full-axis frame, separately from the stricter new axis.
Checkpoint/input drift or failed known-answer check aborts.

For every (line,target), report:

- mean absolute error / source RMS (primary loss; **not DE NMAE**);
- mean squared error / source RMS squared (not jackknife-corrected VCC MSE);
- own-target cosine;
- PDS-like retrieval over all **82 truths** on the fold's fixed 512 genes,
  after joining the two cross-fitted target blocks;
- orthogonal-residual cosine and half-A/half-B reproducibility as descriptive
  diagnostics; source baseline and oracle self-retrieval as instruments.

Uncertainty: 10,000 paired target bootstraps (each target carries all lines)
and separately a crossed line-by-target bootstrap, seed 332. Also report
line-wise changes and leave-one-line-out macro differences. These do not
create independent batches; line/assay batch confounding remains.

All of the following are required to earn a **separately registered larger
effect/count-stage test**, not a submission:

1. Primary absolute-error improvement over source mean, zero, source-only
   TabICL, wrong-recipient, shuffled-label and absolute-error boosting is
   positive with both bootstrap lower bounds >0.
2. Primary median beats same-fit mean in absolute error with target-bootstrap
   lower bound >0. Otherwise the new estimand is not supported.
3. Primary retrieval gain over source mean >=.005 at the point estimate,
   with both lower bounds >=-.005; cosine loss no worse than -.005 at the
   point; squared error no worse than source mean at the point.
4. At least 8 of 12 lines improve primary absolute error versus source-only
   TabICL, and every leave-one-line-out macro difference remains positive.
5. Every integrity/isolation gate passes.

Report every arm whatever the result. If any gate fails, no advancement or
retuning of this construction. An interesting descriptive effect is not a
passed gate. All six competition members could change with these
predictions; **none is priced here**. No count realization/promotion is
authorized by a source-stage point estimate.

## Resource and authority bound

One heavy local process; launch gate >=40% system free, free+inactive >=20 GiB,
swap free >=400 MiB. Max RSS 12 GiB, scratch 4 GiB, biological phase <=2 hours
and per split <=5 minutes; terminal failure/checkpoint on exceeding a bound.
Save split-level predictions/receipts; never rerun a completed biological
split merely because its result is poor. An execution bug needs a disclosed
amendment and preserves prior outputs.

No built-in subagents, Claude worker, paid API/cloud, production edit,
reserve access, package, submission, external message or push.
