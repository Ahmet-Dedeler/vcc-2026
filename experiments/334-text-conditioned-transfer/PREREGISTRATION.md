# 334 — direct primary-description features for source-conditioned transfer

September 9, 2026. Freeze this protocol, code, selected identities and input
hashes before reading response values for this experiment. No parameter
sweep, submissions, reserves, H1 response labels or production writes.

## Question and novelty

Can primary gene-description text improve held-line, held-target response
prediction **on top of same-target source effects**, beyond source-only
features, dimension-matched STRING features and a per-gene diagonal?

087/123 tested target STRING maps, 100/101/104 pairwise sign predictors,
132/161/195 source-conditioned regressors, and 224 a learned transfer
ladder. 224's diagonal was positive. 332's conditional median lost
retrieval and did not beat its recipient-free control. This is a new
information test, not a claim that boosted trees are a new architecture.
The GenePT pretrained vectors from 333 are quarantined and NEVER loaded.

## Fixed data and frame

Six already-used CRISPRi banks, one per distinct line:
K562-gwps, Orion HCT116, Orion HEK293T, RPE1, Nadig Jurkat and Nadig HepG2.
No K562-essential second library in model fitting; no H1 source/label.
These are reused development data, not independent confirmation.

Take 1,024 SHA-256-ordered targets from the six-bank intersection after
excluding A/B/C and the 323 H1 evaluation target list. Response axis:
512 hash-ordered canonical shared genes, excluding the entire source
target intersection and both evaluation panels. Do not choose coordinates
from response outcomes. Exact identities are frozen in `selection.json`.

Read raw `delta.npy` values, with common `target_sum=50000` confirmed from
metadata. No across-target centering, no old precentered response cache
in model fitting. Divide each response row by its own coordinate RMS:
labels are unit-RMS directions, not calibrated log fold changes.
Source features use only eligible donor directions, raw donor row RMS and
donor controls. A held response is never a feature or training donor.

Each outer line has two complementary hash-defined target blocks:
512 train / 512 test. Training examples hold one of the other five lines
out as recipient and use the remaining four as donors. Query each held
target with each of the five four-donor subsets; average predictions.
Thus both training and query use four donors. The baseline average of
the subset means must exactly reproduce the five-donor mean.

For each fit, use 32,768 unique pairs sampled approximately equally over
the five training recipients (max count difference one). Within a recipient
sample uniformly from its 512 training targets x 512 coordinates.
All learned tree arms share samples, responses, splits and hyperparameters.

## Text identity and construction

Use ONLY `NCBI_summary_of_genes.json` from 333's verified archive
(member SHA `3db721b7cfbc35795428c6a7ab4e8a9c59fb5f5ccb8c64af50dd3c97c2d4de2e`).
Exact primary-key lookup; no aliases, stripping keys, case conversion or
manual gene annotations. A nonempty value must begin `Gene Symbol KEY`
with that exact key and token boundary; otherwise treat it as missing
and report it. Remove that identifying prefix, not other biological text.
Require usable descriptions for >=95% selected targets and >=90% genes.
Missing descriptions become zero features; never drop their outcomes.

TF-IDF: lowercase, English stop words, word 1/2-grams, sublinear TF,
min_df=2, max_features=12000, L2 norm. Fit vocabulary and IDF on training
target descriptions plus all selected response-gene descriptions.
Truncated SVD: 16 components, 7 iterations, fixed seed. Fit only those
same documents. Transform held-target text; do not fit on it.
L2-normalize the 16-dimensional representations.

STRING comparator uses the public, response-free `unit` vectors in
100's hash-pinned feature artifact. Refit a 16-component PCA using only
training targets plus response genes. Missing embeddings stay zero;
normalize projected vectors. Do not use the artifact's prefit PCs.

Pair feature blocks are target vector, response-gene vector,
elementwise product and their dot product (49 columns).
Source-only gets 49 zero columns. Main and control fits have the same
tree limits, but zero padding does not create equal information.
Text permutation independently cyclically deranges whole target vectors
and whole response-gene vectors after projection, preserving each
representation bank. This tests correspondence, not a sampled null law.
The additional label control deranges whole training-target response
vectors consistently across recipients, using only training-block labels.

## Fixed models and arms

One `HistGradientBoostingRegressor` residual model per fitted arm:
squared_error, 150 iterations, learning_rate=.05, max_leaf_nodes=15,
max_depth=4, min_samples_leaf=64, l2_regularization=10,
early_stopping=False, fixed random_state. Predict label minus source
coordinate mean, then add source mean back. No recipient identity,
recipient basal feature or recipient-conditioned claim.

Arms: source mean; zero; training-only per-gene diagonal; source-only
residual tree; source+STRING; source+text (primary); same tree with
deranged text pairs; same text tree with deranged training labels.
The diagonal uses all inner training-target pairs, not just the sample:
gain=(sum(P*Y)+median(sum(P*P)))/(sum(P*P)+median(sum(P*P))).
This gives the simple comparator MORE training rows, explicitly.

## Verification and endpoints

Pure tests: recipient/source isolation, no test labels in fitting,
balanced unique sampling, subset-mean identity, whole-vector derangement,
zero/self retrieval, diagonal known answer, direct-text identity checking,
held-text-vocabulary isolation and a synthetic residual learning test.
Before fitting reproduce the five relevant bank arrays in historical
224 S5 by raw loading and its historical centering/float32 conversion.
This historical check is separate and does not preprocess the pilot.
Nadig HepG2 additionally gets a direct-index loader equality check.
Poison held-line responses and prove identical training/query features,
training labels and diagonal. Fail closed on nonfinite data or hash drift.

Report per-target unit-RMS absolute/squared error, own cosine,
full 1,024-target gallery retrieval and decoy cosine. These are NOT VCC
MSE, NMAE or an A/B/C score. No borrowing a board transfer multiplier.
10,000 paired target resamples and crossed line/target resamples;
conditional on fitted models, not a complete training bootstrap.

Advancement requires ALL:

- primary retrieval improvement >=.005 vs source mean;
- retrieval intervals strictly above zero for primary vs source mean,
  source-only, STRING, text-shuffled and diagonal, for BOTH bootstraps;
- cosine point loss no worse than .005 vs source mean and diagonal;
- squared-error point no worse than source mean;
- retrieval vs source-only positive on >=5/6 lines, in every leave-one-line
  average, and separately in lower/higher median-cells-per-perturbation
  halves (metadata-only ordering, fixed before response access);
- all integrity/known-answer gates pass.

Depth split is a coarse cell-number sensitivity check; it is NOT
UMI-depth adjustment. Study/chemistry and line are confounded; no
study-independent or depth-adjusted biological claim is permitted.
Report all six competition members as unpriced and potentially movable.
Passing only earns a separately frozen current-family actual-count test.
Failure closes this fixed construction, not all text representations.

## Resources and stop

One guarded CPU process; six BLAS/OpenMP threads, max 8 GiB RSS,
minimum 40% system free and repository launch gate. Budget 90 minutes.
No test-result inspection or retuning between folds. Save fitted-model,
prediction, input and selection hashes; terminal report after all folds.
At numerical/integrity failure stop and preserve artifacts. A recovery
needs an explicit new freeze, not silent overwriting.
