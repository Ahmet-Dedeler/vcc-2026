# 340 — sparse response-neighborhood transfer

September 9, 2026. Freeze before new response fits or evaluation.
Source-direction pilot only; whole-recipient and target-block exclusion.

## New question, bounded against prior work

Can source responses at FUNCTIONALLY NEIGHBORING OUTPUT GENES predict a
held line's response beyond the source value at that output gene?

334/336 supplied same-coordinate source summaries plus target/gene
annotations; they did not supply the source target's response at nearby
output coordinates. 123 predicted from a target STRING embedding and a
response PCA basis. 108/190 directly propagated signed regulatory priors.
224 fitted an unrestricted gene-space operator and lost retrieval.
This pilot tests a specific sparse neighborhood restriction on such
off-diagonal prediction, with learned edge signs and a diagonal baseline.
It is not a new architecture, information-free renamed rerun, or claim
that greater capacity should help. Source-correlation and random graphs
are mandatory comparators; a graph-name advantage is not assumed.

## Data and graph identity

Reuse 334's six distinct CRISPRi lines, exact 1,024 targets, 512 genes,
unit-RMS raw response frame, and two complementary target blocks.
No H1 response or reserve. Both challenge/H1 panels and the full source
target intersection stay out of the response coordinate axis.

STRING: 100's pinned `string-pca.npz`, ONLY its original normalized
64-dimensional `unit` and availability fields, not its fitted PCs.
Exact gene lookup. For each available response gene select up to eight
OTHER available genes with positive embedding cosine, largest first,
gene-index tie-break. Missing or unsupported neighbors are zero-padded.
Require >=90% of the 512 output genes to have eight neighbors; otherwise
stop at metadata feasibility without response fitting.

This is a k-nearest graph in a released STRING EMBEDDING, not an assertion
that every edge is a documented physical protein interaction.
Freeze graph identities before response access.

Random graph: same per-output degree, uniformly sample distinct other
coordinates using seed 3400909; fixed before response access.
Correlation graph: each fit selects the same per-output degree by largest
absolute Pearson correlation among its inner source-consensus columns,
using TRAINING targets/lines only, no recipient labels; gene-index
tie-break. Zero-variance source columns have correlation zero.

## Model and controls

For each outer held line and 512-target test block:

- Training consists of five inner held-line pairs. Each inner recipient
  contributes all 512 TRAINING targets; its source mean uses the other
  four eligible lines. Outer line absent from every training/source object.
- Fit exactly 334's per-gene diagonal on those pairs; require its query
  prediction and the source mean to reproduce all twelve saved replays.
- For each output gene, input is its own source-consensus coordinate plus
  its eight graph-neighbor coordinates. Response is held-training-line
  direction minus the diagonal prediction.
- Divide each input feature by its training RMS (floor 1e-12), no centering
  or intercept. Solve ridge with coefficient penalty `n_training * I`
  (lambda=1 in mean-square units), prior residual coefficients zero.
- Query the five-donor mean. Because this is a fixed LINEAR predictor, it
  must equal averaging the five four-donor-subset query predictions.
  Baseline diagonal plus fitted residual is the final direction.

Seven arms: source_mean; diagonal; own_refit (same ridge, all neighbor
slots zero); graph (primary); random_graph; correlation_graph;
labels_shuffled (same graph, training recipient responses deranged by
whole target vectors within each recipient using one consistent
permutation; diagonal baseline remains the correctly trained comparator).

No intercept, clipping, output renormalization, blend, grid or early
selection. All arms predict every held target and coordinate.

## Tests and frozen advancement

Pure checks: sparse feature indexing/zero padding; no self/duplicate
neighbors; exact ridge known answer; no target/recipient leakage;
zero residual -> unchanged diagonal; zero input -> zero output;
linearity/subset identity; graph/label derangements; input immutability;
synthetic correlated-neighbor example.

Raw-loader historical checks and all twelve source/diagonal replays
must pass. Poison the outer recipient and verify identical fitted
models/query predictions. Fit every feature scale and correlation graph
inside the training block. Save models, indices and prediction hashes.

Endpoints: per-target full 1,024-gallery retrieval, own/decoy cosine,
unit-RMS squared/absolute error. No VCC-score names or board conversion.
10,000 paired target and crossed line/target bootstraps, seed 3400909,
conditional on fits. All competitors share folds/targets.

Advance the PRIMARY graph only if ALL:

1. Retrieval gain >=.005 over source mean.
2. Both retrieval 95% lower bounds >0 versus source_mean, diagonal,
   own_refit, random_graph, correlation_graph and labels_shuffled.
3. Own-cosine point loss <=.005 versus source mean AND diagonal.
4. Squared-error point no worse than source mean AND diagonal.
5. Retrieval versus diagonal positive on >=5/6 lines, every leave-one-line
   average, and each frozen lower/higher median-cell-number half.
6. All metadata, known-answer, isolation and numerical gates pass.

A pass earns only a separately preregistered actual-count test on the
current family. A failure closes this fixed sparse-neighborhood ridge,
not every possible graph model. Do not tune K/lambda/graph/axis afterward.
All six competition members are unpriced and potentially movable.
These reused data do not constitute independent biological confirmation,
and the cell-number split does not adjust for UMI depth or study.

One serial CPU process, six threads, 8 GiB RSS ceiling, 60-minute total
budget and repository memory launch gate. No model download, paid API,
cloud, worker, reserve, production change, package, or submission.
