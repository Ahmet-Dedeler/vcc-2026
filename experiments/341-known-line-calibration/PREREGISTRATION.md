# 341 — prospective known-line, held-target calibration

September 9, 2026. This is an explicitly DIFFERENT information regime:
the recipient's public responses to other training targets are available.
It is not an unseen-line or cross-study result, does not replace 323, and
does not relax or relabel any prior failed gate.

## Question and closest priors

Does access to the recipient line's own non-evaluation target responses
improve downstream sign prediction beyond an otherwise identical
cross-line model, while retaining target discrimination?

195 measured a held-line boundary for learned corrections. 224 reported a
same-line diagonal as a leaky ceiling, not this held-target sign test.
117/118/141 used same-target Pan integration; 129 used Pan target-neighbor
geometry. 340 tested sparse off-diagonal regression without recipient
labels and failed retrieval. This experiment's NEW information is own-line
training labels; capacity alone is not its hypothesis.

The Claude review's phrase "340 unchanged" was inaccurate for its proposed
full 512-coordinate ridge. Here full-rank ridge is explicitly shared by
the same-line and cross-line arms. Their comparison isolates label regime
within this experiment, not against the old sparse architecture.
Any later application to context A would require an inferred Jurkat
identity, cross-study transfer, and a separate current-family count test.
None is established by this development pilot.

## Fixed frame and leakage boundary

Exactly 334's six banks, 1,024 targets, 512 genes, two complementary
512-target blocks, unit-RMS raw directions and metadata depth halves.
All A/B/C and H1 panel targets/coordinates and the source target
intersection remain excluded as in its frozen selection.
No new response dataset, reserve, H1 reference or challenge outcome.

For recipient r, eligible sources are the OTHER five lines. For each
training recipient s among those five, X_s is the four-source consensus
excluding BOTH r and s, on the TRAINING target block. Stack the five X_s.
They are identical feature rows for every arm below.

- Cross-line labels: the matching s training responses, as in 334.
- Same-line labels: repeat r's TRAINING responses for every X_s.
- Shuffled labels: cyclic derangement of r's whole training-target
  vectors, consistently across all five X_s.

The repeated targets are feature views, not independent experimental
replicates. Never use r's TEST-block response in a feature, scale, label,
operator, baseline fit or selection. Poison it and require identical
predictions. Query with the five-source mean; a fixed linear map of it
must equal the mean of the five four-source-subset predictions.

## Fixed estimator and arms

For any label regime Y, fit the same per-gene diagonal:
g=(sum(X*Y)+median(sum(X*X)))/(sum(X*X)+median(sum(X*X))).
No positivity clipping; do not assume this diagonal cannot flip signs.
The cross-line diagonal must reproduce all twelve 334/340 saved replays.
For that replay use 334's per-recipient accumulation order; the stacked
formula is algebraically identical but can differ at floating precision.

Fit full-coordinate residual ridge to Y-X*g. Input columns divided by
their TRAINING RMS (floor 1e-12), no centering/intercept, coefficient
penalty n_rows*I (lambda=1 in mean-square units). No rank truncation,
parameter search, output normalization or selected blend. The same
feature scale/penalty/512 coordinates are used in every full model.

Arms: source_mean, crossline_diag, crossline_full, same_line_diag,
same_line_full (PRIMARY), labels_shuffled_full. The shuffled arm fits
BOTH its diagonal and residual from shuffled labels.
No wrong-line model whose label secretly appears as one of its input
donors: crossline_full is the matched other-line-label control.

## Endpoints and frozen gates

Reuse full 1,024-gallery retrieval, own/decoy cosine, unit-RMS squared/
absolute error. Add target-macro |truth|-weighted sign accuracy, and
the same endpoint within each target's top 64 |truth| coordinates
(stable gene-index tie breaking). Zero predictions are incorrect for
nonzero truth. Every target must have positive truth mass.
The top-64 set is a response-strength proxy, NOT a reference DE/BH set.

10,000 paired target and crossed line/target bootstraps, seed 3410909,
conditional on fitted models. No pooling all target/gene pairs into
one easier endpoint. Reused development data, not fresh confirmation.

All advancement requirements:

1. Primary full-axis weighted sign gain >=.010 over crossline_diag,
   with BOTH 95% lower bounds >0.
2. The same full-axis sign contrast has both lower bounds >0 versus
   crossline_full, same_line_diag and labels_shuffled_full.
3. Top-64 weighted sign contrasts versus crossline_diag and crossline_full
   both have positive lower bounds under BOTH bootstraps.
4. Retrieval versus crossline_diag and crossline_full has both lower
   bounds >-.005 (a noninferiority guard, not a claimed retrieval gain).
5. On Jurkat alone, target-bootstrap weighted sign lower bounds >0
   versus crossline_diag and crossline_full. This is a prespecified
   recipient check, not an A/B/C forecast.
6. Full-axis sign point contrast versus crossline_diag is positive
   on >=5/6 lines and both frozen cell-number halves.
7. All known-answer, isolation, finiteness and numerical checks pass.

A pass earns only separately registered CROSS-STUDY and actual-count
validation. It cannot clear a package, D/E/F default, or submission.
All six competition members are unpriced and potentially movable.
Failure closes this fixed known-line ridge construction; no retuning.

## Checks and resources

Pure tests: ridge/diagonal known answers, identity residual, linearity,
zero input, consistent label derangement, exact target isolation,
target-macro sign arithmetic/top-64 ties, and input immutability.
Reproduce historical raw loaders and every source-mean/cross-line diagonal
prediction. Independently verify fitted normal equations and endpoint
arithmetic after the run.

One CPU job, six threads, 8 GiB tree RSS ceiling, repository memory gate,
45-minute run budget. No paid API/cloud, worker-driven execution, reserve,
production edit, package or submission. Codex runs and verifies.
