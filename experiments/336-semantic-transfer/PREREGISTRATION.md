# 336 — pretrained biomedical representation versus direct word-frequency text

September 9, 2026. Biological fit is prohibited until 335's terminal
engineering gates pass and this protocol/code/input manifest are committed.
334 remains terminal negative. This is not a parameter recovery for 334:
the new input is a fixed pretrained PubMed semantic representation of the
same source descriptions, constructed locally in 335.

## Matched experiment

Reuse 334's exact six lines, 1,024 targets, 512 genes, two target blocks,
32,768 sampled training pairs, four-donor construction and all raw source
features. No new labels, target selection, context input or model capacity.
This is a repeated-development test; previous outcomes informed the
question, so it is not independent confirmation.

Fit a 16-component full-SVD PCA of 335's 768-dimensional vectors on
training-target plus response-gene rows only, exactly as 334 handled
STRING. L2-normalize projections and keep unavailable rows zero. The
pair block remains target/gene vectors, elementwise product and dot
product: 49 features. Apply the unchanged 334 residual tree and seeds.

New arms: semantic (primary), semantic with independent whole-target
and whole-response-gene cyclic derangements, semantic with the same
whole-training-target label derangement used in 334. Permutations use
334's original seeds, not newly selected favorable shuffles.

Comparators are all eight saved 334 arms, unchanged and hash-verified.
Refit source-only with identical training rows/model seed and query
construction on every split. Require every replayed prediction to match
334 within 1e-12 before accepting the split. Do not replace the old
comparator with a newly tuned one.

## Checks and decision

Require exact canonical row/available-mask/body-hash agreement with 335,
its `ENGINEERING_FEASIBLE_ONLY` result, all registered hashes, and
nonfinite/shape checks. Pure tests verify projected-feature fit exclusion,
missing-vector behavior, permutations and old-model settings. 334's
independent raw-training/held-response audit is hash-bound evidence for
the reused features/labels; this script must never use test labels in a fit.

Report all original endpoints and both 10,000-draw bootstrap types.
All six VCC members remain unpriced and potentially movable.
Primary must satisfy every gate:

- retrieval gain >=.005 vs source mean;
- both target and crossed line/target 95% lower bounds above zero vs
  source mean, source-only, diagonal, STRING, TF-IDF text, and shuffled
  semantic features;
- cosine point no worse than -.005 vs source mean or diagonal;
- squared-error point no worse than source mean;
- primary-source-only retrieval positive on >=5/6 lines, all
  leave-one-line averages, and both original median-cell-number halves;
- all engineering/replay gates pass.

The cell-number split is not a UMI adjustment or an independent-study
test. Whole-target derangements are falsifiers, not a sampled permutation
null distribution. Failure closes this fixed semantic-source construction;
do not switch semantic encoders or retune until a positive appears.
Passing earns only a separately registered actual-count test.

CPU only, six threads, <=8 GiB tree RSS, 45-minute fit budget, repository
memory launch gate. No new reserve, network, paid resource, production
edit, package, worker, public disclosure or submission.
