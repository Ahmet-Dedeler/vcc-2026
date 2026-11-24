# Experiment 138 — DLD1 replicated sparse-sign gating

Status: **implemented and source-safe tested; no development or untouched-test
outcome has been opened by this experiment**.

Date: **2026-08-30**

## Question

Can the reproducible, sparse part of the DLD1 Low1/Low2 response improve the
exact Experiment 132 same-target source-mean comparator without replacing its
strong cross-context geometry?

Experiment 116 found low whole-row Low1/Low2 cosine but very high split
retrieval PDS. Experiment 138 therefore does not add DLD1 as a donor, blend its
response vector, learn target neighbours, or alter the response geometry.
DLD1 supplies only target-specific, split-replicated coordinate signs used to
reweight the incumbent.

## Exact inherited evaluation surface

This experiment imports Experiment 132's implementation and binds its v2 plan:

- development targets: the exact 613 rows in
  `experiments/132-response-decomposition-ridge/work/preregistered-freeze-v2/selection_targets.txt`;
- untouched test targets: the exact 410 rows in
  `experiments/132-response-decomposition-ridge/work/preregistered-freeze-v2/test_targets.txt`;
- response axis: the exact 4,992 genes in
  `experiments/132-response-decomposition-ridge/work/preregistered-freeze-v2/response_genes.txt`;
- outcomes: K562-GWPS, K562-essential, RPE1 and Nadig Jurkat;
- comparator: Experiment 132's raw, uncentered, equal source-context mean with
  strict biological-group exclusion. Both K562 screens are excluded whenever
  either K562 screen is scored.

The 613-target development set contains 140 DLD1-supported targets. The
410-target test contains 94. Every other row is exact incumbent fallback.
All metrics use the complete 613- or 410-target gallery, including fallback
rows.

## DLD1 proposal

The proposal source is the verified Experiment 116 compact artifact:

- `transfer_lfc_low1.npy`;
- `transfer_lfc_low2.npy`;
- `transfer_measured_low1.npy`;
- `transfer_measured_low2.npy`;
- canonical gene and target axes.

For target `p` and response coordinate `g`, a coordinate is actionable for
`k` only when:

1. it is measured in both Low1 and Low2;
2. both values are nonzero and have the same sign;
3. it is among the top `k` absolute values in Low1 and independently among
   the top `k` absolute values in Low2, with response-axis order breaking
   ties.

For incumbent value `h[p,g]` and DLD1 consensus sign `s[p,g]`:

```text
match:     h' = h * exp(+beta)
mismatch:  h' = h * exp(-beta)
```

Non-actionable coordinates are initially unchanged. The complete supported
row is then rescaled to preserve the incumbent L2 norm. Positive
multiplication cannot flip a sign. Targets absent from DLD1 are copied
bit-for-bit and are never included in a target-label null.

The row-L2 restoration is deliberately neutral to both development endpoints:
cosine PDS and signed credit are invariant to positive per-row scale. It keeps
the candidate's overall magnitude bounded for later scientific interpretation;
it is not counted as evidence for this experiment.

Candidate family:

```text
identity
k = 50, 100, 200
beta = 0.025, 0.05
```

Exact ties select the least complex arm: identity first, then smaller `k`,
then smaller `beta`.

## Development procedure

The 613 targets remain in Experiment 132 order and are divided into five
contiguous folds of sizes 123, 123, 123, 122 and 122.

For each held fold:

1. select an arm on the other four folds using an inner-only PDS gallery, so
   held-fold truth cannot influence selection as a decoy;
2. evaluate that frozen fold-specific arm on the held fold using the complete
   613-target gallery;
3. retain signed-credit values only for the matching query rows.

After nested evaluation, one final arm is selected on all 613 development
targets. It is frozen before any 410-target value may be read.

The target-label null deranges DLD1 proposal rows among the DLD1-supported
targets. The response-coordinate null deranges all 4,992 proposal coordinates,
moving Low1, Low2 and both measured masks together. Each null has 49 frozen
draws and repeats the complete five-fold arm-selection procedure.

Uncertainty uses 20,000 target-block bootstrap draws. Blocks use the exact
Experiment 132 construction: connected components of reciprocal five-nearest
DepMap neighbours at cosine at least 0.40.

Development passes only when all frozen checks pass:

- a nonidentity final arm is selected;
- PDS and signed-credit deltas are nonnegative in every dataset and every held
  fold;
- nested macro PDS and signed-credit deltas are each at least `+0.003`;
- the final all-development arm separately clears the same per-dataset and
  macro thresholds;
- both target-block bootstrap lower 95% bounds are strictly positive;
- both observed endpoints exceed both null p95 values;
- maximum positive-gain concentration by biological family and by target
  block is at most `0.5`;
- all sign, norm, finite-value, support and exact-fallback invariants pass.

Any failure freezes a NO-GO and permanently denies this experiment access to
the 410-target test.

## Untouched test

`test` is allowed only after an authenticated development candidate freeze
states that every development gate passed. It uses the final development arm
without refitting or retuning and applies the same per-dataset, macro,
bootstrap, null, concentration and invariant gates to all 410 targets.

A passing test authorizes separate scientific review only. It does not
authorize counts, `cell-eval2`, a challenge build, or submission.

## Commands

Metadata-only audit:

```bash
uv run python experiments/138-dld1-replicate-sign/run.py metadata
```

Create the metadata-only immutable plan:

```bash
uv run python experiments/138-dld1-replicate-sign/run.py plan
```

Source-safe tests:

```bash
uv run python -m unittest \
  experiments/138-dld1-replicate-sign/test_run.py
uv run ruff check experiments/138-dld1-replicate-sign
```

The following commands are intentionally **not run during implementation**:

```bash
uv run python experiments/138-dld1-replicate-sign/run.py freeze

uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/138-dld1-replicate-sign/run.py develop \
  --confirm I_CONFIRM_EXPERIMENT_138_DEVELOPMENT_ONLY

uv run python experiments/138-dld1-replicate-sign/run.py authorize-test \
  --confirm I_REVIEWED_EXPERIMENT_138_DEVELOPMENT_AND_AUTHORIZE_TEST

uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/138-dld1-replicate-sign/run.py test \
  --confirm I_CONFIRM_EXPERIMENT_138_UNTOUCHED_TEST
```

## Hard boundary

- All writes are confined below `experiments/138-dld1-replicate-sign/`.
- No HCASMC, TeloHAEC, Pan panel, Experiment 137 Jurkat-79, challenge outcome,
  `cell-eval2`, count generation, submission build, submission, commit, or
  push is authorized.
- Pan H1 remains only the already-authorized Experiment 132 training context;
  no Pan panel outcome is read.
- Development and test are irreversible, separately locked stages.
- A passing development result does not itself unlock the test. A separate
  post-development authorization must bind the reviewed result and candidate
  freeze before the test command can run.
