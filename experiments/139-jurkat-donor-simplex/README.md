# Experiment 139 — high-MOI Jurkat landmark donor simplex

Status: **implemented and source-safe tested; no Experiment 139 development
or test outcome array has been opened and no authorization exists**.

## Hypothesis

Experiment 107's high-MOI Jurkat estimate is a narrow but target-specific
landmark readout. It may be more useful for choosing among already-allowed
source contexts than as an expression donor.

For each target and held biological group, Experiment 139 compares the
high-MOI Jurkat response with each of Experiment 132's five allowed
source-context vectors on 138 shared exportable landmark genes. Those cosine
similarities determine convex source weights. The Jurkat effect is never
added to, blended with, or substituted for the prediction.

## Frozen inherited surfaces

- Exact Experiment 132 v2 development targets: 613.
- Exact Experiment 132 v2 untouched test targets: 410.
- Exact Experiment 132 v2 response axis: 4,992 genes.
- Outcomes: K562-GWPS, K562-essential, RPE1 and Nadig Jurkat.
- Both K562 screens remain one exclusion group.
- Comparator: the exact equal-weight mean of the five source contexts allowed
  by Experiment 132 for each held biological group.
- High-MOI support: 609/613 development targets and 409/410 test targets.
- Landmark intersection: 138 exportable Experiment 107 response genes on the
  4,992-gene axis.

Unsupported or degenerate targets are bit-exact comparator fallback.

## Candidate family

For target `p`, allowed source `c`, and landmark vectors `j_p` and `x_pc`:

```text
s_pc = cosine(j_p, x_pc)
a_pc = softmax(s_pc / temperature)
w_pc = (1-rho)/5 + rho*a_pc
prediction_p = sum_c w_pc * x_pc
```

Candidates:

```text
identity
rho = 0.25, 0.50
temperature = 0.10, 0.25
```

Weights must be finite, nonnegative and sum to one. Every candidate coordinate
must lie within the coordinatewise minimum and maximum of the five allowed
source vectors. High-MOI values affect only the weights.

## Selection and gates

The 613 targets are divided into five contiguous folds preserving Experiment
132 order. Each held fold uses an arm selected on the other four folds, with
those inner targets also serving as the inner PDS gallery so held truth cannot
influence selection as a decoy. The selected arm is evaluated on held queries
against the complete 613-target gallery. Selection requires nonnegative PDS
and signed credit in every dataset; ties prefer identity, then smaller `rho`,
then larger temperature. One final arm is selected on all development targets
and frozen before test.

There are 49 frozen target-label nulls and 49 landmark-coordinate nulls.
Development nulls repeat complete nested selection. Test nulls use the fixed
development arm without retuning.

Target-block bootstrap uses 20,000 draws and Experiment 132's reciprocal-5NN
DepMap components at cosine 0.40.

Development and test each require:

- nonnegative PDS and signed-credit deltas in every dataset;
- on development, nonnegative deltas in every dataset and held fold;
- macro PDS and signed-credit deltas at least `+0.003`;
- positive bootstrap lower 95% bounds;
- both endpoints above both null p95 values;
- biological-family and target-block positive-gain concentration at most 0.5;
- all simplex, convex-hull and fallback invariants.

Any development failure permanently denies test access.

## Commands

Source-safe:

```bash
uv run python experiments/139-jurkat-donor-simplex/run.py metadata
uv run python experiments/139-jurkat-donor-simplex/run.py plan
uv run python -m unittest experiments/139-jurkat-donor-simplex/test_run.py
uv run ruff check experiments/139-jurkat-donor-simplex
```

Not authorized or run during implementation:

```bash
uv run python experiments/139-jurkat-donor-simplex/run.py freeze

uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/139-jurkat-donor-simplex/run.py develop \
  --confirm I_CONFIRM_EXPERIMENT_139_DEVELOPMENT_ONLY

uv run python experiments/139-jurkat-donor-simplex/run.py authorize-test \
  --confirm I_REVIEWED_EXPERIMENT_139_DEVELOPMENT_AND_AUTHORIZE_TEST

uv run python scripts/run_with_memory_guard.py --max-rss-gib 24 -- \
  uv run python experiments/139-jurkat-donor-simplex/run.py test \
  --confirm I_CONFIRM_EXPERIMENT_139_UNTOUCHED_TEST
```

All writes stay below this directory. HCASMC, TeloHAEC, Pan panel outcomes,
Experiment 137's Jurkat-79 rows, challenge outcomes, `cell-eval2`, counts,
submission construction, submission, commit and push are forbidden.

A passing development result does not itself unlock test access. A separate
post-development authorization must bind the reviewed result and candidate
freeze before the test command can run.
