# Experiment 130 pre-value result

Status: **LOCKED; metadata-only preparation complete**

The TeloHAEC confirmation has not been run. No TeloHAEC `X` or layer
expression value, compact TeloHAEC effect/control value, or Experiment 125
candidate array was opened.

The frozen source has exactly 65 VCC target rows. Every row has at least two
eligible guides and at least 20 eligible cells:

- eligible guides per target: minimum 2, median 8, maximum 18;
- eligible cells per target: minimum 20, median 67, maximum 147;
- eligible control guides: 549;
- all 300 canonical panel response symbols are present on the source axis.

The value gate remains closed because
`CANDIDATE_BINDING_AMENDMENT.json` does not exist and Experiment 125 has no
immutable result directory yet. The later amendment must bind an Experiment
125 `GO`, all passing gates, the frozen modal candidate selected in at least
8/15 folds, `production_log1p` transport, and exact residual/confidence/support
artifacts before any TeloHAEC expression value may be read.

Current verification:

- metadata preflight: passed;
- synthetic fail-closed tests: 6/6 passed;
- Ruff: passed;
- Python compilation: passed;
- biological value access: zero;
- `cell-eval2`, counts, challenge build/submission, commit, and push: not run.
