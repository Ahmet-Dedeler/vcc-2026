# 132 — DepMap response-decomposition Ridge

This experiment tests whether the upstream PCA-50/Ridge-10 model can improve
cross-context perturbation transfer over the exact same-target source mean.

The evaluation is source-only:

- outer outcome datasets: K562-GWPS, K562-essential, RPE1, and Nadig Jurkat;
- K562 screens are one biological exclusion group, so neither is used when
  either K562 screen is scored;
- Pan H1, Orion HCT116, and Orion HEK293T are training-only contexts and are
  never claimed as fresh confirmation outcomes;
- all 300 VCC response coordinates and all gallery-target coordinates are
  removed before fitting or scoring;
- effects are raw, uncentered source deltas;
- metrics are complete-gallery raw PDS and directional signed credit.

The 1,023 common non-panel targets are split by a frozen hash into 613
selection targets and 410 untouched source-test targets. The only candidates
are:

1. `identity_source_mean`: equal mean of the allowed same-target source
   responses;
2. `depmap_pca50_ridge10`: the upstream model family.

Selection and test both require nonnegative deltas on every scored dataset,
macro PDS gain of at least `+0.003`, positive macro signed-credit gain, and
target-label-null control. The untouched test additionally requires positive
target-block-bootstrap lower bounds for both endpoints.

Commands:

```bash
uv run python -m unittest \
  experiments/132-response-decomposition-ridge/test_run.py

uv run python experiments/132-response-decomposition-ridge/run.py plan
uv run python experiments/132-response-decomposition-ridge/run.py freeze
uv run python experiments/132-response-decomposition-ridge/run.py select
uv run python experiments/132-response-decomposition-ridge/run.py evaluate
```

Run `select` and `evaluate` through `scripts/run_with_memory_guard.py`.

No held challenge outcome, `cell-eval2`, counts, challenge artifact,
submission, commit, or push is authorized. A source-only pass authorizes only
separate review.
