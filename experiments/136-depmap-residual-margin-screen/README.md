# 136 — DepMap residual margin screen

This is a development-only screen on Experiment 132's already-consumed
613-target selection split. It asks whether the rejected
`depmap_pca50_ridge10` replacement becomes useful when expressed as a small
residual around the exact same-target source mean.

It evaluates:

1. a fixed global residual-dose grid;
2. a deployment-available per-target rule that chooses the largest dose whose
   PDS rank does not worsen in any observed source context;
3. a stricter rule that preserves every observed source-context nearest-decoy
   cosine margin;
4. an outcome-informed per-target oracle used only as an upper bound on the
   consumed development split.

The Experiment 132 untouched 410-target split remains sealed. This screen can
reject the proposal family, but cannot authorize production or held-test
access.

Run:

```bash
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 8 \
  --min-free-percent 8 \
  --poll-seconds 0.5 \
  -- \
  uv run python experiments/136-depmap-residual-margin-screen/run.py
```

No `cell-eval2`, count generation, challenge build, submission, commit, or
push is authorized.
