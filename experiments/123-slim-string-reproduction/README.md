# 123 — exact SLIM-STRING reproduction

## Question

Can the exact upstream SLIM model learn an unseen-target perturbation response
from the pinned 64-dimensional STRING embeddings, and can that learned
direction improve the existing direct donor-transfer incumbent across unseen
cell contexts?

This is **not** Experiment 087. Experiment 087 fit a low-rank STRING regression
to the residual left by Experiment 064. This experiment reproduces the
upstream SLIM family itself:

```text
Y_hat = G @ W @ P.T + b
```

- `Y` is the training perturbation pseudobulk matrix.
- `b` is the mean training pseudobulk.
- `G` is the PCA score basis of `Y - b`.
- `P` is the pinned STRING target-embedding matrix.
- `W` is the upstream closed-form two-sided ridge solution.

The source-only experiment separately reports:

1. the exact upstream absolute-profile prediction converted back to an effect;
2. the centered SLIM response used as a donor-transfer direction; and
3. a production adaptation that blends that direction with the direct
   donor-transfer incumbent.

Only item 1 is called an exact upstream reproduction. Items 2 and 3 are
explicitly labeled adaptations.

## Leakage boundary

The three source-validation donors are:

- Replogle K562-GWPS;
- Replogle RPE1; and
- Nadig Jurkat.

H1 and both Orion contexts are excluded from every source-only fit, selection,
null, and outcome. They remain external held surfaces.

Before any source effect value is loaded, the experiment freezes:

- a deterministic 600-target STRING-covered non-panel development set;
- 360/120/120 fit/tune/test target partitions;
- one common response axis;
- all target-pairing null assignments;
- all Gaussian-null seeds;
- target-cluster bootstrap assignments;
- the model and candidate grid;
- the exact metric definitions and promotion gates; and
- hashes of every implementation and input file.

Every one of the 300 VCC target coordinates and every coordinate belonging to
any of the 600 development targets is absent from the response axis. Held-donor
truth is raw and uncentered. Every PDS comparison uses the complete frozen
120-target gallery.

## Commands

The commands are intentionally staged:

```bash
uv run python experiments/123-slim-string-reproduction/run.py plan
uv run python -m unittest experiments/123-slim-string-reproduction/test_run.py
uv run ruff check experiments/123-slim-string-reproduction
uv run python experiments/123-slim-string-reproduction/run.py audit
uv run python experiments/123-slim-string-reproduction/run.py select
uv run python experiments/123-slim-string-reproduction/run.py evaluate
uv run python experiments/123-slim-string-reproduction/run.py verify
```

`select` may read fit and tune rows only. It writes a sealed candidate before
`evaluate` is allowed to read any test row. `evaluate` refuses to overwrite an
existing result.

No command invokes `cell-eval2`, creates counts, builds a challenge artifact,
submits, commits, or pushes.

## Promotion

Promotion requires every gate in `preregistered_protocol.json`. A pass
authorizes only preparation of a separately reviewed H1/Orion evaluator. It
does not authorize held evaluation, count generation, a challenge build, or a
submission.

