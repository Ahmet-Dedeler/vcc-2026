# 127 — Panel-only metric ceiling

Read-only analysis of how much the four differential-expression members can
move when only the 299 off-diagonal VCC panel-response coordinates are
mutable.

Every biological surface here is reused development evidence. Nothing in this
experiment is an independent confirmation.

The experiment reports oracle/reference-composition ceilings only. An oracle
may set mutable panel values, calls, signs, and ranking keys using the
already-consumed reference while every non-panel prediction stays fixed.

Only the repaired HCT116 surface has an already-existing exact prediction DE
table. Exact fidelity, Jaccard, and reach counterfactuals are therefore
reported only there. H1 and HEK293T receive exact reference-composition and
oracle NMAE analyses.

Experiment 124's transport is independently rejected and is not used or cited
as a valid mechanism here. Any buildable panel transport is deferred to
Experiment 128.

Pan has already-consumed outcomes, but its compact effect surface has no
scorer-equivalent genome-wide significance table. It is excluded from the
metric ceiling rather than silently promoted into one.

Run:

```bash
uv run python experiments/127-panel-metric-ceiling/freeze.py
uv run python -m unittest experiments/127-panel-metric-ceiling/test_run.py
uv run ruff check experiments/127-panel-metric-ceiling
uv run python experiments/127-panel-metric-ceiling/run.py --analysis
```
