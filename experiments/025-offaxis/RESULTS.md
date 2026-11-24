# 025 — Whitening silently down-weights the genes it cannot reach

**Date:** 2026-08-21 · LODO with a **union** feature axis, so the situation can actually arise.

## The bug

The operator is fitted on genes **every** donor measures — 7,312 of 18,533 on the production
axis. The model predicts on genes **any** donor measures. Whitening multiplies the on-axis
block's norm by ~4.6x, and v5 then rescales the *whole* effect matrix back to its prior median
row norm. Net: **off-axis genes are down-weighted 4.6x relative to on-axis ones**, silently.

At the CPM > 20 gate the scorer tests 7,989 genes in context A: 6,102 on-axis, **1,852
off-axis** (Orion-only), 35 with no donor at all. So it is 23% of the tested universe, and no
previous fold could see it — every LODO run so far scored on the donor *intersection*.

## Result

Feature axis = genes the held-out line shares with **any** predictor; operator axis = the
subset **all** predictors carry. Holding out either Orion line leaves 80% of the scored genes
off-axis (the two Orion lines share a gene axis, so holding out K562 leaves 0% and that fold
is uninformative here).

| | k562 (0% off) | hct116 (80% off) | hek293t (80% off) | macro |
|---|---|---|---|---|
| no whitening | 0.7510 | 0.7123 | 0.7159 | 0.7264 |
| global rescale (what v5 does) | 0.8045 | 0.7658 | 0.7440 | 0.7714 |
| **preserve on-axis norm** | 0.8045 | **0.7744** | **0.7536** | **0.7775** |

**+0.0061 macro**, and the two folds where off-axis genes exist gain +0.0086 and +0.0096 while
the fold where they do not is bit-identical — which is the control this experiment needed.

## Fix

`--whiten-renorm preserve`: rescale only the on-axis block back to its own prior norm and
leave off-axis genes untouched. Production impact should be smaller than the LODO number here
(23% of tested genes off-axis rather than 80%), but the direction is unambiguous and the fix
is free.

**General lesson.** A transform applied to part of a feature space needs its magnitude
restored **on that part**, not globally. Global renormalisation turns a partial transform into
an unintended reweighting of everything it did not touch.
