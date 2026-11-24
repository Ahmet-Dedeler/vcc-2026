# 026 — Re-expressing donor effects at the target's expression level does not help

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations.

## A real mis-specification

Every effect here is `delta = log1p(TS*p_pert) − log1p(TS*p_ctrl)`, TS = 5e4, matching the
scorer. `log1p` is not a log: below `x = TS*p ≈ 1` (**20 CPM**) it is approximately linear, so
`delta` is an absolute difference in proportion; above it, a ratio. The local compression
factor `x/(1+x)` is 0.20 at 5 CPM, 0.50 at 20 CPM, 0.83 at 100 CPM — and the crossover sits in
the middle of the tested gene universe.

So `delta` depends on **that line's own expression level of that gene**. A gene at 500 CPM in
K562 and 20 CPM in the target has its donor delta on the near-ratio branch and its true target
delta half-compressed. Transferring `delta` directly imports the donor's compression into the
target's space; the error is systematic, gene-specific, and largest where most genes live.

The correction is a change of variables, not a heuristic: convert to fold change at the
donor's level, then back at the target's.

```
lfc     = log(expm1(b_donor + delta)) − log(expm1(b_donor))
delta_T = log1p(expm1(b_target) * exp(lfc)) − b_target
```

## Result

| | gamma = 0 | gamma = 0.65 |
|---|---|---|
| **plain (transfer delta directly)** | **0.7210** | **0.7978** |
| half (blend) | 0.7209 | 0.7961 |
| lfc re-expression | 0.7193 | 0.7937 |

It loses, by 0.0017 unwhitened and 0.0041 whitened. Consistent across all three folds.

## Why, and it is a general point

The correction is right about the bias and wrong about what matters. `expm1(b_ctrl)` is tiny
for low-expression genes, and dividing by it to recover a fold change **amplifies noise
exactly where the compression was largest**. The compression `log1p` applies is not only a
readout artefact — it is a variance-stabilising transform, and it is suppressing noise on the
same genes whose bias the correction removes.

So the bias it fixes is cheaper than the variance it releases. **The scorer's own
transformation is doing useful work, and undoing it to "correct" a systematic error is a net
loss** — worth remembering the next time a first-principles correction to the measurement
model looks obviously right.

Third independent failure for the target's basal profile (after donor weighting, exp 010, and
per-gene modulation, exp 017), now including one where it enters as the *parameter of a
change of variables* rather than as a heuristic weight.
