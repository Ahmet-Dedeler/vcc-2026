# 035 — The anatomy of our `nmae`, and why it has no per-gene lever

**Date:** 2026-08-22 · somatic HCT116 reference, 132 real VCC-panel targets, shipped v8 effect
matrix, top-300 genes by |true effect| standing in for the reference's DE set.

`nmae` is our largest genuinely accessible gap — ours +0.065 scaled against Vivai's +0.131,
with `mse` above us closed at our cosine. It had only ever been attacked with global knobs
(alpha, the CPM gate, gamma). This asks which *genes* carry the error.

## The decomposition

Every gene contributes to `nmae` in exactly one of three ways:

| category | genes | share of error | contribution to `nmae` |
|---|---|---|---|
| **silent** (`pred = 0`) — contributes exactly `\|true\|`, a ratio of 1 | 1.7% | 2.2% | 0.0215 |
| **wrong sign** — contributes `\|pred\| + \|true\|`, a ratio **above** 1 | 45.9% | **52.7%** | **0.5120** |
| right sign — contributes `\|\|pred\| − \|true\|\|` | 52.4% | 45.1% | 0.4389 |
| | | | **total 0.9724** |

**Coverage is not the problem** — silent genes are 1.7% of the set and 2.2% of the error, so
the CPM gate and donor coverage are already close to irrelevant here.

**Confident wrong-sign calls are the problem**: 53% of the error, and they are the only
category that can score *worse* than saying nothing.

**Counterfactual: silencing every wrong-sign gene would give `nmae` 0.9090** — past Vivai's
0.9213. So the prize is real and large (+0.105 scaled, +0.018 overall). The question is whether
any signal we have identifies those genes.

## The one signal we have does not separate them enough

Cross-donor sign agreement is the natural candidate, and it is the signal that won experiment
013's ranking test.

| group | genes | sign accuracy | share of `nmae` |
|---|---|---|---|
| donors agree | 12,440 | **55.8%** | 0.2914 |
| donors disagree | 26,490 | **52.1%** | 0.6596 |

**Separation: 3.7 percentage points**, and *both groups are still majority right-signed*.

That is decisive, and the arithmetic is simple. Shrinking a group with sign accuracy `a` by a
small amount changes the error in proportion to `(1−a)·(−1) + a·(+1) = 2a − 1` — it helps only
when `a < 0.5`. At `a = 0.521` the disagreeing group is already **net worth predicting**
(2a − 1 = +0.042), so shrinking it costs more than it saves. At `a = 0.558` the agreeing group
is worse still to shrink.

**There is no subgroup we can profitably silence.** This also explains, quantitatively, why
experiment 029's sign-agreement mask was the worst of its three variants: it silenced a group
that was 52% correct.

## What this closes, and what it points at

`nmae` cannot be improved by per-gene calibration, shrinkage, gating or confidence weighting at
our sign accuracy. Its 0.063 of available headroom is **entirely a function of directional
accuracy on the reference's DE genes**, which is the same binding constraint as `pds`, `mse`
and `jac`.

Combined with the earlier results this closes the last of the "calibrate our way up" family:

| metric | why it is closed |
|---|---|
| `mse` | needs cos ≈ 0.20; we have ~0.13, and whitening trades cosine away deliberately (exp 033) |
| `jac` | 1.31x chance at identifying which genes respond (exp 014) |
| `nmae` | no subgroup is below 50% sign accuracy, so nothing is worth silencing (here) |
| per-perturbation alpha | oracle worth 0.25% squared-error reduction (exp 012) |

**Everything routes back to the direction.** The only two things that have ever moved it are
more donor lines and better geometry — and the geometry family is now exhausted (exp 030f).
