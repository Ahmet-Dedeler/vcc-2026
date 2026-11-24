# 002 — Calibration sweep on the real scorer

**Date:** 2026-08-21 · **Harness:** `h1-fast` (40 perts × 400 cells + 4,000 controls, 2025 H1)
**Scale:** `cell-eval2` 0.16.0, `vcc2026` preset, bundle `h1-fast`. All values **scaled**
(`from_replicate`): 0 = context-mean baseline, 1 = real replicate.

## Round 1 — α sweep at cos θ = 0.5, dispersion = 0.5

| α | pds | mse | fid | reach | jac | nmae | **avg** |
|---|---|---|---|---|---|---|---|
| 0.10 | 1.1113 | 0.0091 | 0.5375 | 0.0230 | 0.0954 | 0.0279 | 0.3007 |
| 0.25 | 1.1113 | **0.1531** | 0.6101 | 0.2275 | 0.1082 | 0.2769 | 0.4145 |
| 0.50 | 1.1113 | 0.0000 | 0.6813 | 0.4569 | 0.1249 | **0.5626** | **0.4895** |
| 1.00 | 1.1113 | 0.0000 | **0.7374** | **0.7582** | **0.1275** | 0.1451 | 0.4799 |

## Two predictions confirmed

**1. `pds_cosine` is exactly scale-invariant.** `1.11127002` at every α, identical to eight
significant figures across a 10× range. Experiment 001 predicted this analytically; the real
scorer reproduces it. The magnitude of the predicted effect is genuinely free for `pds`.

**2. `α*` lands exactly where the theory says.** The synthetic direction is
`d = δ + √3·noise` with `|noise| = |δ|`, so `|d|² = 4|δ|²` and `⟨d,δ⟩ = |δ|²`, giving
`α* = ⟨d,δ⟩/⟨d,d⟩ = 0.25`. Observed `argmax(mse) = 0.25`. Predicted scaled value 0.217 vs
observed 0.153 — the gap is emission noise and the jackknife corrections, which the pure
pseudobulk algebra ignores.

## Correction: the gain is much smaller than I claimed

I estimated tuning α alone would move the leader from 0.134 to 0.176 (+31%). **That estimate was
wrong, because it assumed the other five metrics were indifferent to α. They are not.**

`fid`, `reach` and `jac` all rise monotonically with α — larger effects push more genes past the
Wilcoxon/BH bar, which raises `n_pred`, which lifts `fid`'s yield term and `reach`'s depth. So
shrinking α to serve `mse` costs the DE metrics.

Net on this harness: best single α is **0.5** (avg 0.4895) against α=1.0's 0.4799 — **+0.010, not
+0.04**. `mse` is one sixth of the score and the DE metrics pull the other way.

The lever is real; it is roughly a quarter of the size I first claimed.

## What that implies, and the next test

The DE metrics improve with α only because bigger effects clear the significance bar. **Emission
dispersion reaches the same bar without inflating the pseudobulk effect** — less cell-to-cell
noise means the same effect size is more significant. If that decoupling works, we can hold α at
the `mse` optimum and recover the DE metrics through dispersion.

→ Round 2 (running): α fixed at 0.25, dispersion ∈ {0.1, 0.25, 0.5, 1.0}.

## Caveat that limits round 1

**`pds` raw is 1.0000 — a perfect ranking — at cos θ = 0.5.** On a 40-perturbation panel,
discrimination is trivially easy. The leaderboard's top teams sit at pds raw 0.659–0.784 on a
**300**-perturbation panel, which implies a far lower cos θ than 0.5.

So round 1's cos θ is unrealistically generous and its `pds` column is saturated. Realistic
sweeps must (a) use cos θ in the 0.1–0.3 range and (b) run on a panel large enough that `pds`
is not saturated. `h1-test` (85 perts) is better; the honest version needs ~300.

**Do not read avg = 0.49 as a predicted leaderboard score.** It is the score of an oracle-derived
direction on an easy panel.

## Reproduce

```bash
uv run python experiments/002-calibration-sweep/sweep.py \
  --real data/harness/h1-fast.h5ad --bundle experiments/002-calibration-sweep/fast/bundle \
  --outdir experiments/002-calibration-sweep/fast \
  --cos 0.5 --alpha 0.1 0.25 0.5 1.0 --dispersion 0.5 --jobs 4
```
