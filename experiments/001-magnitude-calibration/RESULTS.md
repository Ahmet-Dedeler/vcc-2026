# 001 — Magnitude calibration is free points

**Date:** 2026-08-21 · **Status:** analytic + numeric, confirmed. Needs end-to-end confirmation
against `cell-eval2` on 2025 H1 once that data lands.

## Claim

`pds_cosine` is exactly invariant to rescaling the predicted effect vector; `expr_mse_unbiased_capped_norm`
is quadratic in it. So predicted-effect magnitude `α` is a free parameter that moves `mse` without
touching `pds`.

## Result

1. **PDS invariance is exact.** `PDS_p` is bit-identical for α ∈ {0.01, 0.1, 0.25, 1, 10, 100}.
   It follows directly from the spec: `d_{p,q} = 1 − cos(δ̂_p, δ_q)`, and cosine ignores `‖δ̂_p‖`.
2. **`min_α mse(α) = 1 − cos²θ`**, attained at `α* = ⟨d,δ⟩/⟨d,d⟩`. Numerics match analysis to 4dp.
3. **`mse(α=1)` values reproduce the leaderboard.** cos 0.5 → 3.00 (Vivai 3.48, Sidechain 3.61,
   joshua 3.96, Jurassic 3.98); cos 0.3 → 10.1 (Gia Khang 6.54, p19 6.97); cos 0.7 → 1.04
   (wombat 1.56). The whole field is sitting on the α=1 curve.

## Why this is airtight

`min_α mse(α) ≤ mse(α=0) ≈ 1.003` — the control-copy value. So **any submission with raw mse > 1
is provably not α-optimised**, with no assumption about the model that produced it. All ten
top-10 teams have raw mse ≥ 1.562.

## Size of the prize

At the leader's implied cos θ ≈ 0.5: raw mse 3.48 → 0.75, scaled 0.000 → 0.251,
overall **0.134 → 0.176 (+31%)**, with `pds` unchanged.

The gain scales with directional accuracy — at cos 0.7 scaled mse is 0.503, at cos 0.9 it is 0.838.
Shrinking a *directionless* prediction gains nothing (control-copy is still 0.000). **α converts
directional accuracy you already have into mse points you are currently throwing away.**

## Caveat / next step

Shrinking α reduces per-gene effect sizes, so fewer genes clear Wilcoxon+BH. That lowers `n_pred`
and triggers `fid`'s yield penalty `min(1, n_pred/n_real)`, and lowers `jac`. The counter-knob is
**NB dispersion** in the emission step: less cell-to-cell noise restores significance at the same α.

→ **Next: joint (α, dispersion) grid against the real scorer on 2025 H1.** That is experiment 002
and it must run before any modelling work.

## Reproduce

```bash
uv run python experiments/001-magnitude-calibration/calibration_math.py
```
