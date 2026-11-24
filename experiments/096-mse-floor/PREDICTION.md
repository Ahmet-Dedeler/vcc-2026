# 096 — pre-registered prediction for the amplitude sweep

Written 2026-08-28 13:40 UTC, before any config in `alpha_grid.json` finished.

## The claim being tested

Our board `expr_mse_unbiased_capped_norm` is **1.3627**. The baseline (context mean response)
is ~0.993, and **predicting the control profile unchanged scores ~1.0** — because the metric
normalises by the measured profile's distance from the control, so a prediction of "no change"
lands at a ratio of about one by construction.

**We are at 1.3627, which is worse than predicting no change at all.** That is not a weak
model; a weak model converges to 1.0. Landing *above* 1.0 requires actively moving away from
the truth, which for a prediction of the form `b_hat = expm1(b_ctrl + alpha*e)` means `alpha`
is past the point where the effect helps.

## The algebra

In log space write the prediction as `b_ctrl + alpha*e` and the truth as `b_ctrl + t`. Then

```
||alpha*e - t||^2  =  alpha^2 ||e||^2  -  2*alpha*<e,t>  +  ||t||^2
```

an upward parabola in `alpha`, minimised at

```
alpha* = <e,t> / ||e||^2 ,      value  ||t||^2 (1 - rho^2),   rho = corr(e, t)
```

Normalising by the metric's denominator (roughly `||t||^2`), this predicts

- `mse_raw(0) ~ 1.0` — predict the control
- `mse_raw(alpha*) ~ 1 - rho^2` — the floor a global amplitude can reach
- `mse_raw(alpha) > 1` whenever `alpha > 2*alpha*`

**So `mse_raw = 1.3627 > 1` implies our production amplitude is more than TWICE optimal.**

## Predictions, with numbers

1. **`mse_raw` is U-shaped in `alpha`**, with an interior minimum. **Falsified if** it is
   monotone over `[0.30, 1.40]`, or flat (range < 0.02), or minimised at the top of the grid.
2. **The minimum sits at `alpha*` in `[0.25, 0.55]`.** Our 08-25 board point `Alpha 0.40`
   scored `0.9640`, so the true optimum is near there. **Falsified if** the argmin is above
   0.8.
3. **The minimum value is `~0.93-0.97`**, i.e. `rho^2 ~ 0.03-0.07`, an effect-truth
   correlation of only **`rho ~ 0.17-0.26`**. This is the uncomfortable part of the
   prediction: it says a global amplitude can buy at most about `+0.025-0.030` scaled, and
   that the leader's `0.8270` (implying `rho ~ 0.42`) is **out of reach by recalibration
   alone** — it needs a better effect direction, not a better scale.
   **Falsified if** the harness minimum is below 0.90.
4. **`reach` and `fid` rise monotonically with `alpha`** over the grid, because a larger
   effect makes more genes clear significance. So there is a genuine crossing and the
   sweep measures its position. **Falsified if** `reach` is flat or falls with `alpha`.
5. **`pds` is nearly flat in `alpha`** — it is a cosine metric and cosine is scale-invariant;
   the only `alpha` dependence is through the emission and the CPM gate.
   **Falsified if** `pds` moves by more than 0.02 scaled across the grid.

## What would change our plans

If (1)-(3) hold, a global amplitude is worth `+0.025` on `mse` and costs `reach`; the
interesting question becomes whether a **non-global** shrinkage (per-gene or per-perturbation,
James-Stein style) beats the global bound, because only that can raise `rho` rather than
merely optimise against it.

If (3) is falsified low — the harness minimum well under 0.90 — then our effect direction is
better than the board implies and the gap is calibration after all, which would be very good
news and would make this the single highest-value lever on the board.

## Known limitation, stated up front

This runs the `031-somatic` model (donor transfer + whitening), **not** the full production
path (PerturbMatch + CD4 + pruned operator). It measures the shape of the `alpha` response and
the position of the optimum, not the production numbers. Harness `n_real` is 92, so `fid`,
`jac` and `reach`'s normalisation are unreadable in magnitude here — see
`docs/RESEARCH-METHOD.md` rule 1. `mse` is the member this harness reads most honestly,
because it is a profile-distance ratio with no `n_real` in it.
