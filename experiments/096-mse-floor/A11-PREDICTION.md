# a11 — pre-registered prediction

Written 2026-08-28 14:55 UTC, while the build was running and before any scoring.

**Arm:** a10 exactly (`--tail-spread 0.15 --overdispersion 16.0`, pruned whitener), with the
single change **`VCC_ALPHA` 1.0796 → 0.838** — a 22.4% cut to the log-space effect amplitude.
One variable.

## Why this and not the larger arm

Experiment 096 measured a bigger intervention (raise the jackknife `claim` to its cap via a
split spread policy, *then* trim α by 6%) worth `+0.0192` gross. It is not what is being
shipped, deliberately:

1. **It moves two things at once.** If it disappointed we could not attribute the miss.
2. **Its call-set cost is genuinely two-sided and unmeasured.** The tail policy churns
   8.7-14% of our calls; depending on whether those new calls are better or worse than chance
   the `fid`+`jac` swing runs from `+0.0028` to `−0.0281`.
3. **This arm tests the load-bearing claim on its own.** The entire 096 frame rests on one
   decomposition. It has passed one out-of-sample test already — it predicted a10's raw `mse`
   at `1.3567 ± 0.01` on a build it never saw, and a10 returned **1.3524**. A second,
   sharper test at a very different amplitude is worth more than the marginal score.

## The mechanism being tested

`expr_mse_unbiased_capped_norm` subtracts a jackknife estimate of sampling noise,
`rho * min(C_hat_p, C_p) + C_p`, from the squared error. Because the quasi emitter repairs
each gene's panel total exactly, the *predictions'* common mode is only 0.84% of their energy
— so as amplitude falls the metric refunds nearly all of it, the curve keeps descending, and
**the optimum is not the least-squares `X/A` (0.066) but the `rho` cliff** at
`s* = sqrt(claim/B)`. For a10's measured `claim = 0.1688` that cliff is `VCC_ALPHA = 0.838`.

**This is the counter-intuitive part, and it is what makes the prediction falsifiable:** a
naive least-squares reading says shrink to ~0.07 of current amplitude; the cliff model says
shrink to 0.838 and stop. The two differ by more than a factor of ten.

## Predictions

Baseline is a10 as scored: overall `0.1541841`, `mse` raw `1.3524148`, `score_mse` `0`,
`pds` `0.7365610`, `nmae` `0.0943384`, `fid` `−0.0278969`, `reach` `0.1235951`,
`jac` `−0.0014931`.

| member | predicted | reasoning |
| --- | ---: | --- |
| `mse` raw | **0.9021** | the cliff value at `claim = 0.1688` |
| `score_mse` | **+0.0903** | `(0.987862879 − 0.9021) / 0.950162566` |
| `reach` | **−0.0214** | 22.4% amplitude cut on the board's own α ladder |
| `nmae` | **−0.0150** | same ladder |
| `pds` | **0.0000 ± 0.005** | `pds_cosine` is a cosine; amplitude is scale-invariant |
| `fid`, `jac` | small, unsigned | no call-set policy change, but `n_pred` moves with α |
| **overall** | **≈ 0.1631** (`+0.0089`) | |

## Falsifiers

1. **`score_mse` returns 0.** The floor is not cleared, the cliff is in the wrong place, and
   the whole 096 frame is wrong rather than merely mis-sized. This is the one that matters.
2. **`mse` raw lands outside `[0.87, 0.94]`.** The decomposition is directionally right but
   quantitatively unusable, and no further arm should be built on its point estimates.
3. **`pds` moves by more than `0.005`.** Amplitude is not scale-invariant through our gate and
   the CPM gate is doing something we have not modelled.
4. **Overall lands below `0.1541841`.** The amplitude cut costs more in `reach`+`nmae` than
   `mse` pays. Note this can happen *even if* predictions 1-3 all hold — the mse model can be
   right and the trade still bad. If so, the correct next move is the split tail policy, which
   buys the same `claim` at a 6% cut instead of 22%.

## What each outcome implies for the next slot

- **As predicted (~0.163):** the model is board-confirmed twice. Ship the tail policy plus a
  6% α cut next, which the same model prices at `+0.0192` gross with a much smaller α cost.
- **`mse` clears but overall falls:** the tail policy becomes the *only* sensible route,
  because it reaches the same `claim` without the amplitude penalty.
- **`mse` stays 0:** stop building on 096's point estimates, and go re-derive the cliff from
  the two-point board solve with a third point.

---

## Amendment, 15:22 UTC — prediction 5 is already in doubt, before submission

The independent amplitude sweep run for this experiment (`alpha_sweep.log`, the `031-somatic`
model scored by the real `cell-eval2` on `orion-hct116-thinned`) returned three of six points
before a11 was ready, and they contradict prediction 5:

| `alpha` | `pds` raw | `pds` scaled |
| ---: | ---: | ---: |
| 0.30 | 0.7797 | 0.618 |
| 0.45 | 0.7942 | 0.650 |
| 0.60 | 0.8017 | 0.667 |

**`pds` rises monotonically with amplitude** — `+0.049` scaled across a 2x change. I predicted
`0.0000 ± 0.005` on the grounds that `pds_cosine` is a cosine and cosine is scale-invariant.
The *statistic* is scale-invariant; **the pipeline that produces it is not.** Counts are
integers, so at lower amplitude more of the effect quantises away and the emitted profile sits
closer to the control, which is genuinely less discriminable. Scale-invariance of the formula
says nothing about a stochastic integer emitter feeding it.

**This is not a reason to withhold the submission, and I am not revising the prediction after
the fact — I am recording that one of its five legs is already weakened.** Specifically:

- The sweep is the `031-somatic` model, not the production pipeline (its `pds` runs 0.62-0.67
  against our board 0.7366), so the slope may not transfer.
- The curve should flatten toward `alpha = 1`; the measured slope is in the steep low-alpha
  region. Extrapolating the flattened slope to a 22.4% cut gives roughly `−0.012` scaled,
  i.e. **`−0.002` overall**.
- If that holds, a11's net is nearer **`+0.007` than `+0.0089`**, and prediction 3's
  `±0.005` band is too tight by about 2x.

**Falsifier 3 stands as written and is now expected to fire.** Recording it beforehand is the
point: if `pds` drops by ~0.012 the mechanism is understood and the cliff model is untouched;
if it drops by much more, amplitude is doing something to the emitter we have not modelled and
intervention 1's 6% cut needs re-pricing before it ships.
