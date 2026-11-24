# v9 board result (2026-08-23)

Entry `80wiD09NGbf67RlP1sBE`, model name `The mitochondria declined to comment`,
published 09:28:26 UTC. **Overall `0.1443973094`, rank 2.**

| metric | v8 | v9 | delta |
| --- | ---: | ---: | ---: |
| pds   | `0.70137` | `0.72718` | **`+0.02581`** |
| mse   | `0.00000` | `0.00000` | `0` |
| nmae  | `0.06518` | `0.07313` | `+0.00795` |
| fid   | `-0.03757` | `-0.03812` | `-0.00055` |
| reach | `0.09460` | `0.10178` | `+0.00717` |
| jac   | `0.00269` | `0.00242` | `-0.00027` |
| **overall** | **`0.13771`** | **`0.14440`** | **`+0.00668`** |

Five of six moved the way h1-test said they would, which is the prediction
b052969 made. The magnitude did not: h1-test projected `+0.0104` avg6 and the
board delivered `+0.00668`. **The h1-test -> board realization ratio is 0.64**,
and that is now the second calibration point. The first (v8) was accurate to
within `0.00002`; this one is not. The difference between them is that v8's
change touched only fid/reach/jac while v9's moved pds and nmae, so the
optimistic reading is that h1-test is well calibrated on the DE metrics and
overstates pds gains by about a third.

Use 0.64 as the discount on any future h1-test pds projection until a third
point says otherwise.

## Standing after v9

Ibrahim Mansour leads at `0.146665`; we trail by `0.00227`.

## Calibration: h1-test is reliable on pds, not on the DE metrics

An earlier version of this file said the opposite. It was wrong, and the error
came from comparing h1-test's own six-metric average against the board's scaled
average. Those use different scaling functions, so the "0.64 realization ratio"
was an artefact, not a property of the harness.

Compared **raw delta to raw delta**, which is the only apples-to-apples
comparison, h1-test looks like this:

| metric | h1-test predicted | board delivered | realization |
| --- | ---: | ---: | ---: |
| pds   | `+0.01030` | `+0.01152` | **1.12x** |
| nmae  | `-0.00580` | `-0.00368` | 0.64x |
| reach | `+0.02680` | `+0.00628` | **0.23x** |
| jac   | `-0.00060` | `-0.00011` | 0.18x |
| fid   | `+0.00040` | `-0.00016` | -0.40x |

h1-test predicts **pds** almost exactly (it slightly under-predicted). It
overstates **reach** by 4x and **jac** by 5x. fid flipped sign but at a
magnitude too small to read anything into.

### How to project a board score from h1-test, correctly

Predict the *raw* delta on h1-test, apply the per-metric realization above,
then multiply by the board's measured raw-to-scaled amplification (from v8->v9):

| metric | raw->scaled amplification |
| --- | ---: |
| pds | `2.24x` |
| nmae | `-2.16x` (raw down = scaled up) |
| fid | `3.46x` |
| reach | `1.14x` |
| jac | `2.59x` |

Do not read h1-test's `avg6` directly. It is computed against h1-test's own
baseline and does not map onto the board.

## The real structure of this competition: transfer cosine

> ## ⚠ RETRACTION, 2026-08-23 (later the same day)
>
> **The transfer-cosine numbers below are unfounded. Do not use them.**
>
> They were derived by applying `min_α mse = 1 − cos²θ` to the board's
> `expr_mse_unbiased_capped_norm` column. That identity holds only when the
> metric is `||pred − true||² / ||true||²`, i.e. when predicting *nothing*
> scores exactly `1.0`. Fitting `mse(α) = c0 − 2αA + α²B` to exp077's three
> alpha points gives **`c0 = 1.3388`**, not `1.0`, so the board column uses a
> different normalization — most likely against a reference prediction rather
> than against zero. Under the fitted model the cosine would be
> `A/√(B·c0) = 0.508`, not the `0.0846` published here.
>
> Both numbers are untrustworthy. Three points fit a three-parameter quadratic
> **exactly**, so the "perfect fit" cited as validation had zero degrees of
> freedom and proved nothing, and the metric is *capped*, which can bend the
> parabola away from a true quadratic.
>
> **Specifically retracted:** our cos θ ≈ `0.0846`; Ibrahim's `cos θ ≥ 0.3008`;
> "3.6x our directional accuracy"; the `1.25x` improvement target; the
> `cos > 0.1314` break-even for low alpha.
>
> **What still stands, being direct observation rather than derivation:**
> our Orion floor is `0.99286` at `α* = 0.392` against a `0.989` threshold, so
> **we miss by `0.0039`**; exp077's low-alpha arms are negative; `pds` is
> scale-invariant in alpha (`+0.00026` across a 2.7x change), so the real trade
> is **nmae-for-mse, not pds-for-mse**; and Ibrahim's `0.9095` is the only
> sub-1.0 raw mse on the board.
>
> Resolving the normalization needs a fourth alpha point — see experiment 088.



`STRATEGY.md:61` says mse is "closed, permanently, for everyone". The second
half is false and the first half is right, for a reason worth stating exactly.

`expr_mse_unbiased_capped_norm` is `||pred - true||^2 / ||true||^2`, so
**predicting no effect at all scores exactly 1.0**, and the baseline to beat is
`0.989`. Every team except one is above 1.0 -- that is, actively worse than
predicting nothing:

| team | raw mse | implied cos(theta) |
| --- | ---: | --- |
| **Ibrahim Mansour** | **`0.9095`** | **`>= 0.3008`** |
| us (v9) | `1.7661` | no bound |
| Cell-tics | `1.8365` | no bound |
| kriagr | `3.7001` | no bound |
| VIVAI | `5.9465` | no bound |
| pop | `6.7378` | no bound |
| Remigiusz Kinas | `10.7532` | no bound |

Since `min_alpha mse = 1 - cos^2(theta)` and any alpha lands at or above that
floor, Ibrahim's `0.9095` forces `cos(theta) >= 0.3008` with no further
assumption. Our measured transfer cosine is `~0.08`.

**He has at least 3.8x our directional accuracy.** That is the whole game, and
it is not a scaling trick -- he predicts which way genes move far better than we
do.

### Why mse really is closed for us right now

At `cos = 0.08` our floor is `0.9936`, which still misses the `0.989` baseline.
So exp063's verdict is correct and no reweighting, blending or alpha change can
score us a single point of mse. STRATEGY's operational conclusion stands.

What is wrong is the word "permanently", and the inference that the metric is
therefore not worth thinking about. The threshold is `cos > 0.105` -- a 1.3x
improvement, not the "roughly doubles" the doc claims. And each further step is
worth a lot:

| cos(theta) | mse floor | scaled mse earned |
| ---: | ---: | ---: |
| `0.08` | `0.9936` | `0` |
| `0.105` | `0.9890` | `~0.002` |
| `0.15` | `0.9775` | `~0.014` |
| `0.20` | `0.9600` | `~0.032` |
| `0.30` | `0.9095` | `~0.082` |

Because `score_avg` is a plain mean of six, the last row is `+0.0137` of
overall on its own.

### The strategic consequence

Our entire approach buys pds by spending cosine -- whitening is documented as
lowering transfer cosine 32% while raising pds. We are the best team on the
board at pds (`0.8278` cosine, first) and we are 5th-worst on raw mse. Ibrahim
is *behind us on pds* (`0.8024`) and wins anyway.

So the ceiling on the current approach is visible: we are near the top of what
pds alone can pay. The next real gain is not another pds lever. It is raising
the underlying transfer cosine, which would improve pds *and* open mse at the
same time.
