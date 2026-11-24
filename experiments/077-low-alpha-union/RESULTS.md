# 077 — low-alpha union sweep

Date: **2026-08-23**. Run `20260823T060938Z`, worker `i-0e4070204e83188fb`
(`r7i.4xlarge`), wall clock `3:41:46`, exit `0`, peak process-tree RSS
`59.92 GiB`. Comparator: Orion HCT116 count-level, union axis, challenge panel
excluded — the same harness as exp065.

## Verdict

**Negative. Do not lower alpha.** Neither arm passes, and the promotion gate
fails on both counts.

| arm | avg6 | delta vs `a=1.0796` |
| --- | ---: | ---: |
| `a=0.35` | `0.11894` | `-0.00270` |
| `a=0.40` | `0.12018` | `-0.00146` |
| baseline `a=1.0796` | `0.12164` | — |

`adjacent_alphas_both_positive: false`, `passes_promotion_gate: false`,
`positive_mse: false` on both. The loss is concentrated in `nmae` (`-0.0289` at
`a=0.40`), partly offset by `fid` (`+0.0209`). pds is essentially flat.

## What it actually bought: the tightest measurement of our distance to mse

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



This is the result worth keeping. Lowering alpha drives raw mse down hard:

| alpha | raw `expr_mse_unbiased_capped_norm` | scaled |
| ---: | ---: | ---: |
| `1.0796` (production) | `2.0537529` | `0` |
| `0.35` | `0.9968821` | `0` |
| `0.40` | **`0.9929697`** | `0` |

The threshold for nonzero scaled mse is `0.989`. **We miss it by `0.00397`.**

Since `min_alpha mse = 1 - cos^2(theta)` and `a=0.40` sits at or near that
floor, this independently measures our transfer cosine:

```
cos(theta) = sqrt(1 - 0.9929697) = 0.0838
```

which confirms `STRATEGY.md:62`'s `~0.08` from a completely different direction.

## The corrected strategic picture

`STRATEGY.md:63` says mse "stays 0 until directional accuracy roughly doubles".
That is too pessimistic. The requirement is:

```
cos needed = sqrt(1 - 0.989) = 0.1049
cos we have                  = 0.0838
factor required              = 1.25x
```

**A 25% improvement in transfer cosine opens the metric**, not a doubling. For
scale, Ibrahim Mansour's board raw mse of `0.9095` forces `cos >= 0.3008`,
which is `3.6x` ours.

The two facts sit together: exp063 and this experiment both confirm that *no
rescaling of our current effect directions can score mse* — we are already at
the floor and the floor is `0.004` too high. The only lever is the direction
itself. See `docs/V9-RESULT.md`.

## Consequence for exp086

Experiment 086 (sparse-responder emission) was preregistered as conditional on
this sweep producing a viable low-alpha operating point. It did not.
**exp086's precondition has failed and it should not be launched.**

## Preserved artifacts

`s3://<your-s3-bucket>/output/exp077-low-alpha-union/20260823T060938Z/`
and `cloud-output/20260823T060938Z/` locally. The two 2.3 GB prediction h5ads
were not downloaded; they are reproducible from the pinned grid.
