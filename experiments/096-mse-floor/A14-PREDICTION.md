# a14 — pre-registered prediction

Written 2026-08-29 17:10 UTC, while building. Both of today's slots are spent; this is for the
next window.

**Arm:** a12 exactly, with `VCC_ALPHA` 1.0108 → **0.982**. One variable.

## Why — a12 shipped 2.9% too hot, and the reason is a real bug in our own model

The artifact audit of `a12-tier8000-cap.vcc` found that `claim` **did** reach the correction cap
(0.245917, **99.75%** of it, binding on 720 of 900 rows). The prediction missed for a different
reason: **the scorer applies `min(claim, budget)` inside each context's own run, and `min` is
concave**, so evaluating it at the panel mean is not the same function. At a12's raised
dispersion, context C's claim (0.265) overran its budget (0.2496) — `rho_C = 0.9418` — which was
invisible to the mean form because the mean claim (0.245917) sat just under the mean budget
(0.245927).

The attribution ladder closes to `+0.00117` of the board:

| raw `mse` | Δ | step |
| ---: | ---: | --- |
| 0.888247 | | as registered |
| 0.888619 | +0.00037 | substituting the artifact's own measured terms |
| **0.910448** | **+0.02183** | **the per-context `min` branch** |
| 0.900182 | −0.01027 | re-solving the anchors with the same fix |
| 0.901349 | +0.00117 | **board** |

**Per-context cliffs sit at `alpha` = 0.982 / 1.008 / 1.042, and the curve turns up at the
LOWEST, not at the mean.** a12 shipped at 1.0108, above context A's cliff.

## Predictions, against a12 (`0.16798741`)

| member | predicted Δ | reasoning |
| --- | ---: | --- |
| `mse` raw | 0.9013 → **0.8813** | the per-context cliff, corrected model |
| `score_mse` | **+0.0210** | |
| `nmae` | +0.0003 to +0.0005 | measured slope −0.011/−0.017 scaled per unit alpha; alpha falls 0.029 |
| `pds` | +0.0003 to +0.0004 | measured slope −0.0138/−0.0116 |
| `fid`, `reach`, `jac` | ~0 | no call-set policy change |
| **overall** | **≈ 0.1716** (`+0.0036`) | |

## Falsifiers

1. **raw `mse` outside `[0.86, 0.90]`.** The corrected model fits a8/a10/a11 to `<=0.0001` and
   a12 to `0.0012`; a miss outside this band means the per-context fix was not the whole story.
2. **`score_mse` below a12's `+0.0912`** — i.e. lowering alpha made the member *worse*, which
   would mean a12 was below the binding cliff rather than above it and the sign of the
   correction is wrong.
3. **`nmae` or `pds` fall.** Both are predicted to rise on the measured amplitude slopes. If
   either falls, the monotone-amplitude model for those two members is wrong and the 096
   invariance retraction needs revisiting.

## What this arm is NOT

It is a calibration fix worth `+0.0036`, not a route to first place. The gap to the leader is
`+0.0432` and lives in `nmae` (34.6%), `reach` (34.6%) and `mse` (30.1%). This closes about 8%
of it. The `nmae` study running under `experiments/098-nmae/` is the one aimed at the gap.
