# The amplitude sweep, completed — the U-shape is real, and the mid-course amendment was the error

Six configs of the `031-somatic` model (donor transfer + whitening) on
`orion-hct116-thinned`, 132 real perturbations, scored by the real `cell-eval2` with the
`vcc2026` preset. Launched before a11 was built; finished after a11 was scored.

| `alpha` | `mse` raw | `mse` scaled | `pds` scaled | `reach` scaled | `avg6` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.30 | 1.0144 | 0 | 0.6180 | **+0.0367** | 0.1316 |
| 0.45 | 1.0037 | 0 | 0.6503 | +0.0243 | 0.1330 |
| **0.60** | **0.9935** | 0 | 0.6667 | +0.0098 | **0.1332** |
| 0.80 | 1.2817 | 0 | 0.6692 | +0.0008 | 0.1328 |
| 1.08 | 1.8841 | 0 | 0.6696 | −0.0077 | 0.1301 |
| 1.40 | 2.7988 | 0 | 0.6675 | −0.0096 | 0.1204 |

## 1. The U-shape is confirmed, and it is violently asymmetric

`mse` raw falls to a minimum near `alpha = 0.60` and then **triples** by 1.40. The
pre-registration (`PREDICTION.md`) predicted exactly this shape — an interior minimum with the
argmin in `[0.25, 0.55]` and a value of `0.93-0.97`. Measured: argmin `~0.60`, value `0.9935`.
Shape right, position right to within one grid step, level about 3% high.

The asymmetry is the operational point. Being 2x above the optimum costs `+1.9` of raw `mse`;
being 2x below costs `+0.02`. **On this member the penalty for over-shooting amplitude is two
orders of magnitude worse than for under-shooting**, which is why thirteen of our submissions
scored exactly 0 on it.

## 2. `pds` saturates above `alpha = 0.6` — which is why a11's `pds` did not fall

`pds` scaled: 0.618 → 0.650 → 0.667 → **0.669 → 0.670 → 0.668**. It rises steeply below 0.6
and is flat to three decimals above it.

Three hours before submitting a11 I saw only the first three points, read a monotone rise, and
amended the prediction to expect `pds` to **fall by ~0.012** at a 22.4% amplitude cut, writing
that the falsifier was "now expected to fire." **a11's `pds` rose by `+0.0033`.**

The completed curve says why, and it is not the reason I gave post-hoc either. I attributed
the miss to the sweep being a different model in a different regime. That is true but
secondary: **the real error was extrapolating a slope from three points that all sat in the
steep part of a curve that flattens.** a11's move — 1.0796 → 0.838 — lies entirely inside the
plateau. The original scale-invariance argument was right, the three-point amendment was
wrong, and the full curve was already running when I made it.

**The pre-registration outperformed the mid-course correction.** That is the whole case for
rule 2 in `docs/RESEARCH-METHOD.md`, demonstrated against myself.

## 3. `reach` moves the OPPOSITE way here from the board — the `n_real` trap again

On this harness `reach` is **best at the lowest amplitude** (`+0.0367` at 0.30, falling to
`−0.0096` at 1.40). a11 cut amplitude and **lost** `reach` on the board (`0.1236 → 0.1122`).
Same axis, opposite sign.

This harness has `n_real = 92` against the challenge's ~4,757, so `reach`'s normalisation is
unreadable here exactly as `docs/RESEARCH-METHOD.md` rule 1 states. Had a11 been chosen on
this sweep's `avg6` column — which peaks at `alpha = 0.60` — we would have cut amplitude far
harder than 0.838 and paid a much larger `reach` bill for `mse` that was already saturated at
the clamp.

**The sweep was worth running and would have been actively harmful to steer by.** It confirms
mechanism and shape; its levels and its `reach` sign do not transfer. The arm that was shipped
came from the `msekit` decomposition against the board's own anchors instead, and predicted
a11's raw `mse` to 0.9%.

## 4. `mse` scaled is 0 in every row

Every config here scores 0 on the member, including the `alpha = 0.60` row at raw `0.9935` —
because the harness's own baseline sits below it. A harness cannot tell you whether you have
cleared the floor on the challenge. Only the board can, and now has.
