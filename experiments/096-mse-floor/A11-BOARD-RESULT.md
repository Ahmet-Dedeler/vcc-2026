# a11 board result — the floor is cleared, and the cliff was where the model said

Entry `oIL3tgNUMocufeha4Kys`, submitted 2026-08-28 15:21 UTC.
**Overall `0.16611503`** against a10's `0.15418408` — **`+0.01193`**, our largest single-arm
gain in the project, and **rank 16 → 7 of 385**.

One variable: `VCC_ALPHA` 1.0796 → 0.838.

| member | a10 | **a11** | Δ | **predicted Δ** |
| --- | ---: | ---: | ---: | ---: |
| `pds` | 0.736561 | 0.739889 | **+0.003328** | 0.0000 ± 0.005 |
| **`mse`** | **0** | **0.082041** | **+0.082041** | **+0.0903** |
| `nmae` | 0.094338 | 0.096982 | +0.002644 | −0.0150 |
| `fid` | −0.027897 | −0.031535 | −0.003638 | small |
| `reach` | 0.123595 | 0.112220 | −0.011375 | −0.0214 |
| `jac` | −0.001493 | −0.002907 | −0.001414 | small |
| **overall** | **0.154184** | **0.166115** | **+0.011931** | **+0.0089** |

Raw: `mse` **1.352415 → 0.909989**. `pds_cosine` 0.831773 → 0.833307. `reach` 0.188226 → 0.178243.

## Against the pre-registered falsifiers

**1. "`score_mse` returns 0" — did not fire.** It returned `+0.0820`. This was the one that
mattered: it would have meant the cliff was in the wrong place and the whole 096 frame wrong
rather than mis-sized.

**2. "raw outside `[0.87, 0.94]`" — did not fire.** Predicted `0.9021`, measured **`0.9100`**,
an error of `+0.9%` on a 33% move of the member. That is the **second** out-of-sample
confirmation of the decomposition: it predicted a10's raw at `1.3567 ± 0.01` (actual
`1.3524`) before a10 was scored, and a11's at `0.9021` (actual `0.9100`) before a11 was built.

**3. "`pds` moves by more than 0.005" — did not fire**, at `+0.0033`. See below.

**4. "overall below a10" — did not fire.** It rose by `+0.0119`, against `+0.0089` predicted.

## The amendment was wrong, and wrong in sign

Three hours before submitting, the independent amplitude sweep returned `pds` rising
monotonically with amplitude (scaled 0.618 → 0.650 → 0.667 across α 0.30 → 0.60), and I
amended the prediction to expect **`pds` to fall by about `0.012`** at a 22.4% amplitude cut,
writing that falsifier 3 was "now expected to fire."

**`pds` rose, by `+0.0033`.** The amendment was wrong in sign, and the original
scale-invariance reasoning it doubted was right.

Why the sweep misled: it is the `031-somatic` model (donor transfer + whitening, no
PerturbMatch, no CD4) and its `pds` runs 0.618-0.667 against our board 0.7366. Those points
sit in a steep low-amplitude regime that our production pipeline is nowhere near — the sweep's
α = 0.60 is below a11's effective operating point, not above it. **Extrapolating a slope
measured on a different model, in a different regime, past the range where it was measured,
was the error.** The amendment is left in `A11-PREDICTION.md` unedited; it is more useful as a
record of a wrong call made in the open than as a corrected sentence.

## Where the +0.0119 came from, against the +0.0089 forecast

The model under-predicted, and both misses are in our favour:

* **`reach` cost half what was budgeted** — `−0.0114` against `−0.0214`. The budget used the
  steepest slope the board supports (the α 0.40 ↔ 1.0796 pair), explicitly flagged as an
  over-estimate because that pair is confounded. It was.
* **`nmae` moved the wrong way and gained** — `+0.0026` against `−0.0150` predicted. Shrinking
  an over-large predicted fold-change toward the truth improves a fold-change error metric;
  the α ladder that predicted a loss was fitted across builds that differ in more than α.

Net DE cost `−0.0138` against `−0.0364` budgeted, and `mse` delivered `+0.0820` against
`+0.0903`.

## What is now established

**The expression-error member is not a wall, and it never needed a better model to move.** We
sat at exactly 0 on it for thirteen submissions — as 376 of 385 teams still do — while the
cause was that our emitted cells were 0.67× as dispersed as the reference's, so we forfeited
most of the sampling-noise credit the metric *offers*, and our amplitude sat above the ρ
cliff. One environment variable was worth `+0.0119`.

The remaining plan is unchanged and now rests on a twice-confirmed model: raise the emitted
dispersion to the `#247` cap and trim α to 1.012 (`a12`, built and checksum-verified, raw
`mse` predicted `0.878`), then the donor-agreement shrinkage.

**The honest ceiling still stands:** calibration alone tops out near raw `0.870` against the
leader's `0.8270`, whose aggregate directional cosine is at least 1.28× ours. This buys rank
3-5. First place needs the direction, not the calibration.
