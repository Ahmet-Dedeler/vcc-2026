# Pre-registered predictions for the over-dispersion arms (a7, then a8)

> The filename says a7 because that arm was written up first. **The arm actually submitted
> is `a8-nb-phi4`** (φ=4), and its prediction is the section headed "a8 prediction (φ=4),
> the one actually submitted". a7 (φ=2) was built, probed at `n_pred` 4,523 and retired
> unused when the finer sweep put the optimum at φ=4.

## a7 — prediction, written before the submission is scored

Arm `a7-nb-phi2` = `a2-prune` (board `0.1503084`) plus `--tail-spread 0.25
--overdispersion 2.0`. One variable from our best board artefact.

Recorded before scoring so the result is a test rather than a story. This project has
mis-predicted the board three times (exp093 `+0.0202 → −0.0029`, exp091 `+0.00088 →
+0.00013`, exp095's ordering counterfactual `+0.00583 → caught before spending`), and every
one of those was reconstructed after the fact.

## The prediction

| member | a2-prune (measured) | **a7 predicted** | Δ predicted |
| --- | ---: | ---: | ---: |
| `pds` raw | 0.831773 | **0.831773** | **exactly 0** |
| `mse` raw | 1.5708 | ~1.4 (clamped, scores 0 either way) | 0 scaled |
| `nmae` raw | 0.944139 | **0.944139 ± 0.0002** | ~0 |
| `reach` raw | 0.170083 | **0.224** | +0.054 |
| `fid` raw | 0.502094 | **0.5135** | +0.0114 |
| `jac` raw | 0.031430 | 0.0318 | +0.0004 |
| **overall** | **0.1503084** | **≈ 0.1671** | **+0.0168** |

`pds` and `nmae` are the load-bearing predictions. The exact-total repair makes each gene's
panel total bit-identical whatever the scatter, so if **`pds` moves at all** the mechanism
is not what I think it is and every number here is void. v14 already demonstrated this on
the board for a full-range `spread` change (`pds` moved by exactly `0.000000`); this is the
same guarantee under a wider scatter.

## Where each number comes from

* **reach**: harness production `0.0838 → 0.1378` under this exact policy, realised through
  the production emitter. Applied **additively** to the board's `0.170083`. Additive rather
  than proportional because proportional (`×1.64`) gives `0.279`, essentially the board's
  all-time maximum of `0.2824`, which is implausible for a one-parameter change.
* **fid**: harness precision `0.5465 → 0.5644` (`+3.3%`) applied to the board's `0.502094`,
  then multiplied by a coverage of `~0.99` at the projected `n_pred` of ~4,880. The coverage
  curve is the one inverted from our own three board points.
* **jac**: harness `+0.0004` raw, taken at face value; it is small either way.

## What would falsify what

| observation | conclusion |
| --- | --- |
| `pds` moves by more than 1e-5 | the exact-total repair does not survive over-dispersion; retract §10 of exp095 |
| `nmae` moves by more than 0.001 | same |
| `reach` moves less than `+0.02` raw | harness reach deltas do not transfer additively; the 4.4× v14 calibration was luck |
| `fid` falls | coverage bites harder than the inversion says, i.e. challenge `n_real` is larger than our three board points imply |
| overall lands near `+0.017` | the realised-emission methodology is calibrated, and `realize.py` becomes the gate for all future emission work |

**The last row is the one worth the slot.** A validated local gate for the emission channel
is worth more than this single submission, because it is what makes the next twenty cheap.

## Addendum — a second calibration point, and why the prediction stays conservative

Two artefacts now have the same policy measured on both surfaces:

| | harness Δreach | board Δreach | ratio |
| --- | ---: | ---: | ---: |
| v14 (`tail 1.0, keep 850`) | +0.0028 | +0.0122 | **4.4×** |
| v15 (`tail 1.0, keep 2000`) | +0.0026 | +0.0104 | **4.0×** |

Consistent to within 10%, from two independent policies. Applied naively to a7's harness
`Δreach` of `+0.0540`, that multiplier predicts `+0.227` — which would put reach at `0.397`
against a board all-time maximum of `0.2824` across 780 entries. **So the multiplier cannot
hold at this magnitude**; it is a small-delta local slope, and the relationship must saturate
somewhere below the ceiling.

The prediction above therefore uses the **additive** transfer (`0.170 + 0.054 = 0.224`),
which is the conservative reading. The honest interval is:

| transfer assumption | board reach | Δoverall from reach |
| --- | ---: | ---: |
| additive (**used above**) | 0.224 | +0.0102 |
| multiplicative, saturating at the board maximum | 0.282 | +0.0211 |

So a7's plausible range is roughly **+0.010 to +0.021 on `reach`**, on top of `fid` and
`jac`. The prediction of `+0.0168` overall sits at the low end deliberately.

## Superseded: φ=4 is the peak, and a7 (φ=2) becomes the fallback

A finer sweep after a7 was already building found φ=2 was **not** the optimum at production
spread:

| φ (spread 0.25) | `n_pred` | reach | Δreach | ΔTOTAL (harness) |
| ---: | ---: | ---: | ---: | ---: |
| 1.0 | 3,528 | 0.1113 | +0.0052 | +0.0101 |
| 2.0 (**a7**) | 3,218 | 0.1378 | +0.0102 | +0.0204 |
| 3.0 | 3,018 | 0.1604 | +0.0144 | +0.0283 |
| **4.0 (a8)** | **2,863** | **0.1675** | **+0.0158** | **+0.0320** |
| 4.0, second seed | 2,880 | 0.1535 | +0.0131 | +0.0298 |
| 5.0 | 2,690 | 0.1418 | +0.0109 | +0.0291 |
| 6.0 | 2,636 | 0.1233 | +0.0075 | +0.0271 |

Reach peaks at φ=4 and falls away on both sides. **Note the seed control**: the two φ=4
draws give reach `0.1675` and `0.1535`, a spread of `0.014` raw against production's
`0.0015`. Over-dispersion adds real emission variance, so the honest φ=4 figure is the mean,
`0.1605`, i.e. `Δreach ≈ +0.0145` overall rather than `+0.0158`.

**a7's probe also calibrated the `n_pred` transfer.** On the real challenge controls φ=2
gives `n_pred` **4,523** (A 4,685 / B 4,878 / C 4,006) against v11's 5,961 — a **24%** cut
where the harness predicted 18%. Applying that same over-shoot to φ=4's harness `−27%`
puts it near **3,870** on the challenge, where the coverage curve inverted from our own
board points still reads ~0.97.

Revised projection for **a8** (φ=4): reach `+0.0145`, fid `+0.0049`, jac `+0.001` →
**≈ +0.020 overall**, against a7's `≈ +0.015`. a7 is built, checksum-verified and held as
the fallback if a8 fails to build in time.

## a8 prediction (φ=4), the one actually submitted

Using a7's **measured** `n_pred` transfer rather than the harness's, and the mean of the two
φ=4 seeds rather than the lucky one:

| member | a2-prune (measured) | **a8 predicted** | Δ |
| --- | ---: | ---: | ---: |
| `pds` raw | 0.831773 | **0.831773** | **exactly 0** |
| `nmae` raw | 0.944139 | **0.944139 ± 0.0005** | ~0 |
| `reach` raw | 0.170083 | 0.2468 | +0.0767 |
| `fid` raw | 0.502094 | 0.5137 | +0.0116 |
| `jac` raw | 0.031430 | 0.0322 | +0.0008 |
| **overall** | **0.1503084** | **≈ 0.1716** | **+0.0213** |

`0.1716` would be **rank 2**, above Ibrahim Mansour's `0.17078`.

Assumptions, each of which can fail independently:

1. **reach transfers additively** (`+0.0767` raw). The two board calibrations say the
   harness *understates* reach deltas ~4×, so additive is the conservative floor; the
   multiplicative reading saturates at the board ceiling of ~0.28.
2. **`n_pred` lands near 3,815** on the challenge — harness `−27%` scaled by the `24/18`
   over-shoot a7's probe measured. The probe on the built artefact checks this before the
   slot is spent.
3. **coverage is ~0.975 there**, from the curve inverted off our own three board points.
4. **precision rises 4.9%**, the harness ratio at φ=4.

If `pds` moves at all, the exact-total repair does not survive φ=4 and everything here is
void — that is the single check that matters most on arrival.

## Footnote: φ=3 and φ=4 are tied once the seed control is applied

The φ=4 single-seed reach of `0.1675` is what made it look like the peak. Its re-seeded
control reads `0.1535`, so the honest φ=4 figure is the mean, **`0.1605`** — and φ=3's
single-seed reading is **`0.1604`**.

| φ | harness `n_pred` | reach (best available estimate) | harness call-set cut |
| ---: | ---: | ---: | ---: |
| 3.0 | 3,018 | 0.1604 *(1 seed)* | −23.3% |
| 4.0 | 2,863 | 0.1605 *(mean of 2)* | −27.2% |

**CORRECTION, measured properly.** The paragraph above compared φ=3's *single* seed against
φ=4's *mean*, which is not a comparison. Re-run with two seeds each:

| φ | seed A | seed B | **mean reach** | mean ΔTOTAL |
| ---: | ---: | ---: | ---: | ---: |
| 3.0 | 0.16044 | 0.15004 | 0.15524 | +0.0272 |
| **4.0** | 0.16751 | 0.15352 | **0.16052** | **+0.0309** |

**φ=4 is genuinely ahead**, by `+0.0037` on the harness total, and the seed-to-seed spread
(`0.0104` at φ=3, `0.0140` at φ=4) is smaller than the gap between them. a8's φ=4 stands on
its own evidence rather than on being the arm that happened to be building.

The `production` control repeated at `−0.00051`, so the noise floor is unchanged.

## a8's probe landed on the projection

Before the slot was spent, `probe_de.py` measured the built a8 artefact against the real
`data/challenge/context_*.h5ad` controls:

| | ctx A | ctx B | ctx C | mean |
| --- | ---: | ---: | ---: | ---: |
| **a8 `n_pred`** | 4,032 | 4,157 | 3,318 | **3,836** |

**Predicted 3,815. Measured 3,836 — 0.5% out.** The chain that produced that number is
worth naming, because it is three transfers deep and each one was measured rather than
assumed: harness `n_pred` at φ=4 (`2,863`, `−27%`) → the over-shoot factor a7's own probe
gave (`24/18`) → the challenge. Coverage at 3,836 reads ~0.975 off the curve inverted from
our three board points, which is the projection's assumption holding.

The decision rule set before the probe was "submit a8 if `n_pred` lands in 3,500–4,200,
fall back to a7 if under 3,000". 3,836 is comfortably inside it, so **a8 is the second
submission** and a7 is retired unused.
