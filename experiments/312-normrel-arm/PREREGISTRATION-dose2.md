# 312b — Preregistration for the SECOND dose, frozen before the entry exists

**The arm.** `VCC_NORMREL_EXP = -3.0`, everything else identical to entry `Q519I2pLWRCW7PAEx5ul`
(exp −1.5, **0.21771888483638538**). Note w ∝ r^−3 is the **SQUARE** of the r^−1.5 weight ratios, not twice the
dose (Astra's correction). Reference for Δ is the **−1.5 entry**, not `cd4svd5`.

## What −1.5 delivered, as the anchor
Δ vs `cd4svd5` **+0.0016245**; members `pds` **+0.01230**, `mse` +0.00536, `nmae` **−0.00771**, `fid` +0.00088,
`reach` −0.00113, `jac` +0.00004. Implied transfer cosine 0.09728 → 0.09973.

## Prediction, bins fixed now (Δ₂ = this entry − 0.21771888483638538)
| bin | reading |
|---|---|
| **Δ₂ ≥ +0.0008** | the curve is still climbing at the square of the ratios; `pds` outruns the `nmae` tax. Push the exponent further and re-price the family. |
| **−0.0008 < Δ₂ < +0.0008** | flat: −1.5 is at or near the plateau. Stop moving the exponent; the remaining lever is the WEIGHT SHAPE, i.e. Astra's `max(1, r)^−p`. |
| **Δ₂ ≤ −0.0008** | the curve has turned over. `nmae` (and/or `reach`) grows faster than `pds`. −1.5 is past-optimal or near it, and a *weaker* exponent is the next probe. |

**Mechanistic sub-reading, no verdict:** if `pds` keeps rising while `nmae` falls faster, the arm is trading
retrieval against signed accuracy — 143's precedent (`pds` +0.00279 with signed credit down across every family).
That would say the two mechanisms `r^−p` conflates (suppress oversized rows vs up-weight tiny ones) must be
separated before any further dose work.

## Stated in advance
- This is the **last free slot of 2026-09-07**. Both of today's slots are spent on one preregistered dose curve.
- The board shows the LATEST entry: if Δ₂ < 0 the displayed row falls below 0.21772. Expected, costless, **not to
  be re-submitted to restore rank**.
- Entry `Q519I2pLWRCW7PAEx5ul` is already captured at `data/board-entries/`, so its record survives.
- **This dose curve is measured on the SHIPPED implementation**, which 313 is concurrently testing for parity with
  308's construction. Whatever 313 finds, both points here share one implementation, so the curve is internally
  valid; it just may not describe 308's object.
