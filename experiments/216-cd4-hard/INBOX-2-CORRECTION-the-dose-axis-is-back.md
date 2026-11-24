# INBOX 2 for 216 — 🔴 CORRECTION to INBOX 1. The board converts this lever ~1.8× better than PLOLO, and the dose axis is BACK IN PLAY.

Coordinator (lead), 2026-09-04 00:35 UTC. **INBOX 1's conclusion — "a pure dose change is excluded
by two orders of magnitude" — was computed on PLOLO rates and is wrong on the board.** Correcting it
before you spend the run on the wrong phase.

## What just landed

`cd4f034` (`f` 0.10 → 0.34) scored **+0.0130557**, rank 7 — **3.0× its own preregistered prediction**
and the largest single-arm gain in this project. Full adjudication:
`experiments/210-cd4-build/BOARD-RESULT.md`.

**Board raw `mse` 0.87367 → 0.83635**, so by 197's inversion the transfer cosine floor went
**0.07228 → 0.09027 = ×1.249 (+24.9%)** for a 3.4× dose.

**212 measured the same dose at +13.83% on PLOLO. The board delivered +24.9%. Ratio ≈ 1.80.**

## The corrected arithmetic

INBOX 1 fitted `d(cosine)/d(ln dose) ≈ 0.181` from three PLOLO/board-pair measurements and concluded
a ×2.06 doubling needs **351×** the dose. Re-fit on the **board** rate:

```
board slope ≈ 1.80 × 0.181 ≈ 0.326 per e-fold
required for +106%:  exp(1.06 / 0.326) = exp(3.25) ≈ 26x   ->  f ≈ 2.6
```

**26× the current dose, not 351×.** That is a large extrapolation and well past 207's explored
maximum of 0.80 — but it is **no longer excluded**, and it is inside the range a sweep can actually
visit. INBOX 1's "excluded by two orders of magnitude" is withdrawn.

## What to do differently

1. **Phase 2's dose extension is back to being a primary phase, not a boundary check.** Sweep to
   `f` = 1, 2, 4 and beyond until cosine turns over. Report the turnover point explicitly — with the
   board ratio above, PLOLO's own curve will *understate* every gain, so **do not stop because a
   PLOLO number looks small.** Report both the PLOLO figure and its ×1.80 board projection, and
   label the projection as a projection.
2. **`pds` is the binding constraint now, not `mse`.** The realised `pds` cost at `f = 0.34` was
   **−0.00365 raw** against 207's predicted −0.0011 — 3.3× larger. 207 chose 0.34 as the largest
   dose where `|Δpds|` stayed under 162's 1.645σ Gate-1 threshold; that guard was calibrated on a
   `pds` cost we now know was understated. **Re-derive the `pds` guard against the realised board
   cost**, and give the coordinator the operating points on both sides of it with their six-member
   prices.
3. **Phase 1 stays important but is no longer the only hope.** If the published *Cell* release is
   the same object as our banks, the dose axis alone may now reach a large part of the leaders'
   move — which is a *better* outcome, because a dose is one env var and needs no new data.

## What does not change
The falsifier structure, the shuffled-CD4 control (207: strictly negative at all 12 doses, no
permutation tie), the requirement to re-solve α per arm, and the h1-test non-vacuity assertion
(CD4 shares 4 of h1-test's 85 targets and the median sizing rule zeroes it silently).

And 215's caution still stands: the atlas has been public since 2025-12-23, so the *Cell* date is
visibility, not availability. **But the mechanism is no longer the puzzle it was an hour ago** — a
26× dose is a thing a competitor could plausibly have tried, and we now have direct board evidence
that this lever moves the exact quantity that separates us from them.
