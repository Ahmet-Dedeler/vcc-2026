# 316 — Preregistration, frozen BEFORE the entry exists

**The arm.** `VCC_NORMREL_EXP=-1.5 VCC_NORMREL_PM=blend`, package `submissions/316-normrel15-pm.vcc`.
Its config differs from the live entry `Q519I2pLWRCW7PAEx5ul` (**0.21771888483638538**) in **exactly one key**,
`normrel_pm` — verified: no other key added, removed or changed. **Reference for Δ is that entry, not `cd4svd5`.**

## What the arm changes, measured not asserted
313 proved the shipped arm is amplitude-changing: delivered non-panel row norms **0.806–2.155**, whole matrix
**1.060**. 316 restores exactly shape-only: whole-matrix norm ratio **1.0000000**, max |delivered-norm ratio − 1|
**4.44e-16**. It is **not** a cleanup — against the shipped arm, max |diff| **2.087**, row cosine median 0.9932,
min 0.8736. All six known-answer rows pass, two **bitwise 0.000000e+00** (unset reproduces the shipped arm; blend
with exponent 0 reproduces `cd4svd5`).

## The prediction
The shipped arm's `nmae` leg was **−0.00771** where the preregistration expected **−0.00095 (a benefit)**. 313's
explanation is that an amplitude-changing object loses `nmae` (AGENTS 8.1d). **If that explanation is right,
restoring shape-only should recover most of that leg.** Full recovery alone would be worth ≈ **+0.00129** overall.

| bin (Δ vs 0.21771888) | reading |
|---|---|
| **Δ ≥ +0.0010** | 313's mechanism is confirmed on the board: the `nmae` loss WAS the accidental amplitude change. Ship this form; re-price the family with the parity ordering. |
| **+0.0002 … +0.0010** | partial recovery — the amplitude change explains some of the leg, not all. |
| **−0.0002 … +0.0002** | neutral: the parity difference does not reach the six-member composite, even though it is a real re-shaping (max |diff| 2.09). |
| **Δ ≤ −0.0002** | **the parity-correct form is WORSE.** The accidental amplitude change was *helping*. That would be the important negative — it would say the +0.0016 we banked is not the weighting mechanism we think it is, and 308's construction is not the thing to chase. |

**Declared secondary, no verdict:** the `nmae` member specifically (shipped +0.10277, incumbent +0.11047 — does it
move back toward the incumbent?); `mse` and `pds`; and the implied transfer cosine (incumbent 0.09728, shipped
0.09973) under §4.1c's validity guard.

## Stated in advance
- The board shows the LATEST entry. If Δ < 0 the displayed row falls below 0.21772. Expected, costless, **not to be
  re-submitted to restore rank**.
- `Q519I2pLWRCW7PAEx5ul` and `1n5LFk6PSEHq15UHENIE` are already captured under `data/board-entries/`.
- **This is one variable.** D5 (weight support) and D7 (PerturbMatch column policy) are frozen per 313 §8 and are
  NOT swept here; 316's arm differs from 313's own `parity_variant` by that freeze alone (max |diff| 0.400).
- Slot: the **first** of 2026-09-08, spent as soon as the allowance resets. The second stays open for the reading.
