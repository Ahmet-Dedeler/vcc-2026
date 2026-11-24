# ADDENDUM to the `cd4f034` preregistration — written 2026-09-03 20:55 UTC, BEFORE submission

The preregistration itself is **not edited**. This addendum records two things that became known
from the *build* (not from any board result), and the classification decision, so both are on the
record before the slot is spent.

## 1. Falsifier 3b FAILED, on a mechanism that is benign and identified

Required: delivered `A` within **1.0%**. Measured: **−2.03%** (median row norm 0.99367×,
`A` 0.97970×). It is real — it survives on quantisation-immune gene sets (0.9682 on the top 2,000
by CPM, where emission noise is 0.00% of `A`).

**Effect-space amplitude is untouched**: the post-preserve median row norm is bit-identical to
sixteen digits. So the divergence is entirely **downstream of the preserve**, in the render. The
build's diagnosis: 207's `amp207.py` renders through 156's `minnorm` mass fix, while production
ships `VCC_MASS_FIX` **unset**. 210 did not re-classify the arm and correctly left that to me.

## 2. Consequence: falsifier 2's range will fire, on the FAVOURABLE side

A lower delivered `A` lowers the `s²A` term, so it makes `mse` **better** than predicted. On 207's
own conversion the extra uncredited favourable `A` channel is **−0.0108 raw**, putting predicted raw
`mse` near **0.8333** — below preregistration falsifier 2's `[0.836, 0.852]`.

**Recorded now, before the result: I expect falsifier 2 to fire low, and that is not a failure of
the hypothesis.** The point prediction that should be judged is **raw `mse` ≈ 0.8333 ± 0.004**, and
the falsifier that still means something is #1: *if raw `mse` does not improve on `rc_callset`'s
0.8738045 at all, the 096/197 identity does not transfer at dose and the CD4 dose family is retired.*

## 3. The classification decision, and it is made on a board point, not on a label

Falsifier 3b's failure removes the *measured* basis for calling this "shape-only", and the two gates
differ materially: **+0.00437 clears +0.0019 but not +0.0089.** So the classification decides the
slot.

**I am classifying it shape-only, on evidence from the same knob.** Rule 13's thresholds are the
realised RMS of the **unpriced** channels (`fid`, `jac`, `reach`) for an arm class. We hold a board
point for a CD4 dose increase on exactly this knob — **a14 (`f` = 0.05) → a15 (`f` = 0.10), a 1.89×
dose step** — and its unpriced channels moved:

| member | a14 → a15, scaled |
|---|---:|
| `fid` | **+0.002189** |
| `reach` | **+0.000826** |
| `jac` | −0.000007 |
| **unpriced sum** | **+0.003008 scaled = +0.000501 overall, favourable** |

**+0.0005 favourable, not ±0.0089.** The +0.0089 threshold was calibrated on arms that moved α or
the emitted dispersion — global amplitude knobs. This arm moves **neither**: `alpha` is 0.982
exactly (207 §6.2 measured α\* = 0.982, so the re-solve is a literal no-op), the emission law is
`rc_callset` unchanged, and the challenge call set is **bit-identical** (`n_pred` 2,599 in A, B and
C, +0.00%), which pins `jac`'s denominator structure and `fid`'s coverage term by construction.

This is the same form of argument the `rc_hybrid` board point forced on rule 13's `reach` toll: a
class-level constant must be established by a board point in that class, not by a label. Here the
board point exists, it is the same knob, and it says the unpriced channel is small and favourable.

**Ship. Predicted +0.00437 conservative / +0.00572 calibrated against +0.0019.**

## 4. What would make me wrong, stated now

If `fid`, `jac` or `reach` moves by more than **0.004 scaled** (preregistration falsifier 5), the
a14→a15 pair did not generalise to a 3.4× dose and the classification was too aggressive. That is
the number to check first when the board returns, before reading the headline.

## 5. Hand-off, not for this slot
The `A` divergence means `amp207.py` and production render through different mass-fix settings.
That is a measurement-infrastructure defect worth fixing — every future amplitude prediction from
that script inherits it — and it belongs to whoever next touches 156's render path. **It does not
affect this build**, which was produced by the production builder with `VCC_MASS_FIX` unset, as a15
and every board arm since have been.
