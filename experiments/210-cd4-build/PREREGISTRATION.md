# PREREGISTRATION — `cd4f034`, written 2026-09-03 18:30 UTC, BEFORE any build

Coordinator (lead). Committed before the builder runs. Not editable after the board result.

## 1. The arm — ONE variable against the current incumbent

```
VCC_CD4_TARGET_F : 0.10  ->  0.34        (3.4x)
emitter          : rc_callset            (unchanged - the incumbent)
everything else  : rc_callset, unchanged
```

Incumbent = `rc_callset`, entry `31oB9yWy8CSNBDEkY4iW`, board **0.1918638367**, rank 11 of 643.

Priced against `rc_hybrid` by experiment 207; this build sits on `rc_callset` instead, which is
`rc_hybrid` plus an **emission-side** call-set union. 207's arm is a **profile** change. The two
channels are orthogonal (`AGENTS.md` §8.1), and `rc_callset` proved on the board that the emission
side leaves `pds`/`mse`/`nmae` bit-identical — so 207's price transfers, and this remains **one
variable against what we actually ship**.

## 2. Why

**This is the first DIRECTION arm in the project's history to clear the gate.** §8.2: `mse`+`pds`+
`nmae` are 135% of the remaining gap to first, and every internal route had closed — 166 (no fifth
line), 184 (adaptation on 50 lines), 193 (the operator), 194 (eleven estimators), 195 (learned
transforms), 203 (hierarchical shrinkage), 204 (donor weighting, decisively), 205 (X-Cell/Pisces do
not exist).

206 then argued the exit is a **lineage-matched donor** and that the one we hold was dosed against
the wrong objective. 207 measured it over a 25-point sweep on PLOLO (4 folds, challenge frame, 744
scored configurations) and **the three maximisers genuinely differ**: cosine at `f = 0.80`
(**+31.5%** over shipped), `q_w` at `f = 0.80`, `pds` at `f = 0.20`. **The shipped `f = 0.10` is
below all three** — because 062 chose the dose by maximising LODO `pds`, the member 197 proved is
least sensitive to direction.

Two of my own premises were wrong and 207 corrected both, which strengthens the case:
- the residual does **not** ship at effective fraction 0.0528 — that is a builder *diagnostic*. a14
  shipped 0.05, **a15 shipped `f = 0.10`**. So 197's measured **+6.8% cosine for a14→a15 was itself
  a 1.89x dose increase**: we already hold a two-point board dose-response curve with positive slope.
- h1-test **cannot price this arm at all** (CD4 shares 4 of its 85 targets; the median-based sizing
  rule then makes the residual identically zero). 200 found it, INBOX-1 flagged it, 207 verified it
  and relocated the guard to `plolo.solve_fraction`.

## 3. Operating point: `f = 0.34`, not the cosine maximiser

Three criteria stated before the arm was chosen, all selecting 0.34: argmax of the measured
six-member price (+0.00437, against +0.00413 at 0.28 and +0.00408 at 0.40); modal
leave-one-fold-out choice; and **the largest dose at which no priceable member is resolvably worse**
(`|Δpds| = 1.23σ` against 162's Gate-1 threshold of 1.645σ; `f = 0.40` fails at 2.43σ).

I am deliberately **not** shipping the cosine maximiser `f = 0.80`. Taking 0.34 forgoes roughly half
the available cosine to keep `pds` inside its own noise floor — this arm's job is to put a **clean,
attributable direction point on the board** and calibrate the cosine→board conversion, not to
maximise one submission.

## 4. Predictions

| member | mechanism | predicted Δ overall vs `rc_callset` |
|---|---|---:|
| `mse` | 096/197 identity, cosine **+18.0%** at this dose, P(>0) = 1.000 | **+0.005206** |
| `pds` | PLOLO, CI spans zero | −0.000685 |
| `nmae` | 188 identity, wald89 pool | −0.000152 |
| `fid`/`reach`/`jac` | ordering precision improves +0.00757; carried at 0 | 0 |
| **total** | | **+0.00437** conservative · **+0.00572** board-calibrated |
| out-of-sample (leave-one-fold-out) | | +0.00340 / +0.00502 |

**Predicted board overall 0.19623 (conservative) to 0.19758 (calibrated).** Gate: **+0.0019**
shape-only — and shape-only is *measured*, not assumed: delivered median row norm 1.0009x,
delivered `A` within 1.0%, and the α re-solve is a **literal no-op** (α\* = 0.982 exactly), so
`VCC_ALPHA` stays at 0.982 and must **not** be re-solved. Clears by **2.3x** in-sample, **1.8x**
out-of-sample.

207's own per-member raw prediction, made against `rc_hybrid`: raw `mse` **0.84413**, raw `pds`
**0.83414**, transfer cosine **0.09350**.

## 5. Falsifiers — pre-committed

1. **raw `mse` does not improve.** It is 84% of the predicted gain and P(>0) was 1.000. If raw `mse`
   is not below `rc_callset`'s **0.8738045**, the 096/197 identity does not transfer at dose and the
   whole CD4 dose family is retired.
2. **raw `mse` outside [0.836, 0.852].** Point prediction 0.84413.
3. **`pds` falls by more than 0.004 raw** (~2.7x the predicted −0.0011). 0.34 was chosen precisely
   as the largest dose where `pds` is not resolvably worse; a real fall means the sweep's `pds`
   noise floor was understated.
4. **`nmae` raw rises by more than 0.002.** Predicted −0.0001.
5. **`fid`, `reach` or `jac` moves by more than 0.004 scaled.** This is a profile arm; 169 falsified
   rule 13's `reach` toll for that class and `rc_hybrid` confirmed the toll belongs to
   profile-*changing emission* arms. A large DE move would mean the call set moved, which it must not.
6. **Config diff is not exactly 2 fields** against `202-rc_callset.h5ad.meta.json`
   (`cd4_target_f` and its deterministic consequence `effects.cd4.fraction`/`solved_fraction`), with
   `operator.sha256` = `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`,
   `effects.cd4.fully_supported_targets` = **276**, `alpha` = **0.982**, and
   `after_all_median_norm` unmoved.

## 6. What the slot buys even if flat

The first board measurement of a **dose** on the one component 197 identified as the only resolvable
cosine gain in eight board arms. Either it converts — and 207's PLOLO cosine instrument (±0.20%
relative) becomes a validated direction instrument for the rest of the project — or it does not, and
the last open direction route closes with a clean number instead of a hope.

**Mechanism worth recording, because it decides D/E/F:** 207 found CD4 must enter **un-whitened** (a
whitened CD4 residual is flat at every dose), the shuffled-CD4 control is **strictly negative at all
12 doses with no permutation tie**, and the gain is **+7.6% to +15.2% on four unrelated non-T
lines** — i.e. **lineage-generic, so it carries to D/E/F** even though the lineages there are
unknown until Oct 22. The mandate's lineage-specificity test was impossible (`nadig_jurkat` covers
0/300).

## 7. Discipline
Joke model name, **empty description** (`AGENTS.md` §14). Submitted at or after the 00:00 UTC slot
reset; both 2026-09-03 slots are spent (`rc_hybrid`, `rc_callset`).
