# Pre-measurement predictions — written BEFORE `cd4f080` finished building

Committed while the known-answer row was still emitting context B, so the numbers below were
fixed before any `f = 0.80` artifact existed on disk. This is not a preregistration (that is
`PREREGISTRATION.md`, binding and untouched); it exists so mandate falsifiers 3 and 4 are read
against a stated expectation rather than explained after the fact.

## Falsifier 3 — delivered amplitude, `cd4f080` against `cd4f034`

210 measured, on the shipped artifacts, `f` 0.10 → 0.34 (Δf = 0.24):

```
median delivered |row| ratio   0.99367   (−0.633%)
A = mean|row|^2 ratio          0.97970   (−2.030%)
```

If the drift is **linear in `f`**, then over 0.34 → 0.80 (Δf = 0.46, i.e. 1.9167× the step):

| quantity | linear-in-f prediction |
|---|---:|
| median row-norm ratio, `cd4f080` / `cd4f034` | **0.98787** (−1.213%) |
| `A` ratio, `cd4f080` / `cd4f034` | **0.96110** (−3.890%) |

The mechanism 210 hypothesised (the implicit renormalisation the production renderer applies,
which is mass-ratio dependent) is **not** linear in `f` — the mass ratio moves with the square
of nothing in particular but the renormaliser divides by a sum that itself grows with the
residual. So a **larger** drop than the linear figure is the expected direction of surprise.
I am predicting linear and recording that super-linear is the plausible miss.

Pass condition adopted here: the ratios must be **below 1** and the drop must **not exceed 3×
the linear prediction** (median ratio ≥ 0.9636, `A` ratio ≥ 0.8833). Beyond that the arm is no
longer a shape-only change in anything but name and the coordinator must be told loudly.

## Falsifier 4 — non-vacuity `‖ΔE‖/‖E‖` against the `f = 0.34` artifact

207 measured `‖delivered(f) − delivered(0)‖ / ‖delivered(0)‖` at **0.1023** (`f = 0.10`) and
**0.3312** (`f = 0.34`). The residual is linear in the fraction, so

```
slope    = (0.3312 − 0.1023) / 0.24  = 0.95375 per unit f
at 0.80  = 0.1023 + 0.95375 × 0.70   = 0.7700
predicted ‖d(0.80) − d(0.34)‖/‖d(0)‖ = 0.7700 − 0.3312 = 0.4388
```

210's analogue of this prediction was 0.2289 and it measured **0.2351**, **+2.7%** on the shipped
frame. Applying the same +2.7%:

| quantity | prediction |
|---|---:|
| `‖ΔE‖/‖E‖`, `cd4f080` vs `cd4f034` | **≈ 0.45** (0.44 linear, 0.451 with 210's +2.7% frame offset) |
| per-block minimum | **> 0.10** (210 measured 0.0706 on a 2.3× smaller step) |
| blocks with displacement exactly 0 | **0 / 900** |

The minimum is the number that matters: 200 and 207 both found the CD4 residual can be
identically zero on a block while every aggregate check passes, because the sizing rule takes a
median over rows. A mean or a median would hide a single silent zero; the minimum cannot.

## Falsifier 5 — challenge `n_pred`

`rc_callset` and `cd4f034` both read **2,599** in all three contexts. Prediction: **2,599 again**
(±0%), because `rc_callset` prescribes `top 400 by |H| ∪ top 2,200 by control CPM` and the call
set's *size* is prescribed, not discovered. The head *membership* must move further than 210's
(min 184 / median 314 of 400) because the dose step is 1.9× larger — predict median head overlap
**below 314**, plausibly ~250–290.

## Meta

`solved_fraction` = 0.80 / 1.0556645369046553 = **0.757817672…**, dose ratio against `cd4f034`'s
0.32207201067578145 = **2.352941×** exactly (= 0.80/0.34).
