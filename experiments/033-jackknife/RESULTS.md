# 033 — The `expr_mse` jackknife credit, measured, and why alpha should go UP not down

**Date:** 2026-08-22 · measured on `data/harness/orion-hct116.h5ad` and a prediction already
on disk.

## What was wrong with "mse is structurally unreachable"

The claim rested on `min_alpha mse = 1 - cos^2`, which needs cos > 0.105 against a measured
~0.08. That bounds the **uncorrected** squared error. The scored numerator is

```
N_p = (1/G_p) ( ||b_hat_p - b_p||^2  -  rho * min(C_hat_p, C_p)  -  C_p )
```

and `C_hat_p` is the delete-one jackknife estimate of **our own submission's** sampling noise.
The metric pays for noise we emit. The cosine bound never applied to that term, and two teams
now score positive on `mse` (Ibrahim Mansour +0.0728 scaled, Immovable Myhal Rock +0.0308).

## Measured

| quantity | value |
|---|---|
| reference `C_p`, median over perturbations | **64.60** |
| reference `C_ctrl` | 0.415 |
| denominator `D_p = ||b_p − b_ctrl||² − C_p − C_ctrl` | 29.52 |
| the reference's own noise as a share of raw effect energy | **68%** |
| **our even emitter's `C_hat_p`** | **7.59** |
| credited `min(C_hat, C_p)`, as a fraction of `D_p` | 0.257 |

**We collect 7.59 of an available 64.60 — about 12%.** (The even emitter is not exactly
zero-variance: integer allocation still leaves a ±1 spread per gene.)

## The part that inverts a standing assumption

`rho` caps the total credit by `B`, the submission's across-perturbation spread — the
anti-gaming term, so you cannot claim credit for noise without also emitting differentiated
signal. Writing `r` for emitted/true effect magnitude, in the regime where `rho` binds the
total credit is `B ~ r^2 ||delta||^2`, and the numerator becomes

```
N  ~  ||delta||^2 (r^2 - 2 r cos + 1)  -  r^2 ||delta||^2  -  C_p
   =  ||delta||^2 (1 - 2 r cos)  -  C_p
```

**The `r^2` terms cancel.** `N` then *decreases* linearly in `r`, so while `rho` is binding,
**a larger alpha is better for `mse`** — the exact opposite of the `alpha* = <d,delta>/<d,d>`
shrinkage argument that has governed this repo since experiment 001, and the opposite of what
every other metric wanted.

The two regimes meet where `B = sum_p min(C_hat_p, C_p)`, i.e. where the predicted effect
energy equals the available credit. Ours is `4.8^2 = 23` against `C_p = 64.6`, so `rho ~ 0.36`
and we are **deep in the binding regime** — alpha should rise, not fall.

## What this predicts, and what tests it

Two changes, both of which the `quasi` emitter can deliver without touching the mean profile:

1. **Raise emitted variance** so `C_hat_p` reaches `C_p` — up to ~8x more credit.
2. **Raise alpha** so `rho` stops binding — our predicted effect norm would need to go from
   4.8 to ~8.0.

Both trade against the DE metrics, and `nmae` hates large alpha (exp 023 measured `nmae` 1.058
at alpha = 1.6 on `h1-fast`). The grid in `experiments/031-somatic/grid_mse.json` sweeps
`spread` from 0.05 to 1.0 and alpha from 0.35 to 2.5 against all six metrics on the somatic
reference, which is the only way to price the trade rather than argue about it.

**Caveat on magnitudes.** This reference carries 150 cells per perturbation where the
competition carries 400, so its `C_p` is inflated by roughly 2.7x and the credit is
correspondingly larger here than on the board. The mechanism, the sign, and the `rho` regime
are what transfer; the size is not.

---

# CORRECTION (2026-08-22, after the sweep): the credit is a de-biasing device, not a bonus

The measurement above is right and the practical conclusion drawn from it was wrong. Emitting
noise to "collect the credit" does not help, and the sweep says so directly:

| config | raw mse |
|---|---|
| alpha 1.08, even emitter (Ĉ ≈ 0) | **2.259** |
| alpha 1.08, spread 0.35 | 2.309 |
| alpha 1.08, spread 0.7 | 2.396 |
| alpha 1.08, spread 1.0 | 2.532 |

More emitted noise makes `mse` monotonically **worse**.

## Why — the term I did not carry through

`||b_hat_p - b_p||^2` is computed from *our emitted cells*, so it already **contains** our
emission noise: `E||b_hat - b||^2 = ||mu - b||^2 + C_hat_p`. The metric then subtracts
`rho * min(C_hat_p, C_p)`. At `rho = 1` and `C_hat <= C_p` the two cancel exactly:

```
N = ||mu - b||^2 + C_hat - C_hat - C_p = ||mu - b||^2 - C_p
```

— identical to what a zero-variance emitter gets. **The jackknife is there to make `mse`
independent of how noisy your emission is.** That is what "unbiased" means in the metric's
name. At `rho < 1` the noise is *under*-compensated and emitting it is strictly harmful, which
is the monotone degradation above.

A minimum-variance emitter is therefore optimal-or-tied for `mse`, and we were not forfeiting
anything.

## The Discord datapoint that confirms the mechanism

A participant submitted the unperturbed controls unchanged and reported raw `mse` **2.0197**,
where the spec says a control-copy should read ~1.0. That gap is exactly `rho`: all 300 of
their predictions were identical, so the across-perturbation spread `B` was zero, `rho = 0`,
and none of their emission noise was compensated. It is a clean external measurement of the
anti-gaming term.

## What survives

`mse`'s only routes are a better cosine or heavy shrinkage, and shrinkage is a losing trade:
`alpha = 0.35` with spread 0.7 reaches raw `mse` 1.030 — still above the 0.989 baseline, so
still scaled 0 — while the board shows the teams that shrink that far (larry, Vu Mai, both at
raw `mse` ~1.003) sit at `reach` 0.090 against our 0.158.

**`mse` stays 0 for us.** The original claim in section 4.1 was right about the conclusion and
wrong about the reason; this correction gets the reason right. Ibrahim Mansour's +0.0728 comes
from a genuinely better cosine (inverting his raw `mse` of 0.919 gives cos ~0.20 against our
~0.13), not from an emission trick.

**But the sweep that killed this found something bigger on the way** — see experiment 034.
