> **INTERIM — the `pds`/`cos`/`nmae` half is measured on two Orion surfaces. `reach`, `fid`,
> `jac` and the exact `mse` are NOT yet priced, and `nmae` has only one surface. Do not ship on
> this document. Experiment 169 completes it.**

# 168 — magnitude sharpening, re-priced with the instrument experiment 030 lacked

**Date:** 2026-09-02 · two held Orion lines, delivered space (gate → production
`mass_consistent_effect` → render), deterministic, no emitter, no `cell-eval2`.

## Why this was reopened

Experiment 030 (2026-08-21) found that a one-parameter power transform of the predicted per-gene
magnitudes, `P → sign(P)·|P|^q` renormalised to the original median row norm, is a **large**
`pds` lever — macro PDS 0.8017 → 0.8415 at q = 2.5, roughly additive with whitening (whitening
alone +0.080, sharpening alone +0.072, both +0.120). **It was never shipped, and its own
write-up says exactly why:**

> "Not shipped, and not resolvable with the tools that produced it. … The `pds` gain is exact and
> large; the `nmae` cost is estimated by an instrument known to be miscalibrated on the axis that
> decides it."

The deciding axis was `nmae`, estimated through experiment 018's proxy — which put the
`nmae`-optimal alpha at ~0.75 against a count-level measurement of 1.08. Three things have since
closed that gap: `debench.Surface.nmae()` computes `nmae` **exactly** (153, 2026-08-31); 156
supplies the renderer, so the scored object is the delivered one; and `mse` was not known to be a
live member until 096 (2026-08-28), a week after 030 was parked.

## Result — the `pds` gain reproduces, and the `nmae` cliff is sharp and late

| q | worst-surface Δpds | Δpds hct116 | Δpds hek293t | Δcos hct | sse* hct | sse* hek | **Δnmae exact (hct)** |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.00 | — | — | — | — | 0.999326 | 0.999803 | — |
| **1.15** | **+0.01145** | +0.01145 | +0.01359 | +11.7% | **0.999148** | **0.999719** | **+0.0036** |
| 1.25 | +0.01799 | +0.01799 | +0.02198 | +19.9% | 0.999003 | 0.999649 | +0.0327 |
| 1.50 | +0.02944 | +0.02944 | +0.03950 | +37.7% | 0.998617 | 0.999447 | +0.2094 |
| 1.75 | +0.03597 | +0.03597 | +0.05419 | +49.9% | 0.998269 | 0.999255 | +0.4004 |
| 2.00 | +0.03985 | +0.03985 | +0.06529 | +58.2% | 0.998017 | 0.999119 | +0.5572 |
| 2.50 | +0.04251 | +0.04251 | +0.08310 | +65.3% | 0.997819 | 0.999990 | +0.9024 |

**`+0.01145` worst-surface raw `pds` at q = 1.15 is the largest local `pds` gain this project has
measured**, it is the same sign and comparable size on two independent cell lines, and `sse*`
(the `mse` proxy at its own best amplitude) improves on both. Sign accuracy on the reference's DE
pool also rises monotonically (0.5417 → 0.5578) and head purity is flat (0.6844 → 0.6785).

**The `nmae` cost is negligible until q ≈ 1.2 and then explodes** — 9x for a 0.10 change in q
between 1.15 and 1.25. That cliff is the whole reason 030 could not decide this, and it is why
the exact instrument mattered: a proxy that mis-sites the cliff by 0.1 in q inverts the verdict.

## Why sharpening should be the RIGHT version of a13's lever

a13 shrank per-gene by **donor agreement** and `nmae` collapsed to raw 1.0185 — past its own
baseline, worse than predicting nothing — while delivering `score_mse` **+0.0548**, the largest
single-member gain in this project's history. Its post-mortem named the mechanism: `nmae` scores
fold-change error over **the reference's** DE genes, and donor agreement is uncorrelated with
which genes those are, so the shrinkage landed squarely on the scored set.

Measured here, AUC for predicting reference-DE membership (HCT116, 127 targets with a ≥10-gene
gate):

| key | mean AUC | median | frac > 0.5 |
| --- | ---: | ---: | ---: |
| **\|predicted effect\|** | **0.5966** | 0.5975 | **0.976** |
| donor agreement (a13's key) | **0.4980** | 0.4919 | 0.425 |
| sqrt(\|eff\| · agreement) | 0.5399 | 0.5417 | 0.874 |

So a13 shrank on a key carrying **no information** about the scored set, and sharpening
concentrates magnitude on a key that does. Experiment 154 independently measured
top-decile-by-\|effect\| sign purity at 0.675 against an overall 0.531, so the magnitude predicts
*correctness* as well as membership — which is the quantity both `mse` and `nmae` are paid in.

A monotone power transform **preserves the \|effect\| ordering exactly**, so `keep_top` and every
\|effect\|-keyed policy are untouched by construction.

## The cliff is NOT simple overshoot — my first explanation was wrong

I predicted the knee would be where our \|lfc\| starts exceeding the reference's, since the docs
record `|h| ≤ |r|` on 99.3% of gated genes at q = 1. Measured on the reference's own DE gate:

| q | frac \|h\|>\|r\| | median \|h\|/\|r\| | q90 \|h\|/\|r\| | mean \|h\| | mean \|r\| |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1.00 | 1.12% | 0.1266 | 0.4319 | 0.0718 | 0.4433 |
| 1.15 | 1.79% | 0.1120 | 0.4497 | 0.0755 | 0.4433 |
| 1.25 | 2.14% | 0.1018 | 0.4548 | 0.0864 | 0.4433 |
| 2.00 | 3.40% | 0.0428 | 0.3931 | 0.2430 | 0.4433 |

The overshoot fraction rises only 1.12% → 3.40%, the **median** ratio *falls* (most genes get
smaller), and mean \|h\| rises 3.4x **toward** mean \|r\| — i.e. by the crude mean-magnitude
criterion, sharpening makes us *better* calibrated while `nmae` gets much worse. So the damage is
concentrated in a heavy tail rather than being a broad miscalibration, and the simple
"overshoot" story is refuted. **The decomposition — Σ|h−r| split by sign-correct/sign-wrong and
by overshoot/undershoot — is left to experiment 169 and must be done before this ships.**

## Status and what is NOT measured

- **`reach`, `fid`, `jac`: unmeasured.** Sharpening changes the emitted magnitudes and therefore
  the emitted p-values, so the composed slot gate (`AGENTS.md` 11.1 rule 13) charges the
  empirical emitted-profile **`reach` toll of −0.0146** until it is priced directly on `h1-test`.
  At q = 1.15 the priced `pds` gain is ≈ `+0.01145/0.45208/6 = +0.00422` overall, so the toll is
  more than half of it. **This arm is not shippable until `reach` is measured.**
- **`nmae` has ONE surface.** HEK293T has no cached reference DE table, so its column is `nan`.
  `h1-test` does have one (157) and is the only board-validated surface.
- **`mse` is proxied by `sse*`, not priced.** 162 showed `mse` is not priceable by a cross-line
  validator at all; it needs 096's `claim` route, and 164 established a15 sits at 99.7% of the
  dispersion cap, so the interaction with `claim` must be computed rather than assumed.
- Absolute values on the Orion pair are not board-comparable (those harnesses carry 26.4% /
  44.1% of their truth's effect energy in a target-independent vector; `AGENTS.md` §6.4). Only
  the relative moves are claimed, and both surfaces agree in sign at every q.

## Files

`sharpen.py` · `work/sharpen.json`.

---

# 168b — the mechanism, and it yields a better arm: HEAD-TARGETED sharpening

Experiment 030 could not explain its own gain ("not distribution matching, not selection"). It
is now identified, and the explanation predicts a strictly better transform.

## The mechanism: our directional signal is almost entirely in the top few percent of |D|

Binning every gated `(target, gene)` pair by the pre-render designed magnitude `|D|`, with the
cross-target mean removed from **both** sides (the Orion harnesses carry 26.4% / 44.1% of their
truth energy in a target-independent vector), and measuring each bin's share of the aggregate
alignment `X = Σ D·T` — which is exactly what the cosine and `mse` are made of:

| | top 5% of \|D\| | top 20% | local slope of E[t\|D] at the top |
| --- | ---: | ---: | ---: |
| HCT116 | **69.8%** of X | 87.0% | +3.11 |
| HEK293T | **73.1%** of X | 85.9% | +3.65 |

The middle bins carry 0–3% of X each and their local slopes are noise (−8.5 to +7.8). So
`E[t|D]` is **strongly convex in log-log**, not a power law: a global fit gives
`|D|^0.787` (HCT116, r = 0.82) and `|D|^0.972` (HEK293T, r = 0.94) — i.e. **q ≈ 1, which is what
we already ship** — while the top ventile's local slope is ≈ 3.

**That resolves the apparent paradox.** The conditional mean is the optimal estimator for
*squared error*; `pds` is a *rank* metric, and a global power law is a crude way to buy the
top-end amplification that discrimination wants. It also explains why 030 found smooth
sharpening beats hard top-K selection: hard truncation estimates each similarity from only K
coordinates and adds variance, which is the standard reason soft weighting beats thresholding in
a matched filter.

**A methodological note, recorded because it nearly fooled me.** My first pass binned the
*rendered* effect without removing the shared mean. It put 6 of 14 bins at mean `|pred|` =
0.00000 (render dust) and reported sign accuracies of **0.9714 and 0.8636** in the smallest bins
— my own mass-fix correction interacting with the HEK293T blob — yielding "derived exponents" of
0.043 and 1.774 from 12/14 and 3/14 bins. Both numbers were garbage and neither was reported as a
result.

## The better arm the mechanism predicts

If the signal is in the head, amplify the head and leave the bulk alone. `top f% × g` multiplies
each row's largest `f%` of gated `|D|` by `g`, then renormalises to the original median row norm.
Exact `nmae` (HCT116), two-surface worst-case `pds`:

| arm | worst Δpds | Δnmae | sse* hct | pds→overall | nmae→overall | **net (pds+nmae)** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| GLOBAL q=1.15 | +0.01145 | +0.0036 | 0.999148 | +0.00422 | −0.00100 | +0.00322 |
| GLOBAL q=1.50 | +0.02944 | **+0.2094** | 0.998617 | +0.01085 | −0.05817 | −0.04732 |
| **top 2% ×2** | **+0.02602** | **+0.0116** | **0.998961** | **+0.00959** | **−0.00323** | **+0.00636** |
| top 3% ×2 | +0.02238 | +0.0089 | 0.999006 | +0.00825 | −0.00246 | +0.00579 |
| top 3% ×3 | +0.03065 | +0.0208 | 0.998899 | +0.01130 | −0.00577 | +0.00553 |
| top 2% ×3 | +0.03452 | +0.0471 | 0.998810 | +0.01273 | −0.01310 | −0.00037 |
| top 1% ×4 | +0.03100 | +0.1442 | 0.998563 | +0.01143 | −0.04005 | −0.02862 |
| top 5% ×2 | +0.01660 | +0.0058 | 0.999091 | +0.00612 | −0.00160 | +0.00452 |

**`top 2% ×2` doubles the best global arm's net.** It buys 2.3x global q=1.15's `pds` for 3.2x its
`nmae`, which is a better trade because the `pds` span is 0.452 and the `nmae` span is ~0.60.
`sign accuracy` on the reference DE pool also rises (+0.0022). Both surfaces agree in sign at
every setting, and the ordering of arms is the same on both.

**Adding `mse`:** `sse*` 0.999326 → 0.998961 is a cosine gain of **+24.2%**. Scaling a23's
*realised* board figure (it delivered `mse` +0.00153 for a +13.6% harness cosine gain) gives
≈ **+0.0027**. Net with the full rule-13 reach toll charged:

```
+0.00959 (pds) − 0.00323 (nmae) + 0.0027 (mse) − 0.00243 (reach toll) ≈ +0.0066 overall
```

3.5x rule 13's +0.0019 shape-only threshold, and it would be our best entry ever
(0.17409 → ~0.1807). If `reach` proves neutral it is ≈ +0.0090.

## The risk that decides it, stated plainly

**`reach` reads the head, and this arm amplifies the head.** a22 is the precedent and the
warning: it aimed the *cell-channel* lever at the top decile and lost `reach` −0.01066, with the
recorded lesson that "the cell-channel lever cannot be aimed at the head without disturbing the
head's order." This arm changes head *magnitudes* rather than head *emission*, and a monotone
change preserves the `|effect|` ordering exactly — but `reach` reads the emitted **p-value**
order, which `docs/CURRENT_STATE.md` records as correlating with `|predicted effect|` at only
**+0.071**. So whether the p-order moves is an open empirical question with a real chance of
going either way, and head purity already dips slightly (−0.0076 at top 2% ×2).

**Do not ship this until `reach` is measured on `h1-test`.** Experiment 169 owns that.

Also unmeasured: `fid`, `jac`, a second surface for `nmae`, and the exact `mse` through 096's
`claim` route (164 established a15 sits at 99.7% of the dispersion cap, so `claim` must be
recomputed rather than assumed).

---

# 168c — the data-derived transform, its control, and a CORRECTION to 168b's mechanism claim

## The transform

168b's arms (`|D|^q`, `top f% x g`) are parametric guesses at the measured `E[t | |D|]` curve, and
their parameters are fitted — which is what will not transfer to D/E/F. So estimate the curve
directly: equal-count bins of `|D|`, `m(|D|) = E[t·sign(D) | |D|]`, clipped at zero and made
isotonic, applied as a per-coordinate multiplier `w = (m/|D|)^tau`. `tau = 1` is the conditional
mean; `tau = 0` is the identity.

**Fitted on one Orion line, scored on the other, both directions** — no truth from the scored
surface enters the transform, and this is also how it would deploy (fit on donor↔donor pairs,
apply to A/B/C).

| tau | worst-surface Δpds | Δnmae (hct) | net pds+nmae → overall |
| ---: | ---: | ---: | ---: |
| 0.25 | +0.00995 | **−0.0007** | +0.00387 |
| **0.50** | **+0.01845** | **+0.0009** | **+0.00654** |
| **0.75** | **+0.02388** | +0.0073 | **+0.00679** |
| 1.00 | +0.02550 | +0.0110 | +0.00635 |
| 1.50 | +0.02805 | +0.0171 | +0.00558 |
| 3.00 | −0.01058 | +0.0061 | collapses |

`tau = 0.5` buys **+0.0185 raw pds at zero `nmae` cost** and improves `sse*` on both surfaces
(0.999326 → 0.999136 and 0.999803 → 0.999723). That is a better trade than any 168b arm.

## The control, and it half-fails

Per AGENTS.md 11.1 rule 10's corollary, the same reweighting machinery was run with the
multiplier vector destroyed in four different ways:

| variant | Δpds hct116 | Δpds hek293t |
| --- | ---: | ---: |
| **FITTED curve** | **+0.01845** | **+0.01874** |
| REVERSED curve | **−0.12838** | **+0.01712** |
| SHUFFLED bins (same multipliers, wrong bins) | −0.02169 | −0.01434 |
| random monotone decreasing | −0.08530 | −0.05361 |
| random monotone increasing | +0.00405 | +0.00145 |
| flat (identity) | 0.00000 | 0.00000 |

**What passes:** the SHUFFLED control is negative on **both** surfaces. The gain therefore
requires the multipliers to be attached to the *correct* `|D|` bins and is not a generic
consequence of reweighting. Random monotone curves are 4–13x smaller or strongly negative.

**What fails:** on HEK293T the **REVERSED** curve is worth +0.01712 against the fitted +0.01874,
and its cosine gain is *larger* (+30.9% vs +17.8%). Two transforms in opposite directions cannot
both be "the better estimator", so at least one HEK293T reading is not measuring alignment with
target-specific structure. The likely cause is that surface's **44.1% target-independent blob**
(`AGENTS.md` §6.4; HCT116 is 26.4%): the curve is fitted with the cross-target mean removed from
both sides, but scoring uses the raw truth, which still contains it. On HCT116, where the blob is
smaller, direction-sensitivity is intact and decisive (fitted +0.018 vs reversed −0.128).

## CORRECTION to 168b

168b concluded that "the conditional mean is optimal for squared error while `pds` is a rank
metric, and a global power law is a crude way to buy the top-end amplification discrimination
wants." **That story is now over-claimed and I am withdrawing the strong form of it.** The
evidence against it is internal to this experiment:

- The fitted conditional mean is **sub-linear** — its multiplier *falls* from 0.51 to 0.098 across
  `|D|`, so it **compresses** the magnitude range. That matches the global power-law fits
  (q = 0.787 / 0.972).
- 168b's head arms **expand** the magnitude range, and also help (`top 2% x2`, +0.02602).
- Compression and expansion cannot both follow from a single "optimal estimator" argument.

What survives is narrower and still useful: (i) the alignment concentration is a solid,
two-surface measurement (top 5% of `|D|` carries 69.8% / 73.1% of `Σ D·T`); (ii) the fitted-curve
gain is real on HCT116 and survives the shuffled control on both; (iii) `pds` is a rank metric
and evidently admits gains from more than one reweighting direction, so **no reweighting arm in
this family should be justified by a mechanism story until it is read on a low-blob surface.**

## Consequence for shipping

**`h1-test` is now the arbiter for this whole family, not just for `reach`.** It carries a 3.6%
target-independent blob against the Orion pair's 26.4% / 44.1%, and it is the only surface
validated against a board `reach` point. Both the head arms of 168b and the isotonic arm here
must be re-read there before either is a candidate. Experiment 169 owns that surface; the
correction and the arms are in `experiments/169-sharpening-six-member/INBOX-from-168.md`.

Nothing in 168 is shippable on the evidence in 168/168b/168c.

---

# 168d — TERMINAL NO-GO. h1-test breaks the tie and kills the family.

`h1_tiebreak.py`. h1-test is the arbiter: **3.6% target-independent blob** against the Orion
pair's 26.4% / 44.1%, the only surface validated against a board `reach` point, and — uniquely —
scoreable through the **production** `union_pruned` operator without leakage, because it is not
an Orion line. Base effect recovered from the frozen 095 artifact
(`out-h1test/.../g0.65_q1.0_a1.08_spread0.25/pred.h5ad`) as `log1p(TS·p) − b_ctrl`, which is
mass-consistent by construction. The transform curve is **percentile-indexed** so it is scale-free
across surfaces shipping different alphas, and it is fitted on the Orion pair only.

## 1. The tie is broken in favour of "direction matters"

| variant, tau = 0.5 | h1-test Δpds | HCT116 | HEK293T |
| --- | ---: | ---: | ---: |
| **FITTED** | **+0.00280** | +0.01845 | +0.01874 |
| REVERSED | **−0.00378** | −0.12838 | +0.01712 |
| SHUFFLED | **−0.01737** | −0.02169 | −0.01434 |

h1-test agrees with HCT116 qualitatively — fitted positive, reversed negative, shuffled worst —
so the transform is **not** a generic reweighting artefact and HEK293T's indifference to
reversal is that surface's problem, not the method's.

## 2. And the magnitude collapses. Sized on h1-test, the family fails.

| arm | Δpds h1-test | Δpds Orion (worst) | Orion overstatement | pds → overall |
| --- | ---: | ---: | ---: | ---: |
| GLOBAL q=1.15 | +0.00112 | +0.01145 | **10.2x** | +0.00041 |
| GLOBAL q=1.50 | +0.00784 | +0.02944 | **3.8x** | +0.00289 |
| isotonic tau=0.5 | +0.00280 | +0.01845 | **6.6x** | +0.00103 |
| **top 2% x2** | **−0.00042** | +0.02602 | **SIGN FLIP** | −0.00015 |
| **top 2% x3** | **−0.00994** | +0.03452 | **SIGN FLIP** | −0.00367 |
| top 3% x2 | −0.00196 | +0.02238 | **SIGN FLIP** | −0.00072 |
| top 5% x2 | +0.00028 | +0.01660 | 59x | +0.00010 |

**Every arm fails rule 13's +0.0019 shape-only threshold**, and the −0.00243 reach toll sinks
them regardless. The one arm with a real `pds` gain on h1-test, GLOBAL q=1.50 (+0.00289 overall),
is the one whose exact `nmae` cost is **+0.2094 raw** — worth about −0.058 overall. There is no
operating point.

**Decision: TERMINAL NO-GO on magnitude sharpening in every form tested — global power,
head-targeted multiplication, and the data-derived isotonic transform.** Experiment 030's
+0.072 macro PDS does not survive rendering onto a low-blob surface. Do not reopen without a
materially different mechanism.

## 3. What I got wrong tonight, twice, and why it is the useful part

1. **168b's mechanism claim** ("the conditional mean is optimal for squared error, `pds` is a
   rank metric") was withdrawn in 168c after my own control showed the fitted curve *compresses*
   while the head arms *expand*, and both appeared to help.
2. **The head-targeted arms I recommended in 168b as "double the best global arm" have the
   WRONG SIGN on the validated surface.** `top 2% x2` reads +0.02602 on Orion and −0.00042 on
   h1-test. I had already flagged that the arm was not shippable until `reach` was measured; the
   number that killed it was `pds`, on a better surface, which I had not measured when I wrote
   the recommendation.

Both errors have the same root: **sizing an arm on the Orion pair.**

## 4. The instrument rule this establishes, measured four independent ways tonight

| evidence | overstatement |
| --- | --- |
| 157: repaired HCT116 vs h1-test on the tier-policy `reach` effect | **4.4x** |
| target-independent blob share: Orion 26.4% / 44.1% vs h1-test 3.6% | — |
| 166: 159's line-pair values are one realisation of a random draw | ±0.001–0.002 |
| **this experiment: Orion vs h1-test on sharpening `pds`** | **3.8–10x, and sign flips** |

**Never size an arm on the Orion pair. Use `h1-test`.** The Orion surfaces remain usable for
arm-to-arm *ordering* on mean-removed statistics, and their controls (SHUFFLED, REVERSED) did
their job here. But a magnitude read off them can be wrong by an order of magnitude and can
invert.

Note this does **not** retroactively invalidate a23: its Orion screen predicted +0.00035/+0.00278
raw `pds` and the board delivered +0.00075, i.e. Orion *under*-stated a small perturbation. The
failure mode is specific to **large reshapings**, which interact with the blob.
