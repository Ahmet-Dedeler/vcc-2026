# 183 — unique-component boost. **NO-GO** (PLOLO +0.00104 raw `pds`, 3× short of the bar; composed net −0.0026 to +0.0017 against a +0.0019 floor)

**Date:** 2026-09-03. Measurement only. No `cell-eval2`, no build, no submission, nothing
written outside this directory. Peak RSS **7.20 GiB** (`plolo183.py`), 2.13 GiB (`h1test183.py`).

## The lever

For each already-whitened predicted-effect row `h`:

```
c  = mean_p h_p                 # centroid of the panel's predicted effects
u  = h - alpha * proj_c(h)      # strip a fraction of the shared DIRECTION
h' = u * (||h|| / ||u||)        # preserve the row's L2 exactly
```

Opened by audit 181's item 3: *"reopen `pds` with one lever that is NOT an un-whitening move —
170's criterion, norm-preserving and raising `X` without `A`, is unexplored outside sign repair."*
This is the smallest operator satisfying the first half of that criterion. `alpha ∈ {0, 0.25,
0.50, 0.75, 1.00}`, `alpha = 0` the known-answer identity.

## Verdict in one line

**It does exactly what it was designed to do and the effect is a fifth of what it needs to be.**
Mean pairwise row cosine falls `+0.00215 → −0.00313`, centroid energy `0.699% → 0.032%`, row
norms preserved to `4.4e-16` — and PLOLO `pds` rises by **`+0.00104` raw** (`+0.00230` scaled,
`+0.00039` on overall), against the `+0.003` **raw** bar this experiment was opened with and at
**1.56×** PLOLO's own leave-one-arm-out sigma (`0.001494` scaled = `0.000675` raw). **It also fails the second half of 170's criterion: it does not raise `X`, it lowers
it** on three of the four surfaces with truth, and it lowers sign accuracy `q` by `−0.0016`.

## 1. Known-answer gates (rule 10's corollary) — all three exact

| # | gate | result |
| --- | --- | --- |
| 1 | `alpha = 0` reproduces 162's own `a23` rows, 4 folds × 3 challenge frames × 3 placements | worst \|Δ`pds`\| = **`0.000e+00`** |
| 2 | `a23 − a15` in the challenge frame = 162's published mass-fix price `+0.004014` scaled | **`+0.004014`** (raw `+0.001815`) |
| 3 | `‖h'‖ / ‖h‖ − 1` over every transform call | **`4.441e-16`** |
| 4 | h1-test `alpha = 0` reproduces 168d/173's base `pds` `0.86134` | `0.861345`, \|d\| **`4.54e-06`** |

Gate 2 is the non-vacuous one: it reproduces an *independently published* PLOLO number to six
decimals, so the folds, the operators, the CD4 solve and the renderer are all the ones 162 used.

## 2. PLOLO, challenge frame — the instrument rule 13 names

Four donor folds over the real 300-target panel, rendered into contexts A/B/C, scored against
the held library's own deep pseudobulk. Deltas vs `alpha = 0`, placement `pre` (§4).

| alpha | k562_gwps | orion_hct116 | orion_hek293t | kolf2 | **MEAN raw** | sd | **MEAN scaled** |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.25 | +0.000541 | +0.000445 | +0.000255 | +0.000162 | **+0.000351** | 0.000150 | +0.000776 |
| 0.50 | +0.000915 | +0.000887 | +0.000395 | +0.000036 | **+0.000558** | 0.000366 | +0.001235 |
| 0.75 | +0.001424 | +0.001237 | +0.000654 | +0.000121 | **+0.000859** | 0.000512 | +0.001900 |
| **1.00** | +0.001845 | +0.001647 | +0.000697 | **−0.000027** | **+0.001040** | 0.000754 | **+0.002301** |

Monotone, 3 of 4 folds positive, but the whole response barely clears **162's own LOO sigma**
(`0.001494` scaled): `alpha ≤ 0.50` sits inside it, and the best point, `alpha = 1.00`, is
**1.56σ** (`+0.000388` overall against `±0.000249`). And **`kolf2` reads zero** — that is the only leak-free fold that uses
the *shipped* `union_pruned` operator and the *shipped* three-source donor mixture, i.e. the
closest available analogue of the real build. Its centroid energy is the same as the others'
(0.73%), so this is not "less fuel"; its base cosine is `0.0076` against `0.022–0.031`, so it is
the fold whose truth is furthest from the donor consensus, and there the lever does nothing.

### The geometry actually moved (4-fold × 3-frame mean, `pre`)

| alpha | `pds` | Δ`pds` | cos-to-truth | Δcos % | **mean pairwise cos** | **centroid energy** | median ‖h‖ | sse\* |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.742112 | — | 0.02188 | — | **+0.00215** | **0.00699** | 3.73002 | 0.999537 |
| 0.25 | 0.742463 | +0.000351 | 0.02187 | −0.05 | −0.00015 | 0.00411 | 3.73000 | 0.999537 |
| 0.50 | 0.742670 | +0.000558 | 0.02186 | −0.13 | −0.00180 | 0.00203 | 3.73000 | 0.999538 |
| 0.75 | 0.742971 | +0.000859 | 0.02183 | −0.24 | −0.00280 | 0.00077 | 3.72998 | 0.999539 |
| **1.00** | 0.743152 | **+0.001040** | 0.02180 | **−0.39** | **−0.00313** | **0.00032** | 3.73001 | 0.999541 |

Median row norm is constant to `1e-5` — `A` and `B` do not move, so `min(B, claim)`, `rho` and
the whole `mse` cliff position (170 §4) are untouched. That part of 170's criterion is met
exactly. **`sse*` rises** (worse) at every `alpha`, and **cos-to-truth falls monotonically**.

## 3. h1-test, delivered space (173's / 168d's exact pipeline)

85 targets, 10,783 gated genes, 17,995 scored. Donors `k562_gwps` + `orion_hek293t`;
`data/donors/h1_2025` is 85/85 oracle leakage and never enters.

| alpha | `pds` post | `pds` re | `pds` pre | Δ`pds` pre | Δcos % | mean pairwise cos | centroid energy |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.861345 | 0.861345 | 0.861345 | — | — | +0.000542 | 0.01743 |
| 0.25 | 0.861625 | 0.861485 | 0.861625 | +0.000280 | −0.17 | −0.004460 | 0.01028 |
| 0.50 | 0.863305 | 0.862605 | 0.862605 | +0.001261 | −0.42 | −0.008067 | 0.00506 |
| 0.75 | 0.864286 | 0.864006 | 0.864286 | +0.002941 | −0.76 | −0.010241 | 0.00191 |
| **1.00** | 0.866667 | 0.864566 | 0.864426 | **+0.003081** | **−1.17** | −0.010959 | 0.00095 |

h1-test *just* clears the `+0.003` bar; **PLOLO misses it by 3×**. Rule 13 says PLOLO is the
priceable `pds` channel and h1-test is not, and audit 181 says the h1-test/board `pds`
discrepancy is a **frame** effect — read `pds` in the challenge frame. So the `+0.00104` is the
number, and `+0.00308` is the smaller-decoy-pool (85 vs ~290) reading of the same move.

## 4. The renderer does NOT eat this transform — and that answers 156 §6

156 measured that a pre-render centroid move is silently undone (`centroid` alone, pre-render,
`+0.00000` on HCT116), and explained Duy's published board failure with it. Three placements
were run here: `post` = transform the delivered object; `re` = push that back through the
renderer; `pre` = transform the effect matrix and render once (**what production would do**).

| alpha | post | re | **pre** | \|re−post\| | \|pre−post\| |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.25 | +0.000370 | +0.000333 | +0.000351 | 0.000037 | 0.000019 |
| 0.50 | +0.000528 | +0.000540 | +0.000558 | 0.000012 | 0.000031 |
| 0.75 | +0.000759 | +0.000798 | +0.000859 | 0.000039 | 0.000100 |
| 1.00 | +0.000758 | +0.000972 | +0.001040 | 0.000214 | 0.000282 |

Row-wise cosine between the shipped `pre` object and the pure `post` object is `0.99988` at
`alpha = 1`. **With a23's `minnorm` mass fix in place, `deliver` is an identity to `1e-12` and
geometry applied before the renderer survives it.** 156 held this over as "the honest
single-variable follow-up if a23 confirms"; a23 confirmed, and this is that measurement. It is
the one durable positive result of this experiment and it applies to *any* future
delivered-space geometry arm: **you may now design in effect space.**

## 5. `mse` — the arm fails 170's actual criterion

`mse − 1 = (s²A − 2sX + Kc − min(s²B, claim)) / sumD`. `A` is preserved by construction
(`rA = 1.0002`), so `X` is the only term that can move. `X` needs withheld truth, so it is
transferred as a ratio on the two held Orion lines (169's `xratio183.py` route; 168d: report
both, never average).

| alpha | rX hct116 | rX hek293t | rA hct116 | rA hek293t |
| ---: | ---: | ---: | ---: | ---: |
| 0.50 | 0.958005 | 1.036744 | 1.000102 | 1.000114 |
| **1.00** | **0.912549** | **1.070642** | 1.000197 | 1.000221 |

**The two surfaces disagree in sign.** Worst-surface `alpha = 1` gives raw `mse` `+0.01356`
(`−0.00237` overall); best-surface gives `−0.00181` overall. Three of the four instruments with
truth (PLOLO's 4-fold cos and sse\*, h1-test's cos, Orion-HCT116's `X`) say alignment **falls**;
only Orion-HEK293T says it rises.

**So the lever satisfies "norm-preserving" and violates "raise `X` not `A`".** 170's `+0.0055`
prize was for a transform that converts the `+16%` of alignment gamma can buy *without* the
amplitude. This one buys `pds` by spending alignment, which is the same currency gamma and sign
repair spend — it is simply a much better exchange rate (§7).

## 6. `nmae` / `reach` / `fid` / `jac` — pinned by `q`, not laundered

Not priceable as levels (rule 9: `n_real`; rule 13: `nmae` 3/6 signs, `reach` 2/6). But the
2026-09-03 02:15 landscape note pins the **mechanism**: at fixed `k`, all three sign-bound
members move only through the DE-pool sign accuracy `q`. Both `q` and `k` are directly
measurable on h1-test's 65,044-pair reference DE pool with no `n_real` normalisation, so they
were measured rather than declared unmeasurable:

| alpha | `q` | Δ`q` | mean \|p\| on pool | r\|p\| | pool energy share |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.00 | 0.53745 | — | 0.050618 | 1.000000 | 0.083507 |
| 0.50 | 0.53708 | −0.00037 | 0.050614 | 0.999909 | 0.083436 |
| **1.00** | **0.53588** | **−0.00157** | 0.050661 | 1.000849 | 0.083310 |

**`k` does not move** (pool magnitude ratio `1.0008`, pool energy share `0.0835 → 0.0833`) — so
the 2026-09-03 design constraint ("must not push DE-pool magnitudes above their true values") is
satisfied, and the arm is clean on that axis. **`q` falls.** With `dnmae/dq = |k−1| − (k+1) =
−1.5` at our `k ≈ 0.75`, `Δq = −0.00157` gives raw `nmae` `+0.00236` = `−0.00065` overall.
`reach` and `fid` read the same `q` and inherit the same sign. Rule 13's `−0.0146` `reach` toll
is **not** applied: this is a profile-class arm and the 181 amendment restricts that toll to
emission-law arms (169 measured `+0.0006` for profile arms).

## 7. All six members, and the composed slot gate (rules 12 + 13)

Shape-only arm ⇒ floor `+0.0019`, no emission toll. `decide183.py`:

| alpha | `pds`→ovr | `mse`→ovr (worst) | `nmae`→ovr | `reach`/`fid`/`jac` | **NET worst** | NET best-`mse` | verdict |
| ---: | ---: | ---: | ---: | :---: | ---: | ---: | :---: |
| 0.25 | +0.000131 | −0.000557 | +0.000050 | ≤ 0, unpriced | **−0.000377** | +0.000688 | VETO |
| 0.50 | +0.000208 | −0.001140 | −0.000154 | ≤ 0, unpriced | **−0.001086** | +0.001051 | VETO |
| 0.75 | +0.000320 | −0.001745 | −0.000438 | ≤ 0, unpriced | **−0.001863** | +0.001350 | VETO |
| **1.00** | **+0.000388** | −0.002373 | −0.000654 | ≤ 0, unpriced | **−0.002640** | +0.001651 | VETO |

`pds` leg = PLOLO scaled × 1.011 / 6. `mse` leg = `−2·X_c·(rX−1)/sumD` at 170 §4's
board-anchored `X_c = 0.017766`, `sumD = 0.229141`, converted at 170's realised factor 0.175.
`nmae` leg = measured Δ`q` through the landscape formula, span 0.60.

**The largest `pds` leg this arm can produce is `+0.000388` overall — 1.56σ of PLOLO's own LOO
noise — and even the branch that believes the one surface where `X` rises tops out at
`+0.00165`, still under the `+0.0019` floor.** Every `alpha` is VETO. **NO-GO. No slot.**

## 8. The mechanism, and why this closes the family rather than just this arm

**The exchange rate is excellent and the fuel tank is empty.**

Define `η = Δpds / (relative Δ cos-to-truth)` — raw `pds` bought per unit of alignment given up.

| lever | surface | η |
| --- | --- | ---: |
| gamma (097 §9h) | h1-test | −0.068 |
| sign repair (173) | h1-test | −0.054 |
| **unique boost** | **h1-test** | **−0.263** |
| **unique boost** | **PLOLO, challenge frame** | **−0.269** |

**The boost is ~4–5× more `pds`-efficient per unit of cosine than any un-whitening move, and the
two independent instruments agree on that number to 2%.** It is a strictly better direction on
the same `pds`↔cos frontier — which is the substantive answer to audit 181's item 3: the
frontier's *slope* is not fixed, and un-whitening moves were sitting on a bad part of it.

It still fails, because the quantity available to trade is tiny. The panel's shared energy after
whitening is **0.699%** and its mean pairwise row cosine is **+0.00215**, against the **+0.0747**
Heisenberg's rank-1 deflation was fixing and the `0.7573 → 0.7961` AgingLab reports. Removing
100% of 0.699% is worth `+0.00104` raw `pds`. **Whitening already spent this lever**, which is
exactly what 156 §6 concluded from the pre-render `centroid` arm (`+0.00000` / `+0.00098`) — and
this experiment shows that conclusion was not an artefact of the renderer eating the transform
(§4), because with the mass fix the renderer does not eat it and the answer is the same.

**Consequence for the reopened `pds` question.** Audit 181 was right that the closure was WEAK
and right to demand a non-un-whitening lever. The lever exists, it behaves as designed, and it
is worth `+0.001` because **the shared-direction reservoir is 35× smaller than the one the
competitors' anti-share numbers come from.** Any future `pds` arm has to find its gain somewhere
other than the panel's own shared structure. 170's `+0.0055` prize — a transform that raises `X`
without raising `A` — is still unclaimed, and this experiment narrows what it can look like:
it must add alignment, not merely redistribute the existing direction, because every
redistribution of the existing direction costs `X` and `q`.

## 9. What was NOT measured, and why

* **`reach` / `fid` / `jac` as levels.** Rule 9 (`n_real`) and rule 13 (2/6 and 3/6 sign hits).
  Their *direction* is pinned by the measured Δ`q` (§6); no emitted-count run was made, so no
  level is quoted. Nothing would change the verdict: they are ≤ 0.
* **`X` on the challenge panel.** It cannot be measured, only solved from a board point (170
  §4). The Orion pair supplies the ratio only, reported per surface and never averaged (168d).
* **Anything outside this one transform.** No emission-law change, no alpha, no donor work, no
  operator change. 182 owns pretrained State.
* **No `cell-eval2`, no build, no submission, no write outside `experiments/183-unique-boost/`.**
  Production is untouched. 162's and 169's files were imported read-only; their cached folds
  and artifacts were reused, none rebuilt.
* **`alpha > 1` (over-correction) was not swept.** The mandate was `{0 … 1}`, the response is
  concave in `alpha` and already flattening between 0.75 and 1.00 on three of four folds, and
  the `mse`/`q` costs are still growing linearly there. It cannot reach `+0.0019`.

## 10. Files

| file | what it does |
| --- | --- |
| `uboost.py` | the transform, the norm assertion, and the shared-energy / hubness readouts |
| `plolo183.py` | 162's four PLOLO folds, challenge frame, three placements, five alphas |
| `h1test183.py` | the h1-test delivered-space surface (173's / 168d's exact pipeline) |
| `xratio183.py` | 096's `X` term on the two held Orion lines (169's route) |
| `deprice183.py` | `q` and `k` on h1-test's 65,044-pair reference DE pool |
| `decide183.py` | rule 13's composed slot gate over all six members |
| `work/plolo183.{json,log}` · `work/xratio183.{json,log}` · `work/decide183.{json,log}` · `work/h1test183.{npz,log,time}` · `work/deprice183.log` | outputs |
