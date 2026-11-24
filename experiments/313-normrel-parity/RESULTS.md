# 313 — Implementation parity: production's per-target weighting is NOT the object 308 measured

**Verdict.** *Production's per-target weighting **IS NOT** the object 308 measured* (max |diff| **2.700** on the
non-panel columns of the pre-whiten effect matrix, against a float32 noise floor of **3.39e-07** — a factor of
**7.9e+06**; median per-target row cosine between the two outputs **0.9857**, 15/20 targets below 0.99, worst
non-degenerate **0.9120**). *The difference is* **three independent things**: (i) the PerturbMatch correction is
applied at the **uniform** donor share `1/|have|` instead of the arm's own per-target share, and **outside** the
norm match; (ii) the weights themselves are computed on a **different vector** — production's norm includes the
panel columns (the target's own knockdown) and excludes the PerturbMatch blend, 308's does the opposite; and
(iii) — pre-existing, not caused by 312 — production applies PerturbMatch to **every official column PerturbMatch
observes**, while the PLOLO harness 308 rides on restricts it to `k562 ∩ PM ∩ non-panel`.

**The consequence that matters for the board reading: 308's arm is exactly shape-only (delivered row-norm ratio
1.000000 for every target, both axes). Production's, as shipped, is not** — per-target delivered norm
**0.9930 … 1.3407** on the full axis, **0.8062 … 2.1549** on the scored non-panel columns, whole-matrix norm
**+6.00 %**. `VCC_NORMREL_EXP=-1.5` was preregistered as *shape-only* with a predicted `nmae` leg of **−0.00095
(a benefit)**; the board returned `score_nmae` **−0.007705**. The shipped object is amplitude-changing, so an
`nmae` loss of the wrong sign is exactly what its shipped form predicts.

Everything below is measured on one shared synthetic fixture with both implementations' **real functions**.
Peak RSS **0.21 GiB**. No bank was loaded, no build was run, nothing outside this directory was written.

---

## 1. The two orderings, precisely

**308 (blend → weight → norm-match).** `sweep308.py:184-185` builds the arm as
`Z.normrel_weights(st["rows"], avail, wc)` followed by `Y.combine(st, W, norm_match=True)`. `st["rows"]` comes from
`lib263.stage` (`lib263.py:123`) via `lib263.prepared_sources` (`lib263.py:61-94`), and that function has **already
blended PerturbMatch into the K562 row**: `lib263.py:82-84` overwrites `eff[rows, shared]` with
`(1-PM_BLEND)·k562 + PM_BLEND·pm` on the `shared = E064.shared_k562_pm_columns(...)` support (`064/run.py:119-136`:
official columns **not in the panel**, present in **both** the K562 bank and the PM bank), for the targets both
carry. The cached uniform base `st["pre"]` is `E062.build_baseline` over those same blended rows, so **the mean
row, the per-target weights and the norm-match reference all see the blended K562 row**. `lib308.normrel_weights`
(`lib308.py:207-221`) then takes `nj = ‖row_j(t)‖` **over `wcols` only** — and the prepared rows have their panel
columns zeroed (`062/run.py:287`, `effects[:, panel_positions] = 0.0`) — divides by that donor's median over the
targets it carries, raises to −1.5 with a floor `max(K, 1e-9)`, and renormalises over the carrying donors.
`lib263.combine` (`lib263.py:160-197`) forms `Σ_j w_j r_j` and, with `norm_match=True`, rescales it to
`‖pre[t]‖` — so **the arm is exactly norm-preserving: shape-only.**

**Production (weight → norm-match → correct).** `build_v9_pre_whiten_effects`
(`scripts/make_submission_pm_cd4.py:1390-1415`) dispatches to `_normrel_pre_whiten_effects` (`MB:1418-1469`) when
`VCC_NORMREL_EXP` is set. That function takes `nv = ‖model.s[name][idx[t]][model.dcols[name]]‖` — the **raw** donor
row over **the donor's whole official-axis coverage, panel columns included**, and **without** the PerturbMatch
blend (`model.s[K562]` is the bank) — divides by the donor's median over the panel targets it carries (`MB:1429-
1436`), raises to the exponent, multiplies by the scalar donor weight, and falls back to the **uniform** weight
when `med` or `nv` is exactly zero (`MB:1452-1457`). It then norm-matches the weighted row to the **PerturbMatch-
free** uniform row `base` (`MB:1464-1467`). Only afterwards does `build_production_effect_matrix` (`MB:2841-2848`)
call `apply_perturbmatch_k562_replacement` (`MB:1472-1593`) with `weights=active_weights` — the **scalar, uniform**
production weights — so `k562_share = active_weights[K562] / total_weight` (`MB:1545`) is `1/|have|` regardless of
what the arm did, and the correction `α · k562_share · effect_scale · blend · (pm_std − incumbent)` (`MB:1559-1565`)
is added **on top of** the norm-matched row, on `pm_target_columns` = **every** official column PerturbMatch
observes (panel columns and columns K562-GWPS does not carry included; the "experiment-061 broad gene policy" of
`MB:1484-1487`). The coordinator's 312 patch (`experiments/312-normrel-arm/apply_patch.py`) touched only the
pre-whiten function, so nothing downstream learned about the per-target weights.

**Algebra.** With `r_j` the prepared donor rows, `c = (r_pm − r_k)|_C`, `b = 0.5`, `u_j = 1/|have|`, `w_j` the
normrel weights, `S_w = Σ w_j r_j`, `S_u = Σ u_j r_j`:

| | row emitted |
|---|---|
| 308 (ORD_A) | `λ_A · (S_w + w_k·b·c)`, `λ_A = ‖S_u + u_k·b·c‖ / ‖S_w + w_k·b·c‖` |
| production (ORD_B) | `λ_B · S_w + u_k·b·c`, `λ_B = ‖S_u‖ / ‖S_w‖` |

`ORD_A − ORD_B = (λ_A − λ_B)·S_w + (λ_A·w_k − u_k)·b·c`. Both terms vanish iff `b = 0` (no PerturbMatch) or
`w = u` (exponent 0). `‖ORD_A‖ = ‖incumbent‖` identically; `‖ORD_B‖` does not.

## 2. The fixture

`work/fixture313.npz`, sha256 `723ce13c0c3b2e63e6d0af80271c90a34f3f06ee2b3130fd0a4509d636de64df`, seed
**3130907** (`fixture313.py`). 20 targets × 60 official genes × 3 donors named exactly `k562_gwps`,
`orion_hct116`, `orion_hek293t` so `MB.donor_weights` accepts the model. Official columns 0–19 are the panel
(target *i*'s own gene is column *i*, carrying a −2.4 knockdown); 20–59 are the operator columns (`wcols`, 40).
Gene masks: k562 covers 0–49, hct116 all 60, hek293t 10–59; PerturbMatch observes 5–54 (50 columns, **15 of them
panel**; `263`'s `shared` support is 30 columns). Cases: **(a)** PM absent for `T003`/`T019`, present for 17 of 20;
**(b)** the three masks above; **(c)** `T014`/`T015` missing from one donor, `T016`/`T017` from two, `T018` missing
from k562 so the PM branch cannot fire; **(d)** `T005` a degenerate hct116 row (set to the mean of that donor's
other rows, so the centred row is ~1e-17 — the division guard) and `T006` a small-but-ordinary hek293t row at
0.05×; **(e)** `T007` a k562 row at 8× its own median.

Three layers are run so an ordering difference is never mistaken for a dtype difference: **L1** production's real
`build_v9_pre_whiten_effects` → `_normrel_pre_whiten_effects` → `apply_perturbmatch_k562_replacement` on a real
`MultiDonorModel`; **L2** the real `lib308.normrel_weights` + `lib263.combine(norm_match=True)`; **L3** a clean-room
float64 statement of each ordering, verified against L1/L2 (KA-4, KA-5). Constants read from the builder:
`ALPHA 1.0796`, `BETA 0.0`, `PM_BLEND 0.5`, `effect_scale 6.817`, `donor_weights` all 1.0.

## 3. Known-answer rows — 9/9 pass, before any conclusion

| row | claim | \|d\| | tol |
|---|---|---|---|
| KA-P1 | production's uniform pre-whiten matrix == `ALPHA` × the harness's uniform mean (non-panel) | 3.395e-07 | 1e-06 |
| KA-1a | **production at exponent 0 == production uniform, bitwise** | **0.000e+00** | 0 |
| KA-1b | 308 at exponent 0 == the cached uniform row | 5.960e-08 | 1e-06 |
| KA-1c | 308's exponent-0 weights are exactly uniform | 0.000e+00 | 1e-15 |
| KA-2a | the clean-room 308 weight rule == `lib308.normrel_weights` | 0.000e+00 | 1e-15 |
| KA-2b | **production's weights, recovered by least squares from its own output, == the production weight rule** | 1.050e-08 | 1e-05 |
| KA-3 | **PerturbMatch absent + identical weights ⇒ the two orderings coincide** (clean float64) | **0.000e+00** | 1e-12 |
| KA-4 | the clean-room ORD_B reproduces production's real output | 5.604e-07 | 1e-04 |
| KA-5 | the clean-room ORD_A reproduces 308's real output | 6.804e-08 | 1e-05 |

KA-1a is exact because production's `(nv/med)**0 * 1.0` returns literally `1.0` for every donor, so `cand` is
`base` bit for bit and `nb/nc` is exactly 1. KA-3 is the mandate's analytic question, answered **yes** and proved
at 0.0e+00: **with PerturbMatch absent and the same weights the two orderings are the same map.** KA-2b is the
required demonstration that production's weight rule is what I say it is — the weights were recovered from
production's own emitted rows, not asserted.

The measured **float32 noise floor** between the two substrates is **3.395e-07** at a median row norm of 3.589
(both implementations cast to float32, at different points). Every divergence below is reported against it.

## 4. Per-case divergence

Non-panel columns, pre-whiten effect units (median row norm 3.589). `ordering` isolates D1+D2 (identical weights,
identical substrate); `max|dW|` is the largest per-donor weight difference between the two weight rules;
`share` is the `k562_share` **recovered from production's own PerturbMatch output**, against the K562 weight each
rule actually assigns.

| target | case | nD | PM | max\|diff\| | rel. row-norm | row cos | ordering only | max\|dW\| | w_k 308 | w_k prod | k562_share used |
|---|---|---|---|---|---|---|---|---|---|---|---|
| T000 | (a) PM present, 3 donors | 3 | Y | 3.969e-01 | +4.29e-02 | 0.99266 | 2.296e-02 | 2.33e-02 | 0.3148 | 0.3071 | **0.3333** |
| T001 | (a) PM present, 3 donors | 3 | Y | 3.149e-01 | +9.21e-02 | 0.98167 | 1.058e-01 | 1.50e-01 | 0.2517 | 0.3854 | **0.3333** |
| T002 | (a) PM present, 3 donors | 3 | Y | 4.650e-01 | +2.57e-02 | 0.97326 | 3.348e-03 | 7.47e-02 | 0.3300 | 0.3508 | **0.3333** |
| T003 | **(a) PM ABSENT**, 3 donors | 3 | – | 1.411e-01 | +8.85e-02 | 0.99958 | **0.000e+00** | 2.08e-02 | 0.2852 | 0.2919 | – |
| T004 | (a) PM present, 3 donors | 3 | Y | 2.254e-01 | −1.08e-02 | 0.99035 | 1.457e-01 | 7.03e-02 | 0.4361 | 0.3658 | **0.3333** |
| T005 | **(d) near-zero donor row** | 3 | Y | 1.315e+00 | −1.76e-01 | **0.05225** | 9.202e-02 | 5.73e-14 | 0.0000 | 0.0000 | 0.3333 |
| T006 | (d) small row, 0.05× | 3 | Y | 3.071e-01 | +1.79e-01 | 0.96085 | 3.018e-01 | 5.83e-02 | 0.0890 | 0.1228 | 0.3333 |
| T007 | **(e) row 8× its own median** | 3 | Y | **2.700e+00** | **+1.16e+00** | 0.91204 | 8.023e-01 | 1.28e-01 | 0.0718 | 0.0288 | 0.3333 |
| T008 | (a) PM present, 3 donors | 3 | Y | 3.357e-01 | +7.70e-02 | 0.98860 | 3.633e-03 | 6.84e-02 | 0.3355 | 0.3320 | 0.3333 |
| T009 | (a) PM present, 3 donors | 3 | Y | 3.354e-01 | +6.99e-02 | 0.99021 | 7.448e-04 | 5.54e-02 | 0.3332 | 0.2778 | 0.3333 |
| T010 | (a) PM present, 3 donors | 3 | Y | 4.204e-01 | +1.30e-02 | 0.98728 | 6.910e-03 | 2.85e-02 | 0.3265 | 0.3550 | 0.3333 |
| T011 | (a) PM present, 3 donors | 3 | Y | 5.044e-01 | −5.46e-04 | 0.98598 | 3.072e-02 | 2.57e-02 | 0.3496 | 0.3714 | 0.3333 |
| T012 | (a) PM present, 3 donors | 3 | Y | 3.145e-01 | +2.20e-02 | 0.98503 | 7.639e-02 | 6.79e-02 | 0.3925 | 0.3961 | 0.3333 |
| T013 | (a) PM present, 3 donors | 3 | Y | 6.321e-01 | +9.07e-02 | 0.96700 | 8.194e-02 | 9.98e-02 | 0.2920 | 0.2675 | 0.3333 |
| T014 | **(c) missing from ONE donor** | 2 | Y | 4.414e-01 | −8.55e-04 | 0.98717 | 2.324e-02 | 2.85e-02 | 0.4864 | 0.5149 | **0.5000** |
| T015 | (c) missing from ONE donor | 2 | Y | 3.600e-01 | +3.01e-03 | 0.98572 | 9.912e-02 | 7.48e-02 | 0.4446 | 0.5194 | **0.5000** |
| T016 | **(c) missing from TWO donors** | 1 | Y | 4.964e-01 | +2.91e-02 | 0.97170 | **1.11e-16** | – | – | 1.0000 | 1.0000 |
| T017 | (c) missing from TWO donors | 1 | Y | 1.275e+00 | +6.51e-02 | 0.93886 | **2.22e-16** | – | – | 1.0000 | 1.0000 |
| T018 | (c) missing from k562, PM inert | 2 | – | 1.237e-01 | −1.25e-02 | 0.99717 | **0.000e+00** | 4.53e-02 | – | – | – |
| T019 | **(a) PM ABSENT**, 3 donors | 3 | – | 1.960e-01 | +8.74e-03 | 0.98569 | **0.000e+00** | 1.13e-01 | 0.3555 | 0.3810 | – |

**Do the per-target weights agree? No.** Over the 18 weighted targets, max |dW| **0.1502**, median |dW|
**0.02851** — i.e. the median donor weight moves by ~2.9 percentage points out of a nominal 0.333 (≈ 9 % relative)
and the worst by 15 points (45 % relative). The two rules disagree even where PerturbMatch is absent (T003, T019).

**`k562_share`, demonstrated.** Recovered from production's own PerturbMatch output on every modified target:
it is `1/|have|` **without exception** — 0.3333 on the 13 three-donor targets, 0.5000 on the two-donor ones,
1.0000 on the single-donor ones — while the arm actually built those rows with K562 weights of 0.0288 … 0.4361.
On `T007` production applies the PerturbMatch correction at 0.3333 to a row it had just built with a K562 weight
of 0.0288: **an 11.6× over-application**. On `T004` it applies 0.3333 where the arm used 0.4361 (0.76×).

## 5. Attribution — and the part that predates 312

| term | what it is | size (max, non-panel, pre-whiten units) |
|---|---|---|
| **total** | production as shipped vs 308 as scored | **2.7000** (7.95e+06 × the float32 floor) |
| D1+D2 | the PerturbMatch **ordering** alone (identical weights, identical substrate) | 0.8662 |
| D3+D5 | the **weight rule** alone (identical ordering, identical substrate) | 1.5584 |
| D7 | the PerturbMatch **column policy** alone (broad official vs `k562∩PM∩non-panel`) | 1.2745 |
| D5 | panel columns inside production's norm | max \|dW\| **0.1502** |
| D3 | the PerturbMatch blend inside 308's norm | max \|dW\| **0.1058** |

**The pre-312 baseline is the load-bearing number.** With *uniform* weights on both sides — i.e. `cd4svd5` itself,
no arm at all — production's incumbent and the 308 harness's cached base already differ by **1.2745** on the
non-panel columns. Replace production's PM column policy with 263's `k562∩PM∩non-panel` and that collapses to
**6.43e-08**, *below the float32 floor*. So **D7 is a pre-existing PLOLO-harness-vs-production discrepancy that
is present in the incumbent and is not caused by 312**: production writes half a PerturbMatch row into official
columns the K562 bank does not carry (where the incumbent value is 0), and the harness does not. Everything else
about the two implementations agrees exactly at uniform weights. `T016`/`T017` in the table are the clean
demonstration — one donor, weight 1.0 on both sides, ordering term 1e-16, and yet 0.50 / 1.27 of divergence,
all of it D7.

**Scale caveat on D5, measured (`scale313.py`, `work/scale313.json`).** The mandated fixture's panel is 33 % of the
axis; production's is 300/18,533 = 1.6 %. Sweeping the same construction over axis width:

| genes | panel % of columns | panel % of row energy | max \|dW\| | p90 \|dW\| | median \|dW\| |
|---|---|---|---|---|---|
| 60 (the fixture) | 33.3 | 61.0 | 1.404e-01 | 8.56e-02 | 4.01e-02 |
| 120 | 16.7 | 39.1 | 5.770e-02 | 3.18e-02 | 1.38e-02 |
| 300 | 6.7 | 18.5 | 2.011e-02 | 1.19e-02 | 4.54e-03 |
| 600 | 3.3 | 10.3 | 1.844e-02 | 7.38e-03 | 3.04e-03 |
| 1200 | 1.7 | 5.5 | 1.044e-02 | 4.57e-03 | 9.79e-04 |
| **1235 (production's 1.62 %)** | **1.62** | **5.2** | **1.058e-02** | **5.12e-03** | **1.31e-03** |

So at production's panel share D5 is real but **secondary**: ≈ 1 point of weight at worst, ≈ 0.4 points at p90,
against D1+D2 and D3 which do not shrink with the axis. D1+D2 (the ordering) and D3 (PerturbMatch inside the norm)
are the terms that carry the divergence in the shipped build.

## 6. Is production's ordering a different-but-defensible estimator, or a bug?

**A bug, on the fixture's evidence — three separate ways, none of them a defensible choice.**

1. **The PerturbMatch coefficient is internally inconsistent.** Production builds the row with per-target donor
   weights and then corrects it as if the row had been built with uniform weights. There is no estimator under
   which "K562 contributes 2.9 % of this row" and "replace 50 % of a 33.3 % K562 share" are both true. On `T007`
   the correction is applied 11.6× too strongly. A defensible weight-first estimator exists — apply the correction
   at `w_k(t)` — and production is not it.
2. **It silently breaks the shape-only property the arm was priced on.** 308's `:nm` makes the row norm identical
   to the incumbent's *by construction*; production restores that norm and then adds an uncontrolled term on top
   of it. Measured (§7): delivered row norms 0.9930 … 1.3407, whole-matrix +6.00 %. That is not an estimator
   choice, it is the norm match failing to do the one thing it is there for.
3. **The division guard protects neither implementation, and they fail differently.** Production falls back to the
   uniform weight only when `nv` is **exactly** 0.0 (`MB:1452`) — unreachable for a row that is merely degenerate:
   `T005`'s centred hct116 row has norm ~1e-17, so `nv > 0.0` is true, the weight is `(1e-17/med)^-1.5 ≈ 1e25`
   with **no cap**, and both rules end up assigning that donor **1.000000** of the weight (measured: `w_308` and
   `w_prod` are both `[0.0, 1.0, 0.0]`). Both implementations then rescale a pure rounding-noise direction back
   to full magnitude — 308 exactly to the incumbent's norm, production to 1.161× it once PerturbMatch is added.
   308's `max(K, 1e-9)` clip (`lib308.py:218`) caps the weight at 3.16e13 but changes nothing here, because the
   weight is renormalised to 1.0 either way. **The two emitted rows have row cosine 0.0523 — unrelated vectors.**
   So case (d) is not a divergence *between* the implementations' guards so much as a shared failure whose
   *output* diverges completely; a real guard would drop such a donor, and neither does.

D7 is a different kind of finding: it is a **deliberate production policy** (`MB:1484-1487`, experiment 061's
broad gene support) that the PLOLO harness has never implemented. It is not a bug in production; it is a fidelity
gap in the instrument, and it is present in every PLOLO reading of every arm on this base.

## 7. Which object does 308's +0.00259 describe, and does this explain the `nmae` sign flip?

**308's `+0.00259 [+0.00000, +0.00527]` describes ORD_A** — blend-first, weights from panel-zeroed rows over
`wcols`, exactly norm-preserving. That object has **never been submitted**. The board point
`Q519I2pLWRCW7PAEx5ul` (0.21771888483638538, Δ +0.0016245) measures **ORD_B**. They differ by 2.700 max on rows
whose median norm is 3.589, with 15/20 rows below cosine 0.99 — so this is not a rounding discrepancy between two
spellings of one estimator, it is two estimators.

**Consequences for the ledger.** BR312's realised factors cannot be used as transfer factors for 308's family:
the "`pds` came in at 2.13× the local raw prediction" compares a board measurement of ORD_B against a PLOLO
prediction for ORD_A. It is not a measurement of a transfer factor at all, and it must not be carried into
`docs/READ-FIRST.md`'s k_pds/k_cos table.

**The `nmae` sign flip — measured, not argued.**

| | delivered row norm, arm / incumbent |
|---|---|
| 308 (ORD_A), full axis | median 1.000000, range **[1.0000000, 1.0000000]** |
| 308 (ORD_A), non-panel | median 1.000000, range **[1.0000000, 1.0000000]** |
| production (ORD_B), full axis | median 1.000000, range **[0.993001, 1.340741]**, 15/20 rows ≠ 1 |
| production (ORD_B), non-panel | median 1.004769, range **[0.806168, 2.154923]** |
| whole-matrix norm ratio | production **1.060008**, 308 **0.9999999996** |

312's preregistration (`experiments/312-normrel-arm/PREREGISTRATION.md:10,14-15`) calls the arm **shape-only** and
predicts an `nmae` leg of **−0.00095, a benefit**. The board returned `score_nmae` **−0.007705** (raw
0.9339736 → 0.9379525, +0.003979), i.e. the opposite sign at 8.1× the magnitude, alongside `score_pds` +0.012302
and `score_mse` +0.005359. **The shipped object is not shape-only**: it moves 15 of 20 fixture rows' delivered
magnitude by up to +34 % (full axis) / +115 % (scored columns) and the panel total by +6.00 %, because the
PerturbMatch correction is added *after* the norm match. `nmae` is a magnitude statistic and the identity
(216/188) prices magnitude changes; an amplitude-changing arm priced against the shape-only gate (+0.0019 rather
than +0.0089) and predicted to *gain* `nmae` is exactly the failure mode the board reported.

This is a **sufficient mechanism, not a proof of the board's number**: the fixture's magnitudes are synthetic,
so +6.00 % is not production's figure, and I did not compute an `nmae` identity leg for the shipped object (that
needs the real rows, which this mandate forbids loading). What is established is that the shipped arm changes
delivered magnitude at all, which the preregistration asserted it could not.

## 8. The parity-correct arm — SPECIFIED, NOT APPLIED

Production is unchanged. `scripts/make_submission_pm_cd4.py` and `scripts/def_pipeline.py` were not touched by
this experiment (`git status` for both is as it was at launch); `Q519I2pLWRCW7PAEx5ul` stands as the result of the
shipped implementation.

**Arm name:** `normrel15-pm` · **Env:** `VCC_NORMREL_PM=blend` (a new knob, default **unset** = today's shipped
weight-first behaviour, so an unset build stays bit-identical to `submissions/312-normrel15`). It is only ever
set **together with** `VCC_NORMREL_EXP=-1.5`; setting it alone must be a hard error, since it changes nothing when
the weighting is off.

**The one function to change:** `_normrel_pre_whiten_effects` (`scripts/make_submission_pm_cd4.py:1418`). Under
`VCC_NORMREL_PM=blend` it takes the prepared `pm` mapping and, per target `t` that both K562-GWPS and
PerturbMatch carry, computes all four steps itself, in ORD_A's order:

```
corr(t)  = alpha * effect_scale * PM_BLEND * (pm_standardized(t) - incumbent_k562(t))   # on pm_target_columns
base_pm  = effect(t, uniform)          + (1/|have|)  * corr(t)      # the blended UNIFORM row = norm reference
cand_pm  = effect(t, normrel_weights)  + w_k562(t)   * corr(t)      # the blended WEIGHTED row
row(t)   = cand_pm * ||base_pm|| / ||cand_pm||                      # norm match LAST -> exactly shape-only
```

with `w_k562(t)` the renormalised normrel weight the same function already computes (`MB:1444-1457`).
**Two supporting edits, both mechanical:** (i) `build_production_effect_matrix` (`MB:2837-2849`) must not apply
`apply_perturbmatch_k562_replacement` a second time to the rows this path already blended — pass the handled-target
set through and keep emitting the same `pm_audit` dict so the postbuild drift check is unchanged; (ii) register
`VCC_NORMREL_PM` in `def_pipeline.py`'s `ARM_BOARD_DELTA` table exactly as 312 registered `VCC_NORMREL_EXP`.
Do **not** "fix" `k562_share` in `apply_perturbmatch_k562_replacement` (`MB:1545`) — that function is also the
incumbent's PerturbMatch path and changing it would move `cd4svd5` and every other arm.

**Two open decisions the arm must settle before it is priced, not during the build:**
- **D5**: whether the weights are taken over the panel-zeroed operator columns (308's object, `wcols`) or over the
  donor's full official coverage (production's, panel included). Parity with 308 means `wcols`. At production's
  panel share this is worth ≤ 0.011 of weight (§5), so it is a second variable — freeze it to 308's choice and say
  so, do not sweep it in the same slot.
- **D7**: production's broad PerturbMatch column policy stays as it is. It is the incumbent's behaviour, and
  changing it would be a second variable. Parity is with 308's *ordering*, on production's *column policy*.

**Expected size of the change, measured on the fixture:** parity-correct vs shipped, max |diff| **1.7297**, median
per-target row cosine **0.9982**, min **0.8736**; and the parity variant's delivered row norm is **1.0000000** for
every target — it restores shape-only exactly. So this is a real re-shaping arm, not a cleanup, and it must be
priced as one.

**Known-answer rows the arm must pass before any build:** (1) `VCC_NORMREL_PM` unset reproduces
`submissions/312-normrel15` byte for byte; (2) `VCC_NORMREL_PM=blend VCC_NORMREL_EXP=0` reproduces `cd4svd5` byte
for byte; (3) the arm's delivered row norm equals the incumbent's for every target to ≤ 1e-6 relative (the
shape-only claim, machine-checked this time rather than asserted); (4) a re-run of 308's `normrel1.5` fold on the
PLOLO harness retrodicts +0.00259 at 0.0e+00.

**Priority note, not a recommendation to build.** The gate for a re-shaping arm is +0.0019 shape-only. 308's local
reading of ORD_A was +0.00259 [+0.00000, +0.00527], priced +0.00104 pessimistic = 0.55× the amplitude gate; the
board has now measured ORD_B at +0.0016245 overall with `pds` +0.012302 and `nmae` −0.007705. Whether the
parity-correct arm is worth a slot depends on re-pricing ORD_A's `nmae` leg on the shipped-scale rows, which
this experiment could not do. That re-pricing should come before the slot.

## 9. Limitations

- The fixture is synthetic. Every **structural** statement (which orderings coincide, which quantities the two
  weight rules read, whether each arm is norm-preserving, what `k562_share` takes) is exact and transfers. Every
  **magnitude** (2.700, +6.00 %, the cosines) is fixture-specific and does not.
- The fixture's panel is 33 % of the axis against production's 1.6 %; §5 sweeps that and reports D5 at production's
  share. D1+D2, D3 and D7 are not panel-share dependent.
- The `nmae` mechanism is established as *present* (the arm is amplitude-changing), not *sized*: no `nmae`
  identity leg was computed for the shipped object, which needs the real rows.
- D7's size in production depends on how many genes PerturbMatch carries that K562-GWPS does not — not measurable
  here without loading the banks. Its existence is established from the code and reproduced on the fixture.
- Everything is measured at exponent −1.5 only, and on `BETA = 0.0` (the shipped value), so the `beta * m_full`
  offset inside the norm match never fires. At `BETA ≠ 0` production's norm match would also rescale that constant
  offset, which 308's does not — a fourth divergence, latent today.
- **Line numbers moved during this run.** `scripts/make_submission_pm_cd4.py` was dirty at launch
  (sha `fd6eb189cd92fea6…`) and became clean at HEAD `0efee8a` (sha `dc024795daf9058a…`) mid-run — a 5-line
  removal *before* line 574 ("Astra's cap-variant patch staged unapplied"). `git diff 7380ad4 HEAD` on that file
  is empty, and `build_v9_pre_whiten_effects`, `_normrel_pre_whiten_effects` and
  `apply_perturbmatch_k562_replacement` are byte-identical across the change. Every number in this report was
  re-run against the live file and reproduced exactly; the line references above are the live ones (all shifted
  −5 from what a reader of the launch-time working copy would see).

## 10. Files

| file | sha256 / note |
|---|---|
| `work/fixture313.npz` | `723ce13c0c3b2e63e6d0af80271c90a34f3f06ee2b3130fd0a4509d636de64df`, seed 3130907 |
| `work/parity313.json` | the full run: KA rows, per-target table, attribution, delivered norms, parity variant |
| `work/parity313.log` | the run transcript |
| `work/scale313.json` | the panel-share sweep |
| `fixture313.py` | the fixture generator |
| `parity313.py` | the three-layer parity harness |
| `scale313.py` | the panel-share sensitivity sweep |

Sources read: `scripts/make_submission_pm_cd4.py` (`dc024795daf9058a…`), `experiments/308-kd-weighted-donors/lib308.py`
(`8bc01663a47271af…`), `experiments/263-donor-selection-arbiters/lib263.py` (`8ced84a480606e1b…`),
`experiments/308-kd-weighted-donors/sweep308.py`, `experiments/062-cd4-post-whiten/run.py`,
`experiments/064-pm-cd4-composition/run.py`, `experiments/061-perturbmatch-hct-count/build_candidate.py`,
`src/vcc2026/multi_donor.py`, `experiments/312-normrel-arm/apply_patch.py`,
`experiments/312-normrel-arm/PREREGISTRATION.md`, `data/board-entries/Q519I2pLWRCW7PAEx5ul.json`,
`data/board-entries/Xv0ayvuoHpsKYWv2iu42.json`. Peak RSS 0.21 GiB. No build, no `vcc`, no `cell-eval2`, no
download, nothing written outside `experiments/313-normrel-parity/`.
