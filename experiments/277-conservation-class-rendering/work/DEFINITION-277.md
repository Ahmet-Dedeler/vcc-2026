# DEFINITION-277 — the low-conservation CLASS and the six renderings, frozen before a single arm was scored

Written 2026-09-06, before `render277.py` composed anything. Everything below is fixed; anything decided
afterwards is labelled POST-HOC in `RESULTS.md`.

## 0. Inputs (path + sha256)

| object | path | sha256 |
|---|---|---|
| conservation map | `experiments/276-conservation-map/work/conservation276.csv` | `e5388e3febe6c44cc97a723cf4d5a4b3c27e9e1e964f131a878bfecf2e3464ea` |
| 276's per-fold outcome arrays (h1 map, donor–donor) | `experiments/276-conservation-map/work/predict276.npz` | `a3a9f28966eb6bcf8f2e1b241d7e1ce8a1f1f464050cc9a8cf91a43545cd3267` |
| 262 per-target KA (k562 / hct116 / hek293t / kolf2) | `experiments/262-where-pds-is-lost/work/diag262-*.npz` | `d4c0a35f…` / `d3fce192…` / `43b5d7d4…` / `7b24c9c6…` |
| 273's Flex-400 instrument | `experiments/273-flex-reference-instrument/work/instrument273.npz` | `43a3ec8adea0473fe6d425dec9ed6026918c09fef18aac0623a745f1113a0d9d` |
| 272's noise spectra (per fold + production) | `experiments/272-noise-corrected-operator/work/spectrum272-*.npz` | `77adfad6…` / `2d41485a…` / `5f260094…` / `158beed9…` |
| PLOLO fold caches | `work/plolo/` (byte copies of 262's, `shasum -c` 7/7 OK) | see `work/plolo/SHA256SUMS.txt` |

## 1. The CLASS (declared once, per fold, before any arm is scored)

For fold *f* ∈ {`k562_gwps`, `orion_hct116`, `orion_hek293t`} (and `kolf2`, reported apart):

* eligible targets = the fold's targets (the row order of `L.stage(f)["targets"]`) that have a **finite**
  `C_lolo_<f>` in `conservation276.csv` — the leave-the-held-line-out conservation score, every source from
  the held line removed (276 §1.3, §4.1). One target per cancer fold has no score and is **never** in any class.
* **CLASS30(f)** = the `n_class = floor(0.30 · n_eligible)` eligible targets with the SMALLEST `C_lolo_<f>`.
  Ties are broken by ascending target symbol, so the set is a deterministic function of the CSV.
* `CLASS20` / `CLASS40` are the same rule at 0.20 / 0.40 (V5 only, reported, never selected on).

For the Flex-400 reference (h1-test, 85 targets): the score is `h1_c_lolo` in `predict276.npz`, which 276
built on h1-test's own targets with the identical construction and with **both H1 banks removed**
(`predict276.py:263`, `L.HELD_LINE_SOURCES["H1"]`) — i.e. it is already the H1-excluded score the mandate
asks for; it is used as-is and cited, not recomputed. Same bottom-30 % rule over the targets with a finite
score.

## 2. The renderings (one variable each: what happens to the CD4 residual, or the whitener, ON THE CLASS)

`st = lib257.L.stage(fold)` gives the whitened base `W`, the raw CD4 residual `cd4` (production's
complete-mean3 rows), the H1-pair union `union` at its shipped carrying-row fraction `frac_x`, and the
fold truths. `d, aud = lib246.deflate(cd4, "svd", 5)` is the shipped deflation. Composition is always
`lib246.compose(st, parts)`: whitened base + ONE merged residual add + ONE global preserve.

* **V0 `cd4svd5`** (the incumbent, and the known-answer row):
  `fr = solve_fraction(W, d, 0.34)`; parts `[(d, fr), (union, frac_x)]`.
* **V1 un-deflated CD4 on the class**: the mixed residual `m1` has
  `mean_effect[i] = cd4.mean_effect[i]` for `i ∈ CLASS30`, `= d.mean_effect[i]` otherwise; `gene_votes` and
  `target_condition_count` are `cd4`'s (identical to `d`'s — `lib246._as_res` copies them).
  **The CD4 fraction is re-solved once for `m1`**: `fr1 = solve_fraction(W, m1, 0.34)`, production's
  all-rows rule (`plolo.solve_fraction`, `sizing="all_rows"`): a single scalar = 0.34 / (median over ALL
  panel rows of ‖m1‖ ÷ median ‖W‖). Parts `[(m1, fr1), (union, frac_x)]`.
* **V2 the mirror**: deflated on the class, un-deflated elsewhere. `m2` = `d` on CLASS30, `cd4` off it;
  fraction re-solved the same way.
* **V3 no CD4 residual on the class**: `m3` = `d` with the class rows zeroed (`mean_effect`, `gene_votes`
  and `target_condition_count` set to 0 there); fraction re-solved the same way.
* For V1/V2/V3 a **fixed-fraction diagnostic** (`V1f` / `V2f` / `V3f`) is also reported: the same mixed
  residual carried at V0's own `fr`, so the re-solve's global up/down-dose is separated from the class edit.
  The RESOLVED arm is the primary (production re-solves).
* **V4 the soft noise floor on the class only** (272 A′, `c = 1`): with `V` = the fold operator's
  eigenvectors, `λ` its eigenvalues and `n` 272's independent noise estimate for that operator
  (`spectrum272-<fold>.npz["n"]`), the incumbent weights are `w_inc = (λ + lam)^(-γ/2)` and the damped
  weights `w_flr = (max(λ, 1·n) + lam)^(-γ/2)`. The mixed whitening applies `w_flr` to the CLASS30 rows and
  `w_inc` to every other row, then **one** global on-axis preserve over the mixed matrix
  (`before/after` medians over all rows, exactly `plolo.Operator.apply`'s scalar). Both residual fractions
  are re-solved against the new base. With CLASS = ∅ or `c = 0` this is `cd4svd5` bit for bit (asserted).
* **V5 class-boundary sensitivity for V1 only**: V1 at CLASS20 and CLASS40. Reported, never selected on.

## 3. Controls for V1 (same construction, different row set)

* **RANDOM**: `n_class` eligible targets drawn uniformly without replacement, seeds 20260906 / 20260907 / 20260908.
* **AGREE**: the bottom 30 % of eligible targets by 262's donor–donor agreement (`<fold>|dd` in
  `predict276.npz` — the pipeline's own build-time proxy). V1 must beat this to justify the map.
* **TOPCONS**: the TOP 30 % by `C_lolo` (must not gain).
* V2 is itself the mirror control.

## 4. What is measured, on every arm

`cos_agg`, macro (`a23_cos`), `pds` from `lib246.score` (= `resid218.score_effects`, the fold instrument);
the per-target `pds` and own cosine (a23, three challenge frames) split by CLASS / non-class and by 262's
group (e: base·truth ≤ 0, marginal (0, 0.05], saturated > 0.05); 188's exact `nmae` identity and the
delivered median row norm (`lib257.Identity`, frame A, wald89 pool); paired target bootstraps against V0
(265's `boot265.boot`, B = 2,000, seed 20260906); leave-top-3-out retention (the fold-mean Δ`pds`
recomputed after deleting the three targets with the largest positive Δ).

## 5. Pricing (READ-FIRST §3, coherent-pessimistic)

Class factors, applied per leg at the end of the measured band that LOWERS the total for that cell's sign:
re-shape of occupied rows `k_cos` 3.76 (band 2.15 macro … 3.76 `cos_agg`) and `k_pds` 1.22 (band 0.85 …
1.22) for V1–V3 (BR249 §2 — the deflation's own class); the operator band for V4; the `nmae` identity leg
×1.5 on a loss (BR237: −0.004783 predicted, −0.00716 realised); an amplitude leg if the delivered norm
moves > 0.2 %.

**CLEARED** = pessimistic ≥ +0.0019 over `cd4svd5` **AND** `pds` > 0 on 3/3 cancer folds with the paired
interval above 0 **AND** leave-top-3-out retention ≥ 50 % **AND** every control (random set, agreement-
selected set, mirror) below half the arm **AND** the Flex-400 reading does not contradict the sign.

## 6. Known-answer rows (the run exits non-zero on any failure)

* **KA-1** V0's endpoints == 257's `sweep257-<fold>.json` row `k5_s1` at 1e-9, and its per-target `pds` /
  cosine == 262's `diag262-<fold>.npz` `cd4svd5|pds` / `|cos` at max|d| **0.0e+00** (bit-identical).
* **KA-2** 273's incumbent rows: `INCUMBENT|pds_eff` mean and `cd4h1:off` reproduced from
  `instrument273.npz` at 0.0 (the Flex-400 read starts from 273's own saved arrays where it can, and
  recomputes the composition otherwise, reporting the difference).
* **KA-3** V4 with CLASS = ∅ reproduces V0 at max|d| 0.0e+00.
* **KA-4** BR249's board point retrodicted from this experiment's own V0-vs-`h1union` PLOLO deltas with
  the published class factors (k_cos 3.76 / k_pds 1.22), against the realised +0.0042634.
