# 300 — PREREGISTRATION: is the whitener's rank-retrieval gain paid by the SPREAD of its training rows?

Written and hashed into `work/PREREGISTRATION.sha256` **before any new operator was fitted and before any new fold
number existed**. Everything read before this hash: `docs/READ-FIRST.md`, `R299`, `R272`, `R286`, 299's scripts
(`run299.py`, `price299.py`, `flex299.py`), `scripts/build_union_whitener.py`, 299's `work/rel/<line>.npz`
(per-row reliabilities, REUSED, not recomputed) and 299's `work/sweep299-<fold>.json` (the `wgt:*` rows carried
here as known answers). No arm below has been run.

## 0. The mechanism under test

R272 §0.10 / R286 §6 / R299 §5 agree on one reading: `pds` is paid by AMPLIFYING directions in which the donor
rows are mostly sampling noise, because the target-specific signal (5–10 % of row energy) is spread over exactly
those directions. If that is right, what the whitener's covariance needs from its training rows is their SPREAD
over many directions, and rows that reproduce (or that align with the bank's leading directions) CONCENTRATE the
covariance and are therefore worth less, not more.

Three predictions, each with a frozen rule (§4):

* **P1 — dose.** The reversed-reliability tilt's Δ`pds` is monotone in the dose inside the shape-only band, and
  the largest shape-only dose beats 299's `d = 1`.
* **P2 — interaction with 062's essential-gene prune.** The prune already removes the MOST reproducible rows there
  are (common-essential knockdowns) and is worth +0.0021 `pds` on the Flex instrument (R299 §6). If the tilt and
  the prune act on one mechanism, tilting ON TOP of the prune gains LESS than tilting on the un-pruned banks.
* **P3 — direct spread manipulation.** Up-weighting rows with LOW leverage on the bank's leading directions —
  reliability untouched — gains at least as much as the reliability tilt, and the two are collinear per target.

If P3 holds the reliability axis was a proxy for spread; if the tilt gains and P3 does not, reliability itself
matters.

## 1. Frozen definitions (no number below is fitted after this hash)

Everything downstream of the whitening operator is `cd4svd5` verbatim (`lib263`/`lib246`/`lib257`, α = 0.982,
γ = 0.65, λ = 0, CD4 `f` = 0.34, `VCC_CD4_DEFLATE_K` = 5, union OFF on PLOLO folds). **One variable per arm: the
weights (or the row set) that teach the operator.** Every arm is fitted at **the incumbent's exact row count and
donor proportions** — all rows of both banks teach every arm — except arm B, whose row set is 062's production
prune, stated as such.

* `u` = the row's rank quantile in [0, 1] of its Spearman–Brown-corrected split-half reliability `R` within its
  own bank (299's `row_weights`, byte-identical code path; a row with no `R` takes the bank's median).
* **Arm A (reversed-reliability tilt):** `w = 1 + d·(0.5 − u)`, `d ∈ {0.5, 1.0, 1.5, 2.0, 2.5}`.
  Σ(0.5 − u) = 0 exactly over a rank quantile, so Σw = n exactly: the row count, donor proportions and coverage
  are the incumbent's on every dose. `d = 1.0` is 299's `wgt:mildrev`; `d = 1.8` is 299's `wgt:steeprev`;
  `d = −1.0` is 299's `wgt:mild`; `d = −1.8` is 299's `wgt:steep`.
* **Arm B (the same tilt on top of 062's prune):** the same `w(d)` at the same doses, on the PRODUCTION row set —
  062's 12,000-rows-per-donor cap (`np.random.default_rng(0)`, `np.sort(rng.choice(n, 12000, False))`, donors
  visited in `("k562_gwps", "orion_hct116", "orion_hek293t")` order, donors under the cap consuming no random
  numbers) followed by dropping every row whose perturbation is in `data/donors/k562_essential/perts.txt`.
  `d = 0` is the prune alone. `u` is re-ranked within the pruned bank.
* **Arm C (direct spread manipulation):** `w = 1 + d·(0.5 − v)`, `d ∈ {1.0, 2.0}`, where `v` is the row's rank
  quantile within its bank of `s_i = ‖P₅₀ z_i‖² / ‖z_i‖²`, the row's energy share on the top 50 eigenvectors of
  **the fold's own cached (un-weighted) operator**, with `z_i` the row centred on the bank's un-weighted column
  mean and scaled by 062's `mult` — i.e. exactly the rows the un-weighted fit sees. `d > 0` up-weights
  LOW-leverage rows. Reliability is not used.
* **Arm D (controls).** (i) *Mirror*: the same weight distribution reflected about 1, i.e. `d → −d`, for A and C
  doses. (ii) *Random-weight*: `w = 1 + d·(0.5 − u_π)` with `u_π` the same `u` vector permuted across the bank's
  rows, seeds 300001 and 300002, at `d = 1.0`.

## 2. Known-answer rows (nothing is priced before they pass; abort if one fails)

| id | claim | tolerance |
|---|---|---|
| KA-0 | `E062.build_baseline` on the rebuilt per-donor rows, float32, == the fold cache's `pre` | 0.0e+00 |
| KA-2a | the un-weighted (`d = 0`) refit's eigenvalues == the cached fold operator's (max rel, λ > 1e-8·max) | ≤ 1e-6 |
| KA-2e | that refit applied to `pre` == the cached operator's base (max row rel) | ≤ 1e-4 |
| KA-3 | `cd4svd5` through this file's plumbing == 263's row on (`cos_agg`, `a23_cos`, `pds`) | ≤ 1e-9 |
| KA-4 | `cd4svd5` per-target `pds` == 263's `pertarget263-<fold>.npz` | 0.0e+00 |
| **KA-W1** | **A at `d = 1.0` per-fold Δ`pds` == 299's `wgt:mildrev`** | **≤ 1e-9** |
| **KA-W2** | **A at `d = −1.0` (mirror) per-fold Δ`pds` == 299's `wgt:mild`** | **≤ 1e-9** |
| **KA-W3** | **A at `d = 1.8` per-fold Δ`pds` == 299's `wgt:steeprev`** | **≤ 1e-9** |
| **KA-B** | **the reproduced 062 prune's kept-row counts == the production operator's 7,784 / 10,807 / 10,803** (R299 §6) | exact |
| R2 | the BR249 retrodiction, `cd4svd5` over `h1union` on 261's POINT set, == 261's row | 0.0e+00 |
| KA-F1 | 273's INCUMBENT on the Flex frame (`flex:prod`) == `instrument273.json` | 0.0e+00 |

KA-B is also the statement of the **residual difference from production**: the PLOLO fold operator is fitted on
the fold's TWO banks, the production operator on THREE, so the row SET is reproduced exactly per bank and the
operator is still a two-bank operator. Anything KA-B does not reproduce is reported as a difference, not hidden.

## 3. Readouts (identical to 299/263/261/265, no new statistic)

Per fold (the three cancer folds `k562_gwps`, `orion_hct116`, `orion_hek293t`): per-target `pds`, per-target own
cosine, `cos_agg`, macro cosine, delivered amplitude, `nmae` identity leg (`nmae_wald89`, `top89_precision`),
262's groups (below-chance / marginal / saturated) and the 0.02-knee crossings, leave-top-3-out retention,
paired target bootstrap (265's `boot`/`boot_folds`, B = 20,000, generator seeded per arm name), 3/3 rule, and
261's price on the three COHERENT scenarios (pessimistic = the end of each leg's measured band that LOWERS the
total for that cell's sign).

## 4. The frozen gate and the frozen prediction rules

**Gate (all four must hold):** pessimistic price ≥ **+0.0019** if the delivered amplitude moves ≤ 0.5 %
(shape-only), else ≥ **+0.0089**; 3/3 folds positive AND the pooled bootstrap interval above 0; leave-top-3-out
retention ≥ 50 %; the arm above BOTH its mirror and its random-weight control by more than half its own Δ`pds`.

* **P1 PASSES** iff the A Δ`pds` is monotone non-decreasing across the five doses AND the largest shape-only dose's
  Δ`pds` exceeds `d = 1.0`'s with the paired (per-target, pooled) bootstrap interval of the DIFFERENCE above 0.
* **P2** is read as `Δ`pds`(B, d*) − Δ`pds`(A, d*)` at `d*` = the best shape-only A dose, with its paired interval
  (the mandate's letter), AND as the interaction `[B(d*) − B(0)] − [A(d*) − A(0)]`, which is the contrast that
  isolates the tilt-on-prune effect from the prune's own effect. Both are reported; the letter is the verdict.
  P2's mechanism prediction is CONFIRMED if the interaction is negative.
* **P3 PASSES** iff Spearman(per-target Δ`pds` of C at `d = 1.0`, of A at `d = 1.0`), pooled over the three folds,
  is positive with a bootstrap interval above 0, AND C's own Δ`pds` at its best dose is ≥ A's at the matched dose.
  C's own gate reading is reported separately.

**Ordering noise.** READ-FIRST §3: the digest-identical ordering-noise SD is **0.004502 overall**, three times
299's whole +0.0014 effect. Every interval quoted in RESULTS.md carries that sentence.

## 5. The Flex reference (273's union-OFF Flex-400 instrument, 85 targets), WITH an interval

`flex:prod` (the production operator) as KA-F1; `flex:full` (the un-pruned matched refit) as the reference; then
the best shape-only A dose and the best C dose. Per-target `pds` is saved for all of them and the arm-vs-reference
difference is given a paired target bootstrap over the 85 targets (B = 20,000). The mandate budgets two refits;
`flex:full` is the third and is required for the reference's per-target values.

## 6. The mechanism row

272's statistic on `orion_hct116`: the share of directions whose independent noise estimate exceeds the
eigenvalue, and the share of the WHITENED PREDICTION ENERGY on them, for the best A dose and the best C dose,
against the incumbent's 76.8 % / 86.8 % and 299's `sel:q0.5` 99.9 %.

## 7. Budget, and what is declared in advance

One operator eigendecomposition (17,942² or 17,806²) per arm per fold, ≈ 5 min alone and ≈ 9–10 min contended;
**at most two fold processes in flight**; peak RSS < 11 GiB per process; `memory_pressure` ≥ 40 % free before each
launch. The 3.5-hour box does not fit every (arm × fold) cell. **The run order is frozen here**, each arm is
flushed to disk as soon as it is measured, and RESULTS.md reports exactly which cells completed and lists every
cell that did not:

1. `cd4svd5`, `STACK` (cached operator, no fit) — KA-3, KA-4, R2
2. `tilt:d1.0` — KA-W1, the primary arm
3. `tilt:d-1.0` — KA-W2, the primary arm's mirror control
4. `tilt:d2.5`, `tilt:d1.5`, `tilt:d0.5`, `tilt:d2.0` — the rest of the A ladder (P1)
5. `prune:d0.0` — the prune alone, KA-B
6. `spread:d1.0` — arm C (P3)
7. `prune:d1.0` — arm B at the primary dose (P2)
8. `spread:d2.0` — arm C's second dose
9. `rand:d1.0:s300001` — the random-weight control
10. `tilt:d1.8` — KA-W3
11. `spread:d-1.0` — arm C's mirror
12. `rand:d1.0:s300002`, `prune:d2.5`, `tilt:d-2.5`, `tilt:d-1.5`, `tilt:d-2.0`, `tilt:d-0.5`, `prune:d1.5`,
    `prune:d0.5`, `prune:d2.0`, `spread:d-2.0` — the remainder, in this order

Rows 1–7 are the minimum that reads P1, P2 and P3 with a control on every fold; the priority order is chosen so
that the gate conditions for the primary arm are complete before any exploratory cell is spent.

## 8. What would make this a candidate

Nothing here is a candidate unless it clears the frozen gate in §4 on all four conditions. If an arm clears, the
result is reported and the coordinator decides; **this experiment does not build, does not prep, does not submit,
and writes nothing outside `experiments/300-spread-tilt/`.**
