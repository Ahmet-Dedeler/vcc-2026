# 299 — PREREGISTRATION: does fitting the whitening operator on RELIABILITY-SELECTED training rows change the held-context transfer?

Frozen 2026-09-06, before any fold number, any arm score and any operator refit exists.
Hash: `work/PREREGISTRATION.sha256` (sha256 of this file), written before `run299.py` is launched.

## 0. What was computed BEFORE this hash, and what was not

Computed before the hash (INPUTS ONLY, no fold, no arm, no operator, no score):

* `rel299.py` rebuilt 267's batch-parity half banks for the three donor lines from 267's exact integer pseudobulk
  sums (`op272.write_half_bank`, verbatim) and wrote, for every bank row present in both halves, the split-half
  cosine, its Spearman–Brown correction, and the row's magnitude (`work/rel/<line>.npz`). These are the SELECTION
  KEYS; they are a deterministic function of data on disk and of the formula fixed in the mandate, and no arm, no
  operator, no fold and no score existed when they were computed.
* The one diagnostic that fixed the definition: the cosine had to be taken over the mapped official columns with
  the 300 PANEL columns EXCLUDED, because that is what `E062.prepare_main_source` (and therefore 271) produces —
  it zeroes the panel columns. With panel columns included the KA against 271 reads max |Δ| 2.1e-02 … 5.1e-02;
  excluded, 1.1e-09 … 2.0e-09 (float32 storage of 271's rows). The panel-EXCLUDED statistic is the frozen one.

Nothing else. No fold was staged, no operator was fitted, no `pds` was computed.

## 1. Known-answer rows (all must pass; a failure aborts the arm they gate)

| id | claim | tolerance |
|---|---|---|
| KA-R | per-row split-half reliability from 267's halves reproduces 271's panel-row reliabilities on the shared rows | ≤ 5e-09 (271's rows are stored float32) — **measured 2.03e-09 / 1.10e-09 / 1.10e-09 for k562_gwps / orion_hct116 / orion_hek293t before this hash** |
| KA-0 | `E062.build_baseline` on the rebuilt per-donor rows, cast float32, equals the fold cache's `pre` | 0.0e+00 |
| KA-1 | transform + preserve through the cached fold operator equals the staged base | 0.0e+00 (bit-identical) |
| KA-2 | the operator refit by `fit299` on the INCUMBENT's row set (every row of both donor banks) reproduces the cached operator: (a) eigenvalues above 1e-8 × max, max relative difference; (b) its whitened base applied to the fold's `pre`, max row relative difference; (c) `cos_agg` / `a23_cos` / `pds` through `cd4svd5`'s construction | (a) ≤ 1e-6, (b) ≤ 1e-4, (c) ≤ 1e-5 — 272's KA-2a/2e/2f tolerances verbatim |
| KA-3 | `cd4svd5` through this file's staging reproduces 263's `pds`, `cos_agg`, `a23_cos` per fold | ≤ 1e-9 |
| KA-4 | `cd4svd5` per-target `pds` bit-identical to 263's `work/pertarget263-<fold>.npz` | 0.0e+00 |
| R2 | the BR249 retrodiction: `cd4svd5` vs `STACK` (h1union) priced on the POINT set reproduces 257/261's +0.00472 | 0.0e+00 against 261's own `compare` on the same fold caches |

KA-2 is run ONCE (on `orion_hct116`, the fold 272 also refit) because it costs a full 5½-minute eigendecomposition
and validates the fit code path, not the fold. Its `full` arm doubles as the incumbent-equivalence row in score terms.

## 2. The selection (frozen)

Per donor line j:

* `r_j(row)` = cos(A_j, B_j) of the row's two batch-parity halves, over the official-axis columns both halves map
  onto **minus the 300 panel columns** (production's own support; 271's definition — see §0).
  `R_j = 2 r / (1 + r)` (Spearman–Brown) is reported; the RANKING is by `r_j` (monotone in `R_j`, so the two rank
  identically).
* `m_j(row)` = ‖(A_j + B_j)/2‖ on the same columns — the magnitude of the same object; the `mag` control's key.
* **Eligible pool** = the line's bank rows that are NOT panel targets and have a defined `r_j` (present in both
  halves with two non-zero half norms). Rows with an undefined `r_j` are in no selected arm.
* **Arm row set for donor j at fraction q** = ⌈q · |pool_j|⌉ rows chosen from the pool + **every panel-target row
  the bank carries** (coverage of the panel is preserved exactly, as in the incumbent).
  * `sel:q` — the TOP-q by `r_j` (most reliable)
  * `rev:q` — the BOTTOM-q by `r_j` (least reliable)
  * `mag:q` — the TOP-q by `m_j` (largest response: the 285 / 288 confound)
  * `rand:q` — a uniform random subset of the pool of the same size (seed = 299000 + fold index)
  * `full` — every row of both banks = the incumbent (KA-2)
* q ∈ {0.5, 0.25} for `sel`; the controls are run at q = 0.5 (the same size as the primary arm they control).
  Ties in `r_j` are broken by the row's name (`np.lexsort`), so every arm is deterministic.

Counts fixed by the rule and the data (stated here, before any fold number):

| line | bank rows | with r | panel rows in bank | pool | q0.5 rows | q0.25 rows |
|---|---:|---:|---:|---:|---:|---:|
| k562_gwps | 9,787 | 9,786 | 269 | 9,517 | 5,028 | 2,649 |
| orion_hct116 | 17,767 | 17,765 | 290 | 17,475 | 9,028 | 4,659 |
| orion_hek293t | 17,977 | 17,976 | 294 | 17,682 | 9,135 | 4,715 |

Per fold (donor proportions, incumbent → q0.5): held `k562_gwps` 17,767 : 17,977 = 0.497 : 0.503 → 9,028 : 9,135 =
0.497 : 0.503; held `orion_hct116` 9,787 : 17,977 = 0.353 : 0.647 → 5,028 : 9,135 = 0.355 : 0.645; held
`orion_hek293t` 9,787 : 17,767 = 0.355 : 0.645 → 5,028 : 9,028 = 0.358 : 0.642.

## 3. The fit and everything downstream

`fit299` is `op272.fit_operator` (= 062 `fit_union_whitener`'s arithmetic) restricted to a row index set per bank:
each bank's selected rows centred on **the selected rows' own column mean** (062's recipe applied to the chosen row
set; identical to the bank mean when the set is every row, so KA-2 is intact), scaled by `ref / raw_median_bank`
with `raw_median_bank` the FULL bank's `prepare_main_source` raw median (a property of the bank, not the selection —
held verbatim), `C = Bᵀ B / N`, `eigh`, eigenvalues floored at max × 1e-12, eigenvectors round-tripped through
float32 exactly as 062 stores them. The only change to 272's code is `scipy.linalg.eigh(..., overwrite_a=True,
driver='evd')` in place of `np.linalg.eigh`, to hold the peak RSS under the mandate's cap; KA-2 tests that path.

Everything downstream is `cd4svd5` verbatim (263's `lib263` / 246's `lib246`): α 0.982, γ 0.65, λ 0, the fold's
own axis and `wcols`, ONE global on-axis preserve, CD4 deflated k = 5 re-solved at f = 0.34 (all-rows rule), the H1
pair union re-solved by the carrying-row rule at g = 0.28, one add, one preserve, panel columns zeroed, scored by
`lib246.score`, `nmae` by 188's exact identity through `lib257.Identity`. **One variable: which rows teach the
operator.**

## 4. Readings

Per arm × fold: `pds`, `cos_agg`, `a23_cos`, `a15_sgn_w`, delivered amplitude, `nmae_wald89`, `top89_precision`,
per-target `pds` and cosine (saved to `work/pertarget299-<fold>.npz`), 262's groups (below-chance / marginal /
saturated by the base·truth cosine), the 0.02-knee crossings, leave-top-3-out retention (261's `retention`), the
paired target bootstrap over the pooled per-target Δ`pds` (265's `boot`, B = 10,000, seed 20260905), and the
six-member price by 261's recipe on the three coherent scenarios (PESSIMISTIC / POINT / OPTIMISTIC), with
k_cos ∈ {1.36, 3.76} and k_pds ∈ {0.85, 1.22} — an operator change is a RE-SHAPE of occupied rows, so the
pessimistic side takes k_cos 3.76 on a cosine loss and k_pds 0.85 on a `pds` gain.

## 5. The gate (frozen, READ-FIRST §3)

An arm CLEARS only if all four hold:

1. **price** PESSIMISTIC ≥ +0.0019 if the delivered amplitude is within 0.5 % of the incumbent's (shape-only), else
   ≥ +0.0089;
2. **3/3** Δ`pds` > 0 on all three cancer folds AND the pooled paired-target bootstrap interval above 0;
3. **retention** leave-top-3-out ≥ 50 %;
4. **controls** Δ`pds`(`sel`) exceeds Δ`pds`(`rand`), Δ`pds`(`mag`) and Δ`pds`(`rev`) by more than half of
   Δ`pds`(`sel`) itself.

Each condition's outcome is stated per arm, pass or fail. A well-supported negative is the expected outcome under
the mechanism evidence on the other side (286's E3 lost with a cleaner covariance) and is a valid result.

## 6. Scope actually attempted, and the order (time box: HARD STOP 21:30Z)

Arms are measured in this priority order so that a stop at any time leaves the primaries read on every fold:
per fold `cd4svd5` (KA) → `STACK` (the R2 retrodiction leg) → `sel:q0.5` → `sel:q0.25` → `rand:q0.5` → `mag:q0.5`
→ `rev:q0.5`; `full` (KA-2) on `orion_hct116` only. Each arm's row is appended to
`work/sweep299-<fold>.json` as soon as it is measured, so a killed process loses only the arm in flight.

**Deviation from the mandate, declared here:** the mandate asks for `rand` at 3 seeds. One full eigendecomposition
costs 5½ minutes and the box allows ≈ 16 of them; `rand` is therefore run with ONE seed per fold (three independent
random draws across the three folds), and extra seeds are added only if the clock allows. §4's Flex re-read (273's
union-OFF instrument) is attempted only if ≥ 40 min remain after the PLOLO folds, per the mandate.

## 7. What would make this a positive result, in one sentence

`sel:q0.5` or `sel:q0.25` gains `pds` on 3/3 cancer folds, above all three controls by more than half its own
Δ`pds`, prices ≥ +0.0019 pessimistic at shape-only amplitude, and keeps ≥ 50 % after its three best targets are
removed — i.e. the whitener's eigenbasis is improved by teaching it with reproducible rows rather than all rows.
