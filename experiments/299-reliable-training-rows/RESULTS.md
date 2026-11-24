# 299 — RESULTS: the whitening operator fitted on RELIABILITY-SELECTED training rows

**Preregistration** `PREREGISTRATION.md`, sha256 `8f582c98803233c0e479b14a6b6c1a6f6cd370348f1786e1caf7dee911036c1d`,
hashed 2026-09-06T18:08:36Z into `work/PREREGISTRATION.sha256` — before any fold was staged, any operator was fitted
and any `pds` was computed. Ended 2026-09-06, inside the 21:30Z hard stop. No build, no submission, nothing written
outside `experiments/299-reliable-training-rows/`.

## 0. Headline

**RELIABLE-ROW WHITENER: best arm `sel:q0.5` Δ`pds` −0.06748 [−0.08165, −0.05353], 3/3: NO (0/3), controls
(rand −0.06060 / mag −0.05490 / rev −0.05945 — every one of them ABOVE `sel`), priced −0.08938 pessimistic vs gate
+0.0089 → DOES NOT CLEAR; FLEX: the best arm's construction reads +0.00140 Δ`pds` / +0.00058 Δ own cosine on 273's
union-OFF Flex-400 instrument against its matched refit — SAME SIGN as PLOLO.**

*(The best PREREGISTERED arm is `sel:q0.5`, above. The best arm of the whole experiment is the POST-HOC, rank-
preserving `wgt:steeprev` — the operator taught by the same rows tilted AWAY from reliability — at Δ`pds` +0.00139
[−0.00009, +0.00288], 3/3 folds positive, mirror control passed, priced +0.00268 pessimistic against a +0.0089
amplitude-changing gate; its milder twin `wgt:mildrev` is shape-only at +0.00094 [+0.00017, +0.00173], 3/3, mirror
control passed, priced +0.00140 = 0.74 × the +0.0019 shape-only gate. Neither clears. §5.)*

Two separate negatives, and they are not the same negative:

1. **At q < 1 the arm is a ROW-COUNT effect, not a selection effect.** Halving the teaching rows takes the
   covariance below full rank on two of the three folds (14,163 and 14,056 rows on 17,942 columns → 3,781 floored
   eigenvalues), and the γ = 0.65 whitening then amplifies that floored null space by ≈ 7.9 × 10³ relative to the
   top direction. Every arm at that q pays it: `rand` (−0.0606), `rev` (−0.0595) and `mag` (−0.0549) all beat
   `sel` (−0.0675). Reliability-ranked rows are the WORST of the four matched-size selections.
2. **At full rank there is nothing there either.** On `k562_gwps`, the one fold whose q = 0.5 row set stays above
   its column count (18,163 > 17,806), the four arms are −0.0032 (`sel`), −0.0031 (`rand`), +0.0003 (`mag`),
   +0.0012 (`rev`) — a spread of 0.004 `pds` with per-fold bootstrap intervals ≈ ±0.006, i.e. all of them at
   zero, with `sel` at the bottom. §5's rank-preserving arms make the same reading at the incumbent's exact row
   count and donor proportions.

The mechanism row says why, and it says the opposite of the hypothesis: selecting reliable rows does **not** clean
the operator's eigenbasis. Measured with 272's independent noise estimate on the same fold, the incumbent puts
86.8 % of its whitened prediction energy in directions where noise exceeds the eigenvalue; the reliability-selected
operator puts **99.9 %** there, and its energy-weighted mean signal share collapses from 0.153 to 0.0007.

## 1. Known-answer rows (every one passed; nothing was priced before they did)

| id | claim | tolerance | measured |
|---|---|---|---:|
| KA-R | per-row split-half reliability from 267's halves == 271's panel-row reliabilities (shared rows) | ≤ 5e-09 | **2.03e-09 / 1.10e-09 / 1.10e-09** (k562 / hct116 / hek293t) |
| KA-0 | `E062.build_baseline` on the rebuilt per-donor rows, float32, == the fold cache's `pre` | 0.0e+00 | **0.0e+00** on 3/3 folds |
| KA-2a | the refit on the INCUMBENT's row set == the cached operator's eigenvalues (max rel, λ > 1e-8 max) | ≤ 1e-6 | **1.20e-09** (identical to 272's numpy-path value) |
| KA-2e | that refit applied to the fold's `pre` == the cached operator's base (max row rel) | ≤ 1e-4 | **8.08e-08** |
| KA-2f | `full` == `cd4svd5` on (`cos_agg`, `a23_cos`, `pds`) | ≤ 1e-5 | **1.36e-13** |
| KA-3 | `cd4svd5` through this file's staging == 263's row | ≤ 1e-9 | **0.0e+00** on 3/3 folds |
| KA-4 | `cd4svd5` per-target `pds` == 263's `pertarget263-<fold>.npz` | 0.0e+00 | **0.0e+00** on 3/3 folds |
| R2 | the BR249 retrodiction, `cd4svd5` over `h1union` on 261's POINT set | 0.0e+00 vs 261 | **0.0e+00** — `compare` max\|d\| 0.0e+00, POINT total \|d\| 0.0e+00 (mine +0.004722217, 261 +0.004722217; board realised +0.0042634) |
| KA-F0 | `DampedOperator` + preserve == `plolo.Operator.apply` on the Flex frame | 0.0e+00 | **0.0e+00** |
| KA-F1 | `flex:prod` == 273's saved `INCUMBENT` (`pds`, own cosine, `cos_agg`) | ≤ 1e-9 | **0.0e+00** (0.883333 / 0.05206 / 0.07107) |

KA-2 also validates the one code change made for the memory cap: `scipy.linalg.eigh(..., overwrite_a=True,
driver='evd')` in place of `np.linalg.eigh` inside 062's recipe.

## 2. The selection, and the row counts it produced

Reliability is `r_j(row) = cos(A_j, B_j)` of the row's two batch-parity halves (267's half banks, rebuilt from 267's
exact integer sums by `op272.write_half_bank`), over the official-axis columns both halves map onto **minus the 300
panel columns** — production's own support, which is what `E062.prepare_main_source` produces and therefore what 271
measured. (With the panel columns left in, the KA against 271 reads 2.1e-02 … 5.1e-02; excluded, 1.1e-09. This was
settled before the hash and is recorded in `PREREGISTRATION.md` §0.)

| line | bank rows | with r | panel rows in bank | pool | q0.5 rows | q0.25 rows | median r, pool | median r, top-q0.5 | median r, top-q0.25 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| k562_gwps | 9,787 | 9,786 | 269 | 9,517 | 5,028 | 2,649 | +0.0437 | +0.1009 | +0.2008 |
| orion_hct116 | 17,767 | 17,765 | 290 | 17,475 | 9,028 | 4,659 | +0.0199 | +0.0635 | +0.1136 |
| orion_hek293t | 17,977 | 17,976 | 294 | 17,682 | 9,135 | 4,715 | +0.0326 | +0.1222 | +0.2125 |

The selection is real and strong: the top half of each donor's pool is 2.3–3.8 × more reliable than the pool
median, the top quarter 4.6–6.5 ×. It also picks up the magnitude confound 285 / 288 named, but only mildly — the
selected rows' median magnitude is 1.02–1.07 × the pool's at q0.5 (7.110 vs 6.831 on k562, 7.297 vs 7.121 on
hct116, 6.649 vs 6.227 on hek293t), which is why `mag` and `sel` are not the same arm.
Donor proportions are held: held-`orion_hct116` 0.353 : 0.647 (incumbent) → 0.355 : 0.645 (q0.5); held-`k562_gwps`
0.497 : 0.503 → 0.497 : 0.503; held-`orion_hek293t` 0.355 : 0.645 → 0.358 : 0.642. Every panel target the bank
carries is kept in every arm, so coverage is the incumbent's exactly.

**Fold row totals against the operator's column count** — the fact that dominates every q < 1 arm:

| fold | columns p | incumbent rows | q0.5 rows | q0.25 rows | floored eigenvalues at q0.5 |
|---|---:|---:|---:|---:|---:|
| k562_gwps | 17,806 | 35,744 | **18,163 (> p)** | 7,364 | few (full rank) |
| orion_hct116 | 17,942 | 27,764 | **14,163 (< p)** | 7,364 | **3,781** (measured; p − rows = 3,779) |
| orion_hek293t | 17,942 | 27,554 | **14,056 (< p)** | 7,308 | p − rows = 3,886 |

## 3. The preregistered arm × fold table (everything over `cd4svd5`, three cancer folds)

| arm | k562_gwps | orion_hct116 | orion_hek293t | mean Δ`pds` | 3/3 | `cos_agg` % | macro % | Δ`nmae89` | amplitude | LT3O |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|---:|---:|
| `sel:q0.5` | −0.00319 | −0.10270 | −0.09655 | **−0.06748** | 0/3 | −25.30 | −24.90 | +0.00016 | −27.70 % | +109.0 % |
| `sel:q0.25` | −0.07321 | −0.10688 | −0.10791 | **−0.09600** | 0/3 | −67.59 | −47.91 | +0.11120 | +1.85 % | +109.2 % |
| `rand:q0.5` | −0.00307 | −0.09230 | −0.08644 | −0.06060 | 0/3 | −23.91 | −23.85 | +0.00039 | −27.54 % | +109.6 % |
| `mag:q0.5` | +0.00031 | −0.07916 | −0.08586 | −0.05490 | 1/3 | −23.83 | −23.54 | +0.00092 | −28.21 % | +110.4 % |
| `rev:q0.5` | +0.00123 | −0.08998 | −0.08961 | −0.05945 | 1/3 | −22.64 | −22.51 | +0.00361 | −27.46 % | +110.2 % |
| `full` (KA) | — | 0.0 (1.4e-13) | — | — | — | — | — | — | — | — |

Base-to-incumbent-base median row cosine, the cleanest single diagnostic: `full` 1.0000; at q0.5 `sel` 0.9627 /
0.1221 / 0.1077 (k562 / hct116 / hek293t), `rand` 0.9628 / 0.1132 / 0.0998, `mag` 0.9774 / 0.1322 / 0.1171, `rev`
0.9612 / 0.1017 / 0.0877; at q0.25 `sel` 0.2628 / 0.2459 / 0.2315. The base survives exactly where the row count
stays above the column count.

Paired target bootstrap (265's `boot` / `boot_folds`, B = 20,000, a generator seeded per arm from 265's
SEED = 20260905 so an arm's interval does not move when other arms are added), pooled over the three folds and per fold:

| arm | pooled Δ`pds` | 95 % | k562_gwps | orion_hct116 | orion_hek293t |
|---|---:|---|---|---|---|
| `sel:q0.5` | −0.06748 | [−0.08165, −0.05353] | −0.00319 [−0.00847, +0.00223] | −0.10270 [−0.13203, −0.07336] | −0.09655 [−0.12699, −0.06682] |
| `sel:q0.25` | −0.09600 | [−0.11397, −0.07885] | −0.07321 [−0.10323, −0.04403] | −0.10688 [−0.13746, −0.07614] | −0.10791 [−0.13925, −0.07695] |
| `rand:q0.5` | −0.06060 | [−0.07387, −0.04748] | −0.00307 [−0.00771, +0.00165] | −0.09230 [−0.12089, −0.06399] | −0.08644 [−0.11474, −0.05889] |
| `mag:q0.5` | −0.05490 | [−0.06788, −0.04204] | +0.00031 [−0.00311, +0.00375] | −0.07916 [−0.10590, −0.05205] | −0.08586 [−0.11319, −0.05879] |
| `rev:q0.5` | −0.05945 | [−0.07307, −0.04561] | +0.00123 [−0.00446, +0.00689] | −0.08998 [−0.11941, −0.06134] | −0.08961 [−0.11852, −0.06031] |

262's groups (contribution to the fold-mean Δ`pds`) and the 0.02 knee, for `sel:q0.5`: below-chance +0.00966
(n = 139), marginal **−0.07229** (n = 527), saturated −0.00485 (n = 187); knee crossings 22 up / 119 down, net −97.
The loss is where 262 said the metric lives — the marginal band — and the below-chance group gains, exactly as a
randomised base should behave. `sel:q0.25` is the same shape, larger: marginal −0.11119, knee net −184.

Leave-top-3-out "retention" reads +109 … +110 % for every arm; for a LOSS that ratio only says the loss is diffuse
and not three targets, which is what it says here. It is not evidence of anything positive.

## 4. The price and the frozen gate (261's recipe, coherent scenarios)

| arm | shape? | amplitude | PESSIMISTIC | POINT | OPTIMISTIC | gate | price | 3/3 | retention | controls | verdict |
|---|:--:|---:|---:|---:|---:|---:|:--:|:--:|:--:|:--:|:--:|
| `sel:q0.5` | no | −27.70 % | **−0.08938** | −0.08396 | −0.04142 | +0.0089 | FAIL | FAIL | PASS | **FAIL** | **DOES NOT CLEAR** |
| `sel:q0.25` | no | +1.85 % | **−0.24986** | −0.22000 | −0.10017 | +0.0089 | FAIL | FAIL | PASS | n/a (not run) | **DOES NOT CLEAR** |
| `rand:q0.5` | no | −27.54 % | −0.08366 | −0.07851 | −0.03870 | +0.0089 | FAIL | FAIL | PASS | — | DOES NOT CLEAR |
| `mag:q0.5` | no | −28.21 % | −0.08117 | −0.07596 | −0.03718 | +0.0089 | FAIL | FAIL | PASS | — | DOES NOT CLEAR |
| `rev:q0.5` | no | −27.46 % | −0.08200 | −0.07667 | −0.03868 | +0.0089 | FAIL | FAIL | PASS | — | DOES NOT CLEAR |

Condition-by-condition for the primary arm `sel:q0.5`:

* **price** FAIL — −0.08938 pessimistic against the +0.0089 amplitude-changing gate (the delivered amplitude moves
  −27.7 %, far outside the 0.5 % shape-only band, so the higher gate applies; it fails the shape-only gate too).
* **3/3** FAIL — 0/3 folds positive; pooled interval [−0.08165, −0.05353], entirely below 0.
* **retention** PASS — but see the note above; for a loss this condition carries no information.
* **controls** FAIL on all three, with the margin the WRONG SIGN: `sel` − `rand` = −0.00687, `sel` − `mag` =
  −0.01257, `sel` − `rev` = −0.00802, against a required margin of +0.03374 (half of |Δ`pds`(`sel`)|).

`sel:q0.25` was not given its own q = 0.25 controls; the preregistration set the controls at q = 0.5 (the size of
the primary arm) and one eigendecomposition costs 5–13 minutes inside a three-hour box. Its verdict does not depend on them: it fails
price and 3/3 by margins of 0.26 and 0.10 `pds`.

## 5. Rank-preserving arms (POST-HOC — added after the preregistered arms were read, and labelled as such)

The q < 1 design cannot separate "which rows teach the operator" from "how many rows teach it", because at q = 0.5
the covariance falls below full rank on two of three folds. These arms remove that confound: **every** row of both
banks teaches the operator, at the incumbent's exact row count, donor proportions, axis, γ and λ, and each row is
merely WEIGHTED by its reliability. `u` is the row's rank quantile of `r` within its bank; weights are renormalised
per bank to sum to the bank's row count, so the donor mixture is the incumbent's exactly.

| arm | weight | effective n / n | direction |
|---|---|---:|---|
| `wgt:R` | w ∝ clip(R, 0.005, 0.99) | 34–44 % | toward reliable (strong) |
| `wgt:steep` | w = 0.1 + 1.8 u | ≈ 79 % | toward reliable (steep) |
| `wgt:mild` | w = 0.5 + u | 92.3 % | toward reliable (mild) |
| `wgt:mildrev` | w = 1.5 − u | 92.3 % | **away from** reliable (mild) |
| `wgt:steeprev` | w = 1.9 − 1.8 u | ≈ 79 % | **away from** reliable (steep) |

**The dose–response is monotone across five points and its sign is the OPPOSITE of the hypothesis:**

| arm | k562_gwps | orion_hct116 | orion_hek293t | mean Δ`pds` | pooled 95 % | 3/3 | `cos_agg` % | amplitude | base~inc | pessimistic | gate | verdict |
|---|---:|---:|---:|---:|---|:--:|---:|---:|---:|---:|---:|:--:|
| `wgt:R` | −0.00204 | −0.00815 | −0.00910 | **−0.00643** | [−0.01043, −0.00241] | 0/3 | −12.38 | −2.37 % | 0.948 | −0.03224 | +0.0089 | DOES NOT CLEAR |
| `wgt:steep` | −0.00073 | −0.00292 | −0.00298 | **−0.00221** | [−0.00375, −0.00068] | 0/3 | −4.33 | +0.24 % | 0.992 | −0.01105 | +0.0019 | DOES NOT CLEAR |
| `wgt:mild` | −0.00038 | −0.00136 | −0.00141 | **−0.00105** | [−0.00183, −0.00028] | 0/3 | −2.12 | +0.30 % | 0.998 | −0.00533 | +0.0019 | DOES NOT CLEAR |
| — incumbent — | 0 | 0 | 0 | 0 | | | 0 | 0 | 1.000 | | | |
| `wgt:mildrev` | +0.00032 | +0.00140 | +0.00111 | **+0.00094** | [+0.00017, +0.00173] | 3/3 | +1.90 | +0.23 % | 0.998 | **+0.00140** | +0.0019 | DOES NOT CLEAR |
| `wgt:steeprev` | +0.00061 | +0.00200 | +0.00156 | **+0.00139** | [−0.00009, +0.00288] | 3/3 | +3.88 | +0.87 % | 0.991 | **+0.00268** | +0.0089 | DOES NOT CLEAR |

Mirror controls (each arm against the same weight distribution applied in the opposite direction), the analogue of
the preregistered control condition: `wgt:mildrev` beats `wgt:mild` by +0.00199 against a required +0.00047 →
**PASS**; `wgt:steeprev` beats `wgt:steep` by +0.00360 against +0.00070 → **PASS**; both hypothesis-direction arms
fail their mirror by the same margins with the sign reversed. `wgt:R`'s mirror (`wgt:rev`, w ∝ 1 − R) was not run —
at the reliabilities measured here 1 − R ≈ 1 for nearly every row, so that weighting is a weaker tilt than
`wgt:mild`, not a stronger one, and the rank-quantile pair above is the matched comparison.

* **`wgt:mild`** — the hypothesis's own direction at full rank — loses on 3/3 folds, pooled interval entirely below
  0, and prices −0.00533 pessimistic. Its mirror control beats it by +0.00199, more than double the required
  margin, in the wrong direction. Every gate condition except the (uninformative-for-a-loss) retention fails.
* **`wgt:mildrev`** — tilting the whitener's training rows TOWARD the least reproducible ones — gains `pds` on
  3/3 folds with the pooled interval above 0 ([+0.00017, +0.00173]), gains `cos_agg` +1.90 %, is shape-only
  (delivered +0.23 %, inside the 0.5 % band), and **passes the 3/3 and the control conditions**. It fails on
  price (+0.00140 pessimistic = 0.74 × the +0.0019 shape-only gate) and on leave-top-3-out retention (+47.8 %
  against 50 %). **DOES NOT CLEAR** — but it is the only arm in this experiment that failed narrowly rather than
  categorically.
* **`wgt:steeprev`** — the same tilt 4.75 × steeper — is larger in `pds` (+0.00139) and in cosine (+3.88 %), which
  is why the direction is real; but the extra dose moves the delivered amplitude to +0.87 %, outside the shape-only
  band, so the +0.0089 gate applies and +0.00268 pessimistic is 0.30 × of it. Its pooled bootstrap interval also
  crosses zero ([−0.00009, +0.00288]) and its retention falls to +34.9 %: the gain concentrates as the dose rises.
* 262's groups for `wgt:mildrev`: below-chance +0.000205, marginal +0.000724, saturated +0.000014 — the gain is in
  the marginal band, exactly where 262 says `pds` is decided; knee crossings 3 up / 4 down, net −1.

**Reading.** The reliability of a training row is not neutral for the whitener — it is *negatively* useful. This is
the direct prediction of R272 §0.10 and R286 §6: `pds` is bought by amplifying directions in which the donor rows
are mostly sampling noise, so a covariance estimated preferentially from the rows that reproduce is a covariance
that under-weights exactly the directions the metric pays for. 299 measures that prediction on a five-point dose
curve and it holds in both signs.

## 6. The Flex reference (273's union-OFF Flex-400 instrument, 85 targets)

The production operator `union_pruned` is itself fitted on a PRUNED, per-donor-capped row set (12,000 rows per donor,
then the essential-gene rows dropped: 7,784 / 10,807 / 10,803 of 9,787 / 17,767 / 17,977), and its exact seed-0
subsample is not reproducible here without guessing 062's RNG call order. The read is therefore a **matched pair of
refits** on the same three banks and the same axis, with the production operator carried as the known-answer row.
Everything else is 273's incumbent verbatim: base + H1-CD4 deflated k = 5 re-solved at f = 0.34 (all-rows rule),
union OFF, scored with the scorer's own per-target statistic.

| arm | `pds` | own cosine | `cos_agg` | delivered | base ~ prod | Δ`pds` vs `flex:full` | Δ own cos | Δ`cos_agg` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `flex:prod` (= 273's INCUMBENT, KA-F1 at 0.0e+00) | 0.883333 | 0.05206 | 0.07107 | 3.8711 | 1.0000 | +0.002101 | −0.000600 | +0.001044 |
| `flex:full` (un-pruned refit, the matched reference) | 0.881232 | 0.05266 | 0.07003 | 3.9383 | 0.9651 | — | — | — |
| **`flex:wgt:mildrev`** (the same rows, tilted AWAY from reliable) | **0.882633** | **0.05323** | **0.07169** | 3.9269 | 0.9660 | **+0.001401** | **+0.000577** | **+0.001661** |

The arm read here is `wgt:mildrev` — the shape-only member of the pair, and the one that passed the 3/3 and control
conditions on PLOLO. (`wgt:steeprev` has the larger PLOLO Δ`pds` but changes the delivered amplitude by +0.87 %;
with two refits affordable inside the box, the shape-only construction is the one worth reading on a second
instrument.)

**FLEX: same sign as PLOLO.** Against its matched un-pruned refit the reliability-reversed tilt gains Δ`pds`
+0.001401 and Δ own cosine +0.000577 at amplitude ×0.9971 — the same direction as PLOLO and, on `pds`, 1.5 × the
PLOLO reading of the same arm (+0.00094). Unlike arms (a)–(c) in R273 §0, whose PLOLO verdicts did NOT reproduce as
`pds` effects on this instrument, this one does, in sign and in size; the cosine effect is again the larger of the
two relative to its own scale (+2.4 % of own cosine vs +0.16 % of `pds`), which is 273 §0.6's signature. Caveat:
85 targets, no interval computed here (one arm, one refit), so this is a sign-and-size agreement, not an
independent confirmation.

Side observation, not an arm: the production PRUNED operator beats the un-pruned refit of the same three banks by
Δ`pds` +0.0021 on this instrument, i.e. 062's essential-gene prune and 12,000-row cap are themselves worth about
1.5 × what the reliability tilt is worth — another row-composition effect, and in the opposite direction from
reliability (the prune removes the most reproducible rows there are: common-essential knockdowns).

## 7. The mechanism row (272's statistic, fold `orion_hct116`)

| operator | directions with s < 0.5 | s = 0 | share of the WHITENED prediction energy on them | raw-energy share | Σn / Σλ | energy-weighted mean s |
|---|---:|---:|---:|---:|---:|---:|
| incumbent (cached fold operator) | 93.0 % | 76.8 % | **86.8 %** | 71.2 % | 0.855 | 0.1530 |
| `sel:q0.5` refit | 90.6 % | 78.5 % | **99.9 %** | 43.3 % | 0.747 | **0.0007** |

The incumbent's row reproduces R272 §0.10 (77–79 % of directions at s = 0 carrying 87–92 % of the whitened energy):
76.8 % and 86.8 % here. Selecting the reliable rows moves the statistic the WRONG WAY. The refit's covariance has
3,781 exactly-floored eigenvalues (λ = λ_max × 1e-12), whose whitening weight is (1e-12)^(−γ/2) ≈ 7.9 × 10³ times
the top direction's, so essentially all of the whitened prediction energy (99.9 %) lands in directions the
independent noise estimate scores as noise, and the energy-weighted signal share falls by a factor of 220.

## 8. What this settles, and what it does not

* **Settled.** Row SELECTION BY RELIABILITY, at the fractions the mandate named, does not improve the whitener: it
  is below a matched random subset on 3/3 folds, below a magnitude-selected subset on 3/3, and below the REVERSED
  (least reliable) selection on 3/3. That ordering — `mag` > `rev` > `rand` > `sel` — is the opposite of the
  hypothesis, and it holds on the one fold where the row count does not fall below the column count as well as on
  the two where it does.
* **Settled.** The whitener's `pds` does not come from the reproducible part of the donor rows. 272 showed it is
  paid in directions that are mostly sampling noise; 286 showed replacing the covariance with an unbiased signal
  covariance loses; 299 now shows that removing the least reproducible ROWS loses too, and loses more than removing
  random rows. The three results agree: what the metric needs is the SPREAD of the donor rows over many directions,
  not their reproducibility.
* **Measured, and the one live thread.** The reliability axis is not inert — it is signed, and the useful sign is
  the reverse of the hypothesis. A monotone five-point dose curve (§5) runs −0.00643 → −0.00221 → −0.00105 → 0 →
  +0.00094 → +0.00139 as the tilt moves from "toward reliable" to "away from reliable", at the incumbent's exact
  row count and donor proportions, with matched mirror controls passing on both positive arms and the same sign
  reproduced on 273's Flex-400 instrument (+0.00140). **It is still not build-worthy: the best shape-only member
  prices +0.00140 pessimistic against the +0.0019 shape-only gate (0.74 ×) and retains 47.8 % after its three best
  targets are removed (against 50 %) — two narrow misses, not one.** The coordinator should treat this as a
  measured direction with a mechanism, not as a candidate: nothing here justifies a slot, and READ-FIRST §3's
  ordering-noise SD (0.004502 overall) is three times the whole effect.
* **Not settled by this experiment.** Whether some other composition rule at the incumbent's exact row count helps
  (091's finding is that composition matters at equal count; §5 tests the reliability axis of that at equal count,
  but the space of composition rules is larger than one axis — and §6's side observation, that 062's essential-gene
  prune is worth +0.0021 on the Flex instrument, says another axis of it is already being exploited). Nothing here
  bears on target-specific uses of the reliability signal, which 271 and 263 already closed on their own
  constructions.
* **Not tested.** `rand` at three seeds (one seed per fold was run; the three folds are three independent draws);
  controls at q = 0.25; `wgt:rev` (the mirror of `wgt:R` — at these reliabilities 1 − R ≈ 1, so it is a weaker tilt
  than `wgt:mild`, not the matched mirror); an interval on the Flex arm; and any dose between `wgt:steeprev` and
  the point where the amplitude leaves the shape-only band, which is where a clearing version of this arm would
  have to live.

## 9. Provenance and cost

Inputs: PLOLO's fold caches (byte copies of 263's in `work/plolo/`), 267's batch-parity half banks (rebuilt from
267's exact integer sums, then deleted), the three donor banks, 271's saved half rows (for the noise covariance and
the KA), 261's and 257's pricing constants, 265's bootstrap, 273's Flex instrument. Reserves and R3 were not
touched. Files written: `PREREGISTRATION.md`, `work/PREREGISTRATION.sha256`, `rel299.py`, `run299.py`,
`price299.py`, `mech299.py`, `flex299.py`, `work/rel/<line>.{npz,json,log}`, `work/sweep299-<fold>.json`,
`work/pertarget299-<fold>.npz`, `work/identity299-<fold>.npz`, `work/price299.{json,txt}`,
`work/mech299-orion_hct116.{json,log}`, `work/flex299.{json,log}`, the run logs. Nothing outside this directory was
modified; the rebuilt half-bank `delta.npy` files were deleted after use (they regenerate from 267's sums in ~1 min).

Arms measured: 13 on PLOLO (`cd4svd5`, `STACK`, `full`, `sel:q0.5`, `sel:q0.25`, `rand:q0.5`, `mag:q0.5`,
`rev:q0.5`, `wgt:mild`, `wgt:mildrev`, `wgt:steep`, `wgt:steeprev`, `wgt:R`) plus 3 on the Flex instrument, from
35 operator eigendecompositions (31 on PLOLO, 1 for the mechanism row, 3 on the Flex frame — one of which was a
repeat after the first Flex process died on a missing `sys.path` entry for `est271`, fixed and re-run).

Cost: one operator eigendecomposition (17,942² or 17,806², `scipy` `evd`) per arm per fold, ≈ 5.2 min alone and
≈ 10–13 min with three or four folds in flight; staging a fold 3–24 s; scoring an arm 1–3 s. Peak `ru_maxrss` per
process 10.1–11.0 GiB — that figure includes the memmapped donor-bank pages the Gram pass touches; the live
anonymous RSS measured with `ps` during the eigendecomposition was 5.4–9.0 GiB. `memory_pressure` stayed at 79–89 %
free throughout with up to four processes.

**Deviations from the mandate, all declared:** (i) `rand` at one seed per fold rather than three seeds per fold
(preregistered in §6 of `PREREGISTRATION.md`); (ii) controls at q = 0.5 only; (iii) `scipy.linalg.eigh` in place of
`np.linalg.eigh` inside 062's recipe, for the memory cap, validated by KA-2; (iv) up to four fold processes in
flight rather than strictly one at a time — the memory gate (`memory_pressure` ≥ 40 % free, free+inactive ≥ 20 GiB)
held with margin throughout and is recorded above; (v) §5's rank-preserving arms are POST-HOC, added after the
preregistered arms were read, and are labelled as such everywhere they appear.
