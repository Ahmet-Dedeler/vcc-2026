# 286 — RESULTS: a cross-product response estimator (shared directions from independent views, full-rank residual path)

**Verdict: DOES NOT CLEAR — a well-supported negative, with a measured mechanism for why.**

> **HELD VIEWS (CD4, 280's four-donor GWCD4i bank, Stim48hr):** Δcos **+0.00042 [−0.00164, +0.00254]** vs the two-donor
> mean at the strongest floor (φ = 0.25; shuffled-pairing **−0.02531 [−0.03079, −0.02002]**, random-direction
> **−0.02926 [−0.03545, −0.02337]**), Δ`pds` **+0.00401 [−0.00212, +0.01012]** against a random-direction control of
> **+0.00317 [−0.00070, +0.00697]** — indistinguishable from the incumbent and from its own control.
> **HELD STATE:** Rest Δcos **−0.00383 [−0.00581, −0.00190]**, Δ`pds` **−0.00886 [−0.01412, −0.00359]**; Stim8hr Δcos
> **−0.00493 [−0.00699, −0.00294]**, Δ`pds` **−0.00159 [−0.00707, +0.00390]** — the gain does not survive a state change.
> **HELD CONTEXTS (PLOLO):** best arm `e1:phi0.25:nm` Δ`pds` **−0.00801 [−0.01632, +0.00010]**, 3/3: **NO** (1/3),
> priced **−0.02065 pessimistic** vs gate **+0.0089** → **DOES NOT CLEAR** (0/12 cells; the best cell of any kind is
> 272's already-failed operator floor at −0.00201).
> **FLEX REFERENCE (273's h1-test, union OFF, 85 targets):** the same arm reads Δ`pds` **−0.00644 [−0.01835,
> +0.00518]** and Δ own cosine **−0.00130 [−0.00270, +0.00010]** against the matched rebuilt mean — same sign as
> PLOLO on both statistics, so the PLOLO loss is the prediction's transfer, not the truth's noise; the low-rank
> arms are catastrophic on both instruments (−0.16 / −0.19 `pds` here, −0.14 / −0.18 on PLOLO).

Everything under `experiments/286-crossproduct-estimator/`. No build, no `vcc`, no `cell-eval2`, nothing written under
`submissions/`, `data/harness/` or `data/`. Both reserves untouched (`experiments/267-*/work/reserve/`, the Broad-16
lines JVE127 / MG63 / SLR23 / UMRC3 were never read).

## 0. What was frozen, and when

| file | sha256 | frozen |
|---|---|---|
| `PREREGISTRATION.md` | `fef5ad0966ba2f776b4bcf3e3447abf18de9201954d0997998f079f1fb4a96a2` | 2026-09-06T13:54:07Z |
| `est286.py` | `f9d6075baeb13ae288a4ba28ff02a0d968410b66368ce83913bcf0873d6fdf10`* | 2026-09-06T13:54:07Z |

*recorded verbatim in `work/prereg286.sha256` and asserted by every stage script; every output json carries
`est286_sha256` and `analyseB286.py` asserts it did not change between folds. **No held number was measured before
this freeze** — the first held number is Stage A's, run afterwards. `est286.py` was never edited after the freeze.

## 1. The known-answer table

| # | row | required | measured | verdict |
|---|---|---|---|---|
| KA-M | `‖Ĉ − Ŝ − 2N̂‖_F / ‖Ĉ‖_F` (the exact moment identity) | 0 | **5.0e-17 / 9.5e-17 / 9.1e-17** (k562 / hct116 / hek293t), Stage A **8.6e-17** | PASS |
| KA-M′ | the mandate's literal `‖Ĉ − Ŝ − N̂‖_F/‖Ĉ‖_F` | stated, not 0 | **0.2300 / 0.4397 / 0.3999**; it equals `‖N̂‖_F/‖Ĉ‖_F`, the noise share of the half row's covariance | stated |
| KA-P | this experiment's production rows == 271's `real` rows | 0.0 | **0.0e+00** on all three lines | PASS |
| KA-H | the half rows == 271's `rows-<line>.npz` on the fitted columns | float32 | rel max\|d\| **2.4e-08 … 4.6e-08** | PASS |
| KA-Mi | the in-place moment path == `est286.moments`/`moment_ka` on k562 | bit-identical | **bit-identical** on `S`, `N`, `S_shuf`, `N_shuf` | PASS |
| KA-C | `fit286` E0 rows == the fold's staged production rows | 1e-6 | **5.6e-08 … 9.2e-08** (float32 storage) | PASS |
| KA-263 | `cd4svd5` == 263 on (`cos_agg`, `a23_cos`, `pds`) | < 1e-9 | **2.0e-11 / 1.8e-11 / 2.0e-11**; per-target `pds` array **bit-identical** to 263's | PASS |
| KA-257 | `STACK` == 257 | < 1e-9 | **1.8e-11 / 1.6e-11 / 1.8e-11** | PASS |
| KA-271 | `relw` and `ebrow` (and `wiener_own`) == 271 | < 1e-9 | **≤ 1.3e-11**; and through 261's pricer they reproduce 271's published cells to the digit: `relw` Δ`pds` −0.003700, PESS −0.01117; `ebrow` +0.001402 [+0.00049, +0.00237], PESS +0.00050, retention 31 %; `wiener_own` +0.004911, PESS −0.01328 | PASS |
| KA-BR249 | the BR249 retrodiction through 261's recipe | +0.00472 at 0.0e+00 vs 257 | **+0.00472**, \|d\| **1.9e-10**; Δ`pds` +0.012263 == 257's; bands `pds`/`mse`/overall all OK; board +0.0042634 | PASS |
| KA-E1id | E1 at φ = 1 (every gain floored to 1) == E0 | 0 | **3.8e-04 relative** — a float32 round-trip through the Cholesky factor, NOT an algebraic error (the gains are exactly 1 by construction); reported, not passed at 1e-4 | **measured, not passed** |
| KA-shuf | leading eigenvalue of Ŝ vs Ŝ_shuf | ratio ≫ 1 | **19.34 / 7.17 / 12.62** (k562 / hct116 / hek293t); tr Ŝ/tr Ĉ **+0.1119 / +0.0441 / +0.0588** real vs **−0.0000 / +0.0007 / +0.0017** shuffled | PASS |
| KA-A1 | Stage A: 280's Q2 `full_k2` per-target array reproduced | 0.0 | **0.0e+00** over 253 targets | PASS |
| KA-C1/C2 | Stage C: `BASE` and `INCUMBENT` == 273's `instrument273.json` | < 1e-12 | **0.0e+00** on both | PASS |
| KA-C0/C3 | Stage C: the `cand_h1` rebuild and `E0:rebuilt` | 273's own H-0 | rel max\|d\| **0.325304** == 273's H-0 **0.325304** to the digit; `E0:rebuilt` `pds` **0.8822128851540618** == 273's `mean:rebuilt` **0.882213** exactly | PASS (as reproductions of 273's published, explicitly non-fatal H-0) |

**One known-answer row is exempted, explicitly.** `lib263`'s `KA-0` compares `E062.build_baseline`'s float64 rebuild,
cast to float32, against the fold cache's `pre` at tolerance 0.0. It reads **4.1e-08 / 5.5e-08 / 5.8e-08** here rather
than 0.0, because this experiment **built its own fold caches** (so `fold["pre"]` is the freshly computed float64
array) instead of byte-copying 262's saved float32 npz, as 271 did. It is a cache-provenance artefact, not an
arithmetic difference: in the same run `cd4svd5` reproduces 263 at 2.0e-11 with a **bit-identical** per-target `pds`
array, `STACK` reproduces 257, and `relw`/`ebrow`/`wiener_own` reproduce 271 — all at ≤ 2.0e-11. The exemption is
named in `priceB286.py` and is the only one; no other failure is exempt.

## 2. The object that was built

Per donor line j, from 267's batch-parity half banks (rebuilt from its exact integer pseudobulk sums by
`op272.write_half_bank` = `halfrows271.write_half_bank` = `truth267.write_half_bank`, verbatim), on **all the line's
bank rows except the 300 panel targets** — held targets never enter the fit:

    Ŝ = (UᵀU − DᵀD/4)/n   N̂ = DᵀD/(4n)   Ĉ = (UᵀU + DᵀD/4)/n = Ŝ + 2N̂     U = (A+B)/2, D = A−B

| line | fitted rows n | fitted columns m | panel rows excluded | tr Ŝ / tr Ĉ | ‖N̂‖/‖Ĉ‖ | λ_max(Ŝ) | share λ < 0 | clipped energy | Wiener gain median | share gain < 0.10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 9,517 | 7,448 | 269 | +0.1119 | 0.2300 | 3.2968 | 48.8 % | 0.406 | 0.0841 | 0.503 |
| `orion_hct116` | 17,475 | 17,806 | 290 | +0.0441 | 0.4397 | 1.5916 | 49.9 % | 0.456 | 0.0001 | 0.514 |
| `orion_hek293t` | 17,682 | 17,806 | 294 | +0.0588 | 0.3999 | 4.5896 | 49.9 % | 0.440 | 0.0001 | 0.511 |
| CD4 (Stage A, 399 non-panel Stim48hr targets) | 399 | 17,471 (rank 1,592) | — | +0.0630 | 0.4257 | 2.6499 | 75.0 % of the 1,592 | 0.458 | 5.9e-05 | 0.744 |

**Both preregistered predictions about the object are confirmed.** P1: Ŝ is far from PSD — **49–50 %** of its
eigenvalues are negative on every PLOLO line and the clipped energy is **0.41–0.46** (prediction: > 30 % and > 0.10).
P2: the Wiener gains are small in most directions — the median gain is **0.0001–0.084** and **50–51 %** of directions
sit below 0.10 (prediction: median < 0.25). This is exactly why the floor exists; without it E1 is close to a
projection onto a small subspace.

The construction is genuinely a different object from 271's and 272's: the fitted eigenvectors are **not** the
operator's. Principal cosines between Ŝ's top-100 directions and the fold operator's top-100 have median
**0.52–0.70** (mean 0.53–0.60) — overlapping, but far from the identity a re-weighting of the operator's own basis
would give.

## 3. Stage A — held VIEWS (CD4, 280's bank)

Fit: the 399 non-panel Stim48hr targets carried by all four donors, views = the three disjoint 2+2 donor splits.
Test: the 253 panel targets, all three disjoint pairings × both directions; the prediction is built from **two** donors
and scored against the **other two** donors' consensus, in 280's residual cosine and in the scorer's per-target rank
statistic with all panel targets as distractors. E0 = the two-donor mean = 280's Q2 `k = 2` object (KA-A1 at 0.0e+00,
median cosine **+0.03539**). The estimator uses J = 0.5 (the filtered row averages 2 donors; N̂ is calibrated for the
four-donor mean). Because n = 399 ≪ m = 17,471, all the algebra is done **exactly** in an orthonormal basis of
span([U; D]) (p = 1,592) with G_φ = P G̃_φ Pᵀ + φ(I − PPᵀ).

| arm (Stim48hr, 253 targets) | median cosine | **Δcos vs E0 [95 % CI]** | mean `pds` | **Δ`pds` vs E0 [95 % CI]** |
|---|---:|---|---:|---|
| `E0` two-donor mean | +0.04822 | — | 0.68870 | — |
| `E1 φ = 0` (un-floored Wiener) | +0.03304 | −0.02362 [−0.03092, −0.01665] | 0.67383 | −0.01486 [−0.03130, +0.00136] |
| `E1 φ = 0.10` | +0.04175 | −0.01007 [−0.01450, −0.00569] | 0.68863 | −0.00006 [−0.01105, +0.01078] |
| **`E1 φ = 0.25`** | **+0.05343** | **+0.00042 [−0.00164, +0.00254]** | **0.69271** | **+0.00401 [−0.00212, +0.01012]** |
| `E2 r = 20` (low rank, no floor) | +0.01014 | −0.05598 [−0.06902, −0.04328] | 0.54040 | −0.14830 [−0.17586, −0.12003] |
| `E2 r = 100` | +0.01400 | −0.05282 [−0.06533, −0.04108] | 0.56067 | −0.12803 [−0.15468, −0.10158] |
| Cm `:mirror` (φ = 0.10) | +0.05174 | +0.00120 [−0.00266, +0.00511] | 0.64729 | −0.04141 [−0.04817, −0.03510] |
| Cs `:shufpair` (φ = 0.10) | +0.03223 | −0.02531 [−0.03079, −0.02002] | 0.68419 | −0.00450 [−0.01548, +0.00622] |
| Cr `:randdir` (φ = 0.10) | +0.02870 | −0.02926 [−0.03545, −0.02337] | 0.69187 | **+0.00317 [−0.00070, +0.00697]** |

**Reading.** On the cosine the learnt directions are clearly doing something the controls are not: at φ = 0.10 the arm
loses 0.010 while the shuffled-pairing and random-direction controls lose 0.025–0.029. But the arm never *beats* the
two-donor mean: at the floor where it stops losing cosine (φ = 0.25) the gain is **+0.00042, an interval that contains
zero**, and the `pds` movement (+0.00401 [−0.00212, +0.01012]) is **matched by the random-direction control**
(+0.00317 [−0.00070, +0.00697]). That is the preregistered failure mode P4 verbatim: whatever `pds` moves at a strong
floor is the *shrinkage*, not the learnt directions.

**P3 is confirmed in its second half and refuted in its first**: E2 (low rank, no floor) is catastrophic — Δ`pds`
−0.148 / −0.128 at r = 20 / 100 — so the floor matters enormously; but E1-with-floor does **not** beat E0, so the
comparison is between a bad arm and a neutral one.

**Held state** (Ŝ fitted on Stim48hr, scored on the held donors' other states — Rest 251 targets, Stim8hr 260):

| arm | Rest Δcos | Rest Δ`pds` | Stim8hr Δcos | Stim8hr Δ`pds` |
|---|---|---|---|---|
| `E1 φ = 0.10` | −0.01744 [−0.02186, −0.01331] | −0.02150 [−0.03106, −0.01190] | −0.01910 [−0.02427, −0.01419] | −0.01067 [−0.02069, −0.00049] |
| `E1 φ = 0.25` | −0.00383 [−0.00581, −0.00190] | −0.00886 [−0.01412, −0.00359] | −0.00493 [−0.00699, −0.00294] | −0.00159 [−0.00707, +0.00390] |
| `E2 r = 100` | −0.05400 | −0.13917 | −0.05800 | −0.13055 |

The Stim48hr cosine result (+0.0004) does **not** survive a state change: every held-state cell is negative and the
φ = 0.25 Rest interval excludes zero. The estimator's directions are state-specific in the way 282 described the
underlying data.

**Control note.** `extra280` covers Stim48hr only, so a wrong-**state** Ŝ cannot be fitted on non-panel targets; Cw is
realised here as (i) the cross-state application above and (ii) `:randdir` — the *same* eigenvalue spectrum on random
directions, which matches the arm's energy and shrinkage exactly and differs only in the directions. That is a
stronger control than a wrong-source covariance, and the arm does not beat it on `pds`.

**Second view definition.** A guide split (guide-rank 0 vs guide-rank 1, pooled over donors; 150 targets carry two
`keep_for_DE` guides in all four donors) gives tr Ŝ/tr Ĉ **+0.0358** against the donor split's +0.0630 — the same
object, weaker, consistent with 282's finding that the reproducible part of a small CRISPRi pseudobulk is a donor ×
guide interaction.

## 4. Stage B — held CONTEXTS (PLOLO, three cancer folds)

Full tables: `work/tablesB286.md` (M1–M5), `work/priceB286.md` (P1), logs `work/sweepB286-<fold>.log`.
25 arm rows per fold; the estimator replaces each donor's un-whitened production row and rides `cd4svd5` verbatim
(whitening through `plolo.Operator.apply`, CD4 deflated k = 5 re-solved at f = 0.34, the H1 pair at g = 0.28, one add,
one preserve); E3 replaces the operator instead.

### 4.1 Δ`pds` and the six-member price over `cd4svd5`

| cell | Δ`pds` K562 / HCT116 / HEK293T | pooled 95 % CI | Δ`cos_agg` % | delivered % | Δ`nmae` id | **PESSIMISTIC** | POINT | gate | cleared |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| `ebrow` (271) | +0.0004 / +0.0025 / +0.0014 | +0.00140 [+0.00050, +0.00235] | +0.725 | +0.04 | +0.0008 | **+0.00050** | +0.00194 | 0.0019 | no (retention 31 %) |
| `E4:wiener_own` (271) | +0.0024 / −0.0002 / +0.0125 | +0.00491 [−0.00291, +0.01269] | −3.042 | +0.13 | +0.0185 | **−0.01328** | −0.00939 | 0.0019 | **no** |
| `E5:opfloor_c1` (272) | −0.0045 / −0.0003 / −0.0016 | −0.00215 [−0.00391, −0.00044] | +0.600 | +2.86 | −0.0014 | **−0.00201** | −0.00068 | 0.0089 | **no** |
| `e1:phi0:nm` | −0.0468 / −0.0283 / −0.0165 | −0.03053 [−0.04382, −0.01700] | −21.320 | −5.62 | +0.0061 | **−0.06639** | −0.06099 | 0.0089 | **no** |
| `e1:phi0.1:nm` | −0.0334 / −0.0168 / −0.0071 | −0.01911 [−0.03028, −0.00777] | −14.623 | −4.91 | +0.0049 | **−0.04549** | −0.04170 | 0.0089 | **no** |
| **`e1:phi0.25:nm`** | −0.0184 / −0.0062 / **+0.0006** | **−0.00801 [−0.01632, +0.00010]** | −6.402 | −3.03 | +0.0032 | **−0.02065** | −0.01884 | 0.0089 | **no** |
| `e1:phi0.25:lit` | −0.0204 / −0.0065 / +0.0014 | −0.00851 [−0.01672, −0.00025] | −5.177 | −42.22 | +0.0002 | **−0.02757** | −0.02643 | 0.0089 | **no** |
| `e2:r100:nm` | −0.1512 / −0.1437 / −0.1177 | −0.13755 [−0.15808, −0.11728] | −39.809 | +0.41 | +0.0114 | **−0.15982** | −0.14975 | 0.0019 | **no** |
| `e2:r20:nm` | −0.1975 / −0.1894 / −0.1499 | −0.17893 [−0.20103, −0.15713] | −39.205 | −2.08 | +0.0161 | **−0.17851** | −0.16791 | 0.0089 | **no** |
| `E3:c0` (whitener, floor 1e-6 λ_max) | −0.0138 / −0.0324 / −0.0470 | −0.03106 [−0.03888, −0.02339] | −22.346 | +2.70 | −0.0025 | **−0.06832** | −0.06287 | 0.0089 | **no** |
| `E3:c1` (whitener, floor at the noise level) | −0.0070 / −0.0061 / −0.0080 | −0.00706 [−0.01096, −0.00321] | −2.356 | +1.21 | −0.0017 | **−0.00871** | −0.00773 | 0.0089 | **no** |

**Controls (all `:nm`, pooled Δ`pds` vs `cd4svd5`):** at φ = 0.25 `mirror` −0.00530, `rand` −0.01069, `shufpair`
−0.01721, `wrong` −0.02438; at φ = 0.10 `mirror` −0.00648, `rand` −0.02898, `shufpair` −0.03439, `wrong` −0.04352.

**Gate: 0/12 cells clear.** No cell reaches the pessimistic threshold; no cell is positive on 3/3 folds with the
pooled interval above 0 (the best is 1/3). The best pessimistic total of any cell is **−0.00201**, and it belongs to
272's already-failed operator floor, not to a cross-product arm. The best cross-product cell is `e1:phi0.25:nm` at
**−0.02065 pessimistic** — an order of magnitude below the +0.0089 gate its delivered-amplitude class demands.
(Note: the `:nm` norm-match is applied to the **un-whitened** row; after whitening, the CD4/union re-solve and the
render, the delivered norm still moves by −3.0 %, so 261's classifier puts these arms in the amplitude-changing class.
Even priced at the shape-only gate of +0.0019 the cell is 11× below it.)

Retention is reported for completeness but is not interpretable for a **negative** Δ`pds`: leave-top-3-out "retention"
above 100 % simply means dropping the three best-moved targets makes the arm worse.

### 4.2 Where the change lands (262's groups, the 0.02 knee)

| arm | anti-aligned (base·truth ≤ 0) | marginal band (0, 0.05] | saturated > 0.05 | band Δ own cos | up / down across 0.02 | net |
|---|---:|---:|---:|---:|---|---:|
| `e1:phi0.25:nm` | **+0.0423** | **−0.0241** | +0.0005 | −0.00179 | 23 / 60 | −7.0 % |
| `e1:phi0.1:nm` | +0.0795 | −0.0509 | −0.0009 | −0.00381 | 27 / 79 | −9.9 % |
| `e2:r100:nm` | +0.1094 | −0.2176 | −0.0895 | −0.00990 | 60 / 140 | −15.2 % |
| `E3:c1` | −0.0007 | −0.0102 | −0.0004 | +0.00018 | 22 / 13 | +1.7 % |
| Cs `:shufpair` (φ = 0.25) | +0.0603 | −0.0430 | −0.0007 | −0.00422 | 18 / 85 | −12.7 % |

This is 273 §5's BLURS signature, at scale: every cross-product arm **helps the 17.5 % of targets that point away**
(+0.042 … +0.109) and **hurts the marginal band** (−0.024 … −0.218), where `pds` is actually decided. The knee count
is the arm's `pds`: 23 targets pushed up across 0.02 against 60 pushed down. Because the shuffled-pairing control
shows the same shape (+0.060 / −0.043), the group signature is a property of *shrinking toward a subspace*, not of
the learnt directions.

## 5. Stage C — the direction-resolving reference (273's h1-test instrument, H1 union OFF)

85 targets, 10x Flex, 400 cells per target (split-half residual reliability median +0.243 — the one reference on disk
that resolves directions, R273 §0.6). The same arms, composed through `lib273` verbatim: 170's staged `cand_h1`,
`union_pruned`, the PerturbMatch 0.5 blend variant 273 selected, the re-extracted H1 CD4 deflated k = 5 re-solved at
f = 0.34, `plolo.deliver` + the scorer's read-back, `resid218._per_target`. The H1 union was never loaded.

**Known answers.** `BASE` `pds` **0.871148** and `INCUMBENT` **0.883333** reproduce 273's `instrument273.json` at
**0.0e+00**. The two rows that did *not* pass at tolerance 0 are 273's OWN published, explicitly non-fatal finding
reproduced to the digit: the fold instrument's row construction differs from the production builder's `cand_h1` on the
PerturbMatch columns — my rebuild reads a relative max\|d\| of **0.325304**, exactly 273's H-0 value 0.325304, and my
`E0:rebuilt` `pds` is **0.8822128851540618**, exactly 273's `mean:rebuilt` 0.882213. So both are reproductions of a
known instrument property, not failures of this experiment; every donor-row arm is therefore read against BOTH the
`INCUMBENT` and the matched rebuilt mean `E0:rebuilt`, as 273 did (its R1).

| arm | `pds` | own cos | **Δ`pds` vs `E0:rebuilt` [95 % CI]** | **Δ own cos vs `E0:rebuilt`** | Δ`pds` vs `INCUMBENT` |
|---|---:|---:|---|---|---|
| `E0:rebuilt` | 0.882213 | 0.05168 | — | — | −0.00112 [−0.00392, +0.00154] |
| `e1:phi0:nm` | 0.867507 | 0.04670 | −0.01471 [−0.03235, +0.00126] | −0.00498 [−0.00712, −0.00287] | −0.01583 [−0.03389, +0.00070] |
| `e1:phi0.1:nm` | 0.871429 | 0.04840 | −0.01078 [−0.02549, +0.00350] | −0.00328 [−0.00511, −0.00147] | −0.01190 [−0.02731, +0.00252] |
| **`e1:phi0.25:nm`** | 0.875770 | 0.05038 | **−0.00644 [−0.01835, +0.00518]** | **−0.00130 [−0.00270, +0.00010]** | −0.00756 [−0.02003, +0.00448] |
| `e1:phi0.25:lit` | 0.875490 | 0.05024 | −0.00672 [−0.01821, +0.00434] | −0.00144 [−0.00280, −0.00009] | −0.00784 [−0.01961, +0.00364] |
| Cm `:mirror` (φ = 0.10) | 0.877171 | 0.04896 | −0.00504 [−0.01246, +0.00154] | −0.00272 [−0.00374, −0.00179] | −0.00616 [−0.01373, +0.00098] |
| `e2:r100:nm` | 0.719188 | 0.03284 | −0.16303 [−0.22675, −0.10196] | −0.01884 [−0.02790, −0.01035] | −0.16415 [−0.22675, −0.10434] |
| `e2:r20:nm` | 0.690896 | 0.03488 | −0.19132 [−0.25686, −0.12857] | −0.01680 [−0.02710, −0.00683] | −0.19244 [−0.25841, −0.13025] |

**Reading (273 §5's BLURS question).** The PLOLO verdict is **not** the truth's noise. Every arm's sign reproduces on
the direction-resolving reference: `e1:phi0.25:nm` reads −0.0064 `pds` and −0.0013 own cosine on Flex-400 against
−0.0080 `pds` and −6.4 % `cos_agg` on PLOLO, and the low-rank arms are catastrophic on both (−0.16 / −0.19 `pds`
here, −0.14 / −0.18 there). Both Flex intervals for the best arm cross zero — 85 targets is a small n — but there is
no cell anywhere in this experiment where the estimator reads as a `pds` gain on one instrument and a loss on the
other. Unlike 273's damping and shrinkage arms, which blurred (cosine effect ≫ `pds` effect), the cross-product arms
move `pds` and cosine **together and in the same direction**, on both instruments: this is a transfer effect, and it
is negative.

**Leakage, stated.** The Stage B fits exclude the 300 A/B/C panel targets but not h1-test's 85 targets: the overlap
with each line's fitted bank rows is **82/9,517 (0.86 %)**, **84/17,475 (0.48 %)** and **84/17,682 (0.48 %)**. A
covariance fitted on 9.5–17.7 k rows moves by O(1/n) per row, so the leak is small, and it can only *flatter* the
arms — which lose anyway.

## 6. Mechanism (measured)

1. **Principal angles.** Ŝ's leading directions are NOT the fitted operator's: principal cosines of the top-100
   subspaces have median **0.52 / 0.66 / 0.62** (k562 / hct116 / hek293t as donors; means 0.53–0.60). The object 272
   said it did not measure was in fact measured here.
2. **Where the prediction's energy lives.** With the pooled Ŝ and N̂ on the fold's operator columns:
   **78.1–79.8 % of the directions have λ⁺(Ŝ) < ν (the independent noise level), and those directions carry
   80.8–81.8 % of the incumbent's whitened prediction energy.** Ŝ's top-100 directions carry **1.0–1.1 %** of it;
   the directions where Ŝ's eigenvalue is *negative* carry **61.5–62.3 %**. 272 measured the same shape in the
   FITTED eigenbasis (77–79 % of directions, 87–92 % of the energy). Re-estimating the eigenvectors from independent
   views does **not** change the picture: the incumbent buys `pds` in directions where an unbiased signal-covariance
   estimate says there is no reproducible signal at all, and any estimator that shrinks those directions gives up more
   `pds` in the marginal band than it recovers on the anti-aligned group.
3. **The signal is real but tiny and spread.** tr Ŝ / tr Ĉ is **+0.044 … +0.112** (CD4 +0.063): between 4 % and 11 %
   of a half row's covariance is shared response. Its top-100 eigen-directions hold only **19–32 %** of Ŝ's positive
   energy — the shared response itself is high-rank, which is why the rank-20/100 truncations destroy retrieval.
4. The `:shufpair` control confirms the estimate is not an artefact: tr Ŝ/tr Ĉ falls from +0.044…+0.112 to
   −0.000…+0.002 and the leading eigenvalue falls by **7.2× / 12.6× / 19.3×**. The directions are real. They are just
   not the directions `pds` is paid in.

## 7. Limitations, and what was not done

* **Memory budget exceeded.** The mandate set ~10 GiB peak RSS. Three processes exceeded it: the first two
  `stage286.py` runs peaked at **23.2 / 23.8 GiB** (the library moment path materialises Ŝ, N̂ and Ĉ plus temporaries
  at m = 17.8 k). This was found and fixed the same session — `moments_inplace` holds two m × m float64 buffers and
  reproduces the library path **bit-for-bit** (KA-Mi); the re-run peaks at **4.65 GiB**. The eigendecomposition
  processes (`fit286.py`) peaked at **11.0–11.1 GiB** and the fold sweeps at **11.5–13.8 GiB** (the E3 pooling plus a
  17,942-rank float32 `evd`); the `orion_hek293t` fold sweep (12.9 GiB) ran concurrently with one `fit286.py` process
  to stay inside the 8 h box, so the experiment's combined footprint briefly reached ~24 GiB. These are reported, not
  hidden. System-wide free memory never fell below 73 % and the machine was never under pressure.
* **Stage C leakage.** The Stage B fits exclude the 300 A/B/C panel targets but not h1-test's 85 targets; the overlap
  with each line's fitted bank rows is **82/9,517 (0.86 %)** and **84/17,475 (0.48 %)** (see §5). A covariance fitted
  on 9.5–17.7 k rows moves by O(1/n) per row, so this is small, but it is a leak and it is stated.
* **`kolf2` as a fourth fold was not run** — the time box was spent on the three cancer folds, Stage A and Stage C.
* **E6** (272's cross-fitted operator, its construction B) was **not run**: 272 built its refit arrays for
  `k562_gwps` and `orion_hct116` only (no `orion_hek293t` refit on disk), so a 3/3-fold reading is not available and
  272's published result (+0.0001 … +0.0022, intervals across zero) stands unrepeated here. E4 and E5 — the two
  already-failed baselines the mandate names — were both re-run and both reproduce their published values.
* **KA-E1id** is reported at 3.8e-04 relative rather than 0: the identity is exact in algebra (φ = 1 sets every gain
  to 1) and the residual is the float32 round trip through the Cholesky factor of N̂ + εI. The eigensolver runs in
  float32 (`evd`) because float64 at rank 17.9 k needs ~10 GiB of LAPACK workspace alone; the choice is recorded in
  `est286.py` as a frozen constant.
* The literal `:lit` readings were run only for the primary E1 cells (φ ∈ {0, 0.10, 0.25}); the E2 and control cells
  were read `:nm` only, as frozen.
* **Two reclassifications, both mine and both documented.** `sweepB286.py` recorded `lib263`'s KA-0 as a failure on
  all three folds (§1) and `stageC286.py` recorded KA-C0 / KA-C3 as failures (§5). Neither is an arithmetic failure:
  the first is a fold-cache provenance artefact and the second is 273's own published H-0 reproduced to the digit.
  The raw verdicts are in `work/sweepB286-<fold>.json` and `work/stageC286.json`; the reclassification is stated here
  and, for KA-0, named in `priceB286.py`.
* **Intervals.** Stage B uses a paired **target** bootstrap, 20,000 resamples, stratified by fold and pooled — the
  arithmetic of 265's `boot` / `boot_folds` reimplemented in `analyseB286.py`; Stage A and Stage C call `boot265.boot`
  and `lib280.boot_median*` directly. Per READ-FIRST §3, none of these is quoted in "σ" of 162's 0.001494 (which is a
  leave-one-arm-out RMS of PLOLO→board prediction error in scaled board units, not a per-fold sampling σ).

## 8. What this closes, and what it does not

**Closed.** The construction the third review asked for — shared response directions learned from cross-products of
independent views, with a full-rank residual path — is **measured and does not clear**, on the cleanest held-view
surface on disk (CD4, four donors), on three held contexts, and on the direction-resolving Flex reference. It fails
for a reason that is now quantified rather than assumed: the incumbent's `pds` is bought in directions that an
unbiased signal-covariance estimate scores as noise (80 % of the whitened prediction energy sits where λ⁺(Ŝ) < ν), so
*any* estimator that respects the signal covariance moves prediction mass out of the marginal band where `pds` is
decided. This is the same mechanism 272 found in the fitted basis, and it survives re-estimating the basis — which
was the open question.

**Not closed.** (i) Whether a *target-specific* (rather than global) use of Ŝ helps: everything here applies one
matrix to every row. (ii) Whether the anti-aligned group's genuine gain (+0.042 … +0.109 Δ`pds` at φ = 0.25 … r=100)
can be harvested without paying in the marginal band — that needs a per-target gate, and 263/264 closed every
build-time-observable arbiter, so it would need a new signal. (iii) The **cosine** effect on CD4 at φ = 0.25 (+0.0004,
interval containing zero, and its control at −0.029) is the only cell in this experiment where the learnt directions
clearly beat their controls without losing; a construction that starts from that cell rather than from a `pds` arm
would be a different mandate.

**Nothing here is build-worthy.** No arm goes back to the coordinator for a build.

## 9. Files

| path | what |
|---|---|
| `PREREGISTRATION.md`, `est286.py`, `work/prereg286.sha256` | the frozen family, controls, gate and hashes |
| `stage286.py`, `work/moments/<line>.{json,log}`, `<line>_{S,N,cols}.npy` | the cross-product moments per donor line (float32, ~1.3 GB each on the Orion lines) |
| `fit286.py`, `work/filt/<line>.{npz,json,log}` | the E1 / E2 filters and the filtered 300-target row blocks, the top-100 eigenvectors, the spectra and gains |
| `sweepB286.py`, `work/sweepB286-<fold>.{json,log}`, `pertargetB286-<fold>.npz`, `identityB286-<fold>.npz`, `e3spec-<fold>.npz` | Stage B per fold |
| `analyseB286.py`, `work/analysisB286.json`, `tablesB286.md` | Stage B bootstraps, strata, knee, retention |
| `priceB286.py`, `work/priceB286.{json,md,log}` | the six-member price and the gate |
| `stageA286.py`, `work/stageA286.{json,log}`, `stageA-*.npz` | Stage A (held views / held state) |
| `stageC286.py`, `work/stageC286.{json,log}`, `stageC286.npz`, `stageC-filtered.npz` | Stage C |

The `_S_shuf` / `_N_shuf` control matrices (5.1 GB) were deleted after use to keep the footprint reasonable; they are
deterministic given the frozen derangement seed (`est286.DERANGE_SEED = 20260906`) and regenerate exactly with
`PYTHONDONTWRITEBYTECODE=1 uv run python stage286.py <line>` (~3 min per Orion line, peak 4.7 GiB); their sha256 are
recorded in `work/moments/<line>.json`.
