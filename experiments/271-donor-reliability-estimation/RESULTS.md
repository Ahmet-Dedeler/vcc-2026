# 271 — DOES NOT CLEAR. Donor noise is real and measurable (each donor's row reproduces its own other half at a median residual cosine of +0.047 / +0.028 / +0.028 for K562-GWPS / HCT116 / HEK293T, i.e. a Spearman–Brown full-depth reliability of **0.09 / 0.05 / 0.05**, with 6 % / 32 % / 34 % of targets at or below zero), but no noise-aware combination of the donor rows moves the transfer: the reliability-weighted mean loses to the equal-weight mean (**−0.0037 `pds` [−0.0103, +0.0028]**, pooled over three cancer folds), the per-gene Wiener filters buy `pds` only where the anti-aligned group is perturbed by chance and pay 4–20× on the `nmae` identity, a13's agreement damping re-implemented loses **−0.0069** (its board sign), and the one estimator with an interval above zero — empirical-Bayes row shrinkage toward the consensus, **+0.0014 [+0.0005, +0.0024]** — prices **+0.00050 pessimistic (0.26× the gate)**, retains 31 % leave-top-3-out and is not separable from its shuffled control (+0.0013 [−0.0001, +0.0027]). 0 of 10 cells clear. Two by-products: the reliability of a donor row is a TARGET property shared across lines (the held line's anti-aligned group is the least reliable tier in every donor), and the production operator anti-correlates a donor's replicate halves in-sample (whitened split-half cosine median −0.13…−0.24, 97–100 % ≤ 0) while an out-of-sample line's halves through the same operator read +0.03…+0.07. Measurement only. No build, no submission.

**Date 2026-09-06 (04:17–05:05 local, 01:17–02:05Z), HEAD `08696a7`.** Deterministic (est271 seeds 20260906; bootstrap seed 20260906). `uv run python` throughout, `PYTHONDONTWRITEBYTECODE=1`. **No build, no `vcc`, no `cell-eval2`, no model name, nothing under `submissions/`, `data/harness/` or `data/`, nothing written outside this directory** (`find -newer work/.marker-start` over `experiments/ data/ scripts/ src/ docs/ submissions/` returns nothing outside `experiments/271-*/`). Fold caches are byte copies of 267's = 263's (`work/plolo/SHA256SUMS.txt`: `c13c1cec…`, `9f3e0452…`, `f9a3cc30…` + the three `bctrl-challenge-*.npy`; `diff` against 267's list empty). **The reserved `orion_hek293t` half-B truths (`experiments/267-…/work/reserve/truthB-orion_hek293t.npz`, `84d3879a…`) were not read**; every arm is scored on the fold's CACHE truth (the deep reference, 162's PLOLO). HEK293T's half-B *bank* rows enter only as that donor's own replicate when HEK293T is a donor (K562 and HCT116 folds), as the mandate directs. **Every known-answer row passes: 3 × 2 on the half rows, 3 × 17 on the scoring path (all at 0.0e+00 except two at ≤ 5e-12), the retrodiction of BR249 through the price recipe at 0.0e+00 against 257** (§2). **Memory: peaks 5.9 GiB (half rows), 5.5–5.7 GiB (sweeps), 6.6 GiB (diagnostic); one line / fold per process; 0.1–0.4 min per line, 0.8 min per fold.** Wall clock 48 min from launch to this file (25 min of measurement).

---

## 0. The answer in nine lines

1. **The donor rows are as noisy as 267 said the references are** (§1, table D1): at their own depth (69–96 cells per half, 140–204 cells full), the batch-parity halves of a donor's row agree at a median residual cosine of **+0.047 (K562-GWPS) / +0.028 (HCT116) / +0.028 (HEK293T)**, p10 +0.008 / −0.05 / −0.12, **share ≤ 0: 6 % / 32 % / 34 %**; Spearman–Brown to the full row: **0.091 / 0.054 / 0.054** at the median. The half-depth banks' median row norm is 1.37–1.39× the full-depth one, the moment-based signature of rows that are ~90 % noise energy. Reliability is a weak function of cell count (ρ +0.08…+0.17) and a strong function of the target: the held line's anti-aligned group (262's (e)) is the least reliable tier in every donor (e.g. K562 rows 0.037 / 0.056 / 0.133 for (e) / marginal / saturated of the HCT116 fold), and the two donors' reliabilities correlate across targets (ρ +0.10…+0.17), as do a donor's and the held reference's (+0.10…+0.16).
2. **Reliability knows a little about which donor is closer to the truth, not enough to act on**: pick-the-more-reliable-donor agrees with 263's whitened oracle on **52 / 61 / 58 %** of two-donor targets (break-even 66–70 %), ρ(R₁−R₂, cos₁−cos₂) +0.02 / +0.26 / +0.34.
3. **(a) the reliability-weighted mean (BLUE weights R/(1−R))** turns that into a near-hard pick — the max weight is > 0.9 on 33–40 % of two-donor targets, R ratio > 3 on 51–59 % — and loses: **−0.0037 [−0.0103, +0.0028]** (−0.0066 / −0.0026 / −0.0019; 0/3 folds), between the reversed weights (−0.0194) / shuffled reliabilities (−0.0107) and the mean. The information is real (relw − reversed +0.0157 [+0.0032, +0.0285]) and the estimator still loses on the marginal band (−0.0143 [−0.0227, −0.0059]; 25 up / 40 down across 0.02) while gaining on the anti-aligned group (+0.035) exactly as a random pick does (+0.067). EB-shrunk reliabilities (κ 0.12–0.89 per donor) halve the loss (−0.0007 [−0.0060, +0.0048]); donor-level scalar weights (K562 0.12 vs Orion 0.04–0.06; the 019 / 142 class) −0.0019 [−0.0063, +0.0025]; reliabilities measured after the operator are degenerate (§1.2) and reproduce the mean.
4. **(b) EB row shrinkage toward the consensus** is the mildest re-weighting (max weight median 0.52, p90 0.58) and the only cell above zero: **+0.0014 [+0.0005, +0.0024]**, 3/3 folds (+0.0004 / +0.0025 / +0.0014), `cos_agg` +0.7 %, `nmae` identity +0.0009; positive on every stratum (≤ +0.0035) but **not separable from its own shuffled control (+0.0013 [−0.0001, +0.0027])**, leave-top-3-out retention **31 %** (K562 −136 %: three targets, DTNBP1 / ZNF783 / ZNF496, carry the fold), 0 up / 6 down across the 0.02 knee. Priced **+0.00050 pessimistic / +0.00194 point / +0.00209 optimistic** — 0.26× / 1.02× / 1.10× the shape-only gate.
5. **(c) per-gene Wiener filters**: the row's own half-difference filter keeps 49–63 % of coordinates (55–57 % of energy) and reads +0.0049 [−0.0029, +0.0127] (its shuffled control +0.0013; the difference +0.0036 [−0.0050, +0.0123]) with the gain entirely in the anti-aligned group (+0.069) against a marginal-band loss (−0.0123) and an `nmae` identity cost of **+0.0185**: priced −0.0133 pessimistic. The gene-pooled filter zeroes 51–85 % of genes (pooled signal variance 0), costs −11.6 % `cos_agg` and **+0.101** on the identity: −0.068 pessimistic. 228's wall, again: magnitude reshaping.
6. **(d) a13 re-implemented** (per-gene agreement damping, γ 0.25; graded fraction 0.549 on two Orion donors and 0.199 on K562 + Orion, matching the builder's 0.545 / 0.206–0.208) does what the board said it did: `cos_agg` **+12.6 %**, `pds` **−0.0069 [−0.0153, +0.0015]**, identity `nmae` **+0.040** (norm-matched) — the board's a13 lost −0.0203 with raw `pds` −0.0288 and `nmae` 1.0185. Priced −0.0227 pessimistic (the point set's +0.0035 is the re-shape k_cos 3.76 flattering a cosine bought by magnitude; the literal amplitude-changing form −0.0088 / +0.0050). Every estimator beats (d) on `pds` (ebrow +0.0083 [+0.0001, +0.0167]); none beats it by the mandate's margin AND clears.
7. **Pricing** (261's recipe imported; retrodiction of BR249 at 0.0e+00 against 257: +0.00472 vs the board's +0.00426): **0 / 10 cells** at the pessimistic gate, 2 / 10 with `pds` 3/3 and the interval above 0 (ebrow, ebrow_eb), 5 / 10 with retention ≥ 50 % (all by being negative), **0 / 10 all conditions**.
8. **In-sample whitening** (§1.2, `diagwht271.py`): the fold operator is fitted on the two donors' full-depth rows; their half rows whitened through it are ANTI-correlated (median −0.13…−0.24, 97–100 % ≤ 0; null pairs 0.000 ± 0.01), while the held line's half rows through the same operator read +0.074 / +0.041 / +0.029 (48 / 32 / 48 % ≤ 0; null sd 0.09–0.64). The amplified directions are the ones in which each fitted row's two halves cancel — the operator has learnt the noise realisation of its own rows. This is a diagnosis of the incumbent operator, not a lever tested here (§7).
9. **The 400-cell question** (§5): on these lines' laws a 400-cell reference reaches a target-specific reliability of only **0.11–0.20** (p = 1) / **0.07–0.13** (p = 0.70); Flex chemistry measured **0.48** at 400 on H1 (159). The per-target observed-cosine ceiling √ρ is **0.23–0.30 on PLOLO's deep truths, 0.27–0.44 for an Orion-like 400-cell reference, 0.69 for an H1-like one**. PLOLO's 0.02 knee is a signal cosine of ≈ 0.07–0.09, which reads 0.05–0.06 on an H1-like A/B/C reference — saturated; the board's marginal band (signal cosine ≈ 0.04) reads ≈ 0.01 on PLOLO, inside the reference's noise. The two marginal bands are different target sets.

**One line: DOES NOT CLEAR — best pessimistic +0.00050 (ebrow; point +0.00194, optimistic +0.00209), `pds` +0.0014 [+0.0005, +0.0024] 3/3 but shuffled control +0.0013 [−0.0001, +0.0027] and retention 31 %; (a) −0.0037, (c) −0.013 / −0.068 pessimistic, (d) −0.0227; 0/10 cells.** The construction "shrink each donor row by its own measured sampling reliability before the mean" is closed at this scope: the reliabilities are real, weakly informative about the donor choice, and their correction does not move the transfer cosine on the marginal band.

---

## 1. Donor-row reliability at the donor's own depth (item 1; `halfrows271.py`, `sweep271.py`, `analyse271.py` table D1)

### 1.1 Construction and known answers

267's exact integer pseudobulk sums per (perturbation, batch-parity half) are re-written as half-banks in the donor-table format (truth267's arithmetic verbatim; 267 deleted the `delta.npy` files) and each half is prepared by the instrument's own path (`E062.prepare_main_source`: the half-bank's own column mean over all 9,786–17,977 perturbations, × `EFFECT_SCALE / raw_median`), then converted to the production row's units by the real bank's raw median. Row = the production un-whitened donor row on the official axis (the K562 row un-blended here; lib263.stage applies the PerturbMatch 0.5 blend downstream exactly as the fold cache was built). Reliability r = cos(row_A, row_B) on the line's mapped columns minus the panel; R = 2r/(1+r) (the full row is the union of the halves); whitened reliability = the same after `plolo.Operator.apply` with the fold's held operator, on the operator columns.

| row | required | K562-GWPS | HCT116 | HEK293T |
|---|---|---|---|---|
| **KA-A** the union-of-halves bank through `prepare_main_source` == the production donor row | float32 | **0.0e+00**, 0 duplicate symbols, raw-median ratio 1.000000000 | **0.00 float32 ulp** after the duplicate-column scalar (4 duplicate symbols enter the real bank's raw median: ratio 0.999875460; min row cosine 1.000000000) | **0.00 ulp** (ratio 0.999923671, 4 duplicates, min row cos 1.000000000) |
| **KA-B** split-half reliability of the half rows == 267 `rel\|batch\|cos_resid` | 1e-5 | **2.4e-09** (269) | **1.6e-09** (290) | **1.7e-09** (294) |
| **KA-C** (per fold) the half-builder's real row == `lib263.stage`'s production row, float32; availability identical | 0 | 0.0e+00 off the 7,075 PerturbMatch columns (blended there, as built) | 0.0e+00 | 0.0e+00 |
| **KA-0 / R1** (per fold) the staged per-donor rows rebuild the cache `pre`; the explicit mean == `cd4svd5` | 0 | 0.0e+00 / bit-identical per-target `pds` | same | same |
| **eqR** equal reliabilities through `weights()` + `combine(norm_match)` == `cd4svd5` on `cos_agg`/macro/`pds` | 1e-6 | 5.0e-12 | 2.9e-12 | 1.9e-12 |
| `cd4svd5` == 263; `STACK` == 257 (incl. `nmae_wald89`); `oracle_wht:all:nm` picks and score == 263; `random:all:s0:nm` == 263 | 1e-9 | all 0.0e+00 | all 0.0e+00 | all 0.0e+00 |

### 1.2 The table (batch-parity halves; SB = Spearman–Brown to the full row; the by-group columns use the HELD fold's 262 group and 267's reference-reliability tercile)

| fold (held) | donor | n | cells / half (median) | cells full | raw mean | raw p10 | raw p25 | **raw median** | raw p75 | raw p90 | **share ≤ 0** | share ≥ 0.05 | **SB median** | SB p10 / p90 | wht median | wht share ≤ 0 | ρ(r, cells) | raw by 262 group e / marg / sat | raw by ref-reliability T1 / T2 / T3 | EB κ (r̄, s_null) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|---|
| `k562_gwps` | `orion_hct116` | 260 | 69 | 152 | +0.035 | −0.048 | −0.012 | **+0.027** | +0.069 | +0.119 | **0.315** | 0.335 | **+0.053** | −0.100 / +0.212 | −0.141 | 0.985 | +0.13 | +0.027 / +0.020 / +0.077 | +0.027 / +0.028 / +0.050 | 0.35 (+0.035, 0.064) |
| `k562_gwps` | `orion_hek293t` | 263 | 94 | 201 | +0.042 | −0.119 | −0.022 | **+0.028** | +0.096 | +0.211 | **0.335** | 0.407 | **+0.054** | −0.269 / +0.348 | −0.141 | 0.992 | +0.15 | +0.051 / +0.028 / +0.070 | +0.026 / +0.046 / +0.054 | 0.12 (+0.042, 0.135) |
| `orion_hct116` | `k562_gwps` | 260 | 80 | 169 | +0.065 | +0.008 | +0.023 | **+0.048** | +0.084 | +0.162 | **0.058** | 0.477 | **+0.091** | +0.015 / +0.279 | −0.142 | 0.973 | +0.16 | +0.037 / +0.056 / +0.133 | +0.056 / +0.061 / +0.079 | 0.87 (+0.065, 0.022) |
| `orion_hct116` | `orion_hek293t` | 289 | 96 | 204 | +0.041 | −0.123 | −0.021 | **+0.028** | +0.097 | +0.212 | **0.343** | 0.412 | **+0.054** | −0.282 / +0.350 | −0.243 | 1.000 | +0.13 | +0.032 / +0.028 / +0.110 | +0.023 / +0.031 / +0.069 | 0.13 (+0.041, 0.135) |
| `orion_hek293t` | `k562_gwps` | 263 | 80 | 169 | +0.065 | +0.008 | +0.023 | **+0.048** | +0.084 | +0.162 | **0.057** | 0.471 | **+0.091** | +0.015 / +0.278 | −0.128 | 0.973 | +0.17 | +0.041 / +0.060 / +0.133 | +0.057 / +0.066 / +0.072 | 0.89 (+0.065, 0.021) |
| `orion_hek293t` | `orion_hct116` | 289 | 70 | 153 | +0.036 | −0.050 | −0.012 | **+0.028** | +0.069 | +0.119 | **0.318** | 0.339 | **+0.054** | −0.104 / +0.213 | −0.241 | 0.997 | +0.08 | +0.019 / +0.026 / +0.124 | +0.024 / +0.022 / +0.062 | 0.43 (+0.036, 0.063) |

Random-cell halves (267's `rel|rand|cos_resid`) read +0.045 / +0.029 / **+0.071** with 6 / 30 / 15 % ≤ 0 — HEK293T's random halves overstate by 2.5× where its 223 batches matter (267 §3 (v)); the batch split is used throughout. The moment-based reading agrees: the half-depth banks' median row norms are ×1.382 (K562), ×1.387–1.393 (HCT116), ×1.370–1.390 (HEK293T) of the full-depth ones, which under ‖half‖² = S + 2N, ‖full‖² = S + N gives S/(S+N) = 0.09 / 0.07 / 0.07 for the full rows (159 §1.4's estimator).

Readings. (i) **The donor rows are the references' noise, seen from the other side**: the K562-GWPS rows (169 cells) are ~9 % signal energy, the Orion rows (150–200 cells) ~5 %; a third of the Orion targets have no measurable direction in either half. (ii) **Reliability is a target property.** Within a donor, r is highest on the HELD line's saturated group and lowest on its anti-aligned group (e) — the group whose truth 267 found unmeasurable — and the two donors' r correlate (ρ +0.14 / +0.10 / +0.17); the weak-response tier is weak in every line (262 §4.3's CD4 DE count said the same). Cell count explains little (ρ ≤ +0.17). (iii) **The sampling error of a single split-half cosine is comparable to the spread across targets** for the Orion donors: the derangement null (cos(A_t, B_t')) has sd 0.063–0.135 against a between-target sd of 0.08–0.14, so the empirical-Bayes shrinkage factor κ is 0.12–0.43 for the Orion rows (0.87–0.89 for K562-GWPS, whose null sd is 0.02). A per-target reliability of an Orion row is mostly noise about noise. (iv) **After the fold's operator the halves ANTI-correlate** (median −0.13…−0.24, 97–100 % ≤ 0) — the diagnostic below.

**In-sample whitening (`diagwht271.py`, `work/diagwht271.json`).** Each fold's operator is fitted on the two donors' full-depth rows. Through it, the donors' half rows read cos(WA, WB) median −0.141 / −0.140 (K562 fold), −0.142 / −0.244 (HCT116 fold), −0.127 / −0.241 (HEK293T fold), 97–100 % ≤ 0, with the deranged null at 0.000 ± 0.008–0.012 (no common component; whitened half-to-full cosine 0.62–0.66 against 0.72 raw). The HELD line's half rows through the same operator — rows the fit never saw — read **+0.074 (K562; raw +0.047), +0.041 (HCT116; raw +0.028), +0.029 (HEK293T; raw +0.028)** with 48 / 32 / 48 % ≤ 0 and a deranged null of sd 0.64 / 0.09 / 0.41. Mechanism: the operator amplifies the low-eigenvalue directions of the fitted covariance; with n/p ≈ 1.6–2 those are the directions in which the fitted rows' own noise realisations are smallest, i.e. in which each row's two halves cancel (A ≈ −B relative to their mean), so the halves anti-correlate there by selection, and any out-of-sample vector projected onto them is dominated by amplified noise (the held line's null sd). Consequence for this mandate: a reliability measured after the operator is meaningless in-sample (every R clips to R_MIN and `relw_wht` reproduces the mean), and the un-whitened reliability is the only usable one. Consequence for the incumbent: its whitened donor rows are, in the amplified directions, clean by construction (their noise was selected small) while any independent replicate of the same line — and the held truth — is not; this is one mechanism behind whitening's −32 % cosine for +8 % `pds` (AGENTS §4.2). Whether cross-fitted whitening (fit on half A, whiten half B, and the mirror, averaged) helps or hurts `pds` is a measurement this experiment did not run (§7).

## 2. The estimators (item 2; `est271.py`, sha256 `fe1da9917a2a…`, frozen before the first fold; formulas in its docstring)

All combinations happen BEFORE whitening — the one place a per-target, per-donor weight enters production (`scripts/make_submission_pm_cd4.py` `build_v9_pre_whiten_effects` → `MultiDonorModel.effect(target, alpha, beta, weights=…)`, line 1380; `donor_weights()` at 1342 is the global `{name: 1.0}` dict) — and every combined row is norm-matched to the cached mean row (263's `:nm`, the direction-only reading), so the arm is shape-only and the residual sizing and the global preserve see the incumbent's row norms. Downstream is `cd4svd5` verbatim (lib263: `plolo.Operator.apply`, CD4 deflated k = 5 re-solved at f = 0.34, the H1-pair union at g = 0.28, one add, one preserve). With two donors every per-target scalar rule is a convex weight w on donor 1:

| cell | rule | what it is |
|---|---|---|
| (e) `cd4svd5` | w = ½ | the incumbent (KA 0.0e+00) |
| (a) `relw` | w_j ∝ R_j/(1−R_j), R = 2r/(1+r) clipped to [0.005, 0.99] | the BLUE of two measurements with noise variances Var(s)(1−R_j)/R_j |
| `relw_eb` | as (a) with r̃_j = r̄_j + κ_j(r_j − r̄_j), κ_j = τ²/(τ² + s²), s² = Var of the deranged-pair cosine | knob-free EB shrinkage of the reliabilities themselves; κ reported per donor per fold |
| `relw_glob` | w_j ∝ R̄_j/(1−R̄_j), the donor's mean R | the scalar donor-weight class (019 / 142 / `VCC_DONOR_WEIGHTS`); K562 R̄ 0.117, Orion 0.041–0.058 → 0.67 : 0.33 |
| `relw_wht` | as (a) with r measured after the operator | degenerate in-sample (§1.2) |
| (b) `ebrow` | row̃_j = m + R_j(row_j − m), then the mean ⇔ w_j = (1 + R_j − R̄)/J | EB shrinkage of each row toward the consensus by its noise share |
| `ebrow_eb` | with r̃ | |
| (c) `wiener_own` | per gene w = max(0, 1 − d²/(4x²)), d = A − B, x the production row; mean; norm-matched | the row's own between-half variance as a diagonal filter |
| `wiener_gene` | σ²_tg = c_g/n_t with c_g = mean_t[n_t d²_tg/4]; s²_g = max(mean_t x² − mean_t σ², 0); w = s²/(s² + σ²) | gene-pooled noise scaled by cell count |
| (d) `a13:nm` / `a13:lit` | α_g = exp(−0.25 · midrank_g/2000) of the cross-donor agreement key (geometric mean of \|row_j(g)\| where all carrying donors agree in sign and are nonzero, else 0; descending midrank over the official axis; rows with < 2 distinct key values untouched), applied to the FINAL composed effect (where production applies it, at render); `:nm` re-normalises each row | the baseline (a13, board −0.0203); graded fraction 0.549 (two Orion donors) / 0.199 (K562 + Orion), the builder's 0.545 / 0.206–0.208 |
| controls | `relw:shuf`, `ebrow:shuf`, `wiener_own:shuf` (the donor's reliability vector / half-difference rows permuted across its targets, seed 20260906); `relw:rev` (w ∝ (1−R)/R); 263's `random:all:s0:nm`; 263's `oracle_wht:all:nm` (the ceiling) | |

Screened: 10 estimator cells + 6 control rows + 4 known-answer rows, on three folds — one frozen list, no dose, no second pass.

## 3. Measurement (item 3; `analyse271.py` → `work/tables271.md`, `work/analysis271.json`)

### 3.1 Every arm vs `cd4svd5` (Δ`pds`, a23 tag, scored on the cache truth; paired target bootstrap, 20,000 draws, stratified by fold; `cd4svd5` = 0.8154 / 0.8017 / 0.7810)

| arm | Δ`pds` K562 / HCT116 / HEK293T | **pooled [95 % CI]** | P(≤ 0) | Δ`cos_agg` % | Δ macro cos % | delivered ratio | Δ`nmae` identity (wald89, 3 folds concat) | Δ top-89 precision | rows changed |
|---|---|---|---:|---:|---:|---:|---:|---:|---|
| (a) `relw` | −0.0066 / −0.0026 / −0.0019 | **−0.0037 [−0.0103, +0.0028]** | 0.866 | −3.23 | −4.01 | 1.011 | +0.0049 | −0.0083 | 259 / 259 / 259 |
| `relw_eb` | −0.0047 / +0.0005 / +0.0022 | −0.0007 [−0.0060, +0.0048] | 0.596 | −0.82 | −2.33 | 1.012 | +0.0037 | −0.0051 | 259 / 259 / 259 |
| `relw_glob` | −0.0020 / −0.0063 / +0.0026 | −0.0019 [−0.0063, +0.0025] | 0.802 | −3.30 | −3.64 | 1.039 | +0.0047 | −0.0052 | 259 / 259 / 259 |
| `relw_wht` | +0.0000 / −0.0004 / +0.0005 | +0.0000 [−0.0004, +0.0005] | 0.409 | +0.14 | −0.25 | 1.000 | +0.0003 | −0.0000 | 255 / 257 / 259 |
| (b) `ebrow` | +0.0004 / +0.0025 / +0.0014 | **+0.0014 [+0.0005, +0.0024]** | 0.001 | +0.73 | +0.45 | 1.000 | +0.0009 | +0.0007 | 259 / 259 / 259 |
| `ebrow_eb` | +0.0001 / +0.0009 / +0.0014 | +0.0008 [+0.0004, +0.0013] | 0.000 | +0.55 | +0.18 | 1.003 | +0.0012 | +0.0008 | 259 / 259 / 259 |
| (c) `wiener_own` | +0.0024 / −0.0002 / +0.0125 | +0.0049 [−0.0029, +0.0127] | 0.105 | −3.04 | −0.48 | 1.001 | **+0.0185** | +0.0025 | 264 / 290 / 293 |
| `wiener_gene` | +0.0064 / +0.0002 / +0.0029 | +0.0031 [−0.0082, +0.0143] | 0.291 | **−11.58** | −8.12 | 1.004 | **+0.1013** | −0.0113 | 264 / 290 / 293 |
| (d) `a13:nm` | −0.0009 / −0.0160 / −0.0039 | **−0.0069 [−0.0153, +0.0015]** | 0.947 | **+12.62** | +14.23 | 1.081 | **+0.0397** | −0.0036 | (final rows) |
| `a13:lit` | +0.0000 / −0.0153 / −0.0030 | −0.0061 [−0.0145, +0.0024] | 0.921 | +10.32 | +15.03 | **0.721** | +0.0021 | −0.0038 | (final rows) |
| `relw:shuf` | −0.0097 / −0.0160 / −0.0063 | −0.0107 [−0.0178, −0.0035] | 0.999 | −8.68 | −8.59 | 1.009 | +0.0024 | −0.0100 | |
| `relw:rev` | −0.0144 / −0.0230 / −0.0209 | −0.0194 [−0.0270, −0.0120] | 1.000 | −11.81 | −10.14 | 0.999 | +0.0013 | −0.0119 | |
| `ebrow:shuf` | −0.0001 / −0.0004 / +0.0008 | +0.0001 [−0.0010, +0.0012] | 0.406 | −0.19 | −0.18 | 1.005 | +0.0013 | −0.0003 | |
| `wiener_own:shuf` | +0.0021 / +0.0022 / −0.0003 | +0.0013 [−0.0054, +0.0082] | 0.352 | −1.21 | −0.52 | 0.984 | +0.0193 | +0.0013 | |
| `random:all:s0:nm` (263) | −0.0248 / −0.0233 / −0.0218 | −0.0233 [−0.0337, −0.0132] | 1.000 | −16.24 | −15.77 | 1.010 | +0.0047 | −0.0177 | |
| `oracle_wht:all:nm` (263, the ceiling) | +0.0490 / +0.0468 / +0.0417 | +0.0459 [+0.0382, +0.0539] | 0.000 | +2.33 | +2.31 | 1.000 | +0.0052 | −0.0064 | |
| `STACK` (h1union; the retrodiction row) | −0.0100 / −0.0094 / −0.0174 | −0.0123 [−0.0190, −0.0055] | 1.000 | +0.50 | +0.84 | 1.005 | +0.0014 | +0.0011 | |

### 3.2 Each estimator against its controls (paired per-target difference, pooled)

| arm | vs its shuffled control | vs `relw:rev` | vs (d) `a13:nm` | vs `random:all:s0:nm` | vs `relw_glob` |
|---|---|---|---|---|---|
| `relw` | +0.0070 [−0.0018, +0.0159] | +0.0157 [+0.0032, +0.0285] | +0.0032 [−0.0072, +0.0137] | +0.0196 [+0.0084, +0.0308] | −0.0018 [−0.0088, +0.0052] |
| `relw_eb` | +0.0100 [+0.0017, +0.0184] | +0.0187 [+0.0077, +0.0301] | +0.0062 [−0.0037, +0.0161] | +0.0226 [+0.0120, +0.0333] | +0.0012 [−0.0053, +0.0077] |
| `relw_glob` | +0.0088 [+0.0015, +0.0163] | +0.0175 [+0.0088, +0.0262] | +0.0050 [−0.0051, +0.0152] | +0.0214 [+0.0111, +0.0316] | — |
| `ebrow` | **+0.0013 [−0.0001, +0.0027]** | +0.0208 [+0.0127, +0.0291] | **+0.0083 [+0.0001, +0.0167]** | +0.0247 [+0.0142, +0.0350] | +0.0033 [−0.0012, +0.0079] |
| `ebrow_eb` | +0.0007 [−0.0004, +0.0019] | +0.0202 [+0.0123, +0.0282] | +0.0077 [−0.0006, +0.0161] | +0.0241 [+0.0140, +0.0343] | +0.0027 [−0.0016, +0.0070] |
| `wiener_own` | +0.0036 [−0.0050, +0.0123] | +0.0243 [+0.0136, +0.0352] | +0.0118 [+0.0006, +0.0232] | +0.0282 [+0.0163, +0.0401] | +0.0068 [−0.0016, +0.0153] |
| `wiener_gene` | +0.0018 [−0.0091, +0.0129] | +0.0226 [+0.0098, +0.0351] | +0.0101 [−0.0035, +0.0236] | +0.0264 [+0.0133, +0.0395] | +0.0050 [−0.0060, +0.0162] |
| `a13:nm` | +0.0038 [−0.0071, +0.0151] | +0.0125 [+0.0007, +0.0242] | — | +0.0164 [+0.0032, +0.0295] | −0.0050 [−0.0153, +0.0052] |

The ordering reversed < shuffled < reliability-weighted < mean (−0.019 < −0.011 < −0.004 < 0) says the reliabilities carry information about the donor choice — and that any per-target re-weighting keyed on them costs more than it earns, because the weights are near-hard picks on targets where the reference cannot adjudicate (§3.3).

### 3.3 By 262 group and by the reference's reliability tercile (pooled; mean Δ`pds` per target [95 % CI]; n = 149 / 564 / 140 and 285 / 283 / 285)

| arm | (e) base·truth ≤ 0 | marginal band (0, 0.05] | saturated > 0.05 | ref-reliability T1 (low) | T2 | T3 (high) |
|---|---|---|---|---|---|---|
| `relw` | +0.0352 [+0.0147, +0.0567] | **−0.0152 [−0.0236, −0.0070]** | −0.0006 | −0.0025 [−0.0144, +0.0093] | −0.0084 [−0.0209, +0.0040] | −0.0002 [−0.0102, +0.0097] |
| `relw_eb` | +0.0276 [+0.0097, +0.0465] | −0.0084 [−0.0152, −0.0018] | −0.0000 | −0.0009 | −0.0030 | +0.0019 [−0.0061, +0.0101] |
| `relw_glob` | +0.0174 [+0.0034, +0.0318] | −0.0073 [−0.0124, −0.0022] | −0.0000 | +0.0046 | −0.0033 | −0.0070 [−0.0148, +0.0005] |
| `ebrow` | +0.0035 [+0.0002, +0.0070] | +0.0012 [+0.0002, +0.0023] | +0.0000 | +0.0005 [−0.0010, +0.0021] | +0.0015 [−0.0002, +0.0034] | +0.0022 [+0.0008, +0.0038] |
| `ebrow_eb` | +0.0022 [+0.0005, +0.0043] | +0.0006 [+0.0002, +0.0011] | +0.0001 | +0.0006 | +0.0008 | +0.0011 [+0.0003, +0.0021] |
| `wiener_own` | **+0.0689 [+0.0434, +0.0949]** | −0.0123 [−0.0217, −0.0032] | +0.0002 | −0.0075 | +0.0134 [+0.0009, +0.0260] | +0.0089 [−0.0010, +0.0191] |
| `wiener_gene` | +0.1164 [+0.0809, +0.1522] | −0.0251 [−0.0381, −0.0124] | −0.0019 | +0.0109 | +0.0116 | −0.0129 |
| `a13:nm` | −0.0012 [−0.0278, +0.0262] | −0.0094 [−0.0197, +0.0010] | −0.0029 | −0.0120 | −0.0054 | −0.0033 |
| `relw:shuf` | +0.0229 [+0.0025, +0.0439] | −0.0216 [−0.0311, −0.0125] | −0.0020 | −0.0062 | −0.0145 | −0.0114 |
| `random:all:s0:nm` | +0.0672 [+0.0362, +0.0985] | −0.0527 [−0.0654, −0.0403] | −0.0053 | −0.0353 | −0.0169 | −0.0176 |
| `oracle_wht:all:nm` | +0.1718 | +0.0242 | −0.0001 | +0.0638 | +0.0439 | +0.0299 |

267's mechanism, once more: on the anti-aligned group ANY change to the row helps (random +0.067, the own-variance filter +0.069, the reliability pick +0.035) and on the marginal band every non-trivial re-weighting loses (−0.008…−0.053). The EB row shrinkage is the exception only because it barely moves the rows.

### 3.4 The marginal band (cd4svd5's own cosine in (0, 0.05], n = 527) — where +0.01 cosine buys +0.095 `pds` (262 §4.2)

| arm | mean Δ own cosine on the band | up across 0.02 | down across 0.02 | net | mean Δ`pds` on the band [95 % CI] | Δ own cosine, all targets |
|---|---:|---:|---:|---:|---|---:|
| `relw` | −0.00089 | 25 | 40 | −2.8 % | −0.0143 [−0.0227, −0.0059] | −0.00137 |
| `relw_eb` | −0.00045 | 13 | 21 | −1.5 % | −0.0098 [−0.0166, −0.0033] | −0.00076 |
| `relw_glob` | −0.00096 | 14 | 33 | −3.6 % | −0.0047 [−0.0102, +0.0005] | −0.00112 |
| `ebrow` | +0.00012 | 0 | 6 | −1.1 % | +0.0009 [−0.0001, +0.0020] | +0.00013 |
| `ebrow_eb` | +0.00004 | 1 | 4 | −0.6 % | +0.0005 [+0.0001, +0.0010] | +0.00005 |
| `wiener_own` | −0.00052 | 25 | 42 | −3.2 % | −0.0126 [−0.0224, −0.0031] | −0.00018 |
| `wiener_gene` | −0.00185 | 41 | 77 | −6.8 % | −0.0263 [−0.0402, −0.0125] | −0.00247 |
| `a13:nm` | **+0.00505** | 70 | 21 | +9.3 % | **−0.0229 [−0.0336, −0.0124]** | +0.00442 |
| `relw:shuf` | −0.00208 | 15 | 47 | −6.1 % | −0.0206 [−0.0305, −0.0110] | −0.00280 |
| `random:all:s0:nm` | −0.00366 | 17 | 70 | −10.1 % | −0.0560 [−0.0689, −0.0436] | −0.00507 |
| `oracle_wht:all:nm` | +0.00173 | 40 | 25 | +2.8 % | +0.0260 [+0.0187, +0.0335] | +0.00067 |

No estimator moves the band's own cosine by more than +0.0001; a13's +0.005 on the band (70 up) comes with −0.023 `pds` there, because the damping moves every row toward the same agreement-weighted shape and the distractors move with it (the cosine gain is not discrimination — AGENTS §4.2's lesson in one row).

### 3.5 Leave-top-3-out retention (261's statistic) and the delivered amplitude

| arm | full Δ`pds` | leave-top-3-out | retained | per fold | top-3 targets |
|---|---:|---:|---:|---|---|
| `ebrow` | +0.00140 | +0.00043 | **31 %** | −136 % / +50 % / +39 % | K562: DTNBP1 +0.085, ZNF783, ZNF496; HCT116: TMF1 +0.152, TIMM50, PDCD4; HEK293T: ZNF219 +0.111, ZNF786, ZNF496 |
| `ebrow_eb` | +0.00082 | +0.00035 | 43 % | −108 % / +57 % / +49 % | the same names |
| `wiener_own` | +0.00491 | +0.00033 | 7 % | −137 % / +2618 % / +71 % | |
| `wiener_gene` | +0.00314 | −0.00321 | −102 % | | |
| `oracle_wht:all:nm` | +0.04586 | +0.04000 | 87 % | 88 / 88 / 85 % | |

Delivered-amplitude ratios (frame A, `lib257.Identity`): every norm-matched cell 1.000–1.012 except `relw_glob` 1.039 and `a13:nm` 1.081 (the damping re-shapes toward the gated genes); `a13:lit` 0.721 (the literal arm removes 28 % of the delivered amplitude — the board's a13 carried it with `--effect-scale`). The `nmae` identity: the (a)/(b) cells +0.0003…+0.0049, the Wiener filters +0.0185 / +0.1013, a13 +0.0397 (`:nm`) / +0.0021 (`:lit`, the amplitude drop hides it: the board realised +0.074 raw against a12's 0.945 → 1.0185, which the identity under-calls here as it did in 098 §"one arm of two").

## 4. Price (item 4; `price271.py` → `work/price271.{json,log,md}`; 261's `price` / `factors` / `scenarios` / `retention` imported, 257's constants)

Retrodiction first: `cd4svd5` vs `STACK` through the recipe at h1union's terms reads Δ`pds` +0.012263, `cos_agg` −0.4617 %, identity −0.00138 — **0.0e+00 against 257's R1** — and prices `pds` +0.00552, `mse` −0.00083, DE-three −0.00035, `nmae` +0.00038 = **+0.00472** (257 +0.00472; BR249 board +0.0042634; all three bands pass; this file's `price()` == `price257.price()` at 0.0e+00).

Scenarios: PESSIMISTIC / POINT / OPTIMISTIC = per leg the band end that lowers / the set that retrodicts BR249 / raises the total for that cell's sign: k_pds {0.85, 1.22}, k_cos {1.36, 3.76}, DE-three {0.0425, 0.24 × 0.776}, identity ×1.5 on a loss and ×0 on a gain (pessimistic). **The board class of a base RE-ESTIMATION is unmeasured**; the re-dose / re-shape bands are the pessimistic reading for the sign of each leg, as the mandate instructs, and a first board point of this class would set its own constants.

| cell | Δ`pds` (3 folds) | pooled 95 % CI | Δ`cos_agg` % | delivered | Δ`nmae` id | `pds` | `mse` | DE-three | `nmae` | **PESSIMISTIC** | POINT | OPTIMISTIC | LT3O | gate / pds / ret / ctl | cleared |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `relw` | −0.00370 (0/3) | [−0.0103, +0.0028] | −3.23 | +1.07 % | +0.0049 | −0.00167 | −0.00439 | −0.00238 | −0.00136 | **−0.01117** | −0.00980 | −0.00422 | — | F / F / — / F | no |
| `relw_eb` | −0.00067 (2/3) | [−0.0060, +0.0048] | −0.82 | +1.19 % | +0.0037 | −0.00030 | −0.00111 | −0.00060 | −0.00103 | **−0.00372** | −0.00304 | −0.00166 | — | F / F / — / F | no |
| `relw_glob` | −0.00189 (1/3) | [−0.0063, +0.0025] | −3.30 | +3.88 % | +0.0047 | −0.00085 | −0.00704 | −0.00244 | −0.00130 | **−0.01298** | −0.01163 | −0.00610 | — | F / F / — / F | no |
| `relw_wht` | +0.00005 (2/3) | [−0.0004, +0.0005] | +0.14 | +0.00 % | +0.0003 | +0.00002 | +0.00019 | +0.00010 | −0.00009 | **−0.00004** | +0.00023 | +0.00026 | — | F / F / F / F | no |
| **`ebrow`** | **+0.00140 (3/3)** | **[+0.0005, +0.0024]** | +0.73 | +0.04 % | +0.0009 | +0.00063 | +0.00101 | +0.00054 | −0.00023 | **+0.00050** | +0.00194 | +0.00209 | 31 % | F / T / F / F | **no** |
| `ebrow_eb` | +0.00082 (3/3) | [+0.0004, +0.0013] | +0.55 | +0.34 % | +0.0012 | +0.00037 | +0.00087 | +0.00040 | −0.00034 | **+0.00017** | +0.00131 | +0.00142 | 43 % | F / T / F / F | no |
| `wiener_own` | +0.00491 (2/3) | [−0.0029, +0.0127] | −3.04 | +0.13 % | +0.0185 | +0.00221 | −0.00422 | −0.00225 | −0.00514 | **−0.01328** | −0.00939 | −0.00464 | 7 % | F / F / F / F | no |
| `wiener_gene` | +0.00314 (3/3) | [−0.0082, +0.0143] | −11.58 | +0.44 % | +0.1013 | +0.00142 | −0.01597 | −0.00855 | −0.02814 | **−0.06821** | −0.05124 | −0.03311 | −102 % | F / F / F / F | no |
| `a13:nm` (d) | −0.00691 (0/3) | [−0.0153, +0.0015] | +12.62 | +8.10 % | +0.0397 | −0.00311 | +0.00828 | +0.00932 | −0.01102 | **−0.02267** | +0.00347 | +0.00711 | — | F / F / — / F | no |
| `a13:lit` (d) | −0.00609 (1/3) | [−0.0145, +0.0024] | +10.32 | −27.95 % | +0.0021 | −0.00274 | +0.00073 | +0.00762 | −0.00058 | **−0.00883** | +0.00503 | +0.00806 | — | F / F / — / F | no |
| `relw:shuf` | −0.01067 | [−0.0178, −0.0035] | −8.68 | | +0.0024 | | | | | −0.02590 | −0.02372 | −0.00863 | | | control |
| `relw:rev` | −0.01941 | [−0.0270, −0.0120] | −11.81 | | +0.0013 | | | | | −0.03689 | −0.03419 | −0.01308 | | | control |
| `ebrow:shuf` | +0.00011 | [−0.0010, +0.0012] | −0.19 | | +0.0013 | | | | | −0.00083 | −0.00058 | −0.00028 | | | control |
| `wiener_own:shuf` | +0.00132 | [−0.0054, +0.0082] | −1.21 | | +0.0193 | | | | | −0.01098 | −0.00786 | −0.00599 | | | control |
| `random:all:s0:nm` | −0.02329 | [−0.0337, −0.0132] | −16.24 | | +0.0047 | | | | | −0.05043 | −0.04631 | −0.01762 | | | control |
| `oracle_wht:all:nm` | +0.04586 | [+0.0382, +0.0539] | +2.33 | | +0.0052 | | | | | +0.01353 | +0.02414 | +0.02464 | 87 % | | ceiling |

(legs at POINT; "gate" = pessimistic ≥ +0.0019, "pds" = 3/3 folds and the pooled interval above 0, "ret" = leave-top-3-out ≥ 50 %, "ctl" = beats its shuffled control AND `a13:nm` by more than their half-widths with the paired interval above 0.)

**THE SCREEN: 10 cells. Pessimistic ≥ +0.0019: 0 / 10. `pds` 3/3 with the interval above 0: 2 / 10 (`ebrow`, `ebrow_eb`). Retention ≥ 50 %: 0 of those two. All conditions: 0 / 10.** Best pessimistic `ebrow` +0.00050 (0.26× gate; point +0.00194 = 1.02×; optimistic +0.00209); its `pds` leg is +0.0006 overall and its cosine leg +0.0010 at k_cos 3.76 — a re-shape constant applied to a +0.7 % cosine move that its shuffled control shares (−0.2 %). The point-scenario positives of `a13:lit` (+0.0050) and `a13:nm` (+0.0035) are k_cos 3.76 on a +10–13 % cosine gain that the board already priced at −0.0203: a13's cosine is a magnitude reshaping (228 §8.1d), and the identity is known to under-call its `nmae` (LOG 08-29 §1; 098 §"one arm of two"); the pessimistic reading (−0.0088 / −0.0227) is the board-consistent one.

## 5. The 400-cell question (item 5; `analyse271.py` table Q; 159 §1.3's laws on 267's batch-half numbers)

| line | cells / half (median) | r at half depth | SB → full donor depth (cells) | power law p = 1.0 → 400 cells | p = 0.70 → 400 | √ρ at donor depth | √ρ at 400 (p = 1.0 / 0.70) |
|---|---:|---:|---:|---:|---:|---:|---:|
| K562-GWPS (10x, Replogle) | 81 | +0.047 | 0.090 (162) | 0.196 | 0.131 | 0.30 | 0.44 / 0.36 |
| HCT116 (X-Atlas/Orion) | 70 | +0.028 | 0.054 (140) | 0.139 | 0.088 | 0.23 | 0.37 / 0.30 |
| HEK293T (Orion) | 95 | +0.028 | 0.055 (190) | 0.109 | 0.074 | 0.23 | 0.33 / 0.27 |
| H1 (Flex, 159 §1.3, MEASURED at 400) | — | — | — | **0.481** | — | — | **0.69** |

Law: 1/ρ − 1 = a · m^(−p) fitted through the half-depth point; 159 measured p 0.70 on H1 and 0.93–1.05 on the Orion lines (its 400-cell extrapolations 0.099–0.123 HCT116 / 0.255–0.273 HEK293T came from random-cell shards; 267 showed HEK293T's random halves read 2.5× its batch halves, so the batch-based 0.07–0.11 here is the conservative HEK293T number). Numbers, and what they imply for A/B/C:

- **Expected reliability of the challenge references** (400 cells, Flex): **0.48** if they behave like the H1 Flex screen that the panel was designed on (the same chemistry, the same depth); **0.07–0.20** if like a 10x/Orion line at 400 cells. The observed per-target cosine of a PERFECT prediction is √ρ: **0.69** (H1-like) vs **0.27–0.44** (Orion-like). PLOLO's deep truths sit at √ρ = **0.23–0.30**.
- **The conversion curve moves with the reference.** 262's knee (per-target `pds` ≈ 1 above an observed cosine of 0.02–0.03 on PLOLO) is a SIGNAL cosine of 0.02/0.30…0.03/0.23 ≈ **0.07–0.13**. On an H1-like A/B/C reference the same targets read 0.05–0.09 — saturated. The board's marginal band (observed 0.02–0.03 on A/B/C) is a signal cosine of ≈ **0.03–0.04**, which on PLOLO reads 0.007–0.012 — inside the references' own noise (30–34 % of Orion targets ≤ 0). **PLOLO's marginal band and the board's are different target sets**, and per-target claims read on PLOLO's band (262 H2, 263, this experiment's 3.4) do not name the targets the board is losing.
- **The PLOLO → board cosine factor has a reference-side explanation.** With the prediction fixed, observed cosine scales as √ρ of the truth: 0.69 / (0.23–0.30) = **2.3–3.0** for an H1-like reference, against the measured absolute ratio board 0.105 / PLOLO 0.025–0.038 = 2.8–4.2 and the class constants k_cos 1.36 (re-dose) … 3.76 (re-shape). Most of the factor is the truth's reliability, not the prediction's; the remainder (and the class dependence) is the prediction side.
- **Ceiling for the incumbent's per-target cosine on A/B/C**: at the board's measured aggregate cosine 0.105 and √ρ 0.69 the signal cosine is ≈ 0.15 (0.24–0.39 if the references are Orion-like), i.e. the 159 §2 reading (production at ~30 % of the context-independent ceiling) is unchanged by the reference side; what changes is that the deficit's location — signal cosines below ≈ 0.04 — cannot be measured on 150–200-cell references at all.

## 6. Where a per-target, per-donor weight would enter production (named, not made)

`scripts/make_submission_pm_cd4.py`: `build_v9_pre_whiten_effects(model, targets, alpha, beta, weights)` (line 1380) calls `MultiDonorModel.effect(target, alpha, beta, weights=active_weights)` per target with ONE dict `{donor: weight}` from `donor_weights()` (line 1342, the literal `{name: 1.0}` since v9; `VCC_DONOR_WEIGHTS` overrides it globally). A per-target table would replace `active_weights` by `weights_by_target[target]` inside that comprehension, renormalised over the carrying donors as `effect` already does — one call site, nothing else touched, `VCC_DONOR_WEIGHTS` unset reproducing every prior build. The K562 slot's PerturbMatch blend (`apply_perturbmatch_k562_replacement`, 0.5 on the shared columns) sits after it and would keep its own equal weights. **Nothing here earns that table**: the only positive cell moves the weights by ±0.02–0.08 around ½ and prices at a quarter of the gate.

## 7. Honest limits

- **One split unit (batch parity), two donors per fold, three cancer lines.** The halves resample batches of the SAME experiment; the estimators were tested with the reliability of a row measured in the row's own library. A donor's reliability against an independent study of the same line is not measured here (159 §3's disattenuation would be the instrument).
- **The K562 row is PerturbMatch-blended in production** (0.5 on 7,075 shared columns); its reliability and half differences come from the un-blended K562 halves, so the (a)/(b) weights and the (c) filters treat the blended row as if it carried K562's noise on those columns. This biases the K562 row's estimated noise UP on the shared columns; a correction could only make the K562 row's weight larger, and the fold where the K562 row is a donor and most reliable (HEK293T held) is where (a) reads +0.0022 — the direction of the bias does not rescue (a) on the other folds.
- **Whitened reliabilities are in-sample** (§1.2) and were declared as an estimator input before that was known; `relw_wht` is degenerate for that reason and reported as such. Cross-fitted whitening (operator on half A, rows from half B, and the mirror) was not run.
- **The (c) filters were frozen with one noise model each** (a one-degree-of-freedom own-variance estimate; a cell-count-scaled gene-pooled estimate). Their failure is on the `nmae` identity and on the marginal band, the same two places 224/228 failed; a better-regularised diagonal filter would have to escape both, and 228's mechanism (matched magnitude cannot collect a reshaping) is not about the estimate's quality.
- **a13's construction is re-implemented, not re-run**: the per-gene agreement damping at the builder's constants (γ 0.25, rank scale 2000, toward zero, degenerate rows untouched) on the composed effect; the graded fractions match the builder's documented values to 0.01, the direction and the `pds` sign match the board, but a13 also shipped an amplitude carrier (`--effect-scale`) and the per-perturbation half was not part of the flag that built it as far as the ledger records; the `:nm` / `:lit` pair brackets that.
- **Multiplicity**: 10 estimator cells × 3 folds, 6 strata, one band statistic, one frozen list; the one interval above zero (`ebrow`, P(≤ 0) 0.001) is reported with its control (P(arm ≤ control) 0.032) and its retention (31 %); no cell was chosen after a first look.
- **The 400-cell extrapolation** rests on 159's power law fitted on 15–100-cell shards (4–5 % error at 200 on H1) applied from a single half-depth point; the p = 0.70 / 1.0 spread is the error bar, and the H1-like / Orion-like dichotomy is a statement about chemistry that only the board's own replicate anchors (`r` = 0.927–0.984 on `pds`, AGENTS §4) can settle — those anchors are consistent with a high-reliability reference.

## 8. Files, memory, provenance

| file | what | peak RSS / time |
|---|---|---|
| `est271.py` (`fe1da9917a2a…`) | the frozen estimators, formulas and controls (sha recorded in every sweep meta; unchanged across folds, asserted by `analyse271.py`) | — |
| `halfrows271.py` (`6e402f345e54…`) → `work/half/rows-<line>.{npz,json,log}` (`8e3957ac…`, `41ecbe77…`, `c59dca6c…`), `work/half/<line>/bank-batch{A,B}/` (delta removed), `work/half/halfrows271.nohup` | the half rows in production units, KA-A / KA-B | 3.1 / 5.9 / 5.6 GiB; 0.1 / 0.4 / 0.3 min |
| `sweep271.py` (`c50104bd57ff…`) → `work/sweep271-<fold>.{json,log,nohup}` (`e80d9eaf…`, `c925f7fd…`, `60875686…`), `pertarget271-<fold>.npz` (`3513fe3b…`, `f7f3e377…`, `6cf19433…`), `identity271-<fold>.npz`, `halfwht271-<fold>.npz` | 20 rows per fold through 263's construction; 17 KA rows per fold; the reliability arrays, weights, per-target `pds` / cosine, the identity arrays | 5.5 / 5.7 / 5.7 GiB; 0.8 min each |
| `analyse271.py` (`613ecaeb4143…`) → `work/analysis271.json` (`f7b11e02…`), `work/tables271.md`, `work/analyse271.log` | tables D1, M1–M5, Q; paired bootstraps (20,000, seed 20260906) | < 1 GiB |
| `price271.py` (`d8a8601302bd…`) → `work/price271.{json,log,md}` (`f70e6e5a…`) | the retrodiction row, the six-member price, the CLEARED rule | < 1 GiB |
| `diagwht271.py` (`c094fc59a550…`) → `work/diagwht271.json` (`fb4de783…`) | in-sample vs out-of-sample whitened split-half reliability | 6.6 GiB; 2 min |
| `work/plolo/` | byte copies of 267's fold caches + `bctrl-challenge-*.npy`, `SHA256SUMS.txt`; `work/.marker-start` the write-audit marker | |

Read-only inputs: `experiments/267-…/work/half/<line>/{sums.npz, acc_batchA.npy, acc_batchB.npy}`, `…/bank-union/`, `…/truth-<line>.npz` (the reliability arrays; the reserve file untouched), `data/donors/{k562_gwps,orion_hct116,orion_hek293t}/`, `data/external/perturbmatch`, `experiments/263-…/{lib263.py, work/sweep263-*, pertarget263-*, picks263-*, cov263-*}`, `experiments/257-…/{lib257.py, price257.py, work/sweep257-*, work/price257.json}`, `experiments/261-…/price261.py`, 162's `plolo.py`, 246 / 218 libraries, 249's `work/amp249.json`, `data/board-entries/*.json`. No builder import for any computation (the builder was read for §6 only); no `vcc`; no `cell-eval2`; no submission; no preregistration (nothing clears).

---

## 9. Verdict

**DOES NOT CLEAR.** The donors' rows carry the same sampling noise 267 found in the references (split-half reliability 0.09 / 0.05 / 0.05 at full depth, a third of the Orion targets without a measurable direction), the reliability is a target property that the anti-aligned group shares across lines, and it knows a little about which donor is closer to the truth (52–61 % against a 66–70 % break-even). Correcting for it does not move the transfer: the reliability-weighted mean loses −0.0037 `pds` (its shuffled and reversed controls −0.011 / −0.019), the per-gene filters buy `pds` only on the anti-aligned group and pay +0.019 / +0.10 on the `nmae` identity, a13's agreement damping loses −0.0069 with a +12.6 % cosine (its board signature), and the mildest estimator — EB row shrinkage toward the consensus — is +0.0014 [+0.0005, +0.0024] but +0.0013 [−0.0001, +0.0027] over its shuffled control, 31 % retained leave-top-3-out, and +0.00050 pessimistic (0.26× the gate). The construction "estimate each donor's response by its own measured reliability before combining" is closed at this scope; what it leaves is a diagnosis (the production operator learns its own rows' noise: in-sample whitened halves anti-correlate at −0.13…−0.24 while an out-of-sample line's read +0.03…+0.07) and a scale (a 400-cell reference is 0.07–0.20 reliable on these lines' laws and 0.48 on Flex; PLOLO's marginal band is not the board's). No build, no submission.
