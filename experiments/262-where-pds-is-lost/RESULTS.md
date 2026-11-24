# 262 — Where `cd4svd5` loses `pds`: on a fold-consistent set of targets whose donor consensus does not point at the held line's response. The loss is NOT diffuse: it lives entirely on targets whose whitened-base cosine to the truth is below 0.05 (53 % of the deficit on the 17.5 % of targets where it is ≤ 0, 47 % on the 66 % where it is in (0, 0.05], 0 % above 0.05); `pds` is a threshold function of the per-target transfer cosine (0.08 at ≤ −0.02, 0.37 at (−0.01, 0], 0.66 at (0, 0.01], 0.86 at (0.01, 0.02], ≥ 0.97 above 0.02). Hub, donor support, truth size (by norm) and magnitude explain nothing (|ρ| ≤ 0.09). The same targets fail whichever cancer line is held out (bottom-decile overlap 4–5× chance; 8 targets in the bottom decile on all three folds; Spearman of per-target `pds` across folds +0.70…+0.76), they are fully supported (2/2 whitened donors, CD4 3/3, mostly union-absent), and their responses do not agree across real lines (held-vs-donor truth cosine ρ +0.63 with `pds`; CD4's own DE count ρ +0.38). The `h1union` → `cd4svd5` step's per-target sign is line-specific (Spearman of Δ across folds −0.10…+0.03); it paid on the marginal band (0.835 → 0.858 on 564 targets) and cost on the anti-aligned group (0.403 → 0.389). Two per-target routes the diagnosis motivates were piloted and read negative: keying the CD4 row weight to observable donor disagreement (`pds` +0.00045, 3/3 but 0.3 σ, at −3.5 % aggregate cosine) and choosing the donor per target with CD4 as arbiter (a coin flip: the arbiter matches the best donor 43–62 % of the time). **Measurement only. No price, no candidate, no preregistration. Verdict line at the bottom (§10).**

**Date 2026-09-05 (03:38–04:08 UTC).** Deterministic. **No build, no `vcc`, no `cell-eval2`, no model name, no description, nothing under `submissions/`, nothing written outside this directory** (`find -newer work/.marker-start -type f` over `experiments/ data/ scripts/ src/ docs/ submissions/` returns nothing outside `experiments/262-where-pds-is-lost/` other than sibling agents' own `claude-run.*` logs; no source file changed; builder not imported). Fold caches are byte copies of 257's (= 261's = 249's = 246's; `work/plolo/SHA256SUMS.txt`: `c13c1cec…`, `9f3e0452…`, `f9a3cc30…`, `d86202a9…`; `bctrl-challenge-{A,B,C}.npy` `f1f04360…`, `4db8bd10…`, `cdebb5ea…`; `cmp` byte-identical). One PLOLO fold per process for the diagnostic (peak RSS **4.81 / 4.92 / 4.94 / 4.86 GiB**, 8–10 s each); the pilot staged the three cancer folds serially in one process (**8.85 GB** peak, 33 s — above the mandate's ~8 GB target, below the 10 GB rule; noted). Memory gate at launch: 75 % free, free+inactive 27 GiB, swap free 1.6 GB, no ≥ 1 GB process live. HEAD `f7bf79a`. `uv run python` throughout.

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** the `cd4svd5` construction (261's `cd4svd5_parts`: CD4 deflated k = 5, RE-SOLVED at f = 0.34, + the H1-pair union at g = 0.28) reproduces 257's `k5_s1` == 261's R1 | 1e−9 on (`cos_agg`, macro, `pds`); per-target arrays bit-identical | **0.0e+00** on the endpoints and **0.0e+00** on the per-target `pds` and cosine arrays vs `pertarget257-<fold>.npz` `k5_s1|*` AND vs `pertarget261-<fold>.npz` `cd4svd5|pds`, on **4/4 folds** (`pds` 0.815384 / 0.801659 / 0.780978 / 0.640294; `cos_agg` 0.037715 / 0.030626 / 0.024857 / 0.013640) | **PASS** |
| **KA-2** `h1union` (the un-deflated stack, `lib246.stack_parts`) reproduces 257's `STACK` | as above | **0.0e+00** endpoints and per-target arrays, 4/4 folds (`pds` 0.805397 / 0.792213 / 0.763623 / 0.636155) | **PASS** |
| **KA-3** this file's own per-target rank statistic, recomputed from the (prediction × truth) cosine matrices it saves, equals the instrument's `_pds_per_target_a23` | < 1e−12 | **0.0** (bit-identical) on 5 arms × 4 folds = 20/20 rows | **PASS** |
| **KA-4** pilot: β = 0 equals `cd4svd5` | 1e−9 | **0.0e+00** on 3/3 cancer folds | **PASS** |
| the `cd4svd5` increment over `h1union`, per fold | 257 R2: +0.00999 / +0.00945 / +0.01736 (mean +0.012262757) | **+0.00999 / +0.00945 / +0.01736** (sums of the per-target Δ; §6) | **PASS** |

---

## 0. The answer in twelve lines

1. **The loss has an address and it is the same address on every fold.** Per-target `pds` is bimodal: 22–30 % of targets sit at exactly 1.0, 15–17 % below 0.5; the bottom decile carries **37–42 %** of the deficit and the bottom quartile **71–79 %** (§3, T1). The bottom deciles of the three cancer folds overlap **9–12** targets against **2.2–2.5 expected** (Jaccard 0.21–0.33); per-target `pds` correlates **+0.70…+0.76** across folds; **8 targets are in the bottom decile on all three** (`CAT`, `LPCAT3`, `SAP30L`, `SNRK`, `TRMT1`, `ZBTB41`, `ZNF32`, `ZNF783`) and 17 on two or more (§3, T2).
2. **One covariate explains it and it is the per-target transfer cosine.** The whitened donor base's cosine to the truth (effect level, before any residual) has Spearman **+0.875** with `pds` pooled (+0.91 / +0.89 / +0.86) and separates the bottom decile with AUC **0.05**; the residual's cosine to the truth +0.61. Everything the mandate named as a candidate cause is near zero: prediction hubness **−0.03**, median cosine to the other truths −0.07, prediction nearest-neighbour cosine **+0.26 (the wrong sign for confusability)**, whitened-source count +0.07, CD4 completeness +0.07, union carriers +0.21, CD4 energy removed by the deflation +0.05, truth ‖Δ‖ **−0.03**, truth nearest-neighbour cosine +0.09, delivered ‖pred‖/‖truth‖ +0.16, the render +0.12 (§4, T4).
3. **`pds` is a threshold function of that cosine** (§4.2): mean `pds` **0.08** at own cosine ≤ −0.02, **0.15** at (−0.02, −0.01], **0.37** at (−0.01, 0], **0.66** at (0, 0.01], **0.86** at (0.01, 0.02], **0.97** at (0.02, 0.03], ≥ **0.99** above 0.03; slope **0.095 `pds` per +0.01 cosine** inside (0, 0.05]. So 57 % of the deficit sits on the 16 % of targets with own cosine ≤ 0 (below chance: a neutral row would score 0.5) and 40 % on the 33 % with cosine in (0, 0.02]; nothing above 0.03. A source arm's `pds` is the count of targets it pushes across ~0.02, which is why 257–261 all found it on three or four targets per fold.
4. **The deficit decomposed** (§5, T5c; targets assigned to the first matching group): **(e) base points away** (base·truth ≤ 0) **49–58 %** of the deficit on 46–52 targets; **(r) the residual breaks a positive base** 12–18 %; **(a) hub / confusable** 3–9 %; **(c) small truth by norm** 3–9 %; **(d) render** 9–12 %; **(b) no base support** 0–3 %; **(o) specific confusion** 7–9 %. By bins of the base cosine: (−1, 0] 17.5 % of targets / `pds` 0.389 / **53.1 %** of the deficit; (0, 0.05] 66 % / 0.858 / **46.8 %**; > 0.05 16 % / ≥ 0.998 / **0.1 %** (T5d).
5. **There is no hub left.** The 10 most-beating truths take **4.3–4.5 %** of all beats against a uniform 3.4–3.7 %; every truth beats something; prediction hubness is uncorrelated with `pds` on `h1union` already (−0.09 / −0.02 / −0.02) and after the deflation (−0.12 / +0.01 / −0.01) (T5a, T5b). The residual-hub story that the deflation fixed was a hub of the *un-deflated CD4 rows over the marginal band* (`h1union` 0.835 < `BASE` 0.853 there), not a truth-side hub.
6. **What the anti-aligned targets are**: fully supported (all 17 recurrent ones: 2/2 whitened donors, CD4 3/3 on 15, union-absent on 15), with **held-vs-donor cross-line truth agreement ≈ 0 or negative** (pooled ρ **+0.63** with `pds`, AUC 0.12) and **donor-donor agreement ≈ 0** (ρ +0.38, AUC 0.37 — the build-time-observable version), and weak responses in an independent line: **CD4's own DE count** (adj p < 0.05, three conditions, median per target) ρ **+0.38** pooled; targets with < 10 CD4 DE genes are 68 % of the panel and carry **78 %** of the deficit; targets with ≥ 100 (9 %) sit at `pds` 0.95–0.98 and carry 1.6 % (§4.3). "Small true effect" IS part of the story once effect size is measured by cross-line reproducibility rather than by the held truth's norm, which is noise-dominated (bottom-decile ‖truth‖ 7.0 vs 6.7).
7. **The components** (§5.3): `BASE` (no residual) 0.781 / 0.782 / 0.774; `CD4ONLY` 0.807 / 0.800 / 0.782; `UNIONONLY` 0.792 / 0.787 / 0.771; `h1union` 0.805 / 0.792 / 0.764; `cd4svd5` 0.815 / 0.802 / 0.781. ρ(`BASE` `pds`, `cd4svd5` `pds`) +0.90…+0.93; **52–61 % of `cd4svd5`'s deficit sits on targets already below 0.5 in `BASE`**; the residual stack lifts 12–17 of ~50 such targets above 0.5 and drops 0–2 of ~150 good ones. The union's whole `pds` is the anti-aligned targets it carries: on group (e) ∧ carried (n 22) `BASE` 0.279 → `UNIONONLY` 0.441 (+0.16 per target); on the marginal band ∧ carried (n 144) 0.877 → 0.878 (nothing).
8. **The `h1union` → `cd4svd5` step per target** (§6): gainers / losers / unchanged **102 / 71 / 96**, **112 / 88 / 90**, **134 / 83 / 77**; sums of gains +6.5 / +8.6 / +10.1 against losses −3.8 / −5.8 / −5.0 `pds`-units. **The per-target Δ is fold-inconsistent**: Spearman across cancer folds **+0.03 / −0.03 / −0.10**; the biggest movers flip sign (`PLCG1` +0.39 / −0.62 / +0.66, `PGLS` +0.50 / −0.53 / +0.09). The one covariate that tracks the loss on every fold is unobservable: the CD4 row's cosine to the held truth BEFORE deflation (ρ(Δ, CD4·truth) **−0.36 / −0.36 / −0.38 / −0.35**) — the projector removed signal that was real for that line. Losers and gainers do not differ on any structural covariate (energy removed AUC 0.47, CD4 norm 0.44, support 0.49, union 0.43); both are drawn from the marginal band (base·truth medians 0.008 / 0.012 vs 0.048 for the unchanged).
9. **Pilot (a)** (§7): keying the deflated CD4 row weight to observable donor-donor disagreement, w = exp(β(0.5 − rank)), re-solved sizing: β = 1 **`pds` +0.00045 (3/3 folds, 0.3 LOO-σ) at −3.46 % aggregate cosine**; shuffled weights −0.00014 at +1.40 %; reverse key −0.00140 (0/3) at +4.19 %; β = 2 −0.00072 at −6.8 %; β = 4 −0.0081 at −20 %. The key's sign is right and its value is negligible: the disagreeing rows are rows where CD4 is unreliable too. **Read negative; not a candidate.**
10. **Pilot (b)** (§7): choosing the donor per target instead of averaging has ORACLE headroom (best single donor's agreement with the held truth > 0.05 on 20–50 % of the marginal band vs 6–26 % for the mean; on group (e) the best donor averages +0.00…+0.015 vs the mean's −0.008…−0.031) but **CD4 as arbiter picks the better donor 43–62 % of the time and lands at the random pick** (group (e): arbiter −0.006 / −0.037 / −0.017 vs random −0.007 / −0.035 / −0.024). No arbiter on disk finds the headroom. **Read negative.**
11. **Observable selectors cannot find the below-chance group**: the best build-time-observable bin (bottom 5 % by min(donor-donor agreement, CD4·base cosine)) has mean `pds` **0.685**; donor-donor ≤ 0 (25 % of targets) 0.707; the oracle group (own cosine ≤ 0) 0.295. Neutralising any observable bin would lose `pds` (and a neutral row scores 0 on `fid` and `reach` by definition, brief §4–5), so the "null the anti-aligned rows" idea is dead on both counts (§7.3).
12. **Diagnosis, ranked**: (1) information — the donor consensus has no transferable direction for ~17 % of the panel (weak, line-specific responses; cross-line agreement ≈ 0), and those targets score BELOW chance, 53 % of the deficit; (2) the marginal band — 66 % of the panel at cosine (0, 0.05], where +0.01 cosine buys +0.095 `pds` per target, 47 % of the deficit, and where the deflation's gain lived; (3) the render, ≤ 0.01 `pds` (156 owns it); (4) hubs, support, truth size by norm, magnitude — nothing. Three hypotheses and their screens in §8; two of them already carry a measured negative from the pilots, the third is the reading rule the next arm must satisfy. `kolf2` (apart, §3.4) shows the same structure with a larger deficit (0.360) and a residual hub on the production operator (ρ(median competitor cosine, `pds`) −0.18, AUC 0.61).

---

## 1. Instrument and construction (`diag262.py`; nothing mathematical is new)

162's PLOLO, three cancer folds (`k562_gwps` 269 targets, donors `orion_hct116` + `orion_hek293t`, operator `held_k562_gwps` `13d988fc…`; `orion_hct116` 290, `k562_gwps` + `orion_hek293t`; `orion_hek293t` 294, `k562_gwps` + `orion_hct116`) and `kolf2` (273, all three donors, the production operator `e2f65ecf…`) reported apart. 246's construction imported through 257 (`lib257` → `lib246`: 218's residuals and `union_rows`, `deflate`, `compose`, `score`); the incumbent is 261's `cd4svd5_parts` — `lib246.deflate(cd4, "svd", 5)`, re-solved by `plolo.solve_fraction` at f = 0.34 (fractions 0.395879 / 0.394269 / 0.397869 / 0.336448 vs the stack's 0.373761 / 0.370782 / 0.376056 / 0.316408; energy removed 0.2022 / 0.2001 / 0.2000 / 0.2046), merged with the H1-pair union at its shipped carrying-row fraction for g = 0.28 (0.223670 / 0.222539 / 0.225247 / 0.189583), one add, one preserve. `h1union` = `lib246.stack_parts`. Three more objects for attribution: `BASE` (the whitened base, no residual), `CD4ONLY` (the deflated CD4 at its resolved fraction, no union), `UNIONONLY` (the union at its shipped fraction, no CD4).

**Per-target `pds`** is the instrument's (`resid218.score_effects`, tag a23: the effect rendered by production's `profile_from_effect` with the minnorm mass fix on each of the three challenge frames A/B/C and read back as the scorer does, the rank statistic per frame averaged over frames) — and it is the scorer's definition (brief §1): for prediction p, the fraction of the other truths q whose cosine to p̂ is below p's own, ties at half. This file recomputes it from the (prediction × truth) cosine matrix it saves (KA-3, bit-identical) and keeps the matrix, so the loss can be read pair by pair: the **beat matrix** B[p, q] = number of frames (0–3) in which truth q ranks above p's own truth for prediction p.

**Covariates per target** (all from data already on disk): which of the fold's whitened sources carry the target (bank `perts.txt`), PerturbMatch blending, CD4 condition count (0–3) and complete-row norm, the per-row energy the k = 5 deflation removed, union carriers (`gse295214_h1`, `h1_2025`) and row norm, the `k562_essential` flag (0 on every fold's targets — measured: `n_essential_on_panel` 0 / 0 / 0 / 0, §8.1g), the whitened base row's norm and cosine to the truth (effect level, scored columns), the merged residual's norm, cosine to the truth and to the base and its share of the row, CD4's and the union's own cosines to the truth and to the base; from the scorer geometry: own cosine, competitor max / median / mean, prediction-prediction and truth-truth nearest-neighbour and mean cosines, delivered ‖pred‖ and ‖truth‖; the effect-level (un-rendered) `pds` per target, so the render's per-target effect is the difference. Two covariates defined outside the pipeline, from the fold caches' raw truths: the **cross-line truth agreement** — for each pair of the four lines, the cosine between the two lines' centred truths for the same target on the intersection of their scored columns (line pairs: mean 0.020–0.034, 26–37 % of targets ≤ 0) — read as *held-vs-donor* (mean over the fold's donors: the transferability of the target's response into the held line) and *donor-donor* (between the fold's donors: observable at build time). And **CD4's own DE count** per target (`adj_p_value.npy` < 0.05 over the non-panel genes, median over the three conditions; panel median 3 genes).

## 2. Known-answer rows

Table at the top. Beyond it: the fold caches validated `legacy` field by field by `fold_cache_status` on every fold; the union's carrying rows 74 / 75 / 76 / 71; CD4 complete rows 251 / 267 / 271 / 253 (counts hist `{0:5,1:4,2:9,3:251}`, `{0:8,1:4,2:11,3:267}`, `{0:8,1:4,2:11,3:271}`, `{0:7,1:3,2:10,3:253}`); whitened-source coverage `{0:5,1:5,2:259}`, `{1:31,2:259}`, `{0:1,1:34,2:259}`, `{1:5,2:28,3:240}` (the fold audits' vote histograms, reproduced from the banks); PerturbMatch rows 0 / 239 / 242 / 228.

## 3. Per-target `pds` of `cd4svd5`

### 3.1 Distribution (T1 in the appendix)

| fold (held) | n | mean | deficit | q10 | q25 | median | at 1.0 | < 0.5 | < 0.9 | bottom-decile share of deficit | bottom-quartile share |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 269 | 0.8154 | 0.1846 | 0.389 | 0.711 | 0.958 | 82 | 43 | 111 | 42.3 % | 78.7 % |
| `orion_hct116` | 290 | 0.8017 | 0.1983 | 0.382 | 0.698 | 0.939 | 77 | 50 | 126 | 40.4 % | 75.8 % |
| `orion_hek293t` | 294 | 0.7810 | 0.2190 | 0.343 | 0.633 | 0.905 | 65 | 51 | 146 | 36.8 % | 70.8 % |
| `kolf2` (apart) | 273 | 0.6403 | 0.3597 | 0.166 | 0.402 | 0.696 | 34 | 90 | 195 | 25.6 % | 55.1 % |

The ceiling used throughout is `pds` = 1 on every target (every prediction nearest its own truth), so the fold deficit is 1 − mean; the board's replicate anchor (0.927–0.984, brief §8) is not reachable on PLOLO's deep pseudobulk truths and is not used.

### 3.2 Bottom and top deciles, and their agreement across folds (T2, T3 in the appendix)

Bottom deciles (27 / 29 / 30 targets): K562 `SAP30L DEAF1 PIN1 TMF1 SNRK SNX4 LPCAT3 VIM ZNF641 SMARCC2 TRMT1 ZNF783 COMMD8 CAT SLC39A1 ZNF32 ZBTB41 ZKSCAN3 TRIM35 FZD2 RNF141 YIPF6 SH2D3A MAP2K2 RAB5A NIPSNAP3A DVL3`; HCT116 `RBM4B SAP30L IFNAR2 ZNF543 ZNF496 ZNF32 PGLS PYCR1 ZNF783 ZBTB41 RAB11FIP5 PLCG1 PHF19 ZNF786 NICN1 DTNBP1 SNRK CAT BCL7C ZNF641 ZNF219 SNN POLD4 TRMT1 PDK1 LPCAT3 TARBP2 RFK IREB2`; HEK293T `NFE2L1 SAP30L RBM4B ZBTB41 IFNAR2 IFI27L2 MYD88 PDK1 ZNF783 HEATR6 RAB11FIP5 SNN TBC1D19 ZNF32 DEAF1 PSMB9 NT5DC1 LPCAT3 PYCR1 TDRD7 ZKSCAN3 CAT SFXN4 VIM TIMM50 FHL3 SNRK LMAN1 RNF141 TRMT1`. Top deciles are the 1.0 block (82 / 77 / 65 targets); 14 targets are in the top decile of all three folds (`RBM18 RBM26 RBM5 RNF10 SMARCA5 SYVN1 TBCK THAP7 UPF3B WBP4 ZDHHC4 ZFP62 ZNF462 ZNF507`).

| pair (259 common targets) | overlap of bottom deciles | expected | Jaccard | Spearman of per-target `pds` |
|---|---:|---:|---:|---:|
| K562 vs HCT116 | **9** | 2.5 | 0.214 | +0.756 |
| K562 vs HEK293T | **12** | 2.5 | 0.308 | +0.705 |
| HCT116 vs HEK293T | **12** | 2.2 | 0.333 | +0.741 |

**The bottom decile is 4–5× more shared than chance and the whole per-target profile correlates +0.70…+0.76 across held lines: the loss is a property of the target, not of the fold.** The 17 targets in the bottom decile on ≥ 2 folds, with the covariates that matter (means over the cancer folds):

| target | `pds` K562 / HCT116 / HEK293T / KOLF2 | base·truth cos | held-vs-donor agreement | donor-donor agreement | CD4 | union | whitened donors |
|---|---|---:|---:|---:|---:|---:|---|
| `CAT` | 0.24 / 0.25 / 0.29 / 0.68 | −0.007 | −0.005 | −0.005 | 3/3 | 0 | 2/2/2 |
| `DEAF1` | 0.06 / 0.43 / 0.23 / 0.38 | −0.020 | −0.085 | −0.085 | 3/3 | 0 | 2/2/2 |
| `IFNAR2` | 0.78 / 0.06 / 0.09 / 0.77 | −0.005 | −0.020 | −0.020 | 3/3 | 0 | 2/2/2 |
| `LPCAT3` | 0.13 / 0.34 / 0.25 / 0.86 | −0.000 | −0.004 | −0.004 | 3/3 | 0 | 2/2/2 |
| `PDK1` | 0.39 / 0.33 / 0.14 / 0.83 | −0.004 | −0.019 | −0.019 | 2/3 | 0 | 2/2/2 |
| `PYCR1` | 0.48 / 0.09 / 0.26 / 0.36 | −0.004 | −0.034 | −0.034 | 3/3 | 0 | 2/2/2 |
| `RBM4B` | 0.93 / 0.03 / 0.04 / 0.67 | −0.006 | −0.087 | −0.087 | 3/3 | 0 | 2/2/2 |
| `RNF141` | 0.31 / 0.67 / 0.33 / 0.67 | +0.002 | +0.022 | +0.022 | 3/3 | 0 | 2/2/2 |
| `SAP30L` | 0.02 / 0.03 / 0.02 / 0.46 | −0.013 | +0.005 | +0.005 | 1/3 | 0 | 2/2/2 |
| `SNRK` | 0.10 / 0.22 / 0.33 / — | −0.010 | −0.003 | −0.003 | 3/3 | 1 | 2/2/2 |
| `TRMT1` | 0.18 / 0.32 / 0.34 / 0.04 | −0.008 | +0.002 | +0.002 | 3/3 | 0 | 2/2/2 |
| `VIM` | 0.15 / 0.82 / 0.30 / 0.22 | −0.007 | −0.021 | −0.021 | 3/3 | 0 | 2/2/2 |
| `ZBTB41` | 0.27 / 0.10 / 0.06 / 0.14 | −0.014 | −0.028 | −0.028 | 3/3 | 0 | 2/2/2 |
| `ZKSCAN3` | 0.29 / 0.53 / 0.27 / 0.92 | −0.003 | −0.019 | −0.019 | 3/3 | 0 | 2/2/2 |
| `ZNF32` | 0.26 / 0.09 / 0.22 / 0.18 | −0.017 | −0.050 | −0.050 | 3/3 | 1 | 2/2/2 |
| `ZNF641` | 0.17 / 0.27 / 0.67 / 0.43 | −0.002 | −0.004 | −0.004 | 3/3 | 0 | 2/2/2 |
| `ZNF783` | 0.20 / 0.10 / 0.15 / 0.39 | −0.009 | −0.025 | −0.025 | 3/3 | 0 | 2/2/2 |

(Held-vs-donor and donor-donor agreement coincide in the mean over the three folds because each line is the held line on one fold and a donor on the other two.) Every one of them is fully supported and none of them has a transferable direction anywhere: the three cancer lines do not agree with each other on these targets' responses (cosines −0.09…+0.02 against a panel mean of +0.02–0.03), and CD4 finds 1–2 DE genes for them (median; the 1.0 block: 12).

### 3.3 The full per-target tables

Appendix T3b/T3t: the bottom decile of every fold with 17 covariates per target, and the top decile's names with its medians against the rest. Bottom-decile medians vs the rest, K562 / HCT116 / HEK293T: own cosine −0.011 / −0.013 / −0.011 vs +0.020 / +0.022 / +0.021; base·truth −0.013 / −0.010 / −0.008 vs +0.019 / +0.020 / +0.018; held-vs-donor agreement −0.016 / −0.026 / −0.021 vs +0.020 / +0.028 / +0.029; ‖truth‖ 7.0 / 7.9 / 6.0 vs 6.7 / 7.0 / 6.2 (no separation); truth NN cosine 0.07 / 0.17 / 0.30 vs 0.08 / 0.13 / 0.26; prediction NN 0.04 everywhere; competitor median 0.000 everywhere; whitened donors 2, CD4 3, union 0 on both sides; ‖CD4‖ 4.1–4.4 vs 4.5; CD4 energy removed 0.09–0.10 on both sides.

### 3.4 `kolf2`, apart

Same structure, larger deficit (0.360): base·truth ρ +0.79, held-vs-donor agreement +0.51, residual·truth +0.73, union·truth +0.80; bottom decile `MBD4 ZNF768 TRMT1 MSANTD4 RABEP2 YIPF6 ZNF714 ABCD1 SHPRH TWF2 DVL3 POLD4 GTF2IRD1 SMYD2 ZNF251 PTOV1 ZNF280C MARCKSL1 TFCP2 CAPRIN2 MTHFSD TMEM209 MYRF ZBTB41 PDE6D SND1 PHF19 SPATA20` — only `TRMT1` and `ZBTB41` shared with the cancer folds' recurrent set (the iPSC line's failures are its own). Unlike the cancer folds, a competitor-side hub is measurable on the production operator: ρ(median cosine to the other truths, `pds`) **−0.18** and ρ(mean, `pds`) −0.34 (AUC 0.61 / 0.71 for the bottom decile) — a shared component that the held operators do not leave, on the one fold that uses the shipped operator and a non-cancer truth; noted, not pursued (the truth side of an iPSC line is the wrong instrument for it).

## 4. The covariates (T4 in the appendix)

### 4.1 Ranking by effect size (three cancer folds; Spearman ρ with `pds` per fold and pooled by within-fold rank; AUC = P(bottom-decile target has the higher value); medians bottom decile vs rest)

| covariate | ρ K562 / HCT116 / HEK293T | **ρ pooled** | **AUC pooled** | bottom / rest median |
|---|---|---:|---:|---|
| base·truth cosine (whitened donor base vs the truth) | +0.913 / +0.885 / +0.855 | **+0.875** | **0.050** | −0.009 / +0.020 |
| held-vs-donor cross-line truth agreement (pipeline-independent) | +0.769 / +0.565 / +0.567 | **+0.630** | **0.121** | −0.014 / +0.027 |
| residual·truth cosine (fr·CD4_svd5 + g·union) | +0.619 / +0.609 / +0.618 | +0.614 | 0.220 | −0.002 / +0.019 |
| deflated CD4 row·truth | +0.577 / +0.574 / +0.553 | +0.566 | 0.242 | −0.002 / +0.016 |
| residual·base cosine | +0.586 / +0.579 / +0.514 | +0.558 | 0.283 | +0.004 / +0.016 |
| union row·truth (carried rows) | +0.569 / +0.625 / +0.448 | +0.519 | 0.209 | +0.003 / +0.030 |
| residual row norm | +0.435 / +0.344 / +0.395 | +0.387 | 0.377 | 1.5 / 1.7 |
| **donor-donor cross-line agreement (observable at build time)** | +0.219 / +0.402 / +0.527 | **+0.381** | **0.371** | +0.005 / +0.024 |
| **CD4's own DE count (independent-line effect strength)** | +0.425 / +0.347 / +0.378 | **+0.381** | **0.35–0.40** | 1–2 / 3 (1.0 block 11–12) |
| residual share of the row | +0.340 / +0.315 / +0.345 | +0.329 | 0.391 | 0.30 / 0.35 |
| CD4 row norm | +0.371 / +0.294 / +0.328 | +0.320 | 0.411 | 4.1 / 4.5 |
| prediction nearest-neighbour cosine (confusability) | +0.216 / +0.295 / +0.264 | **+0.255** | 0.397 | 0.036 / 0.039 |
| union carriers / union row norm | +0.30 / +0.21 / +0.22 | +0.213 / +0.226 | 0.384 | 0 / 0 |
| delivered ‖pred‖ / ‖truth‖ | +0.098 / +0.166 / +0.196 | +0.158 | 0.451 | 0.54 / 0.57 |
| max cosine of the prediction to another truth | +0.190 / +0.065 / +0.080 | +0.116 | 0.460 | 0.032 / 0.034 |
| render effect (delivered − effect-level `pds`) | +0.010 / +0.155 / +0.183 | +0.115 | 0.432 | −0.016 / 0.000 |
| truth nearest-neighbour cosine | +0.304 / −0.039 / +0.009 | +0.086 | 0.498 | 0.119 / 0.124 |
| PerturbMatch row | — / +0.149 / +0.097 | +0.079 | 0.428 | |
| CD4 complete / condition count | +0.09 / +0.04 / +0.08 | +0.074 / +0.066 | 0.397 / 0.408 | 1 / 1 |
| whitened sources carrying the target | +0.120 / +0.104 / +0.086 | +0.068 | 0.439 | 2 / 2 |
| CD4 row energy removed by the deflation | +0.129 / +0.038 / +0.014 | +0.052 | 0.478 | 0.085 / 0.100 |
| `k562_essential` flag | constant 0 | — | — | 0 / 0 |
| truth's mean cosine to the other truths | −0.202 / +0.067 / +0.092 | −0.004 | 0.510 | |
| whitened base row norm | +0.004 / −0.017 / −0.052 | −0.019 | 0.503 | 4.7 / 4.6 |
| truth ‖Δ‖ | +0.070 / −0.054 / −0.087 | **−0.029** | 0.518 | 7.0 / 6.7 |
| prediction's mean cosine to the other predictions (hubness) | −0.115 / +0.006 / −0.008 | **−0.034** | 0.488 | |
| median cosine of the prediction to the other truths | −0.075 / −0.015 / −0.124 | **−0.072** | 0.540 | 0.000 / 0.000 |

Geometry rows that are the loss itself rather than a cause (for the record): own cosine ρ +0.958 (AUC 0.003); own minus best competitor +0.939; `BASE` `pds` +0.900; `CD4ONLY` `pds` +0.973; `UNIONONLY` +0.924; `h1union` +0.945; effect-level `pds` +0.945.

### 4.2 The conversion curve (cancer folds pooled, n = 853): `pds` is a threshold in the per-target cosine

| own cosine (delivered prediction vs its truth) | n | share of targets | mean `pds` | share of the deficit |
|---|---:|---:|---:|---:|
| ≤ −0.02 | 12 | 1.4 % | 0.079 | 6.4 % |
| (−0.02, −0.01] | 33 | 3.9 % | 0.147 | 16.4 % |
| (−0.01, 0] | 94 | 11.0 % | 0.374 | **34.3 %** |
| (0, 0.01] | 143 | 16.8 % | 0.659 | **28.4 %** |
| (0.01, 0.02] | 142 | 16.6 % | 0.860 | 11.6 % |
| (0.02, 0.03] | 132 | 15.5 % | 0.970 | 2.3 % |
| (0.03, 0.05] | 110 | 12.9 % | 0.993 | 0.5 % |
| (0.05, 0.08] | 89 | 10.4 % | 0.998 | 0.1 % |
| > 0.08 | 98 | 11.5 % | 1.000 | 0.0 % |

Least-squares slope inside (0, 0.05]: **9.5 `pds` per unit cosine = +0.095 per +0.01**. The whole deficit sits below cosine 0.03; 57 % of it below 0 (16 % of targets, each scoring below the 0.5 a neutral row would get); 40 % in (0, 0.02]. Above 0.03 the metric is saturated and no cosine gain there buys anything. This is why every source arm's `pds` since the board points was "three or four targets per fold": a source pays only where it carries a target sitting near the threshold, and moves it across.

### 4.3 Effect strength measured in an independent line (CD4's DE count; `work/cd4sig262.log`)

| CD4 DE genes per target (median of 3 conditions) | n | share | mean `pds` | share of the deficit | in group (e) |
|---|---:|---:|---:|---:|---:|
| < 10 | 566 | 68.0 % | 0.769 | **78.2 %** | 18.6 % |
| 10–29 | 111 | 13.3 % | 0.825 | 11.6 % | 18.9 % |
| 30–99 | 82 | 9.9 % | 0.825 | 8.6 % | 20.7 % |
| 100–299 | 40 | 4.8 % | 0.978 | 0.5 % | 0.0 % |
| ≥ 300 | 33 | 4.0 % | 0.946 | 1.1 % | 3.0 % |

ρ(CD4 DE count, `pds`) +0.43 / +0.35 / +0.38 (pooled +0.38); ρ(CD4 DE count, base·truth) +0.41 / +0.30 / +0.34. The panel's "negligible" and "subtle" bins (AGENTS.md §2.1: built by design from a low-depth H1 pilot) are where the loss is; a target with a strong, reproducible response transfers (the ≥ 100-gene tier is at 0.95–0.98). The held truth's norm cannot see this (ρ −0.03) because a deep pseudobulk of a negligible response has the same norm as a real one — its norm is sampling noise.

## 5. Which covariate explains the loss, and the deficit decomposition

### 5.1 The four mechanisms the mandate named, measured

- **(a) the prediction confusable with other predictions / the residual hub** — **no, and nothing remains of it.** Prediction hubness ρ −0.03 (AUC 0.49); the prediction's median cosine to the other truths −0.07 (0.54); its nearest-neighbour cosine to another prediction is *positively* related to `pds` (+0.26: well-predicted targets sit in populated, transferable programmes). Truth-side beat concentration (T5a): the 10 most-beating truths take **4.5 / 4.3 / 4.4 %** of the beats against a uniform 3.7 / 3.4 / 3.4 %; the top 30 take 12–13 %; no truth beats nothing. Before the deflation (`h1union`) the same numbers were 4.7 / 4.6 / 4.6 % and hubness ρ −0.09 / −0.02 / −0.02: the residual-hub story never was a truth-side hub; what the deflation removed was a shared CD4 component that depressed the marginal band (§5.3), and it is gone.
- **(b) low coverage / weak donor support** — **no.** Whitened-source count ρ +0.07, CD4 completeness +0.07, union carriers +0.21 (the union helps, on 22 anti-aligned carried targets, §5.3), PerturbMatch +0.08. Uncovered rows (the mean-response filler) are 5 / 0 / 1 targets and 3.2 / 0.0 / 0.1 % of the deficit. Every recurrent bottom-decile target has 2/2 donors and CD4 3/3.
- **(c) small true effects** — **no by norm (ρ −0.03), YES by reproducibility.** The held truth's norm is noise-dominated; the cross-line agreement of the truths (ρ +0.63) and CD4's DE count (+0.38) are the effect-size measures that work, and they say the anti-aligned group is the panel's weak-response tier: 68 % of targets with < 10 CD4 DE genes carry 78 % of the deficit.
- **(d) magnitude error at correct direction** — **no, by construction and by measurement.** `pds` is invariant to the row's scale; the only magnitude path is the render's non-linearity, which costs −0.0009 / −0.0095 / −0.0046 `pds` (mean, delivered minus effect level; `BASE` −0.007 / −0.013 / −0.004) and is large per target only where the rank is unstable (bottom-decile median −0.016 vs 0.000). ‖pred‖/‖truth‖ ρ +0.16 (AUC 0.45).
- **(e) — the one that explains it**: the donor consensus's direction does not transfer for the target: base·truth ρ +0.875, AUC 0.05; and it does not transfer because the target's response is not shared across lines (held-vs-donor agreement +0.63; donor-donor +0.38; CD4 DE count +0.38).

### 5.2 The partition (T5c; each target assigned to the first group it matches; shares of the fold deficit, ceiling `pds` = 1)

| fold | deficit | (b) no base | **(e) base points away** | (r) residual breaks a positive base | (a) hub quintile | (c) small-truth quintile (by norm) | (d) render > 0.05 | (o) specific confusion | n per group |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `k562_gwps` | 0.1846 | 3.2 % | **57.5 %** | 11.7 % | 3.0 % | 8.7 % | 8.7 % | 7.2 % | 5 / 47 / 14 / 21 / 24 / 19 / 139 (82 at zero deficit) |
| `orion_hct116` | 0.1983 | 0.0 % | **51.5 %** | 16.1 % | 7.8 % | 7.2 % | 10.6 % | 6.8 % | 0 / 46 / 21 / 31 / 27 / 28 / 137 (77) |
| `orion_hek293t` | 0.2190 | 0.1 % | **49.0 %** | 18.1 % | 8.6 % | 2.8 % | 12.3 % | 9.2 % | 1 / 52 / 25 / 28 / 20 / 34 / 134 (65) |
| `kolf2` (apart) | 0.3597 | 0.0 % | **64.0 %** | 10.6 % | 5.6 % | 4.1 % | 11.6 % | 4.1 % | 0 / 101 / 22 / 24 / 22 / 39 / 65 (34) |

Marginal (overlapping) shares are in the appendix (T5c): (e) alone 49–60 %; the hub and small-truth quintiles 17–27 % and 19–22 % each — but those quintiles are one fifth of the panel by construction and carry no more than a random fifth would (a random 20 % of targets carries 20 % of the deficit), which is the AUC ≈ 0.5 statement in another form. The render condition's marginal share (58 %) is the same rank-instability seen from the other side.

By bins of the mechanism covariate and of the two pipeline-independent covariates (T5d, cancer folds pooled):

| bin | base·truth: n / mean `pds` / deficit share | held-vs-donor agreement: n / `pds` / share | donor-donor agreement (observable): n / `pds` / share |
|---|---|---|---|
| ≤ 0 | 149 / 0.389 / **53.1 %** | 178 / 0.564 / **45.3 %** | 210 / 0.707 / 35.9 % |
| (0, 0.05] | 564 / 0.858 / **46.8 %** | 394 / 0.835 / 37.8 % | 348 / 0.808 / 38.9 % |
| (0.05, 0.1] | 108 / 0.998 / 0.1 % | 144 / 0.956 / 3.7 % | 145 / 0.900 / 8.5 % |
| > 0.1 | 32 / 1.000 / 0.0 % | 61 / 0.965 / 1.2 % | 74 / 0.89 / 4.7 % |

Cross-tab, held-vs-donor (rows) × donor-donor (columns), deficit share / n / mean `pds`: held ≤ 0 & donors ≤ 0 **20.1 % / 65 / 0.47**; held ≤ 0 & donors (0, 0.05] 16.7 % / 74 / 0.61; held ≤ 0 & donors > 0.05 **8.4 % / 39 / 0.63** (the donors agree with each other and not with the held line: a line-specific response); held (0, 0.05] & donors ≤ 0 13.9 % / 115 / 0.79; held (0, 0.05] & donors (0, 0.05] 19.8 % / 202 / 0.83; held > 0.05: 4.9 % of the deficit on 205 targets at 0.89–0.99.

### 5.3 The components on the same targets (T7 and `work/extra262.log` §3; cancer folds pooled, mean `pds`)

| rows | n | `BASE` | `CD4ONLY` | `UNIONONLY` | `h1union` | `cd4svd5` | `cd4svd5` − `BASE` |
|---|---:|---:|---:|---:|---:|---:|---:|
| all | 853 | 0.779 | 0.796 | 0.783 | 0.787 | 0.799 | +0.020 |
| group (e), base·truth ≤ 0 | 149 | 0.291 | 0.371 | 0.315 | **0.403** | 0.389 | +0.097 |
| the marginal band, (0, 0.05] | 564 | 0.853 | 0.858 | 0.853 | **0.835** | 0.858 | +0.005 |
| base·truth > 0.05 | 140 | 0.999 | 0.999 | 0.998 | 0.999 | 0.999 | −0.000 |
| group (e) ∧ union carried | 22 | 0.279 | 0.346 | **0.441** | 0.473 | 0.468 | +0.189 |
| group (e) ∧ union absent | 127 | 0.294 | 0.375 | 0.294 | 0.391 | 0.375 | +0.081 |
| marginal band ∧ union carried | 144 | 0.877 | 0.893 | 0.878 | 0.875 | 0.892 | +0.015 |
| marginal band ∧ union absent | 420 | 0.845 | 0.846 | 0.845 | 0.821 | 0.846 | +0.001 |
| donor-donor agreement ≤ 0 | 210 | 0.678 | 0.705 | 0.680 | 0.689 | 0.707 | +0.029 |
| donor-donor agreement > 0.05 | 219 | 0.886 | 0.896 | 0.887 | 0.883 | 0.897 | +0.011 |

Three readings. (i) The residual stack is a modest repair on a base it inherits: 52–61 % of `cd4svd5`'s deficit is on targets already below 0.5 in `BASE`; it lifts 17 / 12 / 17 of 51 / 48 / 55 such targets above 0.5 and drops 0 / 0 / 2 of ~150 targets that were ≥ 0.9. (ii) **The union's entire `pds` is the anti-aligned targets it carries**: +0.16 per target on 22 of them (0.279 → 0.441) and nothing on the 144 marginal ones it carries (0.877 → 0.878) — a direct H1 measurement replaces a wrong direction but cannot sharpen a nearly-right one. (iii) **The deflation's gain lived on the marginal band and its cost on the anti-aligned group**: the un-deflated CD4 rows were a shared component that pushed the marginal targets' own truths down the ranking (`h1union` 0.835 < `BASE` 0.853 on 564 targets; deflated 0.858), while on group (e) the un-deflated CD4 was worth more (0.403 vs 0.389) — the removed component carried line-specific signal there, which is §6's reading.

## 6. The `h1union` → `cd4svd5` step, per target (T6)

| fold | mean Δ`pds` | gainers / losers / unchanged | Σ gains | Σ losses | top-3 share of the net | ρ(Δ, CD4·truth before) | ρ(Δ, base·truth) | ρ(Δ, energy removed) | ρ(Δ, ‖CD4‖) | ρ(Δ, `pds` before) |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | +0.00999 | 102 / 71 / 96 | +6.51 | −3.83 | 55 % | **−0.364** | −0.078 | +0.019 | +0.029 | −0.197 |
| `orion_hct116` | +0.00945 | 112 / 88 / 90 | +8.58 | −5.84 | 55 % | **−0.360** | +0.049 | +0.103 | +0.052 | −0.141 |
| `orion_hek293t` | +0.01736 | 134 / 83 / 77 | +10.10 | −4.99 | 33 % | **−0.382** | −0.049 | +0.157 | +0.154 | −0.238 |
| `kolf2` (apart) | +0.00414 | 113 / 107 / 53 | +7.15 | −6.02 | 76 % | **−0.352** | +0.101 | −0.016 | −0.016 | −0.089 |

- **The per-target Δ does not repeat across held lines**: Spearman of Δ between cancer folds **+0.033 / −0.026 / −0.101** on 260–289 common targets. The biggest movers flip sign: `PLCG1` +0.394 / −0.617 / +0.661 (KOLF2 −0.537), `PGLS` +0.498 / −0.531 / +0.093, `LMAN1` +0.402 / −0.025 / +0.324, `MYD88` −0.157 / +0.255 / −0.265, `SMARCC2` +0.040 / −0.375 / +0.207, `SFXN4` +0.005 / −0.042 / −0.538, `ZNF530` −0.088 / −0.072 / +0.512 (full table in the appendix). Of 27 top movers, 6 keep one sign on all three folds.
- **The losers share exactly one covariate and it is unobservable**: the CD4 row's cosine to the held truth *before* deflation (ρ −0.35…−0.38 on every fold, incl. `kolf2`): the projector removed, from those rows, a component that pointed at that line's response. No structural covariate separates losers from gainers (AUC losers-higher: energy removed 0.47, ‖CD4‖ 0.44, whitened support 0.49, union 0.43, ‖truth‖ 0.49; pooled medians losers / gainers / unchanged: base·truth 0.008 / 0.012 / 0.048; energy removed 0.088 / 0.106 / 0.090; CD4·truth before 0.018 / −0.001 / 0.056; `h1union` `pds` 0.817 / 0.793 / 1.000). Both losers and gainers are the marginal-band targets; the unchanged are the saturated ones.
- **The reading for the next arm**: the deflation is a re-ranking of the marginal band with a line-specific sign per target and a positive mean, not a per-target repair. Any further re-shaping of the CD4 rows (k = 10, direction subsets, softer s — 255/257) re-rolls the same band; its per-target outcome cannot be validated fold by fold, only its mean can, and the mean's transfer (k_pds 1.22 on A/B/C) is a property of those three lines' marginal bands. The direction the next arm must not repeat: a change whose per-target sign is not reproducible across held lines cannot be priced by retention (leave-top-3-out) and cannot be attributed by covariate — 257 §4.3's "undefined retention" and 261 §6's negative retention are this fact seen from the price side.

## 7. Pilots (labelled as such: one variable each, KA row, one control; NOT the screens of §8, NOT priced)

### 7.1 (a) Agreement-keyed CD4 row weight (`pilot262.py`, three cancer folds; `work/pilot262.{json,log}`)

The deflated CD4 row of target p multiplied by w_p = exp(β(0.5 − r_p)), r_p the within-fold rank in [0, 1] of the observable donor-donor agreement (0 = most disagreeing; rows without the covariate w = 1), the object re-solved by `solve_fraction` at f = 0.34 and merged with the union at its shipped fraction (261's construction with one multiplication). β = 0 == `cd4svd5` at 0.0e+00 (3/3). Controls: the same weight multiset on deranged rows (seed 20726); the reverse key −β.

| cell | Δ`pds` vs `cd4svd5` (3-fold mean; per fold) | folds + | rel `cos_agg` | rel macro |
|---|---|---:|---:|---:|
| β = 1 | **+0.000448** (+0.00042 / +0.00008 / +0.00085) | 3/3 | **−3.46 %** | −1.35 % |
| β = 1, shuffled weights | −0.000142 (+0.00032 / +0.00074 / −0.00149) | 2/3 | +1.40 % | +0.17 % |
| β = 1, reverse key | −0.001397 (−0.00133 / −0.00098 / −0.00188) | 0/3 | +4.19 % | +0.80 % |
| β = 2 | −0.000718 | 1/3 | −6.81 % | −2.74 % |
| β = 2, shuffled / reverse | −0.001591 / −0.003164 | 1/3 / 0/3 | +3.16 % / +8.93 % | −0.12 % / +0.49 % |
| β = 4 | −0.008052 | 0/3 | −20.2 % | −6.39 % |
| β = 4, shuffled / reverse | −0.008773 / −0.009018 | 0/3 / 0/3 | −0.03 % / +12.1 % | −3.20 % / −3.93 % |

The key has the right sign (the reverse loses on 3/3) and no value: +0.00045 is 0.3 LOO-σ (0.001494), and it costs 3.5 % of the aggregate cosine — 0.00013 `pds` per −1 % against the deflation's own 0.0266 per −1 % (246 §3.5), 200× worse. Mechanism, read off §5.3: the disagreeing rows are rows where CD4 is unreliable too (CD4·truth −0.002 on the bottom decile), so up-weighting CD4 there adds noise, and down-weighting it on agreeing rows removes cosine that was real. **Negative; the screen in §8.1 is recorded for completeness and is not worth a fold load unless a better observable appears.**

### 7.2 (b) Headroom for choosing the donor per target (`pilot262.py` §b; from the arrays on disk, no new fold)

For each cancer fold: each donor's agreement with the held truth (cross-line truth cosine), the mean of the two (proxy for the base's agreement: Spearman with base·truth **+0.84 / +0.70 / +0.72**), the best single donor (oracle), the donor CD4 would pick (the one whose row CD4 agrees with more — CD4·(donor row) is `cd4_cos_truth` measured on that donor's own fold), and a random pick.

| fold, rows | n | mean of two | **best single (oracle)** | **CD4-arbiter pick** | random pick | > 0.05: mean / best / arbiter | arbiter = best |
|---|---:|---:|---:|---:|---:|---|---:|
| K562, group (e) | 40 | −0.0080 | +0.0091 | −0.0062 | −0.0066 | 0 % / 8 % / 0 % | 62 % |
| HCT116, group (e) | 35 | −0.0310 | +0.0003 | −0.0367 | −0.0346 | 0 % / 3 % / 3 % | 43 % |
| HEK293T, group (e) | 42 | −0.0243 | +0.0153 | −0.0169 | −0.0235 | 0 % / 19 % / 12 % | 62 % |
| K562, marginal band | 141 | +0.0188 | +0.0327 | +0.0187 | +0.0174 | 6 % / 20 % / 9 % | 48 % |
| HCT116, marginal band | 168 | +0.0301 | +0.0618 | +0.0275 | +0.0217 | 24 % / 46 % / 24 % | 52 % |
| HEK293T, marginal band | 171 | +0.0331 | +0.0647 | +0.0320 | +0.0307 | 26 % / 50 % / 27 % | 49 % |

An oracle that knew which donor agrees with the held line would put 20–50 % of the marginal band above the 0.05 saturation point (vs 6–26 % for the mean) and lift group (e) to ≈ 0 … +0.015 — by the conversion curve, roughly +0.1–0.2 `pds` per target there. **An independent line cannot find it: CD4 picks the better donor at coin-flip rate (43–62 %) and lands on the random pick.** Per-target donor selection needs an arbiter that knows the held line; none is on disk (basal-profile similarity was closed as a global weight by 010 and is the only other candidate).

### 7.3 The below-chance group cannot be neutralised (`work/extra262.log` §2)

Targets with own cosine ≤ 0 (16 %) score 0.295 on average; a neutral row (zero effect) scores exactly 0.5 by the tie rule (brief §1), so nulling them would be worth +0.03 `pds` at oracle. No observable selector finds a bin below 0.5: donor-donor ≤ 0 → 0.707 (n 210); ≤ −0.02 → 0.724; CD4·base ≤ 0 → 0.727; residual·base ≤ 0 → 0.704; donor-donor ≤ 0 ∧ CD4·base ≤ 0 ∧ union absent → 0.718; bottom 5 % by min(donor-donor, CD4·base) → **0.685**. Nulling any of these loses `pds`, and a null row scores 0 on `fid` and `reach` by definition (brief §4–5: a submission calling nothing scores 0 whenever the reference found anything). Dead on both counts; recorded so it is not re-proposed.

## 8. Hypotheses (at most three, one variable each; no price, nothing cleared)

**H1 — Dose the CD4 residual by observable donor disagreement.** *Mechanism:* the base is unreliable where the whitened donors disagree with each other (donor-donor agreement ρ +0.38 with `pds`; the residual's lift is +0.029 there vs +0.011 where they agree), so the un-whitened residual should carry more weight there and less where the base is already saturated (`pds` 0.999 above cosine 0.05, where CD4 can only cost cosine). *Covariate group:* (e) and the marginal band, selected by the observable donor-donor agreement. *Builder variable (named, not made):* `VCC_CD4_ROW_WEIGHT_BETA` — w_p = exp(β(0.5 − rank_p)) on the deflated CD4 rows before `solve_cd4_dose`, rank_p from the pairwise cosines of the three production donor rows (un-whitened, operator columns); unset/0 byte-inert. *Screen (≤ 2 h):* PLOLO, three cancer folds one per process, β ∈ {0.5, 1, 2} on `cd4svd5`; KA β = 0 == 257 `k5_s1` at 0.0e+00; controls: deranged weights (218's seed), reverse key, and the same β on the `h1union` base (is it additive to the deflation?); endpoints `pds`, `cos_agg`, macro, 188's `nmae` identity, leave-top-3-out. **The pilot (§7.1) has already read this at β = 1: +0.00045 `pds` (0.3 σ) for −3.5 % cosine, the reverse key −0.0014 — negative. Run the screen only with a different observable; the one on disk is exhausted.**

**H2 — Choose the donor per target where the donors disagree, with an independent arbiter.** *Mechanism:* on the marginal band the best single donor's agreement with the held line is 1.7–2× the two-donor mean's and would put 20–50 % of the band past the 0.05 saturation point; the equal-weight mean of two disagreeing rows is a compromise pointing at neither (§7.2). *Covariate group:* (e) ∧ donor-donor ≤ 0 (65 targets, 20 % of the deficit at `pds` 0.47) and the disagreeing half of the marginal band. *Builder variable (named):* `VCC_BASE_ROW_ARBITER=<source>` with `VCC_BASE_ROW_ARBITER_TAU` — where the donors' pairwise cosine is ≤ τ, replace the pre-whitening mean row by the donor row the arbiter agrees with most; unset = the mean (byte-inert). *Screen (≤ 2 h):* rebuild the two per-source prepared rows per fold (`E062.prepare_main_source`, ~1–2 min per fold), select by arbiter, whiten, compose as `cd4svd5`; KA τ = −1 == `cd4svd5`; controls: random pick, deranged arbiter, and the ORACLE pick (the held truth itself) as the ceiling; endpoints as H1. **The pilot (§7.2) has already read the only arbiter on disk: CD4 matches the best donor 43–62 % of the time and equals the random pick on every group. The screen is worth running only with a new arbiter that knows the held line; the ceiling it would test against is the oracle row.**

**H3 — Read every direction arm on the marginal band, per target, before pricing it.** *Mechanism:* `pds` is a threshold in the per-target transfer cosine (§4.2: +0.095 per +0.01 inside (0, 0.05], zero above 0.03, below chance below 0), so an arm's `pds` is the number of targets it moves across ~0.02, and its sign per target is line-specific when the arm re-shapes rows (§6: Δ Spearman across folds −0.10…+0.03 for the deflation). *Covariate group:* the marginal band (66 % of targets, 47 % of the deficit) vs the saturated block (16 %, 0 %) vs group (e) (17.5 %, 53 %). *Builder variable:* none — this is the instrument rule; the existing per-target arrays every experiment since 246 saves are the input. *Screen (≤ 1 h, no fold load):* for the arms already measured on `cd4svd5` (257's grid, 258's dose and subsets, 259's constructions, 261's projectors — all have `pertarget*.npz`), tabulate Δ`pds` by base·truth bin (this experiment's `work/diag262-<fold>.npz` `cov|base_cos_truth`) and count the targets crossing 0.02 in own cosine; test whether Σ over the marginal band predicts the arm's fold `pds` (it must, by construction) and whether the anti-aligned group's share predicts the fold-to-fold retention (the union arms' gains sat on group (e) ∧ carried, 22 targets; the deflation's on the band). KA: `cd4svd5` − `h1union` per fold reproduces +0.00999 / +0.00945 / +0.01736 from the bins. Controls: the energy-matched random and shuffled rows those experiments saved, which should show no band structure. What it buys: a pricing rule that says in advance whether an arm's `pds` is a hub effect on the band (transfers as a mean, cannot be attributed) or a coverage effect on group (e) (transfers per carried target; its October value is the number of anti-aligned targets the source carries, ~17 % of its coverage at +0.16 each), and therefore which retention condition to apply to it.

No fourth: the "null the anti-aligned rows" idea is dead by §7.3, and the hub, support, truth-norm and magnitude routes are closed by §4–5 at |ρ| ≤ 0.09.

## 9. Honest limits

- **One instrument, three cancer lines.** Every number is 162's PLOLO with deep pseudobulk truths and the held operators; the board's per-target `pds` is not observable, so "the same targets" is a statement about K562 / HCT116 / HEK293T held out, not about A/B/C (JURKAT / HELA / CAL33). The fold-consistency (Spearman +0.70…+0.76) is what licenses the transfer of the diagnosis, not of any number.
- **base·truth cosine is the mechanism, not a cause**: ρ +0.875 is partly the definition of the metric restated. The causal covariates are the pipeline-independent ones — cross-line agreement (+0.63 / +0.38) and CD4's DE count (+0.38) — and they are moderate, so a third of the per-target variance in the anti-aligned group is not explained by anything measured here (line-specific biology, or truth-side noise on a deep pseudobulk of a negligible response).
- **The partition's thresholds are stated, not fitted** (base·truth ≤ 0; residual drop > 0.1; quintiles; render > 0.05) and the group shares move with them; the bin tables (T5d) and the conversion curve are the threshold-free versions and say the same thing.
- **The essentiality flag is uninformative by construction** (0 on every fold's targets, §8.1g); CD4's DE count is the only graded effect-strength measure on disk (the `depmap26q1` extract is a KO table for other lines, not a panel-wide gene-effect matrix) — it is measured in a T cell, not a cancer line.
- **The pilots are pilots**: one seed per control, three folds, β on a rank scale; they carry no price and their negatives are 0.3–1 σ readings whose sign the reverse-key control makes credible. §7.2's arbiter test uses the cross-line truth cosines as stand-ins for the un-whitened donor rows (identical up to the panel-column zeroing and the scored-column subset; the proxy's Spearman with base·truth is +0.70…+0.84).
- **`kolf2`** is reported apart; its residual competitor-side hub on the production operator (ρ −0.18 / −0.34) is a lead about the shipped operator that this experiment cannot follow (an iPSC truth is the wrong instrument for a cancer-panel operator).
- **Peak RSS**: the diagnostic 4.8–4.9 GiB per fold; the pilot staged three folds in one process at 8.85 GB, above the mandate's ~8 GB target.

## 10. Files, memory, provenance

| file | what | peak RSS |
|---|---|---|
| `diag262.py` (`38958b10…`) → `work/diag262-<fold>.{npz,json,log}` | the per-fold measurement: five objects, per-target `pds`/cosine (delivered and effect level), beat matrices, geometry and structural covariates; KA-1/2/3 | 4.81 / 4.92 / 4.94 / 4.86 GiB, 8–10 s |
| `explore262.py` → `work/explore262.log` | the first look (prints) | — |
| `analyse262.py` (`c9b8dfea…`) → `work/analysis262.json`, `work/tables262.md` | distribution, deciles, agreement, covariate ranking, decomposition, the step; the cross-line truth agreement covariate from the fold caches | < 1 GiB |
| `extra262.py` → `work/extra262.{json,log}` | conversion curve, observable selectors, components by group, the cross-tab, top-vs-bottom medians | < 1 GiB |
| `pilot262.py` (`14891cec…`) → `work/pilot262.{json,log}` | §7.1 and §7.2 | 8.85 GB (three folds, one process) |
| `work/cd4sig262.{json,log}` | CD4's DE count per target vs `pds` (§4.3) | < 1 GiB |
| `work/plolo/` | byte copies of 257's fold caches, `SHA256SUMS.txt`; `work/.marker-start` the write-audit marker | |

Nothing outside this directory was created or modified; no builder import; no `vcc`; no `cell-eval2`; no submission.

---

## Verdict

**DIAGNOSIS DELIVERED; NOTHING PRICED, NOTHING CLEARED.** `cd4svd5` loses `pds` on a fold-consistent set of fully-supported targets whose responses do not agree across real lines (the panel's weak-response tier), where the donor consensus points away from the held line's truth and the prediction scores below chance; that group is 17.5 % of the panel and 53 % of the deficit, the marginal band (cosine 0–0.05) is the other 47 %, and above cosine 0.05 the metric is saturated. Hub, support, truth size by norm and magnitude explain nothing. The deflation's gain was a hub effect on the marginal band with a line-specific per-target sign; the union's gain is the anti-aligned targets it happens to carry. The two per-target repairs the diagnosis suggests — dosing by observable disagreement and choosing the donor with an independent arbiter — were piloted and read negative; the loss is information-limited at the selector level as well as at the source level. The rule that survives (H3) is that every direction arm should be read per target on the marginal band and on the anti-aligned group before it is priced, because those two groups transfer by different mechanisms.

---

# Appendix — the per-target tables (`work/tables262.md`, verbatim)

## T1. Per-target `pds` of `cd4svd5`, per fold

| fold (held line) | n | mean `pds` | deficit 1−mean | q10 | q25 | median | q75 | at 1.0 | < 0.5 | < 0.9 | bottom-decile share of deficit | bottom-quartile share |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` (K562) | 269 | 0.8154 | 0.1846 | 0.389 | 0.711 | 0.958 | 1.000 | 82 | 43 | 111 | 42.3 % | 78.7 % |
| `orion_hct116` (HCT116) | 290 | 0.8017 | 0.1983 | 0.382 | 0.698 | 0.939 | 1.000 | 77 | 50 | 126 | 40.4 % | 75.8 % |
| `orion_hek293t` (HEK293T) | 294 | 0.7810 | 0.2190 | 0.343 | 0.633 | 0.905 | 0.998 | 65 | 51 | 146 | 36.8 % | 70.8 % |
| `kolf2` (KOLF2) — apart | 273 | 0.6403 | 0.3597 | 0.166 | 0.402 | 0.696 | 0.952 | 34 | 90 | 195 | 25.6 % | 55.1 % |

## T2. Fold-to-fold agreement of the bottom decile (cancer folds; targets restricted to the 259 common to all three)

| pair | \|A\| | \|B\| | overlap | expected by chance | Jaccard | Spearman of per-target `pds` (all common targets) |
|---|---:|---:|---:|---:|---:|---:|
| K562 vs HCT116 | 27 | 24 | **9** | 2.5 | 0.214 | +0.756 (n = 260) |
| K562 vs HEK293T | 27 | 24 | **12** | 2.5 | 0.308 | +0.705 (n = 263) |
| HCT116 vs HEK293T | 24 | 24 | **12** | 2.2 | 0.333 | +0.741 (n = 289) |

Bottom decile on **all three** cancer folds (8): `CAT`, `LPCAT3`, `SAP30L`, `SNRK`, `TRMT1`, `ZBTB41`, `ZNF32`, `ZNF783`.  On **two or more** (17): `CAT`, `DEAF1`, `IFNAR2`, `LPCAT3`, `PDK1`, `PYCR1`, `RBM4B`, `RNF141`, `SAP30L`, `SNRK`, `TRMT1`, `VIM`, `ZBTB41`, `ZKSCAN3`, `ZNF32`, `ZNF641`, `ZNF783`.  Of those, in `kolf2`'s bottom decile: `TRMT1`, `ZBTB41`.

| target (bottom decile on ≥ 2 cancer folds) | `pds` K562 | `pds` HCT116 | `pds` HEK293T | `pds` KOLF2 | base·truth cos (mean, cancer) | held-vs-donor truth agreement (mean, cancer) | donor-donor agreement (mean) | CD4 | union carriers | whitened sources (cancer folds) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `CAT` | 0.24 | 0.25 | 0.29 | 0.68 | -0.007 | -0.005 | -0.005 | 3/3 | 0 | 2/2/2 |
| `DEAF1` | 0.06 | 0.43 | 0.23 | 0.38 | -0.020 | -0.085 | -0.085 | 3/3 | 0 | 2/2/2 |
| `IFNAR2` | 0.78 | 0.06 | 0.09 | 0.77 | -0.005 | -0.020 | -0.020 | 3/3 | 0 | 2/2/2 |
| `LPCAT3` | 0.13 | 0.34 | 0.25 | 0.86 | -0.000 | -0.004 | -0.004 | 3/3 | 0 | 2/2/2 |
| `PDK1` | 0.39 | 0.33 | 0.14 | 0.83 | -0.004 | -0.019 | -0.019 | 2/3 | 0 | 2/2/2 |
| `PYCR1` | 0.48 | 0.09 | 0.26 | 0.36 | -0.004 | -0.034 | -0.034 | 3/3 | 0 | 2/2/2 |
| `RBM4B` | 0.93 | 0.03 | 0.04 | 0.67 | -0.006 | -0.087 | -0.087 | 3/3 | 0 | 2/2/2 |
| `RNF141` | 0.31 | 0.67 | 0.33 | 0.67 | +0.002 | +0.022 | +0.022 | 3/3 | 0 | 2/2/2 |
| `SAP30L` | 0.02 | 0.03 | 0.02 | 0.46 | -0.013 | +0.005 | +0.005 | 1/3 | 0 | 2/2/2 |
| `SNRK` | 0.10 | 0.22 | 0.33 | — | -0.010 | -0.003 | -0.003 | 3/3 | 1 | 2/2/2 |
| `TRMT1` | 0.18 | 0.32 | 0.34 | 0.04 | -0.008 | +0.002 | +0.002 | 3/3 | 0 | 2/2/2 |
| `VIM` | 0.15 | 0.82 | 0.30 | 0.22 | -0.007 | -0.021 | -0.021 | 3/3 | 0 | 2/2/2 |
| `ZBTB41` | 0.27 | 0.10 | 0.06 | 0.14 | -0.014 | -0.028 | -0.028 | 3/3 | 0 | 2/2/2 |
| `ZKSCAN3` | 0.29 | 0.53 | 0.27 | 0.92 | -0.003 | -0.019 | -0.019 | 3/3 | 0 | 2/2/2 |
| `ZNF32` | 0.26 | 0.09 | 0.22 | 0.18 | -0.017 | -0.050 | -0.050 | 3/3 | 1 | 2/2/2 |
| `ZNF641` | 0.17 | 0.27 | 0.67 | 0.43 | -0.002 | -0.004 | -0.004 | 3/3 | 0 | 2/2/2 |
| `ZNF783` | 0.20 | 0.10 | 0.15 | 0.39 | -0.009 | -0.025 | -0.025 | 3/3 | 0 | 2/2/2 |

Top decile (all at `pds` = 1.0 on every fold; the 1.0 block is 82, 77, 65 targets on the three cancer folds) shared by all three cancer folds: 14 of the 27-target deciles — `RBM18`, `RBM26`, `RBM5`, `RNF10`, `SMARCA5`, `SYVN1`, `TBCK`, `THAP7`, `UPF3B`, `WBP4`, `ZDHHC4`, `ZFP62`, `ZNF462`, `ZNF507`.

## T3b. `k562_gwps` (K562) — bottom decile (27 targets)

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `SAP30L` | 0.02 | 0.00 | 0.02 | -0.023 | -0.016 | — | -0.016 | 0.048 | 7.0 | 0.06 | 0.03 | -0.000 | 2 | 1 | 0 | 0.0 | — |
| `DEAF1` | 0.06 | -0.09 | 0.02 | -0.022 | -0.031 | -0.000 | -0.036 | -0.181 | 8.3 | 0.07 | 0.04 | -0.001 | 2 | 3 | 0 | 5.1 | 0.08 |
| `PIN1` | 0.06 | 0.00 | 0.05 | -0.015 | -0.016 | 0.005 | -0.001 | 0.085 | 13.8 | 0.05 | 0.03 | 0.004 | 2 | 2 | 1 | 0.0 | — |
| `TMF1` | 0.07 | 0.00 | 0.07 | -0.018 | -0.013 | — | -0.020 | 0.263 | 7.0 | 0.07 | 0.03 | 0.004 | 2 | 2 | 0 | 0.0 | — |
| `SNRK` | 0.10 | -0.05 | 0.11 | -0.017 | -0.016 | -0.012 | -0.001 | -0.006 | 5.9 | 0.10 | 0.03 | 0.001 | 2 | 3 | 1 | 5.2 | 0.29 |
| `SNX4` | 0.11 | -0.02 | 0.14 | -0.019 | -0.015 | -0.015 | -0.028 | 0.016 | 9.7 | 0.09 | 0.05 | -0.002 | 2 | 3 | 0 | 5.3 | 0.14 |
| `LPCAT3` | 0.13 | -0.34 | 0.21 | -0.012 | -0.001 | -0.015 | -0.005 | -0.003 | 7.3 | 0.07 | 0.03 | 0.002 | 2 | 3 | 0 | 4.8 | 0.44 |
| `VIM` | 0.15 | 0.01 | 0.07 | -0.012 | -0.020 | 0.021 | -0.026 | -0.011 | 8.5 | 0.05 | 0.05 | 0.002 | 2 | 3 | 0 | 4.4 | 0.04 |
| `ZNF641` | 0.17 | -0.11 | 0.16 | -0.013 | -0.006 | 0.002 | -0.001 | -0.009 | 7.1 | 0.07 | 0.04 | -0.001 | 2 | 3 | 0 | 6.4 | 0.20 |
| `SMARCC2` | 0.17 | 0.04 | 0.08 | -0.013 | -0.014 | 0.003 | -0.007 | -0.028 | 5.1 | 0.05 | 0.04 | 0.000 | 2 | 3 | 0 | 5.7 | 0.20 |
| `TRMT1` | 0.18 | -0.02 | 0.13 | -0.011 | -0.021 | 0.007 | -0.009 | 0.024 | 6.4 | 0.08 | 0.04 | 0.000 | 2 | 3 | 0 | 3.8 | 0.01 |
| `ZNF783` | 0.20 | 0.05 | 0.13 | -0.011 | -0.013 | 0.018 | -0.032 | -0.012 | 6.5 | 0.07 | 0.04 | -0.001 | 2 | 3 | 0 | 3.9 | 0.05 |
| `COMMD8` | 0.21 | 0.09 | 0.15 | -0.009 | -0.016 | 0.014 | 0.017 | 0.040 | 6.9 | 0.10 | 0.04 | 0.001 | 2 | 3 | 0 | 5.1 | 0.13 |
| `CAT` | 0.24 | -0.15 | 0.31 | -0.006 | -0.002 | 0.000 | 0.014 | -0.043 | 5.8 | 0.07 | 0.04 | 0.003 | 2 | 3 | 0 | 3.7 | 0.06 |
| `SLC39A1` | 0.25 | 0.01 | 0.30 | -0.009 | -0.011 | -0.007 | -0.004 | 0.106 | 6.0 | 0.05 | 0.03 | -0.001 | 2 | 3 | 0 | 3.5 | 0.05 |
| `ZNF32` | 0.26 | 0.06 | 0.15 | -0.008 | -0.016 | 0.016 | -0.010 | -0.129 | 8.1 | 0.09 | 0.03 | -0.000 | 2 | 3 | 1 | 4.3 | 0.03 |
| `ZBTB41` | 0.27 | -0.06 | 0.22 | -0.006 | -0.012 | 0.012 | -0.017 | -0.048 | 6.8 | 0.05 | 0.03 | 0.002 | 2 | 3 | 0 | 3.2 | 0.05 |
| `ZKSCAN3` | 0.29 | 0.05 | 0.21 | -0.003 | -0.002 | 0.007 | -0.001 | -0.055 | 8.6 | 0.05 | 0.04 | 0.003 | 2 | 3 | 0 | 5.4 | 0.08 |
| `TRIM35` | 0.31 | 0.01 | 0.30 | -0.008 | -0.007 | -0.003 | -0.020 | 0.073 | 4.5 | 0.07 | 0.03 | -0.002 | 2 | 3 | 0 | 3.7 | 0.05 |
| `FZD2` | 0.31 | 0.00 | 0.31 | -0.006 | -0.008 | — | 0.004 | -0.012 | 6.0 | 0.05 | 0.03 | 0.000 | 2 | 0 | 0 | 0.0 | — |
| `RNF141` | 0.31 | 0.02 | 0.27 | -0.006 | -0.007 | -0.001 | -0.006 | 0.078 | 4.9 | 0.04 | 0.03 | 0.000 | 2 | 3 | 0 | 3.4 | 0.02 |
| `YIPF6` | 0.32 | -0.04 | 0.46 | -0.007 | 0.002 | -0.010 | -0.010 | 0.073 | 9.5 | 0.05 | 0.04 | -0.001 | 2 | 3 | 0 | 3.6 | 0.03 |
| `SH2D3A` | 0.32 | 0.16 | 0.33 | -0.008 | -0.012 | -0.003 | -0.066 | 0.092 | 6.6 | 0.10 | 0.04 | -0.002 | 2 | 3 | 0 | 3.4 | 0.06 |
| `MAP2K2` | 0.36 | -0.02 | 0.59 | -0.004 | 0.012 | -0.032 | 0.014 | 0.108 | 6.4 | 0.06 | 0.04 | 0.001 | 2 | 3 | 0 | 3.7 | 0.06 |
| `RAB5A` | 0.36 | -0.10 | 0.19 | -0.003 | -0.007 | 0.027 | 0.012 | 0.147 | 6.7 | 0.07 | 0.04 | 0.003 | 2 | 3 | 0 | 4.0 | 0.12 |
| `NIPSNAP3A` | 0.37 | 0.00 | 0.35 | -0.003 | -0.007 | 0.004 | -0.011 | 0.137 | 8.0 | 0.08 | 0.03 | 0.001 | 2 | 3 | 0 | 4.8 | 0.11 |
| `DVL3` | 0.38 | 0.19 | 0.40 | -0.002 | 0.007 | -0.006 | 0.005 | -0.008 | 4.6 | 0.08 | 0.03 | 0.002 | 2 | 3 | 0 | 4.0 | 0.18 |
| **median, bottom decile** | 0.24 | 0.00 | 0.19 | -0.009 | -0.012 | 0.001 | -0.007 | 0.016 | 6.8 | 0.07 | 0.04 | 0.000 | 2 | 3 | 0 | 3.9 | 0.06 |
| median, the rest | 0.98 | 0.00 | 0.96 | 0.027 | 0.024 | 0.025 | 0.025 | 0.030 | 6.8 | 0.08 | 0.04 | 0.000 | 2 | 3 | 0 | 4.5 | 0.10 |

## T3t. `k562_gwps` (K562) — top decile (27 targets)

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *27 targets, all at 1.0:* `ZNF628`, `ZNF526`, `ZNF507`, `ZNF462`, `ZNF444`, `ZKSCAN1`, `ZFP62`, `ZDHHC4`, `WBP4`, `UPF3B`, `TRAPPC6A`, `THAP7`, `TEX2`, `TBCK`, `SYVN1`, `STK3`, `SRSF5`, `SND1`, `SMARCA5`, `SLIRP`, `SIN3B`, `SCAF8`, `RSRC1`, `RNF10`, `RBM5`, `RBM26`, `RBM18` |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| **median, top decile** | 1.00 | 0.00 | 1.00 | 0.070 | 0.059 | 0.081 | 0.056 | 0.060 | 6.5 | 0.09 | 0.04 | 0.000 | 2 | 3 | 0 | 4.6 | 0.09 |
| median, the rest | 0.93 | 0.00 | 0.90 | 0.020 | 0.018 | 0.017 | 0.019 | 0.023 | 6.8 | 0.08 | 0.04 | 0.000 | 2 | 3 | 0 | 4.5 | 0.10 |

## T3b. `orion_hct116` (HCT116) — bottom decile (29 targets)

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `RBM4B` | 0.03 | -0.01 | 0.07 | -0.027 | -0.021 | -0.014 | -0.141 | 0.022 | 7.2 | 0.25 | 0.04 | -0.003 | 2 | 3 | 0 | 3.0 | 0.03 |
| `SAP30L` | 0.03 | 0.00 | 0.03 | -0.018 | -0.009 | — | 0.014 | -0.013 | 7.2 | 0.12 | 0.03 | 0.001 | 2 | 1 | 0 | 0.0 | — |
| `IFNAR2` | 0.06 | -0.01 | 0.10 | -0.021 | -0.009 | -0.030 | -0.025 | -0.010 | 5.8 | 0.26 | 0.04 | 0.001 | 2 | 3 | 0 | 5.2 | 0.06 |
| `ZNF543` | 0.06 | -0.02 | 0.12 | -0.016 | -0.011 | -0.018 | -0.027 | -0.002 | 9.3 | 0.28 | 0.05 | -0.002 | 2 | 3 | 0 | 4.0 | 0.08 |
| `ZNF496` | 0.08 | -0.04 | 0.07 | -0.021 | -0.024 | 0.016 | -0.119 | 0.119 | 10.5 | 0.22 | 0.03 | -0.000 | 2 | 3 | 0 | 3.5 | 0.11 |
| `PGLS` | 0.09 | -0.53 | 0.09 | -0.029 | -0.016 | -0.027 | -0.068 | -0.013 | 15.8 | 0.35 | 0.04 | 0.001 | 2 | 3 | 0 | 12.2 | 0.28 |
| `ZNF32` | 0.09 | -0.02 | 0.04 | -0.013 | -0.022 | 0.027 | -0.101 | 0.052 | 12.4 | 0.20 | 0.03 | -0.000 | 2 | 3 | 1 | 4.3 | 0.03 |
| `PYCR1` | 0.09 | 0.04 | 0.57 | -0.020 | 0.006 | -0.051 | -0.046 | -0.008 | 9.7 | 0.38 | 0.03 | 0.001 | 2 | 3 | 0 | 4.7 | 0.05 |
| `ZNF783` | 0.10 | -0.04 | 0.15 | -0.017 | -0.011 | -0.002 | -0.025 | -0.026 | 6.1 | 0.09 | 0.04 | -0.001 | 2 | 3 | 0 | 3.9 | 0.05 |
| `ZBTB41` | 0.10 | 0.00 | 0.11 | -0.017 | -0.010 | -0.015 | -0.037 | -0.008 | 6.6 | 0.32 | 0.04 | 0.000 | 2 | 3 | 0 | 3.2 | 0.06 |
| `RAB11FIP5` | 0.13 | -0.05 | 0.07 | -0.017 | -0.029 | 0.037 | — | — | 14.2 | 0.50 | 0.03 | 0.001 | 1 | 3 | 0 | 4.5 | 0.14 |
| `PLCG1` | 0.14 | -0.62 | 0.41 | -0.017 | -0.002 | -0.023 | -0.013 | 0.028 | 7.9 | 0.23 | 0.07 | 0.002 | 2 | 3 | 0 | 13.8 | 0.77 |
| `PHF19` | 0.16 | -0.11 | 0.28 | -0.034 | -0.011 | -0.039 | — | — | 4.7 | 0.15 | 0.04 | -0.001 | 1 | 3 | 0 | 8.9 | 0.38 |
| `ZNF786` | 0.19 | 0.11 | 0.10 | -0.010 | -0.019 | 0.026 | -0.121 | 0.085 | 10.8 | 0.40 | 0.03 | -0.000 | 2 | 3 | 0 | 4.4 | 0.20 |
| `NICN1` | 0.21 | 0.00 | 0.21 | -0.008 | -0.010 | — | — | — | 7.1 | 0.11 | 0.03 | 0.000 | 1 | 0 | 0 | 0.0 | — |
| `DTNBP1` | 0.21 | 0.10 | 0.10 | -0.013 | -0.015 | -0.000 | -0.113 | 0.078 | 7.8 | 0.30 | 0.04 | -0.000 | 2 | 3 | 0 | 4.1 | 0.03 |
| `SNRK` | 0.22 | 0.01 | 0.14 | -0.015 | -0.007 | -0.007 | 0.002 | -0.012 | 6.9 | 0.07 | 0.04 | 0.001 | 2 | 3 | 1 | 5.2 | 0.30 |
| `CAT` | 0.25 | -0.08 | 0.16 | -0.007 | -0.010 | -0.002 | -0.015 | 0.014 | 4.2 | 0.21 | 0.04 | -0.000 | 2 | 3 | 0 | 3.7 | 0.06 |
| `BCL7C` | 0.26 | -0.34 | 0.24 | -0.007 | -0.010 | -0.007 | 0.008 | 0.031 | 7.8 | 0.10 | 0.03 | 0.000 | 2 | 3 | 0 | 4.8 | 0.19 |
| `ZNF641` | 0.27 | -0.07 | 0.29 | -0.009 | -0.001 | -0.006 | -0.008 | 0.005 | 8.0 | 0.09 | 0.05 | -0.000 | 2 | 3 | 0 | 6.4 | 0.20 |
| `ZNF219` | 0.28 | 0.12 | 0.23 | -0.008 | -0.016 | 0.005 | -0.057 | 0.052 | 13.4 | 0.17 | 0.03 | 0.002 | 2 | 3 | 1 | 7.4 | 0.23 |
| `SNN` | 0.30 | 0.05 | 0.33 | -0.008 | -0.004 | -0.006 | — | — | 9.7 | 0.08 | 0.03 | -0.001 | 1 | 3 | 0 | 3.9 | 0.09 |
| `POLD4` | 0.30 | 0.01 | 0.69 | -0.007 | 0.010 | -0.036 | -0.051 | 0.012 | 4.5 | 0.19 | 0.04 | 0.000 | 2 | 3 | 0 | 3.4 | 0.05 |
| `TRMT1` | 0.32 | -0.05 | 0.50 | -0.004 | -0.003 | -0.012 | 0.002 | 0.002 | 9.3 | 0.07 | 0.04 | 0.000 | 2 | 3 | 0 | 3.8 | 0.02 |
| `PDK1` | 0.33 | 0.00 | 0.33 | -0.004 | -0.004 | — | -0.027 | -0.002 | 8.1 | 0.17 | 0.03 | 0.001 | 2 | 2 | 0 | 0.0 | — |
| `LPCAT3` | 0.34 | -0.02 | 0.25 | -0.004 | -0.002 | 0.003 | -0.008 | 0.002 | 6.7 | 0.07 | 0.04 | 0.001 | 2 | 3 | 0 | 4.8 | 0.44 |
| `TARBP2` | 0.36 | 0.04 | 0.24 | -0.003 | -0.004 | 0.005 | -0.000 | 0.026 | 5.2 | 0.07 | 0.04 | 0.000 | 2 | 3 | 1 | 5.4 | 0.06 |
| `RFK` | 0.37 | -0.06 | 0.28 | -0.004 | -0.009 | 0.018 | — | — | 17.2 | 0.09 | 0.05 | -0.000 | 1 | 3 | 0 | 5.2 | 0.10 |
| `IREB2` | 0.38 | -0.00 | 0.44 | -0.002 | -0.004 | 0.002 | -0.013 | 0.003 | 10.2 | 0.17 | 0.04 | 0.001 | 2 | 3 | 0 | 5.8 | 0.26 |
| **median, bottom decile** | 0.21 | -0.01 | 0.21 | -0.013 | -0.010 | -0.006 | -0.026 | 0.004 | 7.9 | 0.17 | 0.04 | 0.000 | 2 | 3 | 0 | 4.4 | 0.10 |
| median, the rest | 0.97 | 0.00 | 0.96 | 0.022 | 0.020 | 0.017 | 0.028 | 0.023 | 7.0 | 0.13 | 0.04 | -0.000 | 2 | 3 | 0 | 4.5 | 0.10 |

## T3t. `orion_hct116` (HCT116) — top decile (29 targets)

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *29 targets, all at 1.0:* `ZNF507`, `ZNF462`, `ZNF410`, `ZFP62`, `ZDHHC4`, `WBP4`, `UPF3B`, `TRAPPC6A`, `TMEM79`, `TMEM104`, `THAP7`, `TEX2`, `TBCK`, `TAF4`, `SYVN1`, `SSR4`, `SMARCA5`, `SLC25A13`, `SIN3B`, `SCAF8`, `RSRC1`, `RNF2`, `RNF10`, `RBM5`, `RBM26`, `RBM18`, `R3HCC1L`, `PRMT2`, `PRKAB2` |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| **median, top decile** | 1.00 | 0.00 | 1.00 | 0.053 | 0.050 | 0.051 | 0.062 | 0.049 | 6.7 | 0.10 | 0.04 | -0.000 | 2 | 3 | 0 | 4.5 | 0.07 |
| median, the rest | 0.91 | 0.00 | 0.88 | 0.016 | 0.014 | 0.013 | 0.021 | 0.019 | 7.1 | 0.14 | 0.04 | 0.000 | 2 | 3 | 0 | 4.5 | 0.10 |

## T3b. `orion_hek293t` (HEK293T) — bottom decile (30 targets)

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `NFE2L1` | 0.01 | -0.03 | 0.07 | -0.055 | -0.012 | -0.074 | -0.098 | 0.043 | 9.3 | 0.56 | 0.04 | 0.003 | 2 | 3 | 1 | 8.3 | 0.33 |
| `SAP30L` | 0.02 | 0.00 | 0.02 | -0.020 | -0.015 | — | 0.018 | -0.020 | 5.7 | 0.25 | 0.03 | -0.001 | 2 | 1 | 0 | 0.0 | — |
| `RBM4B` | 0.04 | -0.00 | 0.04 | -0.025 | -0.016 | 0.005 | -0.135 | 0.010 | 7.7 | 0.66 | 0.03 | -0.000 | 2 | 3 | 0 | 3.0 | 0.02 |
| `ZBTB41` | 0.06 | -0.02 | 0.03 | -0.013 | -0.019 | 0.023 | -0.028 | -0.027 | 5.1 | 0.08 | 0.03 | -0.000 | 2 | 3 | 0 | 3.2 | 0.06 |
| `IFNAR2` | 0.09 | -0.04 | 0.13 | -0.012 | -0.011 | 0.006 | -0.041 | 0.023 | 4.8 | 0.18 | 0.04 | 0.000 | 2 | 3 | 0 | 5.2 | 0.06 |
| `IFI27L2` | 0.12 | 0.05 | 0.11 | -0.014 | -0.013 | 0.003 | -0.098 | 0.009 | 7.1 | 0.70 | 0.03 | -0.000 | 2 | 3 | 0 | 3.8 | 0.09 |
| `MYD88` | 0.13 | -0.27 | 0.14 | -0.014 | -0.003 | -0.002 | 0.004 | 0.023 | 3.9 | 0.32 | 0.04 | 0.001 | 2 | 3 | 0 | 4.3 | 0.14 |
| `PDK1` | 0.14 | 0.00 | 0.14 | -0.014 | -0.009 | — | -0.022 | -0.012 | 6.0 | 0.12 | 0.04 | -0.002 | 2 | 2 | 0 | 0.0 | — |
| `ZNF783` | 0.15 | -0.02 | 0.28 | -0.012 | -0.004 | -0.019 | -0.019 | -0.038 | 5.8 | 0.52 | 0.03 | -0.000 | 2 | 3 | 0 | 3.9 | 0.05 |
| `HEATR6` | 0.16 | -0.07 | 0.13 | -0.012 | -0.012 | 0.008 | — | — | 6.5 | 0.40 | 0.03 | 0.001 | 1 | 3 | 0 | 3.7 | 0.11 |
| `RAB11FIP5` | 0.16 | 0.01 | 0.18 | -0.035 | -0.045 | -0.020 | — | — | 11.7 | 0.45 | 0.03 | 0.000 | 1 | 3 | 0 | 4.5 | 0.14 |
| `SNN` | 0.20 | 0.01 | 0.18 | -0.011 | -0.010 | 0.010 | — | — | 8.3 | 0.35 | 0.03 | 0.001 | 1 | 3 | 0 | 3.9 | 0.09 |
| `TBC1D19` | 0.20 | -0.02 | 0.34 | -0.010 | -0.004 | -0.018 | — | — | 5.3 | 0.37 | 0.03 | 0.000 | 1 | 3 | 0 | 4.7 | 0.13 |
| `ZNF32` | 0.22 | 0.01 | 0.21 | -0.023 | -0.011 | -0.036 | -0.038 | -0.073 | 9.6 | 0.38 | 0.04 | 0.001 | 2 | 3 | 1 | 4.3 | 0.03 |
| `DEAF1` | 0.23 | 0.04 | 0.11 | -0.011 | -0.019 | 0.008 | -0.120 | -0.014 | 5.2 | 0.55 | 0.04 | -0.001 | 2 | 3 | 0 | 5.1 | 0.09 |
| `PSMB9` | 0.23 | 0.02 | 0.14 | -0.013 | -0.014 | 0.008 | — | — | 3.6 | 0.17 | 0.04 | 0.001 | 1 | 3 | 0 | 5.0 | 0.26 |
| `NT5DC1` | 0.23 | -0.03 | 0.54 | -0.008 | 0.005 | -0.036 | -0.006 | 0.005 | 11.6 | 0.28 | 0.03 | 0.000 | 2 | 3 | 0 | 3.2 | 0.03 |
| `LPCAT3` | 0.25 | -0.05 | 0.36 | -0.010 | 0.002 | -0.024 | -0.000 | -0.012 | 6.0 | 0.13 | 0.04 | -0.001 | 2 | 3 | 0 | 4.8 | 0.44 |
| `PYCR1` | 0.26 | -0.00 | 0.32 | -0.008 | -0.004 | 0.004 | -0.023 | -0.054 | 6.6 | 0.10 | 0.04 | 0.000 | 2 | 3 | 0 | 4.7 | 0.05 |
| `TDRD7` | 0.26 | 0.07 | 0.74 | -0.009 | 0.014 | -0.043 | — | — | 5.8 | 0.49 | 0.03 | -0.000 | 1 | 3 | 0 | 4.9 | 0.14 |
| `ZKSCAN3` | 0.27 | 0.05 | 0.13 | -0.009 | -0.002 | 0.009 | -0.034 | 0.010 | 4.4 | 0.27 | 0.04 | 0.001 | 2 | 3 | 0 | 5.4 | 0.08 |
| `CAT` | 0.29 | -0.05 | 0.18 | -0.014 | -0.009 | -0.001 | -0.014 | 0.013 | 3.0 | 0.09 | 0.04 | -0.000 | 2 | 3 | 0 | 3.7 | 0.06 |
| `SFXN4` | 0.29 | -0.54 | 0.69 | -0.005 | 0.011 | -0.023 | 0.040 | -0.012 | 5.6 | 0.70 | 0.05 | 0.000 | 2 | 3 | 1 | 3.7 | 0.10 |
| `VIM` | 0.30 | 0.04 | 0.35 | -0.006 | -0.007 | -0.010 | -0.032 | -0.001 | 10.9 | 0.34 | 0.03 | -0.001 | 2 | 3 | 0 | 4.4 | 0.04 |
| `TIMM50` | 0.30 | 0.00 | 0.40 | -0.007 | 0.000 | 0.005 | -0.049 | 0.006 | 7.0 | 0.71 | 0.06 | -0.001 | 2 | 2 | 1 | 0.0 | — |
| `FHL3` | 0.33 | -0.32 | 0.54 | -0.005 | -0.008 | -0.020 | 0.006 | 0.011 | 7.5 | 0.10 | 0.04 | -0.001 | 2 | 3 | 0 | 3.5 | 0.13 |
| `SNRK` | 0.33 | 0.02 | 0.18 | -0.010 | -0.009 | -0.006 | -0.009 | 0.010 | 6.1 | 0.07 | 0.04 | 0.002 | 2 | 3 | 1 | 5.2 | 0.30 |
| `LMAN1` | 0.33 | 0.32 | 0.20 | -0.010 | -0.006 | -0.007 | 0.021 | 0.023 | 5.6 | 0.15 | 0.04 | -0.001 | 2 | 3 | 0 | 6.1 | 0.19 |
| `RNF141` | 0.33 | -0.01 | 0.77 | -0.004 | 0.003 | -0.025 | 0.032 | 0.002 | 3.6 | 0.19 | 0.05 | 0.000 | 2 | 3 | 0 | 3.4 | 0.03 |
| `TRMT1` | 0.34 | -0.06 | 0.42 | -0.003 | -0.000 | -0.004 | 0.013 | -0.020 | 7.2 | 0.13 | 0.04 | 0.002 | 2 | 3 | 0 | 3.8 | 0.02 |
| **median, bottom decile** | 0.23 | -0.00 | 0.18 | -0.011 | -0.008 | -0.005 | -0.021 | 0.004 | 6.0 | 0.30 | 0.04 | 0.000 | 2 | 3 | 0 | 4.1 | 0.09 |
| median, the rest | 0.95 | 0.00 | 0.93 | 0.021 | 0.018 | 0.016 | 0.029 | 0.024 | 6.2 | 0.26 | 0.04 | 0.000 | 2 | 3 | 0 | 4.5 | 0.10 |

## T3t. `orion_hek293t` (HEK293T) — top decile (30 targets)

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *30 targets, all at 1.0:* `ZNF526`, `ZNF507`, `ZNF462`, `ZNF444`, `ZFP62`, `ZDHHC4`, `WBP4`, `UPF3B`, `TRIP4`, `TMEM79`, `TMEM104`, `THAP7`, `TBCK`, `TAF4`, `SYVN1`, `SUPT7L`, `SSR4`, `SRSF5`, `SND1`, `SMARCA5`, `SLIRP`, `RNF2`, `RNF10`, `RBM5`, `RBM26`, `RBM18`, `PRMT6`, `PQBP1`, `PHF8`, `PHF23` |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| **median, top decile** | 1.00 | 0.00 | 1.00 | 0.057 | 0.043 | 0.048 | 0.063 | 0.050 | 6.0 | 0.26 | 0.05 | 0.000 | 2 | 3 | 0 | 5.2 | 0.08 |
| median, the rest | 0.87 | 0.00 | 0.87 | 0.014 | 0.013 | 0.010 | 0.021 | 0.019 | 6.2 | 0.26 | 0.04 | 0.000 | 2 | 3 | 0 | 4.4 | 0.10 |

## T3b. `kolf2` (KOLF2) — bottom decile (28 targets) — apart

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `MBD4` | 0.02 | -0.25 | 0.13 | -0.022 | -0.020 | -0.032 | -0.042 | -0.027 | 8.9 | 0.20 | 0.04 | -0.000 | 3 | 3 | 0 | 4.5 | 0.24 |
| `ZNF768` | 0.03 | -0.01 | 0.01 | -0.022 | -0.017 | -0.003 | -0.021 | 0.019 | 10.1 | 0.04 | 0.03 | -0.000 | 3 | 3 | 0 | 5.8 | 0.15 |
| `TRMT1` | 0.04 | -0.03 | 0.19 | -0.014 | -0.005 | -0.016 | 0.014 | 0.002 | 9.3 | 0.11 | 0.04 | -0.000 | 3 | 3 | 0 | 3.8 | 0.02 |
| `MSANTD4` | 0.04 | -0.13 | 0.03 | -0.016 | -0.010 | -0.003 | -0.004 | -0.047 | 7.1 | 0.04 | 0.05 | 0.001 | 3 | 3 | 0 | 5.9 | 0.17 |
| `RABEP2` | 0.05 | 0.01 | 0.25 | -0.016 | -0.004 | -0.041 | -0.003 | 0.019 | 8.5 | 0.10 | 0.04 | 0.000 | 3 | 3 | 0 | 3.5 | 0.09 |
| `YIPF6` | 0.06 | -0.00 | 0.08 | -0.012 | -0.013 | -0.002 | -0.026 | 0.018 | 9.6 | 0.14 | 0.03 | -0.000 | 3 | 3 | 0 | 3.6 | 0.03 |
| `ZNF714` | 0.06 | -0.00 | 0.06 | -0.015 | -0.006 | -0.001 | -0.012 | 0.036 | 7.5 | 0.03 | 0.03 | 0.000 | 3 | 3 | 2 | 3.7 | 0.08 |
| `ABCD1` | 0.06 | 0.00 | 0.06 | -0.015 | -0.002 | — | — | — | 6.5 | 0.05 | 0.04 | 0.000 | 2 | 0 | 0 | 0.0 | — |
| `SHPRH` | 0.07 | -0.03 | 0.11 | -0.014 | -0.010 | -0.010 | -0.035 | 0.073 | 7.5 | 0.04 | 0.03 | -0.000 | 3 | 3 | 0 | 3.7 | 0.08 |
| `TWF2` | 0.07 | -0.01 | 0.08 | -0.013 | -0.001 | -0.009 | -0.005 | -0.000 | 7.8 | 0.06 | 0.04 | -0.000 | 3 | 3 | 1 | 2.9 | 0.02 |
| `DVL3` | 0.08 | 0.01 | 0.12 | -0.012 | -0.004 | -0.004 | 0.005 | 0.001 | 10.3 | 0.07 | 0.03 | -0.000 | 3 | 3 | 0 | 4.0 | 0.18 |
| `POLD4` | 0.09 | 0.03 | 0.12 | -0.011 | -0.009 | -0.009 | -0.009 | -0.030 | 7.8 | 0.09 | 0.04 | -0.000 | 3 | 3 | 0 | 3.4 | 0.04 |
| `GTF2IRD1` | 0.11 | -0.03 | 0.38 | -0.010 | 0.000 | -0.013 | 0.007 | 0.069 | 9.2 | 0.03 | 0.05 | 0.000 | 3 | 3 | 1 | 6.9 | 0.20 |
| `SMYD2` | 0.12 | 0.05 | 0.29 | -0.014 | -0.005 | -0.017 | -0.032 | 0.027 | 8.3 | 0.17 | 0.03 | 0.001 | 3 | 3 | 0 | 5.4 | 0.10 |
| `ZNF251` | 0.12 | 0.02 | 0.14 | -0.009 | -0.006 | 0.002 | 0.001 | -0.030 | 7.5 | 0.10 | 0.03 | -0.000 | 3 | 3 | 0 | 3.2 | 0.09 |
| `PTOV1` | 0.12 | -0.05 | 0.23 | -0.010 | -0.007 | -0.015 | 0.002 | 0.015 | 7.1 | 0.12 | 0.03 | 0.000 | 3 | 3 | 0 | 3.7 | 0.03 |
| `MARCKSL1` | 0.13 | 0.05 | 0.13 | -0.011 | -0.009 | 0.004 | -0.042 | -0.002 | 8.0 | 0.17 | 0.03 | -0.000 | 3 | 3 | 0 | 3.5 | 0.10 |
| `TFCP2` | 0.13 | 0.01 | 0.23 | -0.012 | -0.007 | -0.016 | -0.010 | 0.043 | 10.1 | 0.04 | 0.04 | 0.000 | 3 | 3 | 0 | 4.2 | 0.14 |
| `ZNF280C` | 0.13 | 0.02 | 0.27 | -0.011 | -0.008 | -0.021 | -0.003 | 0.004 | 6.9 | 0.09 | 0.04 | 0.000 | 3 | 3 | 0 | 3.4 | 0.06 |
| `CAPRIN2` | 0.13 | 0.06 | 0.05 | -0.013 | -0.010 | 0.006 | — | — | 5.8 | 0.18 | 0.04 | -0.000 | 2 | 3 | 0 | 7.3 | 0.12 |
| `MTHFSD` | 0.14 | -0.06 | 0.33 | -0.010 | -0.007 | -0.013 | -0.002 | 0.019 | 9.2 | 0.11 | 0.04 | 0.000 | 3 | 3 | 0 | 3.4 | 0.07 |
| `TMEM209` | 0.14 | 0.06 | 0.38 | -0.009 | 0.003 | -0.028 | 0.013 | 0.028 | 11.4 | 0.09 | 0.04 | 0.000 | 3 | 3 | 0 | 4.4 | 0.06 |
| `MYRF` | 0.14 | 0.05 | 0.14 | -0.009 | -0.009 | -0.002 | -0.024 | 0.058 | 9.4 | 0.04 | 0.04 | 0.000 | 3 | 3 | 0 | 4.4 | 0.22 |
| `ZBTB41` | 0.14 | 0.01 | 0.14 | -0.010 | -0.006 | -0.006 | 0.009 | -0.028 | 7.5 | 0.06 | 0.04 | 0.000 | 3 | 3 | 0 | 3.2 | 0.05 |
| `PDE6D` | 0.15 | 0.01 | 0.18 | -0.009 | -0.008 | -0.005 | -0.001 | 0.007 | 8.4 | 0.05 | 0.03 | -0.000 | 3 | 3 | 0 | 5.5 | 0.08 |
| `SND1` | 0.16 | -0.06 | 0.11 | -0.013 | -0.013 | -0.010 | -0.040 | 0.097 | 10.0 | 0.13 | 0.04 | 0.002 | 3 | 3 | 0 | 5.8 | 0.22 |
| `PHF19` | 0.17 | -0.08 | 0.14 | -0.012 | -0.011 | -0.002 | — | — | 7.8 | 0.12 | 0.05 | 0.001 | 2 | 3 | 0 | 8.9 | 0.38 |
| `SPATA20` | 0.17 | 0.00 | 0.18 | -0.013 | -0.004 | -0.001 | — | — | 7.9 | 0.07 | 0.02 | 0.000 | 1 | 3 | 0 | 4.4 | 0.06 |
| **median, bottom decile** | 0.12 | 0.00 | 0.14 | -0.012 | -0.007 | -0.009 | -0.004 | 0.016 | 8.1 | 0.09 | 0.04 | 0.000 | 3 | 3 | 0 | 3.9 | 0.09 |
| median, the rest | 0.74 | 0.00 | 0.67 | 0.006 | 0.005 | 0.009 | 0.004 | 0.023 | 8.7 | 0.07 | 0.04 | -0.000 | 3 | 3 | 0 | 4.5 | 0.10 |

## T3t. `kolf2` (KOLF2) — top decile (28 targets) — apart

| target | pds | Δ vs h1union | base pds | own cos | base·truth | resid·truth | held-vs-donor agree | donor-donor agree | ‖truth‖ | truth NN | pred NN | comp med | src | CD4 | union | ‖CD4‖ | CD4 E rm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *28 targets, all at 1.0:* `ZNF507`, `WBP4`, `TIMM50`, `THAP7`, `TAF4`, `SUPT7L`, `SMARCA5`, `SIN3B`, `SCAF8`, `RBM5`, `RBM26`, `PQBP1`, `PHF8`, `PCGF6`, `PCGF1`, `PAN2`, `MED15`, `KHSRP`, `KHDRBS1`, `HMGXB3`, `GPATCH8`, `GMEB2`, `GIGYF2`, `FOXN2`, `EZH2`, `DHX35`, `DCP2`, `CSNK2A1` |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| **median, top decile** | 1.00 | 0.00 | 1.00 | 0.061 | 0.031 | 0.071 | 0.031 | 0.065 | 9.2 | 0.08 | 0.06 | -0.000 | 3 | 3 | 1 | 6.0 | 0.11 |
| median, the rest | 0.65 | 0.00 | 0.59 | 0.004 | 0.002 | 0.006 | 0.001 | 0.019 | 8.6 | 0.07 | 0.04 | 0.000 | 3 | 3 | 0 | 4.4 | 0.09 |

## T4. Covariate ranking — structural covariates vs per-target `pds` of `cd4svd5` (three cancer folds; `kolf2` apart)

Spearman ρ with `pds` per fold and pooled by within-fold rank; AUC = P(a bottom-decile target has the HIGHER covariate value than a random other target) — 0.5 is no separation, far from 0.5 is separation (below 0.5: the bottom decile has LOWER values); pooled median of the covariate in the bottom decile vs the rest.

| covariate | ρ K562 | ρ HCT116 | ρ HEK293T | **ρ pooled** | AUC K562 | AUC HCT116 | AUC HEK293T | **AUC pooled** | bottom-decile median | rest median | ρ / AUC `kolf2` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| base·truth cos (whitened donor base vs the truth, effect level) | 0.913 | 0.885 | 0.855 | **0.875** | 0.052 | 0.031 | 0.068 | **0.050** | -0.009 | 0.020 | 0.79 / 0.12 |
| cross-line truth agreement: held line vs the fold's donor lines (mean cos) | 0.769 | 0.565 | 0.567 | **0.630** | 0.120 | 0.082 | 0.163 | **0.121** | -0.014 | 0.027 | 0.51 / 0.27 |
| cross-line truth agreement between the fold's donor lines (cos; observable at build time) | 0.219 | 0.402 | 0.527 | **0.381** | 0.473 | 0.393 | 0.236 | **0.371** | 0.005 | 0.024 | 0.27 / 0.40 |
| residual·truth cos (fr·CD4_svd5 + g·union vs the truth) | 0.619 | 0.609 | 0.618 | **0.614** | 0.237 | 0.220 | 0.207 | **0.220** | -0.002 | 0.019 | 0.73 / 0.13 |
| deflated CD4 row·truth cos | 0.577 | 0.574 | 0.553 | **0.566** | 0.254 | 0.224 | 0.250 | **0.242** | -0.002 | 0.016 | 0.65 / 0.15 |
| union row·truth cos (carried rows only) | 0.569 | 0.625 | 0.448 | **0.519** | 0.164 | 0.264 | 0.192 | **0.209** | 0.003 | 0.030 | 0.80 / 0.06 |
| residual·base cos (does the residual agree with the base?) | 0.586 | 0.579 | 0.514 | **0.558** | 0.292 | 0.288 | 0.270 | **0.283** | 0.004 | 0.016 | 0.53 / 0.36 |
| residual row norm (scored cols) | 0.435 | 0.344 | 0.395 | **0.387** | 0.314 | 0.441 | 0.371 | **0.377** | 1.5 | 1.7 | 0.35 / 0.38 |
| residual share of the row, ‖R‖/‖W+R‖ | 0.340 | 0.315 | 0.345 | **0.329** | 0.341 | 0.431 | 0.397 | **0.391** | 0.299 | 0.345 | 0.24 / 0.46 |
| CD4 row norm (complete rows; 0 if absent) | 0.371 | 0.294 | 0.328 | **0.320** | 0.356 | 0.482 | 0.394 | **0.411** | 4.1 | 4.5 | 0.28 / 0.40 |
| CD4 row energy removed by the k=5 deflation | 0.129 | 0.038 | 0.014 | **0.052** | 0.432 | 0.537 | 0.460 | **0.478** | 0.085 | 0.100 | 0.11 / 0.49 |
| union row norm (0 if not carried) | 0.302 | 0.208 | 0.229 | **0.226** | 0.406 | 0.430 | 0.445 | **0.384** | 0.0 | 0.0 | 0.38 / 0.41 |
| union carriers (0/1/2) | 0.296 | 0.206 | 0.219 | **0.213** | 0.407 | 0.431 | 0.447 | **0.384** | 0.000 | 0.000 | 0.37 / 0.42 |
| whitened sources carrying the target (0/1/2; 3 on kolf2) | 0.120 | 0.104 | 0.086 | **0.068** | 0.521 | 0.464 | 0.455 | **0.439** | 2.000 | 2.000 | 0.07 / 0.49 |
| PerturbMatch-blended row (K562 source only) | — | 0.149 | 0.097 | **0.079** | 0.500 | 0.425 | 0.469 | **0.428** | 0.500 | 1.000 | 0.15 / 0.47 |
| CD4 complete (3/3 conditions) | 0.088 | 0.042 | 0.084 | **0.074** | 0.455 | 0.487 | 0.488 | **0.397** | 1.000 | 1.000 | -0.04 / 0.52 |
| CD4 condition count (0–3) | 0.089 | 0.044 | 0.085 | **0.066** | 0.455 | 0.486 | 0.489 | **0.408** | 3.000 | 3.000 | -0.04 / 0.52 |
| whitened base row norm | 0.004 | -0.017 | -0.052 | **-0.019** | 0.469 | 0.544 | 0.494 | **0.503** | 4.7 | 4.6 | 0.06 / 0.42 |
| k562_essential flag | — | — | — | **0.041** | 0.500 | 0.500 | 0.500 | **0.434** | 0.000 | 0.000 | — / 0.50 |
| truth ‖Δ‖ (scored cols) | 0.070 | -0.054 | -0.087 | **-0.029** | 0.488 | 0.594 | 0.470 | **0.518** | 7.0 | 6.7 | 0.09 / 0.40 |
| truth's nearest-neighbour cosine to another truth | 0.304 | -0.039 | 0.009 | **0.086** | 0.330 | 0.615 | 0.537 | **0.498** | 0.119 | 0.124 | 0.12 / 0.55 |
| truth's mean cosine to the other truths | -0.202 | 0.067 | 0.092 | **-0.004** | 0.623 | 0.379 | 0.536 | **0.510** | 0.002 | 0.002 | -0.04 / 0.58 |
| prediction's nearest-neighbour cosine to another prediction | 0.216 | 0.295 | 0.264 | **0.255** | 0.355 | 0.404 | 0.429 | **0.397** | 0.036 | 0.039 | 0.35 / 0.39 |
| prediction's mean cosine to the other predictions (hubness) | -0.115 | 0.006 | -0.008 | **-0.034** | 0.500 | 0.466 | 0.498 | **0.488** | 0.001 | 0.002 | 0.00 / 0.54 |
| median cosine of the prediction to the OTHER truths | -0.075 | -0.015 | -0.124 | **-0.072** | 0.569 | 0.546 | 0.507 | **0.540** | 0.000 | 0.000 | -0.18 / 0.61 |
| max cosine of the prediction to another truth | 0.190 | 0.065 | 0.080 | **0.116** | 0.382 | 0.568 | 0.424 | **0.460** | 0.032 | 0.034 | 0.17 / 0.54 |
| delivered ‖pred‖ / ‖truth‖ (scored cols) | 0.098 | 0.166 | 0.196 | **0.158** | 0.433 | 0.405 | 0.512 | **0.451** | 0.544 | 0.573 | 0.07 / 0.48 |
| render effect on pds: delivered − effect-level | 0.010 | 0.155 | 0.183 | **0.115** | 0.463 | 0.530 | 0.310 | **0.432** | -0.016 | 0.000 | 0.18 / 0.38 |
| own cosine (delivered prediction vs its truth) | 0.973 | 0.966 | 0.962 | **0.958** | 0.002 | 0.002 | 0.004 | **0.003** | -0.011 | 0.023 | 1.00 / 0.00 |
| own cosine minus the best competitor cosine | 0.949 | 0.953 | 0.942 | **0.939** | 0.050 | 0.030 | 0.054 | **0.044** | -0.043 | -0.011 | 0.93 / 0.03 |
| pds of the whitened base alone (no residual) | 0.927 | 0.912 | 0.903 | **0.900** | 0.028 | 0.025 | 0.047 | **0.034** | 0.189 | 0.950 | 0.88 / 0.05 |
| pds of base + deflated CD4 alone (no union) | 0.980 | 0.984 | 0.982 | **0.973** | 0.005 | 0.000 | 0.011 | **0.005** | 0.206 | 0.962 | 0.97 / 0.00 |
| pds of base + union alone (no CD4) | 0.949 | 0.940 | 0.914 | **0.924** | 0.021 | 0.028 | 0.040 | **0.030** | 0.205 | 0.953 | 0.93 / 0.05 |
| pds of h1union (the un-deflated stack) | 0.975 | 0.955 | 0.941 | **0.945** | 0.015 | 0.027 | 0.035 | **0.026** | 0.208 | 0.953 | 0.97 / 0.01 |
| pds at effect level (un-rendered) | 0.971 | 0.964 | 0.958 | **0.945** | 0.017 | 0.009 | 0.025 | **0.017** | 0.198 | 0.978 | 0.92 / 0.05 |

## T5. Deficit decomposition (ceiling: `pds` = 1 on every target, i.e. every prediction nearest its own truth; the fold deficit is 1 − mean `pds`)

### T5a. Where the beats are: prediction side vs truth side

| fold | beats (prediction, truth, frame) | share from the 10 most-beaten predictions | bottom-decile predictions | bottom quartile | share from the 10 most-beating (hub) truths | top-30 hub truths | truths beating nothing | top-10 hub truths |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `k562_gwps` | 39928 | 18.0 % | 42.3 % | 78.7 % | **4.5 %** | 13.0 % | 0 | `RBM18` 197, `BCL7C` 189, `TRMT1` 186, `ZKSCAN1` 181, `PHF23` 177, `LMAN1` 177 |
| `orion_hct116` | 49869 | 16.1 % | 40.4 % | 75.8 % | **4.3 %** | 12.2 % | 0 | `HMGN3` 227, `ZHX3` 226, `ZBTB34` 217, `SRSF5` 212, `DVL3` 211, `PDK3` 209 |
| `orion_hek293t` | 56601 | 14.1 % | 36.8 % | 70.8 % | **4.4 %** | 12.5 % | 0 | `ACLY` 278, `ZNF444` 260, `TMF1` 257, `TBC1D13` 247, `DEAF1` 245, `TWF2` 242 |
| `kolf2` | 80131 | 9.7 % | 25.6 % | 55.1 % | **4.2 %** | 12.3 % | 0 | `BRD1` 351, `PDCD4` 346, `COQ10B` 341, `IFI27L2` 341, `ZNF530` 338, `ADNP` 336 |

A hub truth would beat many predictions; a uniformly spread loss puts every truth at ≈ 1/n of the beats (0.3–0.4 %). The top-10 hub truths take 4–5 % against a uniform 3.4–3.7 %: there is no truth-side hub left. The loss sits on the prediction side.

### T5b. How much of the residual-hub story remains (h1union → cd4svd5)

| fold | mean pairwise prediction cosine (hubness) h1union → cd4svd5 | ρ(prediction hubness, `pds`) h1union → cd4svd5 | ρ(median cos to other truths, `pds`) h1union → cd4svd5 | top-10 hub-truth beat share h1union → cd4svd5 |
|---|---:|---:|---:|---:|
| `k562_gwps` | -0.0001 → -0.0003 | -0.089 → -0.115 | -0.119 → -0.075 | 4.7 % → 4.5 % |
| `orion_hct116` | 0.0026 → 0.0025 | -0.021 → +0.006 | -0.058 → -0.015 | 4.6 % → 4.3 % |
| `orion_hek293t` | 0.0023 → 0.0022 | -0.021 → -0.008 | -0.115 → -0.124 | 4.6 % → 4.4 % |
| `kolf2` | 0.0016 → 0.0014 | -0.003 → +0.005 | -0.234 → -0.181 | 4.2 % → 4.2 % |

### T5c. The partition of the deficit by mechanism (each target assigned to the FIRST group it matches, in this order)

Groups: **(b) no base support** — carried by none of the fold's whitened sources (the base row is the mean-response filler); **(e) non-transferable direction** — the whitened donor base's cosine to the truth ≤ 0 (the donors' consensus points away from this line's response); **(r) residual breaks it** — base·truth > 0 but `pds` fell below the base-only `pds` by more than 0.1 after the residual; **(a) hub / confusable** — the prediction's median cosine to the other truths is in the fold's top quintile; **(c) small true effect** — truth ‖Δ‖ in the fold's bottom quintile; **(d) render / magnitude** — the render moved `pds` by more than 0.05; **(o) other** — none of the above (the prediction points the right way but a few specific truths sit closer).

| fold | deficit | (b) no base | (e) base points away | (r) residual breaks it | (a) hub | (c) small truth | (d) render | (o) other: specific confusion | n targets per group |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `k562_gwps` | 0.1846 | 3.2 % | 57.5 % | 11.7 % | 3.0 % | 8.7 % | 8.7 % | 7.2 % | b 5, e 47, r 14, a 21, c 24, d 19, o 139, zero-deficit 82 |
| `orion_hct116` | 0.1983 | 0.0 % | 51.5 % | 16.1 % | 7.8 % | 7.2 % | 10.6 % | 6.8 % | b 0, e 46, r 21, a 31, c 27, d 28, o 137, zero-deficit 77 |
| `orion_hek293t` | 0.2190 | 0.1 % | 49.0 % | 18.1 % | 8.6 % | 2.8 % | 12.3 % | 9.2 % | b 1, e 52, r 25, a 28, c 20, d 34, o 134, zero-deficit 65 |
| `kolf2` | 0.3597 | 0.0 % | 64.0 % | 10.6 % | 5.6 % | 4.1 % | 11.6 % | 4.1 % | b 0, e 101, r 22, a 24, c 22, d 39, o 65, zero-deficit 34 |

Marginal (overlapping) shares of the deficit, i.e. each condition on its own:

| fold | (b) no base | (e) base·truth ≤ 0 | (r) residual breaks | (a) hub quintile | (c) small-truth quintile | (d) render > 0.05 |
|---|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 3.2 % (n 5) | 60.3 % (n 51) | 11.7 % (n 14) | 18.0 % (n 54) | 21.9 % (n 54) | 58.3 % (n 71) |
| `orion_hct116` | 0.0 % (n 0) | 51.5 % (n 46) | 16.1 % (n 21) | 17.4 % (n 58) | 19.8 % (n 58) | 59.6 % (n 87) |
| `orion_hek293t` | 0.1 % (n 1) | 49.0 % (n 52) | 18.1 % (n 25) | 26.8 % (n 59) | 18.7 % (n 59) | 57.6 % (n 97) |
| `kolf2` | 0.0 % (n 0) | 64.0 % (n 101) | 10.6 % (n 22) | 22.5 % (n 55) | 19.3 % (n 55) | 70.2 % (n 153) |

### T5d. The deficit by bins of base·truth cosine and of held-vs-donor truth agreement (cancer folds pooled)


| base·truth cosine bin | n targets | share of targets | mean `pds` | share of the pooled deficit | deficit per target |
|---|---:|---:|---:|---:|---:|
| (-1, 0] | 149 | 17.5 % | 0.389 | **53.1 %** | 0.611 |
| (0, 0.05] | 564 | 66.1 % | 0.858 | **46.8 %** | 0.142 |
| (0.05, 0.1] | 108 | 12.7 % | 0.998 | **0.1 %** | 0.002 |
| (0.1, 0.2] | 31 | 3.6 % | 1.000 | **0.0 %** | 0.000 |
| (0.2, 1.01] | 1 | 0.1 % | 1.000 | **0.0 %** | 0.000 |

| held-vs-donor truth agreement bin | n targets | share of targets | mean `pds` | share of the pooled deficit | deficit per target |
|---|---:|---:|---:|---:|---:|
| (-1, 0] | 178 | 20.9 % | 0.564 | **45.3 %** | 0.436 |
| (0, 0.05] | 394 | 46.2 % | 0.835 | **37.8 %** | 0.165 |
| (0.05, 0.1] | 144 | 16.9 % | 0.956 | **3.7 %** | 0.044 |
| (0.1, 0.2] | 59 | 6.9 % | 0.964 | **1.2 %** | 0.036 |
| (0.2, 1.01] | 2 | 0.2 % | 1.000 | **0.0 %** | 0.000 |

| donor-donor truth agreement (observable) bin | n targets | share of targets | mean `pds` | share of the pooled deficit | deficit per target |
|---|---:|---:|---:|---:|---:|
| (-1, 0] | 210 | 24.6 % | 0.707 | **35.9 %** | 0.293 |
| (0, 0.05] | 348 | 40.8 % | 0.808 | **38.9 %** | 0.192 |
| (0.05, 0.1] | 145 | 17.0 % | 0.900 | **8.5 %** | 0.100 |
| (0.1, 0.2] | 62 | 7.3 % | 0.903 | **3.5 %** | 0.097 |
| (0.2, 1.01] | 12 | 1.4 % | 0.831 | **1.2 %** | 0.169 |

## T6. The `h1union` → `cd4svd5` step, per target

| fold | mean Δ`pds` | gainers / losers / unchanged | sum of gains | sum of losses | top-3 share of the net gain | ρ(Δ, base·truth) | ρ(Δ, CD4 row energy removed) | ρ(Δ, CD4·truth before) | ρ(Δ, CD4 row norm) | ρ(Δ, `pds` before) | AUC losers-higher: E removed / ‖CD4‖ / base·truth / `pds` before |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `k562_gwps` | +0.00999 | 102 / 71 / 96 | +6.512 | -3.826 | 55 % | -0.078 | +0.019 | -0.364 | +0.029 | -0.197 | 0.50 / 0.44 / 0.30 / 0.34 |
| `orion_hct116` | +0.00945 | 112 / 88 / 90 | +8.582 | -5.843 | 55 % | +0.049 | +0.103 | -0.360 | +0.052 | -0.141 | 0.45 / 0.46 / 0.29 / 0.35 |
| `orion_hek293t` | +0.01736 | 134 / 83 / 77 | +10.096 | -4.993 | 33 % | -0.049 | +0.157 | -0.382 | +0.154 | -0.238 | 0.47 / 0.40 / 0.35 / 0.41 |
| `kolf2` | +0.00414 | 113 / 107 / 53 | +7.148 | -6.018 | 76 % | +0.101 | -0.016 | -0.352 | -0.016 | -0.089 | 0.50 / 0.50 / 0.34 / 0.41 |

Fold-to-fold consistency of the per-target step (Spearman of Δ`pds` across cancer folds, common targets): K562–HCT116 +0.033; K562–HEK293T -0.026; HCT116–HEK293T -0.101.

| biggest movers (top-5 losers and gainers of each cancer fold) | Δ K562 | Δ HCT116 | Δ HEK293T | Δ KOLF2 | same sign on all cancer folds? |
|---|---:|---:|---:|---:|---|
| `PLCG1` | +0.394 | -0.617 | +0.661 | -0.537 | no |
| `CDK8` | -0.389 | — | — | +0.200 | yes |
| `PGLS` | +0.498 | -0.531 | +0.093 | +0.047 | no |
| `PPIP5K2` | +0.260 | +0.722 | +0.041 | +0.017 | yes |
| `LMAN1` | +0.402 | -0.025 | +0.324 | -0.017 | no |
| `MYD88` | -0.157 | +0.255 | -0.265 | -0.036 | no |
| `ZNF530` | -0.088 | -0.072 | +0.512 | +0.163 | no |
| `MBD4` | -0.260 | -0.198 | +0.199 | -0.246 | no |
| `ZNF324` | +0.586 | +0.032 | -0.031 | +0.061 | no |
| `SMARCC2` | +0.040 | -0.375 | +0.207 | +0.067 | no |
| `FHL3` | +0.039 | -0.258 | -0.316 | -0.108 | no |
| `SFXN4` | +0.005 | -0.042 | -0.538 | +0.017 | no |
| `ZNF219` | +0.080 | +0.115 | -0.377 | +0.049 | no |
| `PSMB9` | — | -0.353 | +0.023 | — | no |
| `ZNF768` | +0.042 | +0.360 | +0.150 | -0.015 | yes |
| `ZNF670` | +0.109 | +0.416 | -0.006 | +0.028 | no |
| `SNX4` | -0.019 | +0.000 | +0.506 | -0.048 | no |
| `PI4KB` | -0.106 | -0.053 | +0.337 | — | no |
| `CHKB` | +0.001 | +0.000 | +0.466 | — | no |
| `MYRF` | +0.021 | -0.078 | -0.358 | +0.050 | no |
| `BCKDK` | -0.274 | -0.076 | +0.106 | -0.059 | no |
| `TMCO6` | +0.113 | +0.287 | -0.052 | -0.099 | no |
| `ZNF496` | +0.297 | -0.042 | +0.113 | +0.147 | no |
| `BCL7C` | -0.042 | -0.338 | -0.039 | -0.049 | yes |
| `LPCAT3` | -0.340 | -0.020 | -0.053 | +0.189 | yes |
| `ZNF16` | -0.221 | -0.060 | +0.028 | -0.002 | no |
| `MSANTD4` | +0.001 | +0.286 | +0.015 | -0.134 | yes |

| covariate (cancer folds pooled) | median, losers (Δ < 0) | median, gainers (Δ > 0) | median, unchanged | AUC losers-higher |
|---|---:|---:|---:|---:|
| cov|base_cos_truth | 0.008 | 0.012 | 0.048 | 0.315 |
| cov|cd4_energy_removed | 0.088 | 0.106 | 0.090 | 0.473 |
| cov|cd4_norm | 4.3 | 4.5 | 4.7 | 0.436 |
| cov|cd4_cos_truth | 0.018 | -0.001 | 0.056 | 0.521 |
| cov|cd4_deflated_cos_truth | 0.006 | 0.007 | 0.050 | 0.354 |
| cov|cd4_cos_base | 0.006 | 0.009 | 0.061 | 0.335 |
| h1union|pds | 0.817 | 0.793 | 1.000 | 0.372 |
| h1union|pred_mean | 0.002 | 0.002 | 0.001 | 0.539 |
| agree_held_vs_donors | 0.012 | 0.017 | 0.052 | 0.355 |
| agree_donor_donor | 0.014 | 0.018 | 0.046 | 0.398 |
| cov|union_count | 0.000 | 0.000 | 0.000 | 0.431 |
| cov|n_whitened | 2.000 | 2.000 | 2.000 | 0.487 |
| cd4svd5|truth_norm | 6.7 | 6.8 | 6.8 | 0.487 |
| cov|resid_share | 0.323 | 0.337 | 0.365 | 0.434 |

## T7. Component attribution of the loss (per fold): the same targets, four objects

| fold | mean `pds`: BASE / CD4ONLY / UNIONONLY / h1union / cd4svd5 | ρ(BASE `pds`, cd4svd5 `pds`) | targets with BASE `pds` < 0.5 that cd4svd5 lifts above 0.5 | targets with BASE `pds` ≥ 0.9 that cd4svd5 drops below 0.5 | share of cd4svd5's deficit on targets already < 0.5 in BASE | render: mean `pds` delivered − effect-level (BASE / cd4svd5) |
|---|---|---:|---:|---:|---:|---|
| `k562_gwps` | 0.7811 / 0.8067 / 0.7916 / 0.8054 / 0.8154 | +0.927 | 17 of 51 | 0 of 150 | 61.3 % | -0.0072 / -0.0009 |
| `orion_hct116` | 0.7821 / 0.8004 / 0.7866 / 0.7922 / 0.8017 | +0.912 | 12 of 48 | 0 of 155 | 51.9 % | -0.0128 / -0.0095 |
| `orion_hek293t` | 0.7738 / 0.7816 / 0.7714 / 0.7636 / 0.7810 | +0.903 | 17 of 55 | 2 of 146 | 52.8 % | -0.0035 / -0.0046 |
| `kolf2` | 0.6034 / 0.6204 / 0.6265 / 0.6362 / 0.6403 | +0.884 | 21 of 101 | 0 of 62 | 68.4 % | -0.0127 / -0.0137 |
