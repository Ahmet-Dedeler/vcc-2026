# 255 — `cd4svd10` BUILT: the CD4 residual's top-10 target-shared directions deflated before sizing (`VCC_CD4_DEFLATE_K` 5 → 10 on top of `cd4svd5`). Both known-answer rows pass to the digit; every gate of 249's set passes with `cd4svd5` as the reference. The measured increment over `cd4svd5` is **not a cleared arm** (`pds` +0.00054 on 1/3 PLOLO folds, priced −0.0005 … −0.0002 overall) — the package exists so that branch A of `E249/DECISION-TREE.md` can buy the board's k-response, and for nothing else. **BUILDABLE — verdict line at the bottom (§17).**

**Date 2026-09-04 (19:10–20:00 UTC).** Deterministic. **No `vcc submit`, no `cell-eval2`, no model name, no
description, no slot proposal.** No file edited outside this directory: the builder is git HEAD's
`scripts/make_submission_pm_cd4.py` (`3615d20f…`, 254's commit `32c5219`; `git status` on it EMPTY before and after),
and `VCC_CD4_DEFLATE_K` is 249's flag (`be029c8`), read at import and applied in one place. Packages written under
`submissions/255-*` only. Every ≥ 15 GiB stage ran alone under the mandate's gate (`memory_pressure` ≥ 40 % free AND
`vm_stat` free+inactive ≥ 30 GiB AND zero other builders AND zero `vcc prep` AND the working-tree builder hashed to
HEAD's bytes before each row); peak RSS per stage in §14. The two armed submit scripts (237's pid 74625, 249's pid
77532) slept through everything. **The incumbent was NOT rebuilt**: the flag-off byte identity through these exact
builder bytes is 254 §KA 1 (§4 below).

| row / falsifier | required | measured | verdict |
|---|---|---|---|
| **KA 1** the staged k = 10 arm reproduces 246's `B2_cd4_svd10_resolved` | macro −1.67 % / aggregate −1.54 % / Δ`pds` +0.01281, three cancer folds, to the digit; `energy_removed_frac` ≈ 0.286 | **−1.669499 % / −1.536907 % / +0.012807143**, max\|Δ\| **0.0e+00** vs 246's json rows AND vs `price246.json`'s digits (−0.01669498660858859 / +0.012807142917080442) on 24/24 known-answer rows over 4 folds; `pds` +0.00982 / +0.01503 / +0.01357 (3/3); energy removed **0.2830** (3-fold) / **0.2843** (4-fold; per fold 0.28626 / 0.28149 / 0.28115 / 0.28813); production's own panel 0.2791 over 291 rows; effect matrix vs 246's `compose` max\|Δ\| 0.0 on every fold | **PASS** |
| **KA 2** `VCC_CD4_DEFLATE_K=5` still stages to 249's `a28d5be2…` | `a28d5be2…` | (a) re-hash of 249's staged `effects-svd5.npy`: array sha **`a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11`**, file sha `979959e8…` = 249's record; (b) re-staged through HEAD's builder: **`a28d5be2…`**, fraction 0.34104860341941284; (c) k = 5 through the builder on every PLOLO fold == 246 `B2_cd4_svd5_resolved` and == 249's own `repro249-<fold>.json` rows at 0.0e+00 | **PASS** |
| flag-off byte identity (not rebuilt; cited) | `8df4dadd…` on 12,648,807,800 B through builder `3615d20f…` | **254 §KA 1**: `submissions/254-ka-h1union-g028.h5ad.meta.json` records `output.sha256` `8df4dadd…`, 12,648,807,800 B, builder `3615d20f…`, `run254.sh` 1331.63 s / 20.04 GiB, `def_postbuild` 32/32 (verified off that meta here, §4) | **CITED** |
| **F1** rule-10 audit diff vs `249-cd4svd5-f034.h5ad.meta.json` | exactly `cd4_deflate_k` 5 → 10 + its CD4-frame consequences | **config changed = {`cd4_deflate_k`} (5 → 10), added 0, removed 0; audit changed = the 15 pre-declared keys (249's 11 + the deflation block's `k`, `energy_removed_frac`, `top_k_energy_frac_of_centred`, `basis_u_k_sha256`), added 0, removed 0; `n_carried` 291 pinned; 117/117 pinned scalars (58/58 `effects.h1_union.*`), frozen inputs 0, implementation = the builder file only (`2d63b979…` → `3615d20f…`, the same single difference 254 measured against the same reference), unclassified 0** | **PASS** |
| **F2** profile split by CD4-carried row, scale-corrected | carried rows move in direction; non-carried rows one shared scalar | effect level: 276 carried rows cos ≤ 0.99982 (min 0.8858), 24 non-carried = `cd4svd5` × **0.9966440037807254** to **9.5e−17** (separation 2.0e+14×); delivered: 828 carried blocks cos median 0.99707 (min 0.8739), 72 non-carried cos ≥ 0.99989, norm ratio 0.99630–0.99705, shared scale 0.99666594; scale-corrected displacement **0.02482** (carried min) vs **0.01480** (non-carried max), **CLEAN** 1.7× | **PASS** |
| **F3** non-vacuity | zero-block count = 0 | effect level ‖ΔE‖/‖E‖ **0.096311**, 0/300 rows at zero (min 0.003356); delivered ‖ΔE‖/‖E‖ **0.103137**, per block min 0.007729 / median 0.07194 / max 0.4863, **0/900 at zero**, 0/900 < 1e−6, 900/900 digests differ | **PASS** |
| **F4** challenge `n_pred` within ±2 % of `cd4svd5` | ±2 % | median **2,599 / 2,599 / 2,599** for both packages (40 perturbations × 3 contexts), pooled **+0.00 %** | **PASS** |
| **F5** delivered amplitude off the package | `A_c`/`B_c` within 0.5 % of `cd4svd5`'s (243's stream) | `√(A_arm/A_cd4svd5)` − 1 = **−0.4558 %** (contexts ×0.99081 / 0.99095 / 0.99096), `B` −0.4577 %; **vs `h1union` −0.8626 %** (beyond 249's own bound when read against the ORIGINAL incumbent, §9); predicted board raw `mse` **0.82230** at k_cos 1 (0.82210–0.82348 across 243's classes) = +0.0032 raw over `cd4svd5`'s predicted 0.81908 | **PASS (at 91 % of the bound)** |
| **F6** `vcc prep --dry-run` clean; nnz; max counts/cell | clean; vs 4.75e9; vs 1e6 | ✓ clean (360,000 cells, A/B/C = 120,000 each, targets verified, counts-preserved); nnz **1,577,446,557 = 33.21 %** of cap (`cd4svd5` −0.011 %); max **22,323** per cell (44.8× under); dry-run 6.5 s / 14.92 GiB, prep 15.25 s / 14.96 GiB | **PASS** |

---

## 1. What was built, and what was not changed

**One variable.** `VCC_CD4_DEFLATE_K=10` in place of `cd4svd5`'s 5. `deflate_cd4_residual` (249, 246's
`lib246.deflate(res, "svd", k)` verbatim) centres the CD4 residual over its 291 carried rows and removes the top-10
eigen-directions of `R_c R_cᵀ` before `solve_cd4_dose` sizes it. Nothing in the builder was edited; no flag was added;
every other inert flag (`VCC_CD4_SIZING`/`VCC_CD4_TARGET_G` 247, `VCC_UNION2_*` 248, `VCC_CD4_DEFLATE_BASIS` 251,
`VCC_CD4_PARTIAL_G` 254) is explicitly UNSET by `build255.sh` and `stage255.sh` (`env -u`), and `repro255.py` asserts
at import that none of them is live. Shipped configuration otherwise: `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34
VCC_H1_UNION_G=0.28`, `VCC_UNION_SOURCES` unset (the H1 pair), `--min-cpm 10 --tail-spread 1.0 --keep-top 8000
--keep-spread 0.25 --overdispersion 4.0 --whitener-id union_pruned --emitter rc_callset`, operator `e2f65ecf…`.

**What k = 10 does to the CD4 frame (production, A/B/C; `work/stage/stage-svd10.json`):**

| scalar | `cd4svd5` (k = 5) | **`cd4svd10` (k = 10)** | ratio |
|---|---:|---:|---:|
| `deflation.energy_removed_frac` | 0.19866046506100588 | **0.2790712359118336** | directions 6–10 carry 0.0804 |
| `deflation.top_k_energy_frac_of_centred` | 0.1982359547102955 | **0.27868932324044193** | |
| `deflation.basis_u_k_sha256` | `f06e0262…` (`U_5`) | **`9bbf3d41…`** (`U_10`) | |
| `deflation.n_carried` | 291 | 291 | pinned |
| `residual_median_row_norm_at_unit_fraction` | 4.2162920214222215 | **4.031491147959811** | ×0.95617 |
| `f_per_unit_fraction` | 0.9969253548939966 | **0.953229928812221** | ×0.95617 |
| `solved_fraction` = `effects.cd4.fraction` | 0.34104860341941284 | **0.3566820446181956** | **×1.04584** (PLOLO's frame ×1.0442) |
| merged preserve scalar (one number, written twice) | 0.9115882187626339 | **0.9085289321469312** | ×0.996644 |
| union `solved_fraction` | 0.19207161964444322 | 0.19207161964444322 | untouched |
| `effects.final_sha256` | `a28d5be2…` | **`a2ab1df7f88db52140d97c1f8b717d4a78326ad324f5192a52d0299069f4f9f6`** | |

The deflation audit was re-derived offline from the saved pre-deflation residual through the builder's function at
k = 10 AND at k = 5 and equals both staged records (`effdiff255.py`).

---

## 2. Known-answer row 1 — the staged k = 10 arm reproduces 246 (`repro255.py`, one fold per process)

162's PLOLO with the **production builder's code substituted in** (249's `repro249.py` pattern): the builder's
`deflate_cd4_residual` for `lib246.deflate`, and the whole of `apply_frozen_cd4_residual` (deflate → size → merge →
one preserve) for 246's `compose`, the builder imported with the arm's env (`VCC_CD4_DEFLATE_K=10`). Fold caches are
byte copies of 246's/249's (`work/plolo/`, shas verified). Tolerance 1e−9 on `cos_agg`, `a23_cos`, `pds`.

| fold | KA-A `STACK` vs 246 | KA-B builder deflation k=10 vs `lib246.deflate` | KA-C arm through `apply_frozen_cd4_residual` vs 246 `B2_cd4_svd10_resolved` | solved fraction (stack → k=5 → k=10) | union moved | arm k=10: macro / agg / Δ`pds` |
|---|---|---|---|---|---|---|
| `k562_gwps` | 0.0e+00 | bit-identical, audit identical (264 rows, E rm 0.286264777) | **0.0e+00**; effect matrix vs 246's compose 0.0 | 0.373761 → 0.395879 → **0.413735** (= 246's to 0.0e+00) | 0.0e+00 | −2.83 % / −2.89 % / +0.009821 |
| `orion_hct116` | 0.0e+00 | bit-identical (282 rows, 0.281490551) | **0.0e+00**; 0.0 | 0.370782 → 0.394269 → **0.409793** | 0.0e+00 | −2.04 % / −1.88 % / +0.015030 |
| `orion_hek293t` | 0.0e+00 | bit-identical (286 rows, 0.281145900) | **0.0e+00**; 0.0 | 0.376056 → 0.397869 → **0.416984** | 0.0e+00 | −0.14 % / +0.16 % / +0.013571 |
| `kolf2` (apart) | 0.0e+00 | bit-identical (266 rows, 0.288127475) | **0.0e+00**; 0.0 | 0.316408 → 0.336448 → **0.350322** | 0.0e+00 | −2.82 % / −4.91 % / +0.004893 |

Three cancer folds, against 246's own json rows recomputed the same way and against `price246.json`:

| | this run (builder) | 246 json rows | `price246.json` | \|Δ\| | mandate |
|---|---:|---:|---:|---:|---:|
| rel MACRO cosine | **−0.016694987** | −0.016694987 | −0.016694987 | 0.0e+00 | −1.67 % |
| rel AGGREGATE cosine | **−0.015369070** | −0.015369070 | — | 0.0e+00 | −1.54 % |
| Δ`pds` | **+0.012807143** | +0.012807143 | +0.012807143 | 0.0e+00 | +0.01281 |

`pds` positive on 3/3 cancer folds; `energy_removed_frac` **0.2830** over the three cancer folds, **0.2843** over four
(the mandate's "≈ 0.286" is 246 §3.3's `k562_gwps`-leading table entry; the per-fold values are 0.28626 / 0.28149 /
0.28115 / 0.28813). CD4 solved fraction in PLOLO's frame 0.3735 → 0.4135 (×1.107 over the stack, ×1.0442 over k = 5);
production's frame 0.3221 → 0.3567 (×1.107, ×1.0458). 24/24 known-answer rows at 0.0e+00; peak RSS 4.97–5.21 GiB per
fold, 8.3–8.4 s each.

**Energy-matched random control through the same builder path** (246's `B2_cd4_randE10`, seed `20260904 + 10`,
bisection to svd10's energy): dims 73 / 79 / 78 / 72 (= 246's), the subspace bit-identical to
`lib246.deflate(…, "rand", k′, seed)`; fixed sizing macro −3.837 % / agg −4.622 % / Δ`pds` **−0.003214** (= 246's row to
0.0e+00); re-solved through the builder −3.042 % / −3.401 % / **−0.005256**. Removing the same ~29 % of the energy in
random directions buys nothing and costs 2–3× the cosine: at k = 10 as at k = 5, **the builder applies a projector,
not a norm cut**.

### 2.1 The increment k 5 → 10 — the slot's one variable, measured (KA-C′)

The k = 5 arm was pushed through the SAME builder in the same process (the module global switched to 5) and equals
246's `B2_cd4_svd5_resolved` AND 249's own per-fold rows at 0.0e+00 on every fold, so the increment is read off two
matrices produced by identical code:

| fold | macro (k10 vs k5) | aggregate | Δ`pds` (k10 − k5) | Δq_w |
|---|---:|---:|---:|---:|
| `k562_gwps` | −1.064 % | −1.285 % | **−0.000166** | −0.00024 |
| `orion_hct116` | −0.129 % | +0.074 % | **+0.005584** | −0.00012 |
| `orion_hek293t` | −1.393 % | −1.990 % | **−0.003784** | −0.00149 |
| **mean, 3 cancer folds** | **−0.862 %** | **−1.067 %** | **+0.000544 (1/3 folds, 0.36 PLOLO-σ)** | −0.00062 |
| `kolf2` apart | −0.931 % | −1.488 % | +0.000754 | +0.00041 |

Concentration (246 §3.4's leave-top-3-out on per-target Δ`pds`, mean over cancer folds): the k = 10 arm vs `h1union`
retains **52.5 %** (k = 5: 52.4 %); the increment's retention is undefined in any useful sense (+350 % — its total is
near zero, with 85–103 targets up and 72–99 down per fold). **The k ladder is flat past k = 5 on the A/B/C PLOLO folds:
directions 6–10 remove another 8 % of the residual's energy and buy no `pds` beyond noise, at −1.07 % aggregate
cosine.** That is 246 §3.3's "k = 3–5 is the knee" restated through the builder; 250's k 10–15 plateau is a property
of its disjoint panels, not of A/B/C.

**Price** (`price255.py`; 246 §4's rate card at k_cos 1, `pds` at 0.7 / 1.0 / 2.0×; the two 246 rows reproduce at
0.0e+00):

| arm | macro | Δ`pds` | `mse` | `reach` | `nmae` | `fid` | `jac` | `pds` 0.7 / 1.0 / 2.0× | **TOTAL 0.7 / 1.0 / 2.0×** | vs +0.0019 gate |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| `cd4svd10` vs `h1union` (KA = `B2_cd4_svd10_resolved`) | −1.669 % | +0.012807 (8.6σ) | −0.00055 | −0.00037 | −0.00015 | −0.00013 | −0.00002 | +0.00334 / +0.00477 / +0.00955 | **+0.00214 / +0.00357 / +0.00834** | 1.1× / 1.9× / 4.4× |
| `cd4svd5` vs `h1union` (KA = `B2_cd4_svd5_resolved`) | −0.809 % | +0.012263 (8.2σ) | −0.00027 | −0.00018 | −0.00007 | −0.00006 | −0.00001 | +0.00320 / +0.00457 / +0.00914 | +0.00262 / +0.00399 / +0.00856 | 1.4× / 2.1× / 4.5× |
| **INCREMENT `cd4svd10` vs `cd4svd5`** | **−0.862 %** | **+0.000544 (0.4σ)** | −0.00028 | −0.00019 | −0.00008 | −0.00007 | −0.00001 | +0.00014 / +0.00020 / +0.00041 | **−0.00048 / −0.00042 / −0.00022** | **−0.25× / −0.22× / −0.11×** |

The exact `nmae` identity (six246) exists only for the fixed-sizing `B2_cd4_svd10` row: +0.00006, `e` 0.405 (a
magnitude bound). With the package's own `mse` leg (§9: −0.00057 instead of the card's −0.00028) the increment prices
−0.00077 / −0.00071 / −0.00051. **Not a cleared arm on any reading; what a slot buys is stated in the draft
preregistration (§12).**

---

## 3. Known-answer row 2 — the k = 5 path is untouched (`work/ka2-rehash249.json`, `stage255.sh`)

Nothing in the builder changed since 254 (git log `32c5219..HEAD -- src scripts` is empty; working tree ==
HEAD `3615d20f…`), so the row is cheap and was taken three ways: (a) 249's staged `effects-svd5.npy` re-hashed —
array sha `a28d5be2…`, file sha `979959e8…`, both equal to 249's `stage-svd5.json`; (b) `VCC_CD4_DEFLATE_K=5`
re-staged through HEAD's builder (`stage255.sh` row `svd5`, alone, gated) — **`a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11`**,
fraction 0.34104860341941284, preserve 0.9115882187626339, deflation block identical to 249's (`U_5` `f06e0262…`),
10 s / 20.11 GiB; (c) §2's KA-C′ on every PLOLO fold. The `effects-svd5.npy` this experiment staged is the reference
for §6a.

---

## 4. The flag-off byte identity — cited, not rebuilt

The mandate forbids rebuilding the incumbent and says to cite 254 §KA 1. Verified off the meta it left behind:
`submissions/254-ka-h1union-g028.h5ad.meta.json` records `output.sha256`
**`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6`**, 12,648,807,800 B, built by
`scripts/make_submission_pm_cd4.py` **`3615d20f…`** — byte-identical to the builder this experiment used for every
row — with neither `cd4_deflate_k` nor `cd4_partial_g` in its config (254 RESULTS.md line 16: `run254.sh`, both
flags unset, 1331.63 s / 20.04 GiB, `def_postbuild` 32/32, `moved=[]`). 254 also re-proved the staging identity
`abeb0749…` through the same bytes (its rows `off` and `pre`). Nothing here re-derives the incumbent; §6a reads the
k = 10 matrix against 249's staged `effects-off.npy` only after verifying that file hashes to `abeb0749…`.

---

## 5. Falsifier 1 — the rule-10 audit diff (`metadiff255.py`, sets pre-declared in the source before any diff was read)

Reference: **`submissions/249-cd4svd5-f034.h5ad.meta.json`** (the mandate's incumbent). Every difference, verbatim
(`work/metadiff255-vs249.log`; the `output.*` / `peak_rss_gib` bucket is reported, not asserted):

```
[A-config            ] CHANGED config.cd4_deflate_k                                                 5 -> 10
[B-effects/operator  ] CHANGED effects.cd4.deflation.basis_u_k_sha256                              f06e0262... -> 9bbf3d41...
[B-effects/operator  ] CHANGED effects.cd4.deflation.energy_removed_frac                           0.19866046506100588 -> 0.2790712359118336
[B-effects/operator  ] CHANGED effects.cd4.deflation.k                                             5 -> 10
[B-effects/operator  ] CHANGED effects.cd4.deflation.top_k_energy_frac_of_centred                  0.1982359547102955 -> 0.27868932324044193
[B-effects/operator  ] CHANGED effects.cd4.fraction                                                0.34104860341941284 -> 0.3566820446181956
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm    4.5227571142920135 -> 4.537986580097941
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.fraction                              0.34104860341941284 -> 0.3566820446181956
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.post_add_global_preserve_scale        0.9115882187626339 -> 0.9085289321469312
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.effective_f_at_cd4_fraction             0.04984626774469983 -> 0.04766149644061105
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.f_per_unit_fraction                     0.9969253548939966 -> 0.953229928812221
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction 4.2162920214222215 -> 4.031491147959811
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.solved_fraction                         0.34104860341941284 -> 0.3566820446181956
[B-effects/operator  ] CHANGED effects.final_sha256                                               a28d5be2... -> a2ab1df7...
[B-effects/operator  ] CHANGED effects.h1_union.merged_preserve.candidate_pre_preserve_median_norm 4.5227571142920135 -> 4.537986580097941
[B-effects/operator  ] CHANGED effects.h1_union.merged_preserve.post_add_global_preserve_scale     0.9115882187626339 -> 0.9085289321469312
[D-implementation    ] CHANGED hashes.implementation.scripts/make_submission_pm_cd4.py             2d63b979... -> 3615d20f...

A. CONFIG changed ['cd4_deflate_k'] expected ['cd4_deflate_k'];  added []  removed []
B. AUDIT changed == the 15 pre-declared;  added [];  removed []
   direction: residual_median_row_norm x0.956170 (down), f_per_unit x0.956170 (down), effective_f x0.956170 (down),
              solved_fraction / fraction / post_add_preserve.fraction x1.045839 (up), deflation.k x2 (up),
              energy_removed_frac x1.404766 (up), top_k_energy_frac x1.405838 (up); the merged preserve scalar is ONE number written twice
C. FROZEN INPUTS differences: 0    D. IMPLEMENTATION: the builder file only (expected)    F. UNCLASSIFIED: []
PINNED scalars (117) -- each byte-identical: 117/117, of which 58/58 `effects.h1_union.*` fields
FALSIFIER 1 PASS
```

Pinned by name and byte-identical: `config.alpha` 0.982 · `config.cd4_target_f` 0.34 · `config.h1_union_g` 0.28 ·
`config.h1_union_sources` `['gse295214_h1', 'h1_2025']` · `operator.sha256` `e2f65ecf…` · `operator.axis_sha256`
`edf6ba88…` · `effects.{baseline_pre,candidate_pre,post_whiten}_sha256` · all seven `effects.whitening.*` ·
`effects.cd4.{mode, fully_supported_targets 276, supported_entries 2,712,528, panel_residual_nonzero_entries 0, the
count histogram, deflation.n_carried 291}` · `effects.cd4.residual_frame.{target_f_requested 0.34,
whitened_base_median_row_norm, cd4_fraction_constant, fraction_source}` · **every `effects.h1_union.*` field** except
the two merged-preserve scalars · `frozen_bundle.sha256` · `output.n_rows` 360,000 · the emitter's `rc_callset` triple
and every other `config.*` key. `def_postbuild.py` on the arm (`--allow-config-diff cd4_deflate_k`, reference the
shipped 222 meta): **32/32**, `moved=['cd4_deflate_k']`, arms asserted by value (α, `f`, `g`, sources, emitter), operator
file/axis/metadata pins, CD4 complete targets 276, 360,000 rows.

---

## 6. Falsifier 2 — the profile moved on the CD4-carried rows and nowhere else

### 6a. On the effect matrix — exact (`effdiff255.py`, `effects-svd10.npy` vs this experiment's `effects-svd5.npy`, same builder bytes)

| row class | n | rel. displacement `‖ΔE_i‖/‖E_i‖` min / median / max | cosine(k5, k10) min / max | `‖E_k10 − s·E_k5‖/‖E_k5‖` max |
|---|---:|---|---|---:|
| **CD4-carried** (support > 0) | 276 | 0.019115 / 0.070080 / 0.464572 | 0.8857785 / 0.9998195 | 4.63e−01 (median 7.00e−02) |
| **NOT carried** | 24 | **exactly 0.003356 for all 24** (= \|1 − s\|) | **1.0000000000** | **9.45e−17** |

`s` = 0.9085289321469312 / 0.9115882187626339 = **0.9966440037807254**; the 24 non-carried rows' norm ratio is exactly
`s` (min = max to the last digit). Separation on the scale-corrected displacement **2.0e+14×**. The carried rows hold
99.989 % of ‖ΔE‖²; the five most displaced targets: DCP2 0.465, PGLS 0.430, PLCG1 0.336, PCGF1 0.281, SUPT7L 0.273.
For the record, against `h1union` (249's `effects-off.npy`, verified `abeb0749…`): ‖ΔE‖/‖E‖ 0.169511 (k = 5 was
0.142049), preserve ratio 0.99909, carried rows cos min 0.6606, non-carried scale-corrected max 9.1e−17.

### 6b. On the delivered counts (`profile255.py` = 249's `profile249.py`, carrying set = CD4 support; arm vs `249-cd4svd5-f034.h5ad`)

900 (context, target) blocks, 828 carried (276 targets × 3), 72 not:

| block class | n | digest identical | cosine(arm, inc) min / median / max | norm ratio min / max | `‖ΔE − s·E‖/‖E‖` min / median / max |
|---|---:|---:|---|---|---|
| **carried** | 828 | 0 | 0.87388115 / 0.99706633 / — | 0.8470 / 1.0083 | **0.024818** / 0.076667 / 0.484652 |
| **NOT carried** | 72 | 0 | 0.99988977 / 0.99994485 / — | 0.99629670 / 0.99705431 | 0.006960 / 0.010467 / **0.014798** |

Shared delivered scale read off the non-carried blocks **0.99666594** (effect level 0.99664400). Separation **CLEAN**:
the smallest carried displacement is **1.7×** the largest non-carried one (249: 1.9× under a 0.25 % scalar; here the
scalar is 0.33 % and the emitter's integer re-draw is the same 0.7–1.5 % on those blocks). The exact statement is
§6a's, at 1e−16. **Digest counts: 0/900 identical (0/828 carried, 0/72 non-carried).** Nonzero-structure digests
0/900 identical; panel-total sums 9/900 identical. **PASS.**

---

## 7. Falsifier 3 — non-vacuity. A zero would be the bug; there is none

| level | `‖ΔE‖/‖E‖` | per-block min / p01 / median / max | blocks exactly 0 | blocks < 1e−6 |
|---|---:|---|---:|---:|
| effect matrix (300 targets), vs `cd4svd5` | **0.096311** | 0.003356 / — / — / 0.465 | **0 / 300** | 0 / 300 |
| delivered profile (900 blocks), vs `cd4svd5` | **0.103137** | 0.007729 / 0.009285 / 0.071940 / 0.486349 | **0 / 900** | 0 / 900 |
| delivered, inside the 10-CPM gate | 0.103141 | 0.007710 / — / 0.071675 / 0.486555 | — | — |

Per context: A 0.104733, B 0.102020, C 0.102584. The extra five directions applied on every block in all three
contexts. 210's 1 % per-block floor (the script's own "FALSIFIER 4" print, which reads FAIL) holds on all 828 carried
blocks (min 2.5 % scale-corrected) and is by construction inapplicable to the 72 non-carried blocks (0.33 % preserve
scalar + the emitter's integer re-draw), exactly as 249 §7 recorded; the mandate's criterion (zero-block count = 0)
is met. **PASS.**

For the record: median delivered `|row|` ratio 0.99718, mean 0.99595, `A = mean|row|²` **0.99095**; `rc_callset` head
exposure: top-400 |delivered effect| overlap with `cd4svd5` min 210 / median 375 / max 400, implied call-set jaccard
lower bound min 0.864 / median 0.981; |Δ log2FC| over the DE universe max 0.837, mean 0.0042, signed block-mean
median +0.000005.

---

## 8. Falsifier 4 — challenge `n_pred` (`probe255.py` = 249's `probe249.py`; Wilcoxon + BH against the real `data/challenge/context_{A,B,C}.h5ad` controls, 40 perturbations per context)

| package | A | B | C | pooled median | ρ(−log p, \|H\|) in the call set, A / B / C (a15's key for both, cross-key by design) | call set UP |
|---|---:|---:|---:|---:|---|---|
| `cd4svd5` | 2,599 [2,599, 2,635] | 2,599 [2,599, 2,638] | 2,599 [2,599, 2,660] | 2,599 | +0.9309 / +0.9270 / +0.9239 (mean +0.9273) | 49.7 / 49.8 / 48.0 % |
| **`cd4svd10`** | **2,599** [2,599, 2,637] | **2,599** [2,599, 2,637] | **2,599** [2,599, 2,660] | **2,599** | +0.9273 / +0.9255 / +0.9185 (mean +0.9238) | 49.7 / 49.8 / 47.9 % |

**+0.00 %** against ±2 %. `rc_callset` pins `n_pred` by construction and it held; the `cd4svd5` row reproduces
249 §8's own numbers. **PASS.** (Peak RSS 13.31 GiB, 39.7 s, run alone after the prep.)

---

## 9. Falsifier 5 — delivered amplitude off the package, and the predicted board raw `mse` (`amp255.py`)

Per-block panel totals streamed off all three packages (214's `block_totals` via 243's `totals243.py`), 211's
per-context terms via 243's `terms`. **Two known-answer rows:** `h1union`'s re-streamed (`A_c`, `B_c`, `cap_c`) equal
243's `invert243.json` to **0.0e+00** on all three contexts and 197's identity at those terms reproduces its board raw
`mse` **0.8172630** to 0.0e+00; `cd4svd5`'s re-streamed (`A_c`, `B_c`, `cap_c`, `Kc_c`) equal 249's `amp249.json` to
**0.0e+00**.

| context | `A_c` arm | `A_c` cd4svd5 | ratio | `A_c` h1union | `B_c` arm / cd4svd5 | cap | leg |
|---|---:|---:|---:|---:|---|---:|---|
| A | 0.229505 | 0.231633 | ×0.99081 | 0.233629 | 0.227579 / 0.229698 | 0.233754 | B |
| B | 0.218958 | 0.220958 | ×0.99095 | 0.222708 | 0.217023 / 0.219014 | 0.239678 | B |
| C | 0.228927 | 0.231016 | ×0.99096 | 0.232897 | 0.226948 / 0.229028 | 0.266129 | B |

**Delivered amplitude `mean_c √(A_arm/A_cd4svd5)` − 1 = −0.4558 %** (`B` term −0.4577 %); bound ±0.5 % vs the
incumbent the mandate names: **PASS, at 91 % of the bound.** All three contexts stay on the `B` leg. **Against
`h1union` the drift is −0.8626 %** (`cd4svd5` alone was −0.4087 %): each five directions removed take ~0.45 % off the
delivered amplitude, so read against the ORIGINAL incumbent the k = 10 package is beyond 249's own 0.5 % shape-only
bound. The draft preregistration says this in so many words (§12): the arm is shape-only relative to `cd4svd5` only.

Predicted board raw `mse` from 197's identity at each package's own terms (243 §5.2), with
`cos = cos_h1union × (1 + k_cos × rel_agg)`, `rel_agg` −1.5369 % (k = 10) / −0.4617 % (k = 5) from §2, `cos_h1union`
0.107332:

| k_cos | arm cos | **predicted raw `mse` k = 10** | predicted k = 5 (249's table) | k10 − k5 | vs `h1union`'s 0.8172630 | `mse` leg vs `cd4svd5`, overall |
|---:|---:|---:|---:|---:|---:|---:|
| 0 (amplitude/budget leg alone) | 0.107332 | 0.81903 | 0.81809 | +0.00094 | +0.00176 | −0.00016 |
| 0.94 (new-source class) | 0.105782 | **0.82210** | 0.81902 | +0.00309 | +0.00484 | −0.00054 |
| 1.00 | 0.105683 | **0.82230** | 0.81908 | +0.00322 | +0.00504 | −0.00057 |
| 1.36 (re-dose of occupied rows — this class by 243 §4.4) | 0.105089 | **0.82348** | 0.81943 | +0.00405 | +0.00622 | −0.00071 |

The identity charges the k 5 → 10 step **+0.0032 raw `mse` at k_cos 1** (+0.0009 of it amplitude, the rest the extra
−1.07 % aggregate cosine), an `mse` leg of −0.00057 overall — twice what the rate card charges on the macro cosine.
Draft falsifier 2 fires above +0.0065 raw over `cd4svd5`'s board value.

---

## 10. Falsifier 6 — `vcc prep --dry-run`, nnz, max counts per cell

```
✓ dry run — no .vcc written
  input:  submissions/255-cd4svd10-f034.h5ad   cells: 360000   genes: 18533 (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'   contexts: column 'context' → A=120000, B=120000, C=120000
  targets: verified against the official list   normalization: counts-preserved
```

| cap | limit | measured | margin |
|---|---:|---:|---:|
| stored entries in `.X` | 4,750,000,000 | **1,577,446,557** (`cd4svd5` 1,577,615,064, −0.011 %; `h1union` 1,577,411,397, +0.002 %) | **33.21 %** of cap |
| counts per cell | 1,000,000 | **22,323** (min 17,695; identical to both incumbents', the exact-total emitter) | 44.8× under |
| rows | 360,000 | 360,000 | exact |
| cells per (context, target) block | 400 | 400 on all 900 | exact |

Dry run 6.5 s / 14.92 GiB; the real prep 15.25 s / 14.96 GiB. **PASS.**

---

## 11. The package

| artifact | bytes | sha256 |
|---|---:|---|
| **`submissions/255-cd4svd10-f034.vcc`** | **2,824,243,200** | **`f860632a9858b78065515e5d8cf9d2310b9f4dc6f54bb543b794b3c0a4612799`** |
| **`submissions/255-cd4svd10-f034.h5ad`** | **12,648,807,800** | **`3931e7277de91860a1d63c703a6f186c03d57870a9ba7a3ab2d422275ffab104`** |
| `submissions/255-cd4svd10-f034.h5ad.meta.json` | 25,623 | (`output.sha256` = the h5ad's; `effects.final_sha256` `a2ab1df7…` = the staged `svd10` matrix; builder `3615d20f…` = HEAD, commit `32c5219`) |

`work/build/SHA256SUMS.txt` carries both. Configuration (one variable off `cd4svd5`): `VCC_ALPHA=0.982
VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=10`, every other flag unset, 222's build flags, operator
`e2f65ecf…`. **No model name chosen, no description written, nothing submitted.** §14 applies at submit time.

---

## 12. The preregistration draft (`PREREGISTRATION-DRAFT.md`, for the lead — NOT committed by this run)

Six members vs the `cd4svd5` board point (246 §4's rate card on the increment; `pds` at three transfer factors):
`pds` +0.00014 / +0.00020 / +0.00041 (raw +0.000544 × the realised factor) · `mse` −0.00028 on the card, **−0.00057
from the package stream** (raw 0.82230 vs `cd4svd5`'s 0.81908 at k_cos 1) · `reach` −0.00019 · `nmae` −0.00008
(exact identity +0.00006, `e` 0.405, fixed-sizing row) · `fid` −0.00007 · `jac` −0.00001 · **overall −0.00048 /
−0.00042 / −0.00022 (card), −0.00077 / −0.00071 / −0.00051 (package `mse` leg)**. Five falsifiers in 249's form, read
against the `cd4svd5` row: (1) board raw `pds` below `cd4svd5`'s by more than 0.0040 (the other tail pre-committed
too: above by more than +0.0040 says the knee is beyond k = 5 and a k = 10 pinned basis must be built); (2) raw `mse`
worse than `cd4svd5`'s by more than +0.0065; (3) overall Δ vs `cd4svd5` outside [−0.0030, +0.0020]; (4) the rule-10
diff not exactly `cd4_deflate_k` 5 → 10 + its consequences (measured here: PASS); (5) delivered amplitude vs `cd4svd5`
beyond 0.5 % (measured here: −0.456 %, and −0.863 % vs `h1union`). **The draft says in its own §2 that the increment is
NOT a cleared arm** (negative on the card at every factor, `pds` on 1/3 folds, 0.36σ); what the slot buys — the
board's k-response, which decides October's pinned basis (`V5_abc` or a k = 10 basis) — is stated; the decision is
the lead's, on branch A only, and only if a slot priced at ≈ −0.0005 for information is acceptable.

---

## 13. Process notes, on the record

- Sequence: reads (READ-FIRST v2, AGENTS.md §11.1/§14, 249's scripts and results, 246 §3.3/§4, 243 §5.2, 254 §KA 1);
  KA row 2 re-hash; `repro255.py` × 4 folds; `stage255.sh` (`svd10`, then `svd5`), chained without a pause into
  `run255.sh` (arm build → postbuild → prep dry-run → prep) only after the `svd5` row matched `a28d5be2…`; then
  `effdiff255`, `metadiff255`, `profile255`, `amp255`, `probe255`, `price255`. 19:10–20:00 UTC.
- After the arm build the prep gate held at free+inactive 27.6–28.2 GiB for seven minutes (`memory_pressure` 82 % free)
  — the fresh 12.6 GB h5ad's clean page cache, 247 §3.4's / 251 §5's / 254 §13's symptom. As 254 did, the clean cache
  was dropped with `purge` (→ 39.1 GiB) and the gate opened on the next poll. No process was touched; the two armed
  submit scripts slept through everything (pids 74625 and 77532 still alive at the end).
- `other_builders()` read 0 and `vcc prep` count 0 at every gate; no other agent touched the builder (working tree ==
  HEAD before, during and after; `git status` shows only this experiment's files and other agents' `claude-run.jsonl`).
- `stage255.py` is 249's `stage249.py` byte-for-byte but for the header and module name; the current builder's
  `build_production_effects` takes its 254 additions as optional kwargs, so the call is unchanged (and the `svd5` row's
  `a28d5be2…` is the proof).

---

## 14. Cost and memory (every ≥ 15 GiB stage alone; the gate logged before each)

| stage | wall | peak RSS |
|---|---:|---:|
| KA 2 re-hash (`work/ka2-rehash249.json`) | 1 s | < 0.2 GiB |
| `repro255.py` × 4 folds | 8.1–8.6 s each | 4.97 / 5.19 / 5.21 / 5.12 GiB |
| staging `svd10` / `svd5` | 13 s / 10 s | 17.42 / 20.11 GiB |
| **ARM build** (`build255.sh`, `VCC_CD4_DEFLATE_K=10`) | **1317.78 s** | **20.20 GiB** |
| `def_postbuild.py` (arm) | ~1 s | 0.12 GiB |
| `vcc prep --dry-run` / `vcc prep` | 6.5 s / 15.25 s | 14.92 / 14.96 GiB |
| `effdiff255.py` / `metadiff255.py` / `price255.py` | 1.3 s / < 1 s / < 1 s | 0.64 / — / — GiB |
| `profile255.py` / `amp255.py` / `probe255.py` | 24.7 s / 40.0 s / 39.7 s | 3.59 / 1.18 / 13.31 GiB |

Disk: the arm's `.h5ad` + `.vcc` 15.5 GB under `submissions/`; no known-answer h5ad was built; `work/` holds two 44 MB
effect matrices, two CD4/union residual npz, 270 MB of PLOLO fold-cache copies, four per-target npz and three 28 MB
totals npz (all gitignored).

---

## 15. Files

| file | what it does |
|---|---|
| `repro255.py` | KA row 1 through 162's PLOLO with the production builder's code, one fold per process, + the k = 5 arm through the same builder (KA-C′) and the increment, + the energy-matched random control → `work/repro255-<fold>.{json,log}`, `work/pertarget255-<fold>.npz`; `--summary` → `work/repro255-summary.{json,log}` |
| `stage255.py` / `stage255.sh` | the arm's effect matrix (`svd10`) and KA row 2 (`svd5`), each alone and gated → `work/stage/` |
| `build255.sh` / `run255.sh` | 222's build flags + the one passthrough, every other flag unset; the gated serial sequence (arm build, postbuild, dry-run, prep) → `work/build/` |
| `effdiff255.py` | falsifiers 2/3 at the effect level, k = 10 vs k = 5 (and vs `h1union` for the record), deflation audits re-derived → `work/effdiff255.{json,log}` |
| `metadiff255.py` | falsifier 1 vs `249-cd4svd5-f034.h5ad.meta.json` → `work/metadiff255-vs249.{json,log}` |
| `profile255.py` | falsifiers 2/3 on the delivered counts (249's `profile249.py`) → `work/profile255.{json,log}` |
| `probe255.py` | falsifier 4, challenge `n_pred` (249's `probe249.py`) → `work/probe255.{json,log}`, `work/probe255_lfc.npz` |
| `amp255.py` | falsifier 5 + predicted raw `mse` at each package's own terms + the caps → `work/amp255.{json,log,stdout}`, `work/totals-{h1union,cd4svd5,cd4svd10}.npz` |
| `price255.py` | the six-member price of the arm and of the increment on 246 §4's card → `work/price255.{json,log}` |
| `work/ka2-rehash249.json` | KA row 2 (a): 249's staged k = 5 matrix re-hashed |
| `work/builder.intended.sha256` | HEAD's builder sha, the gate's reference |
| `PREREGISTRATION-DRAFT.md` | §12, for the lead |

---

## 16. Honest limits

- **The increment is measured on one instrument and it is inside its own noise.** `pds` +0.00054 raw is 0.36 of
  PLOLO's LOO σ and changes sign across the three folds; the cosine cost is the only thing the folds agree on. The
  board is the only instrument that can price the k 5 → 10 step, which is the entire reason the package exists.
- **The mandate's `energy_removed_frac ≈ 0.286` is 246's `k562_gwps`-leading entry**, not a mean: the four folds
  read 0.281–0.288 (mean 0.284) and production's own panel 0.279. The k = 10 row itself reproduces to 0.0e+00.
- **Amplitude is shape-only only relative to `cd4svd5`.** Against `h1union` the package sits at −0.86 %, beyond the
  0.5 % bound 249 set against that incumbent; a reader who prices the arm against `h1union` must use the
  amplitude-changing gate (+0.0089), which the k = 10 arm does not clear at 0.7× or 1.0× (+0.0021 / +0.0036).
- **The `mse` prediction rests on 197's identity with the PLOLO aggregate cosine as the only cosine input**, as 249's
  did; the class factor for a re-shape of occupied rows is 243's 1.36 and is not identified for this class.
- **No `cell-eval2`, no h1-test** (oracle-leaked on top of `h1union`), nothing submitted.

---

## 17. Verdict

**BUILDABLE.** `submissions/255-cd4svd10-f034.vcc` (`f860632a…`, 2,824,243,200 B) is built, prepped and gated: both
known-answer rows reproduce to the digit, the rule-10 diff is exactly `cd4_deflate_k` 5 → 10 plus its pre-declared
consequences, the profile moved on the 828 CD4-carried blocks and on nothing else, `n_pred` is unchanged, the
delivered amplitude is within 0.5 % of `cd4svd5`'s, and `vcc prep` is clean. It is NOT a slot candidate on its price
(the increment is −0.0005 … −0.0002 overall on 246's card and its `pds` is inside PLOLO's noise); it is the branch-A
instrument for the board's k-response, to be spent only as `E249/DECISION-TREE.md` pre-decided and only with the
draft preregistration committed first. **No model name, no description, nothing submitted.**
