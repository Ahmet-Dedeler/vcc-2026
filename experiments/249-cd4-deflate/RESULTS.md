# 249 — `cd4svd5` BUILT: the CD4 residual's top-5 target-shared directions deflated before sizing. Flag unset is byte-inert (staging `abeb0749…` and the shipped package `8df4dadd…` reproduced through the changed builder); the staged arm reproduces 246's `B2_cd4_svd5_resolved` to 0.0e+00 on every PLOLO fold; every pre-committed falsifier passes. **BUILDABLE — verdict line at the bottom (§16).**

**Date 2026-09-04.** Preregistration `PREREGISTRATION.md` (13:05 UTC) binding and untouched. **Nothing
submitted, `vcc submit` never invoked, `cell-eval2` never run, no model name chosen, no description
written.** Files edited outside this directory: exactly one, `scripts/make_submission_pm_cd4.py`
(commit **`be029c83d6321737b36cbade1c4464d6007ae41a`**, +133 / −0, file sha256 `2d63b979…`).
Packages written under `submissions/249-*` only. Every heavy stage ran alone, gated on
`memory_pressure` ≥ 40 % free AND `vm_stat` free+inactive ≥ 30 GiB AND zero other builders AND zero
`vcc prep` AND the working-tree builder at this experiment's bytes (§12 says why that last clause
exists); peak RSS per stage in §13.

| row / falsifier | required | measured | verdict |
|---|---|---|---|
| **KA 1** flag unset is byte-inert, staging | `effects_array_sha256` `abeb0749…` (222 `stage-on028.json`, re-verified by 247) | **`abeb0749…`** through the changed builder, **and** through the pre-249 byte copy with the flag set (unknown to it) | **PASS** |
| **KA 1** flag unset is byte-inert, package | `.h5ad` `8df4dadd…` on 12,648,807,800 B | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` on 12,648,807,800 B** (build249.sh, flag unset, 1345 s / 20.00 GiB; `def_postbuild` 32/32, `moved=[]`) | **PASS** |
| **KA 2** the staged arm reproduces 246 | macro −0.81 % / aggregate −0.46 % / Δ`pds` +0.01226, three cancer folds, to the digit | **−0.808675 % / −0.461748 % / +0.012263**, max\|Δ\| **0.0e+00** vs 246's rows on 16/16 fold×arm rows; energy removed 0.2008 (3-fold) / **0.2017** (4-fold = 246's "0.202"); CD4 fraction **0.3735 → 0.3960** (production 0.3221 → 0.3410) | **PASS** |
| **KA 3** energy-matched random control through the same builder path (one seed) | Δ`pds` ≤ +0.001 | **−0.001377** re-solved through the builder (−0.000129 at fixed sizing = 246's `B2_cd4_randE5` to 0.0e+00); 52/54/55 dims, 20.2–21.0 % energy, cosine −2.8 % | **PASS** |
| **F1** rule-10 audit diff | exactly `cd4_deflate_k` + the pre-declared CD4-frame consequences; `operator.sha256`, α 0.982, `cd4_target_f` 0.34, `h1_union_g` 0.28, every `effects.h1_union.*` unchanged | same-builder rebuild as reference: **config +1 (`cd4_deflate_k`), audit changed = the 11 pre-declared, added = the 5 deflation keys, 115/115 pinned (58/58 `effects.h1_union.*`), frozen inputs 0, unclassified 0**; vs the shipped 222 meta the same plus ten post-222 metadata fields 235/229 added to every build since (§5) | **PASS** |
| **F2** per-gene panel-total digest split by CD4-carried row | carried rows move in direction; non-carried rows one shared scalar | effect level: 276 carried rows cos ≤ 0.99985, 24 non-carried = incumbent × 1.0024537 to **9.9e−17**; delivered: 828 carried blocks cos median 0.9943 (min 0.672), 72 non-carried cos ≥ 0.99991, norm ratio 1.0021–1.0028, scale-corrected displacement 0.025 (carried min) vs 0.013 (non-carried max), **CLEAN** | **PASS** |
| **F3** non-vacuity | `‖ΔE‖/‖E‖` per block, zero-block count = 0 | delivered **0.153605**, per block min 0.006617 / median 0.102 / max 0.744, **0/900 at zero**, 0/900 < 1e−6, 900/900 digests differ; effect level 0.142049, 0/300 | **PASS** |
| **F4** challenge `n_pred` within ±2 % | ±2 % | median **2,599 / 2,599 / 2,599** for both packages, pooled **+0.00 %** | **PASS** |
| **F5** delivered amplitude off the package | `A_c`/`B_c` within 0.5 % of `h1union`'s | `√(A_arm/A_inc)` − 1 = **−0.409 %** (contexts ×0.99146 / 0.99214 / 0.99192), `B` −0.410 %; predicted board raw `mse` **0.8191** at k_cos 1 (0.8190–0.8194 across 243's classes) | **PASS** |
| **F6** `vcc prep --dry-run` clean; nnz; max counts/cell | clean; vs 4.75e9; vs 1e6 | ✓ clean (360,000 cells, A/B/C = 120,000 each, targets verified, counts-preserved); nnz **1,577,615,064 = 33.21 %** of cap; max **22,323** per cell (44.8× under) | **PASS** |

---

## 1. The change (`scripts/make_submission_pm_cd4.py`, commit `be029c83`, +133 / −0 on 247's `a17b177`; applied by `patch249.py`)

- **`VCC_CD4_DEFLATE_K`** — integer ≥ 1; unset = inert. `0`, `-3`, `5.0`, `abc` refused at import with a
  clear message. Read next to 247's `VCC_CD4_SIZING` pair.
- **`deflate_cd4_residual(cd4, k, *, basis=None)`** — 246's `lib246.deflate(res, "svd", k)` verbatim.
  On the carried rows (`target_condition_count > 0`: the 276 complete rows plus the 15 partial rows
  whose effects `freeze_complete_mean3` zeroes but whose counts it keeps — exactly the rows 162's
  `cd4_for` hands 246; 291 on A/B/C), subtract the per-gene mean over those rows, re-zero outside the
  gene support, `G = R_c R_cᵀ`, `np.linalg.eigh`, `U_k` = top-k eigenvectors, `R_k = (I − U_k U_kᵀ) R_c`,
  re-zero outside the support. Deterministic, no seed. Returns a `CD4Prepared` with the **same**
  `mode`, votes, counts and condition details (so `effects.cd4.mode`, `supported_entries`,
  `fully_supported_targets` and the count histogram cannot move) plus the audit
  `{k, n_carried, energy_removed_frac, top_k_energy_frac_of_centred, basis_u_k_sha256}`. `basis`
  exists for controls only (KA 3); production never passes it.
- **Where it acts** — one place, `apply_frozen_cd4_residual`, after the panel-coordinate assertions
  and **before** `solve_cd4_dose` (247's router; 099's all-rows rule underneath, **unmodified**), so
  the shipped sizing rule at `f` = 0.34 sizes the object actually added (246's *resolved* variant;
  the fraction rises because the deflated residual's median row norm is lower). The second residual
  (`extra`, the H1 union) is never touched; its `solve_carrying_row_fraction` reads only the whitened
  base and the union, so its solved fraction is bit-identical (§3, every fold; §2.1, production).
- **Meta** — `effects.cd4.deflation` and `config.cd4_deflate_k`, written only when on (the
  `h1_union_g` convention): a default build's config block and audit are byte-identical to every
  prior build's.

**Unit level** (`unit249.py`, two builders in one process — the working tree and the pre-249 byte
copy `work/make_submission_pm_cd4.PRE.py` = git HEAD `24a30055…`): **30/30**. OFF: output array and
both audits bit-identical to pre-249, with and without an H1-style second residual, no `deflation`
key. ON: the deflated residual bit-identical to `lib246.deflate(res, "svd", 5)` with identical
energy numbers; the fraction re-solves upward; the union's frame identical to OFF; rows CD4 does not
carry equal OFF's rows × one scalar to 1e−12; carried rows move in direction; a supplied basis is
applied as a projector (output orthogonal to `U` to 1.7e−15, Pythagoras to 4.9e−17 on a
full-support, no-partial-row residual — with partial rows the post-projection re-zeroing makes the
step non-idempotent, which is 246's construction too and the reason exactness needs full support);
a wrong-shaped basis is refused.

**Regression suites** under the changed builder: 235's hardening tests **9/9** (`work/tests235.log`,
1.07 GiB); the builder's own 57-test suite under both builders **56/57 and 56/57, 0 regressions, 0
fixed** (`suite249.py`; the one failure, `test_panel_exclusion_for_operator_and_cd4_residual`,
predates 222 — 247 §5.1).

---

## 2. Known-answer row 1 — the flag unset is byte-inert

### 2.1 Staging (`stage249.py` = 222's `stage222.py`; `stage249.sh`, each row alone and gated)

| row | builder | env | `effects_array_sha256` | `effects.cd4.fraction` | peak RSS |
|---|---|---|---|---|---:|
| `pre5` | pre-249 byte copy (HEAD `24a30055…`) | `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5` | **`abeb0749…` MATCH** (the flag does not exist there) | 0.32207201067578145 | 19.91 GiB |
| `off` | working tree `2d63b979…` | same, flag unset | **`abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315` MATCH** | 0.32207201067578145 | 19.92 GiB |
| `svd5` | working tree `2d63b979…` | + `VCC_CD4_DEFLATE_K=5` | `a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11` (the arm; = the built package's `effects.final_sha256`) | **0.34104860341941284** | 20.07 GiB |

So the only thing that moves the bytes is the new code responding to the flag. The `svd5` frame:
`residual_median_row_norm_at_unit_fraction` 4.464717385709141 → **4.2162920214222215**,
`f_per_unit_fraction` 1.0556645369046553 → **0.9969253548939966**, `solved_fraction`
0.32207201067578145 → **0.34104860341941284** (×1.05892; PLOLO's frame ×1.0592), the one merged
preserve scalar 0.9093569351340488 → **0.9115882187626339**; the union's `solved_fraction`
0.19207161964444322 and everything else in its block unchanged. `effects.cd4.deflation` =
`{k: 5, n_carried: 291, top_k_energy_frac_of_centred: 0.1982359547102955, energy_removed_frac:
0.19866046506100588, basis_u_k_sha256: f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c}`,
re-derived offline from the saved pre-deflation residual through the builder's function: identical
(`effdiff249.py`).

### 2.2 The full A/B/C package (`build249.sh` = 222's `build222.sh` + the flag passthrough, flag unset)

247's byte-identity build had been deleted after verification (247 §3.4), so it was rebuilt:
`submissions/249-ka-h1union-g028.h5ad` **sha256
`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6`, 12,648,807,800 bytes — MATCH**
(240 row 1 / 247 §3.4), 1344.96 s, peak RSS 20.00 GiB (`work/build/build-ka.log`);
`def_postbuild.py` **32/32**, config vs the 222 incumbent `moved=[]`; its meta records the builder at
`2d63b979…` with no `cd4_deflate_k` and no `deflation` block. The h5ad was then deleted as a byte
duplicate of `submissions/222-h1union-g028.h5ad`; the meta is kept
(`submissions/249-ka-h1union-g028.h5ad.meta.json`, copy in `work/build/`) and is the reference for
§5's exact diff.

---

## 3. Known-answer row 2 — the staged arm reproduces 246 (`repro249.py`, one fold per process)

162's PLOLO with the **production builder's code substituted in** (222's `repro222.py` pattern): the
builder's `deflate_cd4_residual` for `lib246.deflate`, and the whole of `apply_frozen_cd4_residual`
(deflate → size → merge → one preserve) for 246's `compose`, the builder imported with the arm's env
(`VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5`). Fold caches are
byte copies of 246's (`work/plolo/`). Tolerance 1e−9 on `cos_agg`, `a23_cos`, `pds`.

| fold | KA-A `STACK` vs 246 | KA-B builder deflation vs `lib246.deflate` | KA-C arm through `apply_frozen_cd4_residual` vs 246 `B2_cd4_svd5_resolved` | solved fraction (stack → arm) | union fraction moved | arm: macro / agg / Δ`pds` |
|---|---|---|---|---|---|---|
| `k562_gwps` | 0.0e+00 | bit-identical, audit identical (264 rows, E rm 0.202191245) | **0.0e+00**; effect matrix vs 246's compose max\|Δ\| 0.0 | 0.373761 → 0.395879 (= 246's to 0.0e+00) | 0.0e+00 | −1.78 % / −1.63 % / +0.009987 |
| `orion_hct116` | 0.0e+00 | bit-identical (282 rows, 0.200117535) | **0.0e+00**; 0.0 | 0.370782 → 0.394269 | 0.0e+00 | −1.91 % / −1.95 % / +0.009446 |
| `orion_hek293t` | 0.0e+00 | bit-identical (286 rows, 0.199984875) | **0.0e+00**; 0.0 | 0.376056 → 0.397869 | 0.0e+00 | +1.27 % / +2.20 % / +0.017355 |
| `kolf2` (apart) | 0.0e+00 | bit-identical (266 rows, 0.204597745) | **0.0e+00**; 0.0 | 0.316408 → 0.336448 | 0.0e+00 | −1.91 % / −3.47 % / +0.004139 |

Three cancer folds, against 246's own json rows recomputed the same way:

| | this run (builder) | 246 (json) | \|Δ\| | mandate |
|---|---:|---:|---:|---:|
| rel MACRO cosine | **−0.008086749** | −0.008086749 | 0.0e+00 | −0.81 % |
| rel AGGREGATE cosine | **−0.004617483** | −0.004617483 | 0.0e+00 | −0.46 % |
| Δ`pds` | **+0.012262757** | +0.012262757 | 0.0e+00 | +0.01226 |

`pds` positive on 3/3 cancer folds; Δq_w +0.00064. `energy_removed_frac` 0.2008 over the three
cancer folds, **0.2017 over four = 246 §3.3's "0.202"**; production's own panel **0.1987** over 291
carried rows (§2.1). CD4 solved fraction in PLOLO's frame **0.3735 → 0.3960** (246 §3.3: 0.374 →
0.396); production's frame 0.3221 → 0.3410 — the absolute differs because production's whitened base
and the un-held CD4 rows differ, the ratio (×1.0589 vs ×1.0592) does not. 16/16 known-answer rows at
0.0e+00; peak RSS 4.93–5.17 GiB per fold, 5.6–6.0 s each.

---

## 4. Known-answer row 3 — energy-matched random control through the same builder path

246's control, matched (`sweep246b.rand_energy_matched`: seed `20260904 + 5`, bisection on k′ until the
removed energy ≥ svd5's), the random orthonormal basis pushed through `deflate_cd4_residual(…,
basis=U)` and then through `apply_frozen_cd4_residual` itself (the builder's deflation call
substituted by the control basis; sizing, merge and preserve are the builder's own):

| | dims | energy removed | macro | agg | Δ`pds` (3 cancer folds) |
|---|---|---|---:|---:|---:|
| fixed sizing (= 246's `B2_cd4_randE5`) | 52 / 54 / 55 | 0.2097 / 0.2040 / 0.2015 | −2.594 % | −3.660 % | **−0.000129** (246: −0.000129, max\|Δ\| 0.0e+00; the subspace bit-identical to `lib246.deflate(…, "rand", k′, seed)`) |
| **re-solved through the builder** | same | same | −2.017 % | −2.821 % | **−0.001377** (per fold −0.00016 / −0.00541 / +0.00144) |
| the arm — the same energy, the top-5 directions | 5 | 0.2008 | −0.809 % | −0.462 % | **+0.012263** |

Removing the same ~20 % of the energy in random directions buys nothing and costs 3–6× the cosine:
**the builder applies a projector, not a norm cut.** Bound ≤ +0.001; PASS.

---

## 5. Falsifier 1 — the rule-10 audit diff (`metadiff249.py`, sets pre-declared in the source before any diff was read)

**Reference: the same-builder unflagged rebuild** (`submissions/249-ka-h1union-g028.h5ad.meta.json`,
byte-identical package to the incumbent, §2.2), `--same-builder`. Every difference, verbatim
(`work/metadiff249-samebuilder.log`; the `output.*` / `peak_rss_gib` bucket is reported, not asserted):

```
[A-config            ] ADDED   config.cd4_deflate_k = 5
[B-effects/operator  ] ADDED   effects.cd4.deflation.basis_u_k_sha256 = 'f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c'
[B-effects/operator  ] ADDED   effects.cd4.deflation.energy_removed_frac = 0.19866046506100588
[B-effects/operator  ] ADDED   effects.cd4.deflation.k = 5
[B-effects/operator  ] ADDED   effects.cd4.deflation.n_carried = 291
[B-effects/operator  ] ADDED   effects.cd4.deflation.top_k_energy_frac_of_centred = 0.1982359547102955
[B-effects/operator  ] CHANGED effects.cd4.fraction                                                 0.32207201067578145 -> 0.34104860341941284
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm     4.533854576152464   -> 4.5227571142920135
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.fraction                               0.32207201067578145 -> 0.34104860341941284
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.post_add_global_preserve_scale         0.9093569351340488  -> 0.9115882187626339
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.effective_f_at_cd4_fraction              0.05278322684523277 -> 0.04984626774469983
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.f_per_unit_fraction                      1.0556645369046553  -> 0.9969253548939966
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction 4.464717385709141  -> 4.2162920214222215
[B-effects/operator  ] CHANGED effects.cd4.residual_frame.solved_fraction                          0.32207201067578145 -> 0.34104860341941284
[B-effects/operator  ] CHANGED effects.final_sha256                                                abeb0749f8ecc941... -> a28d5be2e9cb3ed1...
[B-effects/operator  ] CHANGED effects.h1_union.merged_preserve.candidate_pre_preserve_median_norm 4.533854576152464   -> 4.5227571142920135
[B-effects/operator  ] CHANGED effects.h1_union.merged_preserve.post_add_global_preserve_scale     0.9093569351340488  -> 0.9115882187626339

A. CONFIG added ['cd4_deflate_k']  expected ['cd4_deflate_k'];  changed []  removed []
B. AUDIT changed == the 11 pre-declared;  added 5 keys, all under `effects.cd4.deflation.`;  unexpected []
   direction: residual_median_row_norm x0.944358 (down), f_per_unit x0.944358 (down), effective_f x0.944358 (down),
              solved_fraction / fraction / post_add_preserve.fraction x1.058920 (up);
              the merged preserve scalar is ONE number written twice: 0.9093569351340488 -> 0.9115882187626339
C. FROZEN INPUTS differences: 0      D. IMPLEMENTATION differences: [] (same builder bytes)      F. UNCLASSIFIED: []
PINNED scalars (115) -- each byte-identical: 115/115, of which 58/58 `effects.h1_union.*` fields
FALSIFIER 1 PASS
```

Pinned by name and byte-identical: `config.alpha` **0.982** · `config.cd4_target_f` **0.34** ·
`config.h1_union_g` **0.28** · `config.h1_union_sources` `['gse295214_h1', 'h1_2025']` ·
`operator.sha256` **`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`** ·
`operator.metadata_sha256` · `operator.columns_sha256` · `operator.axis_sha256` `edf6ba88…` ·
`effects.{baseline_pre,candidate_pre,post_whiten}_sha256` · all seven `effects.whitening.*`
(`after_all_median_norm` 4.229295604454299) · `effects.cd4.{mode 'gene_norm', fully_supported_targets
276, supported_entries 2,712,528, panel_residual_nonzero_entries 0, the count histogram
{0: 9, 1: 4, 2: 11, 3: 276}}` · `effects.cd4.residual_frame.{target_f_requested 0.34,
whitened_base_median_row_norm, cd4_fraction_constant, fraction_source}` ·
`effects.cd4.post_add_preserve.baseline_output_median_norm` · **every `effects.h1_union.*` field**
(the union's `solved_fraction` 0.19207161964444322, `carrying_targets` 77, `target_g` 0.28, source
shas, …) except the two merged-preserve scalars · `frozen_bundle.sha256` · `output.n_rows` 360,000 ·
the emitter's `rc_callset` triple and every other `config.*` key.

**Against the shipped `222-h1union-g028.h5ad.meta.json`** (the mandate's literal reference,
`work/metadiff249-vs222.log`): the same additions and the same 11 changes, 108/108 pinned (51/51
`effects.h1_union.*`), frozen inputs 0, implementation 1 (this file), unclassified 0 — **plus ten
added metadata fields that are not the arm's**: `operator.axis_sha256`,
`operator.axis_sha256_checked_against`, `operator.panel_sha256` (235's axis guard),
`effects.h1_union.source_coverage.{gse295214_h1, h1_2025}`,
`effects.h1_union.source_details[0|1].coverage_frac`, `effects.h1_union.targets_by_source_count.{0,1,2}`
(229's source registry). Measured, not argued: those ten keys are **absent** from the 222 and 229
metas (builders `71cf433d…` / `a8785cab…`) and **present with values identical to the arm's** in
235's and 247's byte-identical `8df4dadd…` rebuilds of `h1union` (builders `715f2b46…` /
`24a30055…`) — 10/10 in each. They are what every build since 235 writes; the 222 meta predates
them. Against 247's rebuild meta (`work/metadiff249-vs247ka1.log`) the only extra difference is the
driver's `config.expect_contexts = 'A,B,C'`, which build222.sh's form does not pass.

`def_postbuild.py` on the arm (`--allow-config-diff cd4_deflate_k`): **32/32**, `moved=['cd4_deflate_k']`.

---

## 6. Falsifier 2 — the profile moved on the CD4-carried rows and nowhere else

### 6a. On the effect matrix — exact (`effdiff249.py`, `effects-svd5.npy` vs `effects-off.npy`)

| row class | n | rel. displacement `‖ΔE_i‖/‖E_i‖` min / median / max | cosine(off, svd5) min / max | `‖E_svd5 − s·E_off‖/‖E_off‖` max |
|---|---:|---|---|---:|
| **CD4-carried** (support > 0) | 276 | 0.017954 / 0.099461 / 0.736029 | 0.6785760 / 0.9998508 | 7.38e−01 (median 9.96e−02) |
| **NOT carried** | 24 | **exactly 0.002454 for all 24** (= \|1 − s\|) | **1.0000000000** | **9.85e−17** |

`s` = 0.9115882187626339 / 0.9093569351340488 = **1.0024536939702957**, the merged preserve ratio; the
24 non-carried rows' norm ratio is exactly `s` (min = max to the last digit). Separation on the
scale-corrected displacement **1.8e+14×**. The carried rows hold 99.997 % of ‖ΔE‖²; the five most
displaced targets: MTF1 0.736, PLCG1 0.633, SMARCA5 0.441, DCP2 0.406, SOCS5 0.406.

### 6b. On the delivered counts (`profile249.py` = 222's `profile222.py`, carrying set = CD4 support)

900 (context, target) blocks, 828 carried (276 targets × 3), 72 not:

| block class | n | digest identical | cosine(arm, inc) min / median / max | norm ratio min / max | `‖ΔE − s·E‖/‖E‖` min / median / max |
|---|---:|---:|---|---|---|
| **carried** | 828 | 0 | 0.67243208 / 0.99425023 / 0.99969185 | 0.5926 / 1.0177 | **0.025027** / 0.107449 / 0.746426 |
| **NOT carried** | 72 | 0 | 0.99991230 / 0.99995901 / 0.99998110 | 1.00207773 / 1.00280166 | 0.006164 / 0.009077 / **0.013277** |

Shared delivered scale read off the non-carried blocks **1.00243568** (effect level 1.00245369).
Separation **CLEAN**: the smallest carried displacement is **1.9×** the largest non-carried one (222:
3.3×, under a 1.2 % scalar; here the scalar is 0.24 %, so the emitter's integer re-draw — 0.6–1.3 %
on those blocks — is a larger share of a smaller number). The exact statement is §6a's, at 1e−16.

**Digest counts, as asked: 0/900 identical (0/828 carried, 0/72 non-carried).** Nonzero-structure
digests 0/900 identical; panel-total sums 10/900 identical. **PASS.**

---

## 7. Falsifier 3 — non-vacuity. A zero would be the bug; there is none

| level | `‖ΔE‖/‖E‖` | per-block min / p01 / median / max | blocks exactly 0 | blocks < 1e−6 |
|---|---:|---|---:|---:|
| effect matrix (300 targets) | **0.142049** | 0.002454 / — / — / 0.736 | **0 / 300** | 0 / 300 |
| delivered profile (900 blocks) | **0.153605** | 0.006617 / 0.007777 / 0.102208 / 0.744456 | **0 / 900** | 0 / 900 |
| delivered, inside the 10-CPM gate | 0.153325 | 0.006620 / — / 0.101279 / 0.744864 | — | — |

Per context: A 0.155525, B 0.152306, C 0.152900. The residual applied on every block, in all three
contexts. 210's stricter 1 % per-block floor (which 222 inherited as its "FALSIFIER 4" print) holds
on all 828 carried blocks (min 2.5 % scale-corrected) and is not applicable to the 72 non-carried
blocks, whose displacement is by construction the 0.24 % preserve scalar plus the emitter's integer
re-draw — the script's own line reads FAIL for that reason and the mandate's criterion (zero-block
count = 0) is met. **PASS.**

For the record: median delivered `|row|` ratio 0.99881, mean 0.99636, `A = mean|row|²` **0.99207**
(= (1 − 0.0041)², §9); `rc_callset` head exposure: top-400 |delivered effect| overlap with the
incumbent min 106 / median 365 / max 400 of 400, implied call-set jaccard lower bound min 0.797 /
median 0.973; |Δ log2FC| over the DE universe max 2.71, mean 0.0064, signed block-mean median
+0.000006.

---

## 8. Falsifier 4 — challenge `n_pred` (`probe249.py` = 222's `probe222.py`; Wilcoxon + BH against the real `data/challenge/context_{A,B,C}.h5ad` controls, 40 perturbations per context)

| package | A | B | C | pooled median | ρ(−log p, \|H\|) in the call set, A / B / C (a15's key for both, cross-key by design) | call set UP |
|---|---:|---:|---:|---:|---|---|
| `h1union` | 2,599 [2,599, 2,634] | 2,599 [2,599, 2,639] | 2,599 [2,599, 2,662] | 2,599 | +0.9333 / +0.9304 / +0.9279 (mean +0.9306) | 49.8 / 50.0 / 48.1 % |
| **`cd4svd5`** | **2,599** [2,599, 2,635] | **2,599** [2,599, 2,638] | **2,599** [2,599, 2,660] | **2,599** | +0.9309 / +0.9270 / +0.9239 (mean +0.9273) | 49.7 / 49.8 / 48.0 % |

**+0.00 %** against ±2 %. The `rc_callset` prescription pins `n_pred` by construction and it held.
**PASS.** (Peak RSS 13.32 GiB — 222's instrument at its known cost, run alone after the last build.)

---

## 9. Falsifier 5 — delivered amplitude off the package, and the predicted board raw `mse` (`amp249.py`)

Per-block panel totals streamed off both packages (214's `block_totals` via 243's `totals243.py`),
211's per-context terms via 243's `terms`. **Known-answer row:** `h1union`'s re-streamed
(`A_c`, `B_c`, `cap_c`) equal 243's `invert243.json` to **0.0e+00** on all three contexts, and 197's
identity at those terms reproduces its board raw `mse` **0.8172630** to 0.0e+00.

| context | `A_c` arm | `A_c` h1union | ratio | `B_c` arm / h1union | cap | leg |
|---|---:|---:|---:|---|---:|---|
| A | 0.231633 | 0.233629 | ×0.99146 | 0.229698 / 0.231677 | 0.233754 | B |
| B | 0.220958 | 0.222708 | ×0.99214 | 0.219014 / 0.220755 | 0.239678 | B |
| C | 0.231016 | 0.232897 | ×0.99192 | 0.229028 / 0.230899 | 0.266129 | B |

**Delivered amplitude `mean_c √(A_arm/A_inc)` − 1 = −0.409 %** (`B` term −0.410 %); bound ±0.5 %:
**PASS**, at 82 % of the bound. 246 had predicted 1.0002 / 1.0033 = −0.31 % from the effect-space
delivered ratio; the package says −0.41 %. All three contexts stay on the `B` leg.

Predicted board raw `mse` from 197's identity **at the arm's own terms** (243 §5.2: never the
floor), with `cos_arm = cos_h1union × (1 + k_cos × rel_agg)`, `rel_agg` = −0.4617 % (§3) and
`cos_h1union` = 0.107332 (243's measured):

| k_cos | arm cos | Δcos | **predicted raw `mse`** | Δ raw | `mse` leg, overall |
|---:|---:|---:|---:|---:|---:|
| 0 (amplitude/budget leg alone) | 0.107332 | 0 | 0.81809 | +0.00083 | −0.00015 |
| 0.94 (new-source class) | 0.106867 | −0.000466 | **0.81902** | +0.00176 | −0.00031 |
| 1.00 | 0.106837 | −0.000496 | **0.81908** | +0.00182 | −0.00032 |
| 1.36 (re-dose of occupied rows — this arm's class by 243 §4.4) | 0.106658 | −0.000674 | **0.81943** | +0.00217 | −0.00038 |

The preregistration's "≈ 0.8183" was k_cos = 1 with `h1union`'s amplitude terms; with the arm's own
terms the −0.41 % amplitude adds +0.00083 raw, so the reading table for the board is the column
above. Preregistration falsifier 2 (raw `mse` worse by more than +0.0030) fires above **0.82026**;
the class prediction sits +0.0008 below it.

---

## 10. Falsifier 6 — `vcc prep --dry-run`, nnz, max counts per cell

```
✓ dry run — no .vcc written
  input:  submissions/249-cd4svd5-f034.h5ad   cells: 360000   genes: 18533 (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'   contexts: column 'context' → A=120000, B=120000, C=120000
  targets: verified against the official list   normalization: counts-preserved
```

| cap | limit | measured | margin |
|---|---:|---:|---:|
| stored entries in `.X` | 4,750,000,000 | **1,577,615,064** (`h1union` 1,577,411,397, +0.013 %) | **33.21 %** of cap |
| counts per cell | 1,000,000 | **22,323** (min 17,695; identical to the incumbent's, the exact-total emitter) | 44.8× under |
| rows | 360,000 | 360,000 | exact |
| cells per (context, target) block | 400 | 400 on all 900 | exact |

Dry run 19.7 s / 21.48 GiB; the real prep 25.2 s / 21.21 GiB. **PASS.**

---

## 11. The package

| artifact | bytes | sha256 |
|---|---:|---|
| **`submissions/249-cd4svd5-f034.vcc`** | **2,824,683,520** | **`6e8b41c6244b32dbe2a9765c7d2f2293b2032c49e59b455982b95bc6bb074533`** |
| **`submissions/249-cd4svd5-f034.h5ad`** | **12,648,807,800** | **`aba2775836df3bcf9619a520809fc7a71f325c0a81a881e6f9eb16c970f65faa`** |
| `submissions/249-cd4svd5-f034.h5ad.meta.json` | 25,633 | (`output.sha256` = the h5ad's; `effects.final_sha256` `a28d5be2…` = the staged `svd5` matrix; builder `2d63b979…` = commit `be029c83`) |
| `submissions/249-ka-h1union-g028.h5ad.meta.json` (KA row 1; its h5ad, a proven byte duplicate of `222-h1union-g028.h5ad`, was deleted) | 25,335 | `output.sha256` `8df4dadd…` |

Configuration (the arm is one variable): `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28
VCC_CD4_DEFLATE_K=5`, `VCC_UNION_SOURCES` unset (the H1 pair), `--min-cpm 10 --tail-spread 1.0
--keep-top 8000 --keep-spread 0.25 --overdispersion 4.0 --whitener-id union_pruned --emitter
rc_callset`, operator `e2f65ecf…`. **No model name chosen, no description written, nothing
submitted.** §14 applies at submit time; the preregistration's §6 places it in the second
2026-09-05 UTC slot, after the armed 237 reading, never in flight together with it.

---

## 12. Concurrency, honestly recorded

244's October rehearsal (`experiments/244-october-rehearsal/`) was running the driver end to end on
238's essential-free bundle while this experiment worked, and two things followed from sharing
`scripts/make_submission_pm_cd4.py` with it:

1. **244's rehearsal build ran on this experiment's modified builder.** Its `[build]` started at
   16:12:45 local, nine seconds after this experiment's last edit to the file; the flag is unset in
   the driver's environment, so 244 exercised the byte-inert path (§1's unit rows, §2.1's `off` row
   and §2.2's package are the proof that the path is inert). 244's meta records the uncommitted
   builder sha `fbdcbaf2…` — this experiment's bytes plus 244's own registry pins — under
   `hashes.implementation`; that is the reason, and it is 244's to note.
2. **244's `pin-artifacts` stage edited the same file mid-run** (`def_pin_panel_artifacts.py --apply`:
   a `union_pruned_def` registry entry and 12 replaced CD4 bank hashes in `FROZEN_INPUT_SHA256`), on
   top of this experiment's change, and restored it afterwards from its own saved copy — whose sha
   (`2d63b979…`) is exactly this experiment's intended builder (`intended249.py` reconstructs it from
   git HEAD + `patch249.py` + the one follow-up edit). With the pins present an A/B/C build would
   have failed the frozen-input check, so every heavy stage here additionally gated on the
   working-tree builder hashing to `2d63b979…` (`work/builder.intended.sha256`), and both build
   scripts' gates also count `vcc prep` processes (an agent's own command line quotes the words, so
   `claude` lines are excluded). No stage of this experiment ran while 244's build or prep was live;
   the staging rows started at 13:35:54Z, 90 s after 244's build exited, and the commit was made
   only after the working tree was verified at `2d63b979…`.

---

## 13. Cost and memory (every ≥ 15 GiB stage alone; the gate logged before each)

| stage | wall | peak RSS |
|---|---:|---:|
| `unit249.py` / `suite249.py` / 235's tests | 2 s / 58 s / ~2 min | 0.21 / 0.32 / 1.07 GiB |
| `repro249.py` × 4 folds | 5.6–6.0 s each | 4.93 / 5.16 / 5.17 / 5.08 GiB |
| staging `pre5` / `off` / `svd5` | 8.8–8.9 s each after inputs | 19.91 / 19.92 / 20.07 GiB |
| ARM build (`build249.sh`, `VCC_CD4_DEFLATE_K=5`) | **1349.8 s** | **20.16 GiB** |
| `def_postbuild.py` (arm / KA) | ~1 s | 0.12 GiB |
| `vcc prep --dry-run` / `vcc prep` | 19.7 s / 25.2 s | 21.48 / 21.21 GiB |
| KA build (`build249.sh`, flag unset) | **1345.0 s** | **20.00 GiB** |
| `effdiff249.py` / `amp249.py` / `profile249.py` / `probe249.py` | 3 s / 15 s / 27 s / ~4 min | 0.3 / 1.11 / 3.61 / 13.32 GiB |

Disk: the arm's `.h5ad` + `.vcc` 15.5 GB under `submissions/`; the KA h5ad (12.6 GB) deleted after
verification; `work/` holds three 44 MB effect matrices, three CD4/union residual npz, 270 MB of
PLOLO fold-cache copies and two 28 MB totals npz (all gitignored).

---

## 14. Files

| file | what it does |
|---|---|
| `patch249.py` | the ONE builder change as four exact-anchor replacements (the working tree, or any byte copy of HEAD) |
| `intended249.py` | HEAD + `patch249.py` + the follow-up edit → `work/make_submission_pm_cd4.249.py` and its sha, the gate's reference |
| `unit249.py` | §1's 30 unit rows: working tree vs the pre-249 byte copy vs `lib246.deflate` → `work/unit249.{json,log}` |
| `suite249.py` | the builder's 57 tests under both builders → `work/suite249.{json,log}`; 235's tests → `work/tests235.log` |
| `repro249.py` | KA rows 2 and 3 through 162's PLOLO, one fold per process → `work/repro249-<fold>.{json,log}`; `--summary` → `work/repro249-summary.{json,log}` |
| `stage249.py` / `stage249.sh` | KA row 1 at staging level + the arm's effect matrix → `work/stage/` |
| `build249.sh` / `run249.sh` | 222's build flags + the flag passthrough; the gated serial sequence (arm build, postbuild, prep, KA build) → `work/build/` |
| `metadiff249.py` | falsifier 1 (`--same-builder`, `--tag`) → `work/metadiff249-{samebuilder,vs222,vs247ka1}.{json,log,stdout}` |
| `effdiff249.py` | falsifiers 2/3 at the effect level → `work/effdiff249.{json,log}` |
| `profile249.py` | falsifiers 2/3 on the delivered counts (222's `profile222.py`) → `work/profile249.{json,log}` |
| `probe249.py` | falsifier 4 (222's `probe222.py`) → `work/probe249.{json,log}` |
| `amp249.py` | falsifier 5 and 6's caps; 197's identity at the arm's own terms → `work/amp249.{json,log}`, `work/totals-{cd4svd5,h1union}.npz` |
| `work/make_submission_pm_cd4.PRE.py` | the git control, HEAD `24a30055…`, taken before the first edit |
| `work/commit-message.txt` | the builder commit's message with the known-answer hashes |

```bash
cd ~/Code/virtual-cell-challenge-2026
git show HEAD:scripts/make_submission_pm_cd4.py > experiments/249-cd4-deflate/work/make_submission_pm_cd4.PRE.py
uv run python experiments/249-cd4-deflate/patch249.py            # + the follow-up edit; intended249.py reproduces both
uv run python experiments/249-cd4-deflate/intended249.py
uv run python experiments/249-cd4-deflate/unit249.py
uv run python experiments/249-cd4-deflate/suite249.py
uv run python experiments/235-october-hardening/tests/test_october_hardening.py
cp experiments/246-direction-hypotheses/work/plolo/*.np* experiments/249-cd4-deflate/work/plolo/
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do uv run python experiments/249-cd4-deflate/repro249.py --fold $f; done
uv run python experiments/249-cd4-deflate/repro249.py --summary
bash experiments/249-cd4-deflate/stage249.sh                       # pre5 / off / svd5, gated
uv run python experiments/249-cd4-deflate/effdiff249.py
bash experiments/249-cd4-deflate/run249.sh                         # arm build, postbuild, prep, KA build, gated
uv run python experiments/249-cd4-deflate/metadiff249.py submissions/249-cd4svd5-f034.h5ad.meta.json submissions/249-ka-h1union-g028.h5ad.meta.json --same-builder --tag metadiff249-samebuilder
uv run python experiments/249-cd4-deflate/metadiff249.py submissions/249-cd4svd5-f034.h5ad.meta.json submissions/222-h1union-g028.h5ad.meta.json --tag metadiff249-vs222
uv run python experiments/249-cd4-deflate/profile249.py submissions/249-cd4svd5-f034.h5ad submissions/222-h1union-g028.h5ad
uv run python experiments/249-cd4-deflate/amp249.py
uv run python experiments/249-cd4-deflate/probe249.py h1union=submissions/222-h1union-g028.h5ad cd4svd5=submissions/249-cd4svd5-f034.h5ad --perts=40
```

---

## 15. Honest limits

- **No board member was measured and no `cell-eval2` was run.** KA row 2 is a reproduction of 246's
  PLOLO instrument through the builder, not new evidence that the `pds` gain transfers; the
  preregistration's §5 says what the slot buys and its falsifiers 1–3 read off the board.
- **The delivered amplitude moved −0.41 %**, inside the 0.5 % bound but at 82 % of it, and larger
  than 246's −0.31 % effect-space estimate. On its own it costs +0.00083 raw `mse` (−0.00015
  overall) before any cosine change; the class prediction for the board's raw `mse` is therefore
  0.8190–0.8194, not the preregistration's 0.8183.
- **The delivered-level carrying/non-carrying separation is 1.9×** (222: 3.3×) because the shared
  scalar is 0.24 % and the emitter's integer re-draw is the larger term on the 72 non-carried
  blocks; the exact statement lives at the effect level (1e−16). Making the delivered statement
  exact would need a third build of the incumbent's effects scaled by 1.00245 — 22 minutes for a
  check §6a already settles.
- **`h1-test` is oracle-leaked** for anything on top of `h1union` (222 §11) and was not used.
- **`probe249.py` peaked at 13.3 GiB**, above the ~10 GB guideline, as 222's did (13.2); it ran alone
  after the last build, with 32 GiB free+inactive.
- **The `nmae` sign is unread here**: 246 §8.1d's exact identity gives +0.00031 as a magnitude
  bound; nothing on disk can sign it (222 §8).
- **k = 5 was not re-optimised** and the `+partial` rider (246 §3.7) was not built: one variable.

---

## 16. Verdict

**BUILDABLE.** `submissions/249-cd4svd5-f034.vcc` (`6e8b41c6…`, 2,824,683,520 B) is the preregistered
`cd4svd5` arm: one variable (`VCC_CD4_DEFLATE_K=5`) on top of the shipped `h1union` configuration,
the builder change committed as `be029c83` with the flag unset proven byte-inert at staging
(`abeb0749…`) and package (`8df4dadd…`) level, the staged arm equal to 246's `B2_cd4_svd5_resolved`
to 0.0e+00 on every fold, the energy-matched random control at Δ`pds` −0.0014, and all six
falsifiers passed with their numbers above (config diff exact; carried rows move, non-carried rows
one scalar; 0/900 zero blocks; `n_pred` +0.00 %; amplitude −0.41 %; prep clean, 33.2 % of the
density cap). Not submitted.
