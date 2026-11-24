# 254 — `cd4svd5+partial` BUILT: CD4's partial-condition rows ride as their own residual (`VCC_CD4_PARTIAL_G`, inert by default) on top of the k = 5 deflation. Flag unset is byte-inert (staging `abeb0749…`, 249's `a28d5be2…` with `K=5` alone, the shipped package `8df4dadd…`); 246's `C2_partial@0.34` and `E_cd4_svd5+partial@0.34` reproduce through the builder to 0.0e+00 on every fold; the package the builder actually carries (the RE-SOLVED deflation + the partial rows) is measured separately and its increment over `cd4svd5` is the partial row itself, **+0.0009 … +0.0018 overall, under the +0.0019 gate** — as 246 said. Every pre-committed falsifier passes. **BUILDABLE — verdict line at the bottom (§17).**

**Date 2026-09-04 (17:00–18:07Z).** Deterministic. **No `vcc submit`, no `cell-eval2`, no model name, no
description, no slot proposal.** One file changed outside this directory, `scripts/make_submission_pm_cd4.py`
(+199 / −10 on 251's `c321f76` bytes `10d0177e…` → `3615d20f…`), committed as **`32c5219`** with the
known-answer hashes in the message (`work/commit-message.txt`). Packages written under `submissions/254-*`
only. Every ≥ 15 GiB stage ran alone under the mandate's gate (`memory_pressure` ≥ 40 % free AND `vm_stat`
free+inactive ≥ 30 GiB AND zero other builders AND zero `vcc prep`, the working-tree builder hashed to this
experiment's intended bytes before each row); peak RSS per stage in §14. 252 and 253 were live throughout
(neither touches the builder; 253 rewrote `docs/READ-FIRST.md` to v2 at 17:00Z, which this run read after the
fact and which changes nothing here).

| row / falsifier | required | measured | verdict |
|---|---|---|---|
| **KA 1** flag unset is byte-inert, staging | `effects_array_sha256` `abeb0749…` (222 `stage-on028.json`) | **`abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315`** through the changed builder, **and** through the pre-254 byte copy with `VCC_CD4_PARTIAL_G=0.34` set (unknown to it); `VCC_CD4_DEFLATE_K=5` alone still gives 249's **`a28d5be2…`** (the deflation path untouched) | **PASS** |
| **KA 1** flag unset is byte-inert, package | `.h5ad` `8df4dadd…` on 12,648,807,800 B (251 rebuilt it at 16:21Z; its verified copy was deleted on match, so it was rebuilt once more through THIS builder) | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` on 12,648,807,800 B — MATCH** (`run254.sh`, both flags unset, 1331.63 s / 20.04 GiB; `def_postbuild` 32/32, `moved=[]`) | **PASS** |
| **KA 2** `VCC_CD4_PARTIAL_G=0.34` alone reproduces 246's `C2_partial@0.34` | macro +0.587 % / `pds` +0.00178, three cancer folds, to the digit | **+0.586691 % / +0.001777069** (aggregate +1.156392 %), max\|Δ\| **0.0e+00** vs 246's rows on every fold; the partial residual == `lib246.cd4_partial` bit for bit on every fold AND on the production panel; its fraction == 246's; CD4's and the union's fractions untouched; effect matrix == 246's `compose` at max\|Δ\| 0.0; leave-top-3-out −6.6 % (246: −6.6 %) | **PASS** |
| **KA 3 (a)** `K=5` + `G=0.34` reproduces 246's `E_cd4_svd5+partial@0.34` | `pds` +0.01405, macro −0.699 %, leave-top-3-out 61.5 % | 246 composed its E arm at the STACK's FIXED CD4 fraction (`sweep246b.py`: `[(d, st["frac_cd4"]), …]`), where the builder's `K=5` path RE-SOLVES the fraction (249's resolved variant = the shipped `cd4svd5`). Through the builder with `solve_cd4_dose` pinned to the un-deflated residual (= the fixed fraction): **macro −0.698529 %, `pds` +0.014045983, leave-top-3-out 61.54 %**, max\|Δ\| **0.0e+00** vs 246 on every fold, effect matrix == 246's `compose` at 0.0 | **PASS (246's construction reproduced to the digit; §3)** |
| **KA 3 (b)** the builder's REAL path (the package's object) | measured, not mandated | resolved deflation + partial: **macro −0.222 %, aggregate +0.702 %, `pds` +0.014046 (3/3), leave-top-3-out 61.7 %**; `K=5` alone through this builder == 246 `B2_cd4_svd5_resolved` at 0.0e+00; **increment over `cd4svd5`: macro +0.586 %, aggregate +1.162 %, `pds` +0.001783 (3/3), leave-top-3-out −6.3 %** | **MEASURED** |
| **F1** rule-10 diff vs `249-cd4svd5-f034.h5ad.meta.json` | exactly `cd4_partial_g` + its consequences | config **+1 (`cd4_partial_g` = 0.34)**, changed 0; audit changed = `effects.final_sha256` + the ONE merged preserve scalar and its norm written twice (0.9115882187626339 → 0.9101933222360135); added = 33 keys, all under `effects.cd4_partial.`; **128/128 pinned** (`config.cd4_deflate_k` 5, `effects.cd4.fraction` 0.34104860341941284, the CD4 frame, the deflation block, `operator.sha256`, α, `cd4_target_f`, `h1_union_g`, every `effects.h1_union.*` field but the preserve pair); frozen inputs 0; implementation = this file; unclassified 0. Same-builder KA rebuild as reference: config +2 (`cd4_deflate_k`, `cd4_partial_g`), audit changed = 249's 11 pre-declared keys, added 38 under `effects.cd4.deflation.` and `effects.cd4_partial.`, 116/116 pinned, implementation 0, PASS | **PASS** |
| **F2** profile split by partial row | the 15 partial rows move in direction; the rest one shared scalar | effect level vs `cd4svd5`: 15 rows cos 0.897–0.994, 285 rows == `cd4svd5` × 0.9984698172947937 to **1.1e−16**, separation 1.0e+15×; delivered: 45 carrying blocks cos median 0.9592 (min 0.880), norm ratio up to 1.179; 855 not-carrying blocks cos ≥ 0.99994, shared scale 0.99847129 (effect level 0.99846982); scale-corrected displacement 0.1258 (carrying min) vs 0.0107 (not-carrying max), **CLEAN**, 11.7× | **PASS** |
| **F3** non-vacuity | zero-block count = 0 | effect level ‖ΔE‖/‖E‖ 0.070779, 0/300 rows at zero; delivered ‖ΔE‖/‖E‖ 0.077454, per block min 0.004012 / median 0.007491 / max 0.535259, **0/900 at zero**, 900/900 digests differ | **PASS** |
| **F4** challenge `n_pred` within ±2 % | ±2 % vs `cd4svd5` | median **2599 / 2599 / 2599** vs `cd4svd5`'s 2599 / 2599 / 2599, pooled **+0.00 %** | **PASS** |
| **F5** delivered amplitude off the package | within 0.5 % of `cd4svd5`'s via 243's stream | `√(A_arm/A_cd4svd5)` − 1 = **+0.227 %** (contexts ×1.00481 / 1.00435 / 1.00446), `B` +0.228 %; vs `h1union` -0.183 %; predicted board raw `mse` **0.81610** at k_cos 1 (-0.00116 vs `h1union`'s 0.8172630, -0.00297 vs `cd4svd5`'s predicted 0.81908) | **PASS** |
| **F6** `vcc prep --dry-run` clean; nnz; max counts/cell | clean; vs 4.75e9; vs 1e6 | ✓ clean (360,000 cells, A/B/C = 120,000 each, targets verified, counts-preserved); nnz **1,577,695,797** = 33.21 % of cap (`cd4svd5` 1,577,615,064, +0.005 %); max **22,323** per cell; prep 23.85 s / 24.51 GiB | **PASS** |

---

## 1. The change (`scripts/make_submission_pm_cd4.py`, commit `32c5219`, +199 / −10 on 251's `c321f76`; `patch254.py`, eleven exact-anchor replacements, `--check` reproduces the working tree from the HEAD byte copy)

- **`VCC_CD4_PARTIAL_G`** — a number in (0, 5); unset or blank = inert. `0`, `5`, `-0.3`, `abc` refused at
  import (`envcheck254.py` **11/11**, incl. the pre-254 byte copy importing with the flag set and carrying no
  attribute, and the flag alongside `VCC_CD4_DEFLATE_K=5` or `VCC_CD4_DEFLATE_BASIS=V5_abc` — both accepted:
  the October composition needs the pinned basis plus this flag).
- **`extract_cd4_partial(prepared)`** — 246's `lib246.cd4_partial(targets, equalise=False)` verbatim on the
  pre-freeze `prepare_cd4_mean` object: rows with `target_condition_count` in {1, 2} keep their effect (already
  the mean over the conditions that carry them, `gene_norm`, at the banks' own rescale — E062's per-condition
  multiplier is a property of the whole bank, so a partial row carries exactly the complete rows' rescale) and
  votes; every other row zeroed; counts kept for the partial rows, 0 elsewhere; `condition_details[0]` = the
  row census and the complete rows' median norm. No re-normalisation (246's `equalise=True` read less `pds`).
- **`cd4_partial_support_audit(partial, target_g=)`** — the meta block's static half (row census, histogram,
  supported entries, the complete rows' median norm, the rule, the sizing rule, `target_g`).
- **`prepare_frozen_cd4(..., return_partial=False)`** — the default call performs exactly the operations it
  always performed; `return_partial=True` (only when the flag is on) also returns the partial residual
  extracted from the SAME prepared object the complete-row residual is frozen from.
- **`apply_frozen_cd4_residual(..., partial=None, partial_audit=None)`** — the partial residual is the LAST
  entry of the accumulation order (union 1, union 2, then the partial rows — 246's `merged_many` order for its
  E arm), sized by `solve_carrying_row_fraction` at `CD4_PARTIAL_G` against the same whitened base (the all-rows
  rule would solve exactly 0.0 for a 15-row residual), ONE add + ONE preserve. Its block comes back through
  `audits_out["cd4_partial"]` (`h1_union`'s fields plus the row census). CD4's own solve, its deflation
  (249/251) and the unions never see it; the only CD4-block scalars it can move are the merged preserve pair.
- **`build_production_effects(..., cd4_partial=, cd4_partial_audit=)`** writes `effects.cd4_partial` only when
  on; **`run()`** prepares it only when on and threads it through BOTH effect builds (the main one and
  `VCC_ALPHA=knee`'s per-context rebuild); **`config.cd4_partial_g`** written only when on (the `h1_union_g`
  convention): a default build's config block and audit are byte-identical to every prior build's.

**Unit level** (`unit254.py`, two builders in one process — the working tree and the pre-254 byte copy
`work/make_submission_pm_cd4.PRE.py` = git HEAD `10d0177e…`): **44/44**. OFF: output array and every audit
bit-identical to pre-254 with no extra, with the H1-style extra, and with a union 2; no `cd4_partial` key
anywhere; `return_partial` defaults to False. ON: `extract_cd4_partial` == lib246's arithmetic bit for bit
(arrays, counts, details); the partial's support disjoint from the frozen object's; the output == 246's
`merged_many` chain (cd4·f + union·f₁ + partial·f_p, votes OR'd, counts max, ONE add, ONE preserve) **bit for
bit** (accumulating the partial BEFORE the union agrees only to 1e−12 — order is a float association, which
is why the builder uses 246's order); CD4's fraction and frame and the union's frame identical to OFF; the
preserve scalar moves and is ONE number in all three blocks; rows the partial does not carry == OFF × one
scalar (≤ 1e−12); partial rows move in direction; with `K=5` as well the CD4 fraction and deflation audit
equal the K5-only run's and the output == `merged_many` at the RE-SOLVED fraction; a partial without a dose,
without `audits_out`, with a panel coordinate, or with no carrying rows is refused (loud, not silent).

**Regression suites** under the changed builder: the builder's 57-test suite **56/56 under both builders, 0
regressions, 0 fixed** (`suite254.py`; the one failure, `test_panel_exclusion_for_operator_and_cd4_residual`,
predates 222 — 247 §5.1); 235's hardening tests **9/9** (`work/tests235.log`).

---

## 2. Known-answer row 1 — the flag unset is byte-inert

### 2.1 Staging (`stage254.py` = 249's `stage249.py` + the partial residual; `stage254.sh`, each row alone and gated; `work/stage/stage254.log`)

| row | builder | env | `effects_array_sha256` | `effects.cd4.fraction` | partial | peak RSS |
|---|---|---|---|---|---|---:|
| `pre` | pre-254 byte copy (HEAD `10d0177e…`) | shipped arms + `VCC_CD4_PARTIAL_G=0.34` | **`abeb0749…` MATCH** (the flag does not exist there) | 0.32207201067578145 | — | 19.95 GiB |
| `off` | working tree `3615d20f…` | shipped arms, flags unset | **`abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315` MATCH** | 0.32207201067578145 | ABSENT | 19.95 GiB |
| `k5` | working tree | + `VCC_CD4_DEFLATE_K=5` | **`a28d5be2e9cb3ed17d11a4c5f40f470893f0c968a1ae87b0b3b358321ead5d11` = 249's staged arm** (the deflation path untouched) | 0.34104860341941284 | ABSENT | 17.39 GiB |
| `partial` | working tree | + `VCC_CD4_PARTIAL_G=0.34` | `192c2aed432c296849327c6b578cd90ac38345e4ede682ca2de47383c5e95e3a` (row 2's object) | 0.32207201067578145 (untouched) | 15 rows, fraction 0.11244119388619626 | 19.99 GiB |
| `composed` | working tree | + `VCC_CD4_DEFLATE_K=5 VCC_CD4_PARTIAL_G=0.34` | **`885f6d8c8f865b3d60ecf3ca078afe01c344b81c4376872b1e82942e6a0ce1a6`** (the arm; = the built package's `effects.final_sha256`) | 0.34104860341941284 (= `k5`'s) | 15 rows, fraction 0.11244119388619626 | 20.15 GiB |

The only thing that moves the bytes is the new code responding to the flag. The `composed` row's
`effects.cd4_partial` block: rows **15** (4 one-condition, 11 two-condition: ANKRD52, EAPP, GMEB2, KMT2A,
PDK1, PIN1, PLAGL2, PRMT7, SAP30L, STK3, TAF4, TIMM50, TMF1, UPF3B, ZKSCAN2), 147,420 supported entries,
carrying-row median norm **12.788556007061318** (the complete rows' 4.522043409698855 — a one-condition row is
not shrunk by averaging), `g_per_unit_fraction` 3.0238028274950555, **`solved_fraction`
0.11244119388619626**, `shipped_median_rule_would_solve` **0.0** (5 % coverage), `entries_outside_cd4_support`
147,420 (all of them: the frozen object has no votes on these rows), `targets_outside_cd4_support` 0 (the
frozen object keeps the counts). The deflation block is byte-identical to `k5`'s (`U_5` `f06e0262…`, energy
removed 0.19866046506100588); the union's `solved_fraction` 0.19207161964444322 in every row; the merged
preserve scalar 0.9115882187626339 (`k5`) → **0.9101933222360135** (`composed`), 0.9093569351340488 (`off`) →
0.908116101919572 (`partial`).

**The production-panel partial residual is 246's** (`partialcheck254.py`, **10/10**): the staged
`partial-composed.npz` == `partial-partial.npz` == `lib246.cd4_partial(300 production targets,
equalise=False)` bit for bit (arrays, counts, census, the complete rows' median norm); the staged complete-row
object == `plolo.cd4_for(production targets)` bit for bit; the partial's solved fraction re-derives from the
saved residual and the frame's recorded base norm (4.229295604454299, the shipped) to 1e−15.

### 2.2 The full A/B/C package (`run254.sh` → `build254.sh` = 222's `build222.sh` + the two passthroughs, both unset, gated, alone; `work/build/run254.log`)

Both flags unset through the changed builder (commit `32c5219`, builder `3615d20f…`), 222's build flags verbatim, strictly serial and alone (the gate opened at free+inactive 39.0 GiB, `memory_pressure` 72 % free, zero other builders, zero preps):

| what | measured | required |
|---|---|---|
| `submissions/254-ka-h1union-g028.h5ad` sha256 / bytes | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` / 12,648,807,800 — MATCH** | `8df4dadd…` / 12,648,807,800 (the shipped `h1union`; 240 row 1, 247 §3.4, 249 §2.2, 248 §2.2, 251 §2.2) |
| `config.cd4_partial_g` / `config.cd4_deflate_k` / `effects.cd4_partial` / `effects.cd4.deflation` | **absent / absent / absent / absent** | absent (the `h1_union_g` convention) |
| `def_postbuild.py` against the 222 incumbent meta | **32/32 checks pass, `moved=[]`** | — |
| `hashes.implementation` for the builder | `3615d20f…` (this commit's bytes) | changed, as 235's / 247's / 249's / 248's / 251's did; the h5ad did not |
| cost | **1331.63 s**, peak RSS **20.04 GiB** | 251: 1329.8 s / 17.30 GiB; 249: 1345 s / 20.00 GiB |

The h5ad was re-hashed and deleted as a verified byte duplicate of `submissions/222-h1union-g028.h5ad`; the meta is kept at `submissions/254-ka-h1union-g028.h5ad.meta.json` and `work/build/`, and is the same-builder reference for §5's second diff. The `.vcc` was not re-prepped: it is a deterministic transform of the h5ad (240 row 1, 235 row 2, 249 §2.2 all took `8df4dadd…` to `f13ce40b…`), so the byte identity rests on the h5ad row above, as 247's, 248's and 251's did.

---

## 3. Known-answer rows 2 and 3 — 246 through the builder on 162's PLOLO (`repro254.py`, one fold × mode per process; `work/repro254-*.{json,log}`, `repro254-summary.log`)

162's PLOLO with the **production builder's code substituted in** (222/249's pattern): `extract_cd4_partial`
for `lib246.cd4_partial`, the builder's `freeze_complete_mean3` for `plolo.cd4_for`, and the WHOLE of
`apply_frozen_cd4_residual` ([deflate →] size → merge → one preserve) for 246's `compose`; the builder
imported with the arm's env (`VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_PARTIAL_G=0.34
[VCC_CD4_DEFLATE_K=5]`, one process per configuration). Fold caches sha-identical to 246's. Tolerance 1e−9 on
`cos_agg`, `a23_cos`, `pds`. **104 known-answer rows over 4 folds × 2 modes, every one at 0.0e+00**; the
leave-top-k-out statistic itself re-derives 246's `analyse246.json` retention for E (0.615372) and C2
(−0.066308) from 246's own per-target arrays at 0.0e+00.

### 3.1 Row 2 — `VCC_CD4_PARTIAL_G=0.34` alone == 246 `C2_partial@0.34`

| fold | partial rows (1 + 2) | fraction (= 246's) | KA-B partial == lib246, frozen == cd4_for | KA-C arm vs 246 | effect vs 246 compose | macro / agg / Δ`pds` |
|---|---|---|---|---|---|---|
| `k562_gwps` | 13 (4 + 9) | 0.130575324940 | bit-identical / bit-identical | **0.0e+00** | 0.0 | +0.588 % / +0.650 % / +0.001096 |
| `orion_hct116` | 15 (4 + 11) | 0.130276983163 | bit-identical | **0.0e+00** | 0.0 | +0.108 % / +1.120 % / +0.000505 |
| `orion_hek293t` | 15 (4 + 11) | 0.131495916103 | bit-identical | **0.0e+00** | 0.0 | +1.064 % / +1.699 % / +0.003730 |
| `kolf2` (apart) | 13 (3 + 10) | 0.110076158101 | bit-identical | **0.0e+00** | 0.0 | +0.017 % / +1.070 % / +0.000332 |

Three cancer folds: **rel macro +0.005866912 (246 +0.005866912), aggregate +0.011563922 (246 same), Δ`pds`
+0.001777069 (246 same), \|Δ\| 0.0e+00 each; `pds` positive 3/3; leave-top-3-out −6.6 % (246 −6.6 %)**.
CD4's fraction and the union's fraction are the stack's on every fold (\|Δ\| 0.0).

### 3.2 Row 3 — `VCC_CD4_DEFLATE_K=5` + `VCC_CD4_PARTIAL_G=0.34` and 246's `E_cd4_svd5+partial@0.34`

**The finding the mandate did not anticipate.** 246 built its E arm as `[(deflated cd4, st["frac_cd4"]),
(union, frac_x), (partial, fr)]` — the deflated residual at the STACK's fraction (its `B2_cd4_svd5`, macro
−1.286 %), plus the partial rows (+0.587 %) = −0.699 % — while the builder's `K=5` path re-solves the fraction
on the deflated object (249 §1; 246's *resolved* variant `B2_cd4_svd5_resolved`, macro −0.809 %), which is the
shipped `cd4svd5`. So "`K=5` + `G=0.34` reproduces E" cannot hold literally: the builder's composition is the
resolved one. Both were measured:

| construction | macro | aggregate | Δ`pds` (3 cancer) | per fold Δ`pds` | LT3O `pds` | vs 246 |
|---|---:|---:|---:|---|---:|---|
| **(a) 246's E, through the builder** — `solve_cd4_dose` pinned to size the UN-deflated residual (= the stack's fixed fraction 0.373761 / 0.370782 / 0.376056), everything else the builder's | **−0.698529 %** | −0.023163 % | **+0.014045983** | +0.010454 / +0.009991 / +0.021693 | **61.54 %** | **0.0e+00** on every fold and on the effect matrices; 246 `E_cd4_svd5+partial@0.34` |
| `K=5` alone (249's arm) through this builder | −0.808675 % | −0.461748 % | +0.012262757 | +0.009987 / +0.009446 / +0.017355 | 56.0 % | 0.0e+00, 246 `B2_cd4_svd5_resolved` (249 KA 2 re-verified) |
| **(b) the builder's real path** — resolved deflation (fraction 0.395879 / 0.394269 / 0.397869) + union + partial; **the package's object** | **−0.222188 %** | **+0.701628 %** | **+0.014045525** | +0.011092 / +0.009957 / +0.021088 | **61.7 %** | effect matrix == lib246 `compose` at the re-solved fraction, max\|Δ\| 0.0 (no 246 row exists) |
| **the increment (b) − `cd4svd5`** — the slot's one variable | **+0.585789 %** (+0.600 / +0.110 / +1.050) | +1.162 % | **+0.001782768** (3/3: +0.001105 / +0.000511 / +0.003733) | | **−6.3 %** | = the partial row (row 2) to 4 decimals: additive |

`kolf2` apart: (a) −2.218 % / +0.004107; (b) −1.888 % / +0.004471; increment +0.018 % / +0.000332. Δq_w: (b)
+0.00046 vs the stack, increment −0.00017. The partial's solved fraction is identical with and without the
deflation on every fold (it reads only the base and the residual); the deflation audit is identical with and
without the partial rows.

**Why (b) is the right object and (a) the right known answer.** The mandate says the composed arm "must run
deflation on the COMPLETE rows exactly as `cd4svd5` does, then add the partial-row residual" — `cd4svd5` is
the resolved variant, so (b) is the mandated composition and the package carries it. 246's E numbers are the
fixed-fraction composition; (a) proves the builder's extract → size → merge → preserve is 246's arithmetic to
the bit. The two differ by 0.48 % macro cosine (the re-solve buys cosine at equal `pds`, 249 §1) and by
5e−7 in `pds`; neither differs in the increment, which is the quantity a slot on top of `cd4svd5` measures.

---

## 4. The price (`price254.py`; 246 §4's route = 218 §7.4's rate card on the macro cosine, `pds` from PLOLO at 0.7× / 1.0× / 2.0×; KA: 246's own three price rows reproduce at 0.0e+00)

| arm | macro | Δ`pds` (σ) | `mse` | `reach` | `nmae` rc / exact | `fid` | `jac` | `pds` 0.7× / 1.0× / 2.0× | **TOTAL 0.7× / 1.0× / 2.0×** | × +0.0019 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `E_cd4_svd5+partial@0.34` (246's construction; KA) | −0.699 % | +0.014046 (9.4) | −0.00023 | −0.00015 | −0.00006 / +0.00017 | −0.00005 | −0.00001 | +0.00367 / +0.00524 / +0.01047 | **+0.00316 / +0.00473 / +0.00997** (= 246 §4) | 1.66 / 2.49 / 5.25 |
| `cd4svd5` = `B2_cd4_svd5_resolved` (KA) | −0.809 % | +0.012263 (8.2) | −0.00027 | −0.00018 | −0.00007 / +0.00031 | −0.00006 | −0.00001 | +0.00320 / +0.00457 / +0.00914 | **+0.00262 / +0.00399 / +0.00856** (= 246 §4) | 1.38 / 2.10 / 4.50 |
| **COMPOSED resolved (the package) vs `h1union`** | −0.222 % | +0.014046 (9.4) | −0.00007 | −0.00005 | −0.00002 | −0.00002 | −0.00000 | +0.00367 / +0.00524 / +0.01047 | **+0.00350 / +0.00508 / +0.01031** | 1.84 / 2.67 / 5.43 |
| **INCREMENT vs `cd4svd5` (the slot's variable)** | +0.586 % | +0.001783 (1.2) | +0.00019 | +0.00013 | +0.00005 / +0.00005 (`e` 0.75) | +0.00004 | +0.00001 | +0.00047 / +0.00066 / +0.00133 | **+0.00089 / +0.00109 / +0.00175** | **0.47 / 0.57 / 0.92** |

The increment is under the shape-only gate at every transfer factor, its `pds` is 1.2 PLOLO σ, its
leave-top-3-out is negative (top-1 carries half), and it sits inside the 0.0045 ordering noise on every
reading — exactly 246 §3.7's "NOT build-worthy alone". The composed arm as a whole clears the gate at every
factor because `cd4svd5` does. Board raw `pds` implied for the increment: **+0.0012 / +0.0018 / +0.0036** on
top of `cd4svd5`'s realised value. The package-level `mse` leg is in §8.

---

## 5. Falsifier 1 — the rule-10 audit diff (`metadiff254.py`, sets pre-declared in the source before any diff was read)

**Against the mandate's incumbent `submissions/249-cd4svd5-f034.h5ad.meta.json`** (`--against cd4svd5`,
`work/metadiff254-cd4svd5.log`; the `output.*` bucket is reported, not asserted). Every non-`output` difference,
verbatim:

```
[A-config            ] ADDED   config.cd4_partial_g = 0.34
[B-effects/operator  ] ADDED   effects.cd4_partial.*  (33 keys: the row census, the frame, the merged preserve, panel_residual_nonzero_entries 0)
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm     4.5227571142920135 -> 4.529688365088246
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.post_add_global_preserve_scale         0.9115882187626339 -> 0.9101933222360135
[B-effects/operator  ] CHANGED effects.final_sha256                                                a28d5be2e9cb3ed1... -> 885f6d8c8f865b3d...
[B-effects/operator  ] CHANGED effects.h1_union.merged_preserve.candidate_pre_preserve_median_norm 4.5227571142920135 -> 4.529688365088246
[B-effects/operator  ] CHANGED effects.h1_union.merged_preserve.post_add_global_preserve_scale     0.9115882187626339 -> 0.9101933222360135
[D-implementation    ] CHANGED hashes.implementation.scripts/make_submission_pm_cd4.py             2d63b979... -> 3615d20f...

A. CONFIG added ['cd4_partial_g'] expected ['cd4_partial_g']; changed [] removed []
B. AUDIT changed == final_sha256 (must) + the preserve pair (may, MOVED); added 33 keys, all under effects.cd4_partial.; unexpected []
   the preserve scalar is ONE number written in all three blocks; every pre-declared partial value OK (rows 15 / 4 / 11, carrying 15,
   target_g 0.34, coverage 0.05, shipped_median_rule_would_solve 0.0, panel entries 0, base norm 4.229295604454299);
   solved_fraction == the staged row's; effects.final_sha256 == the staged composed array 885f6d8c...
C. FROZEN INPUTS differences: 0      D. IMPLEMENTATION: this file only      F. UNCLASSIFIED: []
PINNED scalars (128) -- each byte-identical: 128/128, of which 58/58 `effects.h1_union.*`
FALSIFIER 1 (cd4svd5) PASS
```

Pinned and byte-identical to `cd4svd5`'s: `config.cd4_deflate_k` **5** · `config.alpha` 0.982 ·
`config.cd4_target_f` 0.34 · `config.h1_union_g` 0.28 · `config.h1_union_sources` · `operator.sha256`
`e2f65ecf…` · `operator.axis_sha256` `edf6ba88…` · `effects.{baseline_pre,candidate_pre,post_whiten}_sha256` ·
all seven `effects.whitening.*` · `effects.cd4.fraction` **0.34104860341941284** · the whole CD4 frame ·
**every `effects.cd4.deflation.*` key** (k 5, n_carried 291, energy removed 0.19866046506100588, top-k
0.1982359547102955, `U_5` `f06e0262…`) · every `effects.h1_union.*` field (union fraction
0.19207161964444322, 77 carrying targets, source shas) except the preserve pair · `frozen_bundle.sha256` ·
`output.n_rows` 360,000 · the emitter's `rc_callset` triple and every other `config.*` key.

`def_postbuild.py` on the arm (`--allow-config-diff cd4_deflate_k,cd4_partial_g`, incumbent = the shipped 222
meta as the script requires): **32/32**, `moved=['cd4_deflate_k', 'cd4_partial_g']`.

**Against the same-builder unflagged rebuild** `submissions/254-ka-h1union-g028.h5ad.meta.json` (`--against ka`, `work/metadiff254-ka.log`): config ADDED exactly `['cd4_deflate_k', 'cd4_partial_g']`, changed none; audit CHANGED exactly 249's eleven pre-declared keys (the CD4 frame re-solved on the deflated object, `effects.final_sha256`, the preserve pair 0.9093569351340488 → 0.9101933222360135); ADDED 38 keys, all under `effects.cd4.deflation.` (249's five) and `effects.cd4_partial.` (33); frozen inputs 0; **implementation differences NONE** (same builder bytes); unclassified 0; **116/116 pinned** (58/58 `effects.h1_union.*`). FALSIFIER 1 (ka) PASS.

---

## 6. Falsifiers 2 and 3 — the profile moved on the 15 partial rows and nowhere else

### 6a. On the effect matrix — exact (`effdiff254.py`; the staged matrices, same builder)

| comparison | moving rows | n | rel displacement min / median / max | cosine min / max | NOT moving: n, `‖E_arm − s·E_ref‖/‖E_ref‖` max | separation |
|---|---|---:|---|---|---|---:|
| **composed vs `cd4svd5`** (`k5`) | the 15 partial rows | 15 | 0.112 / 0.278 / 0.490 | 0.8968 / 0.9937 | 285 rows, **1.1e−16** (norm ratio = s = 0.9984698172947937 to the last digit) | **1.0e+15×** |
| `partial` vs `h1union` (`off`) | the 15 partial rows | 15 | 0.112 / 0.278 / 0.490 | 0.8968 / 0.9937 | 285 rows, 1.0e−16 (s = 0.9986354827609096) | 1.1e+15× |
| composed vs `h1union` | 276 deflated + 15 partial | 291 | 0.018 / 0.102 / 0.736 | 0.6786 / 0.9999 | 9 rows, 1.1e−16 (s = 1.0009197566650123) | 1.6e+14× |

Composed vs `cd4svd5`: ‖ΔE‖/‖E‖ **0.070779**, 0/300 rows at zero displacement, the 15 partial rows carry
99.956 % of ‖ΔE‖² (the rest is the 0.153 % preserve scalar). Most displaced: STK3 0.490, GMEB2 0.458, EAPP
0.445, TAF4 0.423, PLAGL2 0.405. CD4's fraction, frame and deflation audit identical between `k5` and
`composed`; the union's fraction identical across `off` / `k5` / `composed`.

### 6b. On the delivered counts (`profile254.py` = 249's `profile249.py` with the carrying set = the 15 partial targets; the arm vs `submissions/249-cd4svd5-f034.h5ad`)

900 (context, target) blocks, 45 carrying (15 targets × 3), 855 not:

| block class | n | digest identical | cosine(arm, cd4svd5) min / median / max | norm ratio min / max | `‖ΔE − s·E‖/‖E‖` min / median / max |
|---|---:|---:|---|---|---|
| **carrying** (the 15 partial targets) | 45 | 0 | 0.88049759 / 0.95920956 / 0.99214462 | 1.00589028 / 1.17852997 | **0.125834** / 0.295324 / 0.535245 |
| **NOT carrying** | 855 | 0 | 0.99994229 / 0.99997338 / 0.99999297 | 0.99797760 / 0.99892332 | 0.003744 / 0.007288 / **0.010731** |

Shared delivered scale read off the non-carrying blocks **0.99847129** (effect level 0.99846982). Separation **CLEAN**: the smallest carrying displacement is **11.7×** the largest non-carrying one (249: 1.9×, 222: 3.3× — here the moving rows are few and large and the shared scalar is 0.15 %, so the emitter's integer re-draw on the 855 non-carrying blocks, 0.4–1.1 %, is well under the 12.6 % floor of the carrying ones). The exact statement is §6a's, at 1e−16.

**Digest counts, as asked: 0/900 identical (0/45 carrying, 0/855 non-carrying).** Non-vacuity: ‖ΔE‖/‖E‖ **0.077454** (per context A 0.078895, B 0.076614, C 0.076789), per block min 0.004012 / p01 0.005103 / median 0.007491 / max 0.535259, **0/900 at zero**, 0/900 < 1e−6. 210's stricter 1 % per-block floor (the script's own "FALSIFIER 4" print) holds on all 45 carrying blocks (min 12.6 % scale-corrected) and, as in 249 §7, cannot hold on the 855 non-carrying blocks whose displacement is by construction the 0.15 % scalar plus the emitter's integer re-draw; the mandate's criterion (zero-block count = 0) is met. **PASS.**

For the record: median delivered |row| ratio 1.00122, mean 1.00205, `A = mean|row|²` ratio 1.00440; `rc_callset` head exposure: top-400 |delivered effect| overlap with `cd4svd5` min 209 / median 398 / max 400 of 400 (the minimum on a partial target), implied call-set jaccard lower bound min 0.8631 / median 0.9985; |Δ log2FC| over the DE universe max 0.637, mean 0.000991, signed block-mean median +0.000007.

---

## 7. Falsifier 4 — challenge `n_pred` (`probe254.py` = 249's `probe249.py`; Wilcoxon + BH against the real `data/challenge/context_{A,B,C}.h5ad` controls, 40 perturbations per context)

| package | A | B | C | pooled median | ρ(−log p, \|H\|) in the call set, A / B / C (a15's key for both, cross-key by design) | call set UP |
|---|---:|---:|---:|---:|---|---|
| `cd4svd5` | 2599 [2599, 2635] | 2599 [2599, 2638] | 2599 [2599, 2660] | 2599 | +0.9309 / +0.9270 / +0.9239 (mean +0.9273) | 49.7 / 49.8 / 48.0 % |
| **`composed`** | 2599 [2599, 2635] | 2599 [2599, 2637] | 2599 [2599, 2660] | 2599 | +0.9274 / +0.9262 / +0.9227 (mean +0.9254) | 49.7 / 49.8 / 47.9 % |

**+0.00 %** against ±2 %. The `rc_callset` prescription pins `n_pred` by construction and it held. **PASS.** (Peak RSS 13.27 GiB, 40.91 s — 222's instrument at its known cost, run alone after the last build.)

---

## 8. Falsifier 5 — delivered amplitude off the package, and the predicted board raw `mse` (`amp254.py`)

Per-block panel totals streamed off all three packages (214's `block_totals` via 243's `totals243.py`), 211's per-context terms via 243's `terms`. **Known-answer rows:** `h1union`'s re-streamed (`A_c`, `B_c`, `cap_c`) equal 243's `invert243.json` to **0.0e+00** on all three contexts and 197's identity at those terms reproduces its board raw `mse` **0.8172630** to 0.0e+00; `cd4svd5`'s re-streamed per-block totals equal 249's `totals-cd4svd5.npz` **bit for bit** and its amplitude vs `h1union` re-reads 249 §9's −0.409 % (-0.4087 %).

| context | `A_c` arm | `A_c` cd4svd5 | ratio | `A_c` h1union | ratio | `B_c` arm / cd4svd5 | cap | leg |
|---|---:|---:|---:|---:|---:|---|---:|---|
| A | 0.232747 | 0.231633 | ×1.00481 | 0.233629 | ×0.99623 | 0.230813 / 0.229698 | 0.233754 | B |
| B | 0.221919 | 0.220958 | ×1.00435 | 0.222708 | ×0.99646 | 0.219973 / 0.219014 | 0.239678 | B |
| C | 0.232048 | 0.231016 | ×1.00446 | 0.232897 | ×0.99635 | 0.230057 / 0.229028 | 0.266129 | B |

**Delivered amplitude `mean_c √(A_arm/A_cd4svd5)` − 1 = +0.2268 %** (`B` term +0.2285 %); bound ±0.5 % vs `cd4svd5`: **PASS**, at 45 % of the bound. Vs `h1union` -0.1829 % — the partial rows put back about half of the amplitude the deflation removed (−0.409 %). All three contexts stay on the `B` leg. The mechanism is 246 §3.7's: the 15 rows carried no residual before, so on them the delivered amplitude rises (norm ratio up to 1.18 on a partial block, §6b) while the preserve scalar takes 0.15 % off everything else.

Predicted board raw `mse` from 197's identity **at the arm's own terms** (243 §5.2: never the floor), with `cos_arm = cos_h1union × (1 + k_cos × rel_agg)`, `rel_agg` = +0.7016 % (composed vs the shipped stack, 3 cancer folds; `cd4svd5` -0.4617 %, the increment +1.1624 %) and `cos_h1union` = 0.107332 (243's measured):

| k_cos | arm cos | **predicted raw `mse`, arm** | Δ vs `h1union` (0.8172630) | predicted `cd4svd5` (249's table) | **arm − `cd4svd5`** | `mse` leg vs `cd4svd5`, overall |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 (amplitude/budget leg alone) | 0.107332 | **0.81761** | +0.00035 | 0.81809 | **-0.00048** | +0.00008 |
| 0.94 (new-source class) | 0.108040 | **0.81619** | -0.00107 | 0.81902 | **-0.00282** | +0.00050 |
| 1.00 | 0.108085 | **0.81610** | -0.00116 | 0.81908 | **-0.00297** | +0.00052 |
| 1.36 (re-dose class) | 0.108357 | **0.81556** | -0.00170 | 0.81943 | **-0.00387** | +0.00068 |

The arm's aggregate cosine is ABOVE `h1union`'s on PLOLO (+0.70 %), so the identity predicts raw `mse` falling below both incumbents: −0.0011 … −0.0017 vs `h1union`, −0.0028 … −0.0039 vs `cd4svd5`'s 249 prediction — a `mse` leg of **+0.0005 … +0.0007 overall vs `cd4svd5`**, larger than the rate card's +0.00019 in §4 (the card prices the MACRO cosine, +0.586 %; the identity reads the aggregate, +1.162 %, which the 15 heavy rows dominate). The preregistration's falsifier 2 fires above `cd4svd5`'s realised raw `mse` + 0.0010.

---

## 9. Falsifier 6 — `vcc prep --dry-run`, nnz, max counts per cell

```
✓ dry run — no .vcc written
  input:  submissions/254-cd4svd5-partial-f034.h5ad   cells: 360000   genes: 18533 (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'   contexts: column 'context' → A=120000, B=120000, C=120000
  targets: verified against the official list   normalization: counts-preserved
```

| cap | limit | measured | margin |
|---|---:|---:|---:|
| stored entries in `.X` | 4,750,000,000 | **1,577,695,797** (`cd4svd5` 1,577,615,064, +0.005 %) | **33.21 %** of cap |
| counts per cell | 1,000,000 | **22,323** (the exact-total emitter) | 44.8× under |
| rows | 360,000 | 360,000 | exact |
| cells per (context, target) block | 400 | 400 on all 900 | exact |

Dry run 14.7 s / 20.45 GiB; the real prep 23.85 s / 24.51 GiB. **PASS.**

---

## 10. The package

| artifact | bytes | sha256 |
|---|---:|---|
| **`submissions/254-cd4svd5-partial-f034.vcc`** | **2,824,908,800** | **`eea5ed1802a2337dc3d69b38bec4e7962835ba15b194c650015049b7ff035bf0`** |
| **`submissions/254-cd4svd5-partial-f034.h5ad`** | **12,648,807,800** | **`a9c3f59a37a678a1d481f924df83b6c68e3f7d8bf8ebe073201bb18a2f07bc8b`** |
| `submissions/254-cd4svd5-partial-f034.h5ad.meta.json` | 27,655 | (`output.sha256` = the h5ad's; `effects.final_sha256` `885f6d8c…` = the staged `composed` matrix; builder `3615d20f…` = commit `32c5219`) |
| `submissions/254-ka-h1union-g028.h5ad.meta.json` (KA row 1; its h5ad, a proven byte duplicate of `222-h1union-g028.h5ad`, was deleted) | 25,335 | `output.sha256` `8df4dadd…` |

Configuration (two flags on top of the shipped `h1union` configuration, one of them `cd4svd5`'s):
`VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 VCC_CD4_DEFLATE_K=5 VCC_CD4_PARTIAL_G=0.34`,
`VCC_UNION_SOURCES` unset (the H1 pair), `--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25
--overdispersion 4.0 --whitener-id union_pruned --emitter rc_callset`, operator `e2f65ecf…`. **No model name
chosen, no description written, nothing submitted, nothing armed.** §14 applies at submit time.

---

## 11. The preregistration draft (`PREREGISTRATION-DRAFT.md`, for the lead to commit before any slot — NOT committed by this run)

Six members vs the `cd4svd5` board point (246 §4's rate card on the increment; `pds` at three transfer
factors): `pds` +0.00047 / +0.00066 / +0.00133 · `mse` +0.00019 · `reach` +0.00013 · `nmae` +0.00005 (exact
identity +0.00005, `e` 0.75, a magnitude bound) · `fid` +0.00004 · `jac` +0.00001 · **overall +0.00089 /
+0.00109 / +0.00175**; board raw `pds` +0.0012 … +0.0036 on top of `cd4svd5`'s. Five falsifiers: (1) board raw
`pds` up by less than +0.0006 over `cd4svd5`'s; (2) raw `mse` worse than `cd4svd5`'s by more than +0.0010;
(3) overall Δ vs `cd4svd5` outside [−0.0010, +0.0040]; (4) the rule-10 diff not exactly `cd4_partial_g` + its
consequences (measured here: PASS); (5) delivered amplitude vs `cd4svd5` beyond 0.5 % (measured here:
+0.227 % vs cd4svd5). **The draft says in its own §2 that the increment is NOT a cleared arm** (0.47× … 0.92× the
gate, inside the ordering noise); what a slot buys is stated, and the decision is the lead's after the
`cd4svd5` point.

---

## 12. October

The driver (`scripts/def_pipeline.py`, untouched) passes `VCC_CD4_PARTIAL_G` through to every builder
subprocess when exported but neither prints nor asserts it — the same visibility gap 251 §7 records for
`VCC_CD4_DEFLATE_K` / `VCC_CD4_DEFLATE_BASIS` / `VCC_UNION2_*`. The flag composes with 247's carrying-row CD4
sizing and with 251's pinned basis (`envcheck254`: `V5_abc` + `G=0.34` imports; the partial rows are never
deflated, so the pinned basis question does not arise for them). On a disjoint panel the partial-row count is
a census quantity (`effects.cd4_partial.rows`, printed by the builder as `[cd4-partial] …`): the carrying-row
rule sizes it at any count ≥ 1, and a panel with **zero** partial rows makes the builder refuse loudly
(`solve_carrying_row_fraction` raises) rather than ship a silent no-op. No runbook row is proposed: the arm is
not cleared.

---

## 13. Process notes, on the record

- 253 rewrote `docs/READ-FIRST.md` (v2, 17:00Z) while this run was underway; read after the fact, nothing here
  depends on it. 252 and 253 were live throughout; neither touches the builder; `other_builders()` read 0 at
  every gate.
- After the arm build the prep gate held at free+inactive 25–27 GiB for two minutes (the fresh 12.6 GB h5ad's
  clean page cache, 247 §3.4's / 251 §5's symptom); as 251 did, the clean cache was dropped with `purge`
  (→ 41.0 GiB, `memory_pressure` 81 % free) and the gate opened on the next poll. No process was touched; the
  two armed submit scripts slept through everything.
- `git status` showed only this experiment's files and the builder modified; the builder was committed with
  its path staged explicitly (never `git add -A`, agents running).

---

## 14. Cost and memory (every ≥ 15 GiB stage alone; the gate logged before each)

| stage | wall | peak RSS |
|---|---:|---:|
| `unit254.py` / `envcheck254.py` / `suite254.py` / 235's tests | 1 s / 12 s / 58 s / ~2 min | 0.20 / — / 0.31 / 1.06 GiB |
| `repro254.py` × 8 (4 folds × 2 modes) | 5–7 s each after staging | 4.63–5.13 GiB |
| staging `pre` / `off` / `k5` / `partial` / `composed` (each alone) | 11–18 s each | 19.95 / 19.95 / 17.39 / 19.99 / 20.15 GiB |
| `partialcheck254.py` / `effdiff254.py` / `price254.py` / `metadiff254.py` | seconds | 0.44 / < 0.5 / — / — GiB |
| ARM build (`build254.sh`, both flags) | **1370.7 s** | **20.25 GiB** |
| `def_postbuild.py` (arm / KA) | ~1 s | 0.12 GiB |
| `vcc prep --dry-run` / `vcc prep` | 14.7 s / 23.85 s | 20.45 / 24.51 GiB |
| KA build (`build254.sh`, both flags unset) | **1331.63 s** | **20.04 GiB** |
| `profile254.py` / `amp254.py` / `probe254.py` | 24.15 s / 40 s / 40.91 s | 3.61 / 1.18 / 13.27 GiB |

Disk: the arm's `.h5ad` + `.vcc` 15.5 GB under `submissions/`; the KA h5ad (12.6 GB) deleted after verification;
`work/` holds five 44 MB effect matrices, the CD4 / partial / union residual npz, 131 MB of PLOLO fold-cache
copies, per-target arrays and three totals npz (all gitignored).

---

## 15. Files

| file | what it does |
|---|---|
| `patch254.py` | the ONE builder change as eleven exact-anchor replacements; `--check <HEAD copy>` reproduces the working tree |
| `intended254.py` | HEAD + `patch254.py` → `work/make_submission_pm_cd4.254.py` and `work/builder.intended.sha256` (`3615d20f…`), the gate's reference |
| `envcheck254.py` | the import-time rules, eleven subprocess imports incl. the pre-254 byte copy → `work/envcheck254.json` |
| `unit254.py` | §1's 44 unit rows: OFF bit-identical, ON == lib246's arithmetic and 246's `merged_many` chain, refusals → `work/unit254.json` |
| `suite254.py` | the builder's 57 tests under both builders → `work/suite254.{json,log}`; 235's tests → `work/tests235.log` |
| `repro254.py` | rows 2 and 3 through 162's PLOLO, one fold × mode per process → `work/repro254-<mode>-<fold>.{json,log}`, `pertarget254-*.npz`; `--summary` → `work/repro254-summary.{json,log}` |
| `stage254.py` / `stage254.sh` | row 1 at staging level + the arm's effect matrix and residuals → `work/stage/` |
| `partialcheck254.py` | the production-panel partial residual == `lib246.cd4_partial`, the frozen object == `plolo.cd4_for` → `work/partialcheck254.json` |
| `effdiff254.py` | falsifiers 2/3 at the effect level, three comparisons → `work/effdiff254.{json,log}` |
| `price254.py` | §4 (246's route; KA against `price246.json`) → `work/price254.{json,log}` |
| `build254.sh` / `run254.sh` | 222's build flags + the two passthroughs; the gated serial sequence (arm build, postbuild, prep, KA build) → `work/build/` |
| `metadiff254.py` | falsifier 1 (`--against cd4svd5` / `--against ka`) → `work/metadiff254-{cd4svd5,ka}.{json,log}` |
| `profile254.py` / `probe254.py` / `amp254.py` | falsifiers 2/3 on the delivered counts, 4, 5 and 6's caps → `work/profile254.json`, `probe254.{json,log}`, `amp254.{json,log}`, `totals-*.npz` |
| `PREREGISTRATION-DRAFT.md` | §11's draft, for the lead |
| `work/make_submission_pm_cd4.PRE.py` | the git control, HEAD `10d0177e…` (commit `c321f76`), taken before the first edit |
| `work/commit-message.txt` | the builder commit's message with the known-answer hashes |

```bash
cd ~/Code/virtual-cell-challenge-2026
git show HEAD:scripts/make_submission_pm_cd4.py > experiments/254-partial-rows/work/make_submission_pm_cd4.PRE.py
uv run python experiments/254-partial-rows/patch254.py; uv run python experiments/254-partial-rows/intended254.py
uv run python experiments/254-partial-rows/envcheck254.py; uv run python experiments/254-partial-rows/unit254.py
uv run python experiments/254-partial-rows/suite254.py; uv run python experiments/235-october-hardening/tests/test_october_hardening.py
cp experiments/246-direction-hypotheses/work/plolo/*.np* experiments/254-partial-rows/work/plolo/
for m in partial composed; do for f in k562_gwps orion_hct116 orion_hek293t kolf2; do uv run python experiments/254-partial-rows/repro254.py --mode $m --fold $f; done; done
uv run python experiments/254-partial-rows/repro254.py --summary; uv run python experiments/254-partial-rows/price254.py
bash experiments/254-partial-rows/stage254.sh                       # pre / off / k5 / partial / composed, gated
uv run python experiments/254-partial-rows/partialcheck254.py; uv run python experiments/254-partial-rows/effdiff254.py
bash experiments/254-partial-rows/run254.sh                         # arm build, postbuild, prep, KA build, gated
uv run python experiments/254-partial-rows/metadiff254.py submissions/254-cd4svd5-partial-f034.h5ad.meta.json submissions/249-cd4svd5-f034.h5ad.meta.json --against cd4svd5
uv run python experiments/254-partial-rows/metadiff254.py submissions/254-cd4svd5-partial-f034.h5ad.meta.json submissions/254-ka-h1union-g028.h5ad.meta.json --against ka
uv run python experiments/254-partial-rows/profile254.py submissions/254-cd4svd5-partial-f034.h5ad submissions/249-cd4svd5-f034.h5ad
uv run python experiments/254-partial-rows/amp254.py
uv run python experiments/254-partial-rows/probe254.py cd4svd5=submissions/249-cd4svd5-f034.h5ad composed=submissions/254-cd4svd5-partial-f034.h5ad --perts=40
```

---

## 16. Honest limits

- **No board member was measured and no `cell-eval2` was run.** Rows 2 and 3 reproduce 246's PLOLO instrument
  through the builder; they are not new evidence that the partial rows' `pds` transfers.
- **The mandate's row 3 digits belong to a construction the builder does not carry** (§3.2). 246's E arm sits
  on the fixed-fraction deflation; the builder's `K=5` is the resolved one. Both were measured; the package
  carries the resolved composition, and the increment over `cd4svd5` is unaffected (additive to 5e−7 in
  `pds`, 0.001 % in macro).
- **The increment is under the gate on every reading** (0.47× … 0.92× at `pds` transfer 0.7 … 2.0×), 1.2 PLOLO
  σ, top-1-concentrated (leave-top-3-out −6.3 %), and inside the 0.0045 ordering noise. The package exists
  because the mandate asked for it before the `cd4svd5` decision; nothing here argues it should ride a slot,
  and the preregistration draft says so.
- **The one-condition rows are heavy**: the carrying-row median norm is 12.79 against the complete rows' 4.52,
  so g = 0.34 delivers 0.34 × the base's typical row norm on rows whose own information is a single condition's
  effect. 246's norm-equalised variant read less `pds` (+0.00132 vs +0.00178) and was not built; g was not
  re-tuned.
- **`h1-test` is oracle-leaked** for anything on top of `h1union` (222 §11) and was not used.
- **`probe254.py` peaks at ~13 GiB**, above the ~10 GB guideline, as 222's and 249's did; it ran alone after
  the last build.
- **The `nmae` sign is unread**: 246's exact identity gives +0.00005 for the partial rows as a magnitude bound.

---

## 17. Verdict

**BUILDABLE.** `submissions/254-cd4svd5-partial-f034.vcc` (`eea5ed18…`, 2,824,908,800 B) is `cd4svd5` plus
one variable (`VCC_CD4_PARTIAL_G=0.34`) on top of the shipped configuration; the builder change is committed
as `32c5219` with the flag unset proven byte-inert at staging (`abeb0749…`; `K=5` alone `a28d5be2…`) and
package (8df4dadd…) level; the partial residual is 246's `cd4_partial` bit for bit on every fold
and on the production panel; 246's `C2_partial@0.34` and `E_cd4_svd5+partial@0.34` reproduce through the
builder at 0.0e+00 (macro +0.587 % / `pds` +0.00178; macro −0.699 % / `pds` +0.01405 / leave-top-3-out
61.5 %); the package's own composition (the resolved deflation + the partial rows) reads macro −0.222 % /
`pds` +0.01405 vs `h1union` and its increment over `cd4svd5` is +0.00089 … +0.00175 overall — **under the
+0.0019 gate**; every falsifier passes with its numbers above (config diff exact; 15 partial rows move, 285
rows one scalar to 1e−16; 0/900 zero blocks; `n_pred` +0.00 %; amplitude +0.227 % vs cd4svd5;
prep clean, 33.21 % of the density cap). Not submitted; the preregistration is a draft for the lead.
