# 222 — `h1union` is **BUILDABLE**. The coverage-aware sizing rule is in, the default path is byte-identical to `cd4f034`, and the production builder reproduces 218 to **0.00e+00**.

**Date:** 2026-09-04 · deterministic · **no `vcc submit`, no `cell-eval2`, no submission.**
Two builds, strictly serial, peak RSS **14.67 GiB** and **9.18 GiB** (`/usr/bin/time -l`), both
under the 22 GB budget. One file changed outside this experiment directory:
`scripts/make_submission_pm_cd4.py`. Nothing under `data/`, `src/`, or any other experiment was
touched.

> ⚠️ **`vcc prep` peaked at 22.33 GiB RSS** (23,976,853,504 B, 22 s). It is a bigger process than the
> build that feeds it and it must never run concurrently with a build on this machine. Recording it
> because nothing in the repo did.

---

## 0. Verdict

```
BUILDABLE.  submissions/222-h1union-g028.vcc
            sha256 f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e
            2,824,785,920 bytes.  All 3 known-answer rows and all 5 falsifiers pass.
```

**Two results the coordinator must read before spending the slot**, both from INBOX-1 and both
adverse to the preregistration's wording:

1. **Exact `nmae` is slightly NEGATIVE, not "up or flat"** — 216's exact 188-identity puts the arm
   at raw `nmae` **+0.00039** over the incumbent, i.e. **−0.000108 overall**. That is far inside
   preregistration falsifier 4's ±0.002 bound, so **the falsifier does not fire** — but the sign of
   the prediction is contradicted. And the same instrument **mis-signs its own family's board
   point** (§8), so read it as a magnitude bound, not as a sign.
2. **The `pds` gain is CONCENTRATED and the addendum's §5 pre-commitment FIRES.** Top-3 targets per
   cancer fold carry **99.8%** of the mean Δ`pds`; leave-top-3-out retains **14.6%**. The *cosine*
   gain is not concentrated (leave-top-3-out retains **72.7%**). Under the pre-committed reading a
   positive board result is **weaker** evidence for 223's stack.

Neither is a stop condition in my mandate. **CD4's solved fraction did not move by one bit**, which
is the stop condition, so the arm is single-variable and the package is ready.

---

## 1. The change

`scripts/make_submission_pm_cd4.py`, one file, behind `VCC_H1_UNION_G` (unset = current behaviour).

| added | what it does |
|---|---|
| `solve_carrying_row_fraction` | **the coverage-aware sizing rule** — `median(‖added_i‖ : ‖added_i‖ > 0)` over the base's all-row median. 218 §7.1's `solve_fraction_carrying`, verbatim. |
| `single_source_residual` | 218's `multi_source_residual` for one bank — `prepare_cd4_mean` generalised through `map_bank_to_axis` (duplicate-symbol safe). `chunk_rows` pinned at 218's **256**. |
| `union_residual` | 218's `union_rows` — per-gene mean over whichever sources carry a coordinate. |
| `prepare_h1_union` | the two H1 banks → the union, panel columns excluded. |
| `apply_frozen_cd4_residual(..., extra=, extra_audit=)` | merges `frac_cd4·cd4 + frac_x·union` into **ONE** add with production's **ONE** global preserve (218 §7.2's `merged()`). |

**`solve_cd4_fraction` is not modified.** CD4 is still sized by the shipped median-over-all-rows
rule. Only the second residual uses the new one.

**The bug, reproduced live inside the production build.** The arm's own meta records
`effects.h1_union.residual_frame.residual_median_row_norm_over_all_rows = **0.0**` and
`shipped_median_rule_would_solve = **0.0**`. The union's carrying-row median is **6.1654** against
the base's **4.2293**, giving `g_per_unit = 1.4578` and a solved fraction of **0.19207** at
`g = 0.28`. Under production's rule that fraction is **exactly zero**. This is 200/207/218's silent
zero, now printed by the builder itself on every run with the flag on.

Union on the challenge panel: **77/300 targets carried, 7 by both sources**, 1,398,013 supported
entries; `gse295214_h1` covers 64, `h1_2025` covers 20.

```bash
experiments/222-h1-union-build/build222.sh submissions/222-ka-cd4f034-f034.h5ad     # flag off
VCC_H1_UNION_G=0.28 experiments/222-h1-union-build/build222.sh submissions/222-h1union-g028.h5ad
```

---

## 2. Known-answer row 1 — the default path is byte-identical to `cd4f034`

Rebuilt with the **changed** builder, `VCC_H1_UNION_G` unset:

| | value |
|---|---|
| `shasum -a 256 submissions/222-ka-cd4f034-f034.h5ad` | **`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`** |
| required | `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a` |
| bytes | **12,648,807,800** (required 12,648,807,800) |
| verdict | **PASS** |

Its `.meta.json` differs from `210-cd4f034.h5ad.meta.json` in **exactly three keys**, all
non-substantive: the builder's own implementation hash (`6cc91501…` → `71cf433d…`), `output.path`,
and `peak_rss_gib` (19.833 → 14.672). `output.sha256`, `output.nnz`, every `effects.*`, every
`operator.*`, every input hash: identical.

## 3. Known-answer row 2 — proved against **git**, not against the working tree

`git show HEAD:scripts/make_submission_pm_cd4.py` was saved to
`work/make_submission_pm_cd4.PRE.py` (sha256 `6cc9150122d5e0689eee64ea89fc831c20c3bb80726be75cee335bcdb49c51a6`, HEAD `ea244a6`) **before the first edit**.
`stage222.py` imports a builder **by path** and stages the production effect matrix without
emitting, so the two builders run in two separate processes on identical inputs:

| tag | builder | `VCC_H1_UNION_G` | `effects.final_sha256` | `effects-*.npy` file sha256 |
|---|---|---|---|---|
| `pre` | `work/make_submission_pm_cd4.PRE.py` (git HEAD `ea244a6`) | unset | `87a4b6d181298f63826dc9caeb9674510a5dbc0d3c388e9e4fbcb359f3cd2de2` | `21396cb3bd7e8884…` |
| `off` | `scripts/make_submission_pm_cd4.py` (changed) | unset | `87a4b6d181298f63826dc9caeb9674510a5dbc0d3c388e9e4fbcb359f3cd2de2` | `21396cb3bd7e8884…` |
| `on028` | `scripts/make_submission_pm_cd4.py` (changed) | **0.28** | `abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315` | — |

**`pre` and `off` are byte-identical.** `off`'s sha equals `210-cd4f034.h5ad.meta.json`'s
`effects.final_sha256`, and `on028`'s equals the arm build's — so the staging harness and the two
real builds compute the same objects.

> ⚠️ **The git control had to come from the saved byte copy, not from `HEAD`.** A coordinator sweep
> commit (`481cd09`) picked up my in-progress edit to `scripts/make_submission_pm_cd4.py` mid-run, so
> `git diff` on that file is now empty and `git show HEAD:` no longer returns the pre-change bytes.
> The copy taken at 04:48 UTC is the control. This is the trap `63400de` records.

**CD4's own solve, `off` → `on028`:** `0.32207201067578145` → `0.32207201067578145`.
**Identical.** The arm is one variable.

## 4. Known-answer row 3 — the production builder reproduces 218 to `0.00e+00`

`repro222.py` runs 162's PLOLO instrument with **the production builder's** `prepare_h1_union` and
`solve_carrying_row_fraction` substituted for 218's own implementations. 4 folds × 6 doses + the
shuffled control + the row split, 75 s, peak RSS 8.27 GiB, **0 failures**.

First, the two implementations produce the **same array**: on every fold,
`max|Δ mean_effect| = 0.000e+00`, `max|Δ votes| = 0`.

| endpoint, 3 cancer folds, `g = 0.28` | 222 measured | 218 recorded | \|Δ\| |
|---|---:|---:|---:|
| relative **MACRO** cosine | **+7.6802815 %** | +7.6802815 % | **0.00e+00** |
| Δ`pds` | **+0.004886738** | +0.004886738 | **0.00e+00** |
| relative aggregate cosine | +10.3065 % | +10.3065 % | 0.00e+00 |

Full dose ladder, three cancer folds (`kolf2` excluded from every headline, as in 218):

| `g` | rel MACRO cos | 218 | Δ`pds` | 218 |
|---:|---:|---:|---:|---:|
| 0.02 | +0.6885 % | +0.6885 % | +0.000587 | +0.000587 |
| 0.05 | +1.6822 % | +1.6822 % | +0.001465 | +0.001465 |
| 0.10 | +3.2316 % | +3.2316 % | +0.002793 | +0.002793 |
| 0.17 | +5.1723 % | +5.1723 % | +0.004028 | +0.004028 |
| **0.28** | **+7.6803 %** | **+7.6803 %** | **+0.004887** | **+0.004887** |
| 0.45 | +10.3980 % | +10.3980 % | +0.004317 | +0.004317 |

**12 of 12 fold × dose cells positive on cosine and `pds`** reproduces (all four folds, doses
0.17/0.28/0.45; per-fold Δ`pds` at 0.28: k562 +0.009516, hct116 +0.002072, hek293t +0.003072,
kolf2 +0.021628).

**Target-shuffled control** (218's falsifier 1), 3 cancer folds: macro **−0.5965 % / −1.4217 % /
−3.1633 %** at `g` 0.17/0.28/0.45, Δ`pds` **−0.002508 / −0.005406 / −0.010357** — strictly negative
on both endpoints at every dose, matching 218 to `0.00e+00`.

**Non-carrying-row bound** (218's falsifier 4, `rowsplit218`'s instrument, `gse295214_h1` alone,
4 folds × 3 doses): largest excursion **0.01778 pp** against 218's reported 0.018 pp — **PASS**,
every `rel_in`/`rel_out` matching 218 to `0.00e+00`. For the **shipped union** the same statistic is
**0.02304 pp**, with carrying rows moving **+11.1 % to +100.3 %**.

---

## 5. Falsifier 1 — the config diff is exactly the pre-declared change

`metadiff222.py`, expected sets declared in the source **before** the diff was read. Full output:
`work/metadiff222.log`. **Every difference is listed, in five buckets.**

### A. CONFIG — 3 additions, 0 changes, 0 removals

```
[A-config            ] ADDED   config.h1_union_g = 0.28
[A-config            ] ADDED   config.h1_union_sizing_rule = 'median over carrying rows'
[A-config            ] ADDED   config.h1_union_sources = ['gse295214_h1', 'h1_2025']
```
pre-declared: `{h1_union_g, h1_union_sizing_rule, h1_union_sources}` — **exact match**.

### B. EFFECTS / OPERATOR — 3 changes, 53 additions all under `effects.h1_union.`

```
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm
                                  cd4f034: 4.479655759355644
                                  h1union: 4.533854576152464
[B-effects/operator  ] CHANGED effects.cd4.post_add_preserve.post_add_global_preserve_scale
                                  cd4f034: 0.9203591354319881
                                  h1union: 0.9093569351340488
[B-effects/operator  ] CHANGED effects.final_sha256
                                  cd4f034: '87a4b6d181298f63826dc9caeb9674510a5dbc0d3c388e9e4fbcb359f3cd2de2'
                                  h1union: 'abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315'
```
pre-declared: exactly those three. **0 unexpected additions, 0 removals.**

The two preserve scalars are the **only** thing in CD4's own block that moves, and they move
because 218's construction is one merged add with one global rescale — `PREREGISTRATION-ADDENDUM.md`
§3 records this in advance as a declared consequence. Ratio **0.909357 / 0.920359 = 0.988046**.
`effects.cd4.post_add_preserve.fraction` is written as CD4's **own** solved fraction (not the
merged add's `1.0`), so a rule-10 diff cannot misread the merge as a dose change.

The new block's substance:

```
effects.h1_union.target_g                                           0.28
effects.h1_union.sizing_rule                    'median over CARRYING rows (coverage-aware)'
effects.h1_union.mode               'union:gene_norm:gse295214_h1+gene_norm:h1_2025'
effects.h1_union.carrying_targets                                     77
effects.h1_union.targets_carried_by_both_sources                       7
effects.h1_union.supported_entries                             1,398,013
effects.h1_union.panel_residual_nonzero_entries                        0
effects.h1_union.merged_add                                         True
effects.h1_union.residual_frame.residual_median_row_norm_over_carrying_rows   6.1654229367116375
effects.h1_union.residual_frame.residual_median_row_norm_over_all_rows        0.0
effects.h1_union.residual_frame.whitened_base_median_row_norm                 4.229295604454299
effects.h1_union.residual_frame.g_per_unit_fraction                           1.4577895501601277
effects.h1_union.residual_frame.solved_fraction                               0.19207161964444322
effects.h1_union.residual_frame.shipped_median_rule_would_solve               0.0      <-- THE BUG
effects.h1_union.residual_frame.coverage_frac                                 0.25666666666666665
effects.h1_union.residual_frame.entries_outside_cd4_support                 740,065
effects.h1_union.residual_frame.targets_outside_cd4_support                       1
effects.h1_union.merged_preserve.post_add_global_preserve_scale               0.9093569351340488
effects.h1_union.sources[0]  gse295214_h1  experiments/115-gse295214-h1/work/compact-donor
     delta.npy   d0fa6555f9a62ce503d2ad55c8bfda954322657b52dcd48a24e29c330e718c7a
     b_ctrl.npy  c88816a57c57c3a088183ae24bca406d8ded374871e3d330aa444fcb13811b4f
     genes.txt   c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723
     perts.txt   7060076272d411d4ee00b898cc26eedcb1eaaac5f5e94e2c60f3b45978320747
effects.h1_union.sources[1]  h1_2025       data/donors/h1_2025
     delta.npy   269628f6eae749857e619bad4556ce9ac8b8810882e7af8a9925e3178eec0cc8
     b_ctrl.npy  bd527b8373d56dab8bfa4d6a1a54560987977e37e22630709fcbbd8f7ce42d7a
     genes.txt   da2b383939d0939e911d273539b136c7a10f0482b2ee509d294e1632d5629822
     perts.txt   977b808e723c73eea5dc4a60cec6e71512387dc520512bf133d3b18fce7287dd
effects.h1_union.source_details[0]  covered_targets 64, mapped_non_panel_genes 18233, duplicates 0,
                                    raw_median_centered_norm 7.084942901352942,
                                    gene_normalized_median_before_rescale 140.2711427079268
effects.h1_union.source_details[1]  covered_targets 20, mapped_non_panel_genes 17777, duplicates 0,
                                    raw_median_centered_norm 5.027205005412602,
                                    gene_normalized_median_before_rescale 130.44886030037927
```

### C. FROZEN INPUTS — **0 differences**. D. IMPLEMENTATION — **1**, the file this arm edits

```
[D-implementation    ] CHANGED hashes.implementation.scripts/make_submission_pm_cd4.py
                                  cd4f034: '6cc9150122d5e0689eee64ea89fc831c20c3bb80726be75cee335bcdb49c51a6'
                                  h1union: '71cf433df2c7075bad4be44975d082a2b4250888175ae435042e06fb8921f43e'
```

### F. UNCLASSIFIED — **empty**. E. build output / runtime — 29 keys, reported not asserted

`output.sha256`, `output.path`, `peak_rss_gib`, `output.nnz` (1,577,729,493 → **1,577,411,397**),
`output.nnz_per_cell`, `output.pct_of_density_cap`, the three `output.contexts.*.stored_entries`,
and 20 `output.rc_callset.*` diagnostics. All are properties of this build, not of the arm.

### The pinned scalars — 55 of 55 byte-identical

`config.alpha` **0.982** · `config.cd4_target_f` **0.34** · `operator.sha256`
**`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`** ·
`operator.metadata_sha256` · `operator.columns_sha256` · `operator.columns` 17,942 ·
`effects.baseline_pre_sha256` · `effects.candidate_pre_sha256` · **`effects.post_whiten_sha256`**
(the whitened base, computed before the add) · all seven `effects.whitening.*` including
`after_all_median_norm` **4.229295604454299** · **`effects.cd4.fraction` 0.32207201067578145** ·
`effects.cd4.residual_frame.{solved_fraction, target_f_requested, f_per_unit_fraction,
residual_median_row_norm_at_unit_fraction, whitened_base_median_row_norm, cd4_fraction_constant,
effective_f_at_cd4_fraction}` · `effects.cd4.{mode, fully_supported_targets 276,
supported_entries 2,712,528, panel_residual_nonzero_entries 0}` ·
`effects.cd4.post_add_preserve.{fraction, baseline_output_median_norm}` · `frozen_bundle.sha256` ·
`config.{emitter, min_cpm, keep_top, keep_spread, tail_spread, overdispersion, spread_map,
effect_scale, emission_seed, rc_key, rc_callset_head, rc_callset_cpm, rc_callset_z_cap,
whitener_id, whiten_gamma, whiten_lambda, cd4_mode, cd4_fraction, cd4_conditions,
cd4_target_policy}` · `output.n_rows` 360,000.

**FALSIFIER 1 — PASS.**

---

## 6. Falsifier 2 — the profile moved on the carrying rows and nowhere else

**The literal form of this falsifier is not achievable for this arm, and the addendum says so in
advance.** 218's construction merges the two residuals into one add followed by production's one
global preserve rescale, so the preserve scalar is a function of both and *every* row is multiplied
by a slightly different constant than `cd4f034`'s (0.988046). Holding CD4's own preserve fixed
instead would let the delivered amplitude float, which is the one thing an α-calibrated build must
not do. So the check is run in the **scale-corrected** form the addendum §3 pre-declares, at two
levels.

### 6a. On the effect matrix — **exact**

`effects-off.npy` vs `effects-on028.npy`, 300 rows, 77 carried:

| row class | n | rel. displacement `‖ΔE_i‖/‖E_i‖` | after dividing out the shared scalar 0.988045753 | cosine(off, on) |
|---|---:|---|---|---|
| **carrying** | 77 | min 0.083599 · median 0.253149 · max 0.493774 | max **4.94e−01**, median 2.53e−01 | 0.8982 – 0.9965 |
| **NOT carrying** | 223 | **exactly 0.011954 for all 223** | max **1.03e−16**, median 9.30e−17 | **1.0000000000** exactly |

`0.011954 = |1 − 0.988046|`. Off the carrying rows the arm's effect matrix is the incumbent's
**times one shared scalar, to machine precision** — 16 orders of magnitude of separation from the
carrying rows.

### 6b. On the delivered counts — the emitter's integer realisation blurs it, but it separates cleanly

`profile222.py` on the two 12.6 GB packages, 900 (context, target) blocks, 231 carrying:

| block class | n | digest identical | cosine(arm, inc) min / median | `‖ΔE − s·E‖/‖E‖` min / median / max |
|---|---:|---:|---|---|
| **carrying** | 231 | 0 | 0.90583963 / 0.97375255 | **0.085295** / 0.235927 / 0.473993 |
| **NOT carrying** | 669 | 0 | 0.99965455 / 0.99983524 | 0.007616 / 0.017937 / **0.025986** |

Shared delivered scale read off the non-carrying blocks: **0.98813355**. Separation **CLEAN** — the
smallest carrying displacement is **3.3×** the largest non-carrying one, and the non-carrying
cosines are 0.9997+ against 0.906–0.997 for the carrying ones. The residual 1.8% on non-carrying
blocks is the quasi/`rc_callset` emitter re-drawing integers under a 1.2% amplitude change, not a
direction change; §6a shows the direction there is identical to 1e−16.

**Digest counts, as asked: 0/900 identical (0/231 carrying, 0/669 non-carrying).**
Nonzero-structure digests: 0/900 identical. Panel-total sums: 5/900 identical.

**FALSIFIER 2 — PASS in the scale-corrected form; the raw digest is 0/900 identical and the reason
is the single shared preserve scalar, which is a declared consequence.**

---

## 7. Falsifiers 3, 4, 5 — non-vacuity, `n_pred`, format

### Falsifier 3 — non-vacuity. **A zero here would be the bug. There is no zero.**

| level | `‖ΔE‖/‖E‖` | per-block min | blocks exactly 0 | blocks < 1e−6 |
|---|---:|---:|---:|---:|
| effect matrix (300 targets) | **0.129030** | **0.011954** | **0 / 300** | 0 / 300 |
| delivered profile (900 blocks) | **0.125080** | **0.014245** | **0 / 900** | 0 / 900 |
| delivered, inside the 10-CPM gate | 0.125134 | 0.014154 | — | — |

Per-block p01 0.016428, median 0.022398, max 0.473538. Per context: A 0.124373, B 0.125403,
C 0.125479. The residual applied, on every block, in all three contexts. **PASS.**

Delivered amplitude, for the record (this is the object prereg falsifier 4 worries about): median
delivered `|row|` ratio **1.00413**, mean **1.00133**, `A = mean|row|²` **1.00301** — a **+0.3 %**
amplitude move, and the *effect*-space median row norm over the operator columns is preserved
exactly by construction.

### Falsifier 4 — challenge `n_pred` within ±2 %

`probe222.py` (210's `probe210.py`, two output filenames changed), Wilcoxon + BH against the real
`data/challenge/context_{A,B,C}.h5ad` control cells, 40 perturbations per context:

| package | A | B | C | pooled median |
|---|---:|---:|---:|---:|
| `cd4f034` | 2,599 | 2,599 | 2,599 | 2,599 |
| **`h1union`** | **2,599** | **2,599** | **2,599** | **2,599** |

**+0.00 %** against a ±2 % bound. **PASS.** (Also: delivered order ρ(−log p, |H|) inside the call
set **+0.9306** against the incumbent's +0.9379; call set 48–50 % UP in both; the `rc_callset`
prescription pins `n_pred` by construction and it held.)

### Falsifier 5 — `vcc prep --dry-run` and the format caps

```
✓ dry run — no .vcc written
  cells: 360000   genes: 18533 (float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts: A=120000, B=120000, C=120000
  targets: verified against the official list
  normalization: counts-preserved
```

| cap | limit | measured | margin |
|---|---:|---:|---:|
| stored entries in `.X` | 4,750,000,000 | **1,577,411,397** | **33.21 %** of cap |
| counts per cell | 1,000,000 | **22,323** (min 17,695) | 44.8× under |
| rows | 360,000 | 360,000 | exact |
| cells per (context, target) block | 400 | 400 on all 900 | exact |

**PASS.** The real package was then written and re-written with `--force`, producing a
**bit-identical** sha256 both times.

---

## 8. INBOX-1 item 1 — the EXACT `nmae`, and why its sign must not be trusted here

`six222.py` = 216's `six216.py` with the ladder cut to three arms and `build_final` replaced by the
production merged construction. 188's exact identity, three pools, `key_precision`, 4 PLOLO folds,
challenge frame A. **Six known-answer values against `216/work/six216.json`, all `0.00e+00`.**

| arm | dlv/a15 | `nmae` wald89 | Δ vs incumbent | overall `nmae` | top-89 precision | Δ prec | `pds` 4f | `pds` 3f |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| a15 `f=0.10` | 1.0000 | 0.98476 | −0.00055 | +0.000152 | 0.58995 | −0.00757 | 0.74214 | 0.78441 |
| `cd4f034` `f=0.34` | 1.0009 | 0.98530 | 0 | 0 | 0.59752 | 0 | 0.74148 | 0.78166 |
| **`h1union` `g=0.28`** | **1.0033** | **0.98569** | **+0.00039** | **−0.000108** | **0.60312** | **+0.00561** | **0.75065** | **0.78655** |

- **`nmae`: raw +0.00039, overall −0.000108.** Preregistration falsifier 4 (`raw nmae rises by more
  than 0.002`) **does not fire** — the measured rise is 5.1× inside the bound. But the
  preregistration predicted "up **or flat**"; on this instrument the member is **slightly down**.
- 🔴 **The instrument mis-signs its own family's board point.** For the same CD4 dose step
  0.10 → 0.34 it gives **−0.000152 overall**; the board delivered **+0.001762** (`210/BOARD-RESULT.md`),
  wrong in sign and 11.6× in magnitude. AGENTS.md §11.1 requires a screen to retrodict its family's
  board point before its numbers are read — **this one fails that retrodiction.** So −0.000108 is a
  bound on magnitude ("this arm cannot move `nmae` much"), not a prediction of sign. On the board's
  own realised rate card (218 §7.4) the same +7.68 % cosine buys `nmae` **+0.00075 overall**.
- **`pds` +0.009166 (4 folds) / +0.004891 (3 cancer folds)** — reproduces 218's +0.004887 through a
  second, independently written scoring path.
- **Ordering precision +0.00561** over the incumbent (which itself gained +0.00757 over a15). The
  sign-accuracy channel that `reach` and `fid` read moves **up**, consistent with 218's `q_w`.
- **Delivered amplitude +0.24 %** over the incumbent. Not a dose increase.

## 9. INBOX-1 item 2 — the concentration check FIRES on `pds`, not on cosine

`conc222.py`, three cancer folds, per-target arrays out of 218's own `score_effects`.

| leave-top-k-out | Δ`pds` retained | top-k mean share of Δ`pds` | Δcos retained |
|---|---:|---:|---:|
| full (no deletion) | +0.004887 (100 %) | — | 100 % |
| **k = 1** | +0.003121 (**63.9 %**) | 42.3 % | **87.3 %** |
| **k = 3** | +0.000716 (**14.6 %**) | **99.8 %** | **72.7 %** |
| k = 5 | −0.000786 (−16.1 %) | 136.4 % | 61.8 % |
| k = 10 | −0.002452 (−50.2 %) | 182.4 % | 43.3 % |

**The addendum §5 pre-commitment fires: the top three carried targets carry 99.8 % of the mean
Δ`pds`, well over half.** Per fold the top-3 shares are 69.7 % / 115.2 % / 114.6 % (over 100 %
because other targets lose rank). Top targets by Δ`pds`: `TIMM50 +0.756`, `C5orf22 +0.646`,
`STT3A +0.382` (k562); `ANKZF1 +0.286`, `TARBP2 +0.227`, `PLAGL2 +0.179` (hct116); `DLAT +0.451`,
`TWF2 +0.345`, `ZNF22 +0.240` (hek293t) — **all nine are carried targets.**

**The cosine gain is a different story and it is not concentrated:** 72.7 % survives deleting the
top three, 43.3 % survives deleting the top ten of ~77 carriers. And **99.8–99.9 % of both Δ`pds`
and Δcos sits on the 26 % of targets the union carries** — there is essentially no spillover, which
is the arm doing exactly what it claims.

The mechanism is worth naming because it changes what the board can tell us: `pds` is a bounded
per-target rank statistic, so a broad direction improvement converts into rank movement only for
the targets sitting near a rank boundary. **A positive board `pds` is therefore weaker evidence for
"a genuinely 4th independent estimate" than the +0.004887 alone suggests**, exactly as
pre-committed. The broad Δcos is the stronger evidence, and it is the channel `mse`/`nmae` read.

## 10. INBOX-1 item 3 — CD4's carrying-row norm ratio, the October risk, measured not changed

`median(‖r_i‖) / q_{0.14}(‖r_i‖)` over CD4's **carrying** rows, per PLOLO fold. **CD4's sizing rule
was not touched.**

| fold | carrying | coverage | median all rows | median carrying | all-rows median sits at | **over-dose now** | **over-dose at the q0.14 draw** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 251 | 0.933 | 4.4677 | 4.5442 | 46.6th pct | **1.0171×** | **1.2436×** |
| `orion_hct116` | 267 | 0.921 | 4.4934 | 4.5570 | 45.7th pct | **1.0142×** | **1.2453×** |
| `orion_hek293t` | 271 | 0.922 | 4.4718 | 4.5252 | 45.8th pct | **1.0120×** | **1.2381×** |
| `kolf2` | 253 | 0.927 | 4.4491 | 4.5172 | 46.2th pct | **1.0153×** | **1.2369×** |

**It is neither ~1.05 nor ~1.6 — it is ≈1.24.** Today the over-dose is **1.2–1.7 %** (a no-op, as
the INBOX expected). Under `docs/DEF-RUNBOOK.md` §5's pessimistic 174/300 draw, where the all-rows
median falls to the ~14th percentile of carrying rows, CD4's solved fraction would rise by
**1.237–1.245×**, taking an intended `f = 0.34` to an effective **f ≈ 0.42** — between 216's priced
`f = 0.40` and `f = 0.55` rungs, i.e. a cost somewhere between roughly −0.0003 and −0.0022 against
the argmax, silently, with the board dark. **This is a real October item, sized. It is not retired.**

## 11. INBOX-1 items 4 and 5

**4 — the h1-test leak.** Acknowledged and obeyed. **No h1-test surface was used anywhere in this
experiment.** Every number here comes from 162's PLOLO instrument in the challenge frame, from the
built packages themselves, or from the real `data/challenge/context_*.h5ad` controls.

**5 — the 10-vs-1 reconciliation.** Both are right; they count different things on different target
sets. 218 §7.2's **10** is on a PLOLO fold's target list, and counts rows where CD4 casts **no vote
at all**. The build's `targets_outside_cd4_support = **1**` is on the **300-target challenge
panel**, and counts rows outside CD4's **complete-mean3** support — the production rule that keeps a
CD4 row only when all three conditions (`cd4_rest`, `cd4_stim8hr`, `cd4_stim48hr`) carry the target,
which covers 276/300 here. On the challenge panel, 76 of the union's 77 carried targets fall inside
CD4's complete-mean3 support and **1** does not; the union also adds **740,065 gene-level entries**
on rows CD4 does reach but on genes it does not.

---

## 12. The package

| artifact | bytes | sha256 |
|---|---:|---|
| **`submissions/222-h1union-g028.vcc`** (the package) | 2,824,785,920 | **`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`** |
| `submissions/222-h1union-g028.h5ad` | 12,648,807,800 | `8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` |
| `submissions/222-h1union-g028.h5ad.meta.json` | 24,797 | — |
| `submissions/222-ka-cd4f034-f034.h5ad` (KA row 1, a proven duplicate of `210-cd4f034.h5ad`; deletable) | 12,648,807,800 | `7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a` |

**No model name chosen, no description written, nothing submitted.** §14 applies at submit time.

## 13. Files

| script | what it does | runtime / peak RSS |
|---|---|---|
| `stage222.py` | stages the production effect matrix from a builder given **by path**, without emitting — KA row 2 | 16–22 s / 19.01 GiB peak (`on028`) |
| `repro222.py` | 162's PLOLO instrument driven by the **production** union + sizing rule — KA row 3, 218's falsifiers 1 and 4 | 75 s / 8.27 GiB |
| `build222.sh` | `build210.sh` + one exported variable | 2085 s / 1775 s |
| `metadiff222.py` | falsifier 1, pre-declared sets | 1 s |
| `profile222.py` | 210's `profile210.py` + the carrying-row split — falsifiers 2 and 3 | 31 s / 3.60 GiB |
| `probe222.py` | 210's `probe210.py` — falsifier 4 (`n_pred`) | 43 s / 13.20 GiB |
| `six222.py` | 216's exact `nmae` identity on this arm — INBOX items 1 and 3 | 20 s / 7.85 GiB |
| `conc222.py` | leave-top-k-out concentration — INBOX item 2 | 12 s / 7.59 GiB |
| `work/make_submission_pm_cd4.PRE.py` | **the git control**, `6cc91501…`, taken before the first edit | — |

```bash
git show HEAD:scripts/make_submission_pm_cd4.py > experiments/222-h1-union-build/work/make_submission_pm_cd4.PRE.py
cd ~/Code/virtual-cell-challenge-2026
export VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34
uv run python experiments/222-h1-union-build/stage222.py --builder experiments/222-h1-union-build/work/make_submission_pm_cd4.PRE.py --tag pre
uv run python experiments/222-h1-union-build/stage222.py --builder scripts/make_submission_pm_cd4.py --tag off
VCC_H1_UNION_G=0.28 uv run python experiments/222-h1-union-build/stage222.py --builder scripts/make_submission_pm_cd4.py --tag on028
uv run python experiments/222-h1-union-build/repro222.py
experiments/222-h1-union-build/build222.sh submissions/222-ka-cd4f034-f034.h5ad
VCC_H1_UNION_G=0.28 experiments/222-h1-union-build/build222.sh submissions/222-h1union-g028.h5ad
uv run python experiments/222-h1-union-build/metadiff222.py submissions/222-h1union-g028.h5ad.meta.json submissions/210-cd4f034.h5ad.meta.json
uv run python experiments/222-h1-union-build/profile222.py submissions/222-h1union-g028.h5ad submissions/210-cd4f034.h5ad
uv run python experiments/222-h1-union-build/probe222.py cd4f034=submissions/210-cd4f034.h5ad h1union=submissions/222-h1union-g028.h5ad --perts=40
uv run python experiments/222-h1-union-build/six222.py
uv run python experiments/222-h1-union-build/conc222.py
uv run vcc prep submissions/222-h1union-g028.h5ad -g data/challenge/gene_names.csv \
  --perts data/challenge/pert_counts.csv -o submissions/222-h1union-g028.vcc     # 22.3 GiB RSS
```

## 14. Honest limits

- **No `cell-eval2` run and no board-side member measured.** `fid`, `jac` and `reach` carry `n_real`
  in a denominator or a union and are unreadable on every harness we hold (rule 9); they are priced
  only through 218's rate card. `mse` is not measured here at all — 197's inversion is the route and
  it needs the board.
- **The `nmae` instrument fails its family's board retrodiction** (§8). This is the sharpest caveat
  on the page and it cuts both ways: the arm's `nmae` could land positive as `cd4f034`'s did.
- **Δ`pds` is concentrated in three targets per fold** (§9). The pre-committed reading fires.
- **The delivered-count separation in §6b is 3.3×, not infinite**, because the emitter re-draws
  integers under a 1.2 % amplitude change. The exact statement lives at the effect-matrix level
  (§6a, 1e−16), and the only way to make the delivered statement exact would be a third build of
  `cd4f034`'s effects scaled by 0.988046 — 35 minutes for a check §6a already settles.
- **`g = 0.28` was not re-optimised here.** It is 218's headline dose and the addendum §4 justifies
  it over LOFO's 0.45 on the `pds` turn-over. I built what was preregistered.
- The `kolf2` PLOLO fold is excluded from every headline, as in 218; it is 3–5× more favourable
  because both sources and that recipient are pluripotent.
