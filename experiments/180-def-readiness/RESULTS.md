# 180 — D/E/F readiness: builder generalised; per-context knee VALIDATES locally and FAILS on the board.

**Board result (2026-09-02 23:28 UTC):** see `BOARD-RESULT.md`. Entry `N3ui4hCkKx9bxExG2pxE`,
overall `0.1740437` (+0.00068 vs a23). Local `mse` price +0.005 scaled; board Δ`mse` **−0.00226**.
**D/E/F default: global alpha, not knee.**

Date: **2026-09-02 → 2026-09-03**. Build + verification only. **No `vcc submit`. No
`cell-eval2`. Nothing written into `submissions/`.** `vcc prep` was run `--dry-run` only,
twice, and wrote nothing.

---

## Verdict in eight lines

1. **The builder is D/E/F-ready.** It reads its contexts from `manifest.json`, builds the
   bundle file list and the frozen-input table from them, and carries a documented
   `--repin-frozen-inputs` / `--repin-panel-support` pair. Rehearsed end to end on a
   synthetic D/E/F bundle through `vcc prep --dry-run`: **PASS**.
2. **The A/B/C build is unchanged under the patch, proven three ways.** The effects audit
   reproduces a23 on **179 compared fields with exactly one difference — the patched
   builder's own file hash** (§2). The D/E/F rehearsal's emitted counts are **byte-identical
   to `a23-massfix.h5ad` on `X/data`, `X/indices` and `X/indptr`**, differing only in the
   three context labels (§4). And the pre-change builder out of git emits a sha256-identical
   artifact at the default emitter through the real writer (§1, unit test).
3. **The knee alpha validates.** Computed from the bundle alone (no board): at 096's own
   `C_p` the knee alphas are **0.98586 / 1.01275 / 1.04568** against the transition plan's
   0.987 / 1.013 / 1.047, and the **unclamped** budget ratio at those alphas reads
   **0.999999 / 0.999800 / 0.999746** (§5). Six independent known-answer rows pass, three of
   them to 0.0e+00.
4. **`mse` prices at `+0.005025` scaled mean = `+0.000838` overall**, reproducing the
   transition plan's `≈ +0.005` / `≈ +0.0008` (§6). The plan's `≈ +0.03 scaled mse` for a
   context at `rho 0.70` also reproduces, at **`+0.02977`** (§8).
5. **`pds` is measured, not assumed: `−0.000076` raw** on PLOLO in the challenge frames
   (`−0.000013` overall), 20× inside 162's LOO sigma, with the instrument's non-vacuity
   established on the same folds (§7).
6. **`nmae` `+0.000193`; `reach` / `fid` / `jac` are UNRESOLVED** — on the same surface a
   re-seeded twin moves `reach` by 3.0× the arm's effect and moves `fid` and `jac` by about
   as much as the arm does (§7).
7. **Rule 13 vetoes the arm on A/B/C**: `predicted = −0.00161` against the amplitude class's
   `+0.0089` floor. Even summing every member optimistically gives `+0.0018`. The knee is a
   **D/E/F insurance policy, not an A/B/C arm** (§9).
8. **Two things I did not achieve**: peak RSS was **18.7–19.8 GiB**, not the mandated 10 GB —
   a production build cannot be done under 10 GB in this codebase and never has been (§11);
   and `VCC_ALPHA=knee` refuses `--emitter resample_mix` (§10).

**Package:** `experiments/180-def-readiness/work/abc-knee.h5ad`
sha256 **`f565ce0f5a38ce78b82759f01408aa560c4438d3bd714db0564ed9800647d1d0`**
19,057,758,808 bytes · 360,000 rows · nnz 2,378,696,726 (50.08% of the density cap)
builder sha256 `82f478a5c6f9960b0a3035cbe7625c4d155d7d7c3e1bfd2070d86ac12b74dec6`.
**Not submitted. The decision to spend a slot is Ahmet's, and rule 13 says do not.**

---

## 1. What changed in the builder

Three files: `scripts/make_submission_pm_cd4.py`, `scripts/build_submission.py` (untouched in
the end — `ChallengeSpec` already read contexts from the manifest), and
`scripts/test_make_submission_pm_cd4.py`. One new data file, `scripts/frozen_inputs.json`.

### Item 1 — the context adapter

| change | before | after |
|---|---|---|
| context guard | `raise unless spec.contexts == ("A","B","C")` | manifest's contexts, cross-checked against the frozen table and `--expect-contexts` |
| bundle file list | six hardcoded names | `bundle_file_names(contexts)` |
| frozen bundle hashes | six literals in source | `scripts/frozen_inputs.json`, `files` keyed `bundle/<name>` |
| perturbation count | `!= 300 -> raise` | `!= frozen_bundle["n_constructs"] -> raise` |
| cells per perturbation | unchecked in `run()` | `!= frozen_bundle["cells_per_pert"] -> raise` |
| PerturbMatch / CD4 coverage | `248` / `276` in source | `frozen_bundle["panel_support"]`, `null` = unpinned |
| re-pin | — | `--repin-frozen-inputs`, `--repin-panel-support`, `--frozen-inputs <path>` |

**Why `scripts/` and not `data/challenge/`.** `.gitignore` line 2 is `data/challenge/*`. A
frozen-input table whose recorded value is not under version control is not a freeze, so the
table lives next to the builder. The A/B/C digests were moved verbatim and the test asserts
all six literally.

**The three-way context cross-check.** `manifest.json` has to be read before anything can be
hashed (it is what names the `context_*.h5ad` files), so it is read untrusted, its context
list is compared against the frozen table's **and** against `--expect-contexts`, and only then
is every bundle file — including `manifest.json` itself — hashed. A tampered manifest that
renames the round aborts before a single byte is read.

**`--repin-frozen-inputs` never builds** and `--repin-panel-support` exits before emission. A
build never re-pins itself. `panel_support` is carried across a re-pin only when
`pert_counts.csv` and `gene_names.csv` are byte-identical, and set to `null` otherwise, because
that pair is a property of (panel × donor tables) and the A/B/C values would abort every D/E/F
build on a new 300-target panel. **This was not in the mandate and is a real blocker it did not
name.**

### Item 2 — `VCC_ALPHA=knee`

`VCC_ALPHA` now takes a number (default `1.0796`, unchanged) or the literal `knee`.
`VCC_ALPHA_REF` (default `1.0796`) is the amplitude the effect matrix is first built at and the
anchor the knee is expressed against:

```
alpha_c = alpha_ref * sqrt( claim_c / B_c(alpha_ref) )
```

**On the mandate's "normalised so that the current A/B/C mean reproduces the shipped scale".**
Read literally as "rescale so the mean of the three knee alphas equals 0.982", that gives
0.954 / 0.979 / 1.012 — which are not the knee, do not put `rho` at 1, and contradict the same
paragraph's validation targets. The reading that survives its own known-answer row is
anchoring: if every context were already on its knee the rule returns `alpha_ref` exactly, so
the shipped scale is reproduced whenever the correction is null. That is what is implemented,
and §5 is the row that decides between the two.

The solve is done **inside the builder from the bundle alone**: `B_c` from our own rendered
prediction, `Ĉ_p` from our own emitted cells, `C_p` from the context's own NTC control cells.
`cell_eval2` is imported lazily and only on the knee path.

The knee is then applied by **rebuilding the effect matrix per context** through
`build_production_effects(..., alpha=alpha_c)` — not by scaling the finished one. The pipeline
happens to be homogeneous of degree 1 in alpha under a15/a23's flags and the built meta
measures that (`homogeneity_residual` **5.9e-08 / 7.3e-08 / 6.7e-08**), but homogeneity is a
property of the current configuration, not of the knob: `VCC_CD4_FRACTION` alone breaks it (a
fixed absolute CD4 residual does not scale with alpha), and §7 shows that difference is
measurable on `pds`.

`effects` and its audit are left exactly as built at `alpha_ref`, so the rule-10 effects diff
against the incumbent stays clean; the per-context matrices are recorded separately under
`meta['knee']` with their own sha256s.

### Gate 1a — the test suite: **57 tests pass** (51 pre-existing + 6 new), peak RSS **0.31 GiB**

| new test | what it asserts, and what it reads if the thing is absent |
|---|---|
| `test_bundle_file_list_follows_the_manifest_and_not_a_hardcoded_abc` | `bundle_file_names(["D","E","F"])`, `read_bundle_contexts`, and `frozen_input_paths` on a D/E/F argument. Absent ⇒ raises, or points a D/E/F build at `context_A.h5ad`. |
| `test_frozen_bundle_table_fails_closed_on_every_shape_error` | 12 malformed tables all rejected (contexts/files desync in both directions, short digest, non-hex digest, bad `panel_support`); `panel_support: null` accepted as a legal state; **all six A/B/C digests asserted literally**. |
| `test_repin_frozen_inputs_moves_the_table_to_a_def_bundle_and_never_builds` | (a) the A/B/C table cannot even verify a D/E/F path set; (b) the re-pin moves contexts and drops `panel_support`; (c) the new table verifies the bundle; (d) one changed byte still aborts; (e) an unchanged panel keeps `panel_support`. Also asserts the **production** table is byte-unchanged. |
| `test_context_effects_default_is_none_and_only_moves_its_own_context` | through the REAL writer: an all-identical `context_effects` map is sha256-identical to `None`; moving ONE context's matrix moves that context's 400-cell blocks and **no other's**; three malformed maps rejected. |
| `test_alpha_knee_env_is_parsed_and_fails_closed` | `knee` / `KNEE` parsed, default still `1.0796` fixed, six malformed env combinations rejected including `VCC_ALPHA_REF` without `knee`. |
| `test_knee_puts_the_unclamped_budget_ratio_on_one_and_is_not_a_no_op` | on a synthetic context placed at `rho_unclamped = 0.80`: the knee scale lands it at 1.000 ± 0.02; over-correcting by the same factor lands **above** the knee where `rho` clamps to 1 and `rho_unclamped` does not; the `min(Ĉ, C_p)` branch collapses `claim` onto the cap. Asserts the start is ≥ 0.15 from 1 so "it lands on 1" is not a statement about the starting point. |

The pre-existing `test_default_emitter_reproduces_the_pre_change_builder_byte_for_byte` — which
resolves the last pre-`--emitter` builder revision **out of git** and compares emitted `.h5ad`
sha256 through the real writer — still passes, and now also exercises the manifest the
synthetic bundle gained.

---

## 2. Gate 1 — the rule-10 bit-identity proof on the real A/B/C run

`effects180.py` drives the production `run()` end to end at a23's exact configuration with
**only `write_submission_from_effects` stubbed out**. Everything the patch touches runs for
real: `load_frozen_bundle`, the three-way context cross-check, `collect_input_hashes` over all
34 frozen inputs, the donor model, PerturbMatch, the pruned operator, the whitening, the CD4
residual, the panel-support assertion.

```
[contexts] ['A', 'B', 'C'] (manifest, frozen table and --expect-contexts all agree)

fields compared: 179       effects.* 54    operator.* 24    (+ config, + hashes)
DIFF hashes.implementation.scripts/make_submission_pm_cd4.py
       built     82f478a5c6f9960b...      incumbent 0852378b7616e719...

effects.final_sha256   22c09f22c749bda82d12a0862bacf25b55e709b7de697513da3e3d11288a9b97
a23  final_sha256      22c09f22c749bda82d12a0862bacf25b55e709b7de697513da3e3d11288a9b97

GATE 1 PASS: 1 difference, 0 outside the patched-file hash.   peak RSS 19.660 GiB
```

The one difference is the file that was patched — and it is the SAME revision that built the
package in §9 (`82f478a5…`), so this gate and the shipped artifact are one builder, not two.
`hashes.inputs` (all 34 frozen inputs) and `config` are identical. The context adapter and the
frozen-inputs split cannot reach the effect matrix.

Operator by path and sha, never by description (rule 11):
`data/donors/whitener_union_pruned.npz`, sha256
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, 17,942 columns,
`all_panel_coordinates_excluded: true`.

---

## 3. The synthetic D/E/F rehearsal — the full log

`defbundle180.py` copies the real controls, renames the files **and** `obs['context']` to D/E/F
through h5py in place (so the job stays under a few hundred MB rather than loading three 110M-
nonzero matrices), and rewrites `manifest.json`.

```
context_A.h5ad -> context_D.h5ad  obs.context ['D']  18400 rows  nnz 109910478
  source sha f22e71968487d3da...   target sha a3b6621302991881...
context_B.h5ad -> context_E.h5ad  obs.context ['E']  18400 rows  nnz 101546777
  source sha 41557555b3febf6d...   target sha c6f58df46cf27d65...
context_C.h5ad -> context_F.h5ad  obs.context ['F']  18400 rows  nnz 108213983
  source sha e090e7c7cfdd0006...   target sha 81af8449ba23e39d...
manifest.json  sha 42c0fffebbff886c...   (panel_id vcc2026-synthetic-def-1)
```

**Step 1 — re-pin the file hashes**, into a throwaway table so a concurrent A/B/C build is not
disturbed (`--frozen-inputs`; the code path is identical, only the file differs):

```
[repin] contexts ['A', 'B', 'C'] -> ['D', 'E', 'F']
[repin]   bundle/context_A.h5ad: f22e7196... -> <absent>      (x3 removed)
[repin]   bundle/context_D.h5ad: <absent> -> a3b66213...      (x3 added)
[repin]   bundle/manifest.json:  2c2b5d42... -> 42c0fffe...
[repin] panel_support: {276, 248} -> {276, 248}
```

`panel_support` was correctly PRESERVED here, because the synthetic bundle reuses the real
`pert_counts.csv` — the rule is "carry it only if the panel files are byte-identical", and they
are. To exercise the other branch it was then nulled by hand and re-pinned from a build:

**Step 2 — `--repin-panel-support`** (effects build, exits before emission):

```
[contexts] ['D', 'E', 'F'] (manifest, frozen table and --expect-contexts all agree)
[repin] panel_support pinned to {'perturbmatch_replacement_targets': 248,
        'cd4_complete_targets': 276}; no submission was built
```

**Step 3 — the full a23 configuration end to end** (`VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10
VCC_MASS_FIX=minnorm`, `--min-cpm 10 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25
--overdispersion 4.0 --whitener-id union_pruned`, `--expect-contexts D,E,F`): 13 minutes, peak
RSS 19.77 GiB, 19,057,758,808 bytes, nnz 2,378,695,978.

**Step 4 — `vcc prep --dry-run`. It FAILED first, and that is the most useful thing in this
section:**

```
vcc prep: Unknown context label(s) in 'context': D, E, F.
Valid labels are A, B, C — reuse the labels from the control files exactly.
```

**`vcc prep` hardcodes A/B/C too** (`vcc/prep.py:53`, `REQUIRED_CONTEXTS = ("A","B","C")`).
There is an undocumented-in-`--help`-tail `--contexts` option (`vcc/cli.py:139`, default
`"A,B,C"`) that fixes it. With it:

```
✓ dry run — no .vcc written
  cells:         360000
  genes:         18533  (encoding: float32)
  contexts:      column 'context' → D=120000, E=120000, F=120000
  targets:       verified against the official list
  normalization: counts-preserved
```

**On Oct 22, `vcc prep` needs `--contexts D,E,F` (or a CLI upgrade), or it rejects the
submission before it is ever uploaded.** This is on the checklist.

---

## 4. Gate 3 — the D/E/F rehearsal delivered a23's counts, relabelled and nothing else

`defidentity180.py`, streaming both 19 GB files:

```
X/data       n 2,378,695,978 vs 2,378,695,978   IDENTICAL   2d6f6dd044d67060
X/indices    n 2,378,695,978 vs 2,378,695,978   IDENTICAL   0e6e4b9c02bff20d
X/indptr     n       360,001 vs       360,001   IDENTICAL   c75ef4c31605a3ea
obs/context      categories ['A','B','C'] -> ['D','E','F'];  codes IDENTICAL
obs/target_gene  300 cats -> 300 cats;                       codes IDENTICAL
var index identical: True (18533 genes)

GATE 3 PASS
```

**What this would read if the generalisation were wrong.** If the adapter had paired context D
with `context_E.h5ad`'s controls, reordered the emission loop, or re-seeded a stream, `X/data`
would move — the emission RNG is a single stream across all 900 blocks and the control profile
is per context. It matches on all three arrays. The `obs/context` row is the other half: if
the labels had NOT moved, the two files would simply be the same file and the identity would
prove nothing.

---

## 5. Gate 2 — the knee alpha, and its six known-answer rows

`knee180.py` runs the **production** knee code (`load_control_frame`, `knee_context_terms`,
`solve_knee_alphas_with_columns`) offline against the effect matrix §2 proved is a23's.

| # | known-answer row | result |
|---|---|---|
| 1 | the effect matrix's sha256 = a23's `effects.final_sha256` | `22c09f22c749bda8…` **exact** |
| 2 | the `C_p` NTC proxy reproduces 096's cached 16-construct vector | **max abs diff 0.0e+00 on all three contexts** |
| 3 | `A`, `B` reproduce 170's `chal_minnorm.json` gamma-0.65 row | `dA = 0.0e+00`, `dB ≤ 2.8e-17` |
| 4 | `claim`, `rho` reproduce 170 **at 096's own `C_p`** | `dclaim ≤ 2.8e-17`, `drho ≤ 2.2e-16` |
| 5 | the knee alphas land on the transition plan's 0.987 / 1.013 / 1.047 | **0.98586 / 1.01275 / 1.04568** (max gap 0.0014) |
| 6 | `X` solved from a23's board point reproduces 170's | `0.017932` vs `0.017932` |

**The `C_p` estimator moved, deliberately, and it is the one judgement call here.** 096 and 170
used the **first 16** ntc constructs; the builder uses **all 46**, because "first 16 in file
order" is not a rule that survives contact with a bundle nobody has seen and 46 has the smaller
standard error. The gate is therefore decomposed: `A` and `B` are our side and reproduce 170
outright, while `claim` and `rho` are re-derived at 096's own `C_p` from the per-row `Ĉ` values
(now recorded in the meta) and reproduce it to 3e-17.

### The per-context knee table

`A`, `B`, `claim`, `Ĉ`, `C_p` measured; `Kc` from the control jackknife; `X` and `sumD`
board-anchored.

| ctx | `A` | `B` | `claim` (all-46) | `claim` (096 `C_p`) | `C_p` all-46 | `C_p` 096-16 | `Ĉ` | cap binds | `Kc` |
|---|---:|---:|---:|---:|---:|---:|---:|:-:|---:|
| A | 0.234830 | 0.233593 | 0.233754 | 0.235434 | 14.4398 | 14.5435 | 17.075 | 24/24 | 0.005084 |
| B | 0.224778 | 0.223624 | 0.237923 | 0.237847 | 14.8057 | 14.7716 | 14.704 | 3/24 | 0.005228 |
| C | 0.234941 | 0.233728 | 0.266129 | 0.265026 | 16.4397 | 16.3715 | 17.723 | 24/24 | 0.005919 |

| ctx | `rho` before | `alpha` (all-46, **shipped**) | `alpha` (096 `C_p`) | plan | ±1 s.e. of `C_p` | `rho` after |
|---|---:|---:|---:|---:|---:|---:|
| A | **0.999313** | 0.982338 | 0.985862 | 0.987 | ±0.0019 | **0.999999** |
| B | **0.939899** | 1.012910 | 1.012747 | 1.013 | ±0.0024 | **0.999800** |
| C | **0.878250** | 1.047858 | 1.045684 | 1.047 | ±0.0044 | **0.999746** |

`rho` here is the **unclamped** ratio `B / claim`. The scorer's own `rho = min(1, B/claim)`
reads 1.0 both **on** the knee and **above** it, so "it reads 1 after the correction" cannot
tell those two apart and is not a usable check. Every check in this feature reads the unclamped
ratio; the unit test in §1 exercises exactly that asymmetry.

**The honest caveat on context A.** A sits essentially on its knee already, and its correction
(0.982 → 0.9823 under all-46, → 0.9859 under 096's 16) is **within 0.2 and 2.0 standard
errors of the `C_p` proxy's own ±0.0019**. A's knee alpha is unresolved at the 0.4% level; B's `+3.1%` and
C's `+6.7%` are 13× and 15× their own standard errors and are resolved. 096 already stated the
proxy is stable to ~3%; `alpha ∝ sqrt(C_p)` halves that to ~0.2% on A/B/C's 46 constructs.

### The verification is not an identity, and the builder's own run proves the wiring

`rho_unclamped` returning to 1.000 would be near-tautological if `B` were exactly quadratic in
the amplitude — it is not, because `profile_from_effect` is nonlinear (log1p/expm1, the CPM
gate, and a23's minimum-norm mass fix). More importantly, the **production build's** own
verification pass re-measures on the matrices it **rebuilt** through
`build_production_effects(alpha=alpha_c)`, not on scaled copies:

```
[knee] A: after rebuild at alpha 0.982338 -> rho 1.0000 (was 0.9993), B 0.233754 vs claim 0.233754
[knee] B: after rebuild at alpha 1.012910 -> rho 0.9998 (was 0.9399), B 0.237898 vs claim 0.237946
[knee] C: after rebuild at alpha 1.047858 -> rho 0.9997 (was 0.8783), B 0.266061 vs claim 0.266129
```

If the per-context alpha never reached the effect matrix, these would read the pre-correction
0.9993 / 0.9399 / 0.8783 back. They do not. The meta additionally records
`frobenius_ratio_to_reference` = 1.000344 / 1.031477 / 1.067065 and the solved CD4 fractions
0.094760 / 0.097709 / 0.101080 against the reference build's 0.094727 — the residual really is
re-solved per context.

---

## 6. `mse`, through 096's claim route

```
board anchor a23 score_mse 0.121581 -> raw 0.872341;  solved X 0.017932  (170 solved 0.017932)
X at the knee 0.018520 (mean kappa 1.032828)
raw a23 0.872341 -> knee 0.867566     d score_mse +0.005025     d overall +0.000838
```

| ctx | `num` before | `num` after | Δ raw `mse` | Δ scaled `mse` |
|---|---:|---:|---:|---:|
| A | +0.006321 | +0.006322 | −0.000050 | **+0.000053** |
| B | +0.006382 | +0.006457 | −0.004582 | **+0.004822** |
| C | +0.007133 | +0.007306 | −0.009693 | **+0.010202** |
| mean | | | −0.004775 | **+0.005025** → **+0.000838 overall** |

**A correction worth recording, because it is the easiest way to get this wrong.** `X` is
`Σ_p w_p ⟨dh_p, D_p⟩` and is **linear** in the delivered effect, so it does not survive an
amplitude change unchanged. Holding it fixed — which the first version of this script did —
turns `+0.000838` into `−0.000064`, i.e. it flips the sign of the whole result. The multiplier
is the per-row ratio `kappa_p = ⟨dh_p(s), dh_p(1)⟩ / ‖dh_p(1)‖²`, measured on the delivered
rows:

| ctx | `kappa` | sd | `kappa / s` max deviation | max `1 − cos(dh(s), dh(1))` full axis | panel-excluded axis |
|---|---:|---:|---:|---:|---:|
| A | 1.000343 | 1.9e-06 | 2.2e-05 | 3.4e-09 | 2.3e-09 |
| B | 1.031375 | 3.4e-04 | 3.3e-03 | 4.5e-05 | 4.7e-05 |
| C | 1.066766 | 1.1e-03 | 7.2e-03 | 2.0e-04 | 1.3e-04 |

`sumD = 0.229141` is 096's board-solved constant; `X` is solved from a23's own board point
(one equation, one unknown) and, as 096 §(e) states, cannot be separated per context — the same
`X` is used for all three and the per-context split above inherits that assumption.

---

## 7. All six members priced (rule 12)

| member | value (overall) | instrument | confidence |
|---|---:|---|---|
| `mse` | **+0.000838** | 096's `claim` route, per context, on the challenge panel | good — six known-answer rows, §5 |
| `pds` | **−0.000013** | **PLOLO**, challenge frames, delivered space, kolf2 fold | good — the rule-13 instrument; 20× inside its own LOO sigma |
| `nmae` | **+0.000193** | h1-test, mean of the three per-context alphas | fair — 30× the twin's noise, but one context |
| `reach` | **UNRESOLVED** (+0.000585 point estimate) | h1-test realised p-order | **poor** — the twin moves `reach` by 3.0× the effect |
| `fid` | **UNRESOLVED** (+0.000135) | h1-test realised call set | poor — twin noise ≈ the effect; rule 9 |
| `jac` | **UNRESOLVED** (+0.000072) | h1-test realised call set | poor — twin noise ≈ the effect; rule 9 |

**No member is unmeasured. Three are measured and unresolved, which is a different and weaker
statement, and it is stated as such.**

### `pds` — PLOLO, and why the screen is not vacuous

`pds` is `pds_cosine` with `exclusion_scope="panel"`: every distance it ranks is a cosine
between our delivered delta row and a reference delta row, in one fixed feature space with all
300 panel targets dropped. A cosine is invariant to a positive per-row scale, so "a per-context
alpha cannot move `pds`" is *close* to an identity — close, because a23's minimum-norm mass fix
makes the renderer non-homogeneous and §6 measures the delivered direction moving by up to
`1 − cos = 1.3e-4`. **A screen that can only ever return 0.000000 proves nothing**, so the
instrument's response was established on the same folds and frames:

* **known-answer**: `a23` re-run reproduces 162's recorded `arms.json` `pds` to **0.0e+00 on
  6/6** (fold, frame) pairs;
* **the instrument moves for real geometry changes**, from 162's own record on these folds:
  `a15 − a23` (the mass fix) `+0.00058…+0.0078`, `a19 − a15` `−0.00063…−0.0032`,
  `a11 − a14` (alpha 0.838 vs 0.982) `+0.00065…+0.0023`.

```
delta pds against a23 (162's LOO sigma on this instrument is 0.00149 raw):
  kolf2      knee_A  chalA +0.000000   knee_B  chalB -0.000108   knee_C  chalC -0.000121
  k562_gwps  knee_A  chalA -0.000014   knee_B  chalB +0.000069   knee_C  chalC +0.000028

  the arm as it would ship — context c at alpha_c in frame chal c:
  kolf2      mean -0.000076 raw pds  = -0.000013 overall
  k562_gwps  mean +0.000028 raw pds  = +0.000005 overall
```

**A mechanism that falls out, and it is worth keeping.** Two deliberately large amplitude arms
at a23's exact configuration (`alpha = 0.90` and `1.10`) also read ≈ 0 (`−0.00039…+0.00018`),
while 162's `a11` vs `a14` — the same alpha channel — reads `−0.0016`. The difference is
`VCC_CD4_FRACTION` vs `VCC_CD4_TARGET_F`: with a solved target `f` the whole pipeline is
homogeneous in alpha and the arm is a pure per-row scale, which `pds_cosine` cannot see; with a
**fixed** absolute CD4 residual, changing alpha changes the residual's *relative* weight, which
is a direction change. **Amplitude is `pds`-free only because a15 moved to `VCC_CD4_TARGET_F`.**

### `nmae` / `reach` / `fid` / `jac` — h1-test, and its limits

**h1-test is one context and its own `rho` is 1.0** — it sits *above* the knee — so it cannot
price `mse` at all (170 §4) and a per-context knee does not exist on it. Only the direction and
rough size of each DE member's response to a pure amplitude change transfers. Each arm is
`170/work/eff/h1_g0.65.npz × (alpha/0.982)` through 170's own `deliver` / `emit_and_score` /
`Surface` stack, unchanged. 170's `mse` known-answer gate (anchored on a real `cell-eval2`
`agg_results.csv`) runs first, and the `s = 1.0` row reproduces 170's gamma-0.65 row on **all
six members to 0.000e+00**.

| arm | `s` | `pds` | `nmae` | `reach` | `fid` | `jac` | `n_pred` |
|---|---:|---:|---:|---:|---:|---:|---:|
| anchor (a23) | 1.000000 | 0.87535 | 0.95988 | 0.17752 | 0.50816 | 0.06425 | 3803 |
| knee A | 1.000344 | 0.87535 | 0.95985 | 0.18478 | 0.50832 | 0.06426 | 3806 |
| knee B | 1.031477 | 0.87465 | 0.95918 | 0.17967 | 0.50843 | 0.06431 | 3876 |
| knee C | 1.067065 | 0.87409 | 0.95853 | 0.17750 | 0.50846 | 0.06463 | 3959 |
| knee mean | 1.032962 | 0.87437 | 0.95912 | 0.18009 | 0.50828 | 0.06435 | 3880 |

**The re-seeded twin, which is what makes three of these unresolved:**

| arm | Δ`pds` | Δ`nmae` | Δ`reach` | Δ`fid` | Δ`jac` |
|---|---:|---:|---:|---:|---:|
| anchor, seed 1 − seed 0 | **+0.00000** | +0.00001 | **−0.00957** | +0.00016 | +0.00006 |
| knee mean, seed 1 − seed 0 | **+0.00000** | −0.00002 | **−0.00504** | +0.00027 | +0.00015 |

`pds` is seed-free on this surface; `nmae` is essentially so. `reach` moves by `0.005–0.010`
raw between seeds against an arm effect of `+0.0031`, and `fid`/`jac` move by as much as the
arm does. Those three are noise on this instrument at this effect size.

h1-test reads `pds` at `−0.00065` raw for the arm as shipped (`−0.00024` overall), 9× larger
than PLOLO's `−0.000076`. **PLOLO is the `pds` instrument** (rule 13; 178 measured h1-test
mis-stating `pds` by ~6× in the other direction on un-whitening arms), so `−0.000013` is the
number carried into the gate and h1-test's is recorded as the disagreement it is.

---

## 8. What the rule is actually for — a D/E/F context off its knee

`forward180.py`, using 096's identity with context C's own measured `A`, `B`, `Kc` and the
A/B/C-anchored `X` and `sumD`, placing one context at each `rho`:

| `rho` | knee `s*` | `alpha*` | Δ raw `mse` | **Δ scaled `mse`** | Δ overall (1 ctx of 3) |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 1.00000 | 0.98200 | 0.00000 | +0.00000 | +0.000000 |
| 0.90 | 1.05409 | 1.03512 | −0.00788 | +0.00829 | +0.000461 |
| 0.80 | 1.11803 | 1.09791 | −0.01715 | +0.01805 | +0.001003 |
| **0.70** | 1.19523 | 1.17371 | −0.02829 | **+0.02977** | +0.001654 |
| 0.50 | 1.41421 | 1.38876 | −0.05953 | +0.06266 | +0.003481 |
| **1.10** | 0.95346 | 0.93630 | −0.08593 | **+0.09043** | +0.005024 |
| 1.25 | 0.89443 | 0.87833 | −0.18854 | +0.19843 | +0.011024 |

The `rho = 0.70` row reproduces the transition plan's `≈ +0.03 scaled mse on that context`.

**The curve is asymmetric, and the high side is where the money is.** Below the knee the
un-refunded part of the amplitude is only `A − B` (0.0012 here, 0.5% of `A`), so being 30%
short costs 0.030 of scaled `mse`; above the knee the whole `2sA` is charged with no refund at
all, so being **10% long costs 0.090**. A global alpha fitted on A/B/C has no way to know which
side a D/E/F context lands on, and the expensive side is the one nobody has been thinking about.

`X` and `sumD` are A/B/C's; a D/E/F context will have its own and neither is knowable before
the round. This is the **shape** of the gain sized on our current alignment, not a forecast.

---

## 9. Rule 13 — the composed slot gate

The knee is an **amplitude-changing** arm, so its floor is `+0.0089`.

```
PLOLO d_pds / 6 * 1.011              -0.000076 / 6 * 1.011  =  -0.0000128
096 d_mse (amplitude moves)                                  =  +0.000838
emitted-profile toll                          -0.0146 / 6    =  -0.002433
                                                             ------------
predicted_overall                                            =  -0.001608

SHIP only if > +0.0089 (amplitude class)                     ->  VETO
```

Even the most generous possible sum — every h1-test point estimate taken at face value despite
being inside its own noise (`mse +0.000838`, `nmae +0.000193`, `reach +0.000585`,
`fid +0.000135`, `jac +0.000072`, `pds −0.000013`) — gives **`+0.00181`**, still 5× below the
floor. The toll is arguably not chargeable here (`reach` *was* priced on h1-test directly, and
the arm changes neither the emission law nor the per-cell library sizes), which is why both
readings are shown; neither clears the bar.

**This is the expected answer and not a disappointment.** The global 0.982 was fitted on A/B/C
board points, so A/B/C are near their knees by construction; §8 is where the rule pays.

### The built package

```
experiments/180-def-readiness/work/abc-knee.h5ad
sha256  f565ce0f5a38ce78b82759f01408aa560c4438d3bd714db0564ed9800647d1d0
19,057,758,808 bytes · 360,000 rows · nnz 2,378,696,726 (50.08% of the cap)
peak RSS 18.67 GiB      builder sha256 82f478a5c6f9960b0a3035cbe7625c4d155d7d7c3e1bfd2070d86ac12b74dec6
```

**Gate 4 — the meta diff against a23 (`metadiff180.py`): exactly three config keys, all
alpha-related.**

```
[config] keys that moved: alpha_mode, alpha_per_context, alpha_ref
[config]   alpha_mode          built 'knee'               incumbent <absent>
[config]   alpha_ref           built 0.982                incumbent <absent>
[config]   alpha_per_context.A built 0.9823377082708273   incumbent <absent>
[config]   alpha_per_context.B built 1.0129099639546535   incumbent <absent>
[config]   alpha_per_context.C built 1.0478576713641228   incumbent <absent>
[effects]  fields 54, differences 0     final_sha256 22c09f22c749bda8… (== a23)
GATE 4 PASS
```

Everything else that moved is outside `config`, stated rather than omitted:

| key | built | a23 | why |
|---|---|---|---|
| `hashes.implementation.scripts/make_submission_pm_cd4.py` | `82f478a5…` | `0852378b…` | the file that was patched; not in `FROZEN_IMPLEMENTATION_SHA256`, and the frozen modules including `donor_transfer.py` are unchanged |
| `output.nnz` | 2,378,696,726 | 2,378,695,978 | **+748**, i.e. +0.00003% — the amplitude change moves a handful of stochastic-rounding outcomes |
| `output.contexts.{A,B,C}.stored_entries` | 816.84 / 771.24 / 790.61 M | 816.82 / 771.23 / 790.65 M | the arm |
| `output.{path, bytes, sha256, nnz_per_cell, pct_of_density_cap}` | — | — | the arm (bytes identical at 19,057,758,808) |
| `frozen_bundle.*` | new block | absent | new: the table's path, sha256, contexts, `panel_support` and the observed pair |
| `knee.*` | new block | absent | new: every per-context term (§5) |
| `peak_rss_gib` | 18.67 | 17.47 | page-cache state, not a code change |

**`vcc prep --dry-run` on the package: PASS** — 360,000 cells, 18,533 genes, A/B/C = 120,000
each, targets verified against the official list, counts preserved. Nothing written.

---

## 10. What remains for Oct 22 — the checklist

| # | item | status | how |
|---|---|---|---|
| 1 | **`vcc prep` needs `--contexts D,E,F`** | **found here, unfixed** | the CLI's `REQUIRED_CONTEXTS` is `("A","B","C")`; pass `--contexts D,E,F`, or re-check after `vcc` is upgraded (`vcc skill install` after every upgrade) |
| 2 | re-pin the frozen inputs | ready, rehearsed | `--repin-frozen-inputs --bundle data/challenge`; review the printed diff |
| 3 | re-pin `panel_support` | ready, rehearsed | `--repin-panel-support` (one effects build, ~4 min); the A/B/C 248 / 276 **will** move on a new panel and would otherwise abort every build |
| 4 | **donor coverage re-audit** | not run — the panel does not exist yet | `uv run python experiments/166-donor-acquisition/census.py` against the new `pert_counts.csv`. K562 covers 269–272/300 of the *current* panel, Orion 290–294; if that drops materially the equal-weight consensus changes character (166) |
| 5 | gene axis | asserted by the builder | `ChallengeSpec.load` rejects ≠ 18,533 or duplicates; `write_submission_from_effects` rejects a control file whose `var_names` differ from the official order, per context |
| 6 | `cells_per_pert` | **now asserted** | read from the manifest and cross-checked against `frozen_inputs.json`; a silent change from 400 aborts |
| 7 | `VCC_ALPHA=knee` with a resampling emitter | **refused, by design** | `knee` measures `Ĉ` through the quasi emitter only and raises on `--emitter resample_mix`; if 177's arm is ever revived, `Ĉ` needs its own probe under that law |
| 8 | `n_real` on D/E/F | unmeasurable (rule 9) | unchanged from the transition plan item 5 |
| 9 | `EXPECTED_PM_REPLACEMENT_TARGETS` / `EXPECTED_CD4_COMPLETE_TARGETS` | left in source as the A/B/C values | the build no longer reads them; they are retained because the pre-existing test asserts them and they are the A/B/C known answer |
| 10 | peak RSS | 18.7–19.8 GiB per build | see §11 — this is not reducible without work nobody has scoped |

**Reference commands, for the day.**

```bash
# 1. re-pin the bundle (never builds; prints old -> new)
uv run python scripts/make_submission_pm_cd4.py --repin-frozen-inputs --bundle data/challenge

# 2. pin the new panel's PerturbMatch / CD4 coverage (effects build, exits before emission)
uv run python scripts/make_submission_pm_cd4.py --repin-panel-support --bundle data/challenge \
    --min-cpm 10.0 --whitener-id union_pruned \
    --whitener data/donors/whitener_union_pruned.npz \
    --whitener-metadata data/donors/whitener_union_pruned.json --max-rss-gib 26.0

# 3. re-audit donor coverage on the new panel
uv run python experiments/166-donor-acquisition/census.py

# 4. build (add VCC_ALPHA=knee VCC_ALPHA_REF=0.982 for the per-context knee)
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm \
uv run python scripts/make_submission_pm_cd4.py --bundle data/challenge \
    --expect-contexts D,E,F -o submissions/<name>.h5ad \
    --min-cpm 10.0 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 \
    --overdispersion 4.0 --whitener-id union_pruned \
    --whitener data/donors/whitener_union_pruned.npz \
    --whitener-metadata data/donors/whitener_union_pruned.json --max-rss-gib 26.0

# 5. validate (NOTE --contexts)
vcc prep submissions/<name>.h5ad -g data/challenge/gene_names.csv \
    --perts data/challenge/pert_counts.csv --contexts D,E,F \
    -o submissions/<name>.vcc --dry-run
```

---

## 11. What was NOT measured, and the constraint I could not meet

* **Peak RSS was 18.7 / 19.7 / 19.8 GiB, not the mandated 10 GB.** The mandate's own deliverable
  section says 20 GB and its constraints section says 10 GB; the two are inconsistent and a
  production build cannot meet the lower one. The floor is the three-donor model load —
  `orion_hct116` and `orion_hek293t` are 2.74 and 2.78 GB of `delta.npy` each — and it has been
  13.74 GiB (179) and 17.47 GiB (a23) on the two most recent production builds by other
  agents. An 18 GiB ceiling aborted the CD4 stage twice, so every heavy run used
  `--max-rss-gib 26.0`. **Every
  heavy job was run strictly one at a time**, and every analysis job (tests 0.31 GiB, `C_p`
  probe, knee solve, h1-test 6.62 GiB, PLOLO, the two 19 GB streaming digests) stayed well
  under 10 GiB. This is stated rather than worked around.
* **No board number was produced or predicted.** The `mse` route is anchored on a23's already-
  scored board point; it forecasts a delta, not a score.
* **`X_c` is not separable per context** (096 §(e)). The per-context `mse` split in §6 holds
  `X` at its pooled value and inherits that assumption; only the panel mean `+0.000838` is
  free of it.
* **`fid` / `jac` / `reach` carry `n_real`** (rule 9). Their h1-test values are reported as
  deltas only, and all three are inside their own re-seed noise.
* **The knee was not measured on a bundle whose panel differs.** The synthetic D/E/F bundle
  reuses the real `pert_counts.csv`, so `panel_support` did not have to move and the donor
  coverage did not change. That is checklist item 4 and it cannot be rehearsed before Oct 22.
* **`--emitter resample_mix` + `VCC_ALPHA=knee` is refused, not supported.** Both switches
  exist; composing them needs a `Ĉ` probe under the resampling law.
* **No arm was submitted, no slot was spent, `cell-eval2` was never run.**

---

## 12. Files

| file | what |
|---|---|
| `cpproxy180.py` | reproduces 096's `C_p` NTC proxy from the bundle and pins the selection rule (`work/cp_proxy_46.json`) |
| `effects180.py` | gate 1 — the real `run()` at a23's config with emission stubbed; 179-field meta comparison (`work/gate1.json`, `work/effects-a23.npy`) |
| `defbundle180.py` | builds the synthetic D/E/F bundle (`work/def-bundle/`, `work/def-bundle-hashes.json`) |
| `defidentity180.py` | gate 3 — the two 19 GB artifacts' `X` arrays and obs categoricals (`work/defidentity180.json`) |
| `knee180.py` | gate 2 — the knee, its six known-answer rows, the delivered geometry and the `mse` price (`work/knee180.json`) |
| `h1alpha180.py` | the DE members on h1-test at the three alphas, plus a re-seeded twin (`work/h1alpha180.json`) |
| `plolo180.py` | `pds` on 162's PLOLO in the challenge frames, with its non-vacuity arms (`work/plolo180.json`) |
| `forward180.py` | what the rule is worth off the knee (`work/forward180.json`) |
| `metadiff180.py` | gate 4 — the built package's meta against a23's (`work/metadiff180.json`) |
| `work/abc-knee.h5ad` | the built package, `f565ce0f…`, **not submitted** |
| `work/def-rehearsal.h5ad` | the D/E/F rehearsal artifact, `6404d006…`, **not submitted** |
| `work/frozen_inputs.def.json` | the throwaway D/E/F frozen table; `scripts/frozen_inputs.json` is byte-unchanged (`a654d66b…`) |
