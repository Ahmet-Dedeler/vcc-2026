# 244 — The Oct-22 sequence, end to end, on the fixed driver, on a panel that is disjoint AND essential-free at the pessimistic coverage draw

**Date 2026-09-04.** Deterministic. **No `vcc submit`, no `cell-eval2`, no model name, no description, nothing built
outside `experiments/244-october-rehearsal/work/`.** `scripts/def_pipeline.py` was **not** modified (the skip-branch
hole 240 named is already closed by 247, §5). No file outside this directory was left changed by this experiment: the one
shared file the day's sequence edits (`scripts/make_submission_pm_cd4.py`, step 7's pin) was byte-restored from the copy
saved at pin time (§7.1 — including the one thing that did not go as planned).

---

## 0. Verdict

```
RUNBOOK EXECUTABLE AS WRITTEN: NOT.  Executable with the section-5 branch taken explicitly.

  All eleven stages ran on 238's ef-bundle (300 targets, 0 in A/B/C, 0 essential, K562 146/300,
  CD4 union 186/300) through the current driver with NOTHING exported, and produced a package
  that carries every shipped arm plus 247's carrying-row CD4 sizing:
      operator re-fit  373df206...  axis f060e926... (== 247's recorded digest for this bundle)
      cd4 complete 182/300 (60.7%), CD4 over-dose factor 1.1701602469240082 (== 247, bit-exact)
      H1 union 28 carrying rows (9.3%), all-rows rule 0.0, carrying rule 0.20198547321112115
      effects.cd4.fraction 0.3417130399758081 at g 0.3443655278649641; cd4_target_f 0.34; h1_union_g 0.28
      def_postbuild.py 32/32; vcc prep --dry-run --contexts D,E,F clean; nnz 33.1% of cap; max 22,323 counts/cell
      stages 1,692.7 s (28.2 min); launch -> prepped artifact 29.7 min wall

  But the runbook's literal command STOPS at step 3 on this coverage (047's CD4 gate, 186 < 290),
  and its section-5 response cites a file (`target_row_map.csv`) that only the RELAXED extraction
  writes, so the order it prescribes -- confirm the mapping, then relax -- cannot be followed as
  written.  The mapping was confirmed off the source's obs column instead (§2.2), the driver was
  re-launched with `--cd4-thresholds relaxed`, and everything after that ran as the runbook says.
  Under the mandate's 30 GiB gate the prep stage refused once (29.5 GiB) and had to be re-launched;
  under the runbook's own 26 GiB default it would have run in the same invocation.  Eight runbook
  lines need changing (§8).  All known-answer rows PASS (§6).
```

---

## 1. The panel and its census

**Bundle: `experiments/238-essential-free-panel/work/ef-bundle`** (238's `efp`, panel id `vcc2026-synthetic-ef-238`,
seed 199, pool = union of the three production donors' targets on the official axis − A/B/C − the 2,021-gene
common-essential list; `pert_counts.csv` sha256 `3d093b33…`, `gene_names.csv` `25bfa667…` = the official axis, controls
byte-identical to 199's `def-bundle`, contexts D/E/F). 235's essential-free panel (`essfree235.py`, seed 0) is a target
list only, not a bundle, so 238's is the one that can be driven end to end. Census, recomputed here (`panel244.py`,
47/47 checks incl. every one of 235's 30 `carry235.json` coverage rows on the three bundles it measured):

| property | ef-bundle (this run) | A/B/C (real) | 199 pessimistic (`def-bundle`) |
|---|---:|---:|---:|
| targets / unique / ∩ A/B/C / ∩ essential | 300 / 300 / **0** / **0** | 300 / 300 / 300 / 0 | 300 / 300 / 0 / **30** |
| `k562_gwps` / `orion_hct116` / `orion_hek293t` | **146** / 293 / 296 | 269 / 290 / 294 | 150 / 291 / 293 |
| carriers 3 / 2 / 1 / 0 | 140 / 155 / **5** / 0 | 259 / 35 / 6 / 0 | 141 / 152 / 7 / 0 |
| CD4 Rest / Stim8hr / Stim48hr, union, **complete** | 183 / 185 / 184, **186**, **182** | 283 / 287 / 284, 291, 276 | 173 / 173 / 168, 174, 167 |
| `panel_support` (pm, cd4) | **133, 182** | 248, 276 | 141, 167 |
| H1 pair `gse295214_h1` / `h1_2025` / union rows | **27 / 4 / 28** (3 by both) | 64 / 20 / 77 | 34 / 3 / — |
| registry: dld1 / igvf / jiang / myllia / nadig_hepg2 / nadig_jurkat / rpe1 / k562_essential | 85 / 27 / 4 / 3 / 9 / 7 / 8 / **0** | 67 / 87 / 9 / 14 / 0 / 0 / 0 / 0 | 73 / 36 / 4 / 4 / 33 / 33 / 32 / 30 |
| operator axis (`bundle_axis`, 235's symbol digest) / n_cols | `f060e926…` / 17,942 | `edf6ba88…` / 17,942 | `c464a2ac…` / 17,942 |

The five single-carrier targets (ELMOD3, KLHDC1, RWDD2B, SMIM19, TTC37; 1.7%) are inside §5's "acceptable up to ~10%".
`def_census.py`'s fail-closed gate (247) passed all six rows on this panel (146/293/296 ≥ 100, union 186 ≥ 120, complete
182 ≥ 100, zero-carrier 0 ≤ 5). Identity on the (relabelled A/B/C) controls: **D = JURKAT, E = HELA, F = CAL33**, 8/8
variants unanimous.

---

## 2. What was run

### 2.1 Run 0 — the runbook's §2 command, literally (`work/run0-asis/`)

```
uv run python scripts/def_pipeline.py --bundle <ef-bundle> --work work/run0-asis --name run0-asis \
    --whitener-id union_pruned_def --emitter rc_callset --frozen-inputs work/run0-asis/frozen_inputs.json \
    --stages preflight,identity,cd4-extract,census
```
preflight **66/66** (2.9 s), identity 9.0 s / 2.74 GiB, CD4 source sha256 `c355f535…` verified, then **stage 3 stops the
run**: `ValueError: only 186 panel targets map across any condition; require 290` (rc 1, 0.8 s) — 047's gate, exactly
238's number. The failed extraction writes **nothing** (no `cd4_compact/`, no `target_row_map.csv`; the extractor writes
that file at `extract_cd4.py:691`, after the gate at `:652`).

### 2.2 Confirming the mapping before relaxing (`mapping244.py` → `work/mapping244.json`)

The runbook's §5 says relax "ONLY after confirming the mapping is right (`target_row_map.csv`)". That file cannot exist
yet, so the check was done on the source itself: `h5py` read of the 16.8 GB atlas's `obs` (no `X`), column
`target_contrast_gene_name` (11,526 symbols; the extractor's own key). **A/B/C maps 291/300 (283 / 287 / 284 per
condition — the known answer). ef-bundle maps 186/300 (183 / 185 / 184); the 114 unmapped targets have no case-variant
hit in the atlas** — they are absent, not mis-spelt. 0.6 s, 0.25 GiB. (The relaxed extraction's `target_row_map.csv` and
`panel_union_missing` list, written afterwards, agree.)

### 2.3 Run 1 — the sequence with §5's pre-decided response (`work/run1/`, launched 13:06:13Z)

```
uv run python scripts/def_pipeline.py --bundle experiments/238-essential-free-panel/work/ef-bundle \
    --work experiments/244-october-rehearsal/work/run1 --name sparrow-in-the-vestibule \
    --whitener-id union_pruned_def --emitter rc_callset \
    --frozen-inputs experiments/244-october-rehearsal/work/run1/frozen_inputs.json \
    --cd4-thresholds relaxed --restore-builder-after --min-free-gib 30
```
**Nothing exported.** The driver printed its defaults — `VCC_ALPHA '0.982' VCC_CD4_TARGET_F '0.34' VCC_H1_UNION_G '0.28'
VCC_CD4_SIZING 'carrying' VCC_CD4_TARGET_G '0.3443655278649641'` (all "driver default"), `VCC_UNION_SOURCES <unset: the
H1 pair>`, emitter `rc_callset`, the shipped flags. Differences from the runbook's §2 line and why: `--frozen-inputs`
(a throwaway copy; the day re-pins `scripts/frozen_inputs.json` and commits it), `--restore-builder-after` (the day
commits step 7's edit), `--cd4-thresholds relaxed` (= the extractor's `--min-targets 0 --min-union-targets 0`, §5's
response), `--min-free-gib 30` (the mandate's gate; the runbook's driver default is 26). Gate at launch: 81% free, 34.4
GiB free+inactive, swap free 1.1 GiB, 0 other builders; 249 was alive but reading.

Ran preflight → … → postbuild without intervention (rc 0 each), then **refused `prep`: "free+inactive memory 29.5 GiB
< 30.0 GiB; refusing"** (rc 1 at 13:34:17Z; the driver's own metric, free+inactive+speculative+purgeable). Because the
restore step sits after `prep`, the builder was left pinned; it was restored by hand from the pin record (§7.1), and
`--stages prep --min-free-gib 30` was re-launched at 30.7 GiB (rc 0, 13:35:52Z). Second invocation:
`work/run1/driver-prep.log`.

---

## 3. The stage table

| # | stage | wall s | peak RSS GiB | rc | runbook (199, 2026-09-03) | look-for, observed |
|---|---|---:|---:|---:|---|---|
| 1 | preflight | 2.8 | 0.31 | 0 | 4.5 / 0.3 | 66/66; contexts D,E,F; 18533 / 300 / 400 / 3×18,400 |
| 2 | identity | 9.2 | 2.75 | 0 | 14.7 / 2.7 | D JURKAT, E HELA, F CAL33, unanimous |
| — | CD4 source sha256 | 5 | — | — | 9–31 | `c355f535…` (page-cached) |
| 3 | cd4-extract (relaxed) | 1.4 | 0.32 | 0 | 2.2 / 0.3 | 183 / 185 / 184, union 186, 10,113 genes; **default gates fire** (run 0) |
| 4 | census | 0.2 | 0.10 | 0 | 0.4 | 140/155/5/0; predicted `panel_support` **133 / 182**; 0 essential; gate 6/6 |
| 5 | whitener | **360.9** | 12.33 (eigh 7.53) | 0 | 437.6 / 12.9 | 17,942 cols; basis 7784+10807+10803 = 29,394; dropped 4,393 (`dropped_in_panel` []); sha **`373df206…`** ≠ `e2f65ecf…`; **axis `f060e926…`** ≠ `edf6ba88…` |
| 6 | repin-frozen | 1.1 | 0.23 | 0 | 2.1 | A/B/C removed, D/E/F added, manifest + pert_counts moved, gene_names unchanged, `panel_support → None` |
| 7 | pin-artifacts | 0.5 | 0.05 | 0 | < 1 | ONE registry entry `union_pruned_def` (`373df206…`, json `3ca6c148…`) + 15 `cd4/…` hashes, 17 lines |
| 8 | repin-panel-support | 9.7 | **19.90** | 0 | 9.2 / 19.8 | `[operator] … MATCHES metadata:whitener_union_pruned_def.json`; pinned **{133, 182}** = step 4's prediction |
| 9 | build | **1290.1** | **20.00** | 0 | 1255.8 / 19.9 (240 on A/B/C: 1302.3 / 20.00) | `[operator] MATCHES`; `[arms]` **6 OK lines** (alpha, cd4_target_f, h1_union_g, h1_union_sources, cd4_sizing, cd4_target_g) |
| 10 | postbuild | 0.9 | 0.12 | 0 | 1 | **32/32**, `moved=['cd4_sizing','cd4_target_g','contexts','expect_contexts','whitener_id']` |
| 11 | prep (`--dry-run --contexts D,E,F`) | 15.9 | 16.25 | 0 (2nd invocation) | 11.8–13.9 / 22–24 | `✓ dry run`, cells 360000, D/E/F = 120000 each, targets verified, counts-preserved |
| | **sum of stages** | **1,692.7 = 28.2 min** | | | **1,748 = 29.1 min** | |
| | launch → prepped artifact, wall | **29.7 min** (13:06:13 → 13:35:52Z, incl. the refused gate and the by-hand restore + re-launch) | | | 29 min | |

Memory gates as observed by the driver: 32.2 GiB before step 8, 38.5 before step 9, **29.5 at step 11 (refused)**,
30.7 at the re-launch. `memory_pressure` stayed 79–81% free throughout; the dip to 26.3 GiB free+inactive during the fit
(memmap pages) recovered on its own. The builder's own peak: 19.998 GiB. **Disk: `work/run1` 13.12 GiB** (h5ad 12.62 GB,
operator 1.29 GB, banks, logs); scratch memmap removed by the fitter; `work/` total 13 GiB.

---

## 4. Verification, with numbers (`verify244.py` → `work/verify244.json`, **52/52**)

**The operator was re-fit for this panel.** `operator.sha256` **`373df20652ec1f323a99137a52415ed74f9929e847cb98f50b9906732d3cc58e`**
(≠ `e2f65ecf…`), metadata `3ca6c148…`, 17,942 columns, `axis_sha256`
**`f060e9262a753dddb28d4a3b80f1ab75d4ca37f8c7de7ef9205f5c156ae0213f`** = the digest recomputed from the bundle by the
driver's `bundle_axis` = the digest 247 recorded for this bundle (`skip247.json`, stale-5) ≠ A/B/C's `edf6ba88…`;
`axis_sha256_checked_against` = `metadata:whitener_union_pruned_def.json` (the fresh json, not the legacy table);
`panel_sha256` `dfe912d1…`; `sha256(file at operator.path)` == the meta's; `columns_sha256` `b327501c…` (= 238's
`efp` operators' index digest, as it must be — same panel). `dropped_in_panel: []`.

**The configuration the day must carry**, read off the written meta: `alpha` 0.982 · `cd4_target_f` **0.34** ·
`h1_union_g` **0.28** · `h1_union_sources` `['gse295214_h1','h1_2025']` · `emitter` `rc_callset` · `whitener_id`
`union_pruned_def` · `cd4_sizing` `carrying` · `cd4_target_g` `0.3443655278649641` · `contexts` D,E,F · `expect_contexts`
`D,E,F` · γ 0.65 · λ 0.0 · `min_cpm` 10 · `keep_top` 8000 · `keep_spread` 0.25 · `tail_spread` 1.0 · `overdispersion`
4.0 · call set 400 / 2200 / 1.9 · seeds 0 / 1000003 / 1000 / 0. The rule-10 diff against `222-h1union-g028.h5ad.meta.json`
moved exactly `['cd4_sizing','cd4_target_g','contexts','expect_contexts','whitener_id']`.

**CD4 at 60.7% coverage** (`effects.cd4`): complete targets **182**, `incomplete_targets_preserved_before_global_scale`
118, `fraction` **0.3417130399758081** = `residual_frame.solved_fraction`; frame: `sizing_rule` carrying,
`target_g_requested` 0.3443655278649641, `carrying_rows` 182, `whitened_base_median_row_norm` 4.3876562671196435,
`residual_median_row_norm_over_carrying_rows` 4.421714683828407, `…over_all_rows` 3.778725773202206,
`g_per_unit_fraction` 1.0077623256324775, `delivered_carrying_row_dose_x_base` 0.3443655278649641,
`legacy_all_rows_rule_would_solve` 0.3947899954530122 (at f 0.34), `shipped_median_rule_would_solve` 0.3998590152352451 (at g).

**The over-dose factor on this panel** — `residual_median_row_norm_over_carrying_rows ÷ …over_all_rows` =
4.421714683828407 / 3.778725773202206 = **1.1701602469240082**, read off the meta AND recomputed from the staged
residual on the banks this run re-extracted (`cd4dose244.py`, 247's construction, 14/14: the run1 banks are
byte-identical to 238's on 12/12 data files, `meta.json` differing only in provenance; A/B/C 1.0128397878381297 and
199's pessimistic 1.196293034774599 reproduce 247 bit-exactly). §5.0's planning number is ~1.24×; **CD4 measures 1.170×
here and 1.196× on 199's draw** — lower, for the reason 247 §3.3 gives (CD4's row-norm distribution is narrower than the
sources 222/223 measured). Consequence in this build: the legacy all-rows rule at f 0.34 would have solved 0.39479 against
the carrying rule's 0.34171 — **×1.1553 = 1.1702 × (0.34 / 0.34437)**, i.e. it would have delivered A/B/C's per-row dose
×1.155 on the 182 carrying rows; with 247's sizing on, the delivered carrying-row dose is exactly A/B/C's (×1.000).

**The H1-union frame at low coverage** (`effects.h1_union`): `source_coverage` {gse295214_h1: 27, h1_2025: 4} = the set
intersection; `targets_by_source_count` {0: 272, 1: 25, 2: 3}; **`carrying_rows` 28** (9.33%); residual median over
carrying rows 6.082337186246015, **over all rows 0.0**; **`shipped_median_rule_would_solve` 0.0** (the pre-222 all-rows
rule is dead below 50% coverage — measured on the real base, not a stub); `g_per_unit_fraction` 1.3862383049068872;
**`solved_fraction` = `applied_fraction` 0.20198547321112115** at g 0.28 (A/B/C: 77 rows, 1.4577895501601277, 0.19207);
`targets_outside_cd4_support` 2, `entries_outside_cd4_support` 261,313; merged-add preserve scale 0.9553649891807288.
The arm shrinks with coverage as §3's row says (28 rows against 77) and does not drop; 235's raise on a zero-coverage
source did not fire (both sources carry > 0). `carry244.py` (235's functions, stub base = the shipped 4.229295604454299)
gives the same 28 rows and 0.19470 / all-rows 0.0, and all nine registry sources that cover this panel size > 0 under the
carrying rule and exactly 0.0 under the all-rows rule; `k562_essential` covers 0 and would raise if selected.

**The emitted artifact:** `work/run1/sparrow-in-the-vestibule.h5ad` sha256
**`47614a44de9d8507e3f7f4511efca8a38c08e872f72e0782cd55451eb803a552`**, 12,623,641,976 B; `output.nnz`
**1,574,288,495 = 33.1% of 4.75e9**, re-counted off the file (h5py stream) to the same number; 360,000 rows, D/E/F 120,000
each (stored entries 536,911,140 / 511,996,301 / 525,381,054; median control library 20,109 / 19,946 / 20,034);
**max counts per cell 22,323** (min 17,695) against the 1,000,000 cap; 0 negative, 0 non-integer, 0 explicit zeros.
`def_postbuild.py` **32/32** (`work/run1/postbuild.json`; the count is 247's, 240's "20/20" predates the operator, row-count
and CD4-sizing rows). `vcc prep --dry-run --contexts D,E,F` clean (§3 row 11). `whitening.after_all_median_norm`
4.3876562671196435 (A/B/C 4.229295604454299; reported, panel-dependent). PerturbMatch replacements 133 =
|panel ∩ k562 ∩ pm|. `peak_rss_gib` 19.998.

---

## 5. The skip-branch hole (240 §7 item 1) — CLOSED by 247, observed here

**The line that closes it:** `scripts/def_pipeline.py:517-523` —
`if wout.exists() and wmeta.exists(): source = verify_existing_operator(wout, wmeta, bundle)` — with
`verify_existing_operator` at `:297-336` (recomputes the bundle's symbol axis from the donors' `genes.txt`, reads
`axis_sha256` from the operator's json, else the builder's `LEGACY_OPERATOR_AXIS_SHA256[sha256(npz)]`, else refuses;
refuses on a mismatch or an `n_cols` mismatch; also run after a fresh fit at `:531`). Committed in `588fbe5` (247 B#18).
So nothing was added to the driver here and 240's row 1 was not required; 247 §3.4 re-ran it through this driver
(`8df4dadd…` MATCH, 1294.4 s). Observed in this run, on this bundle (`work/skip/`):

| case | result |
|---|---|
| `--whitener-out data/donors/whitener_union_pruned.npz` (the A/B/C production operator, json without an axis field) | **rc 1** `REFUSING TO SKIP the whitener stage: the existing operator was fitted for a DIFFERENT PANEL.` recorded (`LEGACY_OPERATOR_AXIS_SHA256[operator sha256]`) `edf6ba88…` vs this bundle `f060e926…`, columns 17942 "(not the check)" |
| `--whitener-out experiments/199-def-runbook/work/whitener_union_pruned_def199.npz` (a refit whose json predates 235) | **rc 1** `… the existing operator carries no axis digest, so the panel it was fitted for cannot be checked.` |
| `--whitener-out work/run1/whitener_union_pruned_def.npz` (the operator this run fitted) | **rc 0** `=== [whitener] exists and its axis matches this bundle (metadata:whitener_union_pruned_def.json); skipping the fit` |

Both refusals happen before any work (0.7 s), which is what 240 and 235 asked for. No diff to report.

---

## 6. Known-answer rows

| row | required | measured | verdict |
|---|---|---|---|
| **235's guard fires on a wrong-panel operator, 4 ways** (`fire244.py` = `fire235.py` with output paths only; `work/fire/`) | `ROW 5 PASS` | PASS_abc loads (`edf6ba88…` via the legacy table); FIRE_def_pessimistic, FIRE_unregistered_operator, FIRE_badaxis all raise; BLIND pre-235 builder accepts the wrong panel silently — **5/5, ROW 5 PASS**, 8.5 s, 16.52 GiB (`ru_maxrss` 17,742,643,200 B; run alone, gate 73% free / 42.2 GiB) | **PASS** |
| **the coverage table reproduces 235's rows for the same bundles** (`panel244.py`, `carry244.py`) | every `carry235.json` row | 30/30 coverage cells on abc / def199 pessimistic / def199 realistic; 7/7 fields of the shipped `h1union` frame bit-exact (77 rows, 6.1654229367116375, 1.4577895501601277, 0.19207161964444322, 0.0, 0.0); every `sized` row 235 recorded on the pessimistic bundle (solved / all-rows / f_per_unit / carrying rows) — **47 known-answer checks, 65/65 overall** | **PASS** |
| 240's row 1 | only if `def_pipeline.py` was touched | not touched (sha `ba0d6533…` before and after); 247 §3.4 holds the last run through this driver | n/a |
| 247's CD4 dose rows (`cd4dose244.py`) | 1.0128397878381297 / 1.196293034774599 / 1.1701602469240082; 276 / 167 / 182 | all bit-exact; the run1 banks give 1.1701602469240082 and are byte-identical to 238's on 12/12 data files | **PASS** |
| 238's bundle and coverage (`panel244.py`) | `pert_counts.csv` `3d093b33…`, `gene_names.csv` `25bfa667…`, 146 / 293 / 296, CD4 183/185/184 union 186 (from `cd4_238.json`) | identical; the extraction reproduced 183/185/184/186 and 10,113 genes | **PASS** |
| 247's ef-bundle axis (`skip247.json` stale-5) | `f060e926…` | the fitter wrote it, the builder matched it, the driver recomputed it | **PASS** |
| the arms guard (240 row 4's control) | no `NO` line | 6 `OK` lines, 0 `NO` | **PASS** |

---

## 7. Two things that did not go as planned, both process, both on the record

### 7.1 A second agent edited the builder under the run (rule-11 hazard, `AGENTS.md` §11.1b)

The byte copy of `scripts/make_submission_pm_cd4.py` taken at 13:02Z was HEAD's (`24a30055…`, commit `a17b177`). Step 7's
pin record (`work/run1/pin-artifacts.json`, 13:12:36Z) shows `builder_sha256_before` **`2d63b979…`** — experiment 249
(`cd4svd5`, whose mandate lets it edit that file) had added its `VCC_CD4_DEFLATE_K` code in the ten minutes between
(`work/run1/builder-249-edit-before-pin.diff`, 172 lines, all additions plus a guarded call site). The flag is read from
the environment, is `None` when unset, and the only new call is behind `if CD4_DEFLATE_K is not None`, so **this build's
behaviour is HEAD's**: the meta carries neither `effects.cd4.deflation` nor `config.cd4_deflate_k` (asserted), and every
number above matches 247's independently-derived values. What it did move: `hashes.implementation` records the pinned
builder as `fbdcbaf2…` (= 249's bytes + the 17-line pin), not `24a30055…` + pin — the rehearsal's recorded builder hash is
not reproducible from any commit. At the end the builder was restored to the pre-pin bytes (`2d63b979…`, 249's working
copy, verified `identical`) so 249's edit was preserved; the builder never changed between the pin and the restore
(hash watched at every poll), so nothing of 249's was lost. `git status --short scripts/` shows only 249's uncommitted
edit. **On the day this must be a rule: step 7 is applied to a CLEAN builder** (§8, edit 3).

### 7.2 The prep gate refused once, and the driver's rehearsal restore is not a `finally`

`prep` was refused at 29.5 GiB under the mandate's 30 GiB gate (the runbook's driver default is 26), and because
`--restore-builder-after` runs only after a successful `prep`, the driver exited with the builder still pinned. It was
restored by hand with `def_pin_panel_artifacts.py --restore --record work/run1/pin-artifacts.json` (which restores the
saved byte copy, never git HEAD — 199's own lesson) and `prep` re-run 95 s later at 30.7 GiB. On the day the pin is
committed rather than restored, so the exposure is rehearsals only; but §7's symptom table should say what to do (§8, edit 6).
247 saw the same gate sit below 30 GiB for 15 min after a build (`247 §3.4`); here it recovered in under two minutes.

---

## 8. RUNBOOK EXECUTABLE AS WRITTEN: **NOT** — the lines that must change (`docs/DEF-RUNBOOK.md` not edited here)

1. **§2 step 3 (line 129) and §5 "CD4 gates fire" (lines 411–413).** On a pessimistic essential-free panel 047's gates
   fire (186 < 290) and **the driver stops** — the §2 command is two invocations on such a day, not one. Replace
   "run step 3 with `--min-targets 0 --min-union-targets 0` ONLY after confirming the mapping is right
   (`target_row_map.csv`)" with: *"confirm the mapping off the source's `obs['target_contrast_gene_name']` (no `X`;
   `experiments/244-october-rehearsal/mapping244.py`: A/B/C must read 291 = 283/287/284, and every unmapped target must
   have no case-variant hit), then re-launch the driver with `--cd4-thresholds relaxed`; `target_row_map.csv` is written
   only by the extraction that succeeds, so it cannot be the pre-check."*
2. **§2 the ⚠️ box (lines 139–147).** "One gap the driver still has, and 235 could not close it" is false since `588fbe5`:
   replace with the closing line (`def_pipeline.py:517-523`, `verify_existing_operator` `:297-336`) and the three observed
   outcomes in §5 above.
3. **§2 step 7 (line 133).** Add: *"`git status --short scripts/make_submission_pm_cd4.py` must be EMPTY before the pin.
   Another agent's uncommitted edit is pinned with it, the meta's `hashes.implementation` becomes unreproducible, and a
   rehearsal restore can clobber it (244 §7.1)."*
4. **§2 step 9 (line 135)** "four `OK` lines" → **six** on D/E/F (`cd4_sizing 'carrying'`, `cd4_target_g
   0.3443655278649641` — 247); the §2 box (lines 87–106) should list `VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641`
   as the day's two extra defaults, off A/B/C only.
5. **§2 step 10 (line 136)** "**20/20 PASS**" → **32/32** on a D/E/F build (31/32 on A/B/C — 247 §1.2), `moved`
   = `['cd4_sizing','cd4_target_g','contexts','expect_contexts','whitener_id']`.
6. **§7 (line 461)** "the driver refuses a heavy stage … wait": add *"if the refusal is at `prep` in a rehearsal the builder
   is still pinned — `def_pin_panel_artifacts.py --restore --record <work>/pin-artifacts.json`, then `--stages prep`"*;
   and note the gate reading that refused (29.5 GiB, 40 s after a 20 GiB build).
7. **§2 line 149** "1,748 s ≈ 29 min" → add 244's measurement: **1,692.7 s of stages (whitener 360.9, build 1290.1, prep 15.9),
   29.7 min launch-to-prepped-artifact including the §5 re-launch**; step 11's "22–24 GiB" → 16.25 GiB measured here
   (235 recorded 27.11: plan for 25+).
8. **§5.-1 (lines 325–354) and §5.0 (356–368).** The rehearsal is no longer stale on the panel half: record that 244 ran the
   H1-union path and the carrying-row rule at the pessimistic essential-free draw (28 carrying rows, all-rows rule 0.0,
   solved 0.20199; CD4 182/300, over-dose **1.170×** measured, 1.196× on 199's draw, against the 1.24× planning number —
   and that with 247's sizing it is not delivered). Add the ef-bundle column to §5's table (line 389–397):
   K562 146, Orion 293/296, carriers 140/155/5/0, CD4 183/185/184 union 186, `panel_support` **133 / 182**, essential **0**.
   Also §1 line 63: the builder hash there (`ab8aa3e`, `6cc9150…`) is three commits stale — HEAD `a17b177` is `24a30055…`.

Everything else in §2–§4 ran as written and printed what the runbook says to look for.

---

## 9. Files

```
panel244.py     the panel, its census, 235/238/247 known-answer rows            -> work/panel244.json (47/47)
mapping244.py   section-5's mapping confirmation off the atlas obs               -> work/mapping244.json
carry244.py     235's carry235 functions on ef-bundle + its rows re-measured     -> work/carry244.json, .log (65/65, 47 KA)
cd4dose244.py   247's cd4dose measure() on the staged residual + the meta        -> work/cd4dose244.json, .log (14/14)
fire244.py      235's fire235.py, output paths only                              -> work/fire/fire244.json, .log (ROW 5 PASS)
verify244.py    the artifact, off the written meta and the h5ad                  -> work/verify244.json, .log (52/52)
work/run0-asis/ the literal command, stopped by 047's gate (driver.log, logs/, timings.json)
work/run1/      the rehearsal: driver.log (main invocation), driver-prep.log, logs/<stage>.log, timings.json,
                timings.main-invocation.json, gate.txt, gate-prep.txt, census/identity/preflight/postbuild.json,
                pin-artifacts.json + .builder-before.py, builder-249-edit-before-pin.diff, frozen_inputs.json,
                whitener_union_pruned_def.{npz,json}, cd4_compact/, sparrow-in-the-vestibule.h5ad (+ .meta.json)
work/skip/      the three skip-branch cases (driver-abcop, driver-def199op, driver-fresh; .log each)
work/pre244/    byte copies of the seven shared scripts at HEAD 86eb077 before anything ran, with sha256.txt
```
The 12.6 GB artifact and the 1.3 GB operator are gitignored (`*.h5ad`, `*.npz`) and reproducible (deterministic build,
6-minute fit); they are kept for re-verification and may be deleted.
