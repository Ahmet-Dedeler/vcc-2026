# 252 — The shared PLOLO harness is GUARDED: 242's #4/#5/#6/#13/#14 and the silent zero now RAISE, every guard shown firing against a byte copy of the pre-252 harness, and every published number that went through it reproduces at 0.0e+00

**Date:** 2026-09-04, 17:00–17:45Z · **no `vcc` command, no build, no `cell-eval2`, no production file.**
Files edited outside this directory: exactly the two the mandate names —
`experiments/162-retrodicting-validator/plolo.py` (commit `3decf80`, sha256 `bf221b21…`, +506/−31) and
`experiments/064-pm-cd4-composition/run.py` (commit `e2f8d71`, sha256 `1caa20d4…`, +42/−3). Pre-252 byte
copies: `work/pre252/plolo.py` (`756903d1…`), `work/pre252/run064.py` (`39d0f5a2…`). Both swaps were
atomic renames (254 was building on the box and imports this harness later in its run). Peak RSS of any
step **8.32 GiB** (233's scorer); everything serial behind `gate252.py` (`memory_pressure` ≥ 40 % free,
free+inactive ≥ 20 GiB, swap ≥ 400 MiB), every step under `/usr/bin/time -l` (`work/ka-status.txt`).

---

## 0. The answer in ten lines

1. **Every arithmetic path an existing experiment used is bit-identical.** 8 known-answer runs, **97 rows,
   worst |Δ| = 0.0e+00** on every one (§2): 216's `f = 0.34/0.10` (8), 218's `compose218` (16, headline
   **+7.6803 % / +0.004887** to the digit), 246's bridge rows through `bridge250` (40, incl. the five named),
   222's `repro222` (0 failures), 249's `repro249` (ROW 2 / ROW 3 PASS, 16 rows), 233's BASE rows on
   abc/real/pess (9 rows × ~50 leaves incl. `nmae`/`fid`/`reach`/`jac`), and the k562_gwps fold **rebuilt
   from scratch** bit-identical to 162's cache on all seven arrays and every audit field.
2. **14 guards shown FIRING against the pre-252 byte copy, 12 PASS twins identical, 0 failures** over the 52
   probe rows (`probes252.py` in `blind` then `fire` mode; `compare252.py`; `work/compare252.md`).
3. **#5:** the three near-identical operators are told apart by name — `union` (6.257), `union_pruned`
   (5.895), 095's 8.036 file **unregistered → raises**; a one-ULP-tampered `held_orion_hct116` (max
   eigenvalue 8.036340, indistinguishable by eye) **raises**; `Operator(union.npz, name="union_pruned")`
   **raises SWAPPED OPERATOR**. Every fold's own operator still loads (§1.1).
4. **#4:** `E064.load_frozen_operator(tampered, expected=sha256_file(tampered))` — the exact call `plolo.py`
   made — **passed pre-252 and now raises**; `expected=""` (which silently skipped the check) now raises (§1.1).
5. **The Oct-22 trap fires at apply time:** 233's pess-panel operator applied on the abc axis (17,942
   columns both, `columns_identical: False`) **ran silently pre-252** (output `3e410ea6…`) and now raises
   `AXIS MISMATCH`; on its own axis it passes with the identical output `253b015c…` (§1.6).
6. **#13:** at 40 % coverage the pre-252 preserve substituted 1.0 and **returned the carrying rows scaled
   ×5.685** with no error; now `ON-AXIS PRESERVE UNDEFINED … 108/269 rows (40.1 % coverage)`. At 60 % both
   harnesses return the identical `764d2426…` (§1.3).
7. **#14:** on an all-zero truth block the pre-252 `endpoints()` did not read 0.0 — it **crashed with
   `ZeroDivisionError` in `sse`** (the silent 0.0 needs an all-zero prediction too); now `sgn`/`sgn_w`
   (and undefined `sse`/`gstar`) are **NaN** and `mean_endpoints()` skips them explicitly with a count (§1.4).
8. **The silent zero:** `telohaec` (58/269, 21.6 %) and `gse295214_h1` (61/269, 22.7 %) on the k562 fold
   at `f = 0.34` **returned fraction 0.0 pre-252 and now raise `SilentZeroError`** with the coverage and
   both medians in the message; `sizing="carrying"` returns 218's rule **bit-identical to
   `compose218.solve_fraction_carrying`** (0.2394513569621213 / 0.27503312984503736);
   `sizing="legacy_all_rows"` returns the zero by name, loudly. CD4 (251/269) is unchanged:
   0.37376132237758763 = 216's published fraction (§1.5).
9. **#6:** a changed source list, operator, `PM_BLEND` or gate-changing `MIN_CPM` **printed `cache hit` and
   restored the stale audit pre-252**; now each **refuses the cache** naming the disagreement, the unchanged
   inputs hit (a legacy cache only after seven recomputed checks), and new caches carry a key digest of
   every input (§1.2).
10. **Call sites whose default-path numbers changed: none.** Call sites that now raise (all measured, §3):
    diagnostic reads of the all-rows per-unit on sub-50 % sources in 233/236/238/223/218 (their applied
    fractions never came from the zero; one named line — `plolo.SIZING_DEFAULT = "legacy_all_rows"` — re-runs
    them unchanged, which is how §2's 233 rows were reproduced), and 233's deliberate trap price
    (`trapprice233.py`, needs `check_axis=False`). **HARNESS GUARDED.**

---

## 1. Item by item

Every FIRE row below is a measurement: the same probe module run against the pre-252 byte copy
(`--harness work/pre252/faux/…`, plolo `756903d1…`, E064 `39d0f5a2…`; `work/probes-blind.json`, and the same
run against the live pre-252 harness before the swap, `work/probes-blind-live-pre252.json`, identical) and
against the guarded harness (`bf221b21…` / `1caa20d4…`; `work/probes-fire.json`). Verdicts: `work/compare252.md`.

### 1.1 #5 / #4 — operator identity (`plolo.Operator`, `E064.load_frozen_operator`)

**Diff.** `plolo.py`: `OPERATOR_REGISTRY` — five operators by NAME, each pinned by npz sha256, metadata-json
sha256, ordered-symbol axis digest (`axis_sha256`, 235's canonical form), 062's fit-time index digest
(`wcols_sha256`) where the json carries one, `n_cols`, max eigenvalue; `_REGISTRY_BY_SHA`; an import-time
assert that the shas agree with `E064.KNOWN_OPERATOR_SHA256`. `Operator.__init__(path, *, name=None)` hashes
the file, resolves it by sha, checks the json sha, `n_cols` and max eigenvalue against the pin, and records
`name` / `sha256` / `axis_sha256` / `provenance`; a file not in the registry is accepted only if its own json
declares `sha256` == the file hash **and** an axis record (`axis_sha256` or `wcols_sha256` — the 233/235
convention; all 20 re-fits under 233/238/244 qualify, `work/probes-fire.log`), otherwise
`OperatorGuardError: UNREGISTERED OPERATOR`; `name=` asserts one specific operator. `apply()` verifies the
axis (§1.6) and prints one provenance line per load (rule 11). `build_fold` passes the registry pin (or the
json's declared sha) as `expected_sha256`, never `sha256_file(op_path)`.
`064/run.py`: `KNOWN_OPERATOR_SHA256` (the same five); `EXPECTED_WHITENER_SHA256` reads from it; the check
now requires a well-formed 64-hex expectation, equality, **and** membership in the table or the metadata's
own `sha256` (a caller-computed hash can no longer satisfy it).

**Independence of the pins.** The three held shas are 062's Aug-23 files (`7cecbf7f…`/`97e486e0…` were already
pinned in 064; `13d988fc…` is the prefix 162's fold audit recorded on Sep 2). The production shas are the
builder's `WHITENER_REGISTRY`. The axis digest `edf6ba88…` of the three 17,942-column folds equals the
builder's `LEGACY_OPERATOR_AXIS_SHA256` value — a cross-implementation check that passed — and the k562
fold's `6eaba848…` was measured over the fold whose index digest equals the json's fit-time
`wcols_sha256` `a000353b…` (verified for all three held operators, `work/probes-fire.log`).

| probe (`work/probes-*.json`) | pre-252 byte copy | 252 |
|---|---|---|
| `Operator(095/out/whitener_union.npz)` max eig 8.036353 | accepted, no name | **raises** UNREGISTERED (sha `5cf3899d…`) |
| `Operator(data/donors/whitener_union.npz)` 6.256702 | accepted, no name | accepted as **`union`** |
| `Operator(data/donors/whitener_union_pruned.npz)` 5.894570 | accepted, no name | accepted as **`union_pruned`** |
| `Operator(TAMPERED held_orion_hct116)` — one float32 ULP moved (`mk_tampered252.py`, sha `5f54b31c…`), max eig 8.036340 | accepted, no name | **raises** UNREGISTERED |
| `Operator(union.npz, name="union_pruned")` | not expressible (TypeError) | **raises** SWAPPED OPERATOR: hashes to `5608e807…` (= `union`), name pins `e2f65ecf…` |
| `Operator(TAMPERED, name="held_orion_hct116")` | not expressible | **raises** SWAPPED OPERATOR |
| `Operator(union_pruned.npz, name="union_pruned")` | not expressible | accepted |
| `Operator(FOLDS[held].op)` ×4 | accepted | accepted: `held_k562_gwps` / `held_orion_hct116` / `held_orion_hek293t` / `union_pruned`, max eigs identical |
| `load_frozen_operator(TAMPERED, expected=sha256_file(TAMPERED))` — plolo's own call form | **passed** (n_wcols 17942, evecs (17942, 17942)) | **raises** "neither pinned … nor declared by its own metadata json" |
| `load_frozen_operator(genuine, expected="")` | **passed** (empty expectation skipped the check) | **raises** "an empty expectation is not a check" |
| `load_frozen_operator(genuine, expected=sha256_file(genuine))` / `(genuine, pinned 7cecbf7f…)` | passed | passed, identical |
| `load_frozen_operator(TAMPERED, pinned 7cecbf7f…)` | raised (the old check works when handed a pin) | raises |

Every existing `Operator(` call site among the harness's callers (`work/operator-call-sites.txt`) passes a
FOLDS operator, a production operator, or a self-described re-fit; all still load.

### 1.2 #6 — `build_fold`'s cache

**Diff.** `fold_inputs(held)` collects every input the cached arrays depend on (held, sources, op path/kind,
op-json sha, `PM_BLEND`, `MIN_CPM`, PerturbMatch and donor paths, panel and official-axis symbol digests,
targets); `fold_cache_key` = sha256 of its sorted json; `fold_cache_status(held)` decides `hit` (keyed cache,
key equal) / `legacy` (a pre-252 cache — no key — accepted **only** after `audit.sources`, `audit.op`, the
registry sha prefix, `PM_BLEND == 0.5` (the value every pre-252 cache was built under), `wcols`, `gate@MIN_CPM`
and `scored` recompute identically from the live inputs; 0.03 s) / `miss` (the first disagreement named) /
`absent`. `build_fold` uses it, prints the decision, writes `cache_key` + `inputs` into new caches, never
rewrites a legacy one, and returns `audit["cache"]`, `op_name`, `cache_key` alongside every pre-existing key.

| probe | pre-252 | 252 |
|---|---|---|
| unchanged inputs | `cache hit` | `cache hit (legacy: audit.sources, audit.op, audit.op_sha256==registry, PM_BLEND==0.5, wcols, gate@MIN_CPM, scored)` |
| `FOLDS[k562_gwps]["sources"] = ("orion_hek293t",)` | **`cache hit`, audit still says `[orion_hct116, orion_hek293t]`** | **refused** (`REBUILD ATTEMPTED` sentinel reached); `fold_cache_status`: `legacy cache sources [...] != live ['orion_hek293t']` |
| `["op"] = held_orion_hek293t.npz` | `cache hit`, stale audit op | **refused**: `legacy cache op … held_k562_gwps.npz != live … held_orion_hek293t.npz` |
| `PM_BLEND = 0.25` | `cache hit` | **refused**: `built at PM_BLEND 0.5; live PM_BLEND 0.25` |
| `MIN_CPM = 20` (the smallest value that changes this fold's gate) | `cache hit` | **refused**: `gate != recomputed at MIN_CPM 20.0` |
| `MIN_CPM = 5` | `cache hit` | still valid — and correctly so: every gene the k562 bank carries is above 5 CPM, the gate is identical (`gate_identical_at_5: true`) |

**The real rebuild** (`ka_cache252.py` §C, `work/ka_cache.json`): the k562_gwps fold rebuilt into an empty
dir at `limit_gib=10` in **4.8 s, peak 2.83 GiB**; `pre`, `truth`, `b_ctrl`, `gate`, `scored`, `wcols`,
`targets` **bit-identical** to 162's Sep-2 cache, every audit field identical; the new file carries the key;
a second call is a keyed `hit`; with the source list changed the keyed cache misses naming
`['donor_paths', 'sources']`. (An earlier unplanned rebuild during the staging smoke test, through a faux
root, gave the same identity: `work/probe/rebuild-identity-faux.json`.) The four caches copied by
216/218/222/246/249 are byte-identical to 162's (shas in `work/`), so every one validates as `legacy`.

### 1.3 #13 — the on-axis preserve (`Operator.apply`)

**Diff.** `before = float(np.median(norm(out[:, wcols])))` — the `or 1.0` is gone; `if not before > 0` raises
`OperatorGuardError: ON-AXIS PRESERVE UNDEFINED … {carrying}/{n} rows ({coverage} coverage)`; same for `after`.
Arithmetic where the median is positive is the same expression (`t * (before / after)`).

| synthetic residual on the k562 fold (rows zeroed by `default_rng(252)`) | pre-252 | 252 |
|---|---|---|
| **40 %** (108/269 rows) | silent: `before_raw 0.0`, output median 0.0, **carrying rows scaled ×5.685** (`08cfe0bc…`) | **raises** `… only 108/269 rows (40.1% coverage) carry anything on it` |
| 60 % (161/269) | `764d2426…`, implied preserve 0.9298 | identical `764d2426…` |
| full `pre` | `apply(full pre)` | identical digest (and every §2 row) |

### 1.4 #14 — `endpoints()` sign credit

**Diff.** `n_live = live.sum()`; if 0, `sgn = sgn_w = NaN`; `sse = NaN` when `a > 0` and `d == 0` (was a
`ZeroDivisionError`), unchanged otherwise (`1 - x²/(a d)` if `a > 0` else 1.0); `gstar` unchanged. New
`mean_endpoints(rows, keys=None)` averages endpoint dicts skipping NaN **explicitly** and returns
`_skipped[key]` / `_n`. Output schema unchanged (no new keys), so every `**e` row layout is preserved.

| | pre-252 | 252 |
|---|---|---|
| all-zero truth, live prediction | **`ZeroDivisionError: float division by zero`** (in `sse`, before the sign line) | `sgn NaN, sgn_w NaN, sse NaN, pds 0.5, cos 0.0` |
| all-zero truth **and** prediction (the silent case, by the old formula `agree.sum()/max(0,1)`) | `sgn 0.0, sgn_w 0.0, sse 1.0` — "every sign wrong" | `sgn NaN, sgn_w NaN, sse 1.0` (smoke test, `work/new/`) |
| `mean_endpoints([real, zero, both-zero])` | — | `sgn 0.5` from the one live row, `_skipped {'sgn': 2, 'sgn_w': 2, 'sse': 1}` |
| real truth (k562 fold) | `sgn 0.496000 sgn_w 0.499716` | identical to every digit |

### 1.5 The silent zero (`solve_fraction`)

**Diff.** `solve_fraction(whitened, cd4, target_f, *, sizing=None)`; `SIZING_RULES = ("all_rows",
"carrying", "legacy_all_rows")`, `SIZING_DEFAULT = "all_rows"`. Same `added`, same medians, same
`per_unit = rn / bn if bn else 0.0`, same `target_f / per_unit` when `per_unit > 0`. Where it is not:
`"all_rows"` **raises `SilentZeroError`** (empty residual, zero base, or — the case — a present residual with
an all-rows median of zero, the message carrying carrying-rows/coverage/both medians and what
`sizing="carrying"` would solve to); `"carrying"` uses `median(‖added_i‖ : ‖added_i‖ > 0)` (218 §7.1, the
builder's `solve_carrying_row_fraction`, 247's October rule) and raises only on an empty residual;
`"legacy_all_rows"` returns `(0.0, per_unit)` printing `[plolo] LEGACY SIZING …`. `sizing=None` reads the
module-level `SIZING_DEFAULT`, the one-line named opt-in for an archive driver.

| k562_gwps fold, `f = 0.34`, whitened at α 0.982 / γ 0.65 | rows | pre-252 (fraction, per_unit) | 252 default | 252 `sizing="carrying"` | `compose218.solve_fraction_carrying` |
|---|---|---|---|---|---|
| **`telohaec`** | 58/269 (21.6 %), median all 0.0, carrying 6.9737 | **(0.0, 0.0)** | **`SilentZeroError`** | (0.2394513569621213, 1.4199126048543733) | identical |
| **`gse295214_h1`** | 61/269 (22.7 %), carrying 6.0715 | **(0.0, 0.0)** | **`SilentZeroError`** | (0.27503312984503736, 1.2362147069030087) | identical |
| CD4 (`cd4_for`) | 251/269 | (0.37376132237758763, 0.9096714390808992) | identical — and = 216's published `fraction` 0.373761322377588 | (0.3673…) | — |
| `sizing="legacy_all_rows"` on the two sources | | | (0.0, 0.0) with the LEGACY line | | |

### 1.6 The Oct-22 trap at apply time (not in the mandate's list; it fell out of #5)

`apply()` verifies the axis it is applied on against the operator's fit record: symbol digest
(`axis_sha256`, 235) or index digest (062's `wcols_sha256`; 233's `axis_sha256` is an index digest — the
pess operator's recorded `d66261d7…` equals the digest of the axis reconstructed from the two source banks
and 199's `def-bundle` panel, `work/probes-fire.json`). Cached per axis; a bare object with no record (216's `RefitOperator`,
which borrows `Operator.apply`) is not checked. `check_axis=False` is the named bypass.

| 233's `pess__held_orion_hct116` (self-described) | pre-252 | 252 |
|---|---|---|
| applied on the **abc** fold's 17,942 columns (`columns_identical: False`) | **ran silently**, output `3e410ea6…` — the −0.0429 trap | **raises `AXIS MISMATCH`** |
| applied on the pess axis | `253b015c…` | identical `253b015c…`, `axis_checked = indices:axis_sha256` |
| `check_axis=False` on the abc axis | — | `3e410ea6…`, bit-identical to the pre-252 silent output (the reproduction path for `trapprice233.py`) |

---

## 2. Regression table — every known-answer row, at 0.0e+00

Each script calls the harness exactly as the original does (its own `score_effects`, constants, fold caches
— byte copies of 162's in `work/plolo/`, or the experiment's own — and the same sequence of `build_fold` /
`Operator` / `apply` / `solve_fraction` / `add_post_whiten_residual`). `worst |Δ|` is the max over every
compared field (endpoints, fractions, per-unit, preserve scalar).

| row | script | what | rows | worst \|Δ\| | peak RSS | s |
|---|---|---|---:|---:|---:|---:|
| 4 | `ka216_252.py` vs `216/work/sweep216.json` | `f0.34`, `f0.1` × 4 folds: `cos_agg`, `pds`, `cos`, `a15_cos`, `a15_pds`, `sgn_w`, `sse`, `fraction`, `f_per_unit`, `preserve_scale` | 8 | **0.0e+00** | 7.22 GiB | 16 |
| 1a | `ka218_252.py` vs `218/work/compose218.json` + `diag218.json` | `cd40.34+none0`, `+h1PAIRCARRY0.28`, `+h1PAIRCARRY_SHUF0.28`, `cd40.1+none0` × 4 folds; headline 3 cancer folds | 16 + 2 | **0.0e+00** | 7.98 GiB | 25 |
| 1b | `ka246_252.py` = `bridge250.stage246/whiten_with/arms` part (a) vs 246's json | `STACK`, `B2_cd4_svd5`, **`B2_cd4_svd5_resolved`**, `B2_cd4_svd{3,10}_resolved`, **`B2_cd4_randE5`**, **`B2_cd4_svd5_SHUFDEF`**, `B3_cd4_svd5_splithalf`, **`C2_partial@0.34`**, `E_cd4_svd5+partial@0.34` × 4 folds | 40 | **0.0e+00** | 4.82–5.08 GiB | 12–16 each |
| 2 | `repro222_copy.py` (byte copy of 222's, `091c223c…`, run from here so its json lands here) | KA-3a macro **+0.076802815** / pds **+0.004886738** (|Δ| 0.00e+00), KA-3c ≤ 0.018 pp, KA-3d shuffled strictly negative; 0 failures | 4 + 12 | **0.0e+00** | 8.26 GiB | 59 |
| 2 | `repro249_copy.py` (byte copy, `7e5722b9…`) × 4 folds + `--summary` | KA-A/B/C/D per fold all `max|d| 0.00e+00`; ROW 2 PASS (arm rel_macro −0.008086748638359342, d_pds +0.012262757492797335 = 246's json); ROW 3 PASS (random control d_pds −0.0014 vs +0.0123) | 16 | **0.0e+00** | 4.93–5.17 GiB | 7–9 each |
| 3 | `ka233_252.py` = `score233.main()` `--panels abc real pess --arms BASE`, with `plolo.SIZING_DEFAULT = "legacy_all_rows"` named | every leaf of each BASE row (49–51 per row: `cos_agg`, `pds`, `cos`, `sgn_w`, `nmae`, `fid`, `reach`, `jac`, sizing frame…) vs `233/work/score233.json` | 9 | **0.0e+00** | 8.32 GiB | 60 |
| 6 | `ka_cache252.py` §C | k562_gwps fold rebuilt from scratch vs 162's cache: 7 arrays + 12 audit fields | 19 | **identical** | 2.83 GiB | 5 |

Headline numbers as reproduced (all equal to the published values to the last digit):

| published number | reproduced |
|---|---|
| 218 `cd40.34+h1PAIRCARRY0.28`, 3 cancer folds, `rel_macro_nokolf` | +0.07680281456197664 (**+7.6803 %**) |
| … `d_pds_nokolf` | +0.004886737506532087 (**+0.004887**) |
| 218 shuffled control at 0.28 | −0.014217197335142632 / −0.005405719499964083 |
| 216 k562 `f0.34` `cos_agg` / `pds` / fraction | 0.034566731076184604 / 0.7958811888512827 / 0.37376132237758763 |
| 216 kolf2 `f0.1` `cos_agg` / `pds` | 0.007814272764670947 / 0.6101280255692021 |
| 233 BASE abc / real / pess (k562) `cos_agg` | 0.03698316970799328 / 0.039059241210066833 / 0.03865781081691754 |
| 233 BASE abc / real / pess (k562) `pds` | 0.8071861880190127 / 0.7139390561129689 / 0.6992841163310962 |
| 249/246 `B2_cd4_svd5_resolved` solved fractions, 4 folds | 0.395879400489 / 0.394269198143 / 0.397869055772 / 0.336447740404 |

`repro222`/`repro249` imported the builder as it sat on disk (`3615d20f…`, 254's in-progress edit with its
flag unset; both scripts passed, as 254's own byte-inertness row requires).

---

## 3. Call sites whose behaviour changed

**On the default paths of every run above: none** (§2 is the measurement). What now raises, by static read of
all 173 referencing files plus a scan of the published JSONs for rows that carried the zero
(`work/probes-*.json`, the scan in the session log):

| call site | what it did | 252 | evidence |
|---|---|---|---|
| `233/score233.py:build_final` (audit reads `plolo.solve_fraction(whitened, h1, 1.0)`), `233/oct233.py:132`, `236/stack236.py:163`, `238/six238.py:179` | read the all-rows per-unit of a sub-50 % source **as a diagnostic** next to the carrying-row solve that was actually applied | the diagnostic read raises under the default; `plolo.SIZING_DEFAULT = "legacy_all_rows"` (one line in a driver) re-runs the script unchanged — §2 row 3 is exactly that | `score233.json` 81/81 rows `h1_per_unit_all_rows == 0.0` (abc included: the union covers 77/300); `oct233.json` 160/160; `stack236.json` 894/900; `six238.json` 340/439 `per_unit_all_rows == 0`; `cd4_per_unit_all_rows > 0` on every panel (0.67–0.94) |
| `218/resid218.py:376` | sized `telohaec` / `gse295214_h1` through the all-rows rule and applied the 0.0 | those 96 arms (+4 `f0`) now raise instead of scoring the baseline again | `resid218.json`: 100/476 rows `audit.fraction == 0.0`, coverage 21.6–22.7 %; 218 §7.1 already classed them vacuous (endpoints bit-identical to the baseline) |
| `223/phase0_repro.py` | the 21 % / 40 % rows exist to show the zero | raise at the first masked row | `phase0_repro.json` 2/4 rows `solved_fraction == 0.0` |
| `233/trapprice233.py:45` | applied the abc operator on the pess/real axis **on purpose** | raises `AXIS MISMATCH`; `apply(..., check_axis=False)` reproduces it bit-for-bit (§1.6) | — |
| any `Operator(path)` on 095/096's `whitener_union.npz`, `whitener_union_full`, `whitener_full`, 099's `whitener_hek293t` | accepted anything of the right shape | unregistered and not self-described → raise | no PLOLO caller passes one (all 60 `Operator(` sites pass FOLDS ops, the two production ops, or 233/238's self-described re-fits) |
| `E064.load_frozen_operator(..., expected_sha256="")` / an unpinned file | skipped / tautology | raises | callers 064, 069, 070, 085 pass pinned tables; 104, 149, 178 pass `sha256_file(path)` on registered files — still load |
| `endpoints()` on an all-zero truth block | `ZeroDivisionError` (or 0.0 with a zero prediction) | NaN | no caller found with such a block |

**Published rows computed through the zero path** (the mandate expected none among 223's 21; these are the
later ones): resid218's 96 `telohaec`/`gse295214_h1` arms (vacuous by 218's own §7.1), 223's two phase-0
demonstration rows, and the **diagnostic** `per_unit_all_rows` fields in 233/236/238 — no applied fraction in
any of them came from the zero. 0 of the §2 numbers depended on it.

---

## 4. What was NOT done, and the residual risks

- **`E062.add_post_whiten_residual`** still carries `float(np.median(...)) or 1.0` (`062/run.py`); 062 is not
  in this mandate's write set. A merged residual whose on-axis median is zero would still get scale
  `1.0/candidate` silently there. Every path into it here goes through `solve_fraction` first, which now
  raises before an absent residual reaches it.
- **Legacy caches are accepted by recomputation, not by key.** Seven checks (§1.2) cover every array except
  `pre`/`truth`, which depend on the sources (checked) and `PM_BLEND` (asserted 0.5). A pre-252 cache built
  under a patched `PM_BLEND` other than 0.5 cannot exist in the archive (no script sets it; grep).
- **`plolo.SIZING_DEFAULT` is a module global.** It is the explicit, named, printed opt-in the mandate asked
  for ("keep the old path reachable under an explicit argument that names it"); it exists so archive scripts
  that read the zero as a diagnostic can be re-run without editing them. New code should pass `sizing=`.
- `Operator()` now costs one sha256 of a 1.3 GB file (~0.3–1.3 s per load, page-cached) and `build_fold`'s
  legacy validation 0.03 s; visible in every `[operator]` line.
- 254 was running throughout (its staging build; it imports this harness through `lib246` for its
  known-answer row 2). Its default-path numbers are covered by §2 row 1b (the same `lib246`/`bridge250`
  construction); its log carried no guard error at the time of writing.
- The 1.3 GB tampered operator is regenerable (`mk_tampered252.py`, seconds) and was deleted after the runs;
  its record is `work/probe/tampered.json`.

---

## 5. Files

| file | what |
|---|---|
| `probes252.py` / `compare252.py` | the guard probes (`--mode blind|fire`, `--harness`) and the verdict table → `work/probes-{blind,fire}.json`, `work/compare252.{md,json}` |
| `mk_tampered252.py` | the swapped operator (one float32 ULP) → `work/probe/tampered.json` |
| `ka216_252.py`, `ka218_252.py`, `ka246_252.py`, `ka233_252.py`, `ka_cache252.py`, `repro222_copy.py`, `repro249_copy.py`, `run_ka252.sh` | the known-answer suite → `work/ka*.json`, `work/repro222.json`, `work/repro249-*.json`, `work/ka-status.txt`, `work/summary252.json` |
| `gate252.py` | the memory gate |
| `work/pre252/` | byte copies of the pre-252 `plolo.py` / `run064.py` and the faux root the blind probes import them from |
| `work/plolo/` | byte copies of 162's four fold caches (gitignored) |
| commits | `e2f8d71` (064 loader), `3decf80` (harness); message bodies in `work/commit-*.txt` |

**HARNESS GUARDED.**
