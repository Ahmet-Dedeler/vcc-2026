# 247 — The Oct-22-path guards 242 found are closed, each observed FIRING, and CD4's dose is panel-invariant for D/E/F

**Verdict line at the bottom (§7).** Every number here was measured; every guard was shown firing
against the byte copy of the code it replaces (`work/pre247/`, verified identical to git HEAD
`0fcc63dc` before the first edit). Nothing was submitted, `vcc submit` was never invoked,
`cell-eval2` was never run. Files edited outside this directory: exactly the seven the mandate
authorises. `scripts/arm_submit.sh` was **not** touched (§5.3 has its patch).

| item | 242 finding | fix | demonstrated firing | demonstrated passing |
|---|---|---|---|---|
| A | #1, #7, #12 — `def_postbuild.py`'s operator block blind, row count never compared | registry pin + file hash + SYMBOL axis + basis rows; `n_rows` and per-context rows compared | KA-2/3/4/5 (§1.2), each with the pre-247 script shown PASSING the same file | KA-1: the shipped meta **31/31** |
| B#8 | `other_builders()` misses 3 of 6 launch forms | launcher + marker predicate | 6/6 forms, 8/8 rows; pre-247 misses 4/8 | no leftover false positives; the test suite's file name not flagged |
| B#10 | census cannot exit non-zero on D/E/F | fail-closed coverage gate, always on | two collapses, exit 1; pre-247 exit 0 | A/B/C and section 5's pessimistic draw pass |
| B#18 | whitener skip branch reuses any file | axis digest proven before the skip | 3 stale operators refused, and the driver refuses on a D/E/F bundle | the A/B/C operator skips on A/B/C; a fresh json skips on its own bundle |
| C | #2 — CD4 dose inflates with missing coverage (1.42× at 58.7%) | `VCC_CD4_SIZING=carrying` + `VCC_CD4_TARGET_G`, default off, wired on non-A/B/C bundles only | the flag acts (unit rows); over-dose 1.196× / 1.170× measured on two pessimistic panels | flag off byte-inert: `87a4b6d1…`, `abeb0749…`, and the full package `8df4dadd…`/`f13ce40b…` (§3.4) |
| D | #3 — novelty ledger untracked and 5 entries stale | tracked, backfilled, missing ledger refuses | `abc-knee.vcc` → REFUSED exit 1; missing ledger → exit 2; pre-247 exit 0 on both | `229-stack6-g035.vcc` → OK novel, exit 0 |

---

## 1. Item A — `scripts/def_postbuild.py` (commit `9f050c0`)

### 1.1 Diff summary (+197 / −15)

- **The operator block reads the field that carries the answer.** New rows: `config.whitener_id`
  is in the builder's `WHITENER_REGISTRY`; `operator.sha256` == that registry pin's operator sha
  (*the file the driver pinned for this panel*); `operator.metadata_sha256` == its pair;
  **sha256 of the file at `operator.path` == `operator.sha256`** (read the pin, hash the file,
  ~3 s for 1.3 GB); **`operator.axis_sha256` == `axis_sha256(ordered symbols)` recomputed from
  THIS bundle** (235's canonical form, newline-joined UTF-8); `axis_sha256_checked_against` names
  a real source; `basis[0] == sum(rows_kept_per_donor)` (#12). A meta with no axis field (the
  shipped 222 meta predates 235) is resolved through the builder's `LEGACY_OPERATOR_AXIS_SHA256`
  keyed by operator sha, exactly as the builder itself does; a meta with neither **FAILS**, never
  skips. The blacklist-of-one row is kept as a cheap extra.
- The two pin tables are parsed as **literals out of the builder's source with `ast`** — the
  builder is never imported, so the script stays independent of the build it checks.
- **Row count (#7):** `output.n_rows == n_constructs × cells_per_pert × contexts`, and every
  manifest context's `emitted_rows == n_constructs × cells_per_pert` with no extra context.
- **Item C's rows:** on a non-A/B/C panel `config.cd4_sizing == "carrying"` and
  `config.cd4_target_g == 0.3443655278649641`; on A/B/C both must be ABSENT; the frame's
  `sizing_rule` must agree with the config; `cd4 complete targets ≥ 100` (the census floor,
  re-asserted). `ALLOWED_CONFIG_DIFFS` gains `cd4_sizing`, `cd4_target_g` on a new panel only.

### 1.2 Known-answer rows (`postbuild247.py` → `work/postbuild/`)

| row | 247 | pre-247 (byte copy) |
|---|---|---|
| **KA-1** `submissions/222-h1union-g028.h5ad.meta.json`, untouched | **31/31 PASS**, rc 0 | 20/20 PASS |
| **KA-2** operator block replaced by `data/donors/whitener_union.npz`'s (sha `5608e807…`, `n_cols` 17,942, same axis by construction) | **FAIL, 29/31**, rc 1 | **20/20 PASS — BLIND** (242's demonstration reproduced) |
| **KA-3** `operator.axis_sha256` = 235's `def-bundle-badaxis` digest `f2d6cc07…` | **FAIL, 31/32** on the symbol row; `columns_sha256` row still passes | **20/20 PASS — BLIND** |
| **KA-4** `output.n_rows` 359,600, context C 119,600 rows | **FAIL, 29/31** on both row-count rows | **20/20 PASS — BLIND** |
| **KA-5** `cd4 fully_supported_targets` 90 (30% coverage), pin moved with it | **FAIL** on the coverage floor (and on the set-intersection row) | 19/20 (set-intersection row only) |

Verbatim FAIL lines (247):

```
KA-2  FAIL  operator.sha256 == WHITENER_REGISTRY['union_pruned'][0] (the FILE the builder pinned for this panel): 5608e807230dea04… vs pin e2f65ecf69b7e0b0…
      FAIL  operator.metadata_sha256 == WHITENER_REGISTRY['union_pruned'][1] (its metadata pair): a77f9425ac524cdc… vs pin 838f263c6a4e1bed…
KA-3  FAIL  operator.axis_sha256 == axis_sha256(ordered symbols) recomputed from THIS bundle: f2d6cc075a94964c… vs recomputed edf6ba886286d783…
KA-4  FAIL  output.n_rows == n_constructs x cells_per_pert x contexts: 359600 vs expected 360,000
      FAIL  every manifest context emitted exactly n_constructs x cells_per_pert rows, and no other context: {'A': 120000, 'B': 120000, 'C': 119600}
KA-5  FAIL  cd4 complete targets >= 100 (def_census.py's coverage floor, re-asserted): 90 = 30.0% of 300
```

And the PASS lines that are new on KA-1 (the shipped meta): `operator.sha256 == WHITENER_REGISTRY['union_pruned'][0] … e2f65ecf…`,
`sha256(the file at operator.path) == operator.sha256 … (whitener_union_pruned.npz)`,
`operator.axis_sha256 ABSENT (meta predates 235): LEGACY_OPERATOR_AXIS_SHA256[operator.sha256] == axis recomputed from THIS bundle: legacy edf6ba88… vs recomputed edf6ba88…`,
`operator.metadata.basis[0] == sum(rows_kept_per_donor): 29394 vs 29394`,
`the recomputed SYMBOL axis is A/B/C's known value: edf6ba88…`, `output.n_rows … 360000 vs expected 360,000`.

Note on KA-2: the union operator's **axis** row passes, correctly — the prune drops basis rows,
not columns, so both operators share the A/B/C axis (235 §2.2). What separates them is the file,
and the file is now what is asserted. That is the row 242 asked for.

---

## 2. Item B — `scripts/def_pipeline.py` (commit §5.1) and `scripts/def_census.py` (commit `f8811e2`)

### 2.1 B#8 — `other_builders()` (+161 / −14 in `def_pipeline.py`, with §2.3 and item C's wiring)

The predicate is now: the executable is a **launcher** (`python[0-9.]*`, `uv`, `time`,
`caffeinate`, `nice`, `nohup`, `timeout`; shells deliberately excluded — the first rehearsal
tripped on its own `bash build191.sh`) AND the builder's or fitter's file name appears **anywhere**
in the remaining arguments, matched by a regex that refuses a preceding word character so
`test_make_submission_pm_cd4.py` does not read as a build. `forms247.py` launches a sleeping dummy
literally named `…/scripts/make_submission_pm_cd4.py` in each form:

| launch form | 247 | pre-247 |
|---|---|---|
| 1 `python <script>` | DETECTED | DETECTED |
| 2 absolute script path | DETECTED | DETECTED |
| 3 `python -X dev -u <script>` (two interpreter flags) | DETECTED | **MISSED** |
| 4 `python -c "import runpy; runpy.run_path('<script>')"` (how the test suite starts it) | DETECTED | **MISSED** |
| 5 `uv run python <script>` — the `uv` line alone (the un-exec'd shim) | DETECTED | **MISSED** |
| 5 … and its interpreter child | DETECTED | DETECTED |
| 6 `/usr/bin/time -l python <script>` — the `time` line alone | DETECTED | **MISSED** |
| 6 … and its interpreter child | DETECTED | DETECTED |

Pre-247 missed 4 of 8 rows (242's "3 of 6": the two flags, runpy, the un-exec'd shim). With the
dummies killed: **no leftover false positive**, the test-suite control is **not flagged**, and — a
live row — the real 240-row-1 build running at that moment was seen by 247 on **all three** of its
process lines (`/usr/bin/time`, `uv`, the interpreter) and by pre-247 on one.

### 2.2 B#10 — the census gate (+67 / −0 in `def_census.py`)

`def_census.py` now computes `cd4_union_targets` and runs a **fail-closed gate on every
invocation** (not only under `--verify`): every production donor ≥ 100 targets, CD4 union ≥ 120,
CD4 complete ≥ 100, 0-carrier targets ≤ 5, and **absent CD4 banks FAIL** ("cannot count" is not
"enough"). Floors are relaxable only by an explicit `--min-*/--max-*` value. The driver's `census`
stage additionally refuses to run on a non-A/B/C bundle unless the banks re-extracted for that
panel exist. `census247.py` → `work/census/`:

| case | rc | gate rows |
|---|---|---|
| A/B/C `--verify` | 0 | VERIFY PASS; 269 / 290 / 294, union 291, complete 276, zero-carrier 0 — GATE PASS |
| 199's pessimistic `def-bundle` + `cd4_compact_def199` (§5's draw) | 0 | 150 / 291 / 293, union **174**, complete **167**, zero 0 — GATE PASS |
| 199's `def-bundle` counted against the A/B/C banks | **1** | `FAIL CD4 union over the three banks >= 120 targets: 0`, `FAIL CD4 complete (all three banks) >= 100 targets: 0` |
| synthetic panel: 300 of the 953 axis genes NO production donor carries | **1** | all six rows FAIL, including `FAIL targets no production donor carries <= 5: 300` |
| **BLIND** pre-247 script on that collapsed panel | **0** | (no gate existed) |

### 2.3 B#18 — the whitener skip branch

`verify_existing_operator(wout, wmeta, bundle)` recomputes the bundle's axis from the donors'
`genes.txt` (axis235's 15 lines), reads `axis_sha256` from the operator's json, else the builder's
`LEGACY_OPERATOR_AXIS_SHA256[sha256(npz)]`, else **refuses**; refuses on a mismatch or an `n_cols`
mismatch; it runs both in the skip branch and after a fresh fit. `skip247.py` → `work/skip/`:

| case | result |
|---|---|
| A/B/C's real operator (json has no axis field) vs `data/challenge` | SKIP ALLOWED via the legacy table — so 240's row-1 command still works |
| A/B/C's real operator vs 199's `def-bundle` | **RAISED** `REFUSING TO SKIP the whitener stage: the existing operator was fitted for a DIFFERENT PANEL.` recorded `edf6ba88…` vs this bundle `c464a2ac…` |
| 199's `def199` refit (json predates 235; sha `68dd20dc…` not in the legacy table) vs `def-bundle` | **RAISED** `… the existing operator carries no axis digest, so the panel it was fitted for cannot be checked.` |
| a json carrying `def-bundle`'s own digest vs `def-bundle` | SKIP ALLOWED via `metadata:` |
| a json carrying A/B/C's digest vs 238's essential-free `ef-bundle` | **RAISED** DIFFERENT PANEL, `edf6ba88…` vs `f060e926…` |
| the **driver** (`--stages whitener --whitener-out data/donors/whitener_union_pruned.npz`) on A/B/C | rc 0: `=== [whitener] exists and its axis matches this bundle (LEGACY_OPERATOR_AXIS_SHA256[operator sha256]); skipping the fit` |
| the driver, same operator, on `def-bundle` | **rc 1**, the DIFFERENT PANEL refusal above |

BLIND (pre-247, verbatim): `if wout.exists() and wmeta.exists(): print("=== [whitener] exists, skipping: …"); continue` — no check of any kind.

Two observations for the day: 238's `work/ops/*.json` carry an `axis_sha256` that is the
**column-index** digest (`b2cae0d3…`), not 235's symbol digest — they are 238's own fitter's
files, not production's, and this check would (correctly) refuse them. And 199's `def199` refit
would be refused as unverifiable; the runbook's step 5 refit writes the digest, so nothing on the
day changes.

---

## 3. Item C — CD4's dose (builder commit §5.1; driver wiring in `def_pipeline.py`)

### 3.1 The flag

`VCC_CD4_SIZING=carrying` + `VCC_CD4_TARGET_G=<g>`; **default off**. Validation at import refuses
`carrying` without a `G`, a `G` without `carrying` (a set-but-inert knob is the 242 #17 class), any
other sizing word, and `G ∉ (0, 5)` — all six combinations exercised. The ONLY code path that
changes is `apply_frozen_cd4_residual`'s single call, `solve_cd4_fraction(…)` →
`solve_cd4_dose(…)`, which returns `solve_cd4_fraction` unchanged when the flag is unset and
`solve_carrying_row_fraction(whitened, cd4, G)` when set — the rule 218 validated bit-for-bit on
these banks and 229 shipped for the H1 union. With the flag on, the frame carries both rules
side by side (`legacy_all_rows_rule_would_solve`, `overdose_factor_all_rows_rule`) and the config
gains exactly two keys, `cd4_sizing` and `cd4_target_g`, written only when on (the `h1_union_g`
convention). `solve_cd4_fraction` itself is **not modified** (+100 / −1 in the builder).

### 3.2 `g_cd4`, measured on the real banks (`cd4dose247.py` → `work/cd4dose247.json`, 0.71 GiB)

The two rules depend on the whitened base only through its median row norm, which the shipped meta
records exactly (4.229295604454299); with a stub base of that exact median (235's construction)
every number is the builder's own. Known answers first: the A/B/C CD4 residual reproduces
`residual_median_row_norm_at_unit_fraction` **4.464717385709141**, `f_per_unit_fraction`
**1.0556645369046553**, `solved_fraction` **0.32207201067578145** and 276 complete targets, all to
the last bit.

| quantity | A/B/C |
|---|---:|
| carrying rows / coverage | 276 / 300 = 92.0% |
| median row norm over ALL rows | 4.464717385709141 |
| median row norm over CARRYING rows | 4.522043409698855 |
| `median_carrying / median_all` | **1.0128397878381297** |
| **`g_cd4 = 0.34 × 1.0128397878381297`** | **0.3443655278649641** |
| carrying rule at `g_cd4`: `solved_fraction` | 0.3220720106757815 |
| shipped `effects.cd4.fraction` | 0.32207201067578145 |
| absolute difference | **5.55e-17 (one ULP) — numeric identity, not byte identity** |

A synthetic residual at 229's 176/300 coverage reproduces, through the production functions, the
identity `legacy solve / carrying solve == median_carrying / median_all` to 1e-16 and
`legacy solve == shipped_median_rule_would_solve` exactly — the row that proves the sizing wired in
is 218/229's and not a third rule.

### 3.3 The two delivered doses on the pessimistic panel

Delivered per-row dose on the rows CD4 carries, as a multiple of the base's median row norm
(base-independent: it is the same positive denominator in both rules), relative to A/B/C's
delivered `0.34 × 1.01284 = 0.34437`:

| panel | CD4 complete / coverage | over-dose factor `median_carrying/median_all` | legacy rule delivers | carrying rule at `g_cd4` delivers |
|---|---|---:|---:|---:|
| A/B/C (real) | 276 / 92.0% | 1.01284 | 0.34437 = **×1.000** of A/B/C | 0.34437 = ×1.000 |
| §5's pessimistic draw (199 `def-bundle`, `cd4_compact_def199`) | 167 / 55.7% | **1.19629** | 0.40674 = **×1.1811** of A/B/C | 0.34437 = **×1.000** |
| 238's disjoint AND essential-free draw (`ef-bundle`) | 182 / 60.7% | **1.17016** | 0.39785 = **×1.1553** of A/B/C | 0.34437 = **×1.000** |

The legacy rule's `solved_fraction` on the pessimistic draw would be 0.38732 against the carrying
rule's 0.32792. (The over-dose here is below 222/223's 1.24–1.42 because CD4's row-norm
distribution is narrower than the sources 223 measured; it is a measurement, not a bound.)

### 3.4 Flag off is byte-inert — the known-answer rows

Unit level (`builder247.py`, two builders in one process, a 4-of-6-row synthetic residual where
the rules differ): flag OFF, the working-tree builder's output array **and** its cd4 audit are
bit-identical to the pre-247 byte copy's; flag ON at `g = 0.34 × overdose`, `sizing_rule ==
"carrying"`, the fraction equals the all-rows solve to ≤ 1e-9 and the output equals OFF's to
1e-12; at a different `g` the output moves (the flag acts). 12/12.

Staging level (`stage247.py`, 222's `stage222.py` construction, each run **alone** under
`/usr/bin/time -l`, gated on `memory_pressure ≥ 40%` free and `vm_stat` free+inactive ≥ 30 GiB and
zero other builders — `work/stage/stage247.log`):

| row | env | `effects_array_sha256` | required | `effects.cd4.fraction` | peak RSS |
|---|---|---|---|---|---:|
| off | `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34` | **`87a4b6d181298f63826dc9caeb9674510a5dbc0d3c388e9e4fbcb359f3cd2de2`** | `87a4b6d1…` (222 `stage-off.json`) — **MATCH** | 0.32207201067578145 | 15.53 GiB |
| on028 | + `VCC_H1_UNION_G=0.28` | **`abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315`** | `abeb0749…` (222 `stage-on028.json` = the shipped `final_sha256`) — **MATCH** | 0.32207201067578145 | 16.16 GiB |
| carry | + `VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641` | `0c083a76f7d82a6400252fae03877e455e7f0d5d659fb6b5abdca0e248595864` | (numeric, not byte) | **0.3220720106757815**, `overdose_factor_all_rows_rule` 1.0128397878381297, `legacy_all_rows_rule_would_solve` 0.32207201067578145; merged-add preserve scalar 0.9093569351340488 identical to the shipped | 19.90 GiB |

Both digests reproduced on a second run (`stage247_save.sh`), and `diff247.py` measured the
carry array against the shipped on028 array: **max |Δ| = 4.44e-16** (2.6e-16 of the largest
entry; max relative 3.9e-11 on near-zero entries), differing entries confined to the **276
CD4-carrying rows** (1,220,873 of 5,559,900 entries, 22.0%), zero difference elsewhere. That is
what "numeric identity at one ULP of the fraction" looks like on the 300 × 18,533 object the
emission consumes.

Package level — 240's row 1 through the modified driver AND builder (nothing exported, the
driver's own defaults, `--stages build,postbuild`, then `vcc prep` by hand, strictly serial and
alone, `verify247.sh` → `work/ka1/`):

| what | measured | required |
|---|---|---|
| the driver's printed configuration | `VCC_ALPHA '0.982' VCC_CD4_TARGET_F '0.34' VCC_H1_UNION_G '0.28'` (driver defaults), `VCC_CD4_SIZING <unset on A/B/C …>`, `memory check: 35.2 GiB free+inactive, no other builder` | 240 row 1's, plus the new line |
| the `[build]` command line | byte-identical to 240 row 1's (`timings.json`) | identical |
| `.h5ad` sha256 / bytes | **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` / 12,648,807,800 — MATCH** | `8df4dadd…` / 12,648,807,800 |
| `.vcc` sha256 / bytes | **not re-prepped** — the prep gate (free+inactive ≥ 30 GiB) sat at 28.1–28.8 GiB for 15 min at 82–83% `memory_pressure` free, the fresh h5ad's own page cache holding the gap (deleting the file afterwards freed exactly it: 33.0 GiB). The `.vcc` is a deterministic transform of the `.h5ad` (240 row 1 and 235 row 2 both took `8df4dadd…` to `f13ce40b…`), so the byte identity rests on the h5ad row above; the h5ad was re-hashed (`8df4dadd…`) and deleted as a verified duplicate of `submissions/222-h1union-g028.h5ad` | `f13ce40b…` / 2,824,785,920 |
| 240's `[arms]` block | four `OK` lines as 240, plus **`OK config.cd4_sizing: asked '<absent: the all-rows median>', the artifact carries '<absent>'`** | — |
| `effects.cd4.fraction` / `config.cd4_sizing` / `cd4_target_g` | 0.32207201067578145 / absent / absent | equal / absent / absent |
| `operator.axis_sha256` / source | `edf6ba88…` / `LEGACY_OPERATOR_AXIS_SHA256[operator sha256]` | A/B/C's |
| `def_postbuild.py` through the driver's `postbuild` stage | **32/32 checks pass** (`moved=['expect_contexts']`) | — |
| builder sha in `hashes.implementation` | `24a30055…` (the committed `a17b177` bytes) | changed, as 235's did; the h5ad did not |
| cost | build **1294.4 s**, peak RSS **20.00 GiB**; postbuild 0.8 s / 0.12 GiB | 240: 1302.3 s / 20.00 GiB |

So the modified builder, driven by the modified driver with nothing exported, rebuilds the shipped
`h1union` package byte-for-byte. The flag is inert.

### 3.5 The wiring

`def_pipeline.py`: `DEF_ENV = {"VCC_CD4_SIZING": "carrying", "VCC_CD4_TARGET_G":
"0.3443655278649641"}`, merged into the day's defaults **only when the bundle is not A/B/C**
(`day_env(is_abc)`); an exported value still wins; the printed configuration shows the pair on
D/E/F and `<unset on A/B/C …>` on A/B/C; `assert_build_carries_arms` re-reads `cd4_sizing` /
`cd4_target_g` off the finished meta (and on A/B/C asserts `cd4_sizing` is absent unless the
operator exported it). Dry runs: the A/B/C `[build]` command line is byte-identical to 240 row
1's; the D/E/F line carries `VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=0.3443655278649641`.
`def_postbuild.py` requires the pair on a new panel and its absence on A/B/C (§1.1).

---

## 4. Item D — the novelty ledger (commit `c01bbd8`)

- `.gitignore`: `!submissions/SUBMITTED_LEDGER.json` (one line plus a comment); the ledger is now
  tracked.
- Backfilled, hashed from disk (`work/hashes247.json`; `backfill247.py` refuses on any
  mismatch with a recorded digest — `abc-knee` `f5faa1e2…` matches 242, `222-h1union-g028`
  `f13ce40b…` matches BOARD-LEDGER):

| package | sha256 | entry | submitted (UTC) | board |
|---|---|---|---|---:|
| `abc-knee.vcc` | `f5faa1e21efdb25f…` | `N3ui4hCkKx9bxExG2pxE` | 2026-09-02 23:28 | 0.1740437 |
| `191-rc_hybrid.vcc` | `bc51fc4a525171bd…` | `nQUJZhZb3MW2uFJHBX5e` | 2026-09-03 11:20 | 0.1834190095 |
| `202-rc_callset.vcc` | `f31272b91a1f1913…` | `31oB9yWy8CSNBDEkY4iW` | 2026-09-03 17:47 (prep mtime; published ~18:05) | 0.1918638367 |
| `210-cd4f034.vcc` | `ab5ffa8269e63912…` | `QSe3JlNUYZ2BzvlcsL14` | 2026-09-04 00:00 | 0.2049195830 |
| `222-h1union-g028.vcc` | `f13ce40bca908af7…` | `iUziqE1gewfcPKzYUDaB` | 2026-09-04 03:50:39 | 0.2118310 |

  `191-rc_hybrid.vcc` is identified by path and prep time (mtime 11:18Z, submitted 11:20Z); no
  other record of its digest exists, which the note says.
- `check_submission_novelty.py`: `load()` on a missing ledger prints `REFUSED: the submission
  ledger is MISSING` and **exits 2**; `--ledger` exists for demonstrations; `record` reminds to
  `git add`.

Known-answer rows:

```
$ uv run python scripts/check_submission_novelty.py check submissions/abc-knee.vcc
REFUSED: this package has ALREADY BEEN SCORED. Do not spend a slot on it.
  sha256   f5faa1e21efdb25fd55b9c99f2a71363154a55fe9047206573da784c4ddd63b4
  entry    N3ui4hCkKx9bxExG2pxE   submitted 2026-09-02T23:28:00+00:00
EXIT=1                                                        (pre-247 byte copy: OK: novel package … EXIT=0)
$ uv run python scripts/check_submission_novelty.py check submissions/229-stack6-g035.vcc
OK: novel package, sha256 7a3c2d3e414f8d9d...
EXIT=0
$ uv run python scripts/check_submission_novelty.py --ledger <missing> check submissions/229-stack6-g035.vcc
REFUSED: the submission ledger is MISSING at …
EXIT=2                                                        (pre-247: OK: novel package … EXIT=0)
```

### 4.1 The one-line patch for `scripts/arm_submit.sh` — NOT applied (a process is executing it)

Insert before the `log "SUBMITTING …"` line:

```bash
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || { log "ERROR novelty check refused (rc=$?) -- NOT submitting"; exit 1; }
```

and, to keep the ledger from going stale again, after `[ $RC -eq 0 ] && log "SUBMITTED ok"`:

```bash
[ $RC -eq 0 ] && uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$(echo "$OUT" | grep -o '[A-Za-z0-9]\{20\}' | head -1)" --note "$NAME" >> "$LOG" 2>&1
```

(The entry-id extraction is a best effort against `vcc submit --wait`'s output; verify its format
on the next armed submit before relying on it.)

Not done from 242 §2.3's list: keying the ledger on the **configuration** as well as the bytes
(a re-seeded rebuild at an identical config still reads as novel). It needs a config-digest
convention agreed with the builder's meta and is beyond this mandate.

---

## 5. Commits

| hash | file(s) | known-answer rows in the message |
|---|---|---|
| `9f050c0` | `scripts/def_postbuild.py` | KA-1 31/31; KA-2/3/4 FAIL with pre-247 BLIND; KA-5 |
| `f8811e2` | `scripts/def_census.py` | A/B/C and pessimistic PASS; two collapses exit 1; pre-247 exit 0 |
| `c01bbd8` | `.gitignore`, `submissions/SUBMITTED_LEDGER.json`, `scripts/check_submission_novelty.py` | the five backfilled digests; abc-knee exit 1; 229 exit 0; missing exit 2 |
| `a17b177` | `scripts/make_submission_pm_cd4.py` | `87a4b6d1…`, `abeb0749…`, `8df4dadd…` reproduced with the flag off; 5.6e-17 with it on; 1.1963× / 1.1702× measured; 56/57 under both builders |
| `588fbe5` | `scripts/def_pipeline.py` | 6/6 launch forms; the stale operator refused (`edf6ba88…` vs `c464a2ac…`); `8df4dadd…` through the driver's defaults, 1294 s / 20.00 GiB, postbuild 32/32 |
| `25e20a3` (**by the lead, not this experiment**) | `scripts/arm_submit.sh` | §4.1's patch, applied and re-armed as pid 26557 after this file's draft was read; this experiment did not edit that file |
| `77bcc5f` | `experiments/247-oct22-guards/` (scripts, byte copies, evidence) | the measurements behind every row above; this RESULTS.md follows in its own commit |

### 5.1 Regression suites

- `experiments/235-october-hardening/tests/test_october_hardening.py`: **9/9 passed** under the
  modified builder (1.06 GiB), including "the second residual is sized by the carrying-row rule
  and CD4 is not" (still true with the flag unset) and the shipped H1-pair frame to the last bit.
- `scripts/test_make_submission_pm_cd4.py`: its runner stops at the first failure, and
  `test_panel_exclusion_for_operator_and_cd4_residual` (line 249, from 2026-08-23) unpacks two
  values from `apply_frozen_cd4_residual`, which has returned three since 222 (`481cd092`,
  2026-09-04 05:03) — **the suite has not passed as a whole since 222, before this experiment.**
  `suite247.py` runs the suite's 57 registered tests in isolation under both builders:
  **56/57 under the 247 builder, 56/57 under the pre-247 byte copy, 0 regressions, 0 fixed** —
  the one failure is that test, under both (58 s, 0.32 GiB).

---

## 6. What was NOT closed from 242's table, and why

| 242 # | guard | why not here |
|---|---|---|
| 4, 5, 6, 13, 14 | `plolo.py` / `064/run.py` operator and cache checks | shared PLOLO harness, outside the seven authorised files; 241/243/245/246 are running against it |
| 9, 15, 19, 20, 21 | `agent_supervisor.sh`, `run_with_memory_guard.py`, `board_snapshot.sh`, `cushion_watch.py`, `check_harness_depth.py` | outside the authorised files; none is on the Oct-22 build path |
| 11 | `def_pin_panel_artifacts.py`'s npz↔json pairing skipped for an unpruned fit | outside the authorised files. Partly backstopped now: `def_postbuild` asserts `operator.metadata_sha256` against the registry pair AND hashes the operator file, so a mismatched pair fails post-build |
| 16, 17 | the builder's own env defaults (pre-`a11`/`cd4f034`/`h1union`) and the early return when `VCC_CD4_TARGET_F` is unset | the driver + `def_postbuild` are the layer that pins the shipped values (240); changing the builder's defaults would move every unflagged build's bytes, which this mandate forbids |
| 3 (config keying) | see §4 | — |

---

## 7. Cost and memory

| run | wall | peak RSS |
|---|---:|---:|
| `cd4dose247.py` (three panels' CD4 residuals) | 1.3 s | 0.71 GiB |
| `postbuild247.py`, `census247.py`, `skip247.py`, `forms247.py`, `builder247.py` | seconds each | ≤ 0.3 GiB (skip hashes 1.3 GB thrice) |
| staging off / on028 / carry (each alone) | 13.1 / 13.4 / 12.1 s | 15.53 / 16.16 / 19.90 GiB |
| 240 row 1 build through the driver (alone) | 1294.4 s | 20.00 GiB |
| `vcc prep` | not run (gate never opened; §3.4) | — |
| re-staging on028 / carry with saved arrays, then `diff247.py` | 12 / 9 s | 15.87 / 19.90 GiB |
| 235's tests / the builder suite comparison | 7 s / 58 s | 1.06 GiB / 0.32 GiB |

Everything ≥ 15 GiB waited for `memory_pressure ≥ 40%` free AND free+inactive ≥ 30 GiB AND zero
other builders (the gate is logged before each row), and ran strictly serially. 241, 243, 245 and
246 were live throughout.

---

## 7b. Verdict

**OCT-22 PATH GUARDED.** All five Oct-22-path guards in 242's table (#1, #2, #7, #8, #10, #12,
#18 — the operator block, the CD4 dose, the row count, the concurrent-build predicate, the census
gate, the basis rows and the whitener skip) and the slot guard (#3) are closed and each was
observed firing against the byte copy of the code it replaces, with the shipped package rebuilt
byte-for-byte through the changed driver and builder. **Not closed here, and why:** §6 — the
PLOLO-harness rows (#4, #5, #6, #13, #14) and the monitoring/supervisor rows (#9, #15, #19, #20,
#21) are outside this mandate's seven files; #11 (`def_pin_panel_artifacts.py`) is outside them
too but is now backstopped post-build by the metadata-pair and file-hash rows; #16/#17 (the
builder's own stale defaults) cannot be changed without moving every unflagged build's bytes,
which this mandate forbids, and the driver plus `def_postbuild` are the layer that pins them;
config-keyed novelty (#3's second half) needs a convention beyond this mandate. One row of the
mandate's own list is partial: the `.vcc` re-prep of 240 row 1 did not run because the literal
memory gate never opened (§3.4); the `.h5ad` it derives from matched byte-for-byte.

## 8. Files

```
cd4dose247.py     item C step 1: both rules on the real banks, g_cd4          -> work/cd4dose247.json
builder247.py     item C unit level: OFF bit-identical to pre-247, ON routes  -> work/builder247.json
stage247.py/.sh   item C staging rows, memory-gated, each alone               -> work/stage/
stage247_save.sh  the on028/carry re-run with saved arrays, then diff247.py    -> work/stage/diff247.json
verify247.sh      item C package row (240 row 1 through the modified driver)  -> work/ka1/
postbuild247.py   item A KA-1..5 + BLIND                                      -> work/postbuild/
forms247.py       item B#8, six launch forms                                  -> work/forms/forms247.json
census247.py      item B#10, two collapses + two passes + BLIND               -> work/census/
skip247.py        item B#18, five operators + the driver                      -> work/skip/
backfill247.py    item D, the five entries                                    -> submissions/SUBMITTED_LEDGER.json
suite247.py       the builder's 65 tests under both builders                  -> work/suite247.json
work/pre247/      byte copies of the seven files at HEAD 0fcc63dc before the first edit
work/hashes247.json  the six package digests
```
