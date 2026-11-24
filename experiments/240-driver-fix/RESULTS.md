# 240 — The Oct-22 driver is **FIXED**: its defaults now rebuild the shipped package byte-for-byte, an operator override wins, and the new guard is observed firing

**Date:** 2026-09-04 · deterministic · **no `vcc submit`, no `cell-eval2`, no model name, no
description, no slot proposed.** Three files changed outside this experiment directory —
`scripts/def_pipeline.py`, `scripts/def_postbuild.py`, `docs/DEF-RUNBOOK.md` — all three
authorised by the mandate. **No shipped number moved and nothing was re-tuned: this is a wiring
fix.** `scripts/make_submission_pm_cd4.py` is untouched — still `715f2b46…` at git `b540148`
(235's state), which is row 5.

---

## 0. Verdict

```
FIXED.

  The driver's defaults, on the A/B/C bundle, with NOTHING exported, rebuild the shipped
  `h1union` package BYTE-FOR-BYTE:
      .h5ad  8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6  12,648,807,800 B
      .vcc   f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e   2,824,785,920 B
  and `def_postbuild.py` passes 20/20 on it against the CURRENT incumbent.

  Before 240 the same command built `rc_callset`'s configuration -- f = 0.10, no H1 union,
  -0.0199671 overall -- and `def_postbuild.py` printed 14/14 PASS on it (row 3, verbatim).

  All 5 mandate known-answer rows PASS.
```

---

## 1. The five known-answer rows

| # | row | required | measured | verdict |
|---|---|---|---|---|
| **1** | **driver defaults reproduce the shipped `h1union`** | `.vcc` `f13ce40b…`, 2,824,785,920 B | **`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`, 2,824,785,920 B**, from an `.h5ad` that is itself `8df4dadd…` on 12,648,807,800 B; config identical to the shipped package **in every field but `expect_contexts`**, and all eight named quantities equal (`effects.cd4.fraction` **0.32207201067578145**, `target_f_requested` **0.34**, `f_per_unit_fraction` **1.0556645369046553**, `fully_supported_targets` **276**, `operator.sha256` `e2f65ecf…`, γ **0.65**, λ **0.0**) | **PASS** |
| **2** | **an operator override WINS** | `VCC_CD4_TARGET_F=0.10` ⇒ the builder reads `f = 0.10` | **the builder reads `f = 0.1`** through the fixed driver; through the **pre-240** driver an operator exporting `0.34` was silently handed **`0.1`** | **PASS** |
| **3** | **the old behaviour, reproduced once** | `f` forced to 0.1, `VCC_H1_UNION_G` absent, `def_postbuild.py` **PASSING** | the pre-240 driver hands the builder **`f = 0.1`, `g = None`** with nothing exported, and the pre-240 `def_postbuild.py` prints **`14/14 checks pass`** on that configuration — quoted verbatim in §4 | **PASS** |
| **4** | **the new assertion FIRES** | a raise on a config missing either arm | **`REFUSING THE BUILD: the artifact does not carry the configuration it was built for`** on `rc_callset`'s config (§5, verbatim), **no raise** on the shipped `h1union` config; `def_postbuild.py` **16/20 FAIL** on the same config, and **still 16/20 FAIL when `--incumbent` is re-pointed at `rc_callset`** | **PASS** |
| **5** | **`solve_cd4_fraction` untouched** | not in the diff; `effects.cd4.fraction` unmoved | `git diff scripts/` contains **0** occurrences of `solve_cd4_fraction`; `scripts/make_submission_pm_cd4.py` is unmodified (`715f2b46…`); the driver-built artifact reports `effects.cd4.fraction` **0.32207201067578145**, equal to the shipped package's | **PASS** |

---

## 2. The defect, in one line each

| | before 240 | after 240 |
|---|---|---|
| the env | `INCUMBENT_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.1"}`, pinned by its own comment to `build191.sh` | `SHIPPED_ENV = {"VCC_ALPHA": "0.982", "VCC_CD4_TARGET_F": "0.34", "VCC_H1_UNION_G": "0.28"}`, pinned to `h1union` / `build222.sh` |
| how it was applied | `env={**os.environ, **(env or {})}` — **override**, the operator loses | `resolve_env` — **default**, the operator wins, and the resolved values are what gets printed, logged and passed |
| the H1 union | never set; would have passed through only on top of `f = 0.10`, a combination never built or scored | set by default, with `VCC_UNION_SOURCES` deliberately unset = the H1 pair (229's row 2) |
| the check | `--incumbent-meta` defaulted to `202-rc_callset` — the config the driver itself produced, so the diff was **empty** | defaults to `222-h1union`; the five arms are asserted **by value** in the tracked script, and the reference is itself checked for carrying them |
| after the build | nothing looked at the artifact | the driver re-reads the arms **off the written meta** and raises otherwise |

**Cost of the defect on the one entry that can win the prize: −0.0199671 overall**
(`cd4f034` +0.0130557 + `h1union` +0.0069114), i.e. board 0.2118310 → 0.1918638, against a gap to
first of +0.01248.

---

## 3. Row 1 in full — the row that proves the fix

```bash
uv run python scripts/def_pipeline.py --bundle data/challenge \
    --work experiments/240-driver-fix/work/ka1 --name 240-ka-h1union --stages build \
    --whitener-out data/donors/whitener_union_pruned.npz --whitener-id union_pruned \
    --cd4-out data/external/cd4t/cloud-compact/cd4_compact --frozen-inputs scripts/frozen_inputs.json
# then, by hand, exactly as runpreps229.sh packaged the shipped arm:
vcc prep …/240-ka-h1union.h5ad -g data/challenge/gene_names.csv \
    --perts data/challenge/pert_counts.csv -o …/240-ka-h1union.vcc
```

**Nothing was exported.** The whitener, the CD4 banks and the frozen table are the production
ones, so the only thing under test is the configuration the driver chooses. The driver printed:

```
the configuration this run will build (driver defaults; an exported value wins):
  VCC_ALPHA            '0.982'    (driver default)
  VCC_CD4_TARGET_F     '0.34'     (driver default)
  VCC_H1_UNION_G       '0.28'     (driver default)
  VCC_UNION_SOURCES    <unset: the H1 pair ['gse295214_h1', 'h1_2025']>
  emitter              'rc_callset'
memory check: 36.2 GiB free+inactive, no other builder
```

and at the end of the build, 240's guard, reading the written meta:

```
[arms] the built config against the configuration asked for:
  OK  config.alpha: asked 0.982, the artifact carries 0.982
  OK  config.cd4_target_f: asked 0.34, the artifact carries 0.34
  OK  config.h1_union_g: asked 0.28, the artifact carries 0.28
  OK  config.h1_union_sources: asked ['gse295214_h1', 'h1_2025'], the artifact carries ['gse295214_h1', 'h1_2025']
```

| what | measured | required |
|---|---|---|
| `.h5ad` sha256 / bytes | `8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6` / 12,648,807,800 | identical (222, reproduced by 229 row 2 and 235 row 2) |
| `.vcc` sha256 / bytes | **`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e` / 2,824,785,920** | **identical** |
| config diff vs `222-h1union-g028.h5ad.meta.json` | **`['expect_contexts']` only** (`'<absent>' → 'A,B,C'`) | meta-only, and the runbook passes it on purpose |
| `effects.cd4.fraction` | 0.32207201067578145 | equal |
| `effects.cd4.residual_frame.target_f_requested` | 0.34 | equal |
| `effects.cd4.residual_frame.f_per_unit_fraction` | 1.0556645369046553 | equal |
| `effects.cd4.fully_supported_targets` | 276 | equal |
| `effects.h1_union.carrying_targets` / `applied_fraction` | 77 / 0.19207161964444322 | equal |
| `operator.sha256` | `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123` | equal |
| γ / λ / α / emitter | 0.65 / 0.0 / 0.982 / `rc_callset` | equal |
| `def_postbuild.py` through the driver's `postbuild` stage | **20/20 checks pass**, `moved=['expect_contexts']` | — |

**Cost:** build **1302.3 s**, peak RSS **20.00 GiB**; `vcc prep` **26 s**, peak RSS **19.72 GiB**;
postbuild 0.2 s / 0.10 GiB. The two heavy stages ran **strictly serially, alone**, with 36.2 GiB
free+inactive at the start and no other builder — they are the artifact-producing stages and
cannot be run under 10 GiB; everything else in this experiment peaked at **0.20 GiB**.
(Artifacts deleted after hashing: they are byte-identical duplicates of
`submissions/222-h1union-g028.*`. The meta, the logs and `work/verify240.json` are kept.)

---

## 4. Row 3, verbatim — the pre-240 driver passing while both arms are missing

**What the pre-240 driver hands the builder**, measured by running its own `Runner` with its own
`INCUMBENT_ENV` and letting the real builder report what it parsed (`rows240.py`, byte copies from
git HEAD verified identical):

```
pre240-asdocumented              operator={}                                              -> builder read alpha=0.982 f=0.1 g=None
pre240-operator-tries-both-arms  operator={'VCC_CD4_TARGET_F':'0.34','VCC_H1_UNION_G':'0.28'} -> builder read alpha=0.982 f=0.1 g=0.28
```

The second line is the trap 239 named: an operator who *knows* about both arms and exports both
still gets `f = 0.1`, and `g = 0.28` then rides on top of a dose combination that has never been
built or scored.

**And the pre-240 `def_postbuild.py` on that configuration** — `submissions/202-rc_callset.h5ad.meta.json`
is the config the pre-240 driver produces, with the pre-240 driver's own arguments
(`--incumbent submissions/202-rc_callset.h5ad.meta.json --allow-config-diff whitener_id`):

```
meta      ~/Code/virtual-cell-challenge-2026/submissions/202-rc_callset.h5ad.meta.json
incumbent ~/Code/virtual-cell-challenge-2026/submissions/202-rc_callset.h5ad.meta.json
bundle    ~/Code/virtual-cell-challenge-2026/data/challenge  panel vcc2026-val-1  contexts ['A', 'B', 'C']
operator  e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123  columns 17942
whitening after_all_median_norm 4.229295604454299 (A/B/C incumbent 4.229295604454299; reported, NOT asserted — it is a panel-dependent quantity)

  PASS  config.contexts == the manifest's: ['A', 'B', 'C'] vs ['A', 'B', 'C']
  PASS  config.targets == n_constructs: 300
  PASS  operator.sha256 is NOT the A/B/C operator (it is panel-dependent; rebuild it): e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123
  PASS  operator.columns == the axis recomputed from THIS panel: 17942 vs 17942
  PASS  operator.columns_sha256 == the axis recomputed from THIS panel: b2cae0d381c3cd70… vs b2cae0d381c3cd70…
  PASS  operator.metadata.n_cols == operator.columns: 17942
  PASS  all_panel_coordinates_excluded: 
  PASS  perturbmatch replacement count == set-intersection prediction: 248 vs 248
  PASS  cd4 complete targets == set-intersection prediction: 276 vs 276
  PASS  cd4 complete targets > 0 (the banks were re-extracted for this panel): 276
  PASS  frozen_bundle.panel_support is pinned and matches: {"cd4_complete_targets": 276, "perturbmatch_replacement_targets": 248}
  PASS  config differences vs 202-rc_callset.h5ad.meta.json are only ['contexts', 'expect_contexts', 'targets', 'whitener_id']: moved=[]
  PASS  nnz under the density cap: 1,578,026,567 = 33.2% of cap
  PASS  row count == n_constructs x cells_per_pert x contexts: expected 360,000 rows

14/14 checks pass
```

**Fourteen of fourteen, on a package worth −0.0199671.** `moved=[]` is the whole story: the check
diffed the artifact against the incumbent the driver was pinned to, so the two agreed by
construction.

> **One caveat, stated because it was measured and not assumed.** This row runs the pre-240 check
> on the A/B/C `rc_callset` meta, not on a D/E/F artifact — building a D/E/F package needs the
> October bundle. On the day the config would additionally carry `contexts`, `whitener_id` and
> `expect_contexts`, all three of which the pre-240 command explicitly allows. **The arm keys are
> the only ones at issue and they are identical in both cases**, so the `moved=[]` above becomes
> `moved=['contexts','whitener_id','expect_contexts']` on the day — still all allowed, still PASS.

---

## 5. Row 2 and row 4 — the override wins, and the guard fires

**Row 2, through the fixed driver** (same measurement method: the builder reports what it parsed):

| operator exported | driver | builder read |
|---|---|---|
| *(nothing)* | pre-240 | α 0.982, **f 0.1**, g **None** |
| `VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28` | pre-240 | α 0.982, **f 0.1** ← silently reset, g 0.28 |
| *(nothing)* | **240** | α 0.982, **f 0.34**, **g 0.28**, sources `['gse295214_h1','h1_2025']` |
| **`VCC_CD4_TARGET_F=0.10`** | **240** | α 0.982, **f 0.1** ← **the operator's value wins**, g 0.28 |
| `VCC_H1_UNION_G=""` | **240** | α 0.982, f 0.34, g **None** — an explicitly empty value is honoured as "arm off", which is how `build222.sh`'s own known-answer row turns it off |

**Row 4, the raise, verbatim** (`assert_build_carries_arms` on `rc_callset`'s config, under the
driver's resolved intent):

```
REFUSING THE BUILD: the artifact does not carry the configuration it was built for.
  OK  config.alpha: asked 0.982, the artifact carries 0.982
  NO  config.cd4_target_f: asked 0.34, the artifact carries 0.1
  NO  config.h1_union_g: asked 0.28, the artifact carries '<absent>'
  NO  config.h1_union_sources: asked ['gse295214_h1', 'h1_2025'], the artifact carries '<absent>'
  meta        …/submissions/202-rc_callset.h5ad.meta.json
  asked for   {'VCC_ALPHA': '0.982', 'VCC_CD4_TARGET_F': '0.34', 'VCC_H1_UNION_G': '0.28'}
  A dropped arm is silent everywhere else -- the artifact builds, `vcc prep` passes,
  and the board is dark on the day. What each arm was measured to be worth:
    VCC_ALPHA: the a11->a15 amplitude series; 211 closed it at 0.982
    VCC_CD4_TARGET_F: +0.0130557 (`cd4f034`)
    VCC_H1_UNION_G: +0.0069114 (`h1union`)
  Do not submit this package. DEF-RUNBOOK section 2 step 9.
```

**The control matters as much as the raise:** the same call on the shipped `h1union` meta does
**not** raise, and the driver's real build (§3) printed four `OK` lines — so the guard is
discriminating, not merely loud.

**And the second layer, `def_postbuild.py`, on the same armless config:**

| reference passed | result |
|---|---|
| default (`222-h1union`) | **FAIL, 16/20** — `arm cd4_target_f`, `arm h1_union_g`, `arm h1_union_sources`, and the rule-10 diff |
| re-pointed at `202-rc_callset` (the pre-240 default) | **FAIL, 16/20** — the three arm checks still fire, **and** *"the `--incumbent` reference IS the current incumbent"* fails with `MISSING ['cd4_target_f','h1_union_g','h1_union_sources']` |
| the shipped `h1union` config (control) | **PASS, 20/20** |

The middle row is the point: **the check can no longer be defeated by choosing an older
reference**, which is exactly how it passed before.

---

## 6. The diff

| file | change |
|---|---|
| `scripts/def_pipeline.py` | `INCUMBENT_ENV` → **`SHIPPED_ENV`** at `h1union`'s values, with `SHIPPED_UNION_SOURCES` and a board-delta table; new **`resolve_env`** (defaults, operator wins, reports which keys the operator supplied); `Runner.run` takes `env_defaults`, resolves them, and prints/logs the **resolved** values plus an `[env] operator-supplied, and it WINS` line; new **`assert_build_carries_arms`**, called after the build; a startup banner of the resolved configuration with an **ARM DRIFT** warning naming the board deltas at risk; `--incumbent-meta` default → `222-h1union`; new `--allow-arm-drift`, forwarded to the postbuild stage; `INCUMBENT_FLAGS` → `SHIPPED_FLAGS` (values unchanged — the flags never moved, only the env did) |
| `scripts/def_postbuild.py` | new **`SHIPPED_ARMS`** table (config key → value, env name, board delta, entry) and `arm_rows`; `--incumbent` default → `222-h1union` with a clear error if it is missing; **five arm checks by value** plus a check that **the reference itself carries the arms**; `--allow-arm-drift`; a loud epilogue when an arm check fails |
| `docs/DEF-RUNBOOK.md` | §2 opens with a box naming the two arms, their board deltas and entries, and the rule that the box and `SHIPPED_ENV` move together; §1's incumbent rows → `h1union` (with `rc_callset` explicitly marked **not** a valid `--incumbent`); §2 steps 8/9/10 → the shipped env, `build222.sh` provenance, the `[arms]` block to look for, the new `--incumbent` and **20/20**; §3's `VCC_CD4_TARGET_F` row → **0.34** with its board evidence, plus **new rows for `VCC_H1_UNION_G` and `VCC_UNION_SOURCES`**; §3.1 records that its own pre-decided rule fired, and adds an `h1union` row and a driver row; §4 gains three known-answer rows (the shipped-config build, **the driver's defaults reproducing it**, and the guard firing); §5.-1 records that the **driver half** of the staleness is closed while the **panel half** is not |

`scripts/make_submission_pm_cd4.py`, `scripts/build_union_whitener.py`, `data/`, `src/` and every
other experiment: **untouched.**

---

## 7. What this does NOT fix — read before the day

1. **The `whitener` stage's skip branch is still unguarded.** `def_pipeline.py` skips the operator
   fit if `<work>/whitener_<id>.npz` and its json already exist, with no check on either, so a
   stale operator left by an earlier run is silently reused. The runbook already names the
   one-line fix (read `axis_sha256` from the json and compare it against the axis recomputed from
   this bundle). Since 235 the **build** catches it 20 minutes later, so this is a time cost, not
   a correctness hole — but it is the last member of the same failure family and it is still open.
   **It was out of this mandate's scope; it should be the next one.**
2. **The panel half of §5.-1's staleness stands.** The H1-union path and the carrying-row median
   have still never run against a low-coverage synthetic D/E/F panel. 240 proves the driver builds
   the right *configuration*; it does not prove the arm behaves at October's coverage.
3. **`submissions/*.meta.json` is untracked** (`git ls-files submissions/` is just `.gitkeep`), so
   the `--incumbent` default points at a file that is not in git. That is the pre-existing
   convention, and it is why the arms are asserted **by value in the tracked script** as well —
   if the meta ever disappears, `def_postbuild.py` fails loudly instead of skipping the check.
4. **This is a wiring fix and nothing else.** No dose, no flag, no seed and no shipped number was
   changed or re-tuned, and no slot is proposed.

---

## 8. Reproduce it

```bash
# rows 2, 3, 4 — 30 s, 0.20 GiB peak, no build
uv run python experiments/240-driver-fix/rows240.py

# row 1 — 22 min build (20.00 GiB) + 26 s prep (19.72 GiB), strictly serial, alone
uv run python scripts/def_pipeline.py --bundle data/challenge \
    --work <tmp> --name <x> --stages build \
    --whitener-out data/donors/whitener_union_pruned.npz --whitener-id union_pruned \
    --cd4-out data/external/cd4t/cloud-compact/cd4_compact
vcc prep <tmp>/<x>.h5ad -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o <tmp>/<x>.vcc
uv run python experiments/240-driver-fix/verify240.py     # digests + config, field by field
```

Artifacts of this experiment: `work/rows240.json` (every case, with the child logs under
`work/probe-*/logs/`), `work/verify240.json`, `work/ka1/logs/build.log`,
`work/ka1/240-ka-h1union.h5ad.meta.json`, `work/ka1/postbuild.json`, `work/ka1/timings.json`,
and the pre-240 byte copies `work/def_{pipeline,postbuild}.PRE240.py` (verified identical to git
HEAD before the change).
