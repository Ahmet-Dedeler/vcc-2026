# 316 — the parity-correct arm `normrel15-pm` (`VCC_NORMREL_PM=blend`): BUILT AND VERIFIED

> **September 9 stage/scope update:** this is the pre-submission build
> handoff. The coordinator later submitted 316 on September 8 at 00:01 UTC:
> **0.21697328748387812**, entry `2K5ppxkowMSj2ugoZqYb` (full capture under
> `data/board-entries/`). The norm identity below is **pre-whitening**;
> experiment 318 §1.2 shows delivered norms do change. The source `.h5ad`
> was removed in the September 8 sweep; package and metadata remain.

**316 (parity-correct `normrel15-pm`): BUILT AND VERIFIED.**

`VCC_NORMREL_PM=blend` implements 313 §8's specification exactly. All three known-answer rows pass
(two of them **bitwise, 0.000000e+00**), the arm is **exactly shape-only on the axis the norm match
controls** (max |delivered-norm ratio − 1| = **4.44e-16** against a 1e-6 tolerance), the single full
build landed with `config.normrel_exp = -1.5` **and** `config.normrel_pm = "blend"`, postbuild is
**33/34 with exactly one FAIL — the expected config-drift line, now naming the new key**, the config
differs from `312-normrel15` in **exactly one key**, prep succeeded, and the novelty check returns
`OK: novel package`.

**Nothing was submitted. `vcc submit` was never run.** The coordinator owns the board reading.

| artifact | value |
|---|---|
| h5ad | `experiments/316-parity-arm/work/build/a-duck-reviews-the-diff.h5ad`, 12,648,807,800 B, sha256 **`5347447bb5ca7d03a534431bb8d9b371978e02e3045549121f13d1470a195d5f`** |
| package | `submissions/316-normrel15-pm.vcc`, 2,824,663,040 B, sha256 `36393703e09b6392…` |
| meta beside it | `submissions/316-normrel15-pm.h5ad.meta.json` (copy of the build's meta) |
| builder, patched | `scripts/make_submission_pm_cd4.py` `d87396e5c9467e38…` (was `dc024795daf9058a…`) |
| driver, patched | `scripts/def_pipeline.py` `b356cae45fff0eef…` (was `33bb0b766c055733…`) |
| revert copies | `work/builder.PRE316.py`, `work/def_pipeline.PRE316.py` (byte copies at launch; **uncommitted**) |
| the patch | `apply_patch316.py` (idempotent), `work/patch316.diff` |

---

## 1. The diff, function by function

Two files, both authorised. `git diff scripts/` is saved verbatim at `work/patch316.diff` (271 lines).
Note one naming correction: 313 §8 calls the second supporting edit's site
`build_production_effect_matrix`; the live function is **`build_production_effects`**
(`scripts/make_submission_pm_cd4.py:2805`). No function of the other name exists.

### (a) module level, `MB:578-594` — the knob
`VCC_NORMREL_PM`, read once at import, default **unset**. Two hard errors, both measured at import:

| condition | measured behaviour |
|---|---|
| `VCC_NORMREL_PM=yes` (any non-`blend`) | `SystemExit: VCC_NORMREL_PM must be 'blend' when set: 'yes'` |
| `VCC_NORMREL_PM=blend` with `VCC_NORMREL_EXP` unset | `SystemExit: VCC_NORMREL_PM is set without VCC_NORMREL_EXP: it changes nothing when the per-target weighting is off, so this is refused rather than silently building the incumbent` |
| `VCC_NORMREL_PM=blend VCC_NORMREL_EXP=-1.5` | imports, `NORMREL_EXP -1.5  NORMREL_PM blend` |

### (b) `build_v9_pre_whiten_effects` (`MB:1407-1440`) — plumbing only
Three new **keyword-only, default-`None`** parameters: `pm`, `operator_columns`, `pm_handled_out`.
They are forwarded to `_normrel_pre_whiten_effects` and are used by nothing else. The
`NORMREL_EXP is None` branch — every non-normrel build, i.e. `cd4svd5` and every arm before 312 — is
textually unchanged.

### (c) `_normrel_pre_whiten_effects` (`MB:1443-1585`) — the arm
Under `VCC_NORMREL_PM=blend`, per target `t` carried by **both** K562-GWPS and PerturbMatch:

```
corr(t)  = alpha * k562_share_or_weight * effect_scale * PM_BLEND * (pm_standardized(t) - incumbent_k562(t))
base_pm  = effect(t, uniform)          + corr(t) @ k562_share      # the blended UNIFORM row = the reference
cand_pm  = effect(t, normrel weights)  + corr(t) @ w_k562(t)       # the blended WEIGHTED row
row(t)   = cand_pm * ||base_pm|| / ||cand_pm||                     # norm match LAST
```

`k562_share` is computed by **`apply_perturbmatch_k562_replacement`'s own expression on its own
objects** (`have` filtered on `active_weights > 0`, `total_weight` summed in `model.names` order),
and `w_k562(t) = w[K562] / sum(w.values())` from the weights the function already computed. The
scalar chain is production's, left to right — `alpha * share * effect_scale * PM_BLEND` **then**
`× drive` — rather than the mandate's algebraically-equivalent `share * (alpha * effect_scale *
PM_BLEND * drive)`. **That spelling is load-bearing and is the one deviation I made**: floating-point
multiplication is not associative, and the production chain is what makes KA row 2 reproduce
`cd4svd5` *bitwise* rather than to the last ulp. It is the same object.

Also under `blend`, **D5 is frozen to 308's support**: the per-target norms `nv` and each donor's
median are taken over `dcols[name] ∩ operator_columns` instead of the donor's whole official
coverage. `production_operator_columns` is panel-free by construction and `load_production_operator`
asserts it, so this is exactly 308's panel-zeroed `wcols` (`E064.operator_columns`: the donor union
minus the panel). **D7 is untouched** — the correction is written on `pm["target_columns"]`,
production's broad experiment-061 policy, unchanged.

Fail-closed: `blend` without `pm`, without `operator_columns`, or on a model with no K562 slot is a
`SystemExit`, not a silent fallback.

### (d) `apply_perturbmatch_k562_replacement` (`MB:1592-…`) — do not blend twice
One new keyword, `already_blended`, default `None`. A target in that set is counted into `support`
and `source_vote_histogram` **exactly as before**, then `continue`s without touching the row.
**`k562_share` is not modified** — the incumbent's PerturbMatch path is bit-identical.
Measured consequence in the real build: the emitted `pm_audit` is **identical to `cd4svd5`'s and to
`312-normrel15`'s, scalars and lists** (248 modified, 21 / 3 / 28 preserved, histogram
`{1: 5, 2: 5, 3: 238}`, 9,373 PM genes).

### (e) `build_production_effects` (`MB:2958-2985`) — threading
Creates `pm_handled: set[str]`, passes `pm=pm, operator_columns=operator_columns,
pm_handled_out=pm_handled` into the pre-whiten call and `already_blended=pm_handled` into the
PerturbMatch call. No other change.

### (f) meta config (`MB:5514`) and (g) `def_pipeline.py:121-125`
`normrel_pm` is written into `config` **only when set** (the `cd4_sizing` / `normrel_exp`
convention), and `VCC_NORMREL_PM` is registered in `ARM_BOARD_DELTA` in the same shape 312 used for
`VCC_NORMREL_EXP`, so the driver announces it rather than treating it as an unexplained key.

---

## 2. The three known-answer rows

Run by `test316.py` on **313's saved fixture** (`work/fixture313.npz`,
`723ce13c0c3b2e63…`, seed 3130907, 20 targets × 60 genes × 3 donors), through the **real production
functions** on a real `MultiDonorModel`, against a **byte copy of the builder as it stood before the
patch** (`work/builder.PRE316.py`, `dc024795daf9058a…`) loaded as an independent module. Output:
`work/test316.json`, `work/test316.log`. Peak RSS **0.22 GiB**.

| row | claim | measured | tol | verdict |
|---|---|---|---|---|
| **KA-1F** | `VCC_NORMREL_PM` **unset**, exp −1.5 == the PRE-316 builder's arm, **bitwise** | **0.000000e+00** | 0 | **PASS** |
| KA-1F′ | the knob unset blends nothing in the pre-whiten path | 0 targets | 0 | PASS |
| **KA-2F** | `VCC_NORMREL_PM=blend` with **exp 0** == the PRE-316 builder's **incumbent**, **bitwise** | **0.000000e+00** | 0 | **PASS** |
| KA-2F′ | the in-place blend leaves `pm_audit` unchanged | identical dict | — | PASS |
| **KA-3** | delivered row norm / incumbent == 1 for **every** target (full axis, relative) | **4.440892e-16** | 1e-6 | **PASS** |
| KA-4F | 313's own `parity_variant` reproduced from its clean-room code | 0.000000e+00 | 1e-9 | PASS |

**6/6.** At exponent 0 the blend path fired on 17/20 fixture targets and still returned the
incumbent bit for bit, which is the check that matters: the same code that re-shapes the row at
−1.5 is provably the identity at 0.

### What rows 1 and 2 are, and are not
The mandate states rows 1 and 2 as h5ad `output.sha256` comparisons and, in the same breath, budgets
**one** full build and says to check them on the fixture first. Each production-scale byte-for-byte
row would need its own ~28-minute build (unset at −1.5, and blend at 0), so **rows 1 and 2 are
discharged on the fixture, at 0.0 bitwise, through the production functions**, not on the real
panel. What *was* verified at production scale instead:

* the built config differs from `submissions/312-normrel15.h5ad.meta.json` in **exactly one key**
  (`normrel_pm: absent → "blend"`) — §4;
* the emitted `pm_audit` is **identical** to `312-normrel15`'s and `249-cd4svd5-f034`'s;
* the unset code path is textually the shipped expression (`work/patch316.diff`: with
  `NORMREL_PM is None`, `norm_cols[name] is model.dcols[name]` and the blend block is skipped).

### KA-3, and the one place parity with 308 is structural rather than exact
The norm match is on the **full row** — production's `MB:1464-1467` and the mandate's formula. 308's
harness rows are panel-**zeroed** (`062/run.py:287`), so for 308 the full row *is* the non-panel row
and both read 1.000000. Production's rows keep the panel columns, so the two are different numbers
here. Measured on the fixture:

| | full axis (what the norm match controls) | non-panel (scored) columns |
|---|---|---|
| **316 parity arm** | median 1.0000000, range **[1.0000000, 1.0000000]** | median 0.997580, range [0.634414, 1.299320] |
| 312 shipped arm | median 1.000000, range [0.993001, 1.340741] | median 1.004769, range [0.806168, 2.154923] |
| whole-matrix ratio | parity **1.0000000** · shipped 1.060008 | — |

313's fixture panel is **33 %** of the axis (61 % of row energy) against production's
**300/18,533 = 1.62 %**, so the non-panel column is the fixture's artefact, not production's number.
`scale316.py` sweeps 313's own construction over axis width, running the real builder functions at
each width (`work/scale316.json`):

| genes | panel % cols | panel % energy | parity non-panel [med, min, max] | shipped non-panel [med, min, max] | parity ‖M‖ | shipped ‖M‖ |
|---:|---:|---:|---|---|---:|---:|
| 60 | 33.33 | 59.47 | 0.997580 0.634414 1.299320 | 1.004769 0.806168 2.154923 | 1.0195 | 1.1603 |
| 120 | 16.67 | 31.27 | 1.001041 0.970171 1.539397 | 1.000000 0.950346 1.706455 | 1.1284 | 1.1738 |
| 300 | 6.67 | 16.76 | 1.000263 0.988161 1.256588 | 1.000063 0.984770 1.700509 | 1.0589 | 1.1802 |
| 600 | 3.33 | 9.58 | 1.000000 0.994513 1.080015 | 1.000925 0.991604 1.482582 | 1.0180 | 1.1208 |
| **1235 (production's 1.62 %)** | **1.62** | **4.79** | **1.000195 0.991720 1.043764** | 1.001391 0.994958 1.464702 | **1.0099** | 1.1220 |

The full-axis delivered norm is **exactly 1.0000000 at every width** for the parity arm, and
[0.996398, 1.406781] for the shipped one at 1235 genes. **At production's panel share the parity arm
removes about 92 % of the arm's amplitude change on the scored columns** (matrix +0.99 % against
+12.20 %) and all of it on the full row. Per 313 §9's rule, the *structure* transfers and the
*magnitudes* are fixture-specific: I did not compute an `nmae` identity leg, which needs the real
rows.

---

## 3. The fixture comparison against 313's `parity313.json`

313's `parity_variant` froze the two second variables the **other** way from this mandate: it took
the weights over the **full axis** (D5 = production's choice) and on the **PM-blended** rows (D3 =
308's choice). Recomputed here from 313's own clean-room functions:

| quantity | 313's `parity313.json` | recomputed by `test316.py` |
|---|---|---|
| max \|diff\| vs the shipped arm | 1.7297e+00 | **1.7297e+00** |
| median row cosine vs shipped | 0.998189 | **0.998189** |
| min row cosine vs shipped | 0.873612 | **0.873612** |

**Reproduced exactly** (KA-4F, 0.0e+00). The arm this mandate specifies — D5 frozen to 308's `wcols`,
weights on the **raw** rows as `MB:1444-1457` already computes them — differs from that variant by
max \|diff\| **3.9951e-01**, median row cosine **0.995875**. Against the shipped arm the built arm
reads max \|diff\| **2.0866e+00**, median row cosine **0.993225**, min **0.873612**; against the
incumbent `cd4svd5`, max \|diff\| **1.7553e+00**. It is a real re-shaping, not a cleanup, and must be
priced as one.

---

## 4. The build

One invocation of the gated two-invocation pattern (`run316.sh`, modelled on `run312.sh`), on the
real A/B/C bundle. `work/build/run316.log`, `work/build/driver.log`, `work/build/logs/`.

```
2026-09-07T20:01:53Z ARMED  VCC_NORMREL_EXP=-1.5  VCC_NORMREL_PM=blend
  builder d87396e5c9467e38…  driver b356cae45fff0eef…  git 1544b5f
  gate: free 83%  free+inactive 32.6 GiB  other builders 0  vcc preps 0
2026-09-07T20:29:11Z invocation 1 rc=1   (the expected postbuild FAIL; artifact written)
  restored scripts/make_submission_pm_cd4.py from pin-artifacts.builder-before.py
    (87b39d57bd14f697… -> d87396e5c9467e38…, identical)
```

**The postbuild line, verbatim:**

```
33/34 checks pass       FAIL lines: 1
  FAIL  config differences vs 249-cd4svd5-f034.h5ad.meta.json are only
        ['contexts', 'expect_contexts', 'targets', 'whitener_id']:
        moved=['expect_contexts', 'normrel_exp', 'normrel_pm', 'whitener_id']
```

This is the same single FAIL 312 got (`33/34`, `moved=['expect_contexts', 'normrel_exp',
'whitener_id']`), with `normrel_pm` added: the guard noticing a deliberate arm. `run316.sh` refuses
to prep if the FAIL count is anything but 1 or if the one FAIL is not the config-drift line.

**The arm landed.** `config.normrel_exp = -1.5`, `config.normrel_pm = "blend"`.

**Config diff against the object this arm is one variable away from:**

| vs | keys differing | |
|---|---|---|
| **`312-normrel15`** | **1** | `normrel_pm: <absent> → 'blend'` |
| `249-cd4svd5-f034` | 4 | `normrel_exp`, `normrel_pm`, plus `expect_contexts` and `whitener_id union_pruned → union_pruned_def` (the same two driver keys 312 moved) |

**Where the change entered, from the meta's own stage digests:**

| digest | 316-normrel15-pm | 312-normrel15 | 249-cd4svd5 |
|---|---|---|---|
| `baseline_pre_sha256` | `92555af13c543f2d` | `9e63576df8436257` | `7678bab9ea9d1b8d` |
| `candidate_pre_sha256` | **`92555af13c543f2d`** | `9b19eff922fc2278` | `0601336fdf1ec5c5` |
| `post_whiten_sha256` | `218e399b4d83df6e` | `f3df0f6426ff1bdc` | `b2426e15d76d9fe5` |
| `final_sha256` | `6c25cdc875b92878` | `ab8997931a48a191` | `a28d5be2e9cb3ed1` |

For 316, **`baseline_pre == candidate_pre`** — `apply_perturbmatch_k562_replacement` changed nothing,
because all 248 both-supported targets were already blended inside the arm. For 312 and for
`cd4svd5` the two digests differ. That is the machine-checkable signature that the correction is
applied **once**, in the right place, and that the double-application was prevented.

`nnz` 1,577,896,070 (312: 1,577,401,690; `cd4svd5`: 1,577,615,064); bytes 12,648,807,800, identical
to both.

**Resources.** whitener 328.3 s / 12.32 GiB · repin-panel-support 10.3 s / 20.14 GiB · **build
1276.6 s / 20.23 GiB** · postbuild 0.7 s. Every heavy step ran serial behind the ≥ 40 % free,
≥ 23 GiB free+inactive, no-other-builder, no-other-prep gate. My own analysis jobs peaked at
**0.22 GiB**; the build and prep are the mandate's explicitly authorised heavy steps.

---

## 5. Prep and novelty

```
✓ prep complete
  input:         .../experiments/316-parity-arm/work/build/a-duck-reviews-the-diff.h5ad
  output:        .../submissions/316-normrel15-pm.vcc
  cells:         360000
  genes:         18533  (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts:      column 'context' → A=120000, B=120000, C=120000
  targets:       verified against the official list
  normalization: counts-preserved
       23.68 real   peak memory footprint 16,046,550,800 B (15.38 GiB max RSS)
```

The driver's `prep` stage runs `vcc prep --dry-run` (validation only, writes nothing), so the
package was written by a real `vcc prep` in the same gated slot.

```
$ uv run python scripts/check_submission_novelty.py check submissions/316-normrel15-pm.vcc
hashing submissions/316-normrel15-pm.vcc (2.82 GB)...
prediction identity: 5347447bb5ca7d03... (output.sha256 of 316-normrel15-pm.h5ad.meta.json (container shape agrees))
content identity: bb293c678c8ebf45... (content of pred.h5ad.zst (X, obs index, target_gene, context, var index))
OK: novel package, sha256 36393703e09b6392...
                                                                              rc=0
```

---

## 6. State left behind, and what the coordinator must decide

* **`scripts/make_submission_pm_cd4.py` and `scripts/def_pipeline.py` are PATCHED AND UNCOMMITTED**
  (`d87396e5c9467e38…`, `b356cae45fff0eef…`). Byte copies of the pre-patch files are at
  `work/builder.PRE316.py` / `work/def_pipeline.PRE316.py` (both `dc024795daf9058a…` /
  `33bb0b766c055733…`, i.e. HEAD `1544b5f`). 312 committed its patch before building; I did not
  commit, per the standing rule about git with agents running. `apply_patch316.py` is idempotent and
  re-derives the patch from the pre-patch text.
* **Nothing else under `submissions/`** was written: only `316-normrel15-pm.vcc` and its sibling
  meta.
* **`vcc submit` was never run**, and no slot was spent.
* **The arm is not priced here.** 313 §8's own note stands: whether it is worth a slot depends on
  re-pricing ORD_A's `nmae` leg on the shipped-scale rows, which this mandate did not ask for and
  which needs the real rows. What this experiment establishes is that the object now buildable **is**
  ORD_A's ordering on production's substrate, exactly norm-preserving on the row the norm match
  controls, at a measured +0.99 % residual amplitude on the scored columns at production's panel
  share against the shipped arm's +12.20 %.

## 7. Limitations

* Rows 1 and 2 are **fixture-scale bitwise**, not production-scale byte-for-byte: two more full
  builds would be needed and the mandate budgets one. §2 lists what replaced them at production
  scale.
* The fixture is synthetic (313 §9). Structural statements — which orderings coincide, that the arm
  is norm-preserving, which targets are blended, that exponent 0 is the identity — are exact and
  transfer. Magnitudes (2.0866, the non-panel ranges, +0.99 % vs +12.20 %) are fixture-specific.
* D5 and D3 are **frozen, not swept**, as instructed: D5 to 308's `wcols`, D3 to production's raw
  rows (the spec's `w_k562(t)` is the weight `MB:1444-1457` already computes), D7 to production's
  broad column policy. The D5 + D3 freeze is worth max \|diff\| 3.9951e-01 on the fixture against
  313's own variant — a second variable that a later slot could test, not this one.
* `BETA = 0.0` throughout, so the `beta * m_full` offset inside the norm match never fires; 313's
  latent fourth divergence at `BETA ≠ 0` is untouched and still latent.
* No `nmae` identity leg, no PLOLO retrodiction, no board reading. 313 §8's KA row (4) — re-running
  308's `normrel1.5` fold to retrodict +0.00259 at 0.0e+00 — was not in this mandate's three rows and
  was not run.

## 8. Files

| file | note |
|---|---|
| `apply_patch316.py` | the patch, idempotent; saves the pre-patch byte copies |
| `test316.py` | the three known-answer rows on 313's fixture (+ 313's variant reproduced) |
| `scale316.py` | the panel-share sweep of the delivered non-panel norm |
| `run316.sh` | the gated build + prep; refuses to prep unless the one FAIL is the config-drift line |
| `work/patch316.diff` | `git diff scripts/`, 271 lines |
| `work/test316.json`, `work/test316.log` | KA rows, delivered norms, the 313 comparison |
| `work/scale316.json`, `work/scale316.log` | the width sweep |
| `work/build/run316.log`, `work/build/driver.log`, `work/build/logs/postbuild.log` | the build |
| `work/build/prep.log` | the prep |
| `work/build/a-duck-reviews-the-diff.h5ad{,.meta.json}` | the artifact, `5347447bb5ca7d03…` |
| `work/builder.PRE316.py`, `work/def_pipeline.PRE316.py`, `work/builder.LAUNCH.py`, `work/def_pipeline.LAUNCH.py` | revert copies |

Sources read: `experiments/313-normrel-parity/{RESULTS.md,parity313.py,fixture313.py,scale313.py,work/parity313.json}`,
`experiments/312-normrel-arm/{apply_patch.py,run312.sh,run312b.sh,work/build/*}`,
`scripts/{make_submission_pm_cd4.py,def_pipeline.py,check_submission_novelty.py}`,
`src/vcc2026/multi_donor.py`, `experiments/308-kd-weighted-donors/{lib308.py,sweep308.py}`,
`experiments/263-donor-selection-arbiters/lib263.py`, `experiments/064-pm-cd4-composition/run.py`,
`experiments/162-retrodicting-validator/plolo.py`, `docs/READ-FIRST.md`, `AGENTS.md` §8.1d/§11.2b,
`submissions/{312-normrel15,249-cd4svd5-f034}.h5ad.meta.json`.
