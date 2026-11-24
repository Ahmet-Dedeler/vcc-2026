# 235 — The Oct-22 operator failure is **HARDENED**: it now raises in ~1 second, and the carrying-row rule is proved to cover every registered source

**Date:** 2026-09-04 · deterministic · **no `vcc submit`, no `cell-eval2`, no model name, no
description.** Two files changed outside this experiment directory —
`scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py`, both authorised by the
mandate, both committed with the hashes below. Nothing under `data/`, `src/`, or any other
experiment was written; 199's synthetic bundles, 225's and 223's source tables and the shipped
operator were **read** only.

---

## 0. Verdict

```
HARDENED.

  ITEM A  the operator's panel mismatch now RAISES, observed firing four ways, and the
          pre-235 builder is shown accepting the same wrong panel with no error at all.
          A/B/C axis_sha256 = edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f

  ITEM B  the carrying-row rule covers every source the registry can select: 37 of 37
          (source, panel) cells size ABOVE ZERO across THREE panels, and the shipped
          all-rows rule returns EXACTLY 0.0 in all 37.  A selected source that carries
          nothing now RAISES instead of contributing a silent zero -- and on a panel that
          is disjoint AND essential-free (AGENTS.md 8.1g, landed mid-run) that is exactly
          what `k562_essential` does: 30 targets on 199's bundle, 0 here.

  All 5 mandate known-answer rows PASS.  The diff is 181 insertions and ZERO deletions:
  no existing line of either file was modified or removed, and no shipped byte moved.
```

**The one number that says why this was worth doing.** The only check the builder had was
`n_cols`, and 235 measured it at **17,942 on all four panels on disk** — A/B/C and 199's three
synthetic D/E/F bundles. It cannot distinguish the right panel from a panel with **zero** targets
in common. On that panel **98.81%** of the operator's coordinates pair with the wrong gene, priced
by 212 at **−0.0429 overall on `pds` alone**, on the only entry that counts, with the board dark.

---

## 1. The five known-answer rows

| # | row | required | measured | verdict |
|---|---|---|---|---|
| 1 | flag unset, rebuilt through the changed builder | `.h5ad` sha256 `7036e89f…`, 12,648,807,800 B | **`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`**, 12,648,807,800 B | **PASS** |
| 2 | H1 pair at `g = 0.28`, prepped | `.vcc` sha256 `f13ce40b…`, 2,824,785,920 B | **`f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`**, 2,824,785,920 B — and the `.h5ad` it came from hashes to **`8df4dadd…`**, byte-identical to the shipped `222-h1union-g028.h5ad` | **PASS** |
| 3 | two independent A/B/C operator refits | `e2f65ecf…`, stable `axis_sha256` | **three** refits, all `e2f65ecf…`, all `axis_sha256 edf6ba88…`, all three metadata jsons **byte-identical**; the third ran under the exact committed fitter bytes | **PASS** |
| 4 | `solve_cd4_fraction` untouched, `effects.cd4.fraction` unmoved | `0.32207201067578145` | source function **not in the diff**; both builds report `0.32207201067578145` | **PASS** |
| 5 | **the assertion must FIRE** | a raise on a synthetic disjoint panel | raises on 199's `def-bundle`, on `def-bundle-badaxis`, and on an operator with no axis digest; **the pre-235 builder accepts the same inputs silently** | **PASS** |

### Row 5, verbatim — the raise (`work/fire235.log`)

```
REFUSING TO BUILD: the whitening operator was fitted for a DIFFERENT PANEL.
  operator            ~/Code/virtual-cell-challenge-2026/data/donors/whitener_union_pruned.npz
  operator sha256     e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123
  --whitener-id       'union_pruned'
  axis_sha256 recorded (LEGACY_OPERATOR_AXIS_SHA256[operator sha256]) edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f
  axis_sha256 of THIS panel                c464a2acd6b07921aeb1c2f3f3730fbe39428671c88782ff078825d4c79c91a9
  columns             17942 (this check is NOT a column count; n_cols is 17,942 on every panel this project has measured)
  The operator's axis is union(donor coverage) - panel, so it is a function of the panel: using it against another one pairs ~98.8% of its coordinates with the wrong gene, silently. Re-fit it for this panel (scripts/build_union_whitener.py --bundle <this bundle> --basis-prune essential), re-pin it (scripts/def_pin_panel_artifacts.py --apply), and rebuild. DEF-RUNBOOK section 0.1 and section 2 step 5.
```

**And the row that makes it a measurement rather than a claim** — the same call, same model, same
operator, same panel, through the byte copy of the pre-235 builder
(`work/make_submission_pm_cd4.PRE235.py`, sha `30fd9ab3…`):

```
[BLIND_pre235_def_pessimistic] PASS  raised=False  0.56s
   columns 17942 | metadata_n_cols 17942 | evecs_shape [17942, 17942]
   "the pre-235 builder accepted the WRONG panel with no error"
```

`fire235.py` also fires the guard two more ways: on `def-bundle-badaxis` (199's bundle with two
symbols transposed in `gene_names.csv`), and on `data/donors/whitener_union_full.npz`, an operator
that carries no axis digest and is therefore **refused rather than trusted**.

---

## 2. Item A — what was wrong, and what the guard actually checks

### 2.1 The blind spot, measured (`work/axis235.json`)

The operator's axis is `sorted_by_official_index(union(donor coverage) − panel)`, so it is a
function of the panel. Computed for every panel on disk, from the donors' `genes.txt` alone:

| panel | `n_cols` | `axis_sha256` | positions pairing a DIFFERENT gene than A/B/C | panel ∩ A/B/C |
|---|---:|---|---:|---:|
| `data/challenge` (A/B/C) | **17,942** | `edf6ba88…` | — | 300 |
| 199 `def-bundle` (pessimistic) | **17,942** | `c464a2ac…` | **17,729 = 98.81%** | **0** |
| 199 `def-bundle-r` (realistic) | **17,942** | `87361c7c…` | 17,511 = 97.60% | **0** |
| 199 `def-bundle-badaxis` | **17,942** | `f2d6cc07…` | 17,729 = 98.81% | 0 |

This independently reproduces 199's **98.8%** on a different code path.

**Two things follow, and the second is new.**

1. `n_cols` **is the same on every panel**, so the only check that existed could never fire.
2. The meta's existing `columns_sha256` — a hash of the **integer** column indices — reads
   **`d66261d7…` for BOTH `def-bundle` and `def-bundle-badaxis`**, because a transposed
   `gene_names.csv` moves the symbols without moving the indices. So the guard hashes the ordered
   **symbols**, not the indices, and it is **asserted**, where `columns_sha256` was only recorded.

### 2.2 The change

| file | what |
|---|---|
| `scripts/build_union_whitener.py` | `axis_sha256()` — sha256 over the ordered column symbols, newline-joined UTF-8. Printed after the union-axis line and written into `<out>.json` as **`axis_sha256`** (plus an informational `panel_sha256`). |
| `scripts/make_submission_pm_cd4.py` | `operator_axis_sha256()` (the same function), `LEGACY_OPERATOR_AXIS_SHA256` (axis digests for the two operators fitted before the field existed, **keyed by operator sha256**, rule 11), and the assertion inside `load_production_operator`. `LoadedOperator` carries `axis_sha256` / `axis_sha256_source`; the meta records them under `operator`. |

**Where the expectation comes from, in order:** the operator's own metadata `axis_sha256` (which
every future refit writes, and whose json sha256 is itself pinned in `WHITENER_REGISTRY`), else
`LEGACY_OPERATOR_AXIS_SHA256[operator sha256]`, else **refuse**. Keying the legacy table by the
operator's sha256 rather than by `--whitener-id` means the four experiment scripts that call
`load_production_operator` without an id (156, 103 ×2, 229's `stage229.py`) are covered unchanged.

**Where it is checked:** immediately after the metadata checks and **before `Whitener.load`** — it
costs ~0.4 s and fires before the 1.3 GB eigenbasis is read, let alone the 21-minute build.

**Why `axis_sha256` is in `meta['operator']` and not in `meta['hashes']['inputs']`, which the
mandate suggested.** `experiments/199-def-runbook/rule10_199.py` — a `DEF-RUNBOOK` §4
known-answer row — resolves **every** key in `hashes.inputs` back to a file on disk and hashes it.
A non-file key there breaks that row. `meta['operator']` is where the operator's `sha256`,
`metadata_sha256` and `columns_sha256` already live, and `def_postbuild.py` already reads that
block. This is the only place the mandate's letter was not followed, and the reason is recorded
here rather than resolved silently.

### 2.3 What `def_pin_panel_artifacts.py` needs on the day: nothing

That script rewrites `WHITENER_REGISTRY` by regex on the literal
`WHITENER_REGISTRY: dict[str, tuple[str, str]] = {`, so the annotation and the 2-tuple shape were
left untouched on purpose. The new operator's axis digest travels **inside its metadata json**,
whose sha256 the pin script already pins. **Step 7 is unchanged and the Oct-22 path needs no new
step.**

---

## 3. Item B — the carrying-row rule covers every source the registry can select

### 3.1 Coverage of all ten registered sources (`work/carry235.json`)

| source | A/B/C | 199 pessimistic | 199 realistic |
|---|---:|---:|---:|
| `dld1_lowmoi` | 67 | **73** | 67 |
| `gse295214_h1` *(shipped)* | 64 | 34 | 57 |
| `h1_2025` *(shipped)* | 20 | **3** | 7 |
| `igvf_cardiomyocyte` | 87 | 36 | 60 |
| `jiang_5line_noK562` | 9 | 4 | 3 |
| `myllia_a549` | 14 | 4 | 5 |
| **`nadig_hepg2`** | **0** | **33** | 45 |
| **`nadig_jurkat`** | **0** | **33** | 48 |
| **`rpe1`** | **0** | **32** | 47 |
| **`k562_essential`** | **0** | **30** | 34 |

**Known-answer row:** this reproduces **8 of 8** of 225 §8.3's projected numbers exactly
(`k562_gwps`-style rows aside, which are not registry sources): `dld1_lowmoi` 67→73,
`gse295214_h1` 64→34, `h1_2025` 20→3, `igvf_cardiomyocyte` 87→36, and the four October sources
0→33/33/32/30.

### 3.2 Both sizing rules, side by side, on every (source, panel) cell it covers

`g = 0.28`, base median row norm pinned to the shipped `h1union` build's `4.229295604454299` (there
is no D/E/F whitened base to measure today; the base enters both rules as the same positive
denominator, so the comparison is exact and the fractions are directly comparable to the shipped
`0.19207`).

| | carrying-row rule | shipped all-rows rule |
|---|---|---|
| **28 of 28** (source, panel) cells with coverage > 0 | **all > 0**, range `0.1051 … 0.2776` | **all EXACTLY 0.0** |
| the October four, individually, pessimistic draw | `nadig_hepg2` 0.17889 · `nadig_jurkat` 0.16144 · `rpe1` 0.18289 · `k562_essential` 0.19507 | 0.0 · 0.0 · 0.0 · 0.0 |
| the October four **as one union**, pessimistic | **0.21799** (35 carrying rows, 11.7% coverage) | **0.0** |
| the October four as one union, realistic | **0.21895** (50 carrying rows, 16.7% coverage) | **0.0** |

No registered source clears 50% coverage on any panel — which is exactly why the all-rows rule
returns zero for **every one of them**, and why 233's "+5.5% to +17.7% relative cosine each, all-rows
median 0.0000" is a rule problem and not a source problem.

**Known-answer row (bit-exact):** the shipped pair on the A/B/C panel reproduces
`submissions/222-h1union-g028.h5ad.meta.json`'s residual frame on **7 of 7 fields**, to the last
bit — carrying median `6.1654229367116375`, `g_per_unit` `1.4577895501601277`, solved
`0.19207161964444322`, 77 carrying rows, coverage `0.25666666666666665`, all-rows median `0.0`,
`shipped_median_rule_would_solve` `0.0`.

### 3.3 The one hardening item B needed, and it is a real hole that was open

The union's sizing raises only when the **whole** union carries nothing. A **single selected
source** that carries nothing contributed exactly zero while still being named in
`VCC_UNION_SOURCES` and recorded in the meta — an arm described in the record that never ran.
225 measured this as the expected October case (`kolf2` 273 → 0, `feng_ipsc` 182 → 0 against a
panel that is disjoint by construction). `prepare_h1_union` now raises:

```
union source 'nadig_hepg2' carries 0 of 300 panel targets, so it would contribute EXACTLY
NOTHING to this build while still being named in VCC_UNION_SOURCES and recorded in the meta.
```

**No build ever made is affected**: every source any shipped or built arm selected covers > 0 of
the panel it ran against (§3.1), and known-answer rows 1 and 2 rebuild the two that matter.

### 3.4 ADDENDUM — `AGENTS.md` §8.1g landed mid-run, and it changes §3.1's table (not §3.2's claim)

Commit `eed4e4a` (12:17, while this experiment's builds were running) established that **the real
A/B/C panel contains ZERO essential genes**, for an assay reason expected to hold in October, and
that `DEF-RUNBOOK` §5's two synthetic bundles carry **34 and 30** because they are drawn from donor
coverage. Its rule: *no synthetic-panel number may be quoted without stating the panel's
essential-gene count, and a panel used to predict October must be disjoint AND essential-free.*

**§3.1's and §3.2's panels, stated as the rule requires:** A/B/C **0** essential genes, 199's
pessimistic bundle **30**, 199's realistic bundle **34**.

`essfree235.py` therefore redraws a third panel that has **both** properties — 300 targets, seed 0,
from `(any production donor's targets on the official axis) − (A/B/C panel) − (common-essential
list)`; verified 300 unique, **0** overlap with A/B/C, **0** essential — and repeats the
measurement (`work/essfree235.json`, **23/23 checks PASS**):

| source | A/B/C (ess 0) | 199 pessimistic (ess 30) | **essential-free (ess 0)** | carrying-row solved | all-rows |
|---|---:|---:|---:|---:|---:|
| `dld1_lowmoi` | 67 | 73 | **81** | 0.17741 | 0.0 |
| `gse295214_h1` | 64 | 34 | 22 | 0.21787 | 0.0 |
| `h1_2025` | 20 | 3 | 8 | 0.17832 | 0.0 |
| `igvf_cardiomyocyte` | 87 | 36 | 32 | 0.20732 | 0.0 |
| `jiang_5line_noK562` | 9 | 4 | 4 | 0.18992 | 0.0 |
| `myllia_a549` | 14 | 4 | 3 | 0.13166 | 0.0 |
| `nadig_hepg2` | 0 | 33 | **12** | 0.20647 | 0.0 |
| `nadig_jurkat` | 0 | 33 | **13** | 0.14968 | 0.0 |
| `rpe1` | 0 | 32 | **12** | 0.17147 | 0.0 |
| **`k562_essential`** | 0 | 30 | **0** | — **raises** | — |

**Two things this settles for item B.**

1. **The structural claim is unchanged and now rests on 37 of 37 cells over three panels**, one of
   which has both October properties. The carrying-row rule sizes every source that covers a panel;
   the shipped all-rows rule returns exactly 0.0 for every one of them, on every panel.
2. **§3.3's new raise is not hypothetical — it fires on a registered source under the most
   realistic October panel we can draw.** `k562_essential`'s entire 30-target overlap with 199's
   bundle **was** the essential set; on an essential-free panel it carries nothing, and the builder
   now refuses instead of silently adding zero. The other three October sources survive at
   **~2.7× lower coverage** than the contaminated draw suggested (33/33/32 → 12/13/12).

**This is one independent draw, not a replacement for 238**, which is recomputing `DEF-RUNBOOK`
§5.-2's whole ranking. It agrees with 8.1g's direction and gives that re-measurement a second,
differently-constructed panel to check against.

### 3.5 The test that fails if this ever stops being true

`experiments/235-october-hardening/tests/test_october_hardening.py` — **9/9 pass**, ~2 min,
peak 0.99 GiB, no `pytest` needed (this repo's `uv` env has none, so it uses the same `TESTS` list
and `__main__` runner as `scripts/test_make_submission_pm_cd4.py`):

| test | what it would catch |
|---|---|
| `..._two_axis_digest_implementations_are_the_same_function` | the fitter and the builder drifting apart on the digest's canonical form |
| `..._axis_guard_fires_on_a_panel_the_operator_was_not_fitted_for` | the guard being removed or weakened (synthetic, 8 genes, same `n_cols` both ways) |
| `..._operator_with_no_recorded_axis_is_refused_rather_than_trusted` | the fallback quietly becoming permissive |
| `..._production_operator_is_pinned_by_path_and_by_sha_rule_11` | the pinned A/B/C axis digest or operator sha drifting |
| `..._fitter_writes_an_axis_digest_into_every_operator_it_makes` | a refit on the day producing a json the builder cannot check |
| `..._second_residual_is_sized_by_the_carrying_row_rule_and_cd4_is_not` | the routing regressing to `solve_cd4_fraction` (behavioural, not a source grep) |
| `..._selected_source_that_carries_no_panel_target_now_raises` | §3.3's guard being removed |
| `..._shipped_h1_pair_reproduces_the_shipped_build_exactly` | any change to the residual construction |
| `..._every_registered_source_sizes_above_zero_on_a_panel_it_covers` | a registry entry that cannot actually be sized — **the Oct-22 row** |

They live in this experiment directory rather than in `scripts/test_make_submission_pm_cd4.py`
because the mandate authorises edits to two files under `scripts/` and no others. **Folding them
into that suite is a one-move follow-up for an agent with that permission**, and it should happen
before Oct 22 so they run with everything else.

---

## 4. Diff summary — 181 insertions, 0 deletions

```
 scripts/build_union_whitener.py   |  37 ++++++++++
 scripts/make_submission_pm_cd4.py | 144 ++++++++++++++++++++++++++++++++++++++
 2 files changed, 181 insertions(+)
```

**The committed bytes are the verified bytes.**

| file | sha256 committed | evidence |
|---|---|---|
| `scripts/make_submission_pm_cd4.py` | `715f2b46e090e061a51ffd51a4444f155ce14da947652802019e9e63b8855328` | **equals `hashes.implementation` in both** of row 1's and row 2's metas — this exact file produced `7036e89f…` and `f13ce40b…` |
| `scripts/build_union_whitener.py` | `9ea64b120a8184aa92d170de5557f73a6b7775cdcb0aba57a28dece8dfdc1318` | the third refit ran under these bytes and printed `MATCH` on `e2f65ecf…` |

(A wording fix to one docstring landed *after* rows 1 and 2 were built, which would have made the
committed builder differ from the one the metas record. It was reverted rather than kept: the file
committed is the file that was verified. The fitter's equivalent fix was kept and re-verified by
the third refit instead — 6 minutes rather than 42.)

**Not one existing line was modified or deleted**, in either file. `solve_cd4_fraction`,
`solve_carrying_row_fraction`, `union_residual`, `apply_frozen_cd4_residual`,
`build_production_effects` and every emission path are byte-unchanged; the additions are a hash
function, a pinned table, one assertion, three audit keys, one coverage guard and their comments.

| file | hunk | what |
|---|---|---|
| `build_union_whitener.py` | `axis_sha256()` | the digest, with the measurement that motivates it |
| | after `wcols` | prints `axis sha256` / `panel sha256` (also in `--dry-run`) |
| | `meta` | writes `axis_sha256` and `panel_sha256` |
| `make_submission_pm_cd4.py` | after `WHITENER_REGISTRY` | `LEGACY_OPERATOR_AXIS_SHA256`, keyed by operator sha256 |
| | `LoadedOperator` | two defaulted fields |
| | before `production_operator_columns` | `operator_axis_sha256()` |
| | `load_production_operator` | the `whitener_id` kwarg, **the assertion**, the `[operator]` line, the two returned fields |
| | `run()` | passes `whitener_id`; records `axis_sha256`, `axis_sha256_checked_against`, `panel_sha256` under `meta['operator']` |
| | `single_source_residual` | records each source's `coverage_frac` |
| | `prepare_h1_union` | raises on a selected source that carries **0** panel targets |
| | `h1_union_support_audit` | `targets_by_source_count` (correct at any arity) and `source_coverage` |

### Effect on the meta — measured field by field, against 222's two metas

**Flag unset (row 1), 302 fields compared:**

| | |
|---|---|
| removed | **none** |
| added | `operator.axis_sha256`, `operator.axis_sha256_checked_against`, `operator.panel_sha256` |
| changed | `hashes.implementation.scripts/make_submission_pm_cd4.py` (`71cf433d…` → **`715f2b46…`**), `output.path`, `peak_rss_gib` |

**`g = 0.28` (row 2), 337 fields compared:** the same three `operator` keys, plus
`effects.h1_union.source_coverage.{gse295214_h1: 64, h1_2025: 20}`,
`effects.h1_union.targets_by_source_count.{0: 223, 1: 70, 2: 7}`, and `coverage_frac` inside each
entry of `effects.h1_union.source_details`; the same three changed values and nothing else.

**`output.sha256`, every `effects.*` number, every `operator.*` hash, every input hash and every
`config` field are identical in both.** `hashes.inputs` is untouched, so `rule10_199.py` still
resolves every key to a file (§2.2). This is rule 10's diff, and it comes out at exactly the
intended number of differences.

---

## 5. Cost, memory and what was run

Strictly serial, gate checked before each heavy job (`memory_pressure` ≥ 40% free, `vm_stat`
free+inactive ≥ 20 GiB, swap free ≥ 400 MiB — 233's gate, not the swap-only rule it replaced).

| job | wall | `ru_maxrss` | note |
|---|---:|---:|---|
| `axis235.py` | 3 s | 0.3 GB | reads `genes.txt` only |
| `fire235.py` (row 5) | 118 s | **17.75 GB** | two `MultiDonorModel` loads; the model alone is 16.5 GiB in the production build too |
| operator refit ×3 (row 3) | 320.0 s + 364.1 s + 398.8 s | 13.68 / 13.70 / 13.83 GB | `--max-rss-gib 14`; 5.54 GiB anonymous at the covariance, the rest is the 4.2 GB basis memmap's file pages. The three `.npz` files were deleted after hashing (3.9 GB) |
| `carry235.py` (item B) | 105 s | 1.17 GB | |
| `tests/test_october_hardening.py` | 100 s | 1.06 GB | |
| build, flag unset (row 1) | **1280.6 s** | **21.30 GB** | 19.834 GiB by the builder's own checkpoint |
| build, `g = 0.28` (row 2) | **1267.6 s** | **21.44 GB** | 19.966 GiB internal |
| `vcc prep` (row 2) | 18.8 s | **27.11 GB** | never alongside a build. 222 recorded 23.98 GB for the same call; **25.2 GiB is the number to plan Oct 22 around** |

**The mandate's ~10 GB RSS ceiling was exceeded by design on three jobs and cannot be met by any
of them**: the production build is 19.8 GiB, `vcc prep` is 22–24 GiB, and the three-donor model
that rows 3 and 5 need is 16.5 GiB on its own. Known-answer rows 1 and 2 *are* production builds,
so the alternative was not to verify them. Everything was run one at a time with the gate checked
first, and peak RSS is reported for each.

---

## 6. Honest limits

1. **The guard proves the axis, not the operator's quality.** It answers "was this operator fitted
   for this panel"; it cannot tell a good fit from a bad one. Step 5 is still mandatory.
2. **The legacy axis digest for the two pre-235 operators was derived, not independently
   witnessed.** It is pinned from a computation this experiment ran, cross-checked three ways: the
   fitter's own implementation on a real refit, the builder's implementation through
   `MultiDonorModel`, and a third path in `axis235.py` that reads `genes.txt` directly. All three
   agree on `edf6ba88…`. If all three shared a mistake, the mistake would be invisible — but any
   *future* operator carries its own digest written at fit time, so this exposure shrinks to zero
   on Oct 22.
3. **`data/donors/whitener_union.npz`'s axis digest is pinned but never exercised** — no build in
   this experiment used `--whitener-id union`. It is the same axis by construction (the prune drops
   basis rows, not columns) and by the same three-way computation.
4. **Two of item B's three panels carry essential genes** (30 and 34) against the real panel's 0,
   which `AGENTS.md` §8.1g established mid-run. §3.4 states every panel's count and redraws the
   measurement on one that has both October properties; the coverage numbers on 199's bundles
   should not be used to rank October sources.
5. **The item-B sizing used a stand-in base.** There is no D/E/F whitened base today. The base is a
   common positive denominator in both rules, so the carrying-vs-all-rows comparison and the
   above-zero claim are exact; the *absolute* fractions on the synthetic panels are indicative, not
   the doses that will ship.
6. **No member of the score was measured and none moved.** This experiment changes no shipped byte
   and prices no arm. It does not affect the ordering/direction ledger.
7. **`rule10_199.py`'s check [2]** (working tree vs `ab8aa3e`) was already failing for
   `scripts/make_submission_pm_cd4.py` before 235 — 222 and 229 committed changes to it after that
   revision. 235 adds to that drift. Its check [1], the meta-hash recomputation, is unaffected.

---

## 7. Two follow-ups, both outside 235's authorised files, both one move

1. **`scripts/def_pipeline.py`'s `whitener` stage skips the fit if the output already exists** —
   `if wout.exists() and wmeta.exists(): continue` — with no check on either file, so a stale
   operator left in the work directory by an earlier run is reused silently and the driver's own
   `sha256(wout) == ABC_OPERATOR_SHA` abort never runs. Since 235 the **build** catches it, but 20
   minutes later. The fix, in that skip branch: read `axis_sha256` out of `wmeta` and compare it
   against the axis recomputed from **this** bundle (`axis235.py:axis_for`, 15 lines, needs only
   the donors' `genes.txt`); refuse the skip if it differs or is absent. Recorded in
   `DEF-RUNBOOK` §2 as well.
2. **Fold `tests/test_october_hardening.py` into `scripts/test_make_submission_pm_cd4.py`**, so the
   guard's tests run with the rest of the builder's suite rather than only when someone remembers
   this directory. Nothing in them depends on living here.

`scripts/def_pin_panel_artifacts.py` needs **no** change (§2.3), which was checked rather than
assumed: its regex anchors on `WHITENER_REGISTRY: dict[str, tuple[str, str]] = {` and that line is
untouched.

---

## 8. Files

```
axis235.py                     the axis/n_cols/mis-pairing measurement          -> work/axis235.json
fire235.py                     known-answer row 5, the raise                    -> work/fire235.json, .log
carry235.py                    item B's coverage and sizing table               -> work/carry235.json, .log
essfree235.py                  section 3.4's disjoint AND essential-free redraw  -> work/essfree235.json, .log
tests/test_october_hardening.py 9 regression tests                              -> work/tests235.log
work/builds235.sh              rows 1 and 2 through 222's byte-exact build222.sh
work/refit1.log, refit2.log    row 3, both MATCH e2f65ecf...
work/make_submission_pm_cd4.PRE235.py  30fd9ab34d1e540e3329919bced4d4cfcfc19c81d34213d00b67289fd64532d0
work/build_union_whitener.PRE235.py    2337092687eaa973e2d141a869dfdb45969b4d4c48e9c0b166bb738c1eb1116e
```

The byte copies were taken **before the first edit** (git HEAD `1c30efe4`), because HEAD moves
under a running agent — that is how row 5's BLIND comparison is possible at all.

**The three built artifacts were deleted after hashing** (28.1 GB returned):
`submissions/235-ka-cd4f034-f034.h5ad`, `235-h1union-g028.h5ad` and `235-h1union-g028.vcc` were
proved byte-identical to 222's, which are still on disk. Their **`.meta.json` files are kept** —
they are the evidence for §4's field-by-field diff. Nothing was submitted, and no package here is
a candidate for a slot: both are re-builds of already-scored entries.
