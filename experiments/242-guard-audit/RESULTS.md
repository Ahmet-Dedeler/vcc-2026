# 242 — Guard audit: every check that can PASS while the thing it protects is WRONG

**Read-only audit. Nothing was fixed, built, submitted or re-measured.** Peak RSS of everything run
here: a `sha256` of one 3.7 GB file and three tiny `uv run python` predicate replays. 241 was live
throughout; no file outside `experiments/242-guard-audit/` was touched.

**Scope covered in full:** all 33 files in `scripts/` (14,355 lines), plus the two shared harnesses —
`experiments/162-retrodicting-validator/plolo.py` (**173** referencing `.py` files) and
`experiments/153-de-bench/` (**157**) — plus `experiments/064-pm-cd4-composition/run.py`, which
`plolo` and eight other experiments import for operator loading.

---

## 0. The headline, and it is not a null

**`scripts/` is NOT clean.** Twenty-one guards in the class were found. **Five of them sit on the
Oct-22 path**, and **three of those are in `def_postbuild.py` and the production builder — the two
files 240 and 235 hardened last week.**

Two results are *demonstrated*, not argued:

| | demonstration |
|---|---|
| **`check_submission_novelty.py` passes a package that has already been scored** | ran it: `check submissions/abc-knee.vcc` → **`OK: novel package, sha256 f5faa1e2…`, exit 0**, on the entry scored **0.1740437 on 2026-09-03** (`docs/BOARD-LEDGER.md` line 32) |
| **`def_pipeline.other_builders()` misses 3 of 6 realistic launch forms** | replayed its predicate verbatim: `MISSED` on two interpreter flags, on `runpy`, and on the un-exec'd `uv` shim |

And a correction to the mandate's framing: **there are five instances on the record today, not four.**
The fifth is already fixed and is the best-documented of them —
`scripts/cloud_score_worker.sh:43-52`, where the spot-interruption EXIT trap wrote `DONE 0` for a run
that produced no output, and *"a launcher that trusts the marker would have reported a passing
known-answer row from a run that never scored anything"* (run `ka-20260904`, 221). That comment is
the best statement of the class anyone here has written and it should be the template.

---

## 1. Ranked findings

Ranked by blast radius: Oct-22 path and shipped-artefact reach first, then shared-harness call-site
count, then everything else.

| # | guard | file:line | signature | blast radius | what it would fail to catch |
|---|---|---|---|---|---|
| **1** | the operator block of the D/E/F post-build falsifier | `scripts/def_postbuild.py:146-166` | **#3 + #1** | **Oct-22 gate**, the only prize-eligible entry | **any wrong operator except one specific file.** All 5 operator checks pass on `data/donors/whitener_union.npz` |
| **2** | `solve_cd4_fraction`'s zero-norm raise | `scripts/make_submission_pm_cd4.py:1509-1548` | **#3 (binary guard on a continuous quantity)** | production; the shipped `cd4f034` arm, **+0.0130557** | CD4 coverage anywhere in **(50%, ~92%)** — dose silently inflated; **1.4222× measured** at 58.7% |
| **3** | `check_submission_novelty.py check` | `scripts/check_submission_novelty.py:46-49, 80-82` | **#1 + #2** | **one submission slot** (the scarcest resource in the project) | **every package scored since 2026-09-02** — the ledger is untracked by git and 5 entries stale. **Demonstrated.** |
| **4** | the frozen-operator sha check | `experiments/064-pm-cd4-composition/run.py:179-183`, called from `plolo.py:309`, `plolo178.py:189`, `104/run.py:548`, `149/build_predictions.py:157` | **#3 (tautology)** | plolo: **173** referencing files | **nothing at all.** The caller passes `sha256_file(op_path)` as the expected value; the check is `x != x` |
| **5** | `plolo.Operator(path)` | `experiments/162-retrodicting-validator/plolo.py:330-336` | **#3 (no check exists)** | ~20 experiments (204, 207, 214, 223, 228, 234, 241…) | **the exact rule-10(b) accident**: 3 operators on disk, near-identical names, all `(17942, 17942)`; swapping them runs silently |
| **6** | `build_fold`'s cache-validity test | `plolo.py:244-251` | **#3** | same 173 files | a changed `FOLDS[held]["sources"]`, `["op"]`, `PM_BLEND` or `MIN_CPM`. Prints `cache hit` either way, **and restores the stale `audit` dict**, so the provenance lies too |
| **7** | "row count == n_constructs × cells_per_pert × contexts" | `scripts/def_postbuild.py:217-219` | **#3 (mislabelled)** | Oct-22 gate | **any row count.** The predicate counts non-empty contexts; `rows` is computed and used **only inside the message string** |
| **8** | `other_builders()` — the concurrent-build / OOM guard | `scripts/def_pipeline.py:195-207` | **#3** | Oct-22 path, ~19–22 GiB stages on a shared 64 GB box | 3 of 6 realistic launch forms (**measured**). A miss prints nothing — identical to "nothing running" |
| **9** | `classify() → RUNNING_NO_RESULT` | `scripts/agent_supervisor.sh:20-34, 45-48` | **#2 (silent default)** | every headless subagent | **an agent killed with no `result` event** — i.e. OOM, the failure AGENTS.md says this box produces *"without a traceback"*. Falls to `*) : ;;` and prints nothing |
| **10** | the `census` stage as a fail-closed pipeline step | `scripts/def_census.py:118-119, 182-215`; `def_pipeline.py:356-361` | **#3** | Oct-22 path | **a donor/CD4 coverage collapse.** `--verify` is passed only when `is_abc`, so on a D/E/F bundle the stage cannot exit non-zero. There is **no minimum-coverage threshold anywhere on the Oct-22 path** |
| **11** | the whitener npz↔json pairing check | `scripts/def_pin_panel_artifacts.py:147` + `build_union_whitener.py:404-418` | **#2** | Oct-22 re-pin | a mismatched pair whenever the fit was **unpruned** — `"sha256"` is written into the metadata only under `--basis-prune != none`, and `meta.get("sha256") is None` skips the check silently |
| **12** | `basis` metadata check | `make_submission_pm_cd4.py:1315` | **#3** | Oct-22 build | `basis[0]`. Only `basis[1]` (=17,942, constant) is asserted; **33,787 vs 29,394 is the sole metadata field separating `union` from `union_pruned`** |
| **13** | the on-axis preserve renormalisation | `plolo.py:340, 343` — `float(np.median(...)) or 1.0` | **#2** | 173 files | a zero median row norm — **the same >50%-zero-rows condition as the known live defect, one function away.** `before/after` silently becomes `1.0/after` |
| **14** | `endpoints()` sign accuracy | `plolo.py:203, 205` — `max(live.sum(), 1)`, `max(wt.sum(), 1e-300)` | **#2** | 173 files | an all-zero truth row-block: reads **sign accuracy 0.0**, indistinguishable from "every sign wrong" |
| **15** | `free_percent()` | `scripts/run_with_memory_guard.py:79-90` | **#2** | every guarded local scorer run | a `memory_pressure` failure or format change → returns `None` → the **system-free half of the guard disables itself**; only `system_free=?%` in the log distinguishes it |
| **16** | the builder's env defaults | `make_submission_pm_cd4.py:493` (`VCC_ALPHA` **1.0796**), `:633` (`CD4_TARGET_F` unset), `:663` (`H1_UNION_G` unset) | **#1 (stale reference)** | any build not driven by `def_pipeline.py` | **the builder's own defaults are the pre-`a11`, pre-`cd4f034`, pre-`h1union` configuration** — the exact `-0.0199671` package 239 §2.1 found. Only the driver + `def_postbuild` restore it |
| **17** | `solve_cd4_fraction`'s early return | `make_submission_pm_cd4.py:1538-1540` | **#2** | production | with `VCC_CD4_TARGET_F` unset (the builder default) the zero-norm raise **is never reached**; `CD4_FRACTION=0.05` is applied to an all-but-empty residual |
| **18** | the `whitener` stage's A/B/C-identity check | `def_pipeline.py:377-384` | **#3** | Oct-22 path | anything, whenever `wout`/`wmeta` already exist — **exactly the case the flag's own help text recommends** (`"pass an existing fit to skip the stage"`). Backstopped by 235's axis guard |
| **19** | `board_snapshot.sh`'s "did it parse" filter | `scripts/board_snapshot.sh:23` — `d.get('num_entries',0)>100` | **#2** | the board series 215 needs | a renamed API field → every snapshot deleted, `nohup … >/dev/null 2>&1` swallows it, and the symptom is identical to "the daemon isn't running". *(Currently healthy: key present, `num_entries` 664, pid 5342 alive.)* |
| **20** | `cushion_watch.osum()` | `scripts/cushion_watch.py:28` — `r.get(m) or 0` | **#2** | rival monitoring | a missing/renamed member silently contributes 0 to the ordering sum — can both mask and fabricate the `RIVAL CLOSING` flag |
| **21** | `check_harness_depth.py` | `scripts/check_harness_depth.py:63` | **#3, already on the record** | every DE number read off a thinned harness | the target-independent blob. **AGENTS.md §6.4 already states it: *"the thinned build passes at ratio 0.9999 and carries it anyway"*** |

**Not in the class, checked and clean:** `verify_expected_hashes` (`make_submission_pm_cd4.py:389-413`
— rejects omissions, additions *and* drift, and runs before any helper is imported);
`load_production_operator`'s 235 axis guard; `assert_build_carries_arms` (240);
`def_postbuild`'s `SHIPPED_ARMS` by-value rows (240); `def_preflight.py` (every check reads a value
from the bundle, and the frozen-table comparison fails safe); `cloud_score_worker.sh`'s `COMPLETED`
marker and its `resolved_de_backend`/`cell_eval2_version` gates (221); `073/run.py`'s
`load_frozen_operator`, which pins **both** shas from a manifest.

---

## 2. The top three, in detail

### 2.1 — `def_postbuild.py:146-166`: the operator block is a blacklist of one file plus two recomputations of itself

```python
op = meta["operator"]
new_panel = list(manifest["contexts"]) != ["A", "B", "C"]
chk("operator.sha256 is NOT the A/B/C operator (it is panel-dependent; rebuild it)",
    (op["sha256"] != ABC_OPERATOR_SHA) if new_panel else True, op["sha256"])      # 155
chk("operator.columns == the axis recomputed from THIS panel",
    op["columns"] == len(expect_cols), ...)                                        # 157
chk("operator.columns_sha256 == the axis recomputed from THIS panel",
    op["columns_sha256"] == expect_cols_sha, ...)                                  # 159
chk("operator.metadata.n_cols == operator.columns",
    op["metadata"].get("n_cols") == op["columns"], ...)                            # 162
chk("all_panel_coordinates_excluded", op["all_panel_coordinates_excluded"], "")    # 164
```

**Ask the mandate's question: what does each print if the operator were the wrong file?**

- `op["columns_sha256"]` is written by the builder at `make_submission_pm_cd4.py:4200` as
  `sha256(operator.columns)`, and `operator.columns` is
  `production_operator_columns(model, panel)` — **recomputed from the donor tables and the panel.**
  `def_postbuild` independently recomputes the same set from `data/donors/*/genes.txt`,
  `gene_names.csv` and `pert_counts.csv`. **Neither side reads the operator file.** So checks 157,
  159 and 162 are a cross-implementation consistency check on *the recomputation*, not on the
  artefact. `n_cols` is 17,942 on all four panels on disk (235). All three print PASS on any
  operator whose `n_cols` is 17,942 — which is all of them.
- Check 155 is a **blacklist with exactly one entry**: `e2f65ecf…` (`whitener_union_pruned.npz`).
  Point the build at `data/donors/whitener_union.npz` (sha `5608e807…`, the *unpruned* A/B/C
  operator, `n_cols` 17,942, identical column axis) and **all five checks PASS.**
- Check 164 is a property of the recomputed columns, not of the file.

**The field that carries the real answer is not read.** `grep -c axis_sha256 scripts/def_postbuild.py`
→ **0**. Same for `def_pipeline.py`. The builder writes `operator.axis_sha256` and
`operator.axis_sha256_checked_against` into the meta at `:4210-4211` precisely so this can be
audited, and the falsifier that exists to catch the silent-operator trap never looks at it. The
file's own docstring even says *"`columns_sha256` in the meta is closer but is only RECORDED, never
asserted"* — that sentence describes the pre-235 builder, and the same sentence is now true of
`def_postbuild` itself.

This matters because the script advertises independence: *"Everything asserted here is recomputed
from the bundle and the donor tables, so the check is independent of the build it is checking."*
For the operator it is not independent — it is **blind**, and it is the last gate before `vcc prep`
on a dark board.

**THE FIX (do not apply while 241 is running).** Add two rows to the operator block:

1. `chk("operator.axis_sha256 == axis_sha256(gene_names[recomputed columns])", …)` — compute the
   digest over the ordered **symbols** with the canonical form
   (`"\n".join(symbols).encode("utf-8")`, no trailing newline; `build_union_whitener.axis_sha256`)
   and compare it to `meta["operator"]["axis_sha256"]`.
2. `chk("operator.axis_sha256_checked_against is a real source", meta["operator"]["axis_sha256_checked_against"] in {"metadata:…", "LEGACY_OPERATOR_AXIS_SHA256[operator sha256]"})` —
   so a meta that carries no axis digest at all **fails** instead of skipping.
   *(Note: `submissions/222-h1union-g028.h5ad.meta.json` predates 235 and carries neither field,
   so the incumbent reference must be exempted by name or the fields read only from `--meta`.)*
3. Replace `op["sha256"] != ABC_OPERATOR_SHA` with a positive assertion: the operator sha must
   equal `WHITENER_REGISTRY[meta["config"]["whitener_id"]][0]`, and the metadata sha its pair.
   A blacklist of one is not a guard; a registry lookup is.
4. Assert `int(op["metadata"]["basis"][0]) == op["metadata"]["rows_kept_per_donor"]` summed, so the
   **row** count is pinned (finding 12).

**Known-answer requirement.** Three rows, all cheap, none needing a build:
- **KA-1 (pass):** `--meta submissions/222-h1union-g028.h5ad.meta.json` on the A/B/C bundle must
  still print the same pass/fail vector it prints today, with the two new rows waived-by-absence.
- **KA-2 (fire):** a copy of that meta with `operator.sha256`/`path` rewritten to
  `whitener_union.npz`'s `5608e807…` must **FAIL** on the new registry row — this is the case that
  currently prints 5/5 PASS, so the demonstration is the whole point.
- **KA-3 (fire):** a copy with `operator.axis_sha256` set to 235's measured `def-bundle-badaxis`
  digest `f2d6cc07…` (`experiments/235-october-hardening/work/fire235.log`) must **FAIL** the symbol
  digest while `columns_sha256` still passes — the two-symbol transposition 235 measured as having an
  *identical* index digest. That row is what proves the new check sees something the old three cannot.

---

### 2.2 — `solve_cd4_fraction`: the raise is a cliff at 50%, but the error is continuous, and it is measured

```python
added        = np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0)
residual_norm = float(np.median(np.linalg.norm(added, axis=1)))   # median over ALL 300 rows
base_norm     = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
f_per_unit    = residual_norm / base_norm if base_norm else 0.0            # :1525
...
if CD4_TARGET_F is None:                       # :1538 — the builder's own default
    return CD4_FRACTION, frame                 #        the raise below is never reached
if f_per_unit <= 0.0:                          # :1541
    raise ValueError("VCC_CD4_TARGET_F was set but the CD4 residual has zero norm; nothing to solve")
solved = CD4_TARGET_F / f_per_unit             # :1545
```

Production **does** raise where the PLOLO mirror silently returns 0 — AGENTS.md is right about that,
and it is the reason this is a different finding from the known live defect. **But the raise fires
only at the cliff.** The median runs over all 300 panel rows, so at coverage `c`:

- `c ≤ 0.5` → median is exactly 0 → **raise** (loud, correct);
- `c > 0.5` → the median lands at relative position `(c − 0.5)/c` **of the carrying rows' norm
  distribution**, not at its median. `f_per_unit` is therefore **understated**, and
  `solved = target_f / f_per_unit` is **over-dosed**, monotonically worse as `c` falls toward 0.5.
  Nothing checks it.

**This is not a hypothetical — production's own meta records the magnitude.** The sister source
`h1_union` runs through the *coverage-aware* `solve_carrying_row_fraction` and therefore reports
**both** rules side by side. From `experiments/229-stack4-build/work/229-stack6-g035-reg10.h5ad.meta.json`:

| quantity | value |
|---|---:|
| carrying rows | **176 / 300 = 58.67% coverage** |
| median row norm over **carrying** rows | 6.439850 |
| median row norm over **all** rows | 4.528221 |
| correct `solved_fraction` (carrying rule) | **0.2298584** |
| `shipped_median_rule_would_solve` (the CD4 rule) | **0.3268952** |
| **over-dose factor** | **1.4222×** |

The field is literally named `shipped_median_rule_would_solve` and its own comment says *"the number
this whole change exists to make visible."* **218 built the fix, 229 shipped it for `h1_union`, and
`make_submission_pm_cd4.py:1562-1563` states the exclusion in terms: *"This function is never applied
to CD4 — `apply_frozen_cd4_residual` still sizes CD4 through `solve_cd4_fraction` unchanged."***

**Why this is an Oct-22 finding and not an A/B/C one.** A/B/C CD4 coverage is
`fully_supported_targets = 276/300 = 92.0%`, where `(0.92−0.5)/0.92 = 0.457` ≈ the median, so
`f_per_unit = 1.0556645` is right and `solved_fraction = 0.3220720` is right. **The Oct-22 panel is
disjoint from A/B/C's by construction** (AGENTS.md §8.3), and 225 measured two registered sources
falling from 273 and 182 panel targets to **exactly zero** on it. The CD4 banks are re-extracted for
the new panel, so their coverage is whatever it is — and **the only downstream check is
`def_postbuild.py:182`, `cd4 complete targets > 0`.** A build at `cd4_complete = 160` (53% coverage)
prints PASS on every check in the repo and ships a CD4 arm dosed several-fold over the `f = 0.34`
that the board paid **+0.0130557** for.

**THE FIX (do not apply while 241 is running).**
1. Route CD4 through `solve_carrying_row_fraction` — the function is already in the file, already
   validated by 218's KA-4 (*"bit-for-bit on the three CD4 banks, max abs diff 0.0e+00, votes
   identical, 4/4 folds"*), already shipped for `h1_union`.
2. Independently of (1), add a **coverage floor** to `def_postbuild.py`: assert
   `cd4_complete_targets / n_constructs ≥ 0.85`, waivable only by an explicit
   `--allow-low-cd4-coverage` that prints the two doses side by side. Raise, do not warn.
3. Add the `f_per_unit`/`carrying` pair to the CD4 `residual_frame` **whichever rule is used**, so
   the discrepancy is always visible in the meta the way `h1_union`'s already is.

**Known-answer requirement.** The bar is exactness on the shipped path:
- **KA-1:** rebuilding `h1union` with (1) applied must reproduce
  `effects.cd4.residual_frame.solved_fraction = 0.32207201067578145` **to all 17 digits**, because at
  92% coverage the two rules coincide to within float noise. If it does not, (1) is not a
  no-op on A/B/C and must not ship. *(222's precedent: the rebuild through the modified builder
  reproduced `7036e89f…` byte-for-byte.)*
- **KA-2:** a synthetic residual masked to 176/300 carrying rows must reproduce 229's pair
  (**0.2298584** carrying, **0.3268952** all-rows, ratio **1.4222**) — that is the row that proves
  the new sizing is the one 218/229 validated and not a third rule.
- **KA-3 (fire):** `def_postbuild` on a meta with `cd4_complete_targets = 160` must **FAIL**; today
  it prints PASS.

---

### 2.3 — `check_submission_novelty.py`: the anti-re-submission guard is untracked and 5 entries stale. Demonstrated.

```python
def load() -> dict:
    if not LEDGER.exists():
        return {}                       # :47-49   a missing ledger makes everything novel
    return json.loads(LEDGER.read_text())
...
prev = ledger.get(digest)
if prev is None:
    print(f"OK: novel package, sha256 {digest[:16]}...")   # :81-83
    return 0
```

`LEDGER = submissions/SUBMITTED_LEDGER.json`. `.gitignore` line: `submissions/*` with a single
exception, `!submissions/.gitkeep`. `git ls-files submissions/` returns **`submissions/.gitkeep`, and
nothing else.**

> **The entire state of the guard that protects our scarcest resource lives in a gitignored
> directory.** A clone, a `git clean -xdf`, or any wipe of `submissions/` disarms it completely, and
> the disarmed output — `OK: novel package` — is **byte-identical** to the healthy output.

And it is already failing, without any of that. The ledger holds **5 entries, newest
`2026-09-02T00:32:04Z` (`a23`)**. `docs/BOARD-LEDGER.md` records at least **five** scored entries
since: `abc-knee` (09-03), `rc_hybrid` `nQUJZhZb3MW2uFJHBX5e` (09-03), `rc_callset`
`31oB9yWy8CSNBDEkY4iW` (09-03), `cd4f034` `QSe3JlNUYZ2BzvlcsL14` (09-04), `h1union`
`iUziqE1gewfcPKzYUDaB` (09-04). **None is in the ledger.** `record` is a manual second step that
nobody has run in two days and nine board points.

**The demonstration, run for this audit:**

```
$ uv run python scripts/check_submission_novelty.py check submissions/abc-knee.vcc
hashing submissions/abc-knee.vcc (3.69 GB)...
OK: novel package, sha256 f5faa1e21efdb25f...
EXIT=0
```

`abc-knee` scored **0.1740437** on 2026-09-03. The guard cleared it for a slot. This is the exact
failure it was written for — its docstring says *"Prose did not stop it, so this is the version that
can"* — and it currently cannot.

**Two independent defects, and the second is worse than the first.** Even fully maintained, the
ledger keys on the `.vcc` sha256. Ahmet's rule is *"never re-submit an already-scored **package**"*,
but the thing that must not be repurchased is the **configuration**: a rebuild at an identical
config with a different emission seed produces a different sha and reads as novel, while its board
result is already known to within the 230 noise floor.

**THE FIX (do not apply while 241 is running).**
1. **Move the ledger out of `submissions/`** — `docs/SUBMITTED_LEDGER.json`, tracked — or add
   `!submissions/SUBMITTED_LEDGER.json` to `.gitignore`. Untracked is the whole problem.
2. **Make `load()` fail closed:** a missing ledger must `raise SystemExit(2)` with "the ledger is
   missing; this check cannot run", never `return {}`. A guard with no state must refuse, not pass.
3. **Backfill the five missing entries** from `docs/BOARD-LEDGER.md` §"entry" rows (each carries the
   entry id and date) plus the packages on disk.
4. **Key on the configuration as well as the bytes.** Record
   `meta["config"]` + `meta["operator"]["sha256"]` + `meta["effects"]` digest alongside the package
   sha, and refuse on a *config* collision with a distinct, louder message. This is the half the
   original design missed and the half that actually implements the standing rule.
5. **Make `record` non-optional** by folding it into whatever wrapper runs `vcc submit`, so the
   ledger cannot silently stop growing again.

**Known-answer requirement.**
- **KA-1 (fire, and it is the whole point):** after the backfill,
  `check submissions/abc-knee.vcc` must **exit 1** and name entry `abc-knee`, board 0.1740437,
  2026-09-03. Today it exits 0. *That transition is the demonstration.*
- **KA-2 (fire):** with the ledger moved aside, `check` on any package must exit 2, not 0.
- **KA-3 (pass):** `check` on a genuinely new package must still exit 0 and print the same line, so
  the change is not a blanket refusal.

---

## 3. The demonstration queue — guards never observed firing

235 and 240 were each required to **show the raise**. That discipline has been applied to exactly
four assertions in this repo. Everything else is untested belief.

### Observed firing (evidence on disk — leave these alone)

| guard | evidence |
|---|---|
| 235's `REFUSING TO BUILD: … DIFFERENT PANEL` | `experiments/235-october-hardening/work/fire235.log` — `FIRE_def_pessimistic`, `FIRE_badaxis`, **and a `BLIND_pre235` control showing the old builder accepted the wrong panel** |
| 235's "operator carries no axis digest" | same log, `FIRE_unregistered_operator` |
| 240's `assert_build_carries_arms` | `experiments/240-driver-fix/` (RESULTS row 4) + `work/probe-pre240-*` |
| 240's `def_postbuild` `SHIPPED_ARMS` rows | 240 reproduced the pre-240 **14/14 PASS on a −0.0199671 package** verbatim |
| `verify_expected_hashes` drift | unit test `test_frozen_input_and_helper_hash_drift_fails_before_use` (`scripts/test_make_submission_pm_cd4.py:807`) |
| `VCC_CD4_TARGET_F` **range** validation | unit test at `:1864` — five bad values. **Note: this covers the `(0,1)` range only, NOT the `f_per_unit <= 0` raise at `:1541`** |
| `check_harness_depth` FAIL | `docs/HARNESS-DEPTH-CONFOUND.md`, observed on `orion-hct116.h5ad` at 0.725 |
| `cloud_score_worker.sh`'s `COMPLETED` marker | run `ka-20260904` — AWS reclaimed the instance 42 min in |

### Never observed firing — this is the queue, in priority order

| # | guard | file:line | why it is on the list |
|---|---|---|---|
| 1 | `solve_cd4_fraction`'s zero-norm raise | `make_submission_pm_cd4.py:1541` | **the shipped arm's sizing.** The only test near it checks the env range, not this branch |
| 2 | `solve_carrying_row_fraction`'s no-carrying-rows raise | `:1594` | the fix 218 built for exactly this failure, never shown firing |
| 3 | `union source … carries 0 of N panel targets` | `:1766` | 235 wrote it against a measured October failure (225: 273→0, 182→0) but did not demonstrate it |
| 4 | `--expect-contexts` / manifest disagreement | `:3895`, `:3901`, `:3923` | the round-identity guard; the runbook passes `--expect-contexts` **on purpose** on the day |
| 5 | `expected {n} perturbations / cells per perturbation (frozen table)` | `:3929`, `:3934` | the only defence against a bundle that is not the round you think it is |
| 6 | `PerturbMatch / complete CD4 mean3 support changed` | `:4181`, `:4186` | the `panel_support` pin; fires only on the `--repin-panel-support` path |
| 7 | `def_postbuild`'s five operator rows | `def_postbuild.py:155-166` | **finding 1** — and the demonstration is what shows they are blind |
| 8 | `def_postbuild`'s "row count" row | `:217-219` | **finding 7** — cannot fire; a demonstration attempt is the proof |
| 9 | `def_pipeline`'s ABC-identity check on the re-fit operator | `def_pipeline.py:383-384` | never fired; and finding 18 says the stage skip bypasses it |
| 10 | `other_builders()` / `free_gib()` refusals | `def_pipeline.py:397-402`, `:429-434` | **finding 8** — measured to miss 3 of 6 launch forms |
| 11 | `def_pin_panel_artifacts` post-write verification | `def_pin_panel_artifacts.py:184-186` | a substring search over the whole file; it cannot tell a digest under the *wrong key* from the right one |
| 12 | `def_pin_panel_artifacts` "records sha256 X but the npz hashes to Y" | `:147` | **finding 11** — skipped whenever the json has no `sha256`, i.e. every unpruned fit |
| 13 | `E064.load_frozen_operator`'s sha check | `064/run.py:179-183` | **finding 4** — at four call sites it is `x != x` and can never fire |
| 14 | `check_submission_novelty`'s REFUSED path | `check_submission_novelty.py:84-99` | **finding 3** — demonstrated *not* firing on a scored package |
| 15 | `agent_supervisor.sh`'s TRANSIENT relaunch | `agent_supervisor.sh:63-77` | has fired in practice (2026-09-03), but the `RUNNING_NO_RESULT` branch (finding 9) has never been exercised at all |

**The cheapest way to work this queue** is 235's own pattern, which is already in the repo and cost
it minutes: one script per guard that constructs the broken input, asserts `raised=True`, records the
**verbatim** message, **and carries a `BLIND_pre` row proving the previous code accepted it.** Rows
7, 8, 13 and 14 need no model load at all — they are meta-json and file-hash manipulations, seconds
each. Row 1 needs one masked residual array.

---

## 4. Two things worth saying about the class itself

**The `n_cols` failure mode recurs because the discriminating field is always the expensive one.**
Every instance here follows the same shape: the guard reads a field that is *cheap and stable*
(`n_cols`, `basis[1]`, a context count, a package sha, the recomputed column set) instead of the
field that actually varies with the thing being protected (the ordered symbol axis, `basis[0]`, the
row count, the configuration digest, the operator file's own contents). Finding 1 and finding 12 are
the *same bug as instance #2 in the mandate*, one field over, in the file written to fix it.

**And the recurring second-order failure is that the fix gets built and then not connected.**
`solve_carrying_row_fraction` exists, is validated bit-for-bit, and ships for `h1_union` — and CD4,
the arm worth +0.0130557, still uses the rule it replaces, by an explicit line of code. `axis_sha256`
is computed, asserted at build time and written into the meta — and the post-build falsifier never
reads it. That is the same pattern 239 named: **diagnosed, recorded, and walked past.** The
demonstration queue in §3 is the cheap instrument that would have caught all of it, because a guard
you have never seen fire is a guard you have never tested.

---

## 5. What was NOT done

No fix was applied. No file outside this directory was modified. Nothing was built, prepped,
submitted or scored. `cell-eval2` was not run. The only writes were this file and the directory that
holds it.
