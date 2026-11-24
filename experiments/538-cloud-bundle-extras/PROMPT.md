You are an engineering agent for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md first
(§6, §6.2 waiter/process rules, §11 method rules, §14), then docs/DEF-RUNBOOK.md, docs/october-decisions.json,
scripts/october_recipe.py, scripts/def_day.sh, scripts/def_day_helpers.py, scripts/def_postbuild.py, scripts/def_census.py,
scripts/def_lodo_key.py, experiments/529-october-passthrough-bundle/PLAN.md, experiments/530-cloud-october-bundle/RESULTS.md
(and its patch_bundle530.py + test_bundle530.py), experiments/528-callprop-block/patch_defday528.py + test_defday528.py, and
experiments/531-cloud-final-phase-tooling/RESULTS.md (its review findings 12, 16, 17, 18 and 20).
This VM has NO gitignored data: write code and tests with synthetic inputs only.

Experiment 538. Write ONLY under experiments/538-cloud-bundle-extras/. Commit to your branch only. Never run `vcc`, never
contact anyone, never touch GitHub issues/PRs, never edit files outside your directory (scripts/, docs/, the builder and the
emitter stay untouched: you deliver PATCH SCRIPTS, not edits).

## Context
The October path (`scripts/def_day.sh`) builds the final D/E/F entry after the test bundle is released on Oct 22; the
leaderboard is dark then, so every knob must be default-off, byte-identical when unset, validated when set, recorded in the
final meta, and checked by `riders_match` and `experiments/524-october-recipe-dr/gate524.py`. Cloud agent 530 delivered five
knobs as `experiments/530-cloud-october-bundle/patch_bundle530.py` (tested 16/16 locally, `--check` clean, NOT applied: it is
applied after 2026-10-05 ~04:00Z, when a dress rehearsal that runs `def_day.sh` has finished). 528 wrote
`patch_defday528.py` (the call-propensity block pass-through). Several items from 531's review are still missing. Your patch
applies AFTER 530's bundle and 528's patch, on the resulting files.

## Task: `experiments/538-cloud-bundle-extras/patch_extras538.py` (+ tests), NOT applied
1. **`DEF_DAY_LODO_HI=<h>`** (531 #12): a fixed LODO key cap, passed to `scripts/def_lodo_key.py` as `--hi <h>`
   (`1 < h ≤ 20`). Needs `DEF_DAY_LODO=1`. Refuse it together with `DEF_DAY_LODO_HI_PCT` (one cap form only). Unset: today's
   behaviour exactly (cap 3 or the percentile form). Record it where the percentile form is recorded (the pinned key's
   name/json), add a resolver switch `lodo_hi` (PENDING) to `docs/october-decisions.json` via your patch, a mapping in
   `scripts/october_recipe.py`, a `riders_match` / `_october530_spec` extension and a gate524 value row.
2. **528's call-propensity block**: fold `experiments/528-callprop-block/patch_defday528.py` into the sequence (apply it as
   is if it applies cleanly after 530's bundle; otherwise re-anchor it in your own patch and prove equivalence with 528's
   tests). Add the resolver switch `callset_block` (PENDING), the october_recipe mapping, and a gate524 row
   (`config.rc_callset_block`: source, sha256, swap, steer).
3. **`scripts/def_postbuild.py --require-config-diff <k1,k2,…>`** (531 #17): fail unless EVERY listed key moved against
   the incumbent AND nothing outside the allow-list moved. `--allow-config-diff` keeps its meaning. Unset: identical output
   and exit codes to today.
4. **Census override** (531 #16): `DEF_DAY_CENSUS_MAX_ZERO=<n>` (default unset = today's floor of 5): raises the census
   stop threshold for zero-carrier targets to n (6 ≤ n ≤ 30), allowed only in a fresh work dir, logged as
   `census override: max zero-carrier <n> (count <c>)`, recorded in the day's env record. Unset: byte-identical.
5. **Repo path from the script** (531 #20): `def_day.sh` derives `REPO` from its own path instead of the hard-coded
   `~/Code/virtual-cell-challenge-2026`, with a guard that the derived path contains `scripts/def_day.sh`. On the
   production Mac the result is the same path; in a git worktree it is the worktree.
6. **One applier**: `apply_october538.py [--check] [--root <repo>]` runs, in order, 530's `patch_bundle530.py`, 528's
   `patch_defday528.py`, then your extras; `--check` reports what would change and exits non-zero on any anchor miss;
   idempotent (refuses to apply twice), atomic writes (temp + os.replace), prints the before/after sha256 of every file it
   touches. Style: `experiments/502-october-de-riders/patch_defday502.py` and 530's patch.
7. **Tests** `test_extras538.py`: on a temp copy of the repo's files (never the live ones), with stubs for anything that
   needs data: (a) every knob unset → the patched `def_day.sh` produces the same exported env and the same command lines as
   the 530+528-patched one (show how you prove "byte-identical when unset": e.g. run the script up to its first data access
   with a stub `uv`/python that records argv, as 530's tests do); (b) each knob set → exactly the intended export/argument
   change; (c) each refusal (ranges, `LODO_HI` with `LODO_HI_PCT`, `LODO_HI` without `LODO=1`, override in a reused dir);
   (d) `def_postbuild --require-config-diff` on synthetic metas: all moved → pass, one missing → fail, an extra → fail;
   (e) `october_recipe.py` resolves the new PENDING switches to no export; (f) `--check` on an unpatched tree is clean and a
   second apply is refused.

## Deliver
`experiments/538-cloud-bundle-extras/` with the patch scripts, the applier, the tests, and `RESULTS.md`: what each knob does
with file:line on the PATCHED copies, the test list with pass counts, the exact commands the coordinator runs on 10-05
(check, apply, test), and anything you found that is wrong or risky in 530's bundle or 528's patch (with evidence).
Commit to branch `cloud/538-bundle-extras`. Be exact; say plainly what you could not test without data.
