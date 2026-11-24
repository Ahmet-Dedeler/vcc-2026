You are an engineering agent for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md first
(§6, §6.2 waiter/process rules, §11 method rules, §14), then docs/DEF-RUNBOOK.md, docs/october-decisions.json,
scripts/october_recipe.py, scripts/def_day.sh, scripts/def_day_helpers.py, experiments/529-october-passthrough-bundle/PLAN.md,
and the two review rounds experiments/519-cloud-october-review/RESULTS.md and experiments/526-cloud-october-review2/RESULTS.md.
This VM has NO gitignored data: write code and tests with synthetic inputs only.

Experiment 530. Write ONLY under experiments/530-cloud-october-bundle/. Commit to your branch only. Never run `vcc`, never
contact anyone, never touch GitHub issues/PRs.

## Task: implement 529's October pass-through bundle as patch files (NOT applied)
Add these `def_day.sh` knobs, each DEFAULT OFF and byte-identical to today's behaviour when unset, each with a resolver
switch in `docs/october-decisions.json` (status PENDING), a mapping in `scripts/october_recipe.py`, a `riders_match` row in
`def_day_helpers.py` (the final meta must carry exactly what was exported), and a value-checking row in
`experiments/524-october-recipe-dr/gate524.py`:
1. `DEF_DAY_RHO_SCALE=<s>` (0.9 ≤ s ≤ 1.1): the re-level loop targets ρ* × s instead of ρ* (find where the re-level target
   0.9793 enters: def_day.sh's re-level loop and experiments/448-def-rehearsal-q175/relevel_iter.py).
2. `DEF_DAY_CD4_STACK=x1|x1.5`: the stack dump's CD4 multiplier (x1 is today's default; experiments/457-cd4-half-def/
   def_day_cd4x15.sh is the one-off variant that built x1.5 on 09-28: fold its difference in as the knob).
3. `DEF_DAY_RC_CALLSET_CPM=<n>` (100 ≤ n ≤ 5000): export `VCC_RC_CALLSET_CPM` on the level builds only (builder knob from 523,
   `scripts/make_submission_pm_cd4.py`; unset = 2200).
4. The Arc-sourced uncalled-gene hook: fold in `experiments/517-flex-programme/patch_defday517.py` and
   `patch_sameline517.py` (`DEF_DAY_MDE_SIGN=hookbank:<bank>:<tau>[:<tau_c4>[:le<N>]]`, `DEF_DAY_HOOK_ESTIMATOR=meanz_ne`).
   526 noted gate524's MDE-sign row parses only `<source>:<tau>`: add the hookbank form (the meta records source `file`).
5. `DEF_DAY_LODO_STRENGTH=<theta>`: 522's October rider (pct × exp(theta · z_src) on targets above the floor; tool
   `experiments/522-october-bulk-key/def_strength_key522.py`), needs DEF_DAY_LODO=1 and the percentile cap.
Rules for the patch: Python patch scripts with string-anchored, idempotent edits and atomic writes (temp + os.replace), in
the style of `experiments/502-october-de-riders/patch_defday502.py`; refuse to apply twice; `--check` mode. 526's fix
(def_day.sh refuses a builder that is not HEAD, a 6-line insertion after the `BSHA0=` line) may or may not be at HEAD when
you start: anchor your edits so they apply either way. Do not touch the builder or the emitter.

## Deliverables
- `patch_bundle530.py` (+ any per-file patch scripts), `test_bundle530.py`: unit tests on copies (bash -n; the resolver for
  every combination of the new switches; riders_match accepts the matching meta and rejects each mismatch; gate524 rows;
  the "unset == today's exports and arguments" proof, e.g. by running def_day.sh's argument/export stage against stubs).
- `RESULTS.md`: what each knob does (file:line), the bit-identity argument per knob, the tests and their output, and the
  coordinator's application procedure: apply after `recipe_dr524.sh` has finished its 10-05 run, then the inert KA (efp
  relaxed through def_day.sh must reproduce `ac49d319…` byte for byte, as 497/502 did) and one rider-on rehearsal.
