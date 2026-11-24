You are a code-review agent for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md
first (§6 mechanics, §6.2 waiter rules, §11 method rules), then docs/DEF-RUNBOOK.md, docs/october-decisions.json,
docs/FINAL-PHASE-FEEDBACK-PLAN.md and experiments/519-cloud-october-review/RESULTS.md (the previous review round and its
five defects, all applied in 2d1fc214). This VM has NO gitignored data (no challenge matrices, donors or harness); you can
read code and write tests with small synthetic inputs.

Experiment 526. Write ONLY under experiments/526-cloud-october-review2/. Commit to your branch only. Never run `vcc`,
never touch GitHub issues/PRs, never contact anyone.

## Scope: everything on the October (D/E/F, Oct 22) path that changed since 2d1fc214
`git log 2d1fc214..HEAD` on: scripts/def_day.sh, scripts/def_day_helpers.py, scripts/make_submission_pm_cd4.py (the
builder), src/vcc2026/rc_emit_k.py (the emitter), scripts/october_recipe.py, scripts/def_postbuild.py,
experiments/524-october-recipe-dr/ (recipe_dr524.sh, gate524.py, chain1006.sh, arm1006.sh). The relevant commits:
8960c102 (496's uncalled-gene sign hook, VCC_DE_MDE_SIGN), 39d38196 (511's tercile tie rules in sharpen_q_by_tercile),
b1162b13 (502: DEF_DAY_DE_MAG / DEF_DAY_MDE_SIGN pass-throughs), 3eb21c57 (523: VCC_RC_CALLSET_CPM), 22cbca95 (524).
Also re-read, with fresh eyes, the parts of def_day.sh / def_pipeline.py / relevel_iter.py / def_sameline.py /
def_lodo_key.py that the Oct 22 run depends on.

## Questions (answer each with file:line evidence and, where possible, a failing synthetic test)
1. Oct 22 failure modes: a D/E/F bundle that differs from A/B/C in any way the code assumes (context labels, number of
   targets or cells per target, gene order, pert_counts.csv / gene_names.csv headers or ordering, ntc_id coverage, a
   context with no same-line bank, a context matched to a bank already in the donor mix, targets absent from every
   donor). Which of these fail closed, which fail silently, which produce a wrong but valid-looking artifact?
2. The resolver + pass-throughs: does `october_recipe.py command` produce exactly the env/args def_day.sh needs for
   every combination of DECIDED/PENDING switches (programme source on/off × prior gate × pct cap × S_P2 × zs × the MDE
   sign hook)? Any combination the builder refuses (e.g. the prior gate without the programme), or accepts but ignores?
3. Bit-identity claims: every new knob claims "unset = bit-identical". Check each by reading the code paths (not the
   claims): any default that changed, any RNG stream consumed differently, any config key recorded differently when
   unset (config keys matter because def_postbuild and the gates diff configs).
4. The 524 chain and gate: can recipe_dr524.sh or gate524.py pass a wrong artifact or refuse a right one? Is its
   builders() gate or its wait loop exposed to the process-name traps in AGENTS §6.2? Does chain1006's READY check
   match the lines arm1006.sh greps for?
5. Restart safety: if the machine reboots mid-def_day on Oct 22, what state is left (pinned builder, partial work dir,
   half-written artifacts), and does the runbook's recovery path restore it correctly?

## Deliverable
`experiments/526-cloud-october-review2/RESULTS.md`: numbered findings, each with severity (would break / would silently
degrade / cosmetic), evidence, a minimal reproduction (synthetic test in your dir where possible), and a proposed patch
(as a file in your dir, not applied). State explicitly what you checked and found correct, so the coordinator knows the
coverage. Do not modify any file outside your directory.
