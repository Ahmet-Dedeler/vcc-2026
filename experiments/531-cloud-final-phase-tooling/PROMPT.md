You are an engineering and operations agent for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master).
Read AGENTS.md first (§3 timeline, §6 mechanics and slot rules, §14), then docs/DEF-RUNBOOK.md, docs/FINAL-PHASE-FEEDBACK-PLAN.md,
docs/october-decisions.json, scripts/october_recipe.py, and experiments/526-cloud-october-review2/RESULTS.md (its question-1
table of D/E/F bundle surprises). This VM has NO gitignored data. Experiment 531: write ONLY under
experiments/531-cloud-final-phase-tooling/; commit to your branch only; never run `vcc`; never contact anyone.

## Tasks
1. **Reader.** `readdef.py`: takes two captured entries (JSON files in the format of data/board-entries/*.json, or entry
   ids) and a row id from FINAL-PHASE-FEEDBACK-PLAN §2, prints the raw and scaled member deltas and that row's bin verdict
   and adoption decision, exactly as the plan states it (bulk rows read exactly; DE rows respect the ~0.009 overall noise
   rule unless prior-only). Tests with synthetic JSON for every row, including the edge cases the plan names.
2. **Per-row commands.** For each row 2–12 of the plan: the exact `october_recipe.py` / `def_day.sh` invocation (env,
   arguments), the expected `def_postbuild.py --allow-config-diff` keys vs the entry it is compared with, and the gate. Where
   a row needs a knob that does not exist yet (ρ* scale, CD4 stack: being built in experiment 530), say so and write the
   command as it will be once 530 lands.
3. **Tabletop exercise, Oct 22 → Nov 5 23:59 UTC.** Walk DEF-RUNBOOK step by step as if it were Oct 22 and produce an
   hour-by-hour plan with gates and decision points, then stress it: a 2.8 GB upload that takes 7.5 h (observed 2026-09-29);
   a reboot mid-build (AGENTS §6.2 recovery order, 526 finding 3); a `vcc` CLI upgrade in the window; each bundle surprise
   in 526's table; members hidden vs returned for our D/E/F entries; the model-name queue and the reserved final name; the
   last-entry rule (the final upload must be the best-supported package, started by Nov 4 20:00 UTC). List every gap as a
   numbered finding (severity, evidence file:line, proposed fix).
Deliverable: `RESULTS.md` with the reader, the command table, the timeline and the findings, plus the code and tests.
