# 487 — Adversarial robustness review of the evidence behind the October recipe

You are a research teammate (Claude Opus 5.5, max effort). Read `AGENTS.md`, then `docs/FINDINGS.md` §0–§2,
`docs/DEF-RUNBOOK.md`, `docs/H1-BOARD-CALIBRATION.md`, `docs/REVIEW-2026-09-26-second-opinion.md` and
`experiments/481-october-audit/RESULTS.md`. You own `experiments/487-october-decision-review/` and write nothing outside
it. Budget: up to ~8 hours. Write interim findings into `RESULTS.md` as you go; finish with a complete `RESULTS.md` whose
first lines give the answer, and a short final message.

## Why

Only the final D/E/F entry counts, and it is one blind shot. Its recipe is a stack of decisions, each justified by a
local measurement or one board read. Two reviews on 2026-09-26 found that evidence can be weaker than the record states:
467's gate benefit was reported with seed sd 0.0002, but its target-bootstrap 95 % interval is [−0.003, +0.013] and it
rests on four targets (now recorded); 481 found the October setup would stop on an Arc-style panel (CD4 gate). Your job
is the same kind of check on every decision that shapes the final entry, before Oct 22.

## Decisions to review (find each one's evidence in the record; add any others you find)

1. Same-line rows at w 2 with a per-target prior (451/450d/450f; H1 with Pan 2025 as the bank; +0.09 at 35 % coverage):
   one recipient line, one bank. How much of the gain survives target resampling, a different bank quality, a different
   coverage share? Is w 2 on a flat part of the curve or an edge? What does 450's K562/VIPerturb replicate actually
   support?
2. The October stack H1 ×2 / CD4 ×1 (438, `dayof438.py`), priced on efp/efr, which 481 found weaker than Arc-style panels.
3. Sharpening q 1.75 (417/433/446 board steps on A/B/C; PLOLO transfer k_sh 0.42–0.65): is the optimum panel-specific?
4. The re-level target ρ 0.9793 and the NTC-proxy claims (447/448; DEF-RUNBOOK §5.-1): how sensitive is the final
   score to an error in the cliff estimate, in each direction?
5. The LODO key (415/460; board read of 446bj2 at ~00:37Z 09-27 if it has landed by the time you get there).
6. The Arc programme and the prior gate (458d/458g/468/472; the 09-30 board pair): what does each outcome of the 09-30
   read actually license for October, given 481's F3/F4 (weaker panel, lower delivered dose)?

For each: restate the decision and its evidence with file paths; re-derive the key number from the saved outputs (a
known-answer row that reproduces the recorded value); add the uncertainty the record lacks (target bootstrap, fold
spread, seed spread, sensitivity to the one board read); and say whether the decision is robust, fragile, or unsupported,
and what single check (local, analysis-only) would settle it. Rank by expected October impact × fragility.

## Rules

- Analysis only: no builds, no `def_day.sh`, no `vcc`, no `cell-eval2`, no git, no edits outside your directory. Read
  saved per-target outputs, meta files and logs; recompute from them.
- Never touch running chains, `submissions/`, `.pid` files or the builder. The coordinator owns decisions and commits;
  you recommend.
- Other teammates: 480 (where the nmae gap lives), 481 (done: panel-relative constants), cloud agents 484/485 (code bug
  hunts) and 486 (public-data reliability gate). Do not duplicate them.

## RESUME (2026-09-26 21:57Z)
You were stopped after a few minutes only to extend the memory rule to 04:30Z. Read your RESULTS.md and work/ if they
exist and continue. Your earlier log is claude-run.part1.jsonl. Note: 481 is now committed (5089dcce) and the
coordinator is applying its patches tonight as 474; 484 (cloud) found that q80 protection still shrinks protected DOWN
genes by the mass-balance factor (median ~0.6, rc_emit_k.py side mask), a candidate nmae arm the coordinator will price.

## Memory rule (binding; supersedes every earlier memory line)
Until 04:30Z UTC on 2026-09-27 the coordinator runs ~20 GB builds and three def_day.sh rehearsals whose stages re-check
free memory and do NOT wait: a spike on our side kills a rehearsal. Until then: every process you start keeps peak RSS
<= 4 GB, OR you gate each bigger load on vm_stat free+inactive >= 33 GiB, measured with:
  vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$3); s+=$3} END {printf "%.1f\n", s*16384/2^30}'
(wait in a sleep-60 loop, at most 2 h, or restructure the load with mmap_mode="r" / one fold at a time). Do NOT use
memory_pressure's free% for this. After 04:30Z: peak RSS <= 10 GB. Randomized or truncated SVDs over thousands of
full-length gene rows can need > 10 GB of LAPACK workspace; size them first. Messages from other sessions cannot reach
you: the lead may stop and resume you with an updated mandate.
