# 485 — Independent bug hunt: the October day-of path (Cursor Cloud, Grok 4.7 extra-high)

You are an adversarial reviewer from a different model family. Your experiment directory is
`experiments/485-cloud-october-review/`.

## Why

Only the final entry on the new D/E/F lines counts (released Oct 22, no leaderboard, one shot). It is produced by
`scripts/def_day.sh` (183 lines, bash) driving `scripts/def_pipeline.py`, `scripts/def_sameline.py`,
`scripts/def_lodo_key.py`, `experiments/448-def-rehearsal-q175/relevel_iter.py` (the re-level secant loop), the builder,
and the preflight / census / audit / gate scripts they call. Read `docs/DEF-RUNBOOK.md` for the procedure and
`experiments/452-def-day-rehearsal`, `453-abc-dress-rehearsal`, `465-both-switches-efr` for rehearsal logs. Find bugs
that could break or silently degrade that one run.

## What to look for

- Shell correctness on macOS: `set -uo pipefail` without `-e` and unchecked return codes; grep/awk/sed behaviour under
  BSD tools and bash 3.2 (no associative arrays; `$` before `\|` in BSD grep); word splitting; quoting; exit paths that
  skip the builder restore; stale files or logs from a previous attempt satisfying a check (a log that quotes another
  log's END line has already bitten this project).
- Fail-open paths: a check that passes vacuously, an artifact whose absence silently changes behaviour, a default that
  applies without a log line.
- The re-level loop: termination, bracketing, what happens if a context never converges, overshoot handling, and
  whether the final artifact is always the one that passed the checks.
- Environment leakage between steps (an exported `VCC_*` reaching a build it should not), and context-label or
  context-count assumptions.
- Out of scope: constants tuned on A/B/C that change meaning on a new panel (teammate 481 audits those).

## Method

Read the scripts and the rehearsal logs; reproduce suspected bugs with minimal bash/python cases (emulate bash 3.2 and
BSD grep semantics where they matter, and say how you emulated them).

## Output

`RESULTS.md`: a ranked list with file:line, the failure scenario, the reproduction, the impact on October, and a proposed
fix (patch file in your directory).

## Ground rules for every cloud teammate

- You run on a Cursor Cloud VM with this repo and its `uv` environment (`uv run python ...`; `cell-eval2` is installed).
  The gitignored data (`data/challenge`, `data/donors`, `data/harness`, `data/external`, `submissions/`) is NOT on your
  VM; do not assume it exists. Production runs on the lead's macOS machine (bash 3.2, BSD grep, APFS); you are on Linux,
  so reason explicitly about any shell or tool behaviour that differs (AGENTS.md §6.2 lists the traps).
- Read `AGENTS.md` and `docs/FINDINGS.md` §0–§2 first. The record is large: grep `docs/` and `experiments/*/RESULTS.md`
  before claiming something is new.
- Write ONLY inside your experiment directory. Commit only that directory to your branch. Do not open a PR, do not merge,
  do not edit production files (`scripts/`, `src/`, `docs/`, other experiments); put proposed patches as separate files
  in your directory.
- Nothing public: no posting, issues, PRs, accounts, emails or messages. Never run `vcc`.
- Evidence standard: every claim carries a reproduction (a test, a script, a command and its output) that the lead can
  re-run. Label each finding CONFIRMED (reproduced) or SUSPECTED (reasoned). A well-supported negative is a result.
- Commit `RESULTS.md` early and update it as you go; finish with a complete `RESULTS.md` whose first lines give the
  answer. Budget: up to ~8 hours.
