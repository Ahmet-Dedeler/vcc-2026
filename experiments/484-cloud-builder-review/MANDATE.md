# 484 — Independent bug hunt: the production builder (Cursor Cloud, Grok 4.7 extra-high)

You are an adversarial reviewer from a different model family. Your experiment directory is
`experiments/484-cloud-builder-review/`.

## Why

Every entry and the final October entry are written by `scripts/make_submission_pm_cd4.py` (7,605 lines; the builder),
with `scripts/build_submission.py` (the format writer) and modules under `src/vcc2026/`. It is heavily tested by
byte-for-byte known-answer rebuilds, but those only prove that a change did not change the output; they cannot catch a
bug that has been in the output all along. Find bugs that silently cost score.

## What to look for (correctness, not style)

- Gene, target or context misalignment: axis orders, panel-gene exclusion, the target's own gene, context mixing.
- Count-level integrity: integer rounding, per-cell totals near the 1,000,000 cap, explicit zeros in CSR, dtype
  overflow, negative counts, empty cells or genes, ties.
- The DE channel's emission (`rc_callset`, carriers, ballast rows, MDE add, zero repair, `VCC_DE_MAG`, sign priors, q80
  protection, the prior gate): code that contradicts its own comment or the metric spec
  (`docs/reference/vcc2026-metrics-BRIEF.md`), edge cases at thresholds, quantiles with ties, NaN or inf paths,
  seeds that make contexts or targets share draws unintentionally.
- Silent defaults: a missing artifact or env var that falls back to something without failing.
- Anything that would behave differently on a new panel or new cell line is teammate 481's territory; note it briefly,
  but spend your time on plain correctness.

## Method

Read the code; build small synthetic inputs (the environment can write a synthetic submission and run `cell-eval2` on
toy data) to reproduce each suspected bug; run the existing tests (`experiments/410-pertscale-revisit/run_builder_tests.py`,
`experiments/468-prior-gate/test_gate468.py`, `experiments/451-sameline-mode/test_sameline451.py`,
`experiments/417-bulk-sharpen/test_*.py`) and note which need data you do not have.

## Output

`RESULTS.md`: a ranked list, each item with file:line, what is wrong, the reproduction, the likely score impact, and a
proposed fix (patch file in your directory).

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
