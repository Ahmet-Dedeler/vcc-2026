# 500 — Independent bug hunt: the bulk channel (pds, mse) (Cursor Cloud, Grok 4.7 extra-high)

You are an adversarial reviewer from a different model family. Your experiment directory is
`experiments/500-cloud-bulk-review/`.

## Why

The two bulk members (pds: cosine discrimination of pseudobulk effects; mse: jackknife-corrected squared error of the
log pseudobulk, capped by a spread-dependent credit) are where our model is strongest (pds rank ~3 on the board). They
come from: donor composition and whitening (`src/vcc2026/whiten.py`, `src/vcc2026/multi_donor.py`,
`src/vcc2026/donor_transfer.py`), the per-target magnitude key (410, `VCC_PERT_SCALE`), bulk sharpening
(`VCC_BULK_SHARPEN_Q`, 417), the stack dump/load (`VCC_BULK_LOAD`), CD4 deflation and the H1 union, the per-context
re-level (`VCC_BULK_CTX_SCALE`; `experiments/417-bulk-sharpen/ntclevel417.py`), and the ballast rows that render the
bulk while leaving the DE channel alone (`VCC_PERT_SCALE_BALLAST`, `ballast_rows` in
`scripts/make_submission_pm_cd4.py`). Earlier reviews tonight covered the DE emitter (484, in master) and the October
scripts (485); the only bulk-side finding so far is that `ntclevel417`'s Cp46 assumes 400-cell control constructs
(481 F2, being patched). Find bugs that silently cost pds or mse.

## What to look for

Gene-axis or target-order mismatches between donor tables and the panel; the panel-gene exclusion (pds drops all 300
panel genes; mse drops the target's own gene); log/linear space mix-ups (log1p of CP50k vs log2 CPM); renormalisation
after sharpening or scaling that changes what the comment says is preserved; the energy-preserving key and clip
boundaries; the stack dump being built with a different configuration than the one that loads it; context-scale factors
applied to the wrong context or twice; the ballast rows' composition (allocation `alloc(S, profile_bulk) - T_carriers`,
clipping, the trim to `k * 999,000`, the jackknife boost) against what pds/mse read; and the re-level's claim/cliff
arithmetic (`docs/reference/vcc2026-metrics-BRIEF.md` §2: mse's ρ credit).

## Method and output

Read the code and the named experiments' RESULTS.md; reproduce each suspected bug on synthetic inputs by executing the
production functions (the builder's unit tests: `experiments/410-pertscale-revisit/run_builder_tests.py`; tests that
need data you do not have, say so). `RESULTS.md`: a ranked list, each with file:line, the reproduction, the likely
member and direction, and a proposed fix (patch file in your directory).

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
