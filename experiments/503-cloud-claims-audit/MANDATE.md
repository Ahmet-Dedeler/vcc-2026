# 503 — Audit the manuscript's claims against the experiment records (Cursor Cloud, Grok 4.7 xhigh)

Your experiment directory is `experiments/503-cloud-claims-audit/`.

## Why

The project's write-up (`experiments/294-writeup/PAPER.md`, claim table `experiments/294-writeup/CLAIMS.md`, ~72 claims
C-1 … C-72) will become the finalists' method description and possibly a paper. On 2026-09-26 alone several claims had
to be corrected after the fact: a "not bit-reproducible" claim (the difference was a file path), a coverage claim read
from a panel-only local copy (kolf2: "0" became 181–251), a "same line across labs ≈ −0.02" claim confounded by
mismatched perturbation sets, a control-noise statement, and a gate benefit reported with seed noise only (its
target-bootstrap interval spans 0). Find the rest before anything leaves the lab.

## Method

For every claim in `CLAIMS.md` (and every quantitative sentence in `PAPER.md` it supports): find the experiment record it
cites (`experiments/<NNN>*/RESULTS.md`, `docs/FINDINGS.md`, `docs/BOARD-LEDGER.md`); check that the number, its unit
(raw vs scaled board member; overall vs member; H1 vs board), its interval or noise statement, its sample size and its
scope (which panel, which targets, which surface) match the source; and check whether a later record corrected or
superseded it (grep the claim's key numbers and experiment IDs across `docs/` and `experiments/*/RESULTS.md`, newest
first). Flag claims that state an effect without the uncertainty unit the source used, and causal wording on
correlational evidence.

## Output

`RESULTS.md`: a table claim → source → verdict (SUPPORTED / OVERSTATED / STALE (superseded) / UNSOURCED / WRONG) with the
exact source line for each, and proposed replacement wording for every non-SUPPORTED claim. Do not edit the paper.

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
