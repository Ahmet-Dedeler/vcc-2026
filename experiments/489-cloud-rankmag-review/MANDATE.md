# 489 — Independent review of the exact-U magnitude pass before it reaches the board (Cursor Cloud, Grok 4.7 xhigh)

Your experiment directory is `experiments/489-cloud-rankmag-review/`.

## Why

`VCC_DE_MAG=dn_k250_zs` (module `src/vcc2026/rank_mag.py`, called from `scripts/make_submission_pm_cd4.py`; history in
`experiments/423-rank-preserving-magnitude` and `experiments/429-zs-production`) moves counts between "carrier" cells
of a DE block so that each gene's arithmetic mean (read by `nmae`) changes while every Wilcoxon rank sum (read by `fid`,
`reach`, `jac`) stays exactly the same. It goes on the leaderboard in the 09-28 slot-2 riders (`446bz`, `446bj2z`,
`446br2z`, see `experiments/456-zs-riders/`) and may enter the final October entry. A 2026-09-26 review
(`experiments/484-cloud-builder-review` on branch `cloud/484-builder-review`) found that the q80-protected down genes it
receives are already shrunk by a mass-balance factor (median ~0.6), and that the pass targets 2.5× the lfc it is
handed. Find any bug in this module and its call site before the board read.

## What to check

- Does every move keep each cell's library, each gene's total, the ballast rows and every gene's U statistic exactly?
  Do the module's asserts actually test that, or could a path skip them?
- Interval logic: "each cell strictly inside its own interval between adjacent control values (evaluator order, count /
  (lib / 1e6))". Ties with control values, zero counts, cells whose library changes, float comparisons.
- Selector semantics (`dn`, `k250`, `zs`) against the comments and 423/429's RESULTS; what "reaches its target" means;
  behaviour when the target is infeasible; interaction with the zero repair and the MDE add that run before it.
- Determinism and seeds; the module sha pin in the builder.

## Method

Read the code and 423/429's RESULTS; build synthetic blocks (small AnnData or arrays) and run the module's real functions
on them; compute the scorer's U / p-values with the installed `cell-eval2` (or scipy's rank-sum with the same tie
handling, and say which) before and after the pass.

## Output

`RESULTS.md`: findings ranked, each with file:line, the reproduction, the likely score effect, and a proposed fix. State
plainly whether you found anything that should stop the 09-28 riders.

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
