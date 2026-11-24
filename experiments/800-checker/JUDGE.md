# The lab's judge (Phase B): how a proposal earns compute

Written 2026-10-15 ~23:10Z by the lab lead ([c61d4c]). The judge is a separate headless run per proposal
(`lab_subagent.sh 8NN-judge-<idea> experiments/800-checker/JUDGE.md` with the proposal path appended), never the idea's own
agent. It reads; it does not compute anything heavier than a lookup. It writes `experiments/<idea>/VERDICT.md` and stops.

## What a proposal must contain (PROPOSAL.md in the idea's dir), or it is returned unread
1. **Mechanism**: why the change should move the board, member by member, in one paragraph a reader can falsify.
2. **Class** (AGENTS 11.1c.1): transfer / scorer mechanics / A/B/C-specific. A/B/C-specific is never funded by the lab.
3. **Prediction before data**: the expected check800 price on BULK-W and/or DE-H1 / DE-CD4, member by member, with bins
   (SHIP / WEAK / NULL / HARM at the checker's rule) and the falsifier (the number that would kill it).
4. **Novelty against the record**: the closed families it is NOT (`docs/FINDINGS.md` §2 and the experiment the closure
   cites), and what is different about this construction. "A new name for the same construction does not count."
5. **What it costs**: cloud $, wall time, and how many configurations it will try (every try is charged by best-of-N).
6. **The day-of path**: what the October day needs to compute it from D/E/F's own inputs (controls, panel, public data),
   and nothing fitted to A/B/C.

## The judge's checks (each one answered in VERDICT.md, with the file and line that settles it)
- Is the mechanism closed already? Search FINDINGS §2, the experiment RESULTS it cites, and `grep -ri` the key words in
  `experiments/*/RESULTS.md`. If closed at its construction, KILL with the citation.
- Does the prediction have a falsifier and bins written before any number? If not, REVISE.
- Is it measurable on the checker's surfaces? If it moves members the surfaces cannot see (DESIGN §4), the proposal must
  say how those are pinned or priced; otherwise REVISE.
- Does it use any seen data as evidence of gain (its own training lines, A/B/C board reads, a sealed reserve)? KILL.
- Is the try budget honest? A family that plans > 10 configurations must say how it will choose without wearing BULK-W
  (nested folds, a held-out subset of folds, or one preregistered config).
- Is the expected gain worth its cost? Below the class noise (bulk direction: readable at ~0.0003; DE same-law re-roll
  0.00025; emission-law changes 0.0045), say so.

## Verdicts
- **GO**: the lab lead grants the cloud budget named in the proposal (≤ $50 per launch, charter caps).
- **REVISE**: list exactly what is missing; the idea agent fixes and resubmits once.
- **KILL**: with the citation; the idea is recorded as closed in the lab's log, with the mechanism.
