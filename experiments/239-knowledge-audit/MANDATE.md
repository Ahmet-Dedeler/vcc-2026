# 239 — KNOWLEDGE AUDIT: find the facts we already measured and are not using. Today this failed twice.

You are a research subagent. Work in `experiments/239-knowledge-audit/`. **Read-only over the repo.
No `vcc`, no `cell-eval2`, no build, no submission, no production change.** Peak RSS ≤ 6 GB — this
is mostly reading. `uv run python`. Two other agents are live (221, 238).

## Why this exists — two failures today, same shape

1. **The essential-gene failure.** `docs/LOG.md` (~line 200) measured that **four independent
   essential-gene screens have ZERO overlap with the 300-target panel** against ~31–36 expected by
   chance each, concluded the panel was *"deliberately built to avoid common-essential genes —
   consistent with the stated selection criterion (>80% on-target knockdown with a measurable,
   non-lethal response)"*, and recorded those datasets as **"dead as direct donors."**
   **225, 233 and 236 then priced those exact four datasets as a +0.0066…+0.0104 October asset**,
   and I wrote all three mandates. **The knowledge was ours and four pieces of work walked past it.**
2. **219's audit** found **three re-opened closures** and its mandate then forbade it from measuring
   any. Two (208's Hamilton, 214's Arm A) were only picked up **today, a month later** — and both
   died, but both produced structural results worth more than the arms.

**Pattern: this project's own record contains decision-relevant facts that current work does not
reach.** With the current-panel search space exhausted (twelve routes closed 2026-09-04), **mining
what we already know is the highest-value activity available.**

## Scope — what to read
`docs/LOG.md`, `AGENTS.md`, `docs/CURRENT_STATE.md`, `docs/BOARD-LEDGER.md`, `docs/DEF-RUNBOOK.md`,
`docs/SLOT-PROGRAM.md`, `docs/STRATEGY-2026-09-04.md`, and **every `experiments/*/RESULTS.md`**
(there are ~200). Prior audits — **181, 198, 219** — did *re-pricing*. **You are doing something
different: a FACT audit.** Read theirs so you do not repeat them.

## What to find — four categories, in priority order

1. **🔴 CONTRADICTIONS.** A measurement in the record that contradicts a belief currently acted on.
   The essential-gene case is the template: a measured fact in the log versus three live mandates.
   **For each: quote both, name the file and line, and say which is right.**
2. **🔴 ORPHANED MEASUREMENTS.** Something measured, recorded, and never used by anything
   downstream — especially measurements made *for* one channel that were never tried in another.
   **This is how 224 found its positive** (189's `theta`, measured in ordering where it bought
   nothing, moved to direction) and how 234 then killed it. **Report the measurement, why it was
   made, and which channel never saw it.**
3. **STALE CLOSURES.** A route closed for a reason that no longer holds — a superseded constant, a
   corrected bar, an instrument since invalidated. **Today alone**: §8.2's direction bar was
   corrected **2.0–2.3×** (219), the `reach` toll was found **mis-keyed** (230), 220's ordering
   bound shrank **65%** (231), and 232's sign exchange rate was corrected **3.2×** (234). **Which
   closures rest on a number that has since moved?**
4. **UNRECONCILED DISAGREEMENTS.** Two experiments measuring the same quantity and disagreeing, with
   no adjudication. Known live one: the PLOLO→board cosine factor is **1.384** (`cd4f034`) vs
   **1.000** (`h1union`) on the matched functional, with a **7.84×** spread over four pairs (231).

## Method — and the discipline that makes this credible
- **Quote, with file and line.** A claim without a locatable source is not a finding.
- **Rank by decision relevance**, not by how interesting it is. *"Would this change something we are
  doing or about to do?"* If no, drop it.
- **Do not re-measure anything.** You are auditing the record. If a finding needs measurement, say
  what would measure it and stop.
- **Apply today's own corrections when you read** — `AGENTS.md` §8.1b (ordering noise **0.004502**,
  2.4× the floor), **§8.1d** (signs pay, magnitudes do not), **§8.1f** (the sign bar is **+0.01537**),
  **§8.1g** (the panel excludes essentials; synthetic panels must satisfy Arc's criterion), **§4.1c**
  (only invert `score_mse` > 0). **A finding that ignores these is noise.**
- **Report the null honestly.** *"The record is consistent and nothing decision-relevant is
  orphaned"* is a real and valuable result, and after today I doubt it — but do not manufacture
  findings to avoid it.

## Deliverable
`RESULTS.md`: a ranked table — **category, the fact, file:line for both sides, and what it changes** —
then the top three in detail with what would settle each. **Propose no slot. Build nothing.
Re-measure nothing.**
