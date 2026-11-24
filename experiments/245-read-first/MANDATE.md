# 245 — READ-FIRST: a ≤ 12 KB map of the project built ONLY from primary sources, with a citation on every number

You are a documentation subagent. Work ONLY in `experiments/245-read-first/`. Read-only everywhere
else — **you do not edit `AGENTS.md`, `docs/CURRENT_STATE.md` or anything outside your directory**;
the lead will move your file into `docs/`. No `vcc submit`, no build, no `cell-eval2`, no model
loads; this is grep, sed and reading. `vcc --json status <entry-id>` (read-only) is allowed.

## Why this exists
`docs/CURRENT_STATE.md` is 297 KB and `AGENTS.md` is 123 KB. Experiment 239 (a fact audit)
found the archive — `experiments/*/RESULTS.md` and `*/BOARD-RESULT.md` — reliable, and the
summary layer at the head of those two files drifting: two red findings were in documents written
that same day. Every new agent is told to read the summaries first. **Produce the short document a
new agent should read instead, and build it from the archive, not from the summaries.**

## Sources you may cite (and the only ones)
`docs/BOARD-LEDGER.md`; every `experiments/*/BOARD-RESULT*.md`; the headline block of every
`experiments/*/RESULTS.md` (239 extracted them to `experiments/239-knowledge-audit/work/headlines.txt`
— use it, then open the RESULTS file for anything you quote); `docs/DEF-RUNBOOK.md` §2–§5 tables;
`AGENTS.md` §2 (hard numbers), §4 (metrics), §6 (mechanics), §7 (rules), §11.1 (method rules),
§14 (naming) — the factual sections, not the dated narrative amendments; `scripts/def_pipeline.py`
(`SHIPPED_ENV`); `submissions/SUBMITTED_LEDGER.json`; the live board via `vcc --json status`.
**Nothing from `docs/CURRENT_STATE.md`, `docs/STRATEGY*.md`, `docs/LOG.md` or `AGENTS.md` §8.x
may be cited as a source of a number.** If a fact exists only there, it is not established and does
not go in.

## Deliverable 1 — `READ-FIRST.md` (≤ 12,000 bytes; measure it)
Sections, in this order, each line carrying `(path:line)` for every number:
1. **Board position** — the current entry, its six raw and scaled members, world rank per member,
   the leader and the gap; the four 2026-09-03/04 arms with their Δ and the one variable each.
2. **The shipped configuration** — the exact environment line and flags the October build must
   carry, and the two arms not in the 09-03 recipe.
3. **Rules that bind every agent** — slots (2/day, reset 00:00Z, buy information never rank,
   one in flight), §14 (joke name, description EMPTY, nothing about method anywhere public),
   one variable per slot, preregister before build, price all six members, known-answer row in
   every screen, never `git add -A` with agents running, memory rules, `uv run python`,
   launch via `scripts/claude_subagent.sh`.
4. **Instruments** — a table: instrument → which members it reads → the board point it was
   validated on → known failures (include: h1-test is oracle-leaked for any arm carrying
   `data/donors/h1_2025`, i.e. every arm on top of `h1union`; the exact `nmae` identity has
   mis-signed the board on both checks; the PLOLO→board cosine factor is not identified — 1.384 vs
   1.000).
5. **Closed families** — one line each, the closing experiment id and its headline number.
6. **Open items** — what is running or preregistered, with owner experiment ids.
7. **October** — the eleven-stage sequence pointer, the guards added (235, 240), the retention
   lower bound (233), the essential-gene fact and the coverage/over-dose numbers (222/223/235).
8. **Where things are** — paths: donors, operators, packages, the two ledgers, the launcher.
End with the line: *"Facts checked on 2026-09-04 against the cited files; anything not cited here
is not established."*

## Deliverable 2 — `work/contradictions.md`
Every place where the head of `AGENTS.md` (§8.x) or `docs/CURRENT_STATE.md`'s live header states
a number or verdict that a cited primary source contradicts or does not support, ranked by
decision relevance, ≤ 15 rows: claim (path:line) · what the primary source says (path:line) ·
what depends on it. Do not resolve disagreements; list both sides.

## Rules
No adjectives, no narrative, no claim without a citation, no number from memory. If two primary
sources disagree, show both. Then write `RESULTS.md` (a short index of what you produced and the
byte count of READ-FIRST.md) and stop.
