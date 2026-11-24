# INBOX — 🔴 TIME-SENSITIVE, your census is being built right now. Record BREADTH, not just overlap.

**The 300-target panel changes on 2026-10-22, and the replacement is disjoint from the current one.**
`DEF-RUNBOOK.md` §0 and §5 measured this on two synthetic panels, *"both 300 targets disjoint from
A/B/C"*:

| quantity | A/B/C today | realistic draw | **pessimistic draw** |
|---|---|---|---|
| K562-gwps coverage | 269 | 300 | **150** |
| CD4 union | 291 | 300 | **174** |
| `panel_support` (pm, cd4) | 248, 276 | 273, 286 | 141, **167** |

**Only the last entry after Oct 22 is prize-eligible.** So every source you are pricing now will be
scored against **300 targets it was not selected for**, and a source chosen because it happens to
overlap today's panel is worth **zero** in the only submission that counts.

## What to change — it is cheap and it is a census field, not a new experiment

For every source in the Phase 1 table, record **breadth alongside overlap**:

1. **Total distinct perturbed targets in the source** (genome-wide? a focused library? how many?).
2. **Its overlap with the union of all donor coverage**, not just with today's 300 — that union is
   what a new panel is drawn from (`def-bundle` / `def-bundle-r` build exactly these draws; reuse
   them rather than inventing a sampling model).
3. **Projected coverage under both `DEF-RUNBOOK` §5 draws**, realistic and pessimistic. This is the
   number that decides whether a source survives October.

## Why this changes your ranking, not just your notes

**A genome-wide source and a focused library can have identical overlap with today's panel and
completely different value in October.** The focused one is fitted to a panel that is about to be
replaced; the genome-wide one is not. **Rank by projected pessimistic-draw coverage, not by today's
overlap** — and where the two disagree, say so explicitly in `RESULTS.md`, because that disagreement
is the single most decision-relevant thing your census can produce.

Note K562-gwps is genome-wide and *still* falls 269 → 150 on the pessimistic draw, so breadth is a
strong mitigation, not immunity. Report the projection; do not assume any source is safe.

## And a second-order consequence worth stating

**This is also an argument for MORE sources rather than better ones.** A stack of many partially
covering sources degrades gracefully when the panel moves; a stack leaning on one well-matched
source does not. If your Phase 4 saturation curve is close to flat between, say, three and five
sources on today's panel, **the extra sources are still worth having** — they are October insurance.
Price that explicitly rather than dropping them for redundancy.

Everything else stands: the 60 GB disk budget (INBOX-1), the concentration test
(INBOX-CONCENTRATION), and the `h1union` baseline and PLOLO-units reporting rule (INBOX-BRANCH-A).
