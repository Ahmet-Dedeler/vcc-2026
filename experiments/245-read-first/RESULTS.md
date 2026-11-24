# 245 — READ-FIRST: a primary-source map of the project

**Date:** 2026-09-04. Documentation only: `grep`, `sed`, reading, one `curl` of the public leaderboard, `vcc --json status`
and `vcc --json whoami` (read-only). **No `vcc submit`, no build, no `cell-eval2`, no model load, nothing written outside
`experiments/245-read-first/`.** Peak RSS < 0.2 GB (`uv run python` for the leaderboard rank computation).

## Produced

| file | what | size |
|---|---|---:|
| `READ-FIRST.md` | the eight-section map, every number carrying `(path:line)`; citation key at the top | **11881 bytes** (limit 12,000) |
| `work/contradictions.md` | 15 rows where `AGENTS.md` §8.x or the `docs/CURRENT_STATE.md` header states something a primary source contradicts or does not support, ranked by decision relevance, both sides shown, nothing resolved | 5301 bytes |
| `work/leaderboard-2026-09-04T1149Z.json` | the live public board at fetch time, 664 published entries | 808,548 bytes |
| `work/ranks.json` | our entry's six raw/scaled members, per-member world rank, top-5, leader gap, derived from that fetch | — |
| `work/vcc-status.log` | `vcc --json status` on our seven recorded entry ids, 11:57Z | — |

## Sources used, and two that sit outside the mandate's list

Cited: `docs/BOARD-LEDGER.md`; all eight `experiments/*/BOARD-RESULT*.md`; headline blocks via
`experiments/239-knowledge-audit/work/headlines.txt`, then the owning `RESULTS.md` for every quoted number;
`docs/DEF-RUNBOOK.md` §2–§5; `AGENTS.md` §2, §4 (incl. §4.1c's validity guard), §6, §7, §11.1, §14;
`scripts/def_pipeline.py` (`SHIPPED_ENV`, `SHIPPED_UNION_SOURCES`, `SHIPPED_FLAGS`, `STAGES`);
`scripts/make_submission_pm_cd4.py` (the union-source registry, for the `h1_2025` path only);
`submissions/SUBMITTED_LEDGER.json`; the live board.

Outside the list, cited only for what they alone establish:
- `AGENTS.md` §11.1b/§11.2 lines 1513, 1524-1528, 1574 — the `git add -A`, memory-gate and launcher rules the mandate
  required in §3; they are process rules, not numbers.
- `experiments/237-gate-validation/submit-armed.log` and `PREREGISTRATION.md`, and the `MANDATE.md` of 241/243/244/246 —
  the only primary evidence for §6 (what is armed, running, or mandated). Nothing from `docs/CURRENT_STATE.md`,
  `docs/STRATEGY*.md`, `docs/LOG.md` or `AGENTS.md` §8.x is cited as a source of a number.

## Three measured facts a reader should not miss

1. **A submission is armed.** `scripts/arm_submit.sh` (pid 6813, alive at 11:57Z) fires `submissions/229-stack6-g035.vcc`
   at 2026-09-05 00:00Z under 237's preregistration; 229 measured that package at −0.00219 and marked it NOT RECOMMENDED.
   `docs/CURRENT_STATE.md:83` names a different next candidate. Both sides are in `work/contradictions.md` row 1.
2. **`vcc status` resolves only the current entry.** The six earlier ids (`cd4f034`, `rc_callset`, `rc_hybrid`, `abc-knee`,
   a15 and its resubmission) return `not_found` on this account (`work/vcc-status.log`). The archive's `BOARD-RESULT.md`
   files are now the only record of those scores.
3. **Two published transfer constants disagree between primary sources**: the h1-test→board `reach` factor (1.290 in
   185's board result and the ledger; 0.769× in 221 on the real scorer) and the PLOLO→board cosine factor (1.384 vs
   1.000, 222's board result). Neither is resolved; 243 is running on the second.

Per-member world rank at 11:49Z (664 entries): `pds` 8, `mse` 6, `nmae` 56, `fid` 106, `reach` 7, `jac` 1; overall 3;
leader `cqawesome` 0.22431, gap +0.01248 (`work/ranks.json`).
