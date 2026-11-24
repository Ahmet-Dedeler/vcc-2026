# 253 — READ-FIRST v2: the cited map rebuilt after today's twelve results

**Date:** 2026-09-04, 16:50–17:10Z. Documentation only: `grep`, `sed`, reading, one `curl` of the public leaderboard,
`vcc --json status` on seven entry ids (read-only), `uv run python` for the rank computation and the citation check.
**No `vcc submit`, no build, no `cell-eval2`, no model load; nothing written outside `experiments/253-read-first-v2/`.**
Peak RSS < 0.2 GB.

## Produced

| file | what | size |
|---|---|---:|
| `READ-FIRST.md` | the eight-section map, every number carrying `(path:line)`; citation key at the top; closing line as mandated | **11963 bytes** (limit 12,000) |
| `work/contradictions.md` | the status of 245's fifteen rows against today's head documents (which are resolved, by correction or by banner, and which are not), then 15 ranked rows, both sides shown, nothing resolved | 8,016 bytes |
| `work/leaderboard-2026-09-04T1653Z.json` | the public board at fetch time, 673 published entries | 821,952 bytes |
| `work/ranks.json` | our six raw/scaled members, per-member world rank, top-5 with model names, leader gap (`work/ranks253.py`) | — |
| `work/vcc-status.log` | `vcc --json status` on our seven recorded entry ids, 16:57Z | — |
| `work/arm-procs.log` | the two armed-submission processes as listed by `ps` at 16:58Z | — |
| `work/checkcites.py`, `work/cite.sh` | the citation checker (176 ranges in the final file resolve to existing lines, 0 bad) and the line printer used to verify every cited line by eye | — |

## What changed against 245's map

- **Board (§1):** rebuilt from a 16:53Z fetch. 673 published entries; we are still rank 3 at 0.21183096. The leader's new
  entry `m0904b` (`cqawesome`, 14:18Z) reads 0.22454757, gap +0.01272 (was 0.22431 / +0.01248). Per-member world ranks moved
  to `pds` 8, `mse` 7, `nmae` 58, `fid` 107, `reach` 8, `jac` 1 (were 8/6/56/106/7/1). A new 4th place (`Hongjian Zhou`,
  0.21092928, 15:33Z) sits 0.0009 behind us; it is in `work/ranks.json` and was cut from the map for bytes.
- **Configuration (§2):** the four inert builder flags, each with its commit, its byte-identity known answer and its October
  status: `VCC_CD4_SIZING`/`VCC_CD4_TARGET_G` (247, `a17b177`, wired on non-A/B/C bundles by the driver), `VCC_UNION2_*`
  (248, `2632b92`, no runbook row, driver-blind, its rider NOT BUILDABLE), `VCC_CD4_DEFLATE_K` (249, `be029c83`, never on
  D/E/F), `VCC_CD4_DEFLATE_BASIS=V5_abc` (251, `c321f76`, the October path if the `cd4svd5` point pays, driver-blind).
- **Instruments (§4):** 243's class factor (k_cos 1.36 re-dose / 0.94 new source) and the measured-cosine rule replace the
  "1.384 vs 1.000 not identified" row; the exact `nmae` identity keeps its two mis-signed board checks and is marked as
  under test by the 237 slot; 252 is named as owner of the PLOLO silent-zero fix.
- **Closed (§5):** the DLD-1 rider (248), 246's hypotheses 2–5, PerturbMatch's re-decision (241).
- **Open (§6):** both armed submissions with their pre-committed readings — 237's three-way Δ reading and 243's raw-`mse`
  bins for tonight's 00:00Z slot; 249's `pds`/overall predictions and five falsifiers for the second slot; 252 and 254 running.
- **October (§7):** RB §0.0's eight rehearsal corrections, §3.1's 241 and 250 updates, 247's guards.
- **Where (§8):** `SUBMITTED_LEDGER.json` now holds 10 entries (five backfilled by 247); the two armed packages are listed.

## Sources used

Cited: `docs/BOARD-LEDGER.md`; all eight `experiments/*/BOARD-RESULT*.md`; `RESULTS.md` files of 157, 162, 164, 166, 184,
193, 194, 195, 197, 203, 204, 205, 209, 211, 213, 214, 216, 218, 220, 221, 223, 224, 225, 226, 227, 228, 229, 230, 231,
232, 233, 234, 238, 239, 240, 241, 243, 244, 246, 247, 248, 249, 250, 251 (older ones via 239's `headlines.txt`, every quoted
line re-read in the file); `docs/DEF-RUNBOOK.md` §0.0, §1–§5; `AGENTS.md` §2, §4, §6, §7, §11.1, §14, plus §11.1b/§11.2 for
the `git add -A`, memory-gate and launcher rules (process rules, as 245 did); `scripts/def_pipeline.py`;
`scripts/make_submission_pm_cd4.py` (flag blocks and the two registries); `submissions/SUBMITTED_LEDGER.json`; the live
board. For §6 only: 237's and 249's `submit-armed.log`, `PREREGISTRATION.md`, `ADDENDUM-243.md`, and the first line of
252's and 254's `MANDATE.md`. Nothing from `docs/CURRENT_STATE.md`, `docs/STRATEGY*.md`, `docs/LOG.md` or `AGENTS.md` §8.x
sources a number in the map; §8.x and the CURRENT_STATE header are read only as the *claims* side of `work/contradictions.md`.

## Three facts a reader should not miss

1. **Two submissions are armed for 2026-09-05**, both processes alive at 16:58Z: pid 26557 fires `229-stack6-g035.vcc` at
   00:00Z (a gate test on a package 229 measured at −0.00219), and pid 80453 fires `249-cd4svd5-f034.vcc` after the first is
   terminal and 00:30Z (preregistered +0.0026…+0.0040). The readings are pre-committed in the cited files; read the first
   one's raw `mse` with the measured cosine, never the floor.
2. **The October deflation decision waits on the second point.** RB §3.1's 250 row says: if board raw `pds` rises against
   0.83572, ship the PINNED basis (`VCC_CD4_DEFLATE_BASIS=V5_abc`) on D/E/F and never a day-of re-fit. The builder path exists
   and is rehearsed (251), but the driver neither prints nor asserts that flag yet (R251 §7).
3. **Of 245's fifteen contradictions, six are resolved by correction, four only by a banner, and five stand**; one new row
   leads the list: `docs/CURRENT_STATE.md:89-90` still names the DLD-1 rider as the alternative for tomorrow, and 248 found it
   NOT BUILDABLE (+0.00216, 1.14× the gate, 25.2% leave-top-3-out).

Per-member world rank at 16:53Z (673 entries): `pds` 8, `mse` 7, `nmae` 58, `fid` 107, `reach` 8, `jac` 1; overall 3;
leader `cqawesome` (`m0904b`) 0.22454757, gap +0.01272 (`work/ranks.json`).
