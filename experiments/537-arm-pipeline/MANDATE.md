# 537 — The arm pipeline for 10-07 → 10-21: find, price and specify the next board arms (worker mandate)

Written 2026-09-30 11:10Z by the coordinator (an interactive Claude Code session). You are a headless research worker
(Claude Opus 5.5, max effort). Read `AGENTS.md` first (above all §4 scoring, §6.4–6.5 harness caveats, §11 method rules,
§14: nothing about our method leaves this machine), then `docs/FINDINGS.md` §0–§2, then `docs/NEXT-SLOT-DECISION.md` from
"2026-09-29 01:10Z" to the end, and `docs/october-decisions.json`.

## Why this matters
Ahmet funds this project to win the Virtual Cell Challenge: only the LAST final entry for contexts D/E/F (Nov 5 23:59Z)
counts, and the leaderboard goes dark on Oct 22 (final-phase members are hidden). Board slots before Oct 22 are
prize-irrelevant except as information. After the slots already armed through 10-06 there are ~30 slots (10-07 → 10-21,
two a day). Today the pipeline of good candidate arms for them is thin (446bj2zmz, 446bj2zcp's follow-ups, K20, the
all-switches production arm, 518's S1 contingency). Your job is to fill it with arms that are worth a slot, each priced,
specified and ranked by what it can add on D/E/F, not on A/B/C.

## Where we stand (verify against the files; do not trust this summary over them)
- Best entry: 446bj2z = 0.2728260 (`data/board-entries/be90smZUDcIV7j5MgQ03.json`): scaled pds 0.8152, mse 0.3804,
  nmae 0.1745, fid 0.0055, reach 0.2329, jac 0.0284.
- Leaders (board 2026-09-30 09:03Z, `data/leaderboard/2026-09-30T0903Z.json`): Illumina AI 0.4192 (pds .884, mse .669,
  nmae .342, fid .167, reach .393, jac .060), cqawesome 0.3276, Testh 0.3061. The coordinator traced their histories across
  our snapshots: all three climbed in many small board-driven steps (cqawesome's pds has been flat at 0.825 since 09-13; its
  gains are mse and DE). Board-tuned parameters cannot be re-tuned on a dark D/E/F board; our October path uses no A/B/C
  board information beyond global choices. Transfer is our edge; protect it.
- Board noise (AGENTS §11.1 rule 13): direction members (pds, mse) reproduce to 16 significant figures, so small bulk gains
  are readable; the ordering channel carries sd ≈ 0.0045 overall between digest-identical arms, so a DE gain must clear
  ~2 sd to be read. Price every candidate in those units.
- Reads already armed (each decides a switch or an arm; see NEXT-SLOT-DECISION and each experiment's PREREGISTRATION):
  10-01 460dr (October path + LODO day-of) and 446bj2zcp (528's call-propensity block, 535); 10-02 479's pair
  (per-target sharpen, `sharpen_pt`); 10-03 495 (LODO key cap 6, `lodo_hi_pct`); 10-04 498/497dr (combined bulk; October
  switches); 10-05 496a (CLEAN4 uncalled-gene hook, `mde_sign`), which also decides 446bj2zmz (517's Arc-sourced hook,
  527); 10-06 524dr (the October recipe dress rehearsal) and 446bj2zk16 (K 1,200) / 446bj2zk20 (K 1,600).

## Who else is working (context; do not interfere)
- The coordinator arms and reads slots, builds arms (`experiments/428-zero-gene-audit/run_arm_generic.sh`), commits. It
  reads your RESULTS.md; it cannot message you (you are headless), so put anything it must act on in RESULTS.md.
- `experiments/535-callprop-arm/rebuild535.sh` is building now (until ~12:00Z 09-30). Automated chains and slot waiters
  run all week; `recipe_dr524.sh` builds the October recipe on 10-05 from 02:00Z (~2 h). Other builds may start whenever a
  read lands. No other research worker is running. Peer sessions are idle.

## Your territory and hard rules
- Write ONLY in `experiments/537-arm-pipeline/`. No git, no `vcc`, no `cell-eval2` CLI, no builds, no edits to
  `scripts/`, `src/`, `docs/`, `data/`, `submissions/` or another experiment's directory. Reading anything is fine.
- Memory and concurrency: keep each process ≤ 8 GB RSS (a real-emitter H1 run is ~7.5 GB). Before any job above 2 GB,
  require `vm_stat` free+inactive ≥ 33 GiB AND no build running:
  `uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))'`
  must print 0. Never run heavy jobs 10-05 01:30Z → 04:30Z. Run heavy jobs serially; report peak RSS
  (`/usr/bin/time -l`).
- Process-name traps (AGENTS §6.2): never start a long-lived command whose arguments contain a chain's or a production
  script's name (`run_arm_generic.sh`, `def_day.sh`, `chain5NN…`, `build527.sh`, `rebuild535.sh`, `arm1001…`). Name your
  own drivers `drive537*.sh` / `*537.py`.
- Harness rules (AGENTS §6.5): h1-test is the only surface validated against the board; never put `data/donors/h1_2025`
  in a donor mix scored on h1-test; audit donor↔surface target overlap before every screen; Orion-thinned harnesses give
  relative direction only.
- Known answers first (rule 10): every screen carries a row whose answer is on record and aborts if it does not reproduce.
  For DE, reproduce one of 528's or 517's real-emitter H1 rows bit for bit (their drivers: `528/run528.py` +
  `drive528*.sh`, `517/…`, `435/run435.py`) before any new number. For the bulk, reproduce one row of 438's instrument or
  534's R6 panels (`534-generalization-audit/RESULTS-B.md`) first.

## Objective
1. **Inventory (RESULTS.md §1).** Every candidate arm or rider on record that is neither on the board nor closed
   (FINDINGS §2 is the closed list): deferred riders, contingencies, items marked untested / candidate / "only if", October-only
   riders that could also run on A/B/C, and combinations the scheduled reads open (for example the hook + propensity block:
   528 §6 priced it +0.0085 on H1 against 446bj2zmz's construction; the all-switches production arm). One row each: what it
   is, source `file:line`, channel (bulk / DE), best available price and instrument, status.
2. **Price (§3).** Price the top candidates on the validated instruments, all six members (or name the pinned ones and the
   mechanism that pins them), with CIs across seeds or folds. Convert to expected BOARD delta with the calibrations on record
   (bulk shape: kolf2's PLOLO fold predicted the board, cancer folds over-read ~2.5×, k_sh 0.4–0.65; DE: 528's 0.72 ×
   (fid + jac) and the other H1→board factors in FINDINGS), and to noise units (ordering sd 0.0045 overall).
3. **Transfer (§4).** For each candidate: does it use A/B/C board information? What do October-like panels (534 R6, 522,
   438) or a mechanism say about D/E/F? A candidate whose value is A/B/C-only goes to the bottom of the list, however large.
4. **Plan (§5).** A slot sequence for 10-07 → 10-21 as a decision tree on the 10-01 … 10-06 reads, with, for the first six
   arms: the exact env diff against 446bj2z (the `VCC_*` knobs the coordinator exports before `run_arm_generic.sh`), the
   expected config-diff keys, the predicted board delta per member, and preregistration-ready bins (the house style is in
   `535-callprop-arm/PREREGISTRATION.md` and `524-october-recipe-dr/PREREGISTRATION.md`). Prefer arms whose read is
   decisive: two slots on one preregistered dose curve beat two unrelated guesses.
5. **Stretch (§6).** If the inventory is thin, propose and price up to three NEW mechanisms (not a new name for a closed
   construction), each with its mechanism, its falsifier and a known-answer row.

## Deliverables and timing
- `RESULTS.md` with §1–§6 above, plus §0: a one-screen summary the coordinator can act on (top arms, their env diffs,
  expected deltas in noise units, bins). A negative (nothing worth a slot beyond what is scheduled) is a valid answer if it
  is measured.
- First complete ranked list (§0, §1, §5 even if some prices are still pending) by **2026-10-02 12:00Z**; keep pricing
  until 2026-10-05 18:00Z; then finish with a final message summarising the numbers.
- Take the time the work needs. Measured numbers over reasoning; say plainly what is unmeasured.
