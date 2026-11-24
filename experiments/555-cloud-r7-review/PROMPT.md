You are an adversarial reviewer for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md
first (§6.2 process rules, §11 method rules, especially §11.1c on sealed reserves), then docs/CONFIRMATION-LEDGER.md's
sections "Reserve R6" / "C10b" and "Reserve R7". This VM has NO gitignored data (no matrices, h5ads, operators or most npz
arrays): review code statically and test with synthetic inputs and stubs only.

Experiment 555. Write ONLY under experiments/555-cloud-r7-review/. Commit to your branch only. Never run `vcc`, never contact
anyone, never touch GitHub issues or PRs, never edit files outside your directory. You deliver findings, repro scripts and
(if needed) patch files that are NOT applied.

## Context
A sealed reserve is a test set nobody has tuned on; it is read once, under a claim preregistered before the first number.
R6 (four Arc-like bulk panels, experiment 534) has one access left, C10b on 2026-10-04 (experiments/551-r6-c10b,
`drive551.py access`). R7 (experiment 552) is the next reserve: R7a = four new Arc-like bulk panels T1–T4 (bundles under
experiments/552-r7-reserve/work/T*/bundle; target lists in work/draws.json), R7b = one cell-level DE surface (Calu-3,
floored; its access script is experiments/553-r7b-harness/access553.py). R7 confirms the FINAL October recipe once, after the
DE switches settle (~10-07…10-10). If any of this machinery is wrong, R7's single reading is wasted or, worse, misleading.

New since 2026-10-02 (commit 21381a77), all by the reserve's owner:
- `experiments/551-r6-c10b/r7a551.py` — R7a's driver (plan / ka-cd4 / record / access). It runs 551's instrument UNCHANGED on
  the T panels (outputs under 551/work/T1..T4, because guard551 confines writes to 551's dir), with FRESH operator fits (C10b
  re-fits R6's operators with --expect-sha256; R7a has no earlier values). It reuses gate551 by patching module globals
  (window_ok, thresholds, GATE_LOG, HEAVY_LOG) in-process.
- `experiments/551-r6-c10b/aggregate_r7a551.py` — R7a's aggregation (PRIMARY = F vs M5, the sum of the post-C10a steps; paired
  target bootstrap seed 552; retention vs R6 and A/B/C), modelled on `aggregate551.py`.
- `experiments/552-r7-reserve/CLAIM-R7-DRAFT.md` — the claim's rules (to be frozen after the DE reads).
- Records: `experiments/551-r6-c10b/work/r7a/plan.json`, `ka-cd4.json`, `test/` (an aggregation smoke test on efp), and
  `experiments/551-r6-c10b/RESULTS.md` (last section).

## Review these, and try hard to break each one
A. **r7a551.py against its models.** Compare `panel()` line by line with `drive551.py panel()` (C10b's) and with
   `experiments/534-generalization-audit/work/partB/drive534.py panel()` + `cd4_534.py` + `axes534.py` (how R6 was built).
   (1) Does any step differ from how 534 built R6 other than the intended ones (fresh fits without --expect-sha256; the axis
   check run in a child process; no --repro)? Flags, scratch dirs, --max-rss-gib, spec.json fields (`synthetic`, `meta`,
   `fold_ops`, `full_op`, `rules`), the order of steps, what is asserted after comp551. (2) The axis check: is
   `axes534.axis_sha` / `tuned_on_axes` / `ABC_AXES` / `DONOR_SETS` used correctly, does computing R6's axes in the child write
   anything (guard534 is imported there), and does the check really prove "own panel, not A/B/C's, not tuned-on, not R6's"?
   (3) The gate patch: list every read of `window_ok`, `thresholds`, `GATE_LOG`, `HEAVY_LOG`, `floor_gib`, `OVERRIDE`,
   `NO_START_AFTER`, `BLACKOUT`, `PHASE` in gate551.py and say whether the patch takes effect at each (module-global lookup at
   call time vs a value bound earlier). Can any path still apply C10b's window or the 10-06 12:00Z "no start" limit? (4)
   Resume and failure: a fit that dies mid-way, a crash between `heavy()` and the state update, operators deleted before every
   comp551 json exists, a re-run after `WINDOW CLOSED`, a second `record`. (5) Sealing: can `plan`, `ka-cd4` or `record` touch
   anything that reveals a T panel's truth (the PLOLO folds of donor data for T's targets)? Can `access` start without the
   claim's sha256 in the ledger? (6) The instrument pin (`instrument_files`): it imports lib551/designs551/price551/gate551 and
   exec's the top of full551/comp551/aggregate551 to list the repo modules they import. Find every module the instrument
   imports LAZILY (inside functions, at run time) in lib551, full551, comp551, price551 and their imports
   (experiments/1xx-4xx modules), and say whether each is covered by the 28 pinned files; if not, that is a finding.
B. **aggregate_r7a551.py.** Against `aggregate551.py` and `price551.py`: (1) the PRIMARY definition for each of the six
   configurations (`designs551.CONFIGS`), including F = M5 (falls back to cumulative M0 → M5) and F = M6; (2) whether the plain
   and the calibrated bootstrap calls see the SAME target draws (same seed, same consumption order) so that rows are paired; (3)
   `n_pos_folds`; (4) retention: the R6 reference (C10a's `results534B.json` steps; C10b's `results551.json` claim, used only
   if its recorded config equals F) and the A/B/C reference (`record534.json`); units on both sides; (5) anything that would
   make a TRANSFERS / FRAGILE / FAILS verdict depend on a choice made after the numbers. Write a synthetic test that builds
   fake Panel inputs (or monkeypatches price551.Panel) and checks the primary sums and the pairing.
C. **The claim's rules (CLAIM-R7-DRAFT.md),** against RESERVE.md (§4 access rule; amendments 1–3) and the ledger: (1) is R7a's
   PRIMARY the right single question for "the final recipe vs the recipe before the switches decided after C10a"? (2) R7b's
   board-likeness criterion (CP's raw jac > 0 and raw fid < 0) and the HARM veto (≤ −2.5 × band_total of the switch-on-recipe
   pair): calibrate them against 553's rehearsal numbers (`experiments/553-r7b-harness/work/draws553.json`, `band553_h1r7.json`,
   RESULTS "the band" and the five-draw table): false-veto rate per switch and for 2–3 switches, power against a K-size harm
   (−0.0017 overall-equivalent, between-draw sd 0.0023), and whether band_total or the between-draw spread is the right noise
   for a single surface. Propose concrete numbers if you disagree. (3) Every point where the draft leaves a choice that could
   be made after seeing an R7 number (that would break the reserve) — list them.
D. **R7b's one-command access** (`access553.py r7bd5000 --owner-go <claim> … --pairs … --ctrl-splits 5 --retro … --halves
   0,1,2,3,4`, command block in 553/RESULTS.md "What it needs"): (1) the `--owner-go` check (what exactly it verifies about the
   claim file; can it run against a claim whose sha256 is not in the ledger?); (2) what happens when a switch-on-recipe pair
   names an arm not emitted on the halves (the halves have "the 6 core arms at seed 0"); (3) band_total's construction in
   band553.py vs what the claim says; (4) the `--retro` ratios: units (raw members) and the reference (`access553_rehearsal5`
   is a five-draw mean). Static review plus stubs; the truth files are not on this VM.

## Output
`experiments/555-cloud-r7-review/RESULTS.md`: findings ranked BLOCKER / SHOULD-FIX / NOTE, each with file:line, the concrete
failure scenario (inputs → wrong output or crash), and the smallest fix. Then a section per task (A–D) with what you checked
and how, including checks that found nothing. Put synthetic tests under `repro/` with the command and its output. Be exact:
quote code, cite lines, and separate what you ran from what you inferred.
