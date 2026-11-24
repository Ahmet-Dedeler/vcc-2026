# 522 — October's per-target bulk key: LODO percentile cap vs the bounded strength key on October-like panels (worker mandate)

Written 2026-09-28 14:45Z by the coordinator. You are a research worker (Claude Opus 5.5, max effort). Read `AGENTS.md`,
then `docs/FINDINGS.md` §0–§2, then `experiments/540-lodo-clip-scan/RESULTS.md` (incl. 540o/540c), `543-quantile-cap`,
`546-october-cap-pertarget`, `541`/`542`/`544` (closed LODO variants), `518-bulk-mse-gap/RESULTS.md` (the strength key,
§3 and §6) and `438`'s October-like panel construction. Your directory: `experiments/522-october-bulk-key/` (write ONLY
there).

## Decision this informs
October's per-target bulk key is the switch `lodo_hi_pct` (`docs/october-decisions.json`, PENDING the 10-03 board read of
495 = cap 6 on A/B/C). The candidate October form is `def_lodo_key.py --hi-pct 96.64697522769474`, which equals cap 6 on
A/B/C but caps at 8–16 on 438's October-like panels (543), where it prices +11.6 … +19.4 % direction over cap 3 (546).
The board can only validate the A/B/C half (cap 6). Whether the heavier October tails are real signal or partly an
artefact of how the October-like panels were constructed (543's own caveat), and whether a BOUNDED estimator of the same
pot is more robust there, cannot be read on the board before Oct 22. 518 found an independent estimator of the same pot
— a donor response-strength key (noise-normalised count of strongly responding genes in each donor's own same-target row)
— that prices like cap 6 on A/B/C folds (+0.0016 … +0.0024 overall; +0.0001 … +0.0003 on top of cap 6) and whose
"src" form needs only the three production donors' panel rows (always on disk for any panel).

## Objective
On 438's October-like panels (efp, efr, real; their cancer folds, lib438, the October stack h2c1, q 1.75, the fold's own
LODO fit), price and stress-test the October key forms:
1. cap 3 (today's October default), fixed cap 6, the percentile cap (p96.647), uncapped;
2. 518's strength key alone (src form; θ as 518 fitted it, and a θ scan) with bounded weights;
3. combinations: LODO pct cap × strength (518's S1/S2 construction), LODO capped at 6 × strength;
4. robustness: (a) the tails' origin — for the targets the percentile cap weights above 6 on October-like panels, do the
   weights come from genuine cross-donor agreement or from panel-construction effects (e.g. targets whose donor rows are
   dominated by one donor, low coverage, essential-like responses)? Compare their held-out truth energies to their
   weights; (b) leave-one-panel-out: choose each form's parameters on two panels and score the third; (c) the A/B/C folds
   as the fourth surface (the known board calibration: shape steps k ≈ 0.4–0.65, composition/key steps ≈ 1.31).

## Known answers first
Reproduce 546's percentile-cap row and 540o's fixed-cap-6 row on at least one October-like panel exactly (|d| 0), and
518's strength-key A/B/C-fold price, before reading anything new.

## Deliverables
- `RESULTS.md`: known answers; a table per panel (efp, efr, real, and A/B/C folds) of relative cosine / literal pds per
  form with folds-positive counts; the tail-origin analysis; the leave-one-panel-out result; and a recommendation for
  October's key form with its expected October value and its worst case across panels.
- If the recommendation differs from "percentile cap": the exact change to `scripts/def_lodo_key.py` or a new day-of
  script (as a file in your dir, NOT applied), with a known-answer design (it must reproduce `dc374281`/495's cap-6
  artifact on A/B/C or state why it cannot), so the coordinator can schedule a rehearsal.

## Rules
- Memory: a production build queue runs until ~02:00Z 2026-09-29 (and 517/521 workers run alongside). Until 02:00Z keep
  EVERY process ≤ 4 GB RSS and check `vm_stat` free+inactive ≥ 33 GiB before any step expected above 2 GB (the PLOLO /
  lib438 folds peaked at ~9 GiB in 544: run those after 02:00Z, or at ≥ 33 GiB free+inactive, strictly one at a time).
  After 02:00Z, ≤ 10 GB per process, serial. Report peak RSS (`/usr/bin/time -l`).
- No builds, no `vcc`, no `cell-eval2` CLI, no git, no edits outside your directory. Never start long-lived shell
  commands whose arguments name production scripts (the build chains watch process lists).
- Known-answer row in every screen; target/fold counts on every claim; negative results are results. Finish with
  RESULTS.md and a final message summarising the numbers.
