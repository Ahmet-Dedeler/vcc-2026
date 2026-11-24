# 518 — The bulk mse gap: where our error lives and the next transferable levers (worker mandate)

Written 2026-09-28 13:15Z by the coordinator. You are a research worker (Claude Opus 5.5, max effort). Read `AGENTS.md`,
then `docs/FINDINGS.md` §0–§2 (the bulk sections, the closed list, the transfer calibrations), then
`docs/H1-BOARD-CALIBRATION.md`. Your directory: `experiments/518-bulk-mse-gap/` (write ONLY there).

## Objective
The bulk members (pds, mse) reproduce on the board to 16 significant figures, so bulk gains are exactly readable. Our
best bulk (446bj2, raw) is pds 0.8674 / mse 0.6264; cqawesome is 0.8747 / 0.5704, Testh 0.8587 / 0.5723, Illumina
0.8823 / 0.4899. Testh has WORSE pds than us and much better mse, so part of the gap is not direction ranking.
Find where our bulk mse error lives and price new transferable levers.

## In flight — do not duplicate (their reads land 10-02 … 10-04)
495 (the LODO key at cap 6, +0.0025 priced), 479's per-target sharpening pair (P2g / S_P2, +0.001 each), 498 (the two
combined), 471c (the h1_2025_full union swap), the October percentile cap (543/546). The re-level (NTC-proxy ρ* 0.9793)
is production. q stops at 1.75.

## Step 1 — decompose (known answers first)
On the instruments that transfer for bulk (the PLOLO cancer folds through lib420 as 540–546 used them, the kolf2 fold
whose calibration predicted 417a, and H1 as a third surface), reproduce one recorded number from 540–546 or 420 exactly,
then split our bulk error into pots:
- per-target magnitude mis-scaling vs direction error (1 − cos²), energy-weighted as mse weights them;
- by target energy decile (mse is dominated by strong targets), by context, and by gene class (expression, panel genes,
  mean-response component);
- the value of a magnitude oracle (per-target optimal scale at our direction) and a direction oracle.
Say which pot is larger on each surface and how the pots map onto the board's mse and pds.

## Step 2 — levers in the larger pot
Price candidates with the class's calibration (shape steps transfer at k ≈ 0.4–0.65 of the cancer folds, composition steps
≈ 1.3× PLOLO; state which class a lever is and why). Examples of directions not on the closed list (check §2 first):
per-target magnitude from predicted response strength beyond the LODO key (e.g. DE-count or cross-donor dispersion as a
key), strong-target-specific treatment (the top energy decile), shrinkage of weak targets' bulk rows toward 0 when the
donors disagree (mse credit vs pds cost), a better mean-response (programme) component in the bulk if it pays on mse more
than it costs on pds under the ballast construction, the ρ* target itself (is 0.9793 the mse-optimal level on the folds?).

## Territory
- YOURS: the bulk (log-pseudobulk) prediction and its levelling.
- NOT yours: the DE channel (workers 516/517), production code, builds, slots (the coordinator).

## Deliverables
- `RESULTS.md`: known answer, the decomposition table per surface, every lever tried with its price in board units
  (with the calibration used and its range), negatives with mechanisms, and a verdict.
- For any lever priced ≥ +0.002 overall on at least two surfaces: a build spec (env knobs or a builder patch file in your
  dir, not applied; the expected config diff vs 446bj2; predicted pds/mse ranges) the coordinator can preregister.

## Rules
- Memory: a production build queue runs until ~02:00Z 2026-09-29. Until then keep EVERY process ≤ 4 GB RSS and check
  `vm_stat` free+inactive ≥ 33 GiB before any step expected above 2 GB (the PLOLO folds peaked at ~9 GiB in 544: run
  those only after 02:00Z or at ≥ 33 GiB free+inactive, one at a time). After 02:00Z, ≤ 10 GB per process, serial.
  Report peak RSS (`/usr/bin/time -l`).
- No builds, no `vcc`, no `cell-eval2`, no git, no edits outside your directory. Never start long-lived shell commands
  whose arguments name production scripts (the build chains watch process lists).
- Known-answer row in every screen; state every price's units; negative results are results. Finish with RESULTS.md and
  a final message summarising the numbers.
