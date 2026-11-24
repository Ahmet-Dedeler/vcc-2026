# 465 — the October day with BOTH new switches on a new panel (efr): PASS

2026-09-26 14:09:41Z → 15:52:38Z (1 h 43 min), `chain465.sh`, after 462l.
Run: `DEF_DAY_LODO=1 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full bash scripts/def_day.sh
experiments/238-essential-free-panel/work/ef-bundle-r … 465efr 1.75 auto pt 4`. It uses 238's essential-free realistic
draw (A/B/C controls relabelled D/E/F, real identity) and 452's efr day operator (`75773d59`, cloned; the whitener stage
verified it and skipped the fit).

## What the new code did on a panel that is not A/B/C's
- **LODO key** (`def_lodo_key.py`, rules def): carried 300/300 (all three donors on every target), median key 0.0130
  (A/B/C: 0.0318). w quantiles 0.5 / 0.5 / 0.5 / 1.0 / 2.20 / 3.0 / 3.0; 39 % at the floor, 16 % at the cap. Pinned as
  `lodo_lo0.5_hi3_def` (`f72f874f…`); the dump was scaled by it.
- **Programme source**: Arc's all-cells programme with efr's 6 overlapping targets excluded (294 rows, coverage 0.975)
  for D, E and F. D kept Nadig Jurkat's 9 covered rows (w auto). Mode ON for 3 contexts with `auto:pt`.
- **Re-level**:
  - lvl0 ρ 0.6732 / 0.6105 / 0.5820.
  - lvl1 overshot F (ρ 1.0925 > 1); that build was deleted as never shippable, and the secant loop recovered.
  - **CONVERGED after 3 levelled builds** at D 1.2198 / E 1.2842 / F 1.2496, ρ 0.9796 / 0.9796 / 0.9782, inside the
    [0.7, 1.4] bound.
  - 452's efr runs converged at similar factors (1.0553 / 1.1735 / 1.2498 when forced).
- **Audits**: zero-gene 0 in 0 blocks. The prep dry run verified 360,000 cells. The builder was restored (`4feae27d`),
  undoing both pins.

## Reading
Both switches compose and run end to end on a new panel within the day's time budget (1 h 43 min, like the
single-switch rehearsals). The October day can set either or both from their board reads without new code. Not
board-scorable: the panel is not A/B/C's.
