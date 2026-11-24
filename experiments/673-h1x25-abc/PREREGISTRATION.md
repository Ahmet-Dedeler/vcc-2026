# 673 — W1: R0's exact recipe with the H1 dose at ×2.5 (def_day / A/B/C), a guard read for 708's Part A

Written 2026-10-20 03:32Z (by the clock) by the coordinator [3c7a7e], before any build, any D/E/F data and any W1 number.

## Why this slot
- 708's Part A may move the stack weights on D/E/F's own folds (438's grid, H1 ×{1 … 3} × CD4 ×{0.75 … 2}). The only weights
  board read on an H1 dose above ×2 is 710's h1s710 (+0.0004824, POSITIVE), on ALLSW's base (the 410-era bulk with CD4 ×1.5),
  not on R0's recipe (LODO percentile key, tercilec, θ 0.25, CD4 ×1). 708's own preregistration notes that dayof438's K3 table
  cannot see a key × weight interaction.
- W1 reads H1 ×2.5 on R0's exact recipe, so the guard 710 gives on ALLSW also holds, or fails, on the recipe that ships. It is one
  change (rule 5). CD4 is left at ×1: 438 found CD4 ×1.5 loses on every October-like panel, and 457dr already reads CD4 ×1.5 on
  the October path.
- It uses an expiring A/B/C slot (10-21 s1, or 10-21 s2) that nothing else needs.

## Arm
- W1 = R-ABC2's resolved recipe (the frozen decisions 04be2af3) + `DEF_DAY_H1_STACK=x2.5` (the dump subshell's
  VCC_H1_UNION_G = 0.70, h1s710's G). No hook. Built on def_day on A/B/C (180's def-bundle) like 695.
- Config diff vs R-ABC2 must be within {bulk_load_sha256 (REAL), bulk_ctx_scale}; the only changed export is DEF_DAY_H1_STACK.
  Frozen tree shas checked before and after.

## Prediction and bins (fixed now)
- pred = R-ABC2's board read 0.2739849 (nnHQfcZHsT3as9EafA3i) + 710's Δoverall vs ALLSW +0.0004824 = **0.2744673**.
- D = W1 − pred. **COMPOSES** |D| < 0.0005; **ABOVE** D ≥ +0.0005; **BELOW** D ≤ −0.0005.
- **BELOW blocks:** no 708 move may raise H1 above ×2 in the final (it extends 710's condition to R0's recipe). COMPOSES or
  ABOVE changes nothing: W1 can only block a move, never create one, and the October weights stay as 708 decides on D/E/F.
- Expected re-roll noise: same emission law, pass 1 redrawn (sd ≈ 0.00025 overall), so ±0.0005 is ~2 sd.

## Class and tries
- Class: transfer (composition weights), used here only as a guard. Tries for the H1 dose on the board: ×1, ×2 (many), ×2.5 on
  ALLSW (710), and this one.
- Reader: `read695.py --w1 <w1 entry> data/board-entries/nnHQfcZHsT3as9EafA3i.json` (to be added and dry-run before the read).
