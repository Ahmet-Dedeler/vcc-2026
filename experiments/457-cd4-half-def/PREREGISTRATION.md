# 457 — the CD4 half of the stack under October rules, on the board (preregistration, 2026-09-26 02:45Z)
Written before 457 is built and before any 09-27/28/29 entry exists.

## The arm
457 = 453 (the October path run blind on A/B/C) with ONE change: the stack dump at H1 ×2 + **CD4 ×1.5**
(`VCC_CD4_TARGET_G=0.5165482917974462`, 432 T6) instead of CD4 ×1. Built with `def_day_cd4x15.sh`, a copy of
`scripts/def_day.sh` (b517b6d2) whose only change is the dump's export. Same bundle (180's relabelled A/B/C), same
q 1.75, same flags, same secant re-level to ρ* 0.9793, relabelled to A/B/C.
- Expected config diff vs 453: `bulk_load` + `bulk_load_sha256` + `bulk_ctx_scale`. The CD4 dose lives in the dump,
  and the re-level follows the spread.

## Why
438 prices the CD4 ×1.5 half on A/B/C's panel under the October composition rules (`abcdef`) at +0.00172 (k = 1,
levelled, 3/3 folds). On three October-like panels it turns negative (−0.0011 … −0.0024), which is why October ships
CD4 ×1. 453 − 446b confounds the CD4 half with the level and the rules; 457 − 453 isolates it on the board. It
calibrates 438's instrument on the axis where its October call depends on panel dependence.

## Prediction (457 − 453, overall board units)
- **+0.0017 (k 1) … +0.0022 (k 1.3, 417e's realised stack rate); band +0.0005 … +0.0035.**
- Carried by pds and mse (bulk). The DE members are pinned up to the ballast's drift (K may rise with the extra spread;
  nmae ≤ +0.0015 raw).

## Bins
- **≥ +0.0010:** the CD4 half transfers as 438 prices it. The instrument is validated on the CD4 axis, and its October
  reading (×1 lowest-regret; ×1.5 loses on October-like panels) stands. October keeps CD4 ×1.
- **|Δ| < 0.0010:** flat on the board, so 438 over-reads CD4 even on A/B/C. October keeps CD4 ×1, with more reason.
- **≤ −0.0010:** ×1.5 loses even on A/B/C. October keeps CD4 ×1, and 438's A/B/C CD4 reading is wrong in sign:
  re-examine its CD4 residual rows before trusting its other axes.
- **≥ +0.0040 (far above 438's price):** 438 under-reads CD4. Re-price ×1.25 as an October hedge before Oct 22.
- A DE member outside nmae +0.0015 / fid 0.001 / jac 0.001 / reach 0.01 (raw): a leak into the carriers. Diagnose.

## Slot
09-29 (rule fixed after the 09-27 reads, before any 09-28 entry). Name: a fresh joke. Description: empty.
