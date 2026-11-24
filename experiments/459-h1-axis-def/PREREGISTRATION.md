# 459 — the H1 axis of the stack under October rules, on the board (preregistration, 2026-09-26 05:10Z)
Written before 459 is built and before any 09-27…09-30 entry exists.

## The arm
459 = 453 (the October path run blind on A/B/C) with ONE change: the stack dump at **H1 ×3** (`VCC_H1_UNION_G=0.84`) instead of
×2, with CD4 ×1. Built with `def_day_h1x3.sh`, a copy of `scripts/def_day.sh` (b517b6d2) whose only change is the dump's
export. Same bundle, q, flags and secant re-level; relabelled to A/B/C.
- Expected config diff vs 453: `bulk_load` + `bulk_load_sha256` + `bulk_ctx_scale`.

## Why
- The October stack's H1 dose is fixed at ×2 in `def_day.sh`. 438's day-of pricer picks the H1 dose per panel: ×3–×5 on
  the `real` draw (+0.008 over ×2), ×2 on `efp`/`efr`.
- On A/B/C under October rules (`abcdef`), 438 prices ×3 against ×2 at CD4 ×1 at only **+0.0004** (levelled k 1, q 1.25
  surface).
- The bulk members read exactly on the board, so even a small step calibrates 438's H1 axis. That is the axis whose
  panel dependence decides whether the October day should run 438's pricer instead of the fixed ×2.

## Prediction (459 − 453, overall board units)
- **+0.0004 (k 1) … +0.0005 (k 1.3); band −0.0005 … +0.0015.** Carried by pds and mse.
- The DE members are pinned up to the ballast's drift (nmae ≤ +0.0015 raw).

## Bins
- **≥ +0.0015:** 438 under-reads the H1 axis on A/B/C. Wire `dayof438.py`'s H1 choice into the October day (it would
  pick ×3–×5 on `real`-like panels).
- **−0.0005 … +0.0015:** as priced. Keep the fixed ×2 (lowest regret).
- **≤ −0.0005:** 438 over-reads H1. Keep ×2, with more reason.
- A DE member outside nmae +0.0015 / fid 0.001 / jac 0.001 / reach 0.01 (raw) means a leak. Diagnose.

## Slot
09-30 slot 2 (the rule is fixed today, before any 09-27 entry). Name: a fresh joke. Description: empty.
