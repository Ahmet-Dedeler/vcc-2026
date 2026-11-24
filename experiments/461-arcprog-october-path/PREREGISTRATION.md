# 461dr — the October path with Arc's denoised H1 programme as every context's block prior, vs 453dr
Preregistration, 2026-09-26 10:37Z: before 461dr is built, before 458d is built or scored, and before any 09-27…09-30
entry exists.

## The arm
- 461dr = 453dr with `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` and sameline_flags `pt` (not `pt:noprog`).
- `def_sameline.py --programme-source h1_2025_full` gives every context Arc's H1 programme from all public cells (458d's
  bank, `26662c52…`, the day's panel excluded). Context D (Jurkat) keeps its 0 same-line rows.
- The bulk is untouched, so the dump and the re-level should reproduce 453dr's.
- Expected config diff vs 453dr: exactly `sameline` (content: programme source and `prog` on). `bulk_ctx_scale` should
  be identical; a difference there is a finding to explain, not a failure.
- It is the October code path's known answer for the programme switch: the twin of 458d (vs 446b) through
  `def_day.sh`.

## When it goes
Only if 458d reads **positive** (fid ≥ +0.0015 vs 446b, 09-30 s1). Then 461dr confirms, before Oct 22, that the day's
switch carries the effect through the October path. With a partial, null or harmful 458d it is never submitted.

## Prediction and bins (461dr − 453dr, raw)
- pds and mse: identical to ≤ 1e-6 (falsifier: the bulk must not move).
- fid, nmae, reach, jac: the same direction as 458d − 446b, at 0.6–1.4× its size. The DE design differs only by CD4 in
  the bulk (453 ×1 vs 446b ×1.5), which the DE channel does not read.
- **Bins:**
  - fid in [0.6, 1.4] × 458d's Δfid, same sign: **the October programme path reproduces the board.** Set
    `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full` with flags `pt` for Oct 22 (DEF-RUNBOOK).
  - Same sign, outside the band: record the ratio. Adopt if Δfid ≥ +0.0010.
  - Opposite sign: a path defect. Keep the switch off until diagnosed.
- Name: a fresh joke at arm time. Description: empty.
