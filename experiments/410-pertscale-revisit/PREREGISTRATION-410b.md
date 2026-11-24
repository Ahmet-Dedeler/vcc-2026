# 410b — UNWELDED per-target bulk magnitude (ballast cells), preregistered 2026-09-22 before its build

## One variable against 410
`VCC_PERT_SCALE_BALLAST=1` on 410's configuration (`VCC_PERT_SCALE=wshare_lo0.3_hi3` on 395). Every block is emitted as
400 − K **carrier** cells from the UNSCALED design — the DE channel exactly as 395 builds it (same profile, spread, key,
call set, CLEAN4 q80 prior, K=3 library classes) — plus K **ballast** cells of ≈995k counts whose per-gene counts put the
block's totals on the carriers' DELIVERED composition times the per-gene design ratio profile_scaled / profile_unscaled
("ratio mode"; genes where that would need negative ballast counts are clipped at the carriers' totals; clipped mass
returned to the ballast is trimmed proportionally so every cell stays ≤ 999,999). K per block = ceil of the ballast mass
the design ratio needs, clipped to 1…16.

Expected config diff vs 410's meta: exactly `pert_scale_ballast` (absent → {"lb": 995000, "kmax": 16}).

## Mechanism
The scorer's pseudobulk (pds, mse) sums raw counts over a block's cells, so ballast cells dominate it; the DE table
(Wilcoxon ranks on per-cell CPM, arithmetic-mean log2FC → fid, reach, jac, nmae) weighs every cell 1/400. 410's welded
arm scales the DE channel's fold changes with the bulk (median target 0.62×), which on H1 cost nmae +0.024 raw and reach
−0.016 whether or not the key assigned the weights (shuffled control: nmae +0.023). Ballast keeps the carriers at 1×.

## H1 evidence before the build (run410b.py; filled in from work/h1 when the ratio arms land)
- bal_ident (design-pinned, w = 1): nmae / fid / jac unchanged (|Δ| ≤ 7e-5); pds −0.0059, reach −0.0125, mse −0.0044 —
  pinning the bulk to the bare design removes the emitter's systematic low-count rounding pattern (emitbias410.py:
  −0.002…−0.0035 log units on 0.5–100-count genes). Ratio mode keeps that pattern. Ratio-mode rows: see RESULTS.md.

## Predictions vs 410 (the welded arm, same bulk)
- pds, mse: within the carriers' realisation noise of 410's values (same bulk design by construction; ratio mode).
- nmae, fid, jac, reach: back to 395's values within board noise (the DE channel is 395's again).
- Therefore vs 395: Δoverall ≈ 410's Δmse/6 + Δpds/6 (the bulk gain alone).

## Bins (primary = the DE members vs 395, secondary = the bulk members vs 410)
| DE members (nmae, fid, reach, jac) vs 395 | reading |
|---|---|
| all within noise (nmae ±0.003 raw, reach ±0.009 raw, fid ±0.003, jac ±0.001) | unwelding works on the board: bulk and DE channel are separately designable |
| nmae / reach recover > half of 410's loss but not all | partial: ballast cells leak into the DE table; size K and the ballast composition |
| no recovery | ballast does not unweld on A/B/C; close the construction |
| bulk members vs 410 | pds within ±0.003 and mse within ±0.005 raw → the construction preserves the bulk; larger → read the jackknife credit first |
Standing: never on rank; the board shows the latest entry.
