# 410c — the UNCAPPED per-target key through ballast, preregistered 2026-09-22 ~19:00Z (before 410b reads, before its build ends)

One variable against 410b: `VCC_PERT_SCALE=wshare_lo0.3_inf` (artifact 94c048e9…, clip(0.3, ∞)) instead of `wshare_lo0.3_hi3`;
ballast (`VCC_PERT_SCALE_BALLAST=1`) and everything else unchanged. Expected config diff vs 410b's meta: exactly
`pert_scale` (hi3 → inf) and `pert_scale_sha256`.

Evidence: PLOLO relative aggregate cosine +34…+39 % uncapped vs +27…+30 % at hi 3 (three cancer folds); H1 through the
shipped ballast construction: mse −0.0913 vs −0.0496 raw, pds +0.0007 vs +0.0056, DE members pinned in both
(board-slope +0.0159 vs +0.0097). Board: 410 (welded, hi 3) moved raw mse −0.0797 — the key transfers.

Predictions vs 410b: raw mse lower by a further 0.02…0.06 (the uncapped/capped ratio on PLOLO and H1 is 1.2…1.8×);
pds within ±0.004 raw of 410b (H1 says slightly lower; 410's weld read −0.0038 vs 395); DE members within noise of 410b.
Bins (Δ overall vs 410b): ≥ +0.003 → dose confirmed, uncapped becomes the base; −0.003…+0.003 → flat, keep the capped
key (smaller DE/ballast footprint); ≤ −0.003 → the tail over-weights on A/B/C; stay capped.
Gate before upload: GO-410c written only if 410b shows the DE members back within noise of 395 (the unwelding holds).

## GO rule for 410c, fixed 19:10Z before 410b reads (Δ are 410b − 395, raw)
GO-410c if ALL: |Δnmae| ≤ 0.010, |Δreach| ≤ 0.012, |Δfid| ≤ 0.003, |Δjac| ≤ 0.001 (the DE channel is back within board
noise — ballast unwelds on A/B/C) AND Δmse ≤ −0.050 (the bulk kept most of 410's −0.080). Otherwise NOGO-410c: the slot
goes to diagnosing the leak (or expires); 410c is not uploaded on a construction that leaks.
