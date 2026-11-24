# 495 — the LODO key at cap 6 on 446bj2 (preregistration, 04:42Z 09-27, before 495 is built)

## The arm
495 = 446bj2 with its stack dumped under `lodo415_lo0.5_hi6` (`work/pert_scale_415pp_lo0.5_hi6.npz`, `dc374281…`) instead
of `lodo415_lo0.5_hi3` (`8c82bef3`), re-levelled by ntclevel to ρ 0.9793.
- The weights are `def_lodo_key.py`'s rule at hi 6 on 8c82bef3's own key vector (540: the median normalisation is inside
  the clip, the energy constant outside).
- 49 targets rise above 3 and 10 reach 6.
- Expected config diff vs 446bj2: exactly `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale`.

## Why
The third session's 540 (preregistered grid: caps 2 / 2.5 / 4 / 6 at floor 0.5). Cap 6 is best on PLOLO and passes all
four rules: direction +2.6 % (3/3 cancer folds), kolf2 +8.3 %. The extension (8 / 10 / 15 / uncapped) does not move the
optimum.

## Prediction (495 − 446bj2, raw; 540's calibration at the key class the board just measured, k_X 1.31, literal pds)
- pds +0.0005 (β-emulated −0.0004);
- mse **−0.0140** (−0.0107 at k 1);
- nmae −0.0009 … +0.0006: the two ballast models disagree (mean K 8.83 → 7.54, fewer ballast cells, vs more displacement);
- fid, reach, jac pinned;
- **overall +0.0019 … +0.0029.**

## Bins (vs 446bj2)
- **Transfers:** overall ≥ +0.0012 and Δmse ≤ −0.006. October then passes `--hi 6` to `def_lodo_key.py` (a
  `DEF_DAY_LODO_HI` switch, after 540o's October-panel check).
- **Partial:** +0.0004 … +0.0012. Keep cap 3 unless a second arm confirms.
- **Fails:** below +0.0004, or Δmse > −0.002. Keep cap 3.
- **DE leak check:** |Δnmae| > 0.002 raw, or fid / jac / reach outside their pinned bands (0.001 / 0.001 / 0.01), means
  diagnose before reading.

## Slot
Not before 10-02. The day's rule is fixed before the previous day's entries exist. Name: a fresh joke. Description:
empty.
