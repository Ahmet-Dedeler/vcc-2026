# 540 — the LODO key's clip around 446bj2's (0.5, 3): preregistration

2026-09-27 ~01:45Z, third session. Written after the A/B/C energies (`abc540.py`, no truth), before any PLOLO number.

## Why
- 446bj2 (the LODO key at clip (0.5, 3) instead of 410's key at (0.3, 3)) gained +0.0022 overall on the board, almost all
  of it mse (−0.0145 raw).
- 492g's retrodiction of that step with the 446bj2-anchored pricer needs k_X 1.31 for mse: a composition-class
  transfer, like 417e's 1.3, not the 0.61 of the sharpening steps. The literal pds × 0.34 predicted its pds
  (+0.00024 vs +0.00021).
- Only floors 0.3 and 0.5 and caps 3 and ∞ were ever priced (415/420). The clip is one continuous knob of a class that
  transfers at full rate, and it has not been scanned.

## Arms (all at q 1.75 on the production stack, the same LODO key vector re-clipped)
Floors 0.4, 0.6, 0.7, 0.85, 1.0 at cap 3. Caps 2, 2.5, 4, 6 at floor 0.5. Combinations (0.7, 2.5) and (0.7, 4).
Base (0.5, 3) = 446bj2.

## Instrument
- PLOLO: `arms540.py`, lib420. The fold's LODO key is recomputed and checked against 418's `V1__PP`. Known answer:
  clip (0.5, 3) == comp420's `stack|V1_5_3+sh1.75`, |d| < 1e-12, or abort.
- Price: `analyse540.py`.
  - 446bj2-anchored cliff model; base and arm each re-levelled at their EV-optimal factors.
  - pds = 0.34 × literal PLOLO Δpds (cancer mean).
  - mse at k_X 1.31, with k 1.0 as the conservative end.
  - Ballast nmae cost by the K model and the displacement model (ballast492g's MDE-era slopes).

## Candidate rule (all must hold; otherwise the clip stays at (0.5, 3))
1. Net calibrated overall ≥ +0.0008 at k_X 1.31 AND ≥ +0.0004 at k 1.0. The worse of the two ballast models is used.
2. Cancer relative cos > 0 on 3/3 folds, and kolf2 relative cos ≥ 0.
3. A/B/C render: context energy moves ≤ 3 % (the re-level refunds it), max per-row ratio ≤ 3.
4. The best arm must beat its nearest neighbour on the grid by less than its own gain over the base. A lone spike on
   the grid is not a candidate.

If a clip passes, it becomes a key-rider candidate for the coordinator, one variable vs 446bj2. It would also be the
October question for `DEF_DAY_LODO=1` (the day-of fit takes `--lo/--hi`), to be re-checked on 438's October-like
panels before it is adopted there.
