# 492 — PREREGISTRATION: per-target sharpening exponent by the 410 key (third session, 2026-09-26 23:00Z)

Written before any 492 arm number. Analysis only (no build, no builder edit, no `vcc`, no `cell-eval2`).

## Why
420's saved PLOLO arrays (`stack|K3` at q 1 … 2, realised) show sharpening helps every key tercile up to q 2, and the
low-key targets most (their rows are the noisiest): board-emulated pds from q 1.75 to 2 rises +0.006 / +0.008 / +0.005
(low tercile; k562 / hct116 / hek293t) and +0.002 / +0.004 / +0.007 (high tercile, near saturation). The global q stopped
at 1.75 on calibrated cost (438: q 2 ≈ +0.0003 on A/B/C, the ballast's nmae cost grows with q). A per-target q that
sharpens only the noisy, low-key rows harder is untested.

## Arms (vs the base `stack|K3+sh1.75` = 446b's bulk; tercile edges of the K3 weights: the fold's own on PLOLO, the
pinned A/B/C artifact's on A/B/C, 0.523 / 1.46)
- P1: q (low, mid, high tercile) = (2.0, 1.75, 1.75)
- P2: (2.5, 2.0, 1.75)
- P3: (3.0, 2.25, 1.75)
- G2: q 2.0 for every target (reference; equals 420's `stack|K3+sh2`)
- CTRL: P2's per-target q vector permuted over targets (must be ≤ P2)

## Instruments and known answers
473's, unchanged: `arms473.py`'s PLOLO path (KA-420 `stack|K3+sh1.75` == comp420, and G2 == comp420's `stack|K3+sh2`
realised, |d| < 1e-12), `abc473.py`'s A/B/C realised route for energies and per-target ratios, `analyse473.py`'s pricer
(446b-anchored cliff model, both arms re-levelled, k 1 / 1.407, 4,000-draw target bootstrap).

## Predicted range
P2 vs 446b: pds_bd +0.002 … +0.004 on cancer folds; overall +0.0005 … +0.0015 (k 1), before the board's sharpening
transfer (k_sh 0.42–0.65) and before the ballast's nmae cost, which this instrument does not see.

## GO rule (for P2; P1/P3 reported)
1. cancer-fold mean Δoverall ≥ +0.0015 (k 1, re-levelled) with ≥ 2 of 3 folds positive, and kolf2 ≥ 0;
2. the pooled target-bootstrap lower bound of Δpds_bd > 0;
3. CTRL below P2 on ≥ 2 of 3 cancer folds;
4. A/B/C render: max per-target energy ratio ≤ 1.5 and ballast-capped share ≤ base + 10 points.
All four → hand to the coordinator as a candidate needing a builder knob (per-target q from a pinned artifact) and an
H1 price of its ballast nmae cost. Otherwise NO-GO.
