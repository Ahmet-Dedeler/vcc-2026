# 541 — a robust LODO key (min / median / mean − 1 se / agreement-weighted) does not beat the mean: NEGATIVE

2026-09-27 04:50Z, third session. Exploratory, one script (`arms541.py`, PLOLO, lib420, the 540 path at q 1.75, realised
16). Peak RSS 5.4 GiB per fold, run at ≥ 29 GiB free+inactive.

## Question
540 found the LODO key's cap is the lever: cap 6 is best on the A/B/C folds, and larger caps lose on k562 (−3.5 %
uncapped). Could a robust aggregate of the per-donor leave-one-out projections tame the few noisy huge keys, and let
the cap rise further? Per target p and carrying donor j: v_j = ⟨s (b^{−j}_p + radd_p), d_{j,p}⟩ / |F_p|². 415's key is
their mean.

## Known answers
- The per-donor terms' mean == lib418.lodo_key to ≤ 4.3e-15 relative.
- The mean key reproduces comp420's `stack|V1_5_3+sh1.75` (cap 3) and 540b's cap 6 / uncapped rows: |d| 0, 4/4 folds.

## Result (relative cos vs the mean key at (0.5, 3); cancer mean (folds +), then per fold k562 / hct116 / hek293t / kolf2)
| key | cap 3 | cap 6 | uncapped |
|---|---|---|---|
| **mean (415 V1)** | 0 | **+2.62 % (3/3)**: +0.82 / +3.84 / +3.20 / +8.34 | +0.91 % (2/3): −3.52 / +4.12 / +2.12 / +11.39 |
| min_j (= mean − 1 se at n = 2) | −1.71 % (0/3) | +1.23 % (2/3) | −0.57 % (2/3) |
| median_j (= mean at n = 2) | 0 | +2.62 % (kolf2 +9.01) | +0.91 % (kolf2 +10.46) |
| mean × agreement (share of v_j > 0) | −0.22 % (0/3) | +2.43 % (3/3) | +0.73 % (2/3) |

- On the cancer folds every carried target has exactly 2 LODO terms (the fold's two remaining cancer donors), so the
  median is the mean and the min is mean − 1 se. With two terms no aggregate can tell a noisy key from a real one.
- The conservative aggregates lose at every cap. Down-weighting targets whose two donors disagree removes real signal
  along with the noise.
- kolf2 has 3 terms for 240 targets. There the median is slightly ahead of the mean at cap 6 (+9.0 vs +8.3 %). That is
  one fold, inside its noise.
- On A/B/C the artifact's targets mostly carry 3 donors (259 of 300; 35 carry 2, 6 carry 1), so production has 3
  terms where the cancer folds have 2. Only kolf2 tests n = 3, weakly. Nothing here supports changing the aggregate.

**Decision:** the mean LODO key stays. The lever is the cap (540: 6), not the aggregation.
