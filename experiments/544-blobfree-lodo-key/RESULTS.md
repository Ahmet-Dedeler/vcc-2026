# 544 — a "blob-free" LODO key (each donor's mean-response direction projected out): hypothesis refuted, gain inside noise

2026-09-27 05:05Z, third session. Exploratory, `arms544.py` (PLOLO, lib420, 541's path). Peak RSS 9.2–9.4 GiB per fold
(it reads each donor's full delta table to form its mean direction), run at ≥ 26 GiB free+inactive after 04:30Z,
inside the ≤ 10 GB rule.

## Hypothesis
- In 540's per-fold detail, only the held-out-K562 fold loses at caps above 6. There the key's two terms are both
  Orion-vs-Orion (hct116, hek293t).
- The Orion lines share a technical component: the mean directions have cos 0.39, and HEK293T's panel rows carry a
  median 7–8 % of their energy along its own mean direction. That component could inflate Orion-vs-Orion terms.
- Test: project each donor's mean-response direction u_j (the mean of all its delta rows, official axis, panel columns
  zeroed) out of both the leave-one-out prediction and the donor row before each term is taken.

## Result (relative cos vs the LODO key at (0.5, 3); per fold k562 / hct116 / hek293t / kolf2)
| key, cap | cancer mean (folds +) | per fold |
|---|---|---|
| LODO, 3 (446bj2) | 0 | 0 / 0 / 0 / 0 |
| LODO, 6 (495) | +2.62 % (3/3) | +0.82 / +3.84 / +3.20 / +8.34 |
| LODO, uncapped | +0.91 % (2/3) | **−3.52** / +4.12 / +2.12 / +11.39 |
| blob-free, 3 | +0.55 % (2/3) | +1.23 / +0.71 / −0.29 / −0.98 |
| blob-free, 6 | +2.99 % (3/3) | +1.75 / +4.16 / +3.05 / +7.67 |
| blob-free, uncapped | +1.08 % (2/3) | **−3.33** / +4.48 / +2.11 / +11.79 |
| blob-free, p96.647 | +2.88 % (3/3) | +1.53 / +4.17 / +2.95 / +7.49 |

- **Refuted.** The held-out-K562 fold still loses uncapped (−3.3 % vs −3.5 %). The shared component does not explain
  why high caps fail there.
- Mechanism: the prediction rows carry essentially nothing along u_j (median energy share 0.000). Projecting u_j out of
  the donor row alone therefore removes little from ⟨P, d⟩. The blob-free key's rank correlation with the LODO key is
  0.97–0.998.
- At cap 6 the blob-free key is +0.37 points better on the cancer mean, but mixed by fold (k562 +0.9, hek293t −0.15,
  kolf2 −0.7) and inside the fold-level noise.

**Decision:** keep the LODO key as defined (495, cap 6; the percentile rule of 543 for October). What makes the
held-out-K562 fold lose at high caps stays unexplained. Its two remaining donors (hct116, hek293t) may simply transfer
worst to K562.
