# 545 — the stack's composition under the cap-6 key: the production mix (CD4 ×1.5, H1 ×2) stays optimal

2026-09-27 05:08Z, third session. Exploratory, `arms545.py` (PLOLO, lib420, 540b's path).
- Grid: CD4 multiplier m_cd4 {1, 1.25, 1.5, 1.75, 2, 2.5} × H1 union g {0.42, 0.56, 0.70, 0.84} (production 1.5 /
  0.56, i.e. `VCC_CD4_TARGET_F` 0.51 and `VCC_H1_UNION_G` 0.56).
- Each composition keyed by the fold's LODO key at (0.5, 6), q 1.75, realised 16. The DE carriers stay PP.
- Known answers: the production mix at cap 3 == comp420 `stack|V1_5_3+sh1.75` and at cap 6 == 540b `V1_0.5_6`, |d| 0
  on 4/4 folds. Peak RSS 5.4–5.6 GiB.

## Question
Composition changes transfer at the full ~1.3 rate (417e; 446bj2's key step). The mix was tuned under older keys
(420/438). Does the optimum move once the key's cap is 6?

## Result (relative cos vs the production mix at cap 6; cancer mean (folds +), then k562 / hct116 / hek293t / kolf2; literal Δpds cancer mean)
| composition | cancer rel cos | per fold | Δpds literal |
|---|---|---|---|
| **CD4 ×1.5, H1 0.56 (production)** | 0 | 0 / 0 / 0 / 0 | 0 |
| CD4 ×1.75, H1 0.70 (best) | +0.56 % (3/3) | +0.56 / +1.10 / +0.02 / +5.25 | −0.0021 |
| CD4 ×1.75, H1 0.84 | +0.47 % (2/3) | +0.44 / +1.49 / −0.53 / +10.9 | −0.0031 |
| CD4 ×2, H1 0.84 | +0.44 % (2/3) | | −0.0047 |
| CD4 ×1.5, H1 0.70 | +0.34 % (2/3) | | −0.0010 |
| CD4 ×1.75, H1 0.56 | +0.17 % (3/3) | +0.23 / +0.18 / +0.10 / −1.28 | −0.0011 |
| less CD4 (×1, ×1.25) or less H1 (0.42) | −0.4 … −2.9 % | | |
| production mix at cap 3 (446bj2) | −2.54 % (0/3) | | −0.0014 |

- Every direction gain from more CD4 or more H1 comes with a literal pds loss (−0.001 … −0.005 on PLOLO).
- The best, CD4 ×1.75 with H1 0.70: +0.56 % direction at the composition class is ≈ −0.003 raw mse, ≈ +0.0005
  overall. Its pds is −0.0021 × 0.34 ≈ −0.0007 raw, ≈ −0.0003 overall. Net ≈ +0.0003, with hek293t flat.
- That is not worth a build or a slot.

**Decision:** the production mix stays under the cap-6 key. The cap is the lever (540 / 495). The composition is
already at its optimum within ~0.0003.
