## 7. Setting (C) — held context AND held family (the honest number; 30 folds; every cancer row predicted by a model that saw neither its line nor its family)

| policy (τ 0.02) | 3-ctx Δ`pds` [95 %] | P(≤ 0) | K562 / HCT116 / HEK293T | acted | captured (lab / 263) | +/3 |
|---|---|---:|---|---:|---|---:|
| keep | 0 | | | 0 | 0 | |
| oracle, labelled / 263 realised | +0.0647 / +0.0459 | | | 0.58 / 1.00 | 100 % / 71 % | 3 / 3 |
| best fixed action (always choose K562) | −0.0134 [−0.0220, −0.0048] | | | 0.59 | −21 % | 0 |
| **ridge** | **+0.0001** [−0.0024, +0.0028] | 0.497 | −0.0010 / −0.0004 / +0.0016 | 0.09 | **+0.1 % / +0.1 %** | 1 |
| ridge, choose only | −0.0003 [−0.0026, +0.0022] | 0.619 | −0.0025 / 0 / +0.0015 | 0.09 | −0.5 % | 1 |
| **GBM** | **−0.0091** [−0.0153, −0.0029] | 0.998 | −0.0165 / −0.0037 / −0.0072 | 0.44 | −14 % / −20 % | 0 |
| GBM, choose only | −0.0094 [−0.0153, −0.0034] | 0.999 | −0.0174 / −0.0025 / −0.0083 | 0.36 | −15 % | 0 |
| ridge / GBM, SHUFFLED FEATURES (3 seeds) | −0.0002 [−0.0005 … +0.0002] / −0.0086 [−0.0136 … −0.0040] | | | | | |
| ridge / GBM, SHUFFLED TARGET LABELS | −0.0006 [−0.0023 … +0.0011] / −0.0076 [−0.0100 … −0.0057] | | | | | |
| random action at ridge's / GBM's rate | −0.0021 [−0.0027 … −0.0016] / −0.0072 [−0.0088 … −0.0064] | | | | | |
| recipient-only | 0 (never acts) | | | | | |

Quality on the held rows: ridge Spearman +0.055 / AUC(u > 0) 0.48 (choose), +0.021 / 0.44 (undeflate); GBM +0.065 / 0.53 and −0.023 / 0.46. **τ curve** (all actions; Δ`pds`, acted): ridge τ 0 −0.0072 (0.44) · 0.005 −0.0034 (0.32) · 0.01 −0.0011 (0.23) · **0.02 +0.0001 (0.09)** · 0.03 −0.0011 (0.03) · 0.05 −0.0006 (0.00) · ≥ 0.075 0; GBM τ 0 −0.0117 (0.70) · 0.005 −0.0128 (0.64) · 0.01 −0.0124 (0.57) · 0.02 −0.0091 (0.44) · 0.03 −0.0063 (0.33) · 0.05 −0.0025 (0.16) · 0.075 +0.0003 (0.05) · **0.10 +0.0008 (0.02, 3/3)** · 0.15 +0.0003 (0.00) · ≥ 0.2 0. The GBM's best cell, +0.0008 at τ 0.10 on 2 % of the targets (about six targets per context), is 1.2 % of the oracle and inside the width of its own random-action range (0.0024); the ridge's +0.0001 sits between its shuffled-feature (−0.0002) and shuffled-label (−0.0006) controls, i.e. it abstains almost everywhere and the few actions it takes are a coin flip (precision at τ 0.02: 0.5 on 26 acted targets). The interval of the ridge policy is [−0.0024, +0.0028]: a real effect of the deflation's size (+0.0055 raw `pds`) is excluded; capturing 10 % of the oracle (+0.0065) is excluded.

Family sensitivity (average-linkage families, sizes 5–63): ridge **−0.0008** [−0.0029, +0.0015] (acted 0.09; −0.0045 / −0.0001 / +0.0024), GBM **−0.0106** [−0.0171, −0.0043] (acted 0.45; −0.0190 / −0.0038 / −0.0091) — the same reading with the unbalanced families (sizes 5–63).

