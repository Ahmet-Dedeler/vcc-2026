## 4. Setting (A) — leave-one-context-out (train on two cancer contexts + KOLF2, test on the third; `work/model269.log`)

Realised Δ`pds` at τ = 0.02, three cancer contexts (paired target bootstrap; captured = share of the labelled oracle +0.0647 / of 263's realised oracle +0.0459):

| policy | 3-ctx Δ`pds` [95 %] | P(≤ 0) | K562 / HCT116 / HEK293T | acted | captured (lab / 263) | +/3 |
|---|---|---:|---|---:|---|---:|
| keep (incumbent) | 0 | — | 0 / 0 / 0 | 0 | 0 | — |
| **oracle, labelled** (ceiling) | **+0.0647** [+0.0573, +0.0724] | 0 | +0.0647 / +0.0659 / +0.0637 | 0.58 | 100 % / 141 % | 3 |
| oracle, 263's realised arm | +0.0459 [+0.0382, +0.0538] | 0 | +0.0490 / +0.0468 / +0.0417 | 1.00 | 71 % / 100 % | 3 |
| always choose K562 | −0.0134 [−0.0220, −0.0048] | — | 0 / −0.0211 / −0.0191 | 0.59 | −21 % | 0 |
| always choose HCT116 | −0.0179 [−0.0248, −0.0109] | — | −0.0221 / 0 / −0.0315 | 0.61 | −28 % | 0 |
| always choose HEK293T | −0.0175 [−0.0257, −0.0092] | — | −0.0221 / −0.0303 / 0 | 0.62 | −27 % | 0 |
| always undeflate | −0.0123 [−0.0190, −0.0055] | — | −0.0100 / −0.0094 / −0.0174 | 1.00 | −19 % | 0 |
| **ridge** | **+0.0002** [−0.0027, +0.0034] | 0.467 | −0.0017 / −0.0005 / +0.0027 | 0.10 | +0.3 % / +0.4 % | 1 |
| ridge, choose only (no undeflate) | +0.0004 [−0.0024, +0.0036] | 0.393 | −0.0017 / +0.0000 / +0.0030 | 0.10 | +0.7 % / +1.0 % | 1 |
| **GBM** | **−0.0080** [−0.0149, −0.0013] | 0.991 | −0.0136 / +0.0009 / −0.0112 | 0.43 | −12 % / −17 % | 1 |
| GBM, choose only | −0.0047 [−0.0108, +0.0014] | 0.935 | −0.0099 / +0.0012 / −0.0055 | 0.33 | −7 % / −10 % | 1 |
| ridge, SHUFFLED FEATURES (3 seeds) | −0.0009 [range −0.0013 … −0.0005] | | | | | |
| GBM, SHUFFLED FEATURES | −0.0079 [−0.0099 … −0.0051] | | | | | |
| ridge, SHUFFLED TARGET LABELS | −0.0002 [−0.0013 … +0.0008] | | | | | |
| GBM, SHUFFLED TARGET LABELS | −0.0089 [−0.0131 … −0.0027] | | | | | |
| random action at ridge's / GBM's rate | −0.0014 / −0.0080 | | | | | |
| recipient-only (best fixed action on training) | 0 (never acts: every fixed action is negative on the training contexts) | | | | | |
| **target-only lookup, 066's same-target form** | **+0.0172** [+0.0093, +0.0253] | 0.000 | +0.0220 / +0.0213 / +0.0082 | 0.51 | **+27 % / +37 %** | **3** |
| target-only lookup at τ 0 | +0.0159 [+0.0076, +0.0243] | 0.000 | +0.0212 / +0.0198 / +0.0068 | 0.67 | +25 % / +35 % | 3 |

Out-of-sample prediction quality on the held rows (reference): ridge Spearman(pred, u) +0.075 on the 1,554 choose rows, AUC(u > 0) 0.49; GBM +0.115, AUC 0.54 / AUC(u > 0.05) 0.55; on the 853 undeflate rows Spearman +0.003 / −0.049 and AUC(u > 0) 0.44 / 0.42 (below chance: the deflation's per-target sign is line-specific, 262 §6, and nothing observable predicts it).

**The τ curve (3-context Δ`pds`, acted fraction; all actions).** Ridge: τ 0 −0.0046 (0.46) · 0.005 −0.0029 (0.34) · 0.01 −0.0026 (0.24) · **0.02 +0.0002 (0.10)** · 0.03 −0.0003 (0.03) · 0.05 −0.0004 (0.00) · ≥ 0.075 0 (never acts). GBM: τ 0 −0.0109 (0.69) · 0.005 −0.0091 (0.62) · 0.01 −0.0085 (0.56) · 0.02 −0.0080 (0.43) · 0.03 −0.0050 (0.29) · 0.05 −0.0008 (0.16) · 0.075 +0.0001 (0.06) · 0.10 +0.0002 (0.02) · 0.15 −0.0004 (0.01) · ≥ 0.2 0. Choose-only GBM peaks at +0.0016 (τ 0.075–0.10, acting on 2–4 % of targets, 1/3 folds). **No margin turns either model positive beyond its shuffled controls**: the ridge's +0.0002 at τ 0.02 is inside its shuffled-feature range's width of zero, and the GBM tracks its own shuffled-feature and shuffled-label controls (−0.008 … −0.009) along the whole curve — what the GBM "learns" is the label marginals (how often a switch loses), and acting on them costs what a random action costs at the same rate (−0.0080).

**The one positive row is the target-only lookup** — and §5 shows it is a property of the harness, not of the targets.

## 5. The target-only lookup: a PLOLO leakage path, measured and closed (`lookup269.py`, `mech269.py`, `xfree269.py`)

066 (HEK293T ↔ HCT116, 064's arms) found that a same-target lookup LOSES −0.007. Here the same lookup over the donor-choice actions WINS +0.017 at 3/3 with a bootstrap interval clear of zero and a shuffled-target permutation p = 0.002. Before reading that as target-level biology ("for target t, donor Z carries the conserved response"), the leakage path specific to leave-one-line-out has to be excluded: **the held line's own response sits as a donor row inside every training fold's donor mean.** For held X with donors (Y, Z), the training label of `choose Z` on fold Y (donors X, Z) is positive when X's row disagrees with Y's truth; the held label of `choose Z` on fold X is positive when Y's row disagrees with X's truth — both carry the pair (x, y)'s disagreement, which is a function of the held truth. A challenge line is a donor nowhere, so this path does not exist for A/B/C.

**T2 — where the lookup's signal comes from (τ 0.02, 3-context mean, shuffled-target permutation null in brackets):**

| training source for the lookup | all actions | choose only | per context (all actions) |
|---|---|---|---|
| all three other contexts (as in §4) | +0.0172 [+0.0092, +0.0253], null −0.0117 ± 0.0032, p 0.002 | +0.0205 [+0.0129, +0.0282] | +0.0220 / +0.0213 / +0.0082 |
| **the other two cancer folds only** (the held line is ONE OF TWO donors in every training mean) | **+0.0333** [+0.0257, +0.0411], null −0.0081, p 0.002 | **+0.0384** [+0.0314, +0.0456] — 59 % of the labelled oracle | +0.0346 / +0.0377 / +0.0277 |
| **KOLF2 only** (the held line is one of three donors, in a three-row mean) | **−0.0090** [−0.0172, −0.0009], null −0.0109, p 0.28 | −0.0083 [−0.0159, −0.0008] | −0.0043 / −0.0056 / −0.0171 |

Per (target, action), the Spearman between the training score and the held label is **+0.50 / +0.54 / +0.53** when the score comes from the other cancer fold (sign agreement 79–81 %) and **+0.03 / +0.06 / +0.04** when it comes from KOLF2 (sign agreement 60 %). The signal lives exactly where the held line's row is half of the training mean.

**T1 — the clean held context.** KOLF2's response is a donor on no other fold, so "KOLF2 held, trained on the three cancer folds" has no such path (an iPSC truth, base 0.640, labelled oracle +0.151):

| lookup on KOLF2 | Δ`pds` [95 %] | P(≤ 0) | acted | precision | captured | shuffled-target null, p |
|---|---|---:|---:|---:|---:|---|
| all actions, τ 0 / 0.01 / 0.02 / 0.05 | +0.0005 / +0.0037 / +0.0012 / −0.0000 (95 % ≈ ±0.018) | 0.47 / 0.34 / 0.45 / 0.50 | 0.67 / 0.54 / 0.48 / 0.37 | 0.46–0.53 | 0–2 % | −0.004 ± 0.008; p 0.19–0.34 |
| choose only, τ 0 / 0.01 / 0.02 / 0.05 | +0.0020 / +0.0028 / −0.0006 / −0.0019 | 0.42 / 0.38 / 0.53 / 0.59 | 0.59 / 0.48 / 0.44 / 0.34 | 0.47–0.52 | −1 … +2 % | p 0.22–0.43 |

Fixed actions on KOLF2 for scale: always choose K562 −0.0073, HCT116 +0.0000, HEK293T −0.0141, undeflate −0.0041. The ridge / GBM routers with KOLF2 held (T4, same features, trained on the three cancer contexts): ridge +0.0003 at τ 0 acting on 5 % (precision 0.21) and never acting above; GBM −0.0084 / −0.0033 / −0.0081 / −0.0041 at τ 0 / 0.02 / 0.05 / 0.10 (precision 0.33–0.44); Spearman(pred, label) +0.03 for both.

**Partial correlation (`mech269.py`).** Controlling the cancer-cancer transfer for the single scalar cos(x_truth, y_row) leaves the Spearman at +0.43 … +0.48 (from +0.50 … +0.57), so the leak is not that one number; the shared object is the held line's whole response inside the training fold's mean, which the next test removes entirely.

**T3 — a HELD-LINE-FREE target-level lookup (`xfree269.py`, 264's matched three-donor construction, `pertarget264m-*` / `pertarget264-kolf2`).** On every cancer fold 264 staged K562, HCT116, HEK293T and KOLF2 (at matched whitened amplitude) and scored the three "drop one donor" rows. For held X with donors (Y, Z, K): the candidate rows are mean(Y,Z) (drop K — the KOLF2-free incumbent-like row), mean(Z,K) (drop Y) and mean(Y,K) (drop Z), and each row's ABSOLUTE per-target `pds` exists on a fold whose truth is not X and whose row contains no X: mean(Z,K) against Y's truth (arm `drop:X` on fold Y), mean(Y,K) against Z's truth (`drop:X` on fold Z), mean(Y,Z) against KOLF2's truth (`drop:X` on the KOLF2 fold). KA: the fold totals reproduce 264 §5 to four decimals (`cd4svd5` 0.8154 / 0.8017 / 0.7810, drop-KOLF2 0.8179 / 0.8015 / 0.7776, whitened oracle drop 0.8623 / 0.8475 / 0.8166). Decision: the drop with the best held-line-free score, acting when it beats the keep row's score by τ; realised utility = 264m's per-target `pds` of that arm on fold X minus the drop-KOLF2 row's; the ceiling of this action set is +0.0526 (best of keep / drop Y / drop Z).

| held-line-free lookup | τ 0 | τ 0.02 | τ 0.05 | τ 0.1 |
|---|---|---|---|---|
| vs the keep score, raw | **−0.0152** [−0.0240, −0.0064] (acted 0.71, precision 0.30) | −0.0137 [−0.0225, −0.0051] | −0.0138 | −0.0110 [−0.0194, −0.0028] |
| vs the keep score, target-centred | −0.0152 | −0.0137 | −0.0138 | −0.0110 |
| best drop only, raw / centred | −0.0158 / −0.0180 | −0.0148 / −0.0163 | −0.0152 / −0.0129 | −0.0153 / −0.0095 |

Every cell is negative with a 95 % interval below zero (per context −0.02 … −0.03 on K562, ≈ 0 … −0.005 on HCT116, −0.01 … −0.02 on HEK293T); precision 0.30–0.36 — a wrong switch costs more than a right one earns (263 §5), and the held-line-free scores pick the right one no better than chance. The transfer statistic itself — Spearman between the held-line-free score difference s(drop Y) − s(drop Z) and the realised difference u(drop Y) − u(drop Z) on the held fold — is **+0.035 / +0.119 / +0.061** (permutation p 0.60 / 0.048 / 0.31; n 259 each), against **+0.21 / +0.20 / +0.29** for the leaky analogue in the same construction (the same actions scored by their utility on a fold where X is a donor) and +0.50 in 263's two-donor arrays. **Reading: target-level "which donor carries this target's conserved response" does not transfer between cancer lines once the held line's own response is removed from the training side. 066's negative stands; the +0.017 of §4 is the harness.** This closes the same-target lookup at its construction and it binds 267/268: any target-level learning across PLOLO folds must exclude the held line from the training folds' donor sets, or it measures the held truth.

## 6. Setting (B) — held target families (10 k-means families of the DepMap co-essentiality profiles; contexts seen)

| policy (τ 0.02) | 3-ctx Δ`pds` [95 %] | P(≤ 0) | K562 / HCT116 / HEK293T | acted | captured (lab / 263) | +/3 |
|---|---|---:|---|---:|---|---:|
| **ridge** | **−0.0005** [−0.0015, +0.0000] | 0.931 | −0.0016 / 0 / 0 | 0.01 | −0.8 % / −1.2 % | 0 |
| ridge, choose only | −0.0004 [−0.0012, +0.0000] | 1.000 | −0.0012 / 0 / 0 | 0.00 | −0.6 % | 0 |
| **GBM** | **−0.0056** [−0.0107, −0.0004] | 0.983 | −0.0046 / −0.0128 / +0.0008 | 0.30 | −8.6 % / −12 % | 1 |
| GBM, choose only | −0.0022 [−0.0066, +0.0021] | 0.843 | −0.0013 / −0.0089 / +0.0036 | 0.17 | −3.4 % | 1 |
| ridge / GBM, shuffled features | −0.0002 [−0.0005 … 0] / −0.0059 [−0.0072 … −0.0034] | | | | | |
| ridge / GBM, shuffled target labels | +0.0003 [+0.0001 … +0.0005] / −0.0055 [−0.0077 … −0.0035] | | | | | |
| random action at ridge's / GBM's rate | +0.0000 / −0.0045 | | | | | |
| recipient-only | 0 (never acts) | | | | | |

Quality: ridge Spearman +0.053 / AUC 0.46 (choose), −0.004 / 0.40 (undeflate); GBM +0.049 / 0.50 and −0.052 / 0.45. τ curves: ridge −0.0009 (0.12) at τ 0 → −0.0005 (0.01) at 0.02 → 0 from 0.05; GBM −0.0105 (0.61) → −0.0056 (0.30) → −0.0034 (0.06) at 0.05 → −0.0005 (0.02) at 0.075 → 0. With the target families held the models act on almost nothing (ridge) or lose at their shuffled controls' level (GBM); the target-only lookup is undefined here (no held target is in training).

