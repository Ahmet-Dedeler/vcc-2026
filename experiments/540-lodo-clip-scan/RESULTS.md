# 540 — the LODO key's clip: its cap of 3 is binding. Cap 6 prices +0.0019 … +0.0029 overall vs 446bj2 (candidate rule passes)

2026-09-27 04:40Z, third session. Preregistered in `PREREGISTRATION.md` (commit `9ebab956`, before any PLOLO number).
Scripts:
- `abc540.py` / `abc540b.py`: A/B/C energies, ballast, level maps.
- `arms540.py` / `arms540b.py`: PLOLO, lib420, the fold's own LODO fit.
- `analyse540.py` / `analyse540b.py`: pricer.
- `headroom540.py`: the side check below.
- Runners `run540.sh` / `run540b.sh`.

540b (caps 8, 10, 15, uncapped; floors 0.3 and 0.4 at cap 6) is an EXPLORATORY extension, written after 540's
preregistered read. Peak RSS 5.3–5.4 GiB per fold, after 04:30Z at ≥ 33 GiB free+inactive.

## Answer
- **The LODO key's cap of 3 costs direction.** Raising it to 4 or 6 improves the aggregate cosine on every fold.
- **Cap 6 (floor 0.5, `V1_0.5_6`) is the preregistered pick, and the exploratory extension confirms it as the peak.**
  - Relative cos +2.62 % (cancer mean, 3/3: k562 +0.82, hct116 +3.84, hek293t +3.20 %), kolf2 +8.34 %.
  - Board-calibrated vs 446bj2: **Δoverall +0.0025 … +0.0029 at the key class's k_X 1.31, +0.0019 … +0.0023 at k 1.**
- Raising the floor above 0.5, or lowering the cap below 3, only loses.
- **Preregistered candidate rule: all four items pass for cap 6.**
  1. Net at k_X 1.31 is +0.0025 with the worse ballast model (≥ +0.0008). Net at k 1 is +0.0019 (≥ +0.0004).
  2. Cancer folds positive 3/3; kolf2 +8.3 %.
  3. A/B/C energy ×1.013 / ×1.016 / ×1.005 (≤ 3 %); max per-row ratio ×2.04 (≤ 3).
  4. Cap 6's gain over its neighbours (cap 4 +0.0016, cap 8 +0.0018) is 0.0006–0.0008, less than its gain over the base.
     The curve is smooth, not a lone spike.

## Answers to the coordinator's two questions (04:42Z)
1. **Priced on the full LODO fit with hi = 6, i.e. weights(key, carried, 0.5, 6.0).**
   - PLOLO recomputes each fold's LODO key vector kv1 (lodo_bases + lodo_key, checked against 418's `V1__PP`), then
     calls `weights(kv1, carried>=2, 0.5, 6.0)`, then `key_scale` (the builder's energy-preserving constant, applied
     after the clip).
   - On A/B/C the arm uses the pinned artifact `8c82bef3`'s raw `key` vector (not its clipped `w`):
     `w = where(carried, clip(key / median(key[carried]), 0.5, 6), 1)`. That is the same function. The median
     normalisation is inside the clip, and the energy constant is outside, in apply_pert_scale.
   - Known answer: re-clipping that key at (0.5, 3) reproduces the artifact's w exactly, and apply_pert_scale on the
     dump's carriers reproduces 446bj2's bulk (|d| 0).
   - So `def_lodo_key.py --rules abc --hi 6` on data/challenge should reproduce this w exactly (its `--ka` already
     reproduces 8c82bef3's key to 1e-12). That is a good KA for 495's artifact.
2. **Cap 6 is the preregistered pick.** The preregistered grid had caps 2, 2.5, 4 and 6 at floor 0.5; cap 6 was best
   and passed the rule. The extension (8, 10, 15, uncapped) did not move the optimum. The cancer folds peak at 6, then
   fall to +1.9 % (8, 2/3), +1.2 % (10) and +0.9 % (uncapped, 2/3). Only kolf2 keeps rising (+11.4 % uncapped).

## Known answers (all pass)
| row | measured |
|---|---|
| clip (0.5, 3) == comp420 `stack\|V1_5_3+sh1.75` (designed and realised 16) | \|d\| 0, 4/4 folds (540 and 540b) |
| fold LODO key: weights(kv1, 0.3, 3) == 418's saved `V1__PP` | bit-identical, 4/4 |
| A/B/C: re-clip of `8c82bef3`'s key at (0.5, 3) == its w; apply_pert_scale(carriers, w) == the dump's bulk | exact, \|d\| 0 |
| A/B/C base emulated ballast K vs 446bj2's meta | median 7, max 16 (meta 7 / 16) |

## The clip scan (vs 446bj2's (0.5, 3); pds = 0.34 × literal PLOLO Δpds; mse re-levelled; ballast by the K and displacement models)
| clip (floor, cap) | rel cos cancer (folds +) | kolf2 | Δpds raw | Δmse raw k_X 1.31 (k 1) | Δnmae raw K / D | **Δoverall net, k_X 1.31 (k 1)** | A/B/C energy | max row |
|---|---|---|---:|---:|---:|---|---|---:|
| 0.4, 3 | −0.03 % (1/3) | +0.31 % | +0.0000 | +0.0002 | −0.0000 / +0.0002 | −0.0001 | ×0.998 | 1.01 |
| 0.6, 3 | −0.04 % (1/3) | −0.46 % | −0.0000 | +0.0002 | +0.0001 / −0.0002 | −0.0001 … +0.0000 | ×1.003 | 1.42 |
| 0.7, 3 | −0.19 % (1/3) | −1.09 % | −0.0000 | +0.0010 | +0.0002 / −0.0004 | −0.0002 … −0.0001 | ×1.006 | 1.91 |
| 0.85, 3 | −0.59 % (0/3) | −2.13 % | −0.0001 | +0.0032 | | −0.0007 … −0.0004 | ×1.012 | 2.73 |
| 1.0, 3 | −1.20 % (0/3) | −3.30 % | −0.0001 | +0.0064 | | −0.0013 … −0.0009 | ×1.018 | 3.65 |
| 0.5, 2 | −3.53 % (0/3) | −6.50 % | −0.0005 | +0.0188 | | −0.0036 … −0.0033 | ×1.008 | 1.65 |
| 0.5, 2.5 | −1.42 % (0/3) | −2.70 % | −0.0002 | +0.0076 | | −0.0015 … −0.0013 | ×1.004 | 1.24 |
| **0.5, 4** | **+1.74 % (3/3)** | **+3.79 %** | +0.0003 | −0.0093 (−0.0071) | −0.0004 / +0.0003 | **+0.0016 … +0.0019 (+0.0013 … +0.0015)** | ×0.995 | 1.30 |
| **0.5, 6** | **+2.62 % (3/3)** | **+8.34 %** | **+0.0005** | **−0.0140 (−0.0107)** | −0.0009 / +0.0006 | **+0.0025 … +0.0029 (+0.0019 … +0.0023)** | **×1.013 / 1.016 / 1.005** | **2.04** |
| 0.7, 2.5 | −1.76 % (0/3) | −4.04 % | −0.0002 | +0.0094 | | −0.0019 … −0.0015 | ×1.012 | 2.36 |
| 0.7, 4 | +1.71 % (3/3) | +2.99 % | +0.0003 | −0.0092 | −0.0003 / −0.0000 | +0.0017 … +0.0018 | ×1.000 | 1.40 |
| *540b:* 0.5, 8 | +1.92 % (2/3) | +10.75 % | +0.0006 | −0.0102 | −0.0012 / +0.0007 | +0.0018 … +0.0024 | ×1.030 / 1.036 / 1.020 | 2.93 |
| *540b:* 0.5, 10 | +1.19 % (2/3) | +11.39 % | +0.0006 | −0.0063 | | +0.0012 … +0.0017 | ×1.03 / 1.037 / 1.018 | 3.71 |
| *540b:* 0.5, 15 and uncapped (identical on these folds) | +0.91 % (2/3) | +11.39 % | +0.0006 | −0.0048 | | +0.0009 … +0.0014 | ×1.03 / 1.037 / 1.018 | 3.71 |
| *540b:* 0.4, 6 | +2.52 % (3/3) | +8.48 % | +0.0005 | −0.0134 | | +0.0023 … +0.0028 | ×1.012 | 2.05 |
| *540b:* 0.3, 6 | +2.41 % (3/3) | +8.58 % | +0.0005 | −0.0129 | | +0.0022 … +0.0027 | ×1.011 | 2.05 |
| *540b:* 0.3, uncapped | +0.69 % (2/3) | +11.55 % | +0.0006 | −0.0036 | | +0.0006 … +0.0012 | ×1.028 | 3.72 |

- Per fold at cap 6: literal Δpds +0.0012 / +0.0011 / +0.0020 (cancer), kolf2 −0.0004.
- What the cap does on A/B/C: the energy constant falls from 0.577 to 0.415, so all rows shrink by 28 %, and the
  top-key 16 % grow up to ×1.44 in effect (×2.04 in realised row energy). Mean ballast K falls from 8.83 to 7.54.
- The two ballast models disagree in sign for cap changes. The K model gives −0.0009 raw nmae (fewer ballast cells),
  the displacement model +0.0006. The range above carries both, ±0.0002 overall.
- The board-emulated pds reads negative for every cap raise (cap 6: −0.0010 on PLOLO). For key changes that emulation
  mis-predicted 446bj2's step 7× (492g). If it were right here, cap 6 would lose 0.0003 overall. It would still price
  ≥ +0.0016.

## Why cap 3 was never revisited
- 415 found an uncapped LODO key stronger on PLOLO (+7 … +10 % vs +3.3 % at cap 3) but kept 3 for caution: "uncapped
  weights reach 10×, which stresses ballast K and the jackknife credit".
- 410c (the uncapped K3 key) was built but never submitted. 414 priced it flat for two reasons:
  - context A's cliff, which the ρ-matched re-level now handles;
  - the board-emulated pds cost, which 446bj2's step showed is wrong for key changes (literal pds × 0.34 predicted that
    step, +0.00024 vs +0.00021).
- Here both are handled: every arm is re-levelled, and pds is literal. The uncapped key is indeed worse than cap 6 on
  the cancer folds. The optimum is interior.

## Side check: how much key headroom is left (`headroom540.py`, 420's saved per-target arrays, no new truth)
- An oracle per-target amplitude (each target's regression coefficient of truth on prediction) would raise the
  aggregate cosine +22 / +24 / +26 % over the LODO key on the cancer folds (+52 % kolf2). The LODO key captured
  ~3–5 points of the ~27 % headroom over K3.
- A first-order analytic rescaling of the saved arrays reproduces the exact V1 → V1_5_3 change with the right sign and
  ~70 % of its size.
- On that analytic route, monotone transforms of the (0.3, 3)-clipped LODO weights buy at most +0.3 %, with mixed fold
  signs: powers 0.5–2, floors 0.3–0.7, and a strength tilt (strength / median)^β with β −0.5 … +0.25.
- Taken with the scan above, the reachable gain is in the key's upper tail (the cap). Reshaping its lower and middle
  range gains nothing. Getting more of the oracle headroom would take new per-target information, not a new transform.

## Hand-off
- The coordinator builds cap 6 today as **495**: 446bj2 + the LODO key at (0.5, 6), a new stack dump, re-levelled, one
  variable vs 446bj2.
- Predicted board row vs 446bj2:
  - Δpds raw +0.0005 (literal; the β emulation says −0.0004);
  - Δmse raw −0.0140 at k_X 1.31 (−0.0107 at k 1);
  - Δnmae raw −0.0009 … +0.0006 (ballast);
  - fid / reach / jac pinned;
  - **Δoverall +0.0019 … +0.0029.**
- Suggested bins:
  - transfers: overall ≥ +0.0012 and Δmse ≤ −0.006;
  - partial: +0.0004 … +0.0012;
  - fails: < +0.0004 or Δmse > −0.002;
  - DE leak check: |Δnmae| > 0.002 raw.
- October (`DEF_DAY_LODO=1`, `def_lodo_key.py --hi`): cap 6 must be re-checked on 438's October-like panels before it
  is adopted there. 540o does that next.
- 479's per-target q (P2g / S_P2) and the cap look orthogonal: every capped target is in the HIGH key tercile at
  q 1.75, so raising the cap changes only HIGH rows' amplitudes. 540c measures the interaction on PLOLO.

## Addendum 540o / 540c (04:50Z, EXPLORATORY): the cap on October-like panels, and whether 479 and 495 stack

### 540o — the cap under the October key on 438's panels (`comp540o.py`, lib438, h2c1, the fold's own LODO fit at (0.5, hi), q 1.75)
Known answer: the LODO key at (0.5, 3), q 1.5 == 438's `h2c1|V1_5_3+sh1.5`, |d| 0 on 9/9 panel/folds. Direction vs
the October default (0.5, 3), mean of three cancer folds (folds positive):

| panel | cap 4 | **cap 6** | cap 8 | uncapped | cap 3 + P2g | cap 3 + S_P2 | **cap 6 + P2g** | **cap 6 + S_P2** |
|---|---|---|---|---|---|---|---|---|
| efp | +5.6 % (3/3) | **+12.7 % (3/3)** | +16.1 % (3/3) | +21.4 % (3/3) | +0.84 % | +1.32 % | +13.4 % | +13.9 % |
| efr | +1.7 % (3/3) | **+5.1 % (3/3)** | +7.7 % (3/3) | +16.0 % (3/3) | +1.18 % | +2.36 % | +6.0 % | +8.0 % |
| real | +3.5 % (3/3) | **+8.2 % (3/3)** | +11.4 % (3/3) | +15.7 % (3/3) | +0.91 % | +0.89 % | +8.8 % | +8.7 % |

- Literal Δpds at cap 6: +0.0013 / +0.0006 / +0.0016. The board-emulated pds reads negative on efr and real (the
  emulation 492g showed wrong for key changes).
- **On October-like panels the cap is worth more than on the A/B/C folds, and every cap from 4 to uncapped helps on
  9/9 folds.** There the gain keeps rising through uncapped, whereas the A/B/C folds peak at 6. Part of the size is the
  lower base cosine on these panels (relative changes are larger), but the sign and the ordering are unambiguous.
- **For `DEF_DAY_LODO=1`:** raise `def_lodo_key.py --hi` from 3. Cap 6 is supported on both surfaces. Anything above 6
  is supported only on the October-like panels (the A/B/C folds lose from 8 upward), so 6 is the robust choice until
  495's board read calibrates the class on A/B/C.
- The per-target arms stack with the cap on these panels too. The combined gain ≈ the product of the parts:
  efp 13.4 vs 13.6 %, efr 6.0 vs 6.3 %, real 8.8 vs 9.2 % for P2g; for S_P2 14.0 vs 14.2 %, 8.0 vs 7.6 %, 8.7 vs 9.2 %.

### 540c — 479 × 495 on the A/B/C folds (`arms540c.py`, `abc540c.py`, `analyse540c.py`)
Known answers:
- cap 3 q 1.75 == comp420; cap 3 + P2g and cap 3 + S_P2 == 492g's rows; cap 6 == 540b's row. All |d| 0 on 4/4 folds.
- The two-step pricer reproduces 492g's gross for P2g and S_P2 on cap 3 (|d| 2e-19 and 2e-17).
- P2g's terciles at cap 6 equal cap 3's exactly (100 / 100 / 100): every capped target is HIGH either way.

| step | rel cos (folds +) | kolf2 | Δoverall net (k_X 0.61 / 0.12) | energy vs 446bj2 |
|---|---|---|---|---|
| P2g on cap 3 (= 479's P2g) | +1.26 % (3/3) | +0.68 % | +0.0012 / +0.0008 | ×0.989 |
| S_P2 on cap 3 (= 479's S_P2) | +1.43 % (3/3) | +1.72 % | +0.0014 / +0.0008 | ×0.94 |
| **P2g on cap 6** (increment over 495) | +0.96 % (3/3) | +0.30 % | **+0.0012 / +0.0009** | ×1.010 (cap 6 + P2g) |
| S_P2 on cap 6 (increment over 495) | +0.75 % (3/3) | +0.99 % | +0.0010 / +0.0008 | ×0.97 (cap 6 + S_P2) |

- **Totals vs 446bj2:** the cap step at the key class plus the increment at the sharpening class (k_X 0.61), with the
  worse ballast model on the cap step.
  - cap 6 + P2g: **+0.0037** (key at k_X 1.31) / **+0.0031** (key at k 1).
  - cap 6 + S_P2: +0.0035 / +0.0029.
- **They stack.**
  - P2g's increment is unchanged on cap 6: less direction, but less ballast growth (low-tercile K +0.93 vs +1.30).
  - S_P2's increment shrinks by a quarter: the cap already reallocates amplitude toward strong, consistent rows,
    which overlaps with strength.
  - With cap 6, P2g becomes the better per-target partner on A/B/C. On October-like panels the two partners are within
    each other's noise.
- Sequencing (the coordinator's call): 479 (P2g vs S_P2 on 446bj2) and 495 (cap 6 on 446bj2) are each one variable vs
  446bj2. The combination is a natural third arm once both have read.

### Per-fold detail for the October cap choice (05:00Z)
| surface / held-out fold | cap 4 | cap 6 | cap 8 | uncapped | uncapped max w |
|---|---|---|---|---|---|
| A/B/C: k562 | +0.90 % | +0.82 % | **−1.00 %** | **−3.52 %** | 11.3 |
| A/B/C: hct116 | +2.28 % | +3.84 % | +4.09 % | +4.12 % | 10.3 |
| A/B/C: hek293t | +2.03 % | +3.20 % | +2.66 % | +2.12 % | 9.6 |
| A/B/C: kolf2 | +3.79 % | +8.34 % | +10.75 % | +11.39 % | 9.8 |
| efp: k562 / hct116 / hek293t | +3.9 / +5.7 / +7.3 % | +8.4 / +13.4 / +16.2 % | +10.9 / +17.2 / +20.3 % | +7.5 / +29.9 / +26.9 % | 21.5 / 36.3 / 27.4 |
| efr: k562 / hct116 / hek293t | +1.7 / +1.1 / +2.3 % | +3.6 / +3.6 / +8.1 % | +4.5 / +6.4 / +12.2 % | +9.9 / +18.4 / +19.5 % | 22.3 / 36.5 / 18.8 |
| real: k562 / hct116 / hek293t | +3.6 / +3.3 / +3.7 % | +8.0 / +8.6 / +8.0 % | +10.9 / +12.3 / +11.0 % | +13.8 / +16.4 / +17.0 % | 23.4 / 32.9 / 18.7 |

- Cap 8 beats cap 6 on all 9 October-like folds. Uncapped beats cap 8 on 8 of 9, the exception being efp with k562
  held out.
- The October-like panels have much heavier key tails: the uncapped max weight is 19–36 there, ~10 on the A/B/C folds.
- The only fold where caps above 6 lose is A/B/C's held-out K562.
  - Hypothesis (untested): there the key's two terms both come from Orion lines (hct116, hek293t; same lab and
    technology), so their shared extremes transfer worst.
  - Production has three terms per target, K562 plus both Orion lines, which resembles the other folds.
- October recommendation: cap 6 now. Cap 8 is the candidate to revisit once 495's board read calibrates the class.
  It wins on every October-like fold and costs only the k562 A/B/C fold.
