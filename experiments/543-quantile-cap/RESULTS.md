# 543 — a PERCENTILE cap for the LODO key: equal to cap 6 on A/B/C, +5 … +6 points more direction on October-like panels

2026-09-27 05:00Z, third session. Exploratory, after 540's preregistered read.
- `arms543.py`: PLOLO, lib420, the 4 A/B/C folds.
- `comp543o.py`: lib438, 438's October-like panels, the October stack h2c1.
- `run543.sh`: runner. Each run ~20 s, peak ~5.5 GiB, at ≥ 35 GiB free+inactive.
- q 1.75, realised 16, floor 0.5.

## Question
- The best fixed cap differs by surface. The A/B/C folds peak at 6 (540). On 438's October-like panels every larger
  cap wins up to uncapped (540o), and their key tails are much heavier (uncapped max w 19–36 vs ~10).
- A cap set at a **percentile of the panel's own carried key / median distribution** adapts to the tail. Does it keep
  cap 6's A/B/C gain and take the October-like gain?

## Known answers (all |d| 0)
- The per-donor mean equals lodo_key.
- Fixed caps 3 and 6 reproduce comp420's and 540b's A/B/C rows (4/4).
- 438's `h2c1|V1_5_3+sh1.5` (9/9).
- Fixed cap 6 reproduces 540o's October rows (9/9).

## Result (relative cos vs the cap-3 key; cancer folds positive)
| cap rule | A/B/C folds: cancer mean (k562 / hct116 / hek293t / kolf2) | cap it lands at on A/B/C folds | efp | efr | real | cap it lands at on October-like panels |
|---|---|---|---|---|---|---|
| fixed 6 (495) | +2.62 % (3/3): +0.82 / +3.84 / +3.20 / +8.34 | 6 | +12.7 % | +5.1 % | +8.2 % | 6 |
| p90 | +1.91 % (3/3) | 3.7–4.4 | +12.1 % | +5.5 % | +9.5 % | 4.8–9.2 |
| p95 | +2.60 % (3/3) | 4.7–6.1 | +18.1 % | +9.5 % | +13.3 % | 6.8–14.0 |
| **p96** | **+2.67 % (3/3)**: +0.93 / +3.90 / +3.19 / +6.16 | **5.0–6.7** | **+18.9 %** | **+10.6 %** | **+13.7 %** | **8.1–15.9** |
| p97.5 | +2.49 % (3/3): +0.47 / +3.98 / +3.02 / +7.49 | 5.6–7.2 | +20.8 % | +12.6 % | +14.7 % | 9.2–20.8 |
| p99 | +1.52 % (2/3): −2.20 / +4.11 / +2.65 / +10.3 | 7.4–9.3 | +21.1 % | +14.1 % | +15.5 % | 12.2–32.0 |

Every October-like fold is positive for every rule. The table shows panel means; per-fold values are in
`work/comp543o-*.json`.

- **p96 ≈ p97 is the sweet spot:**
  - On the A/B/C folds it equals cap 6 in value (lands at 5.0–6.7) and in direction (+2.67 vs +2.62 %).
  - On the October-like panels it adds +6.2 / +5.5 / +5.5 points over cap 6.
  - p99 starts to lose the A/B/C k562 fold, as the fixed caps above 8 did.
- **On the A/B/C production artifact, cap 6 is the 96.647th percentile** of carried key / median (10 of 294 carried
  targets exceed 6). With hi = percentile(key / median over carried, 96.647) the A/B/C weights equal cap 6's to
  4e-15.
  - Correction (coordinator, 2026-09-27, 497): that identity holds at the unrounded P = 96.64697522769474 (hi
    5.999999999999996). At the rounded 96.647, hi = 6.000055 and the 10 capped weights differ from cap 6's by 5.5e-5.
    `def_lodo_key.py --hi-pct` is run with the unrounded P; its KA-pct reproduces `dc374281` to 3e-14.
  - So the percentile rule and 495 are the same arm on A/B/C, and 495's board read calibrates the rule.

## Recommendation (the coordinator's call)
- **October (`DEF_DAY_LODO=1`):** make `def_lodo_key.py`'s cap the 96.647th percentile of the day-of panel's carried
  key / median instead of a fixed number.
  - On A/B/C this reproduces 495 exactly.
  - On heavier-tailed panels it caps higher: +5 … +6 points of direction over cap 6 on every October-like panel, 9/9
    folds.
- Adopt it only after 495 reads "transfers". If 495 fails, neither rule should go to October.
- Caveat: the October-like panels are 438's constructions. Their heavier tails could partly reflect how they were
  built. The percentile rule keeps the A/B/C behaviour fixed, so the October risk is confined to panels whose tails
  really are heavier.
