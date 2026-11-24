# 522 — October's bulk key: keep the percentile cap. The bounded strength key loses on every surface. A θ 0.25 strength rider adds +0.7 … +1.9 % on every October-like panel

2026-09-28, research worker (Claude Opus 5.5). Analysis only: no build, no `vcc`, no `cell-eval2`, no git, nothing written
outside this directory. The rider script `def_strength_key522.py` is **not applied** anywhere.

Surfaces:
- **A/B/C PLOLO folds**: 518's dumps, k562 / hct116 / hek293t, plus kolf2 alongside. This is 540's path: the A/B/C stack,
  q 1.75, realised KMAX 16.
- **438's October-like panels**: efp, efr, real, each on its three cancer folds. This is 546o's path: lib438, the October
  stack h2c1, the fold's own LODO fit, q 1.75, realised 16. The heavy part was run once per fold (`dump522.py`); every arm
  was then re-scored exactly from the dump.

**Units.**
- r is the relative aggregate cosine versus **cap 3**, the October default with `lodo_hi_pct` off.
- A panel's figure is the mean of its three cancer folds; "(n/3)" counts the folds where r > 0.
- "Points" means a difference of two such r values, in percent.
- Board units come from `price522.py`: overall ≈ 0.0710 × k × r × ret + 0.1245 × literal Δpds.
  - The first term is 518's identity (+1 % X = −0.0041 raw mse), with k = 1.31 for key steps and 1 as the low end.
  - ret is 438's per-fold retention: A/B/C 1, efp 0.53–0.78, efr 0.32–0.35, real 0.85–1.01.

## Answer

1. **Every known answer reproduces with |d| = 0** (§1): 546's percentile-cap row and 540o's cap-6 row on all 9
   October-like folds, 438's `h2c1|V1_5_3+sh1.5` on all 9, 518's strength-key A/B/C price, and 540's cap 6 on all 4 A/B/C
   folds.

2. **The percentile cap (P 96.64697522769474) is the robust October form.**
   - Versus cap 3: A/B/C **+2.62 % (3/3)**, identical to cap 6. efp **+19.43 % (3/3)**, efr **+11.58 % (3/3)**, real
     **+14.06 % (3/3)**.
   - Versus cap 6, with a target bootstrap:

     | panel | pct − cap 6 | 95 % interval |
     |---|---:|---|
     | A/B/C | +0.00 % | [−0.39, +0.45] |
     | efp | +5.95 % | [+0.38, +10.09] |
     | efr | +6.16 % | [+2.44, +9.63] |
     | real | +5.41 % | [+2.96, +7.67] |

   - It also holds on the Arc-like re-draws (§4b): essential-like targets removed, the rest re-drawn to A/B/C's strength
     mix. There it is ≥ cap 6 on 8 of 9 folds; the ninth is a tie at −0.07 points. By panel: efp +4.89, efr +3.40, real
     +1.24 points.

3. **Do not raise P or uncap for October.** Every October-like panel prefers the extreme setting: uncapped, P 98–99, or θ 1–2.
   - Leave-one-panel-out picks exactly those whenever the held-out panel is October-like (§5).
   - They lose on the one panel Arc actually built: A/B/C uncapped **+0.91 %** (k562 fold **−3.52 %**), P 99
     **+1.52 %** (k562 −2.20 %).
   - Minimax leave-one-panel-out over all four surfaces lands on P 96, cap 6 and θ 0.25, the A/B/C-calibrated forms.

4. **The bounded strength key does not replace the percentile cap: it loses to it everywhere.**
   - S θ 1 is 518's src fit on cap 3, with weights ≤ 4.9. It gives A/B/C +1.94 %, efp +8.91 %, efr +7.09 %, real +7.09 %.
   - Against the percentile cap (S1 − pct, bootstrap):

     | panel | S1 − pct | 95 % interval |
     |---|---:|---|
     | A/B/C | −0.64 % | [−3.18, +1.80] |
     | efp | −8.59 % | [−13.1, −0.6] |
     | efr | −3.68 % | [−8.35, +2.05] |
     | real | −6.10 % | [−8.56, −3.21] |

   - Its best θ on the October-like panels (θ 2) is worth +0.11 % on A/B/C.
   - Mechanism: a percentile-rank feature is bounded by construction, so it cannot express a heavy tail, and the
     October-like tails are real (point 5).
   - Its value is as a **tilt on top of** the LODO key, not as a substitute (point 6).

5. **Where the October tails come from (§4).**
   - **Two donors agreeing, not one.** In 236 of 245 tail targets (weight > 6 under the percentile cap) across the 9
     October-like folds, both LODO terms are ≥ 3 × the median. 7 are mixed, and only 2 come from a single donor.
   - **The panel mix makes the tail heavy.** Random draws are rich in weak targets:
     - 36–45 % of carried targets sit at the key floor (A/B/C 32–34 %);
     - 24–39 % of keys are ≤ 0 (A/B/C 14–17 %);
     - the src strength median is 1.06–1.19 (A/B/C 1.46);
     - efr's production-form median key is 0.013 against A/B/C's 0.032.

     A lower median inflates key/median for the strong targets.
   - **Reverse test on A/B/C.** Re-drawing A/B/C's own targets to the October-like mix raises its percentile cap from
     6.2–7.0 to 8.7. On A/B/C's own truth the percentile cap then beats cap 6 by **+1.39 points** (+0.02 / +2.19 / +1.95).
     Random A/B/C sub-panels of the same size keep hi at 6.3 and give −0.02.
   - **Essential-like targets.**
     - Where they exist (real: 10–12 in each fold's tail), raising them carries real's tail gain: +3.1 … +7.5 %.
     - With them removed, the percentile cap still beats cap 6 by +4.20 points on real, +7.80 on efp and +6.56 on efr.
     - efr's tail has **none**, and raising it still gains +0.9 … +11.2 %.
   - **The held-out truth backs the tail as a group, not target by target.**
     - Raw held-out truth energy barely follows the key: Spearman 0.0–0.2 on every surface, and tail energies are
       1.0–1.8 × the median target.
     - The truth's preferred weight along our direction (w × u) does follow it: Spearman 0.33–0.57.
     - The median ideal weight of a tail target is 3.7–9.0.
     - About half the tail is over-weighted at the percentile cap. The under-weighted half is worth more, so raising the
       tail from 6 to the percentile weights gains on 9/9 folds (+0.9 … +11.2 %).

6. **A small, consistent October rider:** the percentile cap × exp(0.25 z_src) on targets above the floor, never pushed
   below it (the "fpc" form). z_src is 518's src strength rank.
   - Versus the percentile cap (ratio): A/B/C **+0.05 %** (2/3), efp **+1.35 %** (3/3), efr **+1.81 %** (3/3), real
     **+0.70 %** (3/3).
   - Unchanged under 479's partners. Under S_P2: +0.07 / +1.36 / +1.88 / +0.69. Under P2g: +0.01 / +1.31 / +1.78 /
     +0.69. That is 9/9 October-like folds with either partner.
   - Target bootstrap of the plain form:

     | panel | rider − pct | 95 % interval |
     |---|---:|---|
     | A/B/C | +0.07 % | [−0.35, +0.48] |
     | efp | +1.36 % | [+0.71, +2.28] |
     | efr | +1.81 % | [+0.96, +2.63] |
     | real | +0.71 % | [+0.36, +1.07] |

   - Arc-like re-draws (plain form): +1.58 / +1.22 / +0.72 points, 8/9 folds.
   - Its worst fold is A/B/C hek293t at −0.31 %.
   - The board cannot read it before Oct 22: on A/B/C it is worth ≈ +0.0001 overall.

**Recommendation.**
- **October's key form stays the percentile cap**: `lodo_hi_pct = 96.64697522769474`, gated on 495 exactly as planned. Do not
  raise P, do not uncap, and do not swap in the bounded strength key.
- **Add the θ 0.25 floor-preserving strength rider as a separate October switch** (`lodo_strength_theta = 0.25`), applied after
  the same 495 read and a day-of rehearsal (§8).
  - Its tool and known answers are ready: `def_strength_key522.py`, and `--ka` passes.
  - If no rehearsal fits before Oct 22, ship the percentile cap alone. The rider is worth ≤ +0.001.
- If 495 fails, neither applies. 518's S1 contingency is unchanged (+1.94 % here on the A/B/C folds).

**Expected October value** (versus cap 3, board units at k 1.31 / k 1; the table is in §7):
- Percentile cap:
  - Arc-like re-draws: **+0.0031 … +0.0089** (efr … efp);
  - full October-like panels: +0.0037 … +0.0123;
  - A/B/C calibration: +0.0026 / +0.0021.
- Percentile cap + rider: Arc-like **+0.0035 … +0.0098**; full October-like panels +0.0043 … +0.0130.
- **Worst case across panels.** The worst panel mean is A/B/C random sub-panels: +2.32 % for the percentile cap, about
  +0.0023. The worst single fold is the Arc-like re-draw of real/hek293t: −0.32 % for the percentile cap and −0.36 % with
  the rider (cap 6 −0.62 %, cap 3 by definition 0).

## 1. Known answers (every script aborts on failure)

| row | expected | measured |
|---|---|---|
| strength522 on A/B/C: k562 / hct116 / hek293t / cd4_robust robust strength | 518's `strength518.npz` | bit-identical (4/4 arrays) |
| `dump522` per-donor LODO terms, mean | lodo_key | max relative \|d\| 4.0e-14 (9/9) |
| `dump522` 438's `h2c1\|V1_5_3+sh1.5` (design and realised) | `comp438-<fold>-held.json` | \|d\| 0, 9/9 |
| `dump522` 546o's `base q1.75`, `fixed 6`, `p96.6`; 540o's `V1_0.5_6` | comp546o / comp540o jsons (L.cmp: cos_agg, pds, pds λ 0.5 / 1, render A / B) | \|d\| 0, 9/9 × 4 rows |
| `fold522.Oct` on load: the same three 546o rows, re-scored from the dump | 546o | \|d\| 0 on every load |
| `fold522.Abc` on load (score518.Fold) | arms540's V1_0.5_3 / V1_0.5_6 | PASS, 4/4 |
| `arms522` A/B/C: zsrc | skey518.zfeat(F, 'src') | bit-identical, 4/4 |
| `arms522` A/B/C: S θ 1 on cap 3 (518's strength-key price) | skey518.json `src\|theta1\|on_cap3` r, dpds | \|d\| 0, 4/4 (+2.378 / +1.752 / +1.701 / +3.082 %) |
| `arms522` A/B/C: cap 6 | arms540 `V1_0.5_6` | \|d\| 0, 4/4 |
| `arms522` Oct: cap 6; pct 96.647 | 540o `V1_0.5_6`; 546o `p96.6` | \|d\| 0, 9/9 |
| `sp2_522` Oct: `p96.6 + S_P2`; A/B/C: `cap 6 + S_P2` | 546o; 540c | \|d\| 0, 9/9 and 3/3 |
| `arms522c` Oct: `p96.6 + P2g`; A/B/C: `cap 6 + P2g` | 546o; 540c | \|d\| 0, 9/9 and 3/3 |
| `boot522`: arrays → cos_agg | the exact arms (row_score identity) | < 1e-12, every pair |
| `def_strength_key522 --ka`: z on A/B/C; θ 0 on dc374281; θ 0.25 | abc518.zprod(src); dc374281's arrays; w_cap6 × exp(0.25 z) above the floor | PASS / PASS (arrays identical) / PASS |
| panel means reproduce the record | 540 (A/B/C cap 6 +2.62 %), 543 (A/B/C p96 / p97.5 / p99), 546 (pct: efp +19.43, efr +11.58, real +14.06 %), 540o (cap 6: +12.7 / +5.1 / +8.2 %) | equal |

## 2. Per panel: relative cos vs cap 3, and literal Δpds (full table: `work/analyse522.txt`)

| form | A/B/C folds: mean (n+): k562 / hct116 / hek293t; kolf2 | efp | efr | real | literal Δpds A/B/C / efp / efr / real |
|---|---|---|---|---|---|
| cap 3 (October default; 446bj2's key) | 0 | 0 | 0 | 0 | 0 |
| cap 4 | +1.74 % (3/3): +0.90 / +2.28 / +2.03; +3.79 | +5.64 % (3/3) | +1.69 % (3/3) | +3.51 % (3/3) | +0.0008 / +0.0004 / −0.0000 / +0.0006 |
| cap 6 (495) | +2.62 % (3/3): +0.82 / +3.84 / +3.20; +8.34 | +12.68 % (3/3) | +5.10 % (3/3) | +8.20 % (3/3) | +0.0014 / +0.0013 / +0.0006 / +0.0016 |
| cap 8 | +1.92 % (2/3): −1.00 / +4.09 / +2.66; +10.75 | +16.11 % (3/3) | +7.71 % (3/3) | +11.41 % (3/3) | +0.0017 / +0.0018 / +0.0012 / +0.0022 |
| uncapped | +0.91 % (2/3): **−3.52** / +4.12 / +2.12; +11.39 | +21.44 % (3/3) | +15.95 % (3/3) | +15.72 % (3/3) | +0.0018 / +0.0030 / +0.0025 / +0.0034 |
| pct P 90 | +1.91 % (3/3) | +12.08 % (3/3) | +5.46 % (3/3) | +9.46 % (3/3) | +0.0008 / +0.0012 / +0.0006 / +0.0020 |
| pct P 95 | +2.60 % (3/3) | +18.08 % (3/3) | +9.53 % (3/3) | +13.30 % (3/3) | +0.0014 / +0.0021 / +0.0013 / +0.0026 |
| pct P 96 | +2.67 % (3/3): +0.93 / +3.90 / +3.19; +6.16 | +18.92 % (3/3) | +10.61 % (3/3) | +13.67 % (3/3) | +0.0015 / +0.0022 / +0.0014 / +0.0027 |
| **pct P 96.64697522769474 (the candidate)** | **+2.62 % (3/3): +0.78 / +4.00 / +3.10; +6.67** | **+19.43 % (3/3)** | **+11.58 % (3/3)** | **+14.06 % (3/3)** | +0.0016 / +0.0023 / +0.0016 / +0.0027 |
| pct P 97.5 | +2.49 % (3/3): +0.47 / +3.98 / +3.02; +7.49 | +20.76 % (3/3) | +12.59 % (3/3) | +14.67 % (3/3) | +0.0016 / +0.0025 / +0.0018 / +0.0028 |
| pct P 99 | +1.52 % (2/3): **−2.20** / +4.11 / +2.65; +10.30 | +21.14 % (3/3) | +14.06 % (3/3) | +15.45 % (3/3) | +0.0017 / +0.0028 / +0.0023 / +0.0033 |
| S θ 0.5 (src rank, on cap 3; bounded) | +1.56 % (3/3) | +5.33 % (3/3) | +4.08 % (3/3) | +4.20 % (3/3) | +0.0002 / +0.0001 / −0.0001 / +0.0004 |
| **S θ 1 (518's src fit)** | +1.94 % (3/3): +2.38 / +1.75 / +1.70; +3.08 | +8.91 % (3/3) | +7.09 % (3/3) | +7.09 % (3/3) | +0.0004 / +0.0004 / −0.0002 / +0.0009 |
| S θ 2 | +0.11 % (2/3) | +11.96 % (3/3) | +9.95 % (3/3) | +9.69 % (3/3) | +0.0007 / +0.0011 / +0.0001 / +0.0019 |
| S + CD4 θ 1.5 | +3.24 % (3/3): +3.93 / +2.94 / +2.86; +4.63 | +11.67 % (3/3) | +9.27 % (3/3) | +6.89 % (3/3) | +0.0008 / +0.0010 / −0.0000 / +0.0013 |
| strength alone θ 1 (no LODO key) | −10.30 % (0/3) | −17.66 % (0/3) | −6.84 % (0/3) | −10.94 % (0/3) | −0.0025 / −0.0015 / +0.0022 / −0.0008 |
| cap 6 × S θ 0.25 | +2.78 % (3/3) | +14.85 % (3/3) | +7.34 % (3/3) | +9.68 % (3/3) | +0.0015 / +0.0014 / +0.0006 / +0.0018 |
| cap 6 × S θ 0.5 (518's S2 construction) | +2.64 % (3/3) | +16.53 % (3/3) | +9.25 % (3/3) | +10.84 % (3/3) | +0.0015 / +0.0016 / +0.0006 / +0.0020 |
| **pct × S θ 0.25 (the rider, plain)** | **+2.69 % (3/3): +1.22 / +4.07 / +2.79; +7.47** | **+21.02 % (3/3)** | **+13.53 % (3/3)** | **+14.87 % (3/3)** | +0.0016 / +0.0024 / +0.0017 / +0.0028 |
| pct × S θ 0.5 | +2.47 % (3/3) | +22.21 % (3/3) | +15.11 % (3/3) | +15.44 % (3/3) | +0.0016 / +0.0026 / +0.0017 / +0.0029 |
| pct × S θ 1 | +1.29 % (3/3) | +23.56 % (3/3) | +17.21 % (3/3) | +15.92 % (3/3) | +0.0017 / +0.0028 / +0.0019 / +0.0031 |
| gate: pct where z > 0, cap 6 elsewhere | +2.48 % (3/3) | +18.39 % (3/3) | +13.43 % (3/3) | +14.09 % (3/3) | +0.0016 / +0.0023 / +0.0017 / +0.0027 |
| CTRL: S θ 1 with z permuted | −2.19 % (0/3) | −2.12 % (0/3) | −3.30 % (0/3) | −1.47 % (0/3) | ≈ 0 |

- **The percentile cap's hi.** A/B/C folds 6.15 / 6.52 / 7.03 (kolf2 5.20); efp 8.83 / 17.24 / 9.64; efr 8.81 / 16.84 /
  10.27; real 9.33 / 16.03 / 10.25.
  - The production form uses three donor terms and has lighter tails. On efr's full panel (474's `def_lodo_key` artifact)
    the percentile cap lands at **7.48** (17 targets > 6, max 18.7); on A/B/C it is 6.00 (10 > 6).
  - The fold emulation of the rule caps the same fraction of targets (the top 3.35 %), so it is the right emulation, but
    the production cap on a heavy-tailed panel will be milder than the fold caps above.
- **Exploratory 522b.** Adding CD4's strength to the rider's rank is equivalent within noise: pct × S+CD4 θ 0.25 gives A/B/C
  +2.77 %, efp +20.65 %, efr +13.28 %, real +14.62 % (`arms522b.json`). src-only stays: it has fewer moving parts.

## 3. Bootstrap of the key comparisons (`boot522.py`; targets resampled within each fold, 2,000 draws; the panel mean of the fold ratios)

| comparison | A/B/C | efp | efr | real |
|---|---|---|---|---|
| pct vs cap 3 | +2.62 [+0.46, +5.01] | +19.43 [+9.36, +27.31] | +11.58 [+4.78, +17.88] | +14.06 [+9.80, +17.75] |
| pct vs cap 6 | +0.00 [−0.39, +0.45] | +5.95 [+0.38, +10.09] | +6.16 [+2.44, +9.63] | +5.41 [+2.96, +7.67] |
| S1 vs pct | −0.64 [−3.18, +1.80] | −8.59 [−13.14, −0.58] | −3.68 [−8.35, +2.05] | −6.10 [−8.56, −3.21] |
| pct × S0.25 vs pct | +0.07 [−0.35, +0.48] | +1.36 [+0.71, +2.28] | +1.81 [+0.96, +2.63] | +0.71 [+0.36, +1.07] |
| pct × S0.5 vs pct | −0.14 [−0.97, +0.66] | +2.38 [+1.13, +4.15] | +3.29 [+1.57, +4.94] | +1.21 [+0.53, +1.91] |
| cap 6 × S0.5 vs cap 6 (518's S2 construction on the src form) | +0.03 [−0.79, +0.84] | +3.39 [+2.05, +5.19] | +3.98 [+2.22, +5.49] | +2.43 [+1.47, +3.30] |
| uncapped vs pct | −1.70 [−3.06, +0.11] (P ≤ 0: 0.96) | +1.51 [−2.79, +5.09] | +3.95 [+0.29, +7.49] | +1.47 [−0.13, +3.94] |
| gate vs pct | −0.14 [−0.36, +0.02] | −0.77 [−4.15, +1.53] | +1.73 [−0.07, +3.79] | +0.03 [−0.10, +0.22] |

## 4. Where the October tails come from (4a)

### 4a. Tail anatomy (`tails522.py`)

- **Tail**: carried targets with key/median > 6 that the percentile cap weights above 6.
- **Terms**: the fold's two LODO terms v_j / median. "both" = min ≥ 3; "one donor" = min < 1.
- **Essential-like**: DepMap mean Chronos < −0.5 over 17 lines, or on the K562 essential list.
- **Low coverage**: a source's per-target noise is > 1.5 × its median (few cells).
- **u6**: the truth's local ideal-scale ratio on the cap-6 arm, (x/a)/(X/A). Above 1 means the truth wants more than 6.
- **Raises**: the exact arm "cap 6, with only this class raised to its percentile weight", as relative cos vs cap 6.

| fold | hi | tail | both / mixed / one-donor | ess-like | low cov | median u6 (share > 1) | ideal w median | raise all | raise non-ess | raise ess | raise high-cov |
|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|
| abc/k562 | 6.15 | 10 | n/a | 0.20 | 0.40 | 0.94 (0.50) | 5.6 | −0.04 % | −0.04 % | 0 | +0.06 % |
| abc/hct116 | 6.52 | 11 | n/a | 0.18 | 0.27 | 1.14 (0.73) | 6.8 | +0.15 % | +0.07 % | +0.08 % | −0.16 % |
| abc/hek293t | 7.03 | 15 | n/a | 0.13 | 0.20 | 1.06 (0.60) | 6.3 | −0.10 % | −0.02 % | −0.12 % | −0.12 % |
| efp/k562 | 8.83 | 9 | 9 / 0 / 0 | 0.22 | 0.44 | 1.27 (0.56) | 7.6 | +2.52 % | +3.04 % | −0.82 % | +0.75 % |
| efp/hct116 | 17.24 | 18 | 18 / 0 / 0 | 0.17 | 0.33 | 0.86 (0.44) | 5.2 | +10.43 % | +4.71 % | +4.08 % | −8.19 % |
| efp/hek293t | 9.64 | 11 | 11 / 0 / 0 | 0.27 | 0.45 | 1.21 (0.55) | 7.3 | +4.89 % | +1.93 % | +2.67 % | −2.01 % |
| efr/k562 | 8.81 | 27 | 27 / 0 / 0 | 0 | 0.33 | 0.62 (0.37) | 3.7 | +0.89 % | +0.89 % | 0 | +3.43 % |
| efr/hct116 | 16.84 | 54 | 50 / 4 / 0 | 0 | 0.24 | 0.99 (0.48) | 5.9 | +11.21 % | +11.21 % | 0 | +8.35 % |
| efr/hek293t | 10.27 | 25 | 22 / 2 / 1 | 0 | 0.24 | 1.49 (0.72) | 9.0 | +6.37 % | +6.37 % | 0 | +3.04 % |
| real/k562 | 9.33 | 30 | 30 / 0 / 0 | 0.33 | 0.23 | 1.11 (0.63) | 6.7 | +3.85 % | +0.38 % | +3.14 % | +1.53 % |
| real/hct116 | 16.03 | 44 | 43 / 0 / 1 | 0.27 | 0.23 | 1.25 (0.68) | 7.5 | +7.17 % | −4.80 % | +7.51 % | −0.84 % |
| real/hek293t | 10.25 | 27 | 26 / 1 / 0 | 0.37 | 0.19 | 1.24 (0.59) | 7.4 | +5.21 % | −2.16 % | +6.30 % | −0.22 % |

- **Genuine agreement.** 236 of 245 October-like tail targets have both donor terms ≥ 3 × the median, and only 2 are
  single-donor. The raw donor rows agree only loosely (cos 0.0–0.4): the key is |d| × alignment, and both lines see a
  large, aligned response.
- **Not a coverage artefact.** Low-coverage tail targets help on efp and hurt on efr/k562. There is no consistent sign.
- **Essential-like targets carry real's tail when they are present.**
  - Raising only the non-essential tail, with the essential tail held at 6, costs −4.8 / −2.2 % on two real folds.
  - That compares non-essential targets against essential targets whose true energy is large. The sub-panel test below
    removes the essential targets instead, which is the test that matches Arc's essential-free design.
- **Truth energy versus weight.**
  - Spearman(key/median, raw held-out energy) is 0.0–0.2 on every surface, A/B/C included. The tail's median truth energy
    is 1.0–1.8 × the median target's, against weights of 6–33.
  - Spearman(key/median, the truth's ideal weight along our direction) is 0.33–0.57 on October-like folds and 0.45–0.74
    on A/B/C. The key targets the aligned signal, and that supports tail weights of a median 3.7–9.
  - About half the tail is over-weighted at the percentile cap (u < 1). The under-weighted half is under-weighted by
    more, so the raise gains on 9/9 folds.
  - On A/B/C's k562 fold the top keys over-shoot the K562 truth. ACLY has key/median 11.3 against an ideal of 2.9;
    SMARCA5 9.4 against 4.8; UPF3B 9.2 against 4.4. The two Orion lines respond more strongly to these targets than K562
    does. That is why uncapped loses there.

### 4b. The panel mix (`subpanel522.py`, `sub522sum.py`; 20 draws per mode, everything re-derived on the sub-panel)

Matching variable: the mean robust strength over **the fold's own two sources** (never the held-out line), in quintile bins
of the reference panel's distribution.

| mode | panel | hi | cap 6 | pct | uncapped | S1 | pct × S0.25 | **pct − cap 6** | rider − pct | S1 − pct |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A/B/C re-drawn to the October-like mix | abc | 7.82 | +5.47 | +6.86 | +6.35 | +1.53 | +6.90 | **+1.39** (+0.02 / +2.19 / +1.95) | +0.04 | −5.33 |
| random A/B/C sub-panels, same size | abc | 6.26 | +2.34 | +2.32 | +1.05 | +1.91 | +2.51 | −0.02 | +0.19 | −0.40 |
| re-drawn to A/B/C's mix | efp / efr / real | 9.7 / 11.6 / 11.6 | +7.59 / +5.69 / +6.24 | +10.41 / +9.13 / +9.48 | +10.23 / +11.53 / +9.74 | +6.23 / +5.97 / +4.67 | +11.65 / +10.29 / +9.79 | **+2.82 / +3.44 / +3.24** (9/9 folds > 0) | +1.23 / +1.16 / +0.31 | −4.2 / −3.2 / −4.8 |
| essential-like removed | efp / efr / real | 12.0 / 12.1 / 12.8 | +13.70 / +5.12 / +7.09 | +21.50 / +11.69 / +11.29 | +28.04 / +16.06 / +16.36 | +9.45 / +7.36 / +6.57 | +23.51 / +13.66 / +12.34 | **+7.80 / +6.56 / +4.20** (9/9) | +2.01 / +1.97 / +1.05 | −12.1 / −4.3 / −4.7 |
| **Arc-like: essential-like removed and re-drawn to A/B/C's mix** | efp / efr / real | 9.8 / 10.9 / 9.5 | +9.64 / +5.76 / +3.92 | +14.53 / +9.16 / +5.16 | +16.48 / +11.59 / +9.82 | +7.89 / +6.10 / +4.00 | +16.11 / +10.38 / +5.87 | **+4.89 / +3.40 / +1.24** (8/9; real/hct116 −0.07) | +1.58 / +1.22 / +0.72 (8/9) | −6.6 / −3.1 / −1.2 |

Readings:
- **A heavy tail follows the panel mix.** Give A/B/C's own targets the October-like mix and its percentile cap rises to
  8.7 on two folds, and the higher cap is right on A/B/C truth.
- **Arc-like panels land between the two.**
  - Arc-like re-draws of the October-like panels still carry heavier tails than A/B/C: hi 9.5–10.9.
  - The percentile cap beats cap 6 on them too, by +1.2 … +4.9 points.
  - Uncapped is even better there, but it is A/B/C, the only real Arc panel, that punishes uncapped (§2). The adaptive
    rule is the one form that is right on both kinds of panel.

## 5. Leave-one-panel-out (4b; `analyse522.py`)

Each family's free parameter is chosen on the other panels, then scored on the held-out one.

| family | over the October-like panels (mean criterion): held-out efp / efr / real (choice) | over all four surfaces (minimax criterion): held-out A/B/C / efp / efr / real (choice) | mean, worst (4 surfaces, minimax) |
|---|---|---|---|
| fixed cap | +21.44 / +15.95 / +15.72 (uncapped) | +0.91 (uncapped) / +12.68 / +5.10 / +8.20 (cap 6) | +6.72, +0.91 |
| percentile cap | +21.14 / +12.96 / +15.45 (P 98–99) | +1.52 (P 99) / +18.92 / +10.61 / +13.67 (P 96) | +11.18, +1.52 |
| strength on cap 3 | +11.96 / +9.95 / +9.69 (θ 2) | +0.11 (θ 2) / +8.91 / +7.09 / +7.09 (θ 1) | +5.80, +0.11 |
| cap 6 × strength | +18.62 / +12.08 / +12.30 (θ 1) | +1.63 (θ 1) / +14.85 / +7.34 / +9.68 (θ 0.25) | +8.38, +1.63 |
| pct × strength | +23.56 / +17.21 / +15.92 (θ 1) | +1.29 (θ 1) / +21.02 / +13.53 / +14.87 (θ 0.25) | +12.68, +1.29 |
| (fixed) pct P 96.647 | +19.43 / +11.58 / +14.06 | +2.62 / +19.43 / +11.58 / +14.06 | +11.92, **+2.62** |
| (fixed) pct × S θ 0.25 | +21.02 / +13.53 / +14.87 | +2.69 / +21.02 / +13.53 / +14.87 | **+13.03, +2.69** |

- **The October-like panels always choose the extreme member.** Fitted on the other two, every family picks uncapped,
  P 98–99 or θ 1–2, and does well on the third.
- **A/B/C held out is the counter-example.** Those same choices then lose 1–2.5 points to the A/B/C-calibrated forms.
- **With A/B/C among the training panels, minimax picks the calibrated forms**: P 96, cap 6, θ 0.25.
- **The two fixed candidates have the best worst cases of any form**: P 96.647 (+2.62) and pct × S0.25 (+2.69).

## 6. The rider with 479's partners (`sp2_522.py`, `arms522c.py`; relative to the percentile cap)

The October recipe will carry S_P2 (`stercile`) or P2g (`tercile`); 479's pair reads on 10-02.
- **Plain rider versus P2g.** P2g puts every target at the key's MINIMUM ratio into LOW. A multiplier on the floor block
  would scatter it, and so would a target pushed below the floor. The builder's guard catches only the minimum ratio, and
  its tie rule would send the whole floor block (a third of an October panel) to MID.
- **Hence the "fpc" form.** It multiplies only targets above the floor and clips at 0.5. That touches 0–3 targets per fold.

| | A/B/C (per fold) | efp | efr | real |
|---|---|---|---|---|
| rider (plain) | +0.07 (+0.44 / +0.07 / −0.30) | +1.36 | +1.81 | +0.71 |
| **rider fpc** | +0.05 (+0.41 / +0.05 / −0.31) | +1.35 (3/3) | +1.81 (3/3) | +0.70 (3/3) |
| fpc under S_P2 | +0.07 | +1.36 (3/3) | +1.88 (3/3) | +0.69 (3/3) |
| fpc under P2g (comp546o's rule) | +0.01 | +1.31 (3/3) | +1.78 (3/3) | +0.69 (3/3) |
| (S_P2 on pct; P2g on pct) | +0.67; +0.94 | +0.71; +0.60 | +2.56; +0.83 | +0.24; +0.38 |

The rider is orthogonal to both partners: its increment changes by ≤ 0.1 points. On the October-like panels it is worth
about as much as either partner.

## 7. October value in board units (`price522.py`; versus cap 3; k 1.31 / k 1; retention-scaled)

| surface | cap 6 | **pct** | **pct + rider (plain)** | S1 | uncapped |
|---|---|---|---|---|---|
| A/B/C folds | +0.0026 / +0.0020 | **+0.0026 / +0.0021** | +0.0027 / +0.0021 | +0.0019 / +0.0014 | +0.0011 / +0.0009 |
| efp (full) | +0.0074 / +0.0057 | **+0.0114 / +0.0087** | +0.0123 / +0.0095 | +0.0051 / +0.0039 | +0.0123 / +0.0095 |
| efr (full) | +0.0016 / +0.0013 | **+0.0037 / +0.0029** | +0.0043 / +0.0033 | +0.0022 / +0.0017 | +0.0051 / +0.0040 |
| real (full) | +0.0072 / +0.0056 | **+0.0123 / +0.0095** | +0.0130 / +0.0100 | +0.0062 / +0.0047 | +0.0138 / +0.0107 |
| Arc-like efp | +0.0060 / +0.0046 | **+0.0089 / +0.0068** | +0.0098 / +0.0075 | +0.0046 / +0.0036 | +0.0100 / +0.0077 |
| Arc-like efr | +0.0020 / +0.0016 | **+0.0031 / +0.0024** | +0.0035 / +0.0027 | +0.0019 / +0.0015 | +0.0039 / +0.0031 |
| Arc-like real | +0.0036 / +0.0028 | **+0.0048 / +0.0037** | +0.0055 / +0.0043 | +0.0036 / +0.0028 | +0.0090 / +0.0070 |
| worst panel mean (r) | random A/B/C +2.34 % | random A/B/C **+2.32 %** | random A/B/C +2.51 % | A/B/C to October mix +1.53 % | A/B/C **+0.91 %** |
| worst single fold (r) | −0.62 % | **−0.32 %** | −0.36 % | −0.90 % | **−3.52 %** |

Caveats:
- The October conversion assumes that the A/B/C key-class transfer (k 1.31) holds, and that each fold's retention holds
  on the day.
- The Arc-like rows take the fold's retention unchanged.
- nmae (ballast K) is not modelled. 540 bounded it at ±0.0002 for cap 6.

## 8. The day-of change, if the rider is adopted (`def_strength_key522.py`; NOT applied)

- **What it does.** It post-processes the artifact `def_lodo_key.py` writes, with or without `--hi-pct`:
  w' = max(0.5, w × exp(0.25 z_src)) above the floor, and w' = w at the floor.
  - z_src: for each production donor (K562-GWPS, Orion HCT116, Orion HEK293T), the robust strength of its panel rows
    (518's statistic, with the bundle's panel as the excluded gene set), ranked over the targets that donor carries; then
    the mean over the donors carrying each target.
  - The input is only the donor banks and the bundle's `pert_counts.csv`, both on disk for any panel. It runs in about
    1 s at under 1.5 GB.
- **What stays the same.** It keeps every array of the input (targets, key, carried, n_donors) and replaces `w`. The
  builder reads only `targets` and `w`.
- **Known answers** (`--ka`, all PASS today):
  - KA-z: z on A/B/C == 518's `abc518.zprod(targets, 'src')`, bit for bit.
  - KA-inert: θ 0 on 495's dc374281 → every array identical, so the rider off reproduces 495's cap-6 artifact.
  - KA-arm: θ 0.25 on dc374281 → w == max(0.5, w_cap6 × exp(0.25 z)) above the floor, and min w' == 0.5.
  - On A/B/C the rider's multiplier spans 0.883 … 1.132 (median 1.000), and w' spans 0.5 … 6.79.
- **Why it cannot reproduce dc374281 at θ 0.25.** It changes A/B/C's weights by design. Its A/B/C effect is +0.05 % on
  the PLOLO folds (≈ +0.0001 overall), below any readable board step, so it is an October-only switch, decided from the
  October-like surfaces as 438's H1 × 2 / CD4 × 1 choice was.
- **Proposed wiring** (for the coordinator; `scripts/def_day.sh` is not edited here). Behind `DEF_DAY_LODO_STRENGTH_THETA`,
  and only with `DEF_DAY_LODO_HI_PCT` set: after `def_lodo_key.py` writes `$LK`, run
  `uv run python experiments/522-october-bulk-key/def_strength_key522.py --in "$LK" --bundle "$B" --theta 0.25 --out "$W/pert_scale_lodo_def_lo0.5_pct_s025.npz"`
  and pin that file as `lodo_lo0.5_pct_s025_def` instead of `$LK`.
- **Rehearsal design.**
  1. `--ka` (above).
  2. On a relabelled A/B/C bundle, θ 0 must give arrays identical to its `--hi-pct` input.
  3. On the efr bundle, the day-of rerun (`def_day.sh` with the pass-through on) must produce a build whose meta differs
     from the rider-off build only in `pert_scale` and its sha.
  - Tried here on 474's efr artifact and 497's relabelled-bundle artifact: the multiplier spans 0.899 … 1.133 on efr,
    with donor coverage 300/300/300 (`work/rider_*.json`).

## 9. Memory and deviations

- **Dumps.** The nine `dump522.py` runs peaked at 4.59–5.80 GiB each. They ran one at a time through `run_dumps522.sh`,
  whose gate required free+inactive ≥ 33 GiB on 2 consecutive checks 30 s apart. Each took 15–40 s, 14:43–16:00Z.
  This is the mandate's PLOLO exception.
  - The runner first used 3 checks 60 s apart and never cleared them while a build ramped. It was killed while waiting
    and relaunched with 2 × 30 s. Nothing heavy was running at the kill.
- **Light jobs.** Everything else peaked at 0.07–2.77 GiB: arms 0.4–1.1, sub-panels 1.3–1.6, tails 2.1, sp2 2.2,
  arms522b 2.8. That is under the 4 GB cap. Up to three light jobs ran concurrently beside one dump.
- **Nothing written outside this directory.** Every file newer than `work/.stamp_start` under the upstream experiment
  dirs was checked. The only one was `experiments/414-bulk-reprice/work/cliffcheck414.json`, which the production build
  chain rewrites; it was already modified in git status before this session. No lib438/plolo cache was created: the
  A/B/C surface points plolo at 420's existing read-only cache.
- **The rider changed shape during this experiment.**
  - Preregistered: the plain θ scan in `arms522.py`.
  - Exploratory, found afterwards: the floor-preserving "fpc" form (`arms522c.py`) and the CD4 variant (`arms522b.py`).
    The fpc form exists only for compatibility with P2g; it prices within 0.05 points of the plain form on every fold.

## 10. Files

- **Instruments.**
  - `dump522.py`: heavy, one per October-like fold, with its known answers.
  - `fold522.py`: exact arms from the dumps; sub-panels; the src feature.
  - `strength522.py`: robust strength for any panel, plus DepMap and essentiality. Its KA is 518's arrays.
- **Screens.**
  - `arms522.py`: every form on 13 surfaces.
  - `arms522b.py`: the CD4 rider, exploratory.
  - `arms522c.py`: the rider with P2g / S_P2, and the fpc form.
  - `sp2_522.py`: the rider under S_P2.
- **Analysis.** `analyse522.py` (tables, LOPO, pricing); `boot522.py` (target bootstrap); `tails522.py` (tail anatomy);
  `subpanel522.py` + `sub522sum.py` (the panel-mix test); `price522.py` (October value in board units).
- **Runner.** `run_dumps522.sh`: gated and serial.
- **Day-of candidate.** `def_strength_key522.py`: not applied.
- **Outputs** (`work/`).
  - Results: `analyse522.txt`, `analyse522.json`, `boot522.json`, `tails522.json`, `sub522sum.txt`,
    `subpanel522*.json`, `price522.txt`, `sp2_522.json`, `arms522{b,c}.json`.
  - Per-surface: `arms522-<panel>-<fold>.json` / `.npz`; `strength522.npz`.
  - Logs, including `run_dumps522.log`.
  - The dumps `dump522-*.npz` (70–120 MB each) are local only. Never add them to git: rebuild with `dump522.py` in ~20 s
    each at ≥ 33 GiB free.
