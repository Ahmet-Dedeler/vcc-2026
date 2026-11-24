# 706: where the bulk's remaining error sits by target strength: amplitude vs direction

Written 2026-10-15T05:29Z by the science session [bf003b]. Exploratory analysis on existing instruments (no board slot, nothing built).
Question (coordinator, 10-15): Illumina AI's profile (mse 0.79 scaled at a field-level pds) suggests "better amplitude on
strong targets". Does our per-target amplitude leave the strong tier under-weighted on held-out lines, and if not, what is
left?

## Instruments
- 522's exact fold surfaces (`fold522.load`, its own known answers asserted on every load): 9 October-like folds
  (efp / efr / real × held-out K562-GWPS / Orion HCT116 / Orion HEK293T, D/E/F composition rules) and 4 A/B/C PLOLO folds
  (held-out k562 / hct116 / hek293t / kolf2, A/B/C rules).
- Arm: production's key (percentile cap P 96.647), q 1.75 row sharpening, render, ballast KMAX 16.
- X = cos_agg = Σx / sqrt(Σa Σd) over targets, frame-averaged. 522's board units: overall ≈ 0.0710 · r · ret + 0.1245 · Δpds_lit
  (k 1).
- h1-test (Arc lab, Flex; 435's key rows; the h1_2025 donor never enters) as the lab-matched check.

## 1. Amplitude is near its predictable optimum (`budget706.py` → `work/budget706.json`)
| folds | X | 3-tier oracle by z_src | per-target oracle | k by z_src tier (lo / mid / hi) |
|---|---|---|---|---|
| A/B/C (4) | 0.030–0.066 | +0.2 … +1.5 % | +20 … +43 % | 0.64–0.85 / 0.91–1.09 / 0.99–1.11 |
| real (3) | 0.040–0.070 | +0.4 … +2.3 % | +22 … +26 % | 0.48–0.54 / 0.76–0.92 / 1.03–1.20 |
| efp (3) | 0.023–0.052 | +1.4 … +5.6 % | +25 … +40 % | 0.18–0.49 / 0.43–1.16 / 0.93–1.28 |
| efr (3) | 0.016–0.023 | +1.3 … +21 % | +50 … +84 % | 0.23–0.60 / 0.78–0.92 / 1.10–1.73 |

- The weakest third by source strength is over-weighted on all 13 folds (k < 1); the top third is slightly under-weighted
  on 10/13.
- A perfect 3-tier amplitude rule is worth only +0.2 … +2.3 % X on the A/B/C and real folds, about +0.0001 … +0.0016
  overall. 522's θ 0.25 rider (705) takes part of it.
- The per-target oracle (+20 … +84 %) needs per-target knowledge no public feature carries (522: truth energy follows the
  key at Spearman 0.0–0.2).
- **So the strong targets' remaining gap is direction, not amplitude.** Illumina's edge is not reproducible by
  re-weighting.

## 2. A direction lead: the raw donor mean on the strong tier (`dirswap706.py`, `h1swap706.py`)
- Arm: on the top z_src tertile, a norm-preserving blend of the stack composition toward the source donors' equal-weight
  unit mean U = Σ_j r_j/|r_j| (raw, unwhitened delta rows, panel columns zeroed), then the instrument's own key → sharpen
  → render. Target bootstrap 2,000.
- Controls:
  - fold truths ≈ the held-out line's donor rows (cos 0.99), so U is in the truth's space;
  - the shared across-target component is 0.2–1.8 % of energy in the truth, U and the stack (board-like, ~1 %), and
    centring U changes nothing;
  - removing 522's essential-like targets (DepMap < −0.5 or K562-essential) halves efp's gain and removes real/K562's,
    but efr's and real's Orion folds stand.

| surface | top tier, β 0.25: board units k 1 (CI-low) | folds > 0 | Δpds_lit |
|---|---|---|---|
| efp (3) | +0.0032 (−0.0001) | 3/3 | −0.0001 … −0.0141 |
| efr (3) | +0.0026 (+0.0011) | 3/3 | −0.0080 … +0.0008 |
| real (3) | +0.0043 (+0.0013) | 3/3 | −0.0093 … −0.0044 |
| **A/B/C cancer (3)** | **+0.0010 (−0.0002)** | 3/3 | +0.0012 / −0.0040 / −0.0068 |
| **A/B/C kolf2** | **−0.0019** | — | −0.0020 |
| h1-test (ew-cos, 85 targets) | top +8.5 % [−0.3, +22.6], all +12.8 % [+3.0, +27.7] | | mini-pds −0.009 / −0.026 |

- Larger β or more targets is worse on A/B/C, and kolf2 falls to −7.8 % X at β 0.5. The October-like panels favour more.
  This is the same pattern 522 recorded for caps: random draws prefer the extreme, while the one Arc-built panel does not.
- Mechanism (hypothesis): the strong-tier stack rows carry H1 ×2 and CD4 components that help on Arc-like truths. kolf2
  is stem-cell, like H1, and A/B/C is Arc lab. They do not help on public cancer truths, where the raw donor mean is the
  better direction.
- **Status: a lead, not October material.**
  - Minimax over the four panels gives +0.0003 (A/B/C, 4 folds), about the rider's size.
  - The A/B/C evidence (+1.6 … +2.6 % X on cancer folds, kolf2 negative) is the same strength 702 had before its HARM.
  - Tries: about 50 configurations (3 selections × 3 β × centre / essential-free × 13 folds + h1-test), so best-of-N bias
    applies.
  - **Rules or panel: already answered by 438 §4.4.** Every design prices within 0.00014 on `abcdef` (A/B/C's panel
    under the D/E/F rules) versus `abc`, so the October-like folds differ from A/B/C because the PANEL differs, not the
    day-of rules. D/E/F's panel is Arc-selected the same way A/B/C's is (438 §4.5's caveat: "adequate coverage in public
    datasets", high basal expression). The A/B/C folds are therefore the better proxy for D/E/F's value of this lever:
    +0.0003 … +0.0010.
  - It reopens only with an A/B/C board read of top-tier β 0.25 on ALLSW. The predicted +0.0003 … +0.0013 is readable,
    unlike the rider's +0.0001. A TRANSFERS read there would be on Arc's own panel and truth.

## Files
`budget706.py`, `dirswap706.py` (`--centre`, `--noness`, `abc`), `h1swap706.py`; outputs in `work/` (budget706.json,
dirswap706*.json, h1swap706.json, logs). Peak RSS ≤ 4.1 GB.
