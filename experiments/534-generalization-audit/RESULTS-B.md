# 534 Part B — the sealed mock final (bulk channel): RESULTS

**On reserve R6, four fresh Arc-like panels that no decision has ever read (12 held-out cancer folds), the October bulk
recipe keeps 0.70 [0.62, 0.78] of its A/B/C gain in 438's instrument.** M0 → M6 is +0.0229 [+0.0204, +0.0255]
(levelled route, k = 1), against +0.0327 on A/B/C.
- **All six switches TRANSFER.** None is FRAGILE and none FAILS.
- **The survival is uneven.**
  - The three sharpening steps keep 0.53–0.61 of their A/B/C value.
  - H1 ×2 keeps 0.42.
  - The LODO key keeps only 0.30 [0.10, 0.55]. It is positive on 7 of 12 folds and negative on S4.
  - The percentile cap gains 3.3× its small A/B/C value.
- Both controls read ≤ 0.

Worker for 534 Part B (headless Opus), 2026-09-29, 13:53Z to 17:40Z. It ran `PREREGISTRATION-B.md`
(sha256 `74f108cf…`, unchanged) under `MANDATE-B.md`.
- **Panels priced: all four (S1, S2, S3, S4).** Plus two instrument panels: `efp` for the known answers, and A/B/C's own
  panel as a supplementary exact-design column.
- **Nothing run that the mandate forbids:** no `vcc`, no `cell-eval2`, no git, no builder edit, nothing in `submissions/`.
- **Nothing written outside `experiments/534-generalization-audit/`.** Every Part B process ran under a write guard
  (`guard534.py`) that raises on any write under the repo outside 534. A scan after the run found no modified file in any
  upstream experiment I read from.
- **One deviation, by the owner:** the memory-gate floor was lowered from 33 to 27 GiB (§7).

## 0. Answer

1. **Retention.** Pooled over R6, the cumulative October bulk gain is **+0.02285 [+0.02037, +0.02546]** (overall-
   equivalent, 438's levelled route, k = 1). The record's A/B/C value is +0.03274, so retention is **0.70 [0.62, 0.78]**.
   - Against the exact milestone designs priced on A/B/C's own panel (supplementary, §6) it is **0.68 [0.61, 0.76]**.
   - By panel: S1 +0.0305, S2 +0.0250, S3 +0.0221, S4 +0.0138.
   - Pricing M6 directly against M0 gives +0.0224, close to the sum of steps.
2. **Every milestone TRANSFERS** (§7.1: pooled Δ > 0 and the 95 % lower bound > 0). Retention by step:

   | step | retention |
   |---|---|
   | q 1.25 | 0.55 |
   | H1 ×2 | 0.42 |
   | q 1.5 | 0.53 |
   | q 1.75 | 0.61 |
   | LODO key | **0.30** |
   | percentile cap | **3.31** |

   - The two key switches carry most of the spread between panels (§3).
   - The LODO key reads −0.0006 on S4 (0/3 folds).
   - The percentile cap reads +0.0142 on S1 and +0.0007 on S4.
3. **Tuned-on check (§7.3): NOT OVER-READ.** Only 1 of the 3 October switches reads below half of its real/efp/efr mean.
   - The LODO key: R6 +0.00099 against a half-mean of 0.00147.
   - q 1.75: R6 +0.00282, above its tuned-on mean of +0.00251.
   - The percentile cap: R6 +0.00694, about equal to its tuned-on mean of +0.00703.
4. **Controls (§7.4) pass.**
   - C-a (CD4 ×1.5 on M4) is **−0.00087 [−0.00165, −0.00006]**, so 438's CD4 ×1 decision is confirmed, not flagged.
   - C-b (the LODO key permuted over targets) is **−0.0102 [−0.0116, −0.0085]**, negative on 0 of 12 folds.
5. **Board units (the §6 calibrated transfer).**
   - M0 → M6 projects to **+0.0182 [+0.0157, +0.0206]** on R6, against +0.0234 for the same designs on A/B/C (0.78).
   - Four steps have a direct A/B/C board read: 417a, 433k, 446b and 446bj2, together +0.0161. They project to +0.0078
     on R6, **about half** of the board-read September climb.
   - The percentile cap (495, board read 10-03) projects to +0.0090 on R6, and it is the switch that varies most by panel.
6. **Consequence (§7.5): no switch FAILS or is FRAGILE, so nothing goes to the coordinator as a switch-off candidate.**
   Information only:
   - The LODO key is the weakest October switch off A/B/C. It is the one the tuned-on panels over-read by the half-rule,
     and it is negative on S4.
   - Most of R6's October value sits in the percentile cap, as 522 found on the tuned-on panels.
7. **No surprise under §8's definition.** No step from M3 to M6 flips sign, and cumulative retention is ≥ 0.3.
   - Against §8's prior expectations:
     - the sharpening steps landed as expected (retention ≈ 0.5);
     - H1 ×2 is small, but its lower bound is above 0 (TRANSFERS, where FRAGILE was expected);
     - the LODO key missed its expected retention of ≈ 1 (0.30).
   - I re-ran the price known answer after seeing the result: |d| = 0 again. No instrument code changed between the known
     answers (15:45Z) and the R6 runs (from 15:51Z).

## 1. The step table (PREREGISTRATION-B §4, §6)

Levelled route, k = 1. Each panel's figure is the mean of its three cancer folds, with folds positive in brackets. "R6
pooled" is the mean of the 12 folds, with a paired target bootstrap (4,000 draws, seed 534). Retention is the R6 pooled
value divided by the record's A/B/C value; its interval comes from the bootstrap draws.

| step | adds | A/B/C (record) | real (record) | efp (record) | efr (record) | efp, exact design (this instrument) | S1 | S2 | S3 | S4 | **R6 pooled** [95 % CI] (n+) | verdict | retention [95 % CI] |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|
| M0→M1 | q 1.25 (417a) | +0.01178 | +0.00768 | +0.00626 | +0.00579 | +0.00626 | +0.00648 (3/3) | +0.00611 (3/3) | +0.00780 (3/3) | +0.00569 (3/3) | **+0.00652** [+0.00593, +0.00727] (12/12) | TRANSFERS | 0.55 [0.50, 0.62] |
| M1→M2 | H1 ×2 (417e's H1 half) | +0.00238 | +0.00511 | +0.00255 | +0.00031 | +0.00255 | +0.00220 (3/3) | +0.00094 (2/3) | +0.00063 (2/3) | +0.00019 (1/3) | **+0.00099** [+0.00037, +0.00150] (8/12) | TRANSFERS | 0.42 [0.15, 0.63] |
| M2→M3 | q 1.5 (433k) | +0.00856 | +0.00503 | +0.00391 | +0.00444 | +0.00391 | +0.00404 (3/3) | +0.00485 (3/3) | +0.00500 (3/3) | +0.00443 (3/3) | **+0.00458** [+0.00403, +0.00526] (12/12) | TRANSFERS | 0.53 [0.47, 0.61] |
| M3→M4 | q 1.75 (446b) | +0.00464 ᵃ | +0.00254 ᵃ | +0.00195 ᵃ | +0.00304 ᵃ | +0.00192 | +0.00162 (3/3) | +0.00334 (3/3) | +0.00295 (3/3) | +0.00337 (3/3) | **+0.00282** [+0.00232, +0.00342] (12/12) | TRANSFERS | 0.61 [0.50, 0.74] |
| M4→M5 | LODO key (0.5, 3) (446bj2) | +0.00329 ᵇ | +0.00343 ᵇ | +0.00349 ᵇ | +0.00187 ᵇ | +0.00317 | +0.00198 (3/3) | +0.00053 (2/3) | +0.00206 (2/3) | −0.00060 (0/3) | **+0.00099** [+0.00033, +0.00181] (7/12) | TRANSFERS | **0.30** [0.10, 0.55] |
| M5→M6 | percentile cap P 96.647 (495/522) | +0.0021 ᶜ | +0.0095 ᶜ | +0.0087 ᶜ | +0.0029 ᶜ | +0.00774 | +0.01422 (3/3) | +0.00928 (3/3) | +0.00360 (3/3) | +0.00066 (1/3) | **+0.00694** [+0.00517, +0.00828] (10/12) | TRANSFERS | 3.31 [2.46, 3.94] |
| C-a | CTRL: CD4 ×1.5 on M4 | +0.00170 ᵈ | −0.00160 ᵈ | −0.00111 ᵈ | −0.00236 ᵈ | −0.00087 | −0.00058 (1/3) | +0.00086 (3/3) | −0.00151 (1/3) | −0.00227 (0/3) | **−0.00087** [−0.00165, −0.00006] (5/12) | ≤ 0: pass | — |
| C-b | CTRL: LODO key permuted, on M5 | −0.02225 ᵉ | −0.02351 ᵉ | −0.01489 ᵉ | −0.00476 ᵉ | −0.01729 | −0.01932 (0/3) | −0.00780 (0/3) | −0.00978 (0/3) | −0.00390 (0/3) | **−0.01020** [−0.01160, −0.00850] (0/12) | ≤ 0: pass | — |
| **M0 → M6** | the six steps summed | **+0.03274** | +0.03328 | +0.02686 | +0.01835 | +0.02555 | +0.03054 | +0.02504 | +0.02205 | +0.01375 | **+0.02285** [+0.02037, +0.02546] | — | **0.70** [0.62, 0.78] |

Per-panel 95 % intervals (the same draws; `results534B.json`):

| step | S1 | S2 | S3 | S4 |
|---|---|---|---|---|
| M1→M2 | [+0.00105, +0.00318] | [−0.00052, +0.00186] | [−0.00036, +0.00173] | [−0.00115, +0.00156] |
| M4→M5 | [+0.00046, +0.00383] | [−0.00057, +0.00185] | [+0.00051, +0.00382] | [−0.00187, +0.00085] |
| M5→M6 | [+0.01034, +0.01821] | [+0.00496, +0.01211] | [+0.00161, +0.00551] | [−0.00239, +0.00315] |

The other steps are positive on every panel (lower bounds +0.00057 or more).

**Where the record columns come from** (quoted, not recomputed; `record534.py` → `record534.json`):
- **M0→M1, M1→M2, M2→M3:** the exact milestone designs, read from 438's stored rows through 438's own pricer. M1→M2 equals
  438 §3.1's printed H1 ×2 half. The M0→M1 A/B/C board read is +0.0070 (417a).
- **ᵃ M3→M4:** 438 §3.4, the stack's q 1.75 row minus its q 1.5 row. This is the CD4 ×1.5 analogue: 438's grid never
  had h2c1 at q 1.75.
- **ᵇ M4→M5:** 438 §3.5, "LODO V1 clip (0.5, 3), q 1.5 (vs K3 q 1.5)", on the stack at q 1.5.
- **ᶜ M5→M6:** 522 §7, k 1. That is 522's own unlevelled route, against cap 3, scaled by each fold's retention.
- **ᵈ C-a:** 438 §3.1, "CD4 ×1.5 on top of H1 ×2", at q 1.25.
- **ᵉ C-b:** 438's CTRL: the LODO key (0.3, 3) permuted, on PP at q 1.25, unlevelled.
- **The analogues are fair.** The exact designs priced on `efp` by this instrument (column 7) sit next to each record
  analogue:

  | step | efp exact | efp record analogue |
  |---|---:|---:|
  | M3→M4 | +0.00192 | +0.00195 |
  | M4→M5 | +0.00317 | +0.00349 |
  | M5→M6 | +0.00774 | +0.0087 (522's route) |
  | C-a | −0.00087 | −0.00111 |

  On A/B/C, the exact designs give M3→M4 +0.00502, M4→M5 +0.00344 and M5→M6 +0.00234 (§6).

## 2. Known answers on `efp` (PREREGISTRATION-B §5), run 15:45–15:47Z, before the draw

My instrument is `lib534.py`, a byte copy of `lib438.py` (sha256 `44aec6bd…`). It differs only in the write guard, the
docstring and HERE/WORK/ROOT, so every cache lives in `work/partB/`. 438's panel-independent `bankstats/` were cloned in
and are sha-identical. Scoring scripts: `comp534.py` (folds) and `full534.py` (full panel), both comp438/full438 restricted
to the preregistered designs. Pricer: `price534.py`, which is price438.Panel.

| known answer | compared against | result |
|---|---|---|
| KA-FOLD | my fold caches, built from scratch, vs 438's `efp/fold-<held>-held.npz` (targets, pre, truth, b_ctrl, gate, scored, wcols) | bit-identical, 3/3 folds |
| KA-438 | h1c1\|none, M0 `h1c1\|K3`, M1, M2, M3 and C-a `h2c1.5\|K3+sh1.75`, design and realised 16, plus each levelled twin, vs `comp438-<held>-held.json` (cos_agg, pds, pds λ 0.5/1, render A/B per frame, frac_capped, K_mean) | \|d\| = 0 on 23 rows × 3 folds |
| KA-493 (M4) | `h2c1\|K3+sh1.75` vs 493's `base q1.75` (lib438 path; design and realised 16) | \|d\| = 0, 3/3 |
| KA-546 (M5) | `h2c1\|V1_5_3+sh1.75` vs 546o's `base q1.75` (design and realised 16), and vs 522's `cap3` | \|d\| = 0, 3/3 |
| KA-522 (M6) | `h2c1\|V1pct+sh1.75` vs 522's `pctPROD` (cos_agg, pds, pds λ 0.5, K_mean, render A/B, the fold's hi, w_max) | \|d\| = 0, 3/3 |
| KA-COMP | full-panel PP constants vs 432's ef_b1 build meta | worst rel 8.93e-06 (tolerance 5e-5; 438's own value) |
| KA-KEY | full-panel 410 key vs the builder's artifact `c3e429dd` | w identical |
| KA-CLAIMS | ntclevel claims from the bundle's NTC cells vs 438's efp claims | \|d\| = 0 |
| KA-E438 | realised A/B, anchored AB / AB_level, ρ / ρ_level and level factors of 6 designs vs 438's `full438.json` / `level438.json` | worst rel 0 |
| KA-P | `price534.Panel` vs 438's `price438.Panel` on 438's jsons: 6 pairs × levelled / unlevelled × k 1 / 1.407 × 3 folds | 72 prices, \|d\| = 0; H1 ×2 half +0.002545904997, as in 438 |
| frames | bctrl D/E/F and the A/B/C frames vs 438's caches | sha256 identical |
| re-run after reading R6 | `ka_price534.py` | \|d\| = 0 |

- **M4 and M5 were checked against stored lib438 rows, not rows stored by 438 itself.** The preregistration says M0–M5
  reproduce "438's stored values", but 438's grid never had h2c1 at q 1.75. So M4 and M5 were checked against the lib438
  rows that 493 and 546o stored (both computed through lib438), and M5 also against 522's `cap3`. Their levelled twins have
  no stored reference anywhere; they come from the same code, whose other levelled rows reproduce.
- **The supplementary A/B/C panel passed its own known answers.**
  - Full panel: KA-COMP 8.87e-06 against 410b's meta, KA-KEY `d7f76a22` identical, KA-E438 0.
  - Folds: 438's abc rows and 493's abc M4 row reproduce with |d| = 0 on 25 rows × 3 folds, and the fold caches are
    bit-identical.
  - Its percentile caps land at 6.15 / 6.52 / 7.03 on the folds and 6.00 on the full panel, which are 522's A/B/C values.

## 3. R6: the draws and what the instrument sees on each panel

Drawn by `draw534.py` at 2026-09-29T15:50:51Z (`draws.json`), in the order S1, S2, S3, S4.
- **Exclusions at draw time:** 9 `pert_counts.csv` files (7 distinct panels, exactly the preregistered list) and 1,829
  excluded targets.
- **Essential list:** 2,021 genes, sha256 `b1e52795…`.
- **Pools before exclusion:** union 18,544; ∩ CD4 atlas 7,738.
- **Asserted for each panel:** 300 unique targets, 0 essential, 0 overlap with any excluded panel or earlier draw, all on
  the gene axis.
- **Bundles** follow `bundle238.py`'s recipe. Each relabelled control file is byte-identical to 199's and is stored as an
  APFS clone.
- **Operators:** 16 fits by `scripts/build_union_whitener.py` with 238's flags (`--basis-prune essential`, cap 12,000,
  seed 0, `--max-rss-gib 14`).
  - Every fold and full-panel operator passed `plolo.Operator._check_axis` on its own panel, against the builder's own
    `symbols:axis_sha256`.
  - `axes534.py` asserted that each operator's axis equals the builder's axis recomputed for that panel. It also asserted
    that the axis differs from A/B/C's (`edf6ba88…`, `6eaba848…`) and from all 16 tuned-on axes (9 files × 4 donor sets).
- **CD4:** `extract_cd4.py` per panel. The coverage gates passed on the efr-rule panels. On S2 and S4 they aborted, as on
  efp, and the relaxed run was used (DEF-RUNBOOK §5).

| panel | rule, seed | pool → essential-free | sha256 | K562 / HCT116 / HEK293T | H1-union rows | CD4 complete | percentile cap: full (folds) | PLOLO retention ret_f |
|---|---|---|---|---|---:|---:|---|---|
| S1 | efr, 534101 | 6,411 → 5,599 | `93237a3a…` | 300 / 300 / 300 | 31 | 294 | 13.07 (14.47 / 18.41 / 17.80) | 0.669 / 0.646 / 0.750 |
| S2 | efp, 534201 | 15,751 → 14,016 | `835d4b32…` | 137 / 291 / 293 | 13 | 172 | 11.88 (11.35 / 10.03 / 11.36) | 0.486 / 0.311 / 0.303 |
| S3 | efr, 534102 | 5,996 → 5,184 | `d9112014…` | 300 / 300 / 300 | 42 | 290 | 7.19 (11.40 / 9.99 / 8.20) | 0.507 / 0.472 / 0.505 |
| S4 | efp, 534202 | 15,151 → 13,416 | `88a7802a…` | 126 / 293 / 293 | 23 | 181 | 8.22 (8.03 / 9.94 / 12.70) | 0.296 / 0.208 / 0.211 |
| *efp (KA)* | | | | 146 / 293 / 296 | 28 | 182 | 12.43 (8.83 / 17.24 / 9.64) | 0.775 / 0.614 / 0.528 |
| *A/B/C (suppl.)* | | | | 269 / 290 / 294 | 77 | 276 | 6.00 (6.15 / 6.52 / 7.03) | 1 |

**Why the panels differ** (`extras534.json`: each step decomposed into relative aggregate cosine d_rel and
board-emulated Δpds):
- **The sharpening steps are paid mostly through pds.** Δpds_bd is +0.005 to +0.015 on every panel. Their direction
  gain shrinks with q: d_rel is +7 to +13 % at q 1.25, +3 to +9 % at q 1.5, and −1 to +5 % at q 1.75. They are
  panel-robust: 12/12 folds each.
- **The LODO key and the percentile cap are paid through direction (d_rel).** Their board value scales with the panel's
  PLOLO retention ret_f.
  - LODO key d_rel: S1 +4.6 %, S2 +2.5 %, S3 +5.3 %, **S4 −1.9 %**.
  - Percentile cap d_rel: S1 +33.0 %, S2 +32.4 %, S3 +9.3 %, **S4 +0.0 %**.
  - For comparison: A/B/C +3.4 % and efp +19.2 % (levelled). 522 measured +2.62 % and +19.43 % against cap 3 on its
    unlevelled rows.
- **S4 is the panel with the weakest donor transfer** (ret 0.21–0.30, the lowest measured on any panel).
  - The LODO key loses there on both direction (d_rel −1.9 %) and pds (Δpds_bd −0.0010).
  - The percentile cap adds no direction there, only +0.0011 of pds.
  - The panel's cumulative falls to +0.0138.

## 4. Controls, tuned-on check, surprise check (§7.3, §7.4, §8)

- **C-a, CD4 ×1.5 on M4: −0.00087 [−0.00165, −0.00006], positive on 5/12 folds.** Its interval sits below 0, so 438's CD4
  ×1 decision is not flagged.
  - S2 alone reads +0.00086 (3/3).
  - On A/B/C's own panel the same step is +0.00126 (supplementary), the A/B/C-only pattern 438 found.
- **C-b, the LODO key permuted over targets on M5: −0.0102 [−0.0116, −0.0085], 0/12 folds.** The key's per-target
  information is real on every panel.
- **Tuned-on check: 1 of the 3 in-scope switches is below half its real/efp/efr mean, so the panels are NOT OVER-READ**
  (the rule needs 2).

  | switch | R6 pooled | mean of real/efp/efr (record) | half | below half? |
  |---|---:|---:|---:|---|
  | M3→M4, q 1.75 | +0.00282 | +0.00251 | 0.00126 | no |
  | M4→M5, LODO key | +0.00099 | +0.00293 | 0.00147 | **yes** |
  | M5→M6, percentile cap | +0.00694 | +0.00703 | 0.00352 | no |

  With efp's exact-design value in place of its analogue, the LODO row is unchanged: the half-mean is 0.00141.
- **Surprise (§8):** no sign flip on M3–M6 pooled, and cumulative retention is 0.70 (≥ 0.3).

## 5. The calibrated board-unit projection (§6)

Each step is priced at its calibrated transfer, on the same draws:
- shape steps: k_sh 0.65 on X and 0.42 on pds;
- composition: k 1.3;
- keys: k 1.31.

The A/B/C column prices the exact designs on A/B/C's own panel through the same calibration (supplementary). Board reads
are quoted from FINDINGS §0.

| step | S1 | S2 | S3 | S4 | **R6 pooled** [95 % CI] | A/B/C, same instrument | A/B/C board read |
|---|---:|---:|---:|---:|---|---:|---:|
| M0→M1 | +0.00330 | +0.00290 | +0.00376 | +0.00280 | **+0.00319** [+0.00292, +0.00357] | +0.00617 | +0.0070 (417a) |
| M1→M2 | +0.00288 | +0.00133 | +0.00103 | +0.00027 | **+0.00138** [+0.00063, +0.00199] | +0.00333 | — (no read of the H1 half alone) |
| M2→M3 | +0.00190 | +0.00225 | +0.00222 | +0.00219 | **+0.00214** [+0.00188, +0.00249] | +0.00441 | +0.0046 (433k) |
| M3→M4 | +0.00048 | +0.00147 | +0.00104 | +0.00159 | **+0.00115** [+0.00091, +0.00145] | +0.00231 | +0.0023 (446b) |
| M4→M5 | +0.00264 | +0.00068 | +0.00258 | −0.00067 | **+0.00130** [+0.00051, +0.00230] | +0.00405 | +0.0022 (446bj2) |
| M5→M6 | +0.01902 | +0.01181 | +0.00458 | +0.00075 | **+0.00904** [+0.00679, +0.01069] | +0.00314 | pending (495, 10-03) |
| **M0 → M6** | +0.03022 | +0.02044 | +0.01521 | +0.00693 | **+0.01820** [+0.01566, +0.02061] | +0.02342 | +0.0161 (the four reads) |

- **Against A/B/C in the same instrument,** R6 keeps 0.78 of the calibrated cumulative.
- **Against the board.** The four board-read steps (417a, 433k, 446b, 446bj2) read +0.0161 on A/B/C and project to +0.0078
  on R6 (0.48).
  - The instrument reproduces the three sharpening reads on A/B/C (+0.0062 / +0.0044 / +0.0023 against +0.0070 /
    +0.0046 / +0.0023).
  - It over-reads the LODO key's read 1.8× (+0.0041 against +0.0022).
- **The percentile cap is the largest term on R6 and the least certain:**
  - about three times its A/B/C value in the same instrument;
  - ≈ 0 on S4;
  - its A/B/C board read is still pending.

## 6. Secondary readings (not the preregistered verdict route)

| step | primary (levelled k = 1) | unlevelled route | levelled, extrapolated to ρ* ᶠ | 522's route (M5→M6 only), k 1 | exact A/B/C (suppl.) | retention vs exact A/B/C |
|---|---:|---:|---:|---:|---:|---|
| M0→M1 | +0.00652 | +0.00724 | +0.00656 | | +0.01178 | 0.55 [0.50, 0.62] |
| M1→M2 | +0.00099 | +0.00081 | +0.00099 | | +0.00238 | 0.42 [0.15, 0.63] |
| M2→M3 | +0.00458 | +0.00605 | +0.00467 | | +0.00856 | 0.53 [0.47, 0.61] |
| M3→M4 | +0.00282 | +0.00377 | +0.00297 | | +0.00502 | 0.56 [0.46, 0.68] |
| M4→M5 | +0.00099 | +0.00136 | +0.00102 | | +0.00344 | 0.29 [0.10, 0.53] |
| M5→M6 | +0.00694 | +0.00500 | +0.00735 | +0.00728 (r +18.9 %) | +0.00234 | 2.97 [2.21, 3.55] |
| C-a | −0.00087 | −0.00066 | −0.00085 | | +0.00126 | — |
| C-b | −0.01020 | −0.00929 | −0.01033 | | −0.02691 | — |
| **M0 → M6** | **+0.02285** | | **+0.02357** | | **+0.03351** | **0.68 [0.61, 0.76]** |

- **ᶠ The levelling route falls short on some panels.** 438's levelled route applies one ntclevel step (438 §2.3). On the
  efr-rule panels the factors reach 1.28 (S1) and 1.39 (S3), and one step lands short of ρ* 0.979.
  - Landed ρ, M6 (`results534B.json` `landed_rho`): S1 0.938 / 0.930 / 0.916, S3 0.932 / 0.920 / 0.916.
  - Landed ρ, M4: S1 0.956 / 0.951 / 0.944, S3 0.950 / 0.943 / 0.938.
  - On the efp-rule panels: 0.958–0.983.
  - October's re-level is a converging secant loop (448). The ρ* column extrapolates every design's energy to ρ*
    (A, B × t², amplitude × t, PLOLO statistics held fixed).
  - It moves no verdict, and it adds +0.0007 to the cumulative (retention 0.72). The primary route slightly
    under-prices the q 1.75 and percentile-cap steps on low-energy panels.
- **Unlevelled route:** every step keeps its sign. The percentile cap reads lower (+0.0050) because it lowers the
  realised render energy: M6's ρ sits 3–11 % below M5's on every panel (`logs/full534-S*.log`). The unlevelled route
  charges that as distance from the cliff, and the re-level refunds it.
- **Exact A/B/C (supplementary, not preregistered):** the same instrument on A/B/C's own panel, with 438's abc inputs:
  - rules abc (CD4 all-rows sizing, SVD-5 deflation);
  - 187's claims, sd 0.002;
  - 062's held-out fold operators and the production operator `e2f65ecf…`.

  It replaces the record's analogues with exact designs. Retention moves from 0.70 to 0.68.

## 7. Deviations and notes

1. **The memory gate's floor, by the owner of 534.**
   - What happened:
     - PREREGISTRATION-B §9 requires `vm_stat` free+inactive ≥ 33 GiB before any heavy start.
     - The gate stayed closed from 14:16Z to 15:45Z, at 25.3–29.3 GiB.
     - memory_pressure read 78–82 % free, no builder or chain was running, and swap was unused.
   - Why: macOS keeps active ≈ inactive, so free+inactive ≈ (free + pageable)/2. That left ~26–30 GiB with the desktop
     apps open.
   - At 15:40Z I messaged the owner and added a read-only hook to the gate: `GATE-OVERRIDE.json`, which only the owner
     writes.
   - At 15:44Z the owner wrote a 27 GiB floor. It is kept by their `swapguard534.sh`, which returns the floor to 33 if
     swap free drops below 400 MiB. It is recorded here as their deviation.
   - Every other gate condition was unchanged: memory_pressure ≥ 40 %, no builder or chain, strictly serial, and the
     coordinator's 23:45Z and 00:30–06:00Z windows.
   - All 29 heavy starts opened at 27.1–35.9 GiB (`gate.log`). A 30 s sampler (`memwatch534.log`, 314 samples) read
     21.1–32.0 GiB, median 26.2.
2. **Operator fits at `--max-rss-gib 14`, authorised by the mandate** above the launcher's ~10 GB default. There were 16
   fits of 310–332 s each, with ru_maxrss 10.42–13.28 GiB.
3. **One scoring process per panel for the three folds**, so a panel's scoring needs one gate opening. The arithmetic is
   unchanged: the known answers ran this way.
4. **The supplementary A/B/C run's third fold** tripped 062's internal 8 GiB guard (`MemoryError`). Peak RSS is cumulative
   across folds in one process. The fold was re-run alone at 5.74 GiB, and all three abc folds' known answers pass. The
   R6 scoring processes peaked at 7.05–7.75 GiB.
5. **C-b's energy convention is 420/438's.** A permuted key keeps its unpermuted design's energies, so C-b uses M5's
   energies and level factors.
6. **Donor hygiene.**
   - Each fold's operator and keys are fitted on the other two lines only. The held-out line never enters them.
   - The panel's genes are off every operator axis (asserted).
   - h1_2025 appears only as one of the two H1-union residual sources on cancer folds, and no H1 surface is in scope.
   - The operator basis contains the source lines' rows for the panel's targets, as in 438 and as on the October day.
     Nothing from the held-out truth enters.
7. **Disk.**
   - After each panel was priced, its four operators (4.8 GiB) were deleted, so the coordinator's builds keep their
     headroom: 38 GiB free at the end.
   - The fits are deterministic (seed 0). Their sha256 values and exact commands are in `state534.json`, for C10b's
     re-fit.
   - The bundles cost ~0 disk (APFS clones).

## 8. Memory and timing (`heavy.log`, `state534.json`)

| phase | when (UTC) | heavy processes | wall | peak RSS (ru_maxrss) |
|---|---|---|---|---|
| gate closed at 33 GiB | 14:16–15:45 | — | — | — |
| efp known answers | 15:45–15:47 | 2 | 135 s | 5.27 / 7.56 GiB |
| draw + bundles | 15:50:51 | light | 5.5 s | — |
| S1 | 15:51–16:15 | 4 fits + 2 scoring | 24 min | fits ≤ 13.28, scoring ≤ 7.63 GiB |
| S2 | 16:15–16:38 | 4 + 2 | 23 min | ≤ 12.78 / 7.75 GiB |
| S3 | 16:38–17:02 | 4 + 2 | 24 min | ≤ 12.33 / 7.65 GiB |
| S4 | 17:02–17:26 | 4 + 2 | 24 min | ≤ 12.78 / 7.74 GiB |
| A/B/C supplementary | 17:26–17:29 | 3 (1 failed, re-run) | 2.4 min | ≤ 8.00 GiB |
| aggregate (4,000 × 2 bootstrap) | 17:30 | light | 20 s | 0.12 GiB |

In all: 29 heavy processes, 5,955 s of heavy wall clock, strictly serial.

## 9. Files (`work/partB/`)

- **Instrument:**
  - `guard534.py`: the write guard.
  - `lib534.py`: lib438's copy.
  - `designs534.py`: the preregistered designs and steps.
  - `comp534.py`: folds.
  - `full534.py`: full panel.
  - `price534.py`: the pricer and bootstrap.
  - `gate534.py`: the gate and heavy runner.
  - `drive534.py`: the ordered, resumable driver.
- **Steps:**
  - `precheck534.py`
  - `ka_price534.py`
  - `draw534.py`
  - `cd4_534.py`
  - `axes534.py`
  - `record534.py`
  - `aggregate534.py`
  - `census534.py`
  - `extras534.py`
- **Outputs:**
  - `results534B.json`: every number above plus the rule readings.
  - `tables534B.md`
  - `draws.json`
  - `record534.json`
  - `census534.json`
  - `extras534.json`
  - `state534.json`
  - `heavy.log`, `gate.log`, `memwatch534.log`, `logs/`
  - `<panel>/{comp534-*.json, arrays534-*.npz, full534.json, level534.json, keys534-full.npz, fold-*-held.npz, bundle/,
    cd4/, ops/*.json, axes534.json}` for efp, S1–S4 and abc.
- **Ledger draft for the owner:** `ledger-C10a.md`. It is not appended to `docs/CONFIRMATION-LEDGER.md`.
- **Reproduce:** `uv run python experiments/534-generalization-audit/work/partB/drive534.py all`. It is gated and skips
  done steps. Then run `aggregate534.py`, `extras534.py` and `census534.py`.
