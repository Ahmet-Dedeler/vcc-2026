# 438: the September bulk stack re-priced on October-like panels. CD4 ×1.5 does not survive the panel change, H1 ×2 does, q 1.5 and the LODO key keep their sign; one command prices it on Oct 22

Status: **COMPLETE** (2026-09-23, 16:30Z–22:10Z).
- 420's PLOLO pricer, re-wrapped so that the panel is an argument.
- Five panels:
  - `abc`, the real A/B/C panel;
  - `abcdef`, A/B/C's panel under the D/E/F composition rules;
  - `real`, 199's realistic draw, which carries 34 essential targets;
  - `efp`, 238's essential-free pessimistic draw, K562 146/300;
  - `efr`, 238's essential-free realistic draw, K562 300/300.
- No `vcc`, no `cell-eval2`, no builder run, no git. Nothing written outside `experiments/438-october-stack-repricing/`.
- **Memory:** every job ran strictly serially behind the memory gate. The gate never read below 78 % free / 26 GiB
  free+inactive. Max resident set per process was ≤ 8.80 GB, with **two exceptions, both in the first hour, both
  fixed:**
  1. The very first run streamed four donor banks in one process to build the panel-independent bank statistics. Its
     max RSS was **11.19 GB**, mostly file-backed memmap pages. 062's own memory guard stopped it at 10.4 GiB. Those
     statistics are now cached (`work/bankstats/`) and never recomputed.
  2. The first complete A/B/C fold held all 42 compositions in RAM and peaked at **10.09 GB**. It was re-run after
     `lib438.Finals` (an LRU of compositions); every later fold is ≤ 7.3 GB.
  - Both are deviations from the mandate's 10 GB cap, reported here as AGENTS 9 requires.
  - **A third, smaller deviation from "strictly serial":** the eight one-fold `diag438.py` probes (≤ 8.8 GB each, run one
    at a time) overlapped the first ~5 min of the day-of rehearsal's full-panel stage (8.7 GB). The gate read
    ≥ 28 GiB free+inactive throughout; the rehearsal's first-stage timing includes that contention.

## 0. The answer

1. **(b) Known answers: every one reproduces exactly before any October number is read** (§1).
   - The mandate's four stack rows (`PP|K3+sh`, `b:h1x2|K3+sh`, `a:cd4x1.5|K3+sh`, the stack) and `PP|K3` reproduce
     420's `comp420-<fold>.json` with **|d| = 0** on all four folds. The comparison covers design and ballast-realised
     KMAX 16/24/32, cos_agg, pds, board-emulated pds at λ 0.5 and 1, and render A/B per frame.
   - These are also bit-identical, through the panel-agnostic code: the fold arrays, donor rows, PP composition, and K3
     and LODO keys. Energies match to 2e-16.
   - The pricer reproduces 420's r417a one-step prices to 4e-17.
2. **(d) CD4 ×1.5 does not hold on any October-like panel; it is the half that breaks.** Levelled k = 1, CD4 ×1.5 on
   top of H1 ×2 is worth:
   - A/B/C **+0.0017** (3/3);
   - `real` **−0.0016** (1/3);
   - `efp` **−0.0011** (0/3);
   - `efr` **−0.0024** (0/3, CI [−0.0037, −0.0010]).

   It survives leave-top-3-out, the Arc-like exclusion and the X-retention sensitivity (§4.1, §4.3). The October CD4 optimum is
   ×1 (`real`, `efp`) or ×0.75 (`efr`).
   - **A panel effect, not a rules effect.** A/B/C's panel under the D/E/F rules prices within 0.00014 of A/B/C on every
     levelled design, median 0.00004 (§4.4).
   - **Not the pinned deflation basis either.** Re-fitting it on the panel changes nothing (§4.5).
   - **Mechanism.** On the draws, a CD4 residual row points at the truth 0.52–0.79× as well as the donor base does,
     against 0.89× on A/B/C. More of it dilutes the base and costs pds.
   - **Without a day-of measurement,** CD4 ×1 is the lowest-regret dose across the four panels (max 0.0017, mean 0.0006);
     ×1.5 is the worst (0.0031 / 0.0015).
3. **(d) The H1 ×2 half survives the H1 union's coverage collapse** (77 carrying rows on A/B/C → 59 / 28 / 35). Levelled
   k = 1 vs the 417a analogue:
   - A/B/C +0.0024;
   - `real` +0.0051, inflated: 78–96 % of its cosine gain sits on the panel's 34 essential targets;
   - `efp` +0.0026 (3/3);
   - `efr` +0.0003 (1/3, flat).

   It is fragile on `efp`: 70 % of it comes from five near-essential targets (TAF11, TBPL1, SNRPB2, PRCC, ACAT2). It is
   never negative. More H1 than ×2 is panel-dependent:
   - robustly better on `real` (×2 → ×3 +0.0036; +0.0023 without its essential targets);
   - better but not robust on `efp` (+0.0012 → −0.0003 without its top-3);
   - ≤ 0 on `efr` (−0.0002);
   - A/B/C: +0.0005, not significant.

   **Keep ×2 as the default.** The day-of rule moves to ×2.5–×3 only when the move survives leave-top-3-out on the day's
   panel. It did on `real` (×3) and on A/B/C (×2.5), not on `efp` or `efr`.
4. **(d) q and the key keep their sign; q's value roughly halves.**
   - **q.** The q 1.5 step on the stack is +0.0041 to +0.0049 (k = 1), or +0.0019 to +0.0022 at 417a's calibrated
     transfer, 3/3 on every October panel. A/B/C: +0.0085 / +0.0044.
   - **Key.** The LODO key at clip (0.5, 3) vs the 410 key at (0.3, 3), on the stack: +0.0031 (`efp`), +0.0032 (`real`),
     +0.0017 (`efr`). A/B/C: +0.0027.
5. **What keeping A/B/C's weights costs, levelled k = 1.**
   - Against H1 ×2, CD4 ×1 (the day-of rule's October default and its choice on `efp` and `efr`): `efp` −0.0011 (3/3),
     `real` −0.0016 (2/3), `efr` −0.0024 (3/3).
   - Against the rule's choice on `real` (H1 ×3, CD4 ×1): −0.0051.
   - On `efr` the A/B/C stack is worse than no stack at all: **−0.0021 vs the 417a analogue**, 0/3, CI
     [−0.0037, −0.0005].
6. **Instrument checks that make those numbers usable (§2).**
   - **The composition follows the October driver's rules.** The D/E/F constants reproduce the builder's own meta on
     432's builds to ≤ 2.6e-5, and the re-fitted keys equal 432's artifacts bit for bit.
   - **The cross-panel energy route predicts 432's measured October render energies.** Error: +0.2 to +0.6 % on the 417a
     analogue (`real` b0), +0.9 to +1.2 % on the stack (`real` b1), +1.8 to +2.0 % on the stack (`efp` ef_b1).
   - **The levelled route used for October barely depends on that energy anchor.** Moving the anchor by the measured
     2 % moves levelled prices by ≤ 0.00005. Correcting the stack's 0.6 % step error moves them by ≤ 0.00017.
     Unlevelled prices move by up to 0.0010, on `efp`, whose context D sits past its cliff.
7. **(e) The day-of command** (§5):

       uv run python experiments/438-october-stack-repricing/dayof438.py --bundle $B --op $OP --cd4-root $CD4 --tag def_$N --board-q <q>

   - It runs after 432 §8 step 1 and before step 6.
   - **Rehearsed end to end on `efp`:** 1,650 s (27.5 min) wall clock, max RSS 8.72 GB, strictly serial. It reproduced
     the chain runs exactly: energies and PLOLO statistics |d| = 0.
   - On a 300/300 panel expect ~35 min, since each fold takes ~10 min.
   - It prints the re-priced stack table, the levelled H1 × CD4 surfaces and a RECOMMENDATION line, which is **step 6's
     dump recipe**: `VCC_H1_UNION_G=<0.28 × h>` [+ `VCC_CD4_TARGET_G=<0.3443655278649641 × c>` when c ≠ 1]. It also
     prints a q veto flag. q itself stays what the A/B/C board settles.
   - **The rule** (§5): the incumbent is H1 ×2, CD4 ×1, the lowest-regret point. It moves only if a grid point beats that
     by ≥ 0.0010, on 3/3 folds, and after leave-top-3-out.
   - **What the rule chose on the five panels:**
     - A/B/C, and A/B/C under DEF rules: H1 ×2.5–×3, **CD4 ×1.5**. It recovers A/B/C's answer.
     - `real`: H1 ×3, CD4 ×1.
     - `efp`, `efr`: **H1 ×2, CD4 ×1**.
     - These choices are identical under held-out or day-of operators and at q 1.25 or 1.5.
   - **If the command cannot run on the day,** dump H1 ×2 alone (`VCC_H1_UNION_G=0.56`, CD4 at the DEF_ENV default), not
     432's `VCC_CD4_TARGET_G=0.5165…`.

## 1. Known answers (every script aborts on failure; all run before any October number)

| # | row | expected | measured |
|---|---|---|---|
| KA-ROWS | panel-agnostic donor rows from cached bank statistics == 415's `full-<donor>.npz` (3 donors + PerturbMatch) | bit-identical | pass, 4 folds |
| KA-FOLD | my fold builder (targets, pre, truth, b_ctrl, gate, scored, wcols) == 420's cached `fold-<fold>.npz` | bit-identical | pass, 4/4 folds |
| KA-DONORS | `Donors438` (rows, availability, normrel medians, PM correction) == `sens415.Donors` | bit-identical | pass, 4/4 folds |
| KA-PP | `compose(rules abc, m = 1)` == `sens415.compose_pp` (final rows + merged residual) | bit-identical | pass, 4/4 folds |
| KA-KEYS | fold K3 and LODO-V1 weights == 418's `W3` / `V1__PP` | bit-identical | pass, 4/4 folds |
| KA-418 | `PP\|none` design == `comp418-<fold>.json` | \|d\| = 0 | pass, 4/4 folds |
| **KA-420 (the mandate's (b))** | `PP\|K3`, `PP\|K3+sh1.25`, `h1x2\|K3+sh1.25`, `cd4x1.5\|K3+sh1.25`, `stack\|K3+sh1.25` == `comp420-<fold>.json`: design + realised KMAX 16/24/32, cos_agg, pds, pds λ 0.5/1, render A/B per frame, capped fraction | \|d\| = 0 | **\|d\| = 0 on all 5 rows × 4 folds** (k562, hct116, hek293t, kolf2) |
| KA-REAL | my row scorer == `lib418.score` on every KA row | 1e-12 | pass |
| KA-P | my pricer on my A/B/C jsons == `analyse420.json` r417a one-step prices (6 rows × 4 folds × k 1 / 1.407) | 1e-12 | **4.2e-17** |
| KA-COMP (full panel) | PP constants (whitening medians, operator preserve, CD4/H1/partial fractions, post-add preserve, deflation energy) == the builder's meta on the same bundle | ≤ 5e-5 rel | A/B/C 8.87e-6 (vs 410b; 418's own KA value), `real` 2.57e-5 (vs 432 b0), `efp` 8.93e-6 (vs 432 ef_b1); deflation energy exact on all |
| KA-KEY (full panel) | the 410 key re-fitted on the panel's 3-donor operator == the artifact the builder used | identical w | A/B/C `d7f76a22`, `real` `4f1f0209`, `efp` `c3e429dd`: identical; A/B/C LODO-V1 == 415pp key (1e-15; w 2e-14) |
| KA-E420 | A/B/C design + realised energies of 6 rows == `abc420.json` | 1e-12 rel | 2.2e-16 |
| KA-CLAIMS | ntclevel claims from the synthetic bundles' own NTC cells == 432's ntclevel output | equal | 0.233754 / 0.239678 / 0.266129 on all three measured builds |
| KA-RESUME | a resumed run (grid extension H1 ×4/×5) re-builds the 417a-analogue row == the stored one | \|d\| = 0 | pass on every resumed fold and full panel |
| KA-BOOT | the bootstrap's identity draw == the json statistics | 1e-12 | pass on every bootstrapped row |
| CTRL | the LODO key permuted over targets (PP, q 1.25), unlevelled | negative | −0.0223 (A/B/C, 420: −0.0222), −0.0235 (`real`), −0.0149 (`efp`), −0.0048 (`efr`); 0/3 everywhere |

## 2. The instrument

### 2.1 What follows the panel (`lib438.bind`)
- **The panel** comes from the bundle's `pert_counts.csv`, and the gene axis is asserted against `gene_names.csv`.
- **Frames** (b_ctrl, gate) come from the bundle's own control cells. For a synthetic bundle they are asserted equal to
  A/B/C's, because its D/E/F controls are A/B/C relabelled.
- **PLOLO folds.** Targets = panel ∩ held-library coverage. Each fold uses the held-out operator fitted for that
  (panel, fold), from 233 for `real` and from 238 for `efp`/`efr`. Its axis is asserted through
  `plolo.Operator._check_axis`, so the Oct-22 trap raises.
- **Donor rows, normrel medians and the PM correction** use `Donors438`, which is sens415's arithmetic for any target
  list.
- **CD4 banks** are the panel's own extraction. The H1 union is masked to the panel's columns.
- **Keys.** The 410 key and 415's LODO-V1 are re-fitted per fold and on the full panel.
- **Every module global the old code read at call time is re-pointed:** `plolo.PANEL/PANEL_SET/CD4_ROOT/WORK`,
  `sens415.*`, `resid218.PANEL_SET`, `arms414.FRAMES`.

### 2.2 The composition rules follow the October driver
- On A/B/C (`ABC_ENV`): CD4 is sized on all rows at f 0.34 × m, and deflated by SVD-5 re-fitted on the rows. This is
  420 exactly.
- Everywhere else (`DEF_ENV`): CD4 is sized on carrying rows at g = 0.3443655278649641 × m (`VCC_CD4_TARGET_G`), and
  deflated by the pinned A/B/C gene basis `V5_abc`.
- CD4 ×1.5 on the day is `VCC_CD4_TARGET_G=0.5165482917974462` (432 T6).
- H1 union and partial rows are carrying-sized in both.
- `abcdef` is A/B/C's panel under the D/E/F rules. It separates "the rules changed" from "the panel changed".

### 2.3 Pricing (420's r417a route, with the panel's quantities)
- overall = (2.2138 Δpds_bd − 1.051 Δmse_raw)/6, with 187's cliff identity and 414's board-emulated pds at λ 0.5.
- k = 1 and 1.407 for the direction transfer.
- The four panel-dependent inputs:
  - **Energies.** An arm's A_c/B_c = 417a's measured A/B/C energy × (this panel's full-panel realised render energy of
    the arm ÷ A/B/C's of 417a). This is 420's anchor carried across panels.
  - **Claims.** 187's on A/B/C. On a D/E/F bundle, ntclevel417's f = 1 claim (Cp46 × Σw from the bundle's own NTC
    cells), which is what the day's re-level uses.
  - **X.** 410b's board-solved alignment × ret_f, where ret_f = cos_agg(PP|K3) on this panel's fold f ÷ the same on
    A/B/C's fold f. This is 233/238's retention assumption; ret_f is 1.00 on A/B/C, 0.85–1.01 on `real`, 0.53–0.78 on
    `efp` and 0.32–0.35 on `efr`. Pricing with ret = 1 moves no decision (§4.3).
  - **β (board-emulated pds).** 414's A/B/C artifact pattern, zeroed on the panel, at the same share of the fold's PP|K3
    row energy as on A/B/C.
- **LEVELLED prices, the D/E/F default** (432: the re-level is mandatory).
  - Each design gets ntclevel's factor s_c = sqrt(ρ* claim_c / B_c), with ρ* = exp(2 sd Φ⁻¹(0.148)): 0.9793 at sd 0.01
    (D/E/F), 0.9958 at sd 0.002 (A/B/C).
  - Its pds/cos are re-scored on PLOLO at that factor, and its energies re-rendered at that factor.
  - One ntclevel step lands within ~0.4 % of ρ\* for factors ≤ 1.07 and up to 1.2 % short at factors ≥ 1.09 (A/B/C's q 1.5
    stack). That is the render non-linearity 432 measured, and it is inside the price.
- **Unlevelled prices** are 420's primary route, reported for continuity.

### 2.4 Validating the cross-panel energy route on 432's measured October builds

ρ = B / claim per context; predicted by the route above with no artifact of the panel, measured by `ntclevel417` on 432's
builds (`work/<panel>/full438.json` "validation").

| build | design | measured ρ D/E/F | predicted ρ D/E/F | error in B |
|---|---|---|---|---|
| `real` b0 (433h − LOAD − CTX) | 417a analogue | 0.9420 / 0.8716 / 0.8201 | 0.9475 / 0.8738 / 0.8220 | +0.58 / +0.25 / +0.24 % |
| `real` b1 (+ LOAD stack) | stack | 0.9916 / 0.9137 / 0.8586 | 1.0036 / 0.9216 / 0.8660 | +1.22 / +0.87 / +0.86 % |
| `efp` ef_b1 (LOAD stack) | stack | 1.0712 / 0.9925 / 0.9351 | 1.0924 / 1.0109 / 0.9519 | +1.97 / +1.86 / +1.80 % |

- **The route over-states the stack's energy step by about 0.6 %.** The b1/b0 ratio is predicted 1.059 / 1.055 /
  1.054; measured 1.053 / 1.048 / 1.047.
- **On `efp` it sits a further ~1.3 % high.**
- **The unlevelled October prices therefore carry a ~1–2 % energy uncertainty.** On `efp` that sits past the D cliff,
  which is why the levelled route is primary.
- **In the levelled route every arm is placed at the same ρ\*,** so a common offset cancels. Measured on `real` and
  `efp`:
  - a uniform −2 % anchor shift moves levelled prices by ≤ 0.00005, unlevelled by ≤ 0.0010;
  - placing every non-PP design 0.6 % higher at level, the measured stack-step error, moves levelled prices by
    ≤ 0.00017 and leaves every CD4 row unchanged at the 5th decimal.

### 2.5 `efr` has no 3-donor operator: the fold-route energy surrogate

A 3-donor fit peaks at ~13 GB, over this mandate's cap. Two surrogates were tested on the three panels that do have
the real operator (`surr438.py`, `work/surr438.json`):

| surrogate | anchored B error (mean) | arm/417a ratio error (max) | price \|d\| max, unlevelled | levelled |
|---|---|---|---|---|
| same-axis held-out-Orion operator, s1 (abc / real / efp) | −20 / −47 / −56 % | 19 / 66 / 77 % | 0.0040 / 0.0070 / 0.0068 | 0.0007 / 0.0019 / 0.0016 |
| same-axis held-out-Orion operator, s2 | −48 / −59 / −69 % | 81 / 115 / 131 % | 0.0096 / 0.0091 / 0.0080 | 0.0031 / 0.0032 / 0.0030 |
| **fold route** (each design's realised energy averaged over the three cancer folds) | −0.4 / −0.6 / +19 % | 2.8 / 5.9 / 9.2 % | **0.0006 / 0.0009** / 0.0113 | **0.0006 / 0.0011** / 0.0057 |

- **The fold route is exact enough wherever every fold carries the whole panel,** as on A/B/C (269–294) and `real`
  (300). It fails on `efp`, whose K562 fold carries 146 targets.
- **`efr`'s folds carry 300/300,** so `efr` is priced on the fold route, carrying `real`'s error bar (≤ 0.0011
  levelled).
- **`efr`'s levels are fold-route factors** (`foldlevel438.py`).
- **The fold route keeps `real`'s and A/B/C's surface orderings.** CD4 ×1 > ×1.5 and the H1 axis ordering are identical;
  only magnitudes shift, by ≤ 0.0011.

## 3. Tables (c)

Every entry is Δoverall vs that panel's own 417a analogue, `h1c1|K3+sh1.25`: PP composition under the panel's rules,
the panel's re-fitted 410 key clip (0.3, 3), q 1.25. The unit is the cancer-fold mean at k = 1; "n+" is folds
positive; brackets give the 95 % bootstrap CI of the mean (targets resampled within folds, 2,000 draws).

Folds: k562 / hct116 / hek293t; `efp`'s k562 fold carries 146 targets. Coverage (full panel):

| panel | K562 / HCT116 / HEK293T | H1 union rows (gse295214_h1 / h1_2025) | CD4 complete rows | V5_abc captured share of centred CD4 energy |
|---|---|---|---|---|
| abc | 269 / 290 / 294 | 77 (64 / 20) | 276 | 0.198 (own SVD-5 under ABC rules) |
| real | 300 / 300 / 300 | 59 (57 / 7) | 286 | 0.117 |
| efp | 146 / 293 / 296 | 28 (27 / 4) | 182 | 0.133 |
| efr | 300 / 300 / 300 | 35 (30 / 6) | 290 | 0.125 |

### 3.1 The stack and its halves, LEVELLED (the D/E/F default)

| | abc | abcdef | real | efp | efr |
|---|---|---|---|---|---|
| H1 x2 half | +0.00238 3/3 [+0.0011, +0.0035] | +0.00239 3/3 [+0.0011, +0.0035] | +0.00511 3/3 [+0.0036, +0.0065] | +0.00255 3/3 [+0.0011, +0.0040] | +0.00031 1/3 [-0.0008, +0.0013] |
| CD4 x1.5 half | +0.00183 3/3 [+0.0005, +0.0030] | +0.00189 3/3 [+0.0006, +0.0031] | -0.00149 1/3 [-0.0028, -0.0000] | -0.00096 0/3 [-0.0023, +0.0003] | -0.00242 0/3 [-0.0037, -0.0011] |
| **stack (H1 x2 + CD4 x1.5)** | **+0.00406** 3/3 [+0.0024, +0.0055] | +0.00409 3/3 [+0.0025, +0.0055] | +0.00351 3/3 [+0.0017, +0.0053] | +0.00143 2/3 [-0.0006, +0.0034] | **-0.00205** 0/3 [-0.0037, -0.0005] |
| stack, k = 1.407 | +0.00606 3/3 | +0.00611 3/3 | +0.00549 3/3 | +0.00250 2/3 | -0.00187 0/3 |
| CD4 x1.5 on top of H1 x2 | +0.00170 3/3 [+0.0004, +0.0029] | +0.00172 3/3 [+0.0004, +0.0029] | -0.00160 1/3 [-0.0029, -0.0002] | -0.00111 0/3 [-0.0024, +0.0002] | -0.00236 0/3 [-0.0037, -0.0010] |
| H1 x2 on top of CD4 x1.5 | +0.00226 3/3 [+0.0011, +0.0034] | +0.00223 3/3 [+0.0010, +0.0033] | +0.00500 3/3 [+0.0036, +0.0063] | +0.00239 3/3 [+0.0010, +0.0037] | +0.00037 2/3 [-0.0007, +0.0013] |

Per fold (k562 / hct116 / hek293t), levelled k = 1:

| | abc | real | efp | efr |
|---|---|---|---|---|
| H1 x2 half | +0.0036 / +0.0029 / +0.0007 | +0.0074 / +0.0038 / +0.0041 | +0.0048 / +0.0028 / +0.0000 | +0.0018 / -0.0009 / -0.0000 |
| CD4 x1.5 half | +0.0032 / +0.0019 / +0.0005 | +0.0008 / -0.0020 / -0.0032 | -0.0018 / -0.0003 / -0.0008 | -0.0021 / -0.0019 / -0.0033 |
| stack | +0.0063 / +0.0047 / +0.0012 | +0.0078 / +0.0019 / +0.0009 | +0.0024 / +0.0025 / -0.0006 | -0.0003 / -0.0027 / -0.0031 |
| CD4 x1.5 on top of H1 x2 | +0.0027 / +0.0019 / +0.0005 | +0.0004 / -0.0019 / -0.0033 | -0.0024 / -0.0004 / -0.0006 | -0.0021 / -0.0019 / -0.0031 |

### 3.2 The same, UNLEVELLED (420's primary route)

| | abc | abcdef | real | efp | efr |
|---|---|---|---|---|---|
| H1 x2 half | +0.00254 3/3 [+0.0013, +0.0037] | +0.00256 3/3 | +0.00587 3/3 [+0.0044, +0.0072] | +0.00026 2/3 [-0.0012, +0.0017] | +0.00014 1/3 [-0.0009, +0.0011] |
| CD4 x1.5 half | +0.00172 3/3 [+0.0004, +0.0029] | +0.00179 3/3 | -0.00097 1/3 [-0.0023, +0.0005] | -0.00019 2/3 [-0.0015, +0.0011] | -0.00257 0/3 [-0.0038, -0.0013] |
| **stack** | **+0.00369** 3/3 [+0.0021, +0.0051] | +0.00380 3/3 | +0.00437 3/3 [+0.0026, +0.0062] | +0.00042 2/3 [-0.0016, +0.0024] | -0.00230 0/3 [-0.0038, -0.0009] |
| stack, k = 1.407 | +0.00563 3/3 | +0.00575 3/3 | +0.00623 3/3 | +0.00150 2/3 | -0.00211 0/3 |
| CD4 x1.5 on top of H1 x2 | +0.00117 3/3 [-0.0001, +0.0024] | +0.00125 3/3 | -0.00150 1/3 [-0.0028, -0.0001] | +0.00016 2/3 [-0.0012, +0.0014] | -0.00244 0/3 [-0.0038, -0.0012] |
| CONTROL: LODO key permuted (PP) | -0.02225 0/3 | -0.02225 0/3 | -0.02351 0/3 | -0.01489 0/3 | -0.00476 0/3 |

Per fold, unlevelled k = 1: A/B/C's stack +0.0059 / +0.0044 / +0.0009 (420: identical); `efp` +0.0014 / +0.0015 /
−0.0016; `efr` −0.0006 / −0.0030 / −0.0033.

`efp`'s unlevelled route has context D past its cliff (ρ 1.08 for the 417a analogue, 1.09 for the stack; 432 measured
1.07 for the stack), so every energy-raising step pays the cliff there. That is why unlevelled `efp` reads the H1 half
at +0.0003 and the levelled route at +0.0026.

### 3.3 Stack neighbours, vs A/B/C's stack (H1 ×2, CD4 ×1.5, q 1.25), LEVELLED k = 1

| | abc | abcdef | real | efp | efr |
|---|---|---|---|---|---|
| H1 x0, CD4 x1.5 | -0.00803 0/3 | -0.00802 0/3 | -0.01272 0/3 | -0.00628 0/3 | -0.00196 1/3 |
| H1 x1, CD4 x1.5 | -0.00229 0/3 | -0.00226 0/3 | -0.00510 0/3 | -0.00243 0/3 | -0.00037 1/3 |
| H1 x1.5, CD4 x1.5 | -0.00085 0/3 | -0.00085 0/3 | -0.00224 0/3 | -0.00113 0/3 | -0.00011 2/3 |
| H1 x2.5, CD4 x1.5 | +0.00048 2/3 | +0.00045 2/3 | +0.00193 3/3 | +0.00085 2/3 | -0.00008 1/3 |
| H1 x3, CD4 x1.5 | +0.00062 2/3 | +0.00058 2/3 | +0.00350 3/3 | +0.00138 2/3 | -0.00014 1/3 |
| H1 x4, CD4 x1.5 | +0.00012 2/3 | +0.00005 2/3 | +0.00573 3/3 | +0.00213 2/3 | -0.00040 1/3 |
| H1 x5, CD4 x1.5 | -0.00113 2/3 | -0.00121 2/3 | +0.00699 3/3 | +0.00244 2/3 | -0.00066 1/3 |
| H1 x2, CD4 x0.75 | -0.00385 0/3 | -0.00392 0/3 | +0.00151 2/3 | +0.00074 2/3 | **+0.00306** 3/3 |
| **H1 x2, CD4 x1** | -0.00172 0/3 | -0.00174 0/3 | **+0.00160** 2/3 | **+0.00111** 3/3 | **+0.00236** 3/3 |
| H1 x2, CD4 x1.25 | -0.00058 1/3 | -0.00061 1/3 | +0.00106 3/3 | +0.00072 3/3 | +0.00128 3/3 |
| H1 x2, CD4 x1.75 | -0.00005 2/3 | +0.00000 1/3 | -0.00142 0/3 | -0.00097 0/3 | -0.00154 0/3 |
| H1 x2, CD4 x2 | -0.00042 1/3 | -0.00038 1/3 | -0.00302 0/3 | -0.00202 0/3 | -0.00312 0/3 |

The unlevelled neighbour table (`work/summary438.md`) has the same orderings on the CD4 axis. On the H1 axis, A/B/C's
×2.5 / ×3 read +0.0004 / +0.0005 at 2/3 unlevelled, 418's "H1 peaks at ×2" (it pays the energy cliff unlevelled).

### 3.4 Sharpening q on the stack (q step vs the stack at q 1.25)

| | abc | abcdef | real | efp | efr |
|---|---|---|---|---|---|
| q 1.35, levelled k = 1 | +0.00388 3/3 | +0.00387 3/3 | +0.00230 3/3 | +0.00194 3/3 | +0.00194 3/3 |
| **q 1.5, levelled k = 1** | **+0.00847** 3/3 | +0.00850 3/3 | +0.00492 3/3 | +0.00412 3/3 | +0.00454 3/3 |
| q 1.75, levelled k = 1 | +0.01311 3/3 | +0.01319 3/3 | +0.00746 3/3 | +0.00607 3/3 | +0.00758 3/3 |
| q 2, levelled k = 1 | +0.01486 3/3 | +0.01490 3/3 | +0.00781 3/3 | +0.00655 3/3 | +0.00889 3/3 |
| **q 1.5, levelled, 417a-calibrated (k 0.65, k_pds 0.42)** | **+0.00436** 3/3 | +0.00437 3/3 | +0.00221 3/3 | +0.00189 3/3 | +0.00211 3/3 |
| q 1.75, levelled, calibrated | +0.00648 3/3 | +0.00651 3/3 | +0.00294 3/3 | +0.00258 3/3 | +0.00339 3/3 |
| q 2, levelled, calibrated | +0.00680 3/3 | +0.00682 3/3 | +0.00245 3/3 | +0.00242 3/3 | +0.00377 3/3 |
| q 1.5, unlevelled k = 1 | +0.00627 3/3 | +0.00628 3/3 | +0.00309 3/3 | +0.00688 3/3 | +0.00431 3/3 |
| q 1.5, unlevelled, calibrated | +0.00212 3/3 | +0.00212 3/3 | +0.00023 2/3 | +0.00460 3/3 | +0.00176 3/3 |

The calibrated transfer is 417a's board read (FINDINGS 0: k_sh 0.65 from mse, 0.42 from pds). Levelled, q 1.5's step on
October panels is 0.49–0.58× A/B/C's at k = 1 and 0.43–0.51× at the calibrated transfer. It is positive 3/3 and
survives leave-top-3-out everywhere (§4.1).

Unlevelled, `efp` reads it larger (+0.0069) because sharpening lowers energy, and `efp`'s D context sits past its cliff.
The re-level refunds that on the day, so the levelled rows are the ones to read.

### 3.5 Keys on the stack (vs the stack with the 410 key clip (0.3, 3)), LEVELLED k = 1

| | abc | abcdef | real | efp | efr |
|---|---|---|---|---|---|
| **LODO V1 clip (0.5, 3), q 1.25** | +0.00274 3/3 [+0.0016, +0.0041] | +0.00272 3/3 | +0.00318 2/3 [+0.0014, +0.0060] | +0.00306 3/3 [+0.0016, +0.0050] | +0.00169 3/3 [+0.0002, +0.0034] |
| LODO V1 clip (0.3, 3), q 1.25 | +0.00189 3/3 | +0.00186 3/3 | +0.00308 3/3 | +0.00269 3/3 | +0.00106 2/3 |
| 410 key clip (0.5, 3), q 1.25 | +0.00096 3/3 | +0.00095 3/3 | +0.00017 2/3 | +0.00039 2/3 | +0.00111 3/3 |
| LODO V1 clip (0.5, 3), q 1.5 (vs K3 q 1.5) | +0.00329 3/3 | +0.00326 3/3 | +0.00343 2/3 | +0.00349 3/3 | +0.00187 3/3 |
| LODO V1 clip (0.5, 3) on PP, q 1.25 | +0.00292 3/3 | +0.00294 3/3 | +0.00274 3/3 | +0.00254 3/3 | +0.00198 3/3 |

Unlevelled, the LODO (0.5, 3) step on the stack is +0.0028 (A/B/C; 420: +0.0031 on PP) and +0.0032 (`real`).

### 3.6 The H1 × CD4 surface (K3), LEVELLED k = 1: argmax, and what A/B/C's weights cost

| panel | q | argmax (H1, CD4) | argmax | A/B/C weights (×2, ×1.5) | argmax − A/B/C [paired 95 % CI] (n+) | best H1 at CD4 ×1.5 | best CD4 at H1 ×2 |
|---|---|---|---:|---:|---|---|---|
| abc | 1.25 | x3, x1.75 | +0.00475 | +0.00406 | +0.00073 [-0.0003, +0.0016] (2/3) | x3 | x1.5 |
| abc | 1.5 | x3, x1.75 | +0.01304 | +0.01246 | +0.00064 [-0.0004, +0.0016] (2/3) | x3 | x1.5 |
| abcdef | 1.25 | x3, x1.75 | +0.00479 | +0.00409 | +0.00074 [-0.0003, +0.0016] (2/3) | x3 | x1.5 |
| real | 1.25 | x5, x1 | +0.01150 | +0.00351 | +0.00816 [+0.0050, +0.0112] (3/3) | x5 | x1 |
| real | 1.5 | x5, x1 | +0.01690 | +0.00841 | +0.00870 [+0.0056, +0.0119] (3/3) | x5 | x1 |
| efp | 1.25 | x5, x1 | +0.00458 | +0.00143 | +0.00321 [-0.0004, +0.0066] (2/3) | x5 | x1 |
| efp | 1.5 | x5, x1.25 | +0.00878 | +0.00554 | +0.00332 [+0.0001, +0.0066] (2/3) | x5 | x1 |
| efr | 1.25 | x2, x0.75 | +0.00102 | -0.00205 | +0.00306 [+0.0010, +0.0052] (3/3) | x2 | x0.75 |
| efr | 1.5 | x2.5, x0.75 | +0.00522 | +0.00248 | +0.00275 [+0.0005, +0.0050] (3/3) | x2.5 | x0.75 |

- The grid is H1 ×{0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5} × CD4 ×{0.75, 1, 1.25, 1.5, 1.75, 2} at q 1.25 and 1.5.
- The full surfaces per panel are in `work/recommend438-<panel>.txt`.
- The `real`/`efp` argmax sits at the grid's H1 edge (×5). §4.2 explains why that is not taken.

### 3.7 Energies (ρ = B/claim per context slot; unlevelled; the anchored route)

| panel | 417a analogue | H1 x2 half | CD4 x1.5 half | stack | stack q 1.5 | stack after one ntclevel step |
|---|---|---|---|---|---|---|
| abc | 0.991/0.922/0.880 | 1.001/0.933/0.890 | 0.994/0.919/0.878 | 0.987/0.915/0.873 | 0.914/0.842/0.807 | 0.996/0.993/0.992 |
| abcdef | 0.990/0.921/0.879 | 1.001/0.933/0.891 | 0.994/0.919/0.877 | 0.988/0.915/0.874 | 0.914/0.843/0.807 | 0.996/0.993/0.992 |
| real | 0.947/0.874/0.822 | 0.988/0.912/0.858 | 0.973/0.892/0.839 | 1.004/0.922/0.866 | 0.914/0.840/0.786 | 0.980/0.977/0.974 |
| efp | **1.080**/1.004/0.946 | 1.103/1.026/0.966 | 1.067/0.986/0.929 | **1.092**/1.011/0.952 | 1.045/0.968/0.909 | 0.981/0.980/0.979 |
| efr (fold route) | 0.843/0.778/0.739 | 0.828/0.764/0.726 | 0.828/0.759/0.721 | 0.820/0.752/0.715 | 0.747/0.713/0.671 | 0.973/0.974/0.972 |

- **The A/B/C 417a analogue's ρ reproduces 420's measured anchor exactly.**
- **On `efp` both the 417a analogue and the stack sit past context D's cliff unlevelled.** 432 measured 1.071 for the
  stack. The re-level is mandatory there, as 432 found.
- **On `efr` everything sits far below its cliffs.** The re-level factors are 1.08–1.17, inside the builder's
  [0.8, 1.2].

## 4. (d) The answer, in detail

### 4.1 Robustness (`robust438.py`, levelled k = 1; `work/robust438_all.md`)

For each step: the price on all targets; without each fold's top-3 contributing targets (236's falsifier, the union
over folds); without the four essential screens' 2,389 targets (the Arc-like panel, 436); and without the K562-essential
list.

| step | panel | all | minus top-3 | Arc-like | minus K562-essential |
|---|---|---|---|---|---|
| CD4 x1.5 on top of H1 x2 | abc | +0.00170 3/3 | +0.00128 2/3 | +0.00170 | +0.00170 |
| | real | -0.00160 1/3 | -0.00197 0/3 | -0.00161 | -0.00172 |
| | efp | -0.00111 0/3 | -0.00189 0/3 | -0.00092 | -0.00111 |
| | efr | -0.00236 0/3 | -0.00278 0/3 | -0.00242 | -0.00236 |
| CD4 x1.5 → x1 at H1 x2 | efp | +0.00111 3/3 | +0.00069 2/3 | +0.00092 | +0.00111 |
| | efr | +0.00236 3/3 | +0.00203 3/3 | +0.00241 | +0.00236 |
| H1 x2 half | abc | +0.00238 3/3 | +0.00143 3/3 (SIN3B, SMARCA5, TAF4, UPF3B) | +0.00238 | +0.00238 |
| | real | +0.00511 3/3 | +0.00382 3/3 | +0.00319 | +0.00316 |
| | efp | +0.00255 3/3 | **+0.00074** 2/3 (ACAT2, PRCC, SNRPB2, TAF11, TBPL1) | +0.00254 | +0.00255 |
| | efr | +0.00031 1/3 | -0.00042 1/3 | +0.00093 | +0.00031 |
| H1 x2 → x3 (CD4 x1) | abc | +0.00047 2/3 | -0.00034 1/3 | – | – |
| | real | +0.00355 3/3 | +0.00276 3/3 | +0.00263 | +0.00231 |
| | efp | +0.00124 2/3 | **-0.00032** 1/3 | +0.00145 | +0.00124 |
| | efr | -0.00021 1/3 | -0.00070 1/3 | +0.00045 | -0.00021 |
| H1 x3 → x5 (CD4 x1) | abc | -0.00227 0/3 | – | – | – |
| | real | +0.00307 3/3 | +0.00228 3/3 | +0.00238 | +0.00169 |
| | efp | +0.00087 2/3 | -0.00162 1/3 | +0.00137 | +0.00087 |
| | efr | -0.00058 0/3 | -0.00133 0/3 | +0.00052 | -0.00058 |
| q 1.25 → 1.5 on the stack | abc / real / efp / efr | +0.0085 / +0.0049 / +0.0041 / +0.0045, all 3/3 | +0.0087 / +0.0047 / +0.0036 / +0.0041 | – / +0.0045 / +0.0053 / +0.0045 | – |
| key K3 → LODO (0.5, 3) on the stack | abc / real / efp / efr | +0.0027 / +0.0032 / +0.0031 / +0.0017 | +0.0021 / +0.0021 / +0.0021 / +0.0010 | – / +0.0041 / +0.0054 / +0.0017 | – |

**What carries H1 on `real`: essential targets (`ess438.py`).**
- On `real`, 78–96 % of the H1 ×2 half's summed Δx sits on the panel's 34 K562-essential targets (11 % of targets).
  With them removed its relative cosine gain falls from +0.08–0.10 to +0.04–0.07 on every fold.
- This is 238's finding for the October sources, now seen in H1: essential knockdowns transfer across lines, H1 hESC
  included.
- Arc's panel carries none, so `real` overstates H1.
- On `efp` the carriers are near-essential genes the K562 list misses: TAF11, TBPL1, SNRPB2. TAF11 and TBPL1 are in the
  four-screen union. SNRPB2 is spliceosome.

### 4.2 The direct answer

- **CD4 ×1.5 (417e's second half): does not hold.**
  - Levelled, it costs 0.0011 (`efp`), 0.0016 (`real`) and 0.0024 (`efr`) against CD4 ×1 on top of H1 ×2. It is
    positive on 0/3, 1/3 and 0/3 folds.
  - It gets worse with the top-3 targets removed, and is unchanged on the Arc-like panels and with ret = 1.
  - Its cost is almost entirely pds: Δpds_bd −0.0037 (`efp`), −0.0065 (`efr`) for the half, with Δrel cos ≈ 0. On A/B/C
    it cost −0.0007 pds and bought +3 % cosine.
  - **The optimum moves down, to CD4 ×1** (`real`, `efp`; ×0.75–×1.25 within 0.0005) **or ×0.75** (`efr`).
  - The move is a property of the panel. On A/B/C's own panel the D/E/F rules keep ×1.5 (§4.4).
  - The mechanism is in §4.5.
- **H1 ×2 (the first half): holds.**
  - **Does it survive the coverage collapse? Yes.** Levelled vs H1 ×1 it is +0.0026 (`efp`, 28 rows), +0.0051 (`real`,
    59 rows) and +0.0003 (`efr`, 35 rows), against +0.0024 on A/B/C (77 rows).
  - **Coverage lowers the number of rows H1 touches, not the value per row.** The H1 residual is carrying-sized, so ×2
    doubles the dose on whatever rows it has.
  - **Its gain is thinner and more concentrated on October panels:** `efp` loses 70 % of it to leave-top-3-out, and
    `efr` is flat.
  - **It is never significantly negative,** so it stays.
  - **More H1 is a panel-dependent October move, not a default.**
    - Levelled, the H1 axis keeps rising to ×5 on `real` and `efp`, the grid edge.
    - On `real` the rise is robust: ×2 → ×3 is +0.0036 (3/3); +0.0028 without each fold's top-3; +0.0023 without the
      essential targets. The ×2 *half* itself is essential-inflated there.
    - On `efp` the ×2 → ×3 and ×3 → ×5 steps turn negative without their top-3 targets.
    - On `efr` both steps are ≤ 0.
    - A/B/C's own levelled optimum is ×3 (+0.0006 over ×2, 2/3, not significant).
  - **Keep ×2 as the default.** The day-of rule (§5) moves to ×2.5–×3 only when the move survives leave-top-3-out on the
    day's panel. It chose ×3 on `real` and ×2.5 on A/B/C, and kept ×2 on `efp` and `efr`.
- **q 1.25 → 1.5: holds in sign, about half the size.**
  - Levelled k = 1: +0.0041 to +0.0049 on October panels vs +0.0085 on A/B/C.
  - At 417a's calibrated transfer: +0.0019 to +0.0022 vs +0.0044.
  - It is 3/3 everywhere and robust.
  - Whatever the A/B/C board decides for q (433k is q 1.5) carries to October. The October gain is about half.
- **Key clip (0.3, 3) vs 420's LODO key at clip (0.5, 3): the LODO key's advantage holds.** It is +0.0031 (`efp`), +0.0032
  (`real`) and +0.0017 (`efr`), against +0.0027 on A/B/C, and keeps a positive mean after leave-top-3-out (+0.0010 to
  +0.0021). It is not wired for October: it would need its own day-of artifact and pin. `full438` already computes the
  panel's LODO key on the day's operator.
- **What the October stack loses by keeping A/B/C's weights (H1 ×2, CD4 ×1.5)**, levelled k = 1:

  | | vs H1 ×2, CD4 ×1 (the rule's October default) | vs the rule's choice | vs the surface argmax (H1 ≤ ×3) |
  |---|---|---|---|
  | `efp` | −0.0011 (3/3) | −0.0011 (H1 ×2, CD4 ×1) | −0.0023 (×3, ×1) |
  | `real` | −0.0016 (2/3) | −0.0051 (H1 ×3, CD4 ×1) | −0.0051 (×3, ×1) |
  | `efr` | −0.0024 (3/3) | −0.0024 (H1 ×2, CD4 ×1) | −0.0031 (×2, ×0.75) |

  - On `efr` the A/B/C stack is **below the 417a analogue** (−0.0021, 0/3).
  - H1 ×2 alone (the rule's choice there) is +0.0003, and the surface's best point (×2, ×0.75) +0.0010.

### 4.3 Sensitivity to the X-retention assumption (ret = 1: the board's alignment term not scaled by the panel's PLOLO cosine)

| step (levelled k = 1) | real (ret / ret = 1) | efp | efr |
|---|---|---|---|
| H1 x2 half | +0.00511 / +0.00551 | +0.00255 / +0.00381 | +0.00031 / +0.00147 |
| CD4 x1.5 on top of H1 x2 | -0.00160 / -0.00166 | -0.00111 / -0.00105 | -0.00236 / -0.00266 |
| stack | +0.00351 / +0.00385 | +0.00143 / +0.00276 | -0.00205 / -0.00119 |
| key LODO (0.5, 3) on the stack | +0.00318 / +0.00332 | +0.00306 / +0.00465 | – |

- Removing the retention scaling raises every direction gain.
- The CD4 conclusion is unchanged to the fourth decimal: CD4's cost is pds, not X.
- H1 ×2 becomes clearly positive on `efr` (+0.0015). (`work/summary438_noret.md`.)

### 4.4 Rules or panel? `abcdef`: A/B/C's panel under the D/E/F composition rules

- The rules are carrying-row CD4 sizing at g 0.3443655 and the pinned V5_abc deflation.
- Every design on `abcdef` prices within 0.00014 of `abc` levelled: median 0.00004 over 153 designs; unlevelled
  within 0.00030, median 0.00007.
  - Stack +0.00409 vs +0.00406; CD4 ×1.5 on top of H1 ×2 +0.00172 vs +0.00170.
  - The argmax (×3, ×1.75) is the same, and so is the rule's choice (H1 ×2.5–×3 with CD4 ×1.5).
- On A/B/C's own rows V5_abc *is* the SVD-5 basis (captured 0.1982 = own top-5 0.1982). The carrying rule solves CD4 at
  f 0.3432 vs all-rows 0.3461.
- **So everything that moved on `real`/`efp`/`efr` moved because the panel changed.** The day-of rules are not the
  cause.

### 4.5 Why CD4 ×1.5 fails on October-like panels: not the pinned basis, but weaker per-row transfer of the CD4 residual

**Hypothesis 1, rejected: the pinned V5_abc basis removes too little of the shared T-cell hub.**
- V5_abc captures 0.12–0.13 of the October panels' centred CD4 energy, against 0.20 on A/B/C (250 measured 0.10–0.14).
- `efpsvd` is `efp` with the SVD-5 deflation re-fitted on its own rows, 250's rejected day-of re-fit. It removes 23.2 %
  of the energy vs the pinned basis's 13.6 %.
- CD4 ×1.5 on top of H1 ×2 is then **−0.00135** levelled (−0.0025 / +0.0001 / −0.0017), against −0.00111 with the
  pinned basis. Its pds cost is unchanged (Δpds_bd −0.0037 either way).
- The render's shared-mean share (A−B)/A does not rise with CD4 on any panel.

**What is measured instead (`diag438.py`): the CD4 residual points at the held-line truth less well, per carrying row,
than the donor base does.** Mean row cosine over the fold's scored genes (three folds k562 / hct116 / hek293t):

| panel | CD4 residual (deflated) | whitened donor base | H1 union | CD4 / base | H1 / base | CD4 rows with cos > 0 |
|---|---|---|---|---:|---:|---|
| abc | +0.0299 / +0.0207 / +0.0159 | +0.0298 / +0.0238 / +0.0211 | +0.0601 / +0.0457 / +0.0347 | **0.89** | 1.88 | 0.75 / 0.71 / 0.68 |
| real | +0.0217 / +0.0111 / +0.0092 | +0.0219 / +0.0143 / +0.0115 | +0.0944 / +0.0412 / +0.0394 | 0.88 | 3.67 | 0.67 / 0.62 / 0.59 |
| efp | +0.0139 / +0.0090 / +0.0060 | +0.0175 / +0.0105 / +0.0085 | +0.0728 / +0.0585 / +0.0201 | 0.79 | 4.14 | 0.56 / 0.62 / 0.58 |
| efr | +0.0068 / +0.0061 / +0.0030 | +0.0117 / +0.0104 / +0.0082 | +0.0369 / +0.0179 / +0.0155 | **0.52** | 2.32 | 0.53 / 0.59 / 0.56 |

- **On A/B/C a CD4 row is almost as good a pointer as the donor consensus** (0.89×). Adding 50 % more of it buys cosine
  (+3 %) at a small pds cost.
- **On the essential-free draws it is 0.52–0.79× the base,** and barely better than a coin flip on sign (53–62 % of rows
  positive). More of it dilutes the base and costs retrieval.
- **H1 rows stay 2–4× better than the base on every panel,** which is why H1 ×2 survives.
- **The caveat this raises.** Arc selects its panel for "adequate coverage in public datasets" and high basal
  expression. A/B/C's CD4 transfer may therefore be a selection effect that D/E/F shares, and the random draws may
  understate it. That is the case for pricing on the day rather than carrying either answer.
  - The day-of table is exactly that measurement.
  - `diag438.py` prints the CD4/base ratio in one minute per fold.

**Regret of each CD4 dose at H1 ×2 (levelled k = 1, q 1.25; shortfall vs the best CD4 dose on that panel):**

| CD4 | abc | real | efp | efr | max | mean |
|---|---:|---:|---:|---:|---:|---:|
| ×0.75 | 0.0038 | 0.0001 | 0.0004 | 0 | 0.0038 | 0.0011 |
| **×1** | 0.0017 | 0 | 0 | 0.0007 | **0.0017** | **0.0006** |
| ×1.25 | 0.0006 | 0.0005 | 0.0004 | 0.0018 | 0.0018 | 0.0008 |
| ×1.5 (417e) | 0 | 0.0016 | 0.0011 | 0.0031 | 0.0031 | 0.0015 |

At q 1.5 the picture is the same: max regret ×1 0.0016, ×1.25 0.0014, ×1.5 0.0027. **With no day-of measurement, CD4 ×1
(the DEF_ENV default) is the lowest-regret dose; ×1.5 is the worst of the four.**

## 5. (e) The day-of command — text for DEF-RUNBOOK (not applied)

> ### 5.x (438) Re-price the bulk stack's WEIGHTS on the day's panel — one command, between 432 §8 steps 1 and 6
>
> The September stack (417e) was priced on A/B/C. On three October-like panels its CD4 ×1.5 half LOSES (−0.0011 to
> −0.0024 levelled vs CD4 ×1) and its H1 ×2 half holds (438 §4). The weights must be re-priced on the day's panel, not
> carried.
>
> **When:** after 432 §8 step 1 (the light stages + `whitener`: the day's operator `$OP` and CD4 banks `$CD4` exist),
> before step 6 (the stack dump). It needs no key artifact, no pin and no build.
>
> ```bash
> uv run python experiments/438-october-stack-repricing/dayof438.py --bundle $B --op $OP --cd4-root $CD4 --tag def_$N \
>     --board-q <the q the A/B/C board settled: 1.25 or 1.5>
> #   27.5 min rehearsed (efp: full panel 4.5 min incl. the NTC claims; folds 3.4 / 9.7 / 9.9 min); ~35 min on a 300/300
> #   panel. <= 8.8 GB per process, strictly serial behind the memory gate.
> #   Output: work/def_$N/{full438.json,level438.json,comp438-*-prod.json,recommend438-prod.json}
> ```
>
> **It prints:** the stack / halves / q / key table (per fold, mean, 95 % CI, unlevelled and levelled), the levelled
> H1 × CD4 surface at q 1.25 and 1.5, the movers that pass the rule, and one line:
>
>     RECOMMENDATION (def_$N): H1 x<h>, CD4 x<c>, q <q> (the board's), key K3 (410 re-fitted on the day's operator)
>       dump recipe (432 section 8 step 6): VCC_H1_UNION_G=<0.28 h> [VCC_CD4_TARGET_G=<0.3443655278649641 c>] VCC_BULK_DUMP=$DUMP
>
> **What it changes in 432 §8:** only step 6's dump overrides (today `VCC_H1_UNION_G=0.56
> VCC_CD4_TARGET_G=0.5165482917974462`). Paste the printed override(s) instead. Check the dump log for `[h1-union] g=<0.28 h>` and
> `[cd4] sizing=carrying g=<0.344… c>`. When c = 1 the CD4 override is the DEF_ENV default; drop it rather than
> re-exporting it. Nothing else moves: the key (step 2) stays the driver's `wshare_lo0.3_hi3_def`, q (step 5) stays the
> board's, and the re-level (steps 8–9) is unchanged and still mandatory. If it prints `VETO q 1.5`, build at q 1.25.
>
> **The rule (fixed before Oct 22; `recommend438.py`):** levelled route, k = 1, at the board's q (`--board-q`), paired on
> the day's three folds.
> - The incumbent is H1 ×2, CD4 ×1: 417e's H1 half alone, the lowest-regret point over A/B/C and three October-like
>   draws (438 §4.5). 417e's CD4 ×1.5 is the highest-regret CD4 dose there.
> - Candidates are H1 ×{1, 1.5, 2, 2.5, 3} × CD4 ×{0.75 … 2}. H1 is capped at ×3.
> - A candidate moves only if it gains ≥ +0.0010 on the mean, on 3/3 folds, and keeps a positive mean after
>   leave-top-3-out. The most robust mover wins.
> - q 1.5 is vetoed only if its 417a-calibrated step is ≤ 0 on any fold.
>
> **What it chose on the rehearsal panels** (identical under held-out and day-of fold operators, and at q 1.25 and 1.5):
> - A/B/C: H1 ×2.5, **CD4 ×1.5**. On A/B/C's own panel it moves back to 417e's CD4 dose.
> - A/B/C under DEF rules: ×2.5–×3, ×1.5.
> - `real`: H1 ×3, CD4 ×1.
> - `efp`: **H1 ×2, CD4 ×1**.
> - `efr`: **H1 ×2, CD4 ×1**.
> - No q veto anywhere.
>
> **If it cannot run:** dump **H1 ×2 alone** (`VCC_H1_UNION_G=0.56`, CD4 at the DEF_ENV default). That is the
> incumbent: it beats A/B/C's stack by +0.0011 / +0.0016 / +0.0024 levelled on the three draws, and trails it by 0.0017
> on A/B/C.
>
> **Rehearsal (2026-09-23, 238's `efp` bundle, `--fold-ops prod`):** full panel 267 s at 8.72 GB; folds 201 / 584 /
> 593 s at ≤ 6.6 GB; total 1,650 s. It reproduced the chain runs bit for bit.
> - `work/logs/dayof_efp.log` is that run's printout, made under the first version of the rule.
> - `work/logs/dayof_efp_final.log` is the same outputs re-read with `--skip-runs` under the final rule:
>   `RECOMMENDATION: H1 x2, CD4 x1`, dump recipe `VCC_H1_UNION_G=0.56` with no CD4 override.
>
> **Fold operators:** default `--fold-ops prod` uses the day's 3-donor operator on every PLOLO fold, so it needs no
> extra fits. `--fold-ops held --held-ops <dir>` uses three held-out fits (`build_union_whitener.py --donors <the other
> two> --basis-prune essential`, ~5 min / 13 GB each, serial). The shortcut's measured error is in 438 §6.

## 6. The day-of shortcut: the 3-donor operator on the PLOLO folds

`--fold-ops prod` puts the day's own 3-donor operator on every PLOLO fold. It is leaky (the held library is in its basis)
but needs no extra operator fit. It was run on `efp` and `abc`: same folds, designs and level factors
(`prodcmp438.py`, `work/prodcmp438.md`). Levelled k = 1:

| panel | step | held-out fold operators | 3-donor operator on the folds | d | n+ held / prod |
|---|---|---|---|---:|---|
| efp | H1 x2 half | +0.00255 | +0.00267 | +0.00013 | 3/3 / 3/3 |
| efp | CD4 x1.5 half | -0.00096 | -0.00105 | -0.00009 | 0/3 / 1/3 |
| efp | stack | +0.00143 | +0.00153 | +0.00010 | 2/3 / 2/3 |
| efp | CD4 x1.5 on top of H1 x2 | -0.00111 | -0.00114 | -0.00003 | 0/3 / 1/3 |
| efp | H1 x2 → x3 (CD4 x1) | +0.00124 | +0.00137 | +0.00013 | 2/3 / 2/3 |
| efp | q 1.25 → 1.5 on the stack | +0.00412 | +0.00533 | +0.00121 | 3/3 / 3/3 |
| efp | key K3 → LODO (0.5, 3) on the stack | +0.00306 | +0.00191 | -0.00115 | 3/3 / 3/3 |
| abc | H1 x2 half | +0.00238 | +0.00232 | -0.00006 | 3/3 / 3/3 |
| abc | CD4 x1.5 half | +0.00183 | +0.00183 | -0.00000 | 3/3 / 3/3 |
| abc | stack | +0.00406 | +0.00398 | -0.00008 | 3/3 / 3/3 |
| abc | CD4 x1.5 on top of H1 x2 | +0.00170 | +0.00169 | -0.00001 | 3/3 / 3/3 |
| abc | q 1.25 → 1.5 on the stack | +0.00847 | +0.00888 | +0.00041 | 3/3 / 3/3 |
| abc | key K3 → LODO (0.5, 3) on the stack | +0.00274 | +0.00293 | +0.00018 | 3/3 / 3/3 |

- **Composition steps are unchanged by the shortcut to ≤ 0.00013**, well inside the rule's 0.0010 margin.
- **The rule's choice (§5's final rule) is identical under both operator sets, at q 1.25 and 1.5:**
  - `efp`: H1 ×2, CD4 ×1;
  - `abc`: H1 ×2.5, CD4 ×1.5.
  - A first version of the rule, anchored on A/B/C's weights, split on `efp` under the shortcut (×2.5 with CD4 ×1 vs
    ×1.25). That is one reason the rule is anchored on the lowest-regret point instead.
- **The key and q steps move by ±0.0012 on `efp`.** The key is fitted with the operator, so a leaky operator leaks into
  the LODO key. Neither is decided by the day-of rule: q follows the board, and the LODO key is not wired.
- **So `prod` is the day-of default.** Use `held` (three extra fits, ~15 min, 13 GB each) only if the day's printout
  sits within 0.0005 of the margin.

## 7. Limits

- **Synthetic panels.** `efp`/`efr` are random draws from donor-coverage pools minus the K562-essential list (238).
  - Arc selects from an H1 pilot screen. Arc's A/B/C carries 77 H1-union rows; the draws carry 28–35.
  - Arc's panel also carries none of the four essential screens' targets (436); the draws carry 9.
  - The Arc-like columns remove those 9. The H1 conclusions are the least certain part of this file.
  - The CD4 conclusion holds on all three draws, under every exclusion and sensitivity tried. It rests on the draws'
    CD4 residual transferring worse per row than A/B/C's does (0.52–0.79× the base vs 0.89×).
  - If Arc's selection (public-data coverage, high basal expression) gives D/E/F A/B/C's CD4 quality, ×1.5 could be
    right again. The day-of rule re-measures that on the day, and on A/B/C's panel it does move back to ×1.5.
- **The controls are A/B/C's.** The synthetic D/E/F controls are A/B/C relabelled.
  - Frames, claims and the energy anchor therefore transfer trivially here.
  - On Oct 22 the frames and claims come from the day's cells, and the anchor's A/B/C measured/emulated ratio is
    assumed to transfer to new lines. The levelled route does not depend on it.
- **Truth-side constants.**
  - Kc and sumD are A/B/C's. Kc cancels; sumD scales every mse difference by the same factor. 238 measured similar
    truth norms across panels.
  - X is scaled by PLOLO retention. ret = 1 is reported and changes no conclusion.
- **Transfer constants are A/B/C's.**
  - k = 1 / 1.407 for composition and the key.
  - 417a's calibrated k_sh for q.
  - 418's history constants were fitted on one board step per family.
- **`efr` is priced on the fold-route energies,** error bar ≤ 0.0011 levelled (§2.5), with fold-route re-level factors.
- **The energy route over-states the stack's energy step by ~0.6 %** (§2.4). This matters only for unlevelled prices
  near a cliff (`efp` D).
- **The H1 dose is extrapolated beyond ×1.**
  - There is no board point for H1 above the shipped ×1 (418's calibration). 417e's ×2 is on the board on 09-24.
  - The day-of rule caps H1 at ×3.
- **The builder's bank statistics are cached in `work/bankstats/`.** They are panel-independent. Recomputing them streams
  all four banks in one process, 11 GB file-backed (§ Status). Keep the cache for Oct 22.

## 8. Files and reproduction

```
cd experiments/438-october-stack-repricing          # PYTHONDONTWRITEBYTECODE=1; everything writes into ./work only
bash chain438b.sh full:<panel> comp:<panel>:<fold>:held:work/<panel>/level438.json ...   # serial, memory-gated
#   full438.py <panel> [--tag T --op OP] [--resume]   the full panel: KA-COMP/KEY/E420, energies, rho, level factors
#   comp438.py <panel> <fold> [--fold-ops held|prod] [--level F] [--resume]   PLOLO fold: KAs (A/B/C), grid, levelled twins
bash efr438.sh                                      # efr: folds -> foldlevel438.py (fold-route factors) -> levelled twins
uv run python price438.py                           # KA-P
uv run python summary438.py --panels abc,abcdef,real,efp,efr --full efr=efr_s1 --fold-energy efr   # section 3
uv run python summary438.py ... --exclude-oct4 --out summary438_arclike | --no-ret --out summary438_noret
uv run python robust438.py --panels abc,real,efp,efr                    # section 4.1
uv run python recommend438.py <panel> [--full efr_s1 --fold-energy]     # the rule; work/recommend438-<panel>.txt
uv run python surr438.py                                                # section 2.5
uv run python ess438.py <panel> <arm> <ref> L                           # essential / top-3 decomposition
uv run python prodcmp438.py efp abc                                     # section 6 (the day-of operator shortcut)
uv run python diag438.py <panel> <fold>                                 # section 4.5 (per-row transfer of each residual)
bash chain438b.sh full:efpsvd comp:efpsvd:<fold>:held:work/efpsvd/level438.json   # section 4.5 (own-panel deflation probe)
uv run python dayof438.py --bundle B --op OP --cd4-root CD4 [--tag T] [--board-q 1.5]   # section 5 (dayof_rehearsal438.sh)
```

- **Libraries:**
  - `lib438.py`: the panel registry and `bind`, the fold builder, `Donors438`, the composition rules, the scorer, β,
    claims and `Finals`;
  - `price438.py`: `Panel`, the pricer, the bootstrap, `fold_energies`, `october_union`.
- **Outputs:**
  - `work/<panel>/{full438.json, level438.json, comp438-<fold>-<ops>.json, arrays438-<fold>-<ops>.npz, fold-*.npz}`;
  - `work/summary438*.{md,json}`, `work/robust438*.md`, `work/surr438.json`, `work/recommend438-*.txt`;
  - `work/chain438.log` (every job with its gate reading and max RSS), `work/logs/`.
- **Imported read-only:** 162 (plolo), 218, 228, 410, 169, 414, 415, 418, 420 machinery; 233/238/244/432 operators,
  CD4 banks, bundles, keys and metas; `data/donors`, `data/challenge`.
