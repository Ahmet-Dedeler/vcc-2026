# 814 PROPOSAL: draw rule v2 ("arc814") for the lab's held-out panels

Written 2026-10-16 ~11:45Z by the 814 worker. Class: **instrument**. It prices no candidate. The judge reviews this before
any v2 panel is drawn or any fold is built. The evidence is in `RESULTS.md` (sections cited as R§).

## 1. Why a v2 is needed
- On the same held lines, the same arms and the same recipe, A/B/C's real panel sits at the bottom of the arc800 draws for
  key and direction changes: 14/21 bottom ranks, mean percentile 0.19 (panelrank800).
- LABC (A/B/C's panel on the L-fold recipe) is just as low (0.194, 14/21; R§3). So the cause is the targets, not the fold
  rules.
- D/E/F's panel is drawn by Arc. A fold set that over-reads key changes ×3.8 and direction changes ×4.7 misprices every
  candidate of those classes for the final.

## 2. What A/B/C's panel has that arc800 draws do not (R§2, R§3)
| feature | A/B/C vs arc800 null | moves the gains? |
|---|---|---|
| targets of the RPE1/Jurkat essential-gene screens | 0 vs 16.9 ± 3.7 (pct 0.000) | **yes**: direction 0.10 → 0.39, all cells 0.194 → 0.357 (5/21 bottom). Same-count controls do less: top truth 0.276, random 0.207 |
| responsiveness in the source cancer banks (robust DE count) | median 1.42 vs 1.05 ± 0.04 (pct 1.000) | **partly, after the essential cut**: key 0.17 → 0.24 (resp5) / 0.26 (prop10) |
| knockdown efficiency in public screens | weak-KD share K562 4 % vs 12 %; median own-gene lfc pct 0.000 | no (alone or after the cut) |
| public coverage (Pan H1, VIPerturb, CD4 guide KD, TeloHAEC) | pct 0.985–1.000 | not tested alone; inside prop10 |

## 3. The rule
"arc814" = draw801's arc800 rule with changes A and B only:
- **Pool (A).** arc800's pool, built exactly as draw801 builds it: the official axis ∩ the three cancer donors, minus the
  K562 essential list, minus every pert_counts.csv on disk (L1–L8, R6, R7, efp/efr/real, …) as recorded at draw time,
  minus targets below A/B/C's 5th-percentile basal CPM.
  - **New in v2:** minus every target of `data/donors/rpe1/perts.txt` ∪ `data/donors/nadig_jurkat/perts.txt` (sha
    `4377d2dc…`, `f5a5c1bb…`).
  - Mechanism: Arc's cells are FACS-sorted about 14 days after the guide, with ≥ 80 % knockdown and 400 of about 500 cells
    per target. Fitness genes do not survive that selection. arc800 removed only K562's list.
- **Strata.** Unchanged: coverage class × CD4 Rest strength bin, with A/B/C's quotas after draw801's moves (the L panels'
  cell counts).
- **Within-cell weights (B).** A target is drawn with probability proportional to w(decile of s). s is a frozen L2
  logistic score for A/B/C vs the pool (`flags814.py`; coefficients, means and sds in `work/flags/prop10.json`
  `1ec590e2…`; deciles in `prop10.tsv` `66412560…`).
  - The six features: source-bank responsiveness, weak knockdown, VIPerturb-K562 coverage, Pan 2025 H1 coverage, CD4
    guide-level KD, and VCC-2025 H1 control expression.
  - w = A/B/C share / pool share per decile: 0.00, 0.17, 0.40, 0.40, 0.77, 0.90, 0.87, 1.03, 2.00, 3.46.
  - Mechanism: Arc binned strength on an H1 pilot and sampled for diverse response modes. It also required clean
    knockdown and "adequate representation in public datasets". The CD4 Rest bins do not carry responsiveness in
    proliferating lines (A/B/C: 43 % in the pool's top responsiveness quintile, 7 % in the bottom).
- **Draw.** draw801's per-cell procedure (`rng.choice` over the sorted cell members), with `p = w / Σw` in each cell, for
  six INDEPENDENT panels with seeds 814101–814106. Each panel comes from the full v2 pool. Panels overlap each other,
  expected 52.6 of 300 (36–64), and overlap L1–L6 by about 183.
  - Disjoint v2 panels are infeasible. The pool left after L1–L8 cannot fill C3|Q4 (7 members for a quota of 47).
    Sequential weighted draws exhaust the responsive targets by panel 3 (R§4).
  - **KA for the draw script:** the six target lists must reproduce the sha256 that `sim814.py --indep6` recorded for V2p
    (`work/sim2/sim814.json` → `rules.V2p.indep6.panels`):
    - V1 `546e658f…`
    - V2 `7743ca9e…`
    - V3 `7ba4fa00…`
    - V4 `6659c5a1…`
    - V5 `ba540186…`
    - V6 `5939c5bc…`
  - **Draw-time asserts:** 0 K562 essential; 0 essential-screen; 0 overlap with any excluded panel, L7/L8 included;
    every cell's count equals its quota.

## 4. Simulation: A/B/C under v2 (R§4; 1,000 draws, seeds 814000+)
- **Primary statistics.** These were named before the run. A/B/C's percentile under arc800 → v2:
  - essential screens 0.000 → 0.500;
  - responsiveness median 1.000 → 0.936, q90 1.000 → 0.915;
  - Pan H1 coverage 1.000 → 0.372;
  - CD4 guide KD 1.000 → 0.526;
  - VIPerturb 0.997 → 0.626;
  - K562 own-gene KD 0.000 → 0.051, HCT116 own-gene KD 0.071 → 0.989.
- **Against the 0.1–0.9 target: 4 of 8 inside, 2 borderline (0.92 / 0.94), 2 outside.** The two outside are the KD pair,
  in opposite directions: A/B/C's knockdown is stronger in K562 and weaker in HCT116 than in v2 draws. No public KD proxy
  reproduces Arc's own-sgRNA KD selection. I count this as a known miss, not a fix to chase. KD alone does not move the
  gains (R§3).
- **All 226 statistics:** A/B/C has 71 outside [0.1, 0.9] under v2, against 90 under arc800. A typical draw has 29–44
  under its own rule. v2 does not make A/B/C typical on every feature. What it removes are the differences that move the
  gains.
- **Preview** (not a validation). Thinning the L folds to v2's composition puts ABCp at 0.412 with 2/21 bottom. On LABC it
  is 0.412 with 2/21 (key 0.31, direction 0.40, shape 0.52). That is the best of 12 preview configurations, so it is
  biased upward. The essential-screen cut alone reads 0.357 / 5 (LABC 0.303 / 5).

## 5. Validation plan (the bar is the lead's, written in the mandate before any v2 fold exists)
1. **Draw.** One script (801 or this lane, after GO) draws V1–V6 by §3. It asserts the KA shas, writes bundles the way
   draw801 does (pert_counts.csv, A/B/C's gene_names, relabelled control clones), and records every sha at draw time.
   Local, seconds.
2. **Folds (801's recipe, unchanged).** Per panel: CD4 banks (`cd4_801.py`), the day-recipe operator fit per held line
   (`build_union_whitener.py --bundle <panel> --basis-prune essential`), and dumps for K562 GWPS, Orion HCT116 and Orion
   HEK293T. That is 18 folds, built on EC2 the way 801's run g built them ($1.98 for 18). Use the same engine KAs: per-donor
   mean == LODO key, and PT0 ×2 = 0.
3. **Arms (unchanged scripts, same arm sets as panelrank800).** ka804's L arm set (495 cap6 / cap3, 479p2g, 479sp2, 705,
   707), 810's sgu1 and 811's c0_l25, all on the 18 v2 folds. About $1.5 at 804/810/811's per-run costs.
4. **Read.** panelrank800 with the v2 folds in place of L1–L6 (one root list changed). ABCp and LABC are both reported,
   and overread800 is re-run against v2.
   - **Bar (the lead's): A/B/C's mean percentile ≥ 0.35 and bottom count ≤ 6/21, on ABCp.** LABC is reported beside it and
     must not contradict it (LABC < 0.30 or > 9 bottom → PARTIAL).
   - Also report the key / direction / shape class percentiles and the new over-read factors with their CIs.
5. **Decision.**
   - **PASS:** v2 becomes the lab's WORKING draw rule for key- and direction-class pricing. BULK-W reads of those classes
     move to the v2 folds, and the ×3.8 / ×4.7 over-read factors are replaced by v2's.
   - **FAIL:** arc800 stays. Key and direction candidates keep ABCp-based calibration. A second rule (for example the
     essential cut alone, whose preview reads 0.357 / 5) needs the judge and counts as try 2 of this family.
6. **Wear.** The v2 folds share about 60 % of their targets with L1–L6, so part of L's wear carries over. Until the bar is
   read, no new candidate is priced on v2 folds. The seven validation arms are known answers already priced on L.

## 6. Cost, class, tries
- **Cost:** ≤ $5 (folds about $2, arms about $1.5, one failure margin). Lane spend so far: $0.041 of the $5 lane budget.
  Building needs the judge's GO and fits the lane budget, or it gets its own line.
- **Class:** instrument. It changes how the lab draws held-out panels and prices nothing. D/E/F's day-of path is not
  touched: D/E/F's own panel is Arc's.
- **Tries:** the validation is try 1 of the v2-validation family. The design used 5 simulated rules and 15 + 4 preview
  configurations (R§5); they are recorded so the validation read is not mistaken for a first look.
- **Free prediction for 10-22:** D/E/F's panel holds ≤ 3 RPE1 ∪ Jurkat essential-screen targets, and its responsiveness
  median is above the arc800 null's 95th percentile. If it fails, Arc drew D/E/F differently, and v2's premise has to be
  revisited before any v2 read is used for the final.
