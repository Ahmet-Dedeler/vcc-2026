# 814: why A/B/C's panel does not look like the lab's arc800 draws (instrument lane)

Worker run, 2026-10-16, 10:52Z → 11:45Z. Mandate: `MANDATE.md`. Class: **instrument**. No candidate was priced. The draw
rule proposal is in `PROPOSAL.md`.

## 0. Answer in four lines
1. **A/B/C's panel has none of the targets from the RPE1/Jurkat essential-gene screens. An arc800 draw is expected to hold
   16.9 ± 3.7** (1,000 redraws, percentile 0.000). L1–L6 hold 14–24 each. arc800 removes only the K562 essential list.
2. **Those ~17 targets per fold carry most of the direction gains and part of the key gains on the L folds.** I took them out
   of the L folds and left A/B/C's panel uncut. ABCp's rank among L1–L6 then moves from a mean percentile of 0.194 with
   14/21 bottom ranks to **0.357 with 5/21**. Direction alone goes from 0.10 (8/9 bottom) to 0.39 (2/9). On the
   rules-matched comparator LABC it moves from 0.194 (14/21) to 0.303 (5/21).
3. **This is about which targets they are, not about how strong they are.** Removing the same number of the strongest truths
   in each fold gives 0.276 (8/21). Removing the strongest predictions gives 0.241 (12/21), and removing random targets
   gives 0.207 (14/21).
4. **A second axis exists: responsiveness in the source banks.** A/B/C's targets respond in K562/HCT116/HEK293T far more
   often than CD4-bin-matched pool targets: 26 % are silent in at least one bank, against 44 %. Matching this axis on top
   of the essential-screen cut moves the key class from 0.17 to 0.24–0.33. Knockdown efficiency, public coverage (Pan H1,
   VIPerturb, CD4 guide KD) and held-line basal expression also differ, but none of them moves the gains on its own.

Every number below comes from one of the files in §6. Where several configurations were tried, the try counts are in §5.

## 1. Inputs and known answers
- **Pool reconstruction** (`feat814.py`, EC2 run `20261016T110536Z-feat814a`). It rebuilds draw801's arc800 pool from the
  pert_counts files that draw801 recorded at draw time (`801/work/draws.json`, sha `8c8e67d2…`). Each file was read only to
  exclude its targets.
  - **KA PASS:** n_excluded 4373, basal floor 13.412203349132618 (bit-equal), CD4 strength edges, pool 5125, every stratum
    count, and L1–L6 redrawn in order with seeds 800101–800106 reproduce their recorded sha256.
  - The null pool is pool0 minus the sealed L7/L8 lists, used only as an exclusion: 4,525 targets. L7 and L8 were never
    redrawn, featurised or written.
  - Redraws use seeds 814000–814999.
- **Feature table:** `work/run1/features814.csv.gz` (sha `bdee8c2f…`), one row per gene of the null pool ∪ A/B/C ∪ L1–L6
  (4,825 genes). Source shas are in `work/run1/feat814.json` (`364b4f73…`). Sources:
  - banks: K562 GWPS, Orion HCT116, Orion HEK293T, VCC-2025 H1 (`h1_2025_full`), Pan 2025 H1 (115's compact donor);
  - CD4: GWCD4i DE table and guide-KD table;
  - DepMap 17-line subset (360);
  - A/B/C control CPM;
  - membership lists: RPE1, Jurkat, VIPerturb, PerturbMatch, TeloHAEC, Feng iPSC fitness / non-fitness / targeted.
  - Per bank: coverage, basal CPM, own-gene log effect (knockdown), off-target ‖δ‖ on genes > 10 CPM, and a robust DE
    count. The count is 518's statistic with a bank-wide gene scale from a fixed 2,000-row sample, so it does not depend
    on the panel. z_src was not used: it is a within-panel rank, so its panel distribution is uniform by construction.
- **Gain arrays** (`gains814.py`, local). These are the standard-format npz that panelrank800 reads.
  - **KA PASS on every configuration:** the all-target fold values equal `800/work/panelrank800.json` (sha `95db88fd…`),
    max |diff| 0.0.
  - LABC comparator (801 follow-up 4, `801/work/ka-labc`, `810/work/labc`, `811/work/run-m811labc`): the L-fold KA holds.
- **Rule simulation** (`sim814.py`, EC2). **KA:** R0 and R1 reproduce feat814's null statistics bit-equal.
  - The first launch (`sim814a`) **failed this KA and aborted**: `read_csv`'s default float parser does not round-trip.
  - Fixed with `float_precision="round_trip"`, and the A/B/C CD4 edges are now taken as an ndarray. The KA then passed.
- **Truth strength per held line:** the `d` arrays in those npz files (‖δ_truth‖² on each fold's scored genes).

## 2. Where A/B/C sits: percentile among 1,000 arc800 redraws (R0)
`work/run1/feat814.json`; full table in `work/sim2/sim814.json` (`b2cff6d4…`). Under R0, an arc800 draw (L1–L6) has 29–44
of the 226 panel statistics outside [0.1, 0.9]. **A/B/C has 90.**

| feature (panel statistic) | A/B/C | null mean ± sd | A/B/C pct | L1–L6 pct |
|---|---|---|---|---|
| **in RPE1 ∪ Jurkat essential screens (count)** | **0** | 16.93 ± 3.70 | **0.000** | 0.22–0.97 |
| RPE1 screen / Jurkat screen | 0 / 0 | 15.4 / 16.0 | 0.000 / 0.000 | |
| **source-bank responsiveness, median** (mean over the 3 production banks of log1p robust DE count) | 1.421 | 1.048 ± 0.039 | **1.000** | 0.11–0.89 |
| same, q90 / q10 | 2.266 / 0.597 | 1.987 / 0.417 | 1.000 / 1.000 | |
| robust DE count, median: K562 / HCT116 / HEK293T | 3 / 4 / 3 | 2.0 / 2.1 / 1.7 | 0.995 / 1.000 / 1.000 | |
| silent (0 robust DE) in ≥ 1 production bank, share | 0.263 | L 0.437, pool 0.434 | | |
| **own-gene log effect in K562, median (knockdown)** | −0.746 | −0.639 ± 0.030 | **0.000** | 0.22–0.86 |
| estimated KD < 0.5, share: K562 / HCT116 | 0.040 / 0.103 | L 0.121 / 0.225 | | |
| CD4 guide-level significant KD (count) | 294 | 273.4 ± 3.7 | 1.000 | 0.00–0.83 |
| Pan 2025 H1 screen coverage (count) | 64 | 37.9 ± 5.2 | 1.000 | 0.17–0.80 |
| VIPerturb K562 coverage | 284 | 271.2 ± 4.8 | 0.997 | |
| TeloHAEC coverage | 65 | 50.5 ± 6.3 | 0.985 | |
| VCC-2025 H1 coverage | 25 | 1.9 ± 1.3 | 1.000 | structural: the pool excludes the DE-H1 surfaces' 211 targets |
| CD4 Rest n_downstream, median (the arc800 stratifier) | 4 | 3.75 | 0.625 | matched by construction |
| CD4 max over Rest/Stim8/Stim48 n_downstream, median | 14 | 7.7 | 1.000 | |
| DepMap mean Chronos, median / min Chronos, median | −0.059 / −0.272 | −0.051 / −0.252 | 0.146 / 0.053 | |
| basal CPM of the target in K562 / HCT116 / HEK293T controls, median | 35.4 / 28.1 / 36.2 | 30.4 / 29.3 / 37.0 | 0.989 / 0.270 / 0.379 | |
| share below Arc's floor (13.41 CPM) in held-line controls: K562 / HCT / HEK | 0.146 / 0.283 / 0.156 | L 0.131 / 0.228 / 0.178 | | not the difference (agrees with 708) |
| VCC-2025 H1 control CPM, q05 / share < 5 CPM | 6.6 / 0.027 | L 3.2 / 0.077 | | modest |

On the folds themselves, using the `d` arrays and the cap3 base:
- **Truth strength sits inside the L spread** (801 §14 found the same). Median ‖δ_T‖: ABCp 6.79 / 7.09 / 6.15 (K562 / HCT /
  HEK); L 6.30–6.65 / 7.03–7.31 / 6.09–6.37. Top-30 share of Σd: ABCp 0.292 / 0.310 / 0.309; L 0.304–0.363 / 0.293–0.341 /
  0.281–0.340.
- **A/B/C's base predictions point the right way about twice as often.** Median per-target cosine (base vs truth): ABCp
  0.040 / 0.028 / 0.022; L 0.014–0.019 / 0.010–0.012 / 0.008–0.011.
- **On Orion, the predicted-to-true strength rank correlation is lower on A/B/C.** Spearman ρ(a, d): ABCp 0.27 / 0.13 /
  0.15; L 0.12–0.29 / 0.21–0.30 / 0.15–0.24.
- On the L folds the essential-screen targets have stronger truths than the other targets: mean ‖δ_T‖ 9.8 vs 7.0 on K562,
  8.7 vs 7.7 on HCT116, 7.3 vs 6.8 on HEK293T.

**The 1–3 extreme features:**
1. RPE1/Jurkat essential-screen membership (0 against 17; it drives the gains).
2. Source-bank responsiveness (A/B/C targets are rarely silent in the cancer banks; 43 % of A/B/C sits in the pool's top
   responsiveness quintile).
3. Knockdown efficiency in public screens (A/B/C has a lighter weak-KD tail).

## 3. What explains which targets carry the gains
Method: the restricted panelrank (`gains814.py`). A flagged set is removed from the L folds, or thinned to A/B/C's share.
Arc's panel is never cut. Then ABCp is re-ranked among L1–L6 on the same seven arm pairs × three held lines as
panelrank800 (21 cells). This is a post-hoc preview on the folds that found the problem. It is not a validation (§4 of the
proposal has the validation).

| L folds changed by | mean pct | bottom | key (6) | direction (9) | shape (6) |
|---|---|---|---|---|---|
| nothing (panelrank800) | 0.194 | 14/21 | 0.07, 6/6 | 0.10, 8/9 | 0.50, 0/6 |
| **essential screens removed** | **0.357** | **5/21** | 0.17, 3/6 | **0.39, 2/9** | 0.50, 0/6 |
| control: same count, top truth d removed | 0.276 | 8/21 | 0.07, 6/6 | 0.33, 2/9 | 0.40 |
| control: same count, top prediction a removed | 0.241 | 12/21 | 0.07, 6/6 | 0.25, 6/9 | 0.40 |
| control: same count, random (20 repeats) | 0.207 | 14/21 | 0.07, 6/6 | 0.10, 8/9 | 0.50 |
| weak-KD removed / matched (keep 0.5) | 0.207 / 0.207 | 13 / 14 | 0.07 | 0.15 / 0.10 | |
| silent matched (keep 0.6) | 0.201 | 13/21 | 0.07 | 0.12 | |
| responsiveness quintile matched (resp5) | 0.248 | 10/21 | 0.10 | 0.18 | |
| propensity decile matched (prop10) | 0.282 | 9/21 | 0.14 | 0.18 | |
| ess removed + weak-KD matched | 0.330 | 6/21 | 0.14 | 0.36 | |
| ess removed + silent matched | 0.323 | 6/21 | 0.14 | 0.40 | |
| ess removed + kd2 matched | 0.337 | 6/21 | 0.14 | 0.37 | |
| ess removed + resp5 matched | 0.364 | 3/21 | 0.24, 2/6 | 0.39, 1/9 | 0.45 |
| ess removed + resp5 × kd2 matched | 0.364 | 4/21 | 0.21 | 0.42 | 0.43 |
| **ess removed + prop10 matched** | **0.412** | **2/21** | 0.26, 1/6 | 0.42, 1/9 | 0.55 |
| *LABC comparator:* nothing / ess / ess + resp5 / ess + prop10 | 0.194 / 0.303 / 0.357 / 0.412 | 14 / 5 / 3 / 2 | 0.07 / 0.17 / 0.33 / 0.31 | 0.10 / 0.31 / 0.34 / 0.40 | |

Per arm, the essential-screen targets carry the L-fold direction gain on K562. Mean L fold value, all targets → without
them:
- 811 c0_l25: +0.00563 → +0.00006;
- sgu1: +0.00476 → −0.00048;
- 707: +0.00223 → +0.00149.

On the Orion lines the cut changes the L means by −39 % … +8 %. 811 HCT116 is −39 %; 707 and sgu1 on HEK293T rise
by 6–8 %. The key arm (495) falls from +0.00832 to +0.00526 on K562 (ABCp
+0.00132). 705 (the z_src rider) stays at the bottom in 3/3 cells under every rule tried. Its residual is a key-class gap
that none of these features closes. 801 §14 points the same way: A/B/C's LODO key has a light tail (cap hi ≈ 6 against
8–20), which is a property of how its targets look in the SOURCE donors. The responsiveness axis is the part of that
property I can measure, and it moves key from 0.17 to 0.24–0.33.

**Mechanism (my reading).** Arc keeps targets whose cells survive about 14 days to a FACS sort, with ≥ 80 % knockdown and
400 of about 500 cells per target. Fitness genes drop out. A/B/C carries 0 of the K562 list, which arc800 also excludes,
and 0 of the RPE1/Jurkat lists, which arc800 does not.
- Those genes respond strongly and coherently across lines (ribosome, proteasome, mitochondrial programmes). Donors
  predict their direction well (per-target cosine 0.05–0.12 against the L median of 0.008–0.019), so direction changes
  and key boosts earn on them.
- A/B/C's panel has none of them. Its gains have to come from typical targets, which are mostly responsive in the cancer
  banks but individually less predictable. That is why direction and key changes over-read on L folds (×4.7 / ×3.8 vs
  ABCp), while shape changes, which re-weight within the panel, are exchangeable (×1.16).

## 4. Draw rule simulations (A/B/C's percentile under candidate v2 rules)
1,000 redraws each, seeds 814000+. Files: `work/sim2/sim814.json` (R0, R1, V2r, V2rk, V2p, plus fresh-pool and sequential
runs) and `work/sim3/sim814.json` (`89387e1b…`, V2pc). The eight PRIMARY statistics were named in `sim814.py`'s docstring
before the first run.

| rule | pool | ess | resp med | resp q90 | K562 KD | HCT KD | Pan H1 | CD4 guide KD | VIPerturb | A/B/C stats outside [0.1, 0.9] of 226 |
|---|---|---|---|---|---|---|---|---|---|---|
| R0 arc800 | 4525 | 0.000 | 1.000 | 1.000 | 0.000 | 0.071 | 1.000 | 1.000 | 0.997 | 90 |
| R1 = R0 − ess screens | 4355 | 0.500 | 1.000 | 1.000 | 0.000 | 0.060 | 1.000 | 1.000 | 0.996 | 95 |
| V2r = R1 + resp5 weights | 4355 | 0.500 | 0.888 | 0.995 | 0.000 | 0.190 | 0.999 | 1.000 | 0.997 | 73 |
| V2rk = R1 + resp5 × kd2 | 4355 | 0.500 | 0.904 | 0.998 | 0.014 | 0.867 | 1.000 | 1.000 | 0.997 | 64 |
| **V2p = R1 + prop10 weights** | 4355 | 0.500 | 0.936 | 0.915 | 0.051 | 0.989 | **0.372** | **0.526** | **0.626** | 71 |
| V2pc = R1 + continuous odds (capped) | 4355 | 0.500 | 0.982 | 0.688 | 0.064 | 0.974 | 0.292 | 0.633 | 0.663 | 78 |

- **prop10:** the pool decile of an L2 logistic score for A/B/C vs the pool (`flags814.py`, `work/flags/prop10.json`
  `1ec590e2…`). Six features, frozen before the fit; out-of-fold AUC 0.739. Coefficients on standardised features:
  responsiveness +0.50, weak-KD −0.29, VIPerturb +0.12, Pan H1 +0.20, CD4 guide KD +0.87, H1 control CPM +0.11.
  A/B/C's decile counts are 0, 5, 12, 12, 23, 27, 26, 31, 60, 104.
- **V2p** puts 3 more of the 8 primary statistics inside [0.1, 0.9]: Pan H1, CD4 guide KD and VIPerturb.
  - Responsiveness is still borderline (0.94 / 0.92).
  - The two KD statistics fall outside in OPPOSITE directions: A/B/C has a stronger K562 knockdown and a weaker HCT116
    knockdown than v2 draws. Arc selected on its own sgRNAs in its own six lines, and public screens use other guides, so
    no public KD proxy captures that selection.
  - Under V2p the L panels have 69–118 statistics outside [0.1, 0.9]. A/B/C's 71 is at the low end of that range but above
    the 29–44 a typical draw has under its own rule. **v2 does not make A/B/C fully typical on features; it removes the
    differences that move the gains (§3).**
- **Capacity.**
  - The fresh pool (null pool minus L1–L6) cannot fill C3|Q4 (7 members for a quota of 47) or C3|Q5 (33 for 50).
  - Six disjoint v2 panels are infeasible under the weighted rules. In the sequential simulation the responsiveness
    percentile falls to 0.00 by panel 3, and V6 is infeasible.
  - So the validation must use six INDEPENDENT v2 draws. Under V2p they overlap each other by 52.6 targets on average
    (36–64 of 300) and overlap L1–L6 by 180–187 targets. Each independent V2p draw is typical on the primary statistics
    within draw noise (`indep6` in `sim2/sim814.json`; panel sha256s recorded, no lists written).

## 5. Tries, cost, rules kept
- **Tries.** 15 restricted-panelrank configurations on ABCp: 12 candidate flags or matches and 3 controls (table §3). Plus
  4 LABC re-reads of 4 of them. Five draw rules simulated: R1, V2r, V2rk, V2p, V2pc. R0 is the KA. The best preview
  (ess + prop10, 0.412) is the best of 12, so its preview number carries the best-of-N bias. The essential-screen cut
  alone (0.357) was the first try, named from the file comparison before any gain array was read.
- **Cloud: $0.041 in total.** Every run was on-demand and terminated, and `labcloud.sh status` shows no running lab
  instance (11:37Z).
  - `feat814a` (r7i.xlarge, $0.016)
  - `sim814a` (KA FAIL, aborted, $0.006)
  - `sim814b` (KA PASS; the fresh-pool step stopped on short cells, by design, $0.003)
  - `sim814c` (seq6 crashed on an exhausted weighted cell, $0.006)
  - `sim814d` (DONE 0, $0.006)
  - `sim814e` (DONE 0, $0.004)
- **Local work:** file reads and sub-second aggregations on arrays under 50 MB (gains814, flags814, small tables, an n ≤ 20
  sim check). There was no gate602 job, no local job above 100 MB, and nothing under data/ was written. The coordinator's
  WINDOW hold was in force for the whole run, so every computation that streams banks or context files, and every
  1,000-draw simulation, ran on EC2.
- I read the sealed L7/L8 and production's R6/R7 panel files only to exclude their targets (draw801's rule). No list from
  them was written or printed. The sealed lists were never redrawn or featurised.

## 6. Files (sha256)
| file | sha256 |
|---|---|
| `feat814.py` | `0984ae4eee3a22eb956d3cdfe0beed8d7aa74002a40edef339739eb4223e4767` |
| `gains814.py` | `0cf99a9598112f4e6d9d3662a773faf0b54a6f1de91108e6612a4a58bc09f532` |
| `flags814.py` | `2828c3e3f5253377c066bb6137dae5af68a18d7331f15a8ed45af48b880644b8` |
| `sim814.py` | `ceaf96990cf67a4f1dd569bc13ee746333b7b48b5e592c80e4424866c112985a` (reproduces sim2's V2p indep6 panel shas; sim2 itself ran `b5968122…`, before V2pc was appended) |
| `work/run1/features814.csv.gz` | `bdee8c2fedfab527c774bcc7648ff597706ec3f03af2c68986697fefe2ba8cd5` |
| `work/run1/feat814.json` | `364b4f73d99ae4f099d7343e78f107053d033de20b81c95cb736e38f9497a666` |
| `work/run1/null814-R0.npz` / `-R1.npz` | `3ded9bb3…` / `4dc437a9…` |
| `work/flags/essscreen.txt` | `83ae10bf5249cdd504aae63db5fa53fc0a8af7794508d8039edb1b2c0016a869` |
| `work/flags/prop10.tsv` / `prop10.json` | `66412560b0a68cf8e7ebaab0618fef7f310c116d017486a0faaacf6a575425b3` / `1ec590e21ab3e54929411016d962b2398b47647decaf2fb9147317e3290041de` |
| `work/flags/resp5.tsv`, `resp5kd.tsv`, `kd2.tsv`, `propc.tsv` | `050dd8e7…`, `afdbb120…`, `464f7b6a…`, `b0b90368…` |
| `work/sim2/sim814.json` / `work/sim3/sim814.json` | `b2cff6d4…` / `89387e1b…` |
| `work/gains814-*.json` | one per configuration in §3; each records its flag, keep share, exclusion, match, comparator and KA |
| inputs: `data/donors/rpe1/perts.txt` / `nadig_jurkat/perts.txt` | `4377d2dc…` / `f5a5c1bb…` |
| inputs: `801/work/draws.json`, `801/draw801.py`, `800/work/panelrank800.json` | `8c8e67d2…`, `c4b8507e…`, `95db88fd…` |

## 7. What is left
- The validation in `PROPOSAL.md` §4 (judge first): six V2p draws, then 801's fold builds, the seven known arms, and
  panelrank800 against the bar.
- A prediction to check on 10-22, costing nothing: **D/E/F's panel will hold ≤ 3 targets of the RPE1 ∪ Jurkat essential
  screens** (A/B/C 0; arc800 expectation 17). Its responsiveness median (the §2 statistic) should sit above the arc800
  null's 95th percentile. If both fail, Arc drew D/E/F differently and this lane's rule does not apply.
- The sealed LZ folds (L7, L8) were drawn by arc800, so they carry the essential-screen share. That can be checked only at
  their one access, by counting their essential-screen targets before the price is read.
- 705's key-rider gap is unexplained by every feature tried here.

## 8. Addendum 2026-10-16 ~11:55Z (after VERDICT-1, REVISE; text only)
- `PROPOSAL-2.md` replaces `PROPOSAL.md`. Part A is presented as FINDINGS §8.1g's recorded rule, which arc800 broke.
- **HepG2 gap closed from the source files** (obs categories only). Replogle RPE1 and Nadig HepG2 / Jurkat share one
  essential library of 2,393 targets. HepG2 has CDCA8 and Jurkat has INCENP beyond the shared set; the union is 2,394
  (`work/flags/esslib.txt`, `24fffbfe…`).
- The donor `perts.txt` lists used in this RESULTS are filtered to ≥ 20 cells and miss 43 library targets. 7 of those are
  in the null pool (CDK7, CEBPG, CMTR2, POLR3A, POT1, RBM22, ZNF830), so the corrected pool is 4,348 (`5eb2568e…`), not
  4,355.
- The numbers in §2–§4 used the donor lists. Those 7 genes are ≤ 0.2 % of the pool, so no conclusion here changes.

## 9. The frozen 10-22 reader: `predict814.py` (2026-10-17 10:25Z; the lead's request; VERDICT-2 item 4)
- **`predict814.py` sha256 `3c7e2658a32e2a8d67e90d17ededd191c3afe0042477eaa9565a70abfbd7660c`.** It does not change after
  this commit. On 10-22 the lead runs `uv run python experiments/814-panel-rule/predict814.py <D/E/F pert_counts.csv>` once.
  It writes a json and prints both verdicts.
- **P-A:** the count of panel targets in `work/flags/esslib.txt` (`24fffbfe…`). SUPPORT ≤ 3, PARTIAL 4–7, FAIL ≥ 8.
- **P-B:** the count of panel targets in prop10's pool deciles 9–10, over a denominator held at 300. SUPPORT ≥ 130,
  PARTIAL 96–129, FAIL ≤ 95.
  - **Score.** A target's score is the 10-16 fit's own decision value where it exists (`work/flags/propc.tsv`,
    `b0b90368…`, 4,825 genes). Otherwise it is the frozen formula on `work/predict/feat6_allgenes.csv.gz` (`4d6a3320…`,
    20,111 genes).
  - **Constants.** The coefficients, intercept, means and sds are written into the script as literals. The nine decile
    edges are literals too: −5.1669, −3.9834, −3.4409, −3.1340, −2.9277, −2.7474, −2.5519, −2.3318, −1.9699 (full
    precision in the script).
  - **Missing-feature rule (VERDICT-2).** A target counts as outside deciles 9–10 if it is absent from the table, or if
    its responsiveness or its H1 control expression is missing.
- **Same-panel branch:** if the panel's target set equals A/B/C's, both verdicts are VOID. The counts are still reported.
- **Built by `build814p.py`** (`d9e20244…`):
  - `--constants`, local, under a second.
    - Recovers the full-precision coefficients from propc's scores; max error 6.2e-15, and they round to prop10.json's.
    - Cuts the edges over the 4,355-gene fit pool (`e3b6b04c…`).
    - KA: all 4,825 decile labels equal `prop10.tsv`, and A/B/C has 164 in deciles 9–10.
    - One gene, OSGEPL1, sits exactly on the median edge. A rebuilt score would flip it, which is why the reader uses
      propc's frozen scores wherever they exist. Output: `work/predict/constants814.json` (`cc1955f0…`).
  - `--table`, EC2 run `20261017T102050Z-tab814p`, r7i.xlarge, $0.014, DONE 0, terminated.
    - KA PASS: the six features are bit-equal to `features814.csv.gz` on its 4,825 genes.
- **KA PASS** (`predict814.py --ka`): A/B/C gives P-A 0 and P-B 164/300, with deciles 0, 5, 12, 12, 23, 27, 26, 31, 60,
  104. Both verdicts come out VOID, as the same-panel rule requires.
- **Reader checks on existing panels.** These are not predictions; they show the reader separates arc800-style panels:
  - L1 (arc800): P-A 24 → FAIL; P-B 87 → FAIL.
  - 238's efp: P-A 9 → FAIL; P-B 48 → FAIL.
