# 480 — Our nmae gap is information-limited. H1 is a calibrated instrument for DE-magnitude levers. The cancer folds disagree only through their power and a same-study donor. No class-based magnitude rule beats the MDE

Status: **COMPLETE** (worker, Opus 5.5). Written 2026-09-26 21:14–21:56Z after the 21:20Z redirect, then resumed at
21:57Z (memory-rule extension only) for an independent audit and three additions (§1, §3.1, §5.1, §5.2; finished 22:12Z). Analysis
only: no builds, no vcc, no cell-eval2 (the CLI or the library), no git. Every process peaked at ≤ 3.0 GB RSS; the
resume-2 scripts peaked at ≤ 0.41 GB (memory rule: ≤ 4 GB per process until 04:30Z; `/usr/bin/time -l` on each
step). Nothing was written outside `experiments/480-nmae-instrument/`.

## 0. Answer

1. **H-instrument has no support, and the premise was a misattribution by the lead.** The "K562 response means, scale1"
   description is SQ W's 09-07 entry (`SourceResponse-GammaPoisson-v2`, nmae **0.9818**), not the pack's later
   0.878–0.880 (undescribed). `BOARD-INTEL-2026-09-26.md` §5 attached it to the wrong numbers. The one disclosed pair
   agrees with H1: the board has pack-style K562 at scale 1 worse than ours by **+0.044 raw**, and H1 says
   **+0.018 … +0.070** (§2, §4.4). The range depends on the donor-to-lfc conversion and on whether our base or the
   shipped MDE arm stands in for our 09-07 construction.
2. **Second instrument (cancer folds, §4): the ranking of the donor constructions is identical on all six surfaces.**
   Two folds seemed to invert the board pair (HCT116 +0.006, K562 +0.024). Both inversions have mechanisms, and
   neither is the cell line:
   - **Power.** The folds have 150 cells per target and ~12k UMI. Thinning h1-test to that power reproduces the flip
     (−0.010 → +0.002). It raises the L1-optimal donor scale from 0.82 to 1.41, and the gain from adding donor
     magnitude to ours from −0.003 to −0.022 (dose-response in §4.3).
   - **A same-study donor.** On Orion HCT116, whose power sits in A/B/C's band, the cross-study donor alone (K562)
     gives **−0.0103 [−0.0175, −0.0041]** and k 0.77. That matches H1's −0.0102 and k 0.82 almost exactly. The
     inversion comes only from the same-study sibling (Orion HEK293T: +0.0052, k 1.48), which shares the Orion
     technical blob.
3. **h1-test has A/B/C's nmae difficulty, measured with the scorer's own anchor (§4.4).** Its split-half replicate
   nmae is **0.383**; 323's package-built anchor is 0.379, and A/B/C's are **0.369–0.431**. The K562 fold (0.517) and
   h1-test at fold power (0.517) lie outside that band. A/B/C also score 70–87 % of targets (≥ 10 gate genes), against
   h1-test's 68 %. So **434's closure stands on two instruments** once the second is read under board-like conditions
   (cross-study donors, A/B/C power).
4. **Where the gap lives (§5).** From our 0.9427 on H1:
   - oracle magnitudes at our signs recover −0.0985. 74 % of that is on genes we already call, 41 % on our top 200, and
     45 % on the 17 small-gate targets;
   - oracle signs at our magnitudes recover −0.254. 65 % is on uncalled genes (our accuracy there 0.54), 71 % on
     low/mid-expressed genes;
   - cis neighbours hold 0.9 % of the gate mass within 100 kb;
   - design strength does not separate the pots.

   Each +0.01 of gate-wide sign accuracy (iid) is worth about −0.006 … −0.008 raw nmae (§5.2), so chill's 0.835
   needs ~+0.08 accuracy. **The gap is sign information on weak responses, not magnitude.** The same split holds on
   every surface whose replicate anchor is in A/B/C's band: the sign pot is 2.6× the magnitude pot for ours on H1
   (1.8× for the donor analog), 2.0× on Orion HCT116, and all of the pot on Orion HEK293T (donor signs at chance). It
   flips only at fold power (anchor 0.517), where the gate keeps just the strong genes (§5.1).
5. **No cheap class-based magnitude rule beats the shipped MDE (§6).** Expression class, target strength, cis
   position, predicted direction, call status, donor magnitude and donor agreement, alone or crossed (up to 24
   parameters):
   - held out over H1 targets, every Δ is within ±0.003, and every interval spans 0 or is worse;
   - leave one surface out, no class structure transfers between lines;
   - the one rule that "transfers" among cancer folds (direction × expression) learns the Orion expression-direction
     artefact (plus a low-power detection asymmetry, §5.1) and costs H1 +0.029.

   **The class is closed.**
6. **The one construction a power-uncertain instrument recommends (§4.6, not a slot on its own).** Add 0.2 × the panel
   donors' same-target consensus log2 ratio to the DE channel's delivered lfc:
   - on H1 it is −0.0030 raw at A/B/C power (target-CV −0.0025 [−0.0061, +0.0009]);
   - it is −0.007 at the 250-cell / 11k-UMI powers that the replicate anchor cannot exclude;
   - it is never worse at any power in A/B/C's band.

   About 80 % of the gain sits on cells whose sign the rider keeps: it is agreement-weighted magnitude, larger where
   the donors agree with us and smaller where they disagree. Board value, if nmae prices like a pure magnitude lever
   (×1 … ×2.3): about **+0.001 … +0.005 overall**. The ~20 % from flipped signs is not priceable on H1 (calibration
   rules 3–4). It is a rider for a DE-magnitude build, not a lever to chase.
7. **The board's members cannot contradict this, and do not (§3.1, §5.2).** The nmae leaders' fid sits only +0.007 …
   +0.035 above ours (chill 0.530, cqawesome 0.538; ours 0.515), but fid is blind to gate signs at our call-set size:
   we call a median of 2,566 genes per target against n_real 89, so an iid +0.08 gate-wide sign gain moves H1 fid by
   only +0.007 (known answer: fid 0.5089291 reproduced exactly). Only three steps on the whole board moved nmae with
   pds, mse, fid and jac bit-identical. The largest are cqawesome's two (−0.0206, −0.0069). They can come from
   magnitudes of called genes or from anything about uncalled genes, and the board cannot say which.
8. **Audit (§1).** Two scripts that share no code with the analysis re-derived every headline number: 19 H1 rows and
   15 fold rows, bootstrap intervals included. All reproduce. The only fix is rounding: oracle signs give 0.68916
   (−0.2535), not 0.6891 (−0.2536).

## 1. Known answers (480's own code)

| row | expected | measured |
|---|---:|---:|
| 426 `base_s0_bps` (434's base), `ka480.py` | 0.95237 | **0.95237** |
| 434's analytic MDE row (c 0.1, cap 0.1, min count 0.3, our sign) | −0.0098 | **−0.00982** |
| emitted coarse MDE (435 `prodem_mde0.1cap0.1min0.3coarse`, the shipped DE construction on H1) | — | 0.94270 (−0.00966) |
| `gates480.py` on Orion HCT116 thinned vs 098's independent `refde_orion.npz` | identical | max \|Δq\| 0, max \|Δlfc\| 7e-12 (block sizes 4096 and 1024) |
| `gates480.py` on full h1-test vs 323's `reference.npz` | identical | same 10,779 genes; max \|Δlfc\| 2e-7; 1 of 916k cells flips at q = 0.050001 vs 0.049994; a_base nmae 0.952349 vs 0.952366 |
| zero prediction on every surface (formula check) | 1 exactly | 1.0000 on all six surfaces |
| `replicate480.py` nmae replicate anchor on h1-test vs 323's package-built anchor (`h1-anchors.json`) | 0.3791 | 0.3833 (5 splits, different seeds; split sd ≈ 0.02) |
| `decomp480.py` class contributions sum to the total | exact | asserted to 1e-9 for every class variable |
| `rules480.py` analytic shipped MDE (m 1, c 0.1) | 0.94255 (434) | 0.94254 |
| `fid480.py` fid of the shipped H1 arm, BRIEF §4 re-implemented (no cell_eval2 import) vs 435's harness value (cell_eval2's `dispatch_de_metrics`) | 0.5089291 | **0.5089291** |

**Audit (resume 2).** `audit480.py` (H1) and `audit480b.py` (folds) share no code with `common480.py` / `fold480.py`.
Each rebuilds the gate from the reference table, the nmae formula from BRIEF §7, and the donor conversion from
`scripts/build_donor_table.py`'s definition (log1p(5e4 · P / ΣP)). Then it recomputes the numbers this file quotes.

| audit rows | reproduce | notes |
|---|---:|---|
| H1: scored targets, zero/oracle checks, KA base, shipped coarse, the three §5 oracles, the §2 donor rows ((b) pc1 and naive20, (c) ×3, (d), (e)), the §2 (b) − (a) Δ, the §4.6 rider at k 0.1 and 0.2, the §5.2 gate accuracy | 19 / 19 | one rounding correction: oracle signs → 0.68916 (−0.2535), previously quoted 0.6891 (−0.2536); fixed in §5. Peak RSS 0.22 GB |
| folds (HCT116, HEK293T, K562): scored targets, mde_cons, cons@1, a_analog, a_analog − cons@1 with its 4,000-draw target-bootstrap interval | 15 / 15 | intervals identical to §4.1's to 4 decimals. Peak RSS 0.41 GB |

## 2. H1 constructions (`h1con480.py`; 58 scored targets, 38 panel-like)

Donor rows are the exact inverse of the donor table's log1p-CP50k pseudobulk: log2 of the pseudobulk-fraction ratio,
with a pseudocount in CPM (pc). `naive20` = delta/ln 2, which is a 20-CPM pseudocount. Scale 1 = pack-style.

| construction | pc0 | pc1 | pc5 | naive20 |
|---|---:|---:|---:|---:|
| (a) ours, shipped coarse MDE | 0.9427 | | | |
| (a) ours, base (426 `base_s0_bps`) | 0.9524 | | | |
| (b) K562-GWPS same target, scale 1 | 1.0123 | 1.0052 | 0.9879 | 0.9705 |
| (b) at fitted global scale k (in-sample) | 0.9686 (k 0.48) | 0.9693 (0.49) | 0.9684 (0.57) | 0.9677 (0.78) |
| (c) Orion HCT116, scale 1 | 1.1397 | 1.0487 | 1.0001 | 0.9814 |
| (c) Orion HEK293T, scale 1 | 1.0884 | 1.0171 | 0.9913 | 0.9863 |
| (c) consensus of 3 (mean of available), scale 1 | 1.0571 | 0.9632 | 0.9536 | 0.9602 |
| (e) family-balanced (K562 vs mean of the two Orion lines), scale 1 | 1.0557 | 0.9611 | 0.9510 | 0.9564 |
| (d) our sign × K562 magnitude, scale 1 | 1.0280 | 1.0208 | 1.0028 | 0.9830 |

Δ vs (a) coarse, target bootstrap (4,000 draws, unit = scored target), all 58 | 38 panel-like:

| construction | all | panel-like |
|---|---|---|
| (b) K562 scale 1, pc1 | +0.0625 [+0.0376, +0.0884] | +0.0673 [+0.0314, +0.1056] |
| (b) K562 scale 1, naive20 | +0.0278 [+0.0145, +0.0415] | +0.0363 [+0.0177, +0.0554] |
| (c) consensus scale 1, pc1 | +0.0205 [+0.0046, +0.0366] | +0.0247 [+0.0033, +0.0476] |
| (e) family-balanced, pc1 | +0.0184 [+0.0014, +0.0352] | +0.0216 [−0.0011, +0.0457] |

Every panel-donor construction is worse than ours on H1, as 434 found with a learned model.

## 3. Board calibration pair (from `data/leaderboard/*.json`, 176 snapshots, deduplicated by entry id)

| entry | date | nmae | pds | mse | description |
|---|---|---:|---:|---:|---|
| SQ W `NTC-GammaPoisson-v1` | 09-07 08:37 | 1.0098 | 0.496 | 1.029 | NTC-only null |
| **SQ W `SourceResponse-GammaPoisson-v2`** | 09-07 13:20 | **0.9818** | 0.736 | 5.498 | K562 response means, scale1 (the quoted text) |
| SQ W `sq-0908c` | 09-08 12:56 | 0.9485 | 0.774 | 0.950 | none |
| SQ W `sq-0912/13` | 09-12/13 | 0.894 / 0.890 | 0.805 | 0.806 | none |
| pack best (`hh-v24`, `M21`, `run_2026_16`) | 09-22/23 | 0.878–0.880 | 0.843–0.845 | 0.789–0.806 | none |
| **ours** "the cat sat on the modem" | 09-07 17:58 | **0.9380** | 0.856 | 0.817 | — |
| ours 446b | 09-26 | 0.8960 | 0.867 | 0.641 | — |

The misattribution (§0.1) was the lead's, in `BOARD-INTEL-2026-09-26.md` §5. The pack's 0.98 → 0.88 path after 09-07 is
undisclosed. BOARD-INTEL §5 now carries the correction (commit `16eba04c`, checked 2026-09-26 22:05Z).

### 3.1 What the board's members say about where rivals' nmae comes from (`board480.py`; context, not a test)

11,608 scored entries from all snapshots, deduplicated by entry id. Each team's latest entry, raw members:

| team | overall | nmae | fid | reach | jac | pds | mse |
|---|---:|---:|---:|---:|---:|---:|---:|
| Illumina AI | 0.3579 | 0.8397 | 0.5494 | 0.4217 | 0.0382 | 0.8898 | 0.5324 |
| cqawesome | 0.3183 | 0.8268 | 0.5379 | 0.3373 | 0.0373 | 0.8727 | 0.6060 |
| Testh | 0.3014 | 0.8503 | 0.5281 | 0.3292 | 0.0363 | 0.8607 | 0.6000 |
| CellWM | 0.2798 | 0.8677 | 0.5255 | 0.3100 | 0.0451 | 0.8597 | 0.6873 |
| nailong king | 0.2747 | 0.8512 | 0.5216 | 0.2887 | 0.0389 | 0.8641 | 0.6980 |
| chill | 0.2708 | 0.8351 | 0.5303 | 0.3076 | 0.0368 | 0.8436 | 0.7488 |
| pack (Vivian Xu / Ava Brown / HY H / SQ W) | 0.243–0.249 | 0.890–0.924 | 0.528 | 0.286–0.294 | 0.034 | 0.846–0.853 | 0.739–0.762 |
| **ours (446b)** | 0.2690 | 0.8960 | 0.5146 | 0.2793 | 0.0408 | 0.8672 | 0.6409 |

- **The nmae leaders' sign members sit close to ours.** chill and cqawesome are 0.061–0.069 better on nmae. Their fid
  is only +0.016 / +0.023 higher, reach +0.03 / +0.06, and jac is slightly lower. That does **not** mean their signs
  are no better: at our call-set size fid barely reads gate signs (§5.2: an iid +0.08 gate-wide sign gain moves H1
  fid by +0.007). Jac ≈ ours suggests call sets of our size, so the same holds for them.
- **nmae-only steps.** Only three consecutive-entry steps on the whole board kept pds, mse, fid and jac bit-identical
  while nmae moved by ≥ 0.001: cqawesome 09-18 `m0920d → m0921a` **−0.0206** (reach +0.0004), cqawesome 09-22
  `m0923s → m0924kc` **−0.0069** (reach +0.0021), and a +0.0023 elsewhere. Such a step leaves the call set and the
  signs of called genes untouched. It can change called genes' magnitudes, or anything about uncalled genes' means.
  The board cannot say which.
  - For scale: H1 says magnitude beyond the MDE is worth ≤ 0.002 raw at our signs (§6), ≤ 0.005 on the board at ×2.3.
    So cqawesome's −0.0206 is either its first MDE-like fix, or information about uncalled genes. A magnitude lever
    pays in proportion to (2a − 1): the same MDE is worth ~1.8× more at a = 0.62 than at our 0.568.
  - chill has no such step: 0.835 from its first entry.
- Testh's two DE-only steps (pds/mse bit-identical: −0.0165 on 09-20, −0.0087 on 09-25) moved fid and reach too.

## 4. Cancer folds as a second instrument, with a power control (`fold480.py`, `h1thin480.py`, `k562harness480.py`)

### 4.1 Surfaces, constructions and the calibration test

**Surfaces.** Truth gates are built with the scorer's rules (Wilcoxon vs the surface's own 18,400 controls, BH per
target over genes > 5 CPM in control, lfc = log2 ratio of arithmetic mean CPM; `gates480.py`, known answers in §1).
The held-out line's own donor table is never used on its fold.

| surface | truth cells / target | median ctrl UMI | scored targets (≥ 10 gate genes) | median gate | donors | donor-consensus gate sign acc (\|l\|-weighted) |
|---|---:|---:|---:|---:|---|---:|
| h1 (h1-test) | 400 | 20,190 | 58 / 85 | 434 | K562, HCT116, HEK293T | 0.579 |
| h1_c150 (150 random cells of each target) | 150 | 20,190 | 48 / 85 | 236 | same | 0.597 |
| h1_thin (150 cells + binomial thinning 0.6) | 150 | 12,114 | 45 / 85 | 134 | same | 0.610 |
| Orion HCT116 thinned harness | 150 | 12,174 | 127 / 132 | 70 | K562 + HEK293T (same-study sibling) | 0.583 |
| Orion HEK293T thinned harness | 150 | 14,262 | 125 / 132 | 154 | K562 + HCT116 (same-study sibling) | 0.513 |
| K562-GWPS panel cells (267 targets, 150 cells) | 150 | 11,328 | 71 / 267 | 28 | the two Orion lines (cross-study) | 0.660 |

A/B/C: 400 cells per target, median library ~20k (AGENTS §2): **h1-test has A/B/C's power; the folds do not.**

**Constructions** (identical code everywhere; sign from the available donors' consensus `cons`, pc 1 CPM):
`mde_cons` = sign(cons) × the shipped MDE law min(0.1/√cnt, 0.1) for cnt ≥ 0.3; `cons@1` = pack-style scale 1;
`cons@k` = in-sample L1-optimal global scale; `cons@kLOSO` = k fitted on the other lines; `mde+cons@k` = both;
`a_analog` = the H1-calibrated construction mde_cons + 0.575·cons (H1's in-sample k), applied unchanged; `small_base` =
0.502·cons (0.502 = median |our base| / median |cons| on H1's gate). (e) family-balanced equals `cons` on every fold
(two donors; on the K562 fold both are Xaira).

| construction (nmae) | h1 | h1_c150 | h1_thin | HCT116 | HEK293T | K562 |
|---|---:|---:|---:|---:|---:|---:|
| zero (known answer) | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| mde_cons | 0.9804 | 0.9825 | 0.9838 | 0.9844 | 0.9971 | 0.9753 |
| cons@1 (pack-style) | 0.9632 | 0.9495 | 0.9522 | 0.9356 | 0.9905 | 0.8731 |
| cons@k | 0.9598 | 0.9474 | 0.9474 | 0.9328 | 0.9863 | 0.8422 |
| cons@kLOSO | 0.9739 | 0.9474 | 0.9487 | 0.9344 | 0.9997 | 0.8832 |
| mde+cons@k | 0.9530 | 0.9381 | 0.9388 | 0.9242 | 0.9851 | 0.8334 |
| a_analog | 0.9530 | 0.9499 | 0.9537 | 0.9417 | 0.9851 | 0.8967 |
| **(a) ours, shipped coarse** | **0.9427** | **0.9403** | **0.9422** | — | — | — |
| (a) ours, base | 0.9524 | 0.9519 | 0.9541 | — | — | — |
| in-sample k of cons@k / of mde+cons@k | 0.82 / 0.575 | 1.21 / 1.03 | 1.41 / 1.20 | 1.23 / 1.10 | 0.63 / 0.575 | 1.74 / 1.60 |
| ceiling: sign(cons) × oracle \|l\| | 0.7936 | 0.7850 | 0.7653 | 0.8485 | 0.9851 | 0.6602 |

The ordering among donor constructions is the same on all six surfaces (mde+cons@k < cons@k < cons@1 < mde_cons,
single donors mostly between). What changes with the surface is the preferred donor scale.

**The calibration test** (the board's one disclosed pair: pack-style scale 1 loses to our construction by +0.044 raw).
Δ = H1-calibrated construction − pack-style scale 1; negative agrees with the board. Target bootstrap, 4,000 draws,
unit = scored target of that surface:

| surface | a_analog − cons@1 | real (a) coarse − cons@1 | agrees with the board? |
|---|---|---|---|
| h1 (A/B/C power) | −0.0102 [−0.0191, −0.0014] | **−0.0205 [−0.0366, −0.0046]** | yes |
| h1_c150 | +0.0005 [−0.0067, +0.0077] | −0.0092 [−0.0251, +0.0060] | a_analog: no (tie) |
| h1_thin (fold power) | +0.0015 [−0.0072, +0.0103] | −0.0101 [−0.0283, +0.0075] | a_analog: no (tie) |
| Orion HCT116 | +0.0061 [+0.0014, +0.0104] | — | **inverts** |
| Orion HEK293T | −0.0054 [−0.0090, −0.0020] | — | yes (donor signs near chance, 0.51) |
| K562-GWPS | +0.0236 [+0.0137, +0.0335] | — | **inverts** |

### 4.2 Adding donor magnitude to ours: the power control

**Adding the donor consensus to our shipped construction** (`a_coarse + k·cons`, k in-sample; the most favourable
reading):

| surface | k | Δ nmae, all scored | Δ, panel-like targets (no essential-screen coverage) |
|---|---:|---|---|
| h1 (A/B/C power) | 0.175 | −0.0031 [−0.0067, +0.0002] | −0.0042 [−0.0090, +0.0004] (38) |
| h1 at full power, restricted to h1_thin's 45 targets (`extras480.py aplus45`) | 0.10 | −0.0012 | |
| h1_c150 | 0.675 | −0.0196 [−0.0332, −0.0067] | −0.0256 [−0.0430, −0.0070] (29) |
| h1_thin | 0.80 | −0.0224 [−0.0415, −0.0044] | −0.0290 [−0.0516, −0.0072] (27) |

**Reading.**
- The same line, the same donors and the same predictions give a 7× larger donor-magnitude gain when only the truth's
  power drops. The restriction to the same 45 targets keeps the full-power gain at −0.0012, so this is the gate's gene
  content, not target selection. At low power only large-|l| genes pass BH; there the donor signs are better (0.579 →
  0.610) and the donors' magnitudes are the right size.
- The K562 fold and h1_thin sit below A/B/C's power (replicate anchor 0.517 against 0.369–0.431, §4.4), so they
  over-read donor-magnitude levers. The K562 fold over-reads most (cross-study, sign accuracy 0.66, median gate 28).
  The Orion folds' anchors sit at the band's upper edge. Their inversion comes from the same-study donor (§4.5).
- The MDE law's own value on a small base (0.502·cons + MDE vs 0.502·cons): h1 −0.0138, h1_c150 −0.0168, h1_thin
  −0.0157, HCT116 −0.0154, HEK293T −0.0020, K562 −0.0240. The shipped add on our real base: −0.0097 (h1), −0.0116,
  −0.0119. The board read 433a's MDE at −0.0168 (1.76× H1).

### 4.3 Power dose-response (`power480.py`; same predictions, only the truth's power changes)

| surface | scored | median gate (scored) | donor-cons acc | k (cons@k) | a_analog − cons@1 [95 %] |
|---|---:|---:|---:|---:|---|
| h1 400 cells, 20k UMI (h1-test) | 58 | 434 | 0.579 | 0.82 | −0.0102 [−0.0191, −0.0014] |
| h1 250c 20k | 52 | 385 | 0.586 | 0.92 | −0.0074 [−0.0152, +0.0001] |
| h1 400c 11k | 55 | 259 | 0.584 | 0.99 | −0.0034 [−0.0110, +0.0042] |
| h1 250c 11k | 49 | 198 | 0.596 | 1.10 | −0.0021 [−0.0091, +0.0049] |
| h1 150c 20k | 48 | 236 | 0.597 | 1.21 | +0.0005 [−0.0067, +0.0077] |
| h1 100c 20k | 45 | 144 | 0.609 | 1.44 | +0.0044 [−0.0036, +0.0126] |
| h1 150c 11k | 46 | 125 | 0.608 | 1.48 | +0.0039 [−0.0037, +0.0117] |
| K562 fold 150c 11k (267 targets) | 71 | 28 | 0.660 | 1.73 | +0.0236 [+0.0137, +0.0335] |
| K562 fold 150c, the 55 targets with ≥ 250 cells | 12 | 22 | 0.693 | 1.75 | +0.0339 [+0.0183, +0.0517] |
| K562 fold 250c, same 55 targets | 14 | 48 | 0.677 | 1.29 | +0.0113 [−0.0106, +0.0328] |
| HCT116 fold 150c 12k | 127 | 70 | 0.583 | 1.23 | +0.0061 [+0.0014, +0.0104] |
| HEK293T fold 150c 14k | 125 | 154 | 0.513 | 0.63 | −0.0054 [−0.0090, −0.0020] |

The same monotone law holds within H1 (cells and depth) and within K562 (150 → 250 cells on the same targets): less
power, smaller gates of larger-|l| genes, higher donor sign accuracy on the gate, larger optimal donor scale.

### 4.4 Which power is A/B/C's? The scorer's replicate anchor (`replicate480.py`, `power480b.py`)

`replicate480.py` reimplements cell-eval2's `lfc_nmae_ref` raw estimator. It uses disjoint halves of floor(n/2) cells
per group, controls split too, each half against its own controls, and a half's filtered or non-finite lfc read as 0.
The gate and denominator come from the full table. Five splits. **Known answer:** h1-test 0.383 vs 323's
package-built anchor 0.379 (`323/work/h1-anchors.json`; splits 0.349–0.409 there, 0.363–0.412 here; different seeds).

| surface | replicate nmae anchor r | scored / targets | in A/B/C's band 0.369–0.431? |
|---|---:|---:|---|
| **h1-test (400c 20k)** | **0.383** | 58 / 85 (68 %) | **yes** |
| h1 250c 20k | 0.381 | 52 / 85 (61 %) | yes |
| h1 400c 11k | 0.405 | 55 / 85 (65 %) | yes |
| h1 150c 20k | 0.437 | 48 / 85 (56 %) | just above |
| h1 150c 12k (h1_thin) | 0.517 | 45 / 85 (53 %) | no |
| Orion HCT116 fold | 0.427 | 127 / 132 | yes |
| Orion HEK293T fold | 0.440 | 125 / 132 | just above |
| K562 fold | 0.517 | 71 / 267 | no |
| **A/B/C (metrics BRIEF §8)** | **0.369–0.431** | **209–261 / 300 (70–87 %)** | — |

- A/B/C's anchor and scored fraction match h1-test at full power. Its scored fraction (70–87 %) is at or above h1-test's
  68 %, and every thinned level scores fewer targets.
- The anchor cannot exclude H1 at 250 cells or at 400 cells / 11k (0.381, 0.405). The board's own n_real estimates
  bracket the same range: audit 181 gives a median of 90–110 (h1-test's all-target median is 89), and audit 192 gives
  31–45 (h1 250c: 44; 400c 11k: 47).
- **The board pair is consistent with that range** (`power480b.py`, pack-style K562 at scale 1 minus ours):

| H1 level | all-target median gate | our real MDE add | K562@1 (pc1) − ours | K562@1 (naive20) − ours | cons@1 − ours | ours + k·cons (k in-sample) |
|---|---:|---|---|---|---|---|
| 400c 20k (h1-test) | 89 | −0.0097 [−0.0135, −0.0057] | +0.0625 [+0.0376, +0.0884] | +0.0278 [+0.0145, +0.0415] | +0.0205 | −0.0031 (0.18) |
| 250c 20k | 44 | −0.0108 | +0.0383 [+0.0196, +0.0588] | +0.0182 | +0.0093 | −0.0096 (0.40) |
| 400c 11k | 47 | −0.0106 | +0.0346 [+0.0163, +0.0558] | +0.0177 | +0.0089 | −0.0098 (0.43) |
| 150c 20k | 23 | −0.0116 | +0.0181 [+0.0000, +0.0351] | +0.0180 | +0.0092 | −0.0196 (0.68) |
| 150c 11k | 16 | −0.0133 | +0.0082 [−0.0109, +0.0272] | +0.0194 | +0.0057 | −0.0321 (0.88) |

  The board's +0.044 falls between h1-test (+0.0625) and the 250-cell level (+0.038) under the pc1 conversion. It is
  far from the fold-power levels (+0.008 … +0.018). The board's 1.76× on the MDE is not explained by power: the MDE add
  grows only to −0.011 … −0.013 at reduced power.

### 4.5 The Orion inversions come from the same-study donor (`extras480.py crossstudy`)

| fold | donor for cons | donor sign acc | k (cons@k) | a_analog − cons@1 [95 %] |
|---|---|---:|---:|---|
| HCT116 (anchor 0.427, in band) | K562-GWPS (cross-study) | 0.549 | 0.77 | **−0.0103 [−0.0175, −0.0041]** |
| HCT116 | Orion HEK293T (same study) | 0.588 | 1.48 | +0.0052 [+0.0013, +0.0091] |
| HEK293T (anchor 0.440) | K562-GWPS (cross-study) | 0.517 | 0.50 | −0.0129 [−0.0181, −0.0086] |
| HEK293T | Orion HCT116 (same study) | 0.515 | 0.74 | −0.0028 [−0.0060, +0.0006] |
| h1-test (reference) | three panel donors (all cross-study) | 0.579 | 0.82 | −0.0102 [−0.0191, −0.0014] |

A/B/C have no same-study donor. Under that condition, at in-band power, the cancer fold and H1 give the same number.
The Orion same-study sibling shares the Orion technical blob (§5.1): its magnitudes "predict" the fold's truth for
technical reasons.

### 4.6 The donor-consensus rider, priced at every A/B/C-compatible power (`rider480.py`)

Rider: ours (shipped coarse) + k × X. X is either `cons` (mixed: it can flip our sign where |k·cons| > |ours|) or
s·|cons| (magnitude only, our sign kept). Target-CV k: 4/5 of targets, 50 partitions. Target bootstrap, 4,000 draws.

| H1 level | kind | fixed k 0.1 / 0.2 / 0.3 / 0.4 | target-CV Δ [95 %] (mean k) | panel-like |
|---|---|---|---|---|
| 400c 20k (h1-test) | mixed | −0.0024 / **−0.0030** / −0.0016 / +0.0015 | −0.0025 [−0.0061, +0.0009] (0.18) | −0.0035 [−0.0081, +0.0011] |
| 400c 20k | magnitude-only | −0.0013 / −0.0006 / +0.0022 / +0.0069 | −0.0006 [−0.0028, +0.0016] (0.12) | −0.0019 [−0.0047, +0.0010] |
| 250c 20k | mixed | −0.0042 / **−0.0073** / −0.0089 / −0.0096 | −0.0086 [−0.0155, −0.0017] (0.40) | −0.0098 [−0.0183, −0.0007] |
| 250c 20k | magnitude-only | −0.0031 / −0.0050 / −0.0055 / −0.0048 | −0.0047 [−0.0094, +0.0002] (0.30) | −0.0066 [−0.0126, −0.0002] |
| 400c 11k | mixed | −0.0040 / **−0.0070** / −0.0089 / −0.0098 | −0.0087 [−0.0158, −0.0015] (0.42) | −0.0112 [−0.0210, −0.0006] |
| 400c 11k | magnitude-only | −0.0036 / −0.0062 / −0.0075 / −0.0075 | −0.0066 [−0.0128, −0.0007] (0.35) | −0.0106 [−0.0190, −0.0023] |
| 150c 20k (just above the band) | mixed | −0.0050 / −0.0096 / −0.0134 / −0.0162 | −0.0183 [−0.0311, −0.0051] (0.67) | −0.0241 [−0.0418, −0.0064] |

- **k 0.2 (mixed) is the robust choice.** It is at the h1-test optimum (−0.0030), and it gains −0.007 at the two
  reduced-power levels the anchor cannot exclude. It is never worse at any level. Larger k pays only if A/B/C sit at
  the low-power end.
- At h1-test's power the always-boost magnitude-only form is null (−0.0006). The mixed form's gain is
  agreement-weighted magnitude: 80 % on cells whose sign it keeps, 20 % on sign flips (class table below). The flip
  part is sign-like, and H1 has mis-signed nmae for DE-sign/carrier changes twice (`docs/H1-BOARD-CALIBRATION.md`
  rules 3–4).
- Board translation (each −0.01 raw ≈ +0.0028 overall):
  - pure-magnitude reading (rule 1, ×1 … ×2.3) for the kept-sign 80 %: −0.003 … −0.007 raw on H1 → **+0.001 …
    +0.005 overall**;
  - the flip part (≈ 20 %): unpriceable, with a board-sized downside.
- **Where the rider's gain sits** (k 0.2, mixed; class contributions to Δ, `extras480.py riderclass`):

  | class | h1-test (−0.0030) | 250c 20k (−0.0073) |
  |---|---|---|
  | cells whose sign the rider flips / keeps | −0.0006 / −0.0024 | −0.0006 / −0.0067 |
  | donors agree / disagree with our sign | −0.0019 / −0.0011 | −0.0062 / −0.0011 |
  | truth down / up | −0.0086 / +0.0056 | −0.0105 / +0.0033 |
  | we say down / up | −0.0040 / +0.0009 | −0.0061 / −0.0012 |
  | control count low / mid / high tertile | −0.0033 / +0.0002 / +0.0001 | −0.0049 / −0.0021 / −0.0002 |
  | called / uncalled | −0.0011 / −0.0020 | −0.0035 / −0.0038 |
  | cis < 100 kb / rest | −0.0003 / −0.0027 | −0.0004 / −0.0069 |

  About 80 % of the gain is on cells whose sign the rider keeps. It is agreement-weighted magnitude: the magnitude
  grows where the donors agree with us and shrinks where they disagree. That is why the always-boost magnitude-only form
  is null. It leans down (the donor consensus is down-biased like CLEAN4) and sits in low-expressed genes. This is 434's
  sign-confidence idea as one linear term; 434's GBM with the panel donors on panel-like targets did not find it
  (+0.0029).
- **Robustness** (`extras480.py riderrobust`, fixed k, no fitting), at h1-test / 250c 20k:
  - the conversion barely matters: k 0.2 with the three-donor consensus gives −0.0024 / −0.0077 (pc 0.5), −0.0030 /
    −0.0073 (pc 1), −0.0030 / −0.0059 (pc 5) and −0.0023 / −0.0042 (pc 20). At k 0.1 with pc 1 it is −0.0024
    [−0.0044, −0.0006] on h1-test;
  - **K562-GWPS carries it.** K562 alone gives −0.0036 / −0.0078 at k 0.2, pc 1. The two Orion lines alone give
    −0.0001 / −0.0044.
- This is 427e/430's "blend" measured again: 430 found the panel-donor blend's increment over the MDE to be −0.0014 on
  38 targets. The new information is only that its value grows as the reference's power falls. It is not a new lever.

## 5. Where our remaining nmae gap lives (step 2, `decomp480.py`, H1, our shipped construction, 58 targets)

nmae is additive over any partition of a target's gate cells, so each class's contribution below adds to the total
exactly (asserted in the script). Per class: **mass** = share of the gate's |l| (the zero prediction's contribution);
**ours** = our contribution; **acc** = |l|-weighted sign accuracy of our sign; **rec** = how much of ours an oracle
recovers: **mag** = oracle magnitudes at our signs (434's bound, total 0.844), **sign** = oracle signs at our
magnitudes, **abst** = abstain exactly where our sign is wrong; **best c** = the in-sample best class-constant MDE add
(shipped 0.1) and its gain.

Totals: ours 0.9427. Oracle magnitudes at our signs → 0.8442 (−0.0985). Oracle signs at our magnitudes → 0.6892
(−0.2535). Abstain where wrong → 0.8132 (−0.1295). Both → 0.

| class | mass | ours | acc | rec mag | rec sign | rec abst | best c (gain) |
|---|---:|---:|---:|---:|---:|---:|---|
| **called by us (q < 0.05)** | 0.345 | 0.303 | 0.652 | **0.0725** | 0.0877 | 0.0460 | 0.09 (−0.0000) |
| **uncalled** | 0.655 | 0.640 | 0.540 | 0.0260 | **0.1659** | 0.0835 | 0.13 (−0.0006) |
| truth up (oracle class) | 0.473 | 0.463 | **0.481** | −0.0544 | 0.1404 | 0.0712 | 0.00 (−0.0081) |
| truth down (oracle class) | 0.527 | 0.480 | 0.643 | 0.1529 | 0.1131 | 0.0583 | 0.20 (−0.0090) |
| we say up | 0.378 | 0.355 | 0.534 | 0.0278 | 0.1131 | 0.0583 | 0.13 (−0.0004) |
| we say down | 0.622 | 0.588 | 0.593 | 0.0707 | 0.1404 | 0.0712 | 0.10 (0.0000) |
| control count low (< 0.87/cell) | 0.472 | 0.447 | 0.551 | 0.0361 | 0.0857 | 0.0433 | 0.16 (−0.0012) |
| mid | 0.271 | 0.255 | 0.570 | 0.0253 | 0.0952 | 0.0483 | 0.10 (0.0000) |
| high (≥ 2.60/cell) | 0.257 | 0.241 | 0.613 | 0.0371 | 0.0726 | 0.0379 | 0.08 (0.0000) |
| cis: TSS < 10 kb of the target | 0.004 | 0.002 | 1.000 | 0.0023 | 0.0000 | 0.0000 | 0.60 (−0.0005) |
| cis: 10–100 kb | 0.005 | 0.004 | 0.651 | 0.0027 | 0.0001 | 0.0001 | 0.21 (0.0000) |
| 0.1–1 Mb | 0.006 | 0.005 | 0.599 | 0.0034 | 0.0007 | 0.0003 | 0.20 (−0.0001) |
| same chromosome > 1 Mb | 0.048 | 0.044 | 0.574 | 0.0048 | 0.0121 | 0.0062 | 0.12 (0.0000) |
| other chromosome | 0.919 | 0.870 | 0.569 | 0.0847 | 0.2364 | 0.1208 | 0.12 (−0.0002) |
| design strength (467) 0–1 | 0.224 | 0.213 | 0.554 | 0.0156 | 0.0582 | 0.0293 | 0.13 (−0.0002) |
| 2–5 | 0.448 | 0.419 | 0.548 | 0.0587 | 0.1113 | 0.0565 | 0.12 (−0.0003) |
| ≥ 6 | 0.328 | 0.311 | 0.578 | 0.0242 | 0.0841 | 0.0436 | 0.10 (0.0000) |
| gate size n_real 10–99 (17 targets; oracle) | 0.293 | 0.265 | 0.606 | **0.0442** | 0.0523 | 0.0264 | 0.15 (−0.0005) |
| 100–499 | 0.224 | 0.215 | 0.544 | 0.0121 | 0.0653 | 0.0333 | 0.09 (0.0000) |
| 500–1,999 | 0.259 | 0.248 | 0.555 | 0.0193 | 0.0741 | 0.0379 | 0.13 (−0.0003) |
| ≥ 2,000 | 0.224 | 0.215 | 0.574 | 0.0229 | 0.0618 | 0.0319 | 0.10 (0.0000) |
| our rank 0–199 | 0.081 | 0.054 | 0.860 | **0.0404** | 0.0096 | 0.0058 | 0.07 (−0.0002) |
| 200–999 | 0.098 | 0.087 | 0.661 | 0.0222 | 0.0353 | 0.0187 | 0.06 (−0.0003) |
| 1,000–2,999 | 0.190 | 0.186 | 0.579 | 0.0116 | 0.0523 | 0.0263 | 0.16 (−0.0006) |
| ≥ 3,000 | 0.631 | 0.616 | 0.539 | 0.0243 | **0.1563** | 0.0787 | 0.13 (−0.0008) |
| donor consensus agrees with our sign | 0.709 | 0.651 | 0.605 | **0.1240** | 0.1874 | 0.0962 | 0.13 (−0.0007) |
| disagrees | 0.284 | 0.284 | **0.482** | −0.0244 | 0.0646 | 0.0326 | 0.00 (−0.0009) |
| donor \|cons\| top tertile | 0.420 | 0.377 | 0.606 | 0.0868 | 0.0919 | 0.0474 | 0.12 (−0.0002) |

**Where the pots sit.**
- **Magnitude (−0.0985 at our signs).** 74 % is on genes we call (35 % of the mass), 41 % on our top-200 genes (8 % of
  the mass), 45 % on the 17 small-gate targets, and 88 % in the top tertile of donor |cons|. Our called genes carry
  median |l̂| 0.071 against |l| 0.170. This pot needs per-gene effect sizes inside our confident set. No class-constant
  rule captures it (best in-sample class gains ≤ 0.0016; §6). The panel donors carry little of it at A/B/C's power: the
  rider in §4.6 recovers −0.003 of the −0.0985.
- **Sign (−0.254 at our magnitudes).** 65 % is on uncalled genes and 71 % on low/mid-expressed genes, where our accuracy
  is 0.54 and |l| is large (median 0.34–0.47, set by the MDE). This pot needs the recipient line's direction for weak
  responses: the programme family (458d / 458g, board test 09-30), not magnitude.
- **Cis neighbours are negligible**: 1.5 % of the mass within 1 Mb, 0.4 % within 10 kb (425's finding, now with the
  100 kb window the mandate asked for: 0.9 % within 100 kb). They are already right (acc 1.0 < 10 kb).
- **Target strength does not separate the pots**: design-strength classes have sign accuracy 0.55–0.58 and the same
  best c. n_real (an oracle) separates magnitude (small gates) from sign (large gates), which is 434's finding.
- **Our down-bias shows**: on truth-up genes our sign accuracy is 0.481 (below chance), so oracle magnitudes at our
  signs would lose there (−0.054). Knowing the truth direction would let the MDE add go to 0 on up genes and 0.2 on down
  genes (−0.017 in-sample), but our predicted direction does not separate them (we say up: 0.534; down: 0.593).
- H1's truth leans down with expression (P(up) 0.51 in the lowest control-count decile → 0.41 in the highest;
  `extras480.py updown`). Our signs already lean harder (P(we say up) 0.38 → 0.26).

### 5.1 The same decomposition on the folds

The folds have no "ours", so `decomp480b.py` decomposes the a_analog construction (the H1-calibrated donor analog,
§4.1) on every surface. Totals are summarised by `pots480.py`. Pots are the nmae each oracle would recover; the split
columns are shares of the pot.

| surface (construction) | replicate anchor | nmae | magnitude pot | sign pot | abstain pot | mag pot on truth-down | mag pot in gates < 100 | sign pot in low/mid expression |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| h1 (ours, shipped; §5) | 0.383 | 0.9427 | 0.0985 | **0.2535** | 0.1295 | 155 % (up genes −55 %) | 45 % | 71 % |
| h1 (a_analog) | 0.383 | 0.9530 | 0.1600 | **0.2922** | 0.1532 | 88 % | 50 % | 70 % |
| Orion HCT116 (a_analog) | 0.427 | 0.9417 | 0.0932 | **0.1870** | 0.0936 | 78 % | 56 % | 77 % |
| Orion HEK293T (a_analog; donor signs 0.51) | 0.440 | 0.9851 | −0.0002 | **0.2573** | 0.1292 | — | — | 75 % |
| h1_thin (150 cells, 12k UMI) | 0.517 | 0.9537 | **0.1963** | 0.1655 | 0.0842 | 94 % | 60 % | 66 % |
| K562 fold (150 cells, 11k UMI) | 0.517 | 0.8967 | **0.2641** | 0.0924 | 0.0463 | 84 % | 97 % | 75 % |

- **At A/B/C-band power the pot is sign on every surface.** It is 1.8–2.6× the magnitude pot on H1 and 2.0× on HCT116.
  On HEK293T, where donor signs are at chance, it is all of the pot. The sign pot sits 66–77 % in low/mid-expressed
  genes everywhere.
- **At fold power the balance flips** (h1_thin 0.196 magnitude vs 0.166 sign; K562 0.264 vs 0.092). The gate keeps only
  large-|l| genes, where donor signs are better, so magnitude becomes the pot. This is §4's power effect in pot form. It
  is also why the best class-constant add at fold power is 3–5× H1's (c 0.3–0.6 vs 0.1–0.15).
- The magnitude pot sits on truth-down genes on every surface, because donor and CLEAN4 signs are down-biased.

**The truths' expression-direction profile** (`artefact480.py`; share of gate cells whose truth is up, by the gene's
control expression, gate-pooled):

| surface | cells × depth | lowest decile | tertiles (low / mid / high) | highest decile |
|---|---|---:|---|---:|
| h1-test | 400 × 20k | 0.51 | 0.50 / 0.49 / 0.44 | 0.41 |
| h1 250c / 400c 11k | 250 × 20k / 400 × 11k | 0.55 / 0.56 | 0.52 / 0.49 / 0.43 ; 0.51 / 0.47 / 0.43 | 0.40 / 0.40 |
| h1 150c / 100c | 150 × 20k / 100 × 20k | 0.61 / 0.68 | 0.55 / 0.49 / 0.43 ; 0.59 / 0.50 / 0.43 | 0.38 / 0.39 |
| h1_thin (150c, 12k) | 150 × 12k | 0.71 | 0.59 / 0.48 / 0.42 | 0.38 |
| K562 fold | 150 × 11k | 0.90 | 0.73 / 0.51 / 0.47 | 0.47 |
| Orion HCT116 | 150 × 12k | 0.85 | 0.70 / 0.31 / 0.21 | **0.09** |
| Orion HEK293T | 150 × 14k | 0.87 | 0.76 / 0.44 / 0.12 | **0.08** |

Two artefacts, separable by the thinned H1 gates:
- **The up-excess at low expression is partly power.** Decreases of low-expressed genes are hard to call with few
  cells, so their gate fills with increases. H1's lowest decile goes 0.51 → 0.71 as power falls to fold level; the
  K562 fold reaches 0.90.
- **The down-excess at high expression is Orion's own.** Its highest decile is 0.08–0.09 up, against 0.37–0.41 on every
  H1 power level and 0.47 on K562. This is the "technical blob" the harness caveat warns about.

An expression × direction rule learnt on the Orion folds (§6) therefore learns one artefact A/B/C lack entirely (the
Orion blob), and one they have only at H1's level (the power asymmetry).

### 5.2 What sign information is worth (`signvalue480.py`, H1, our shipped construction)

Raise our |l|-weighted gate sign accuracy in a class by +Δa: flip random wrong-signed cells to right, iid within the
class, 20 draws; draw sd ≤ 0.002. Keep our magnitudes, or re-fit the global MDE c (in-sample). This is 434f's iid
assumption, so treat it as a price for the information, not a lever.

| class (pooled \|l\| share) | our acc | +0.02 | +0.05 | +0.10 | +0.05 with c re-fit |
|---|---:|---:|---:|---:|---:|
| all gate genes (1.00) | 0.568 | −0.0126 | −0.0321 | −0.0641 | −0.0395 (c 0.17) |
| uncalled (0.75) | 0.540 | −0.0079 | −0.0196 | −0.0390 | −0.0235 (c 0.15) |
| called (0.25) | 0.652 | −0.0053 | −0.0137 | −0.0271 | −0.0159 (c 0.14) |

- Each +0.01 of gate-wide accuracy is worth about −0.0065 raw nmae, or −0.008 with magnitudes re-fitted. chill's
  0.835 against our 0.896 (board) is −0.061, roughly +0.08 of accuracy, which agrees with 434f's map.
- The table prices the analytic construction (426 base + MDE, 0.94254). The same flips on the emitted coarse lfc
  (`fid480.py`, 0.94270) give −0.0116 / −0.0294 / −0.0470 at +0.02 / +0.05 / +0.08, which is −0.0059 per +0.01. A flip
  gains 2·min(|l̂|, |l|), and the emitted arm's wrong-signed cells carry slightly smaller |l̂|.
- **fid does not see this accuracy** (`fid480.py`; known answer: our H1 fid 0.5089291, identical to 435's
  harness value). We call a median of **2,566** genes per target against a median n_real of **89**. The yield is ~1
  (0.972), so fid ≈ sign precision over all our calls (0.525), and most calls lie outside the reference's gate. The
  called share of the gate's |l| is 0.25.

  | iid gate-wide sign gain | Δnmae (emitted arm) | Δfid, flips on gate cells only | Δfid, same flip rate on our non-gate calls too |
  |---|---:|---:|---:|
  | +0.02 | −0.0116 | +0.0019 | +0.0215 |
  | +0.05 | −0.0294 | +0.0046 | +0.0537 |
  | +0.08 | −0.0470 | **+0.0073** | +0.0859 |

  The last column is an upper bound: non-gate calls have reference lfcs dominated by noise, so real information
  cannot move their signs much. A rival whose call sets are our size (its jac ≈ ours says so, §3.1) could carry +0.08
  of gate sign accuracy and show only ~+0.007 of fid. **The leaders' fid (+0.007 … +0.035 over ours) therefore neither
  confirms nor rules out 434f's reading.** Before this check, a fid close to ours could have been read as "their signs
  are no better, so their edge is magnitude". That inference is invalid.
- The information that would buy it:
  - the recipient line's direction for weak responses in low/mid-expressed genes, the uncalled 65 % of the sign pot.
    That is the programme family (458d / 458g, board test 09-30), not magnitude;
  - per-gene effect sizes inside our confident calls, for the magnitude pot. The panel donors do not carry these at
    A/B/C's power (§4.6).

## 6. A cheap class-based magnitude rule? None beats the shipped MDE (step 3)

**Rule shapes** (oracle-free classes only): c-only `base + s·c_C·min(1/√cnt, 1)`; m+c `m_C·base + s·c_C·…`. Shipped =
(m 1, c 0.1), whose analytic form reproduces 434's row (0.94254; the emitted coarse arm reads 0.94270).

**H1 (A/B/C power), fitted on 4/5 of the targets, tested on the held-out fifth** (`cv480.py`; 5-fold over targets × 50
random partitions; per-target held-out Δ averaged over partitions; target bootstrap, 4,000 draws, unit = scored target):

| class variable | shape | in-sample Δ (bound) | held-out Δ vs shipped [95 %] | panel-like (38) |
|---|---|---:|---|---|
| none (one global class) | c | −0.0003 | −0.0002 [−0.0011, +0.0007] | −0.0005 [−0.0017, +0.0008] |
| none | m+c | −0.0007 | +0.0001 [−0.0016, +0.0018] | +0.0006 [−0.0019, +0.0030] |
| expression tertile | c | −0.0012 | −0.0003 [−0.0015, +0.0008] | −0.0009 [−0.0025, +0.0006] |
| expression tertile | m+c | −0.0039 | −0.0008 [−0.0042, +0.0031] | −0.0009 [−0.0056, +0.0048] |
| cis < 100 kb vs rest | c | −0.0006 | −0.0001 [−0.0011, +0.0010] | −0.0004 [−0.0019, +0.0010] |
| cis < 100 kb vs rest | m+c | −0.0019 | +0.0002 [−0.0023, +0.0031] | +0.0010 [−0.0028, +0.0051] |
| target strength tertile (donor consensus) | c | −0.0011 | −0.0002 [−0.0016, +0.0012] | −0.0006 [−0.0026, +0.0013] |
| target strength tertile | m+c | −0.0016 | **+0.0026 [+0.0002, +0.0050]** | +0.0023 [−0.0011, +0.0056] |
| predicted direction (up/down) | c | −0.0004 | **+0.0016 [+0.0005, +0.0028]** | +0.0009 [−0.0004, +0.0021] |
| predicted direction | m+c | −0.0030 | +0.0018 [−0.0014, +0.0050] | +0.0001 [−0.0041, +0.0041] |
| called by us | c | −0.0006 | +0.0003 [−0.0007, +0.0012] | −0.0006 [−0.0019, +0.0007] |
| called by us | m+c | −0.0016 | +0.0007 [−0.0020, +0.0038] | +0.0006 [−0.0033, +0.0055] |
| donor \|cons\| tertile | c | −0.0003 | +0.0006 [−0.0003, +0.0014] | +0.0003 [−0.0009, +0.0015] |
| donor \|cons\| tertile | m+c | −0.0031 | −0.0005 [−0.0045, +0.0043] | +0.0011 [−0.0043, +0.0081] |
| donor agrees with our sign | c | −0.0016 | −0.0005 [−0.0033, +0.0024] | +0.0005 [−0.0033, +0.0044] |
| donor agrees with our sign | m+c | −0.0021 | +0.0009 [−0.0029, +0.0055] | +0.0038 [−0.0018, +0.0102] |
| direction × agreement (6) | m+c | −0.0068 | −0.0018 [−0.0058, +0.0021] | −0.0025 [−0.0076, +0.0031] |
| direction × expression (6) | m+c | −0.0074 | −0.0016 [−0.0059, +0.0026] | −0.0038 [−0.0091, +0.0018] |
| called × direction × agreement (12) | m+c | −0.0101 | −0.0016 [−0.0061, +0.0030] | −0.0021 [−0.0084, +0.0047] |

(c-only rows of the 6- and 12-class variables: in-sample −0.0023 … −0.0047, held-out −0.0004 … +0.0016, all spanning 0.)

**Across surfaces, leave one surface out** (`loso480.py`; the a_analog on every surface, MDE add power-normalised by
√(400/n_cells) so c means the same at 400 and 150 cells; per-class c fitted on the other three surfaces):
- Every single-variable class rule is within ±0.007 of a single global c on every held-out surface, with mixed signs
  (expression, predicted direction, donor |cons|, target strength, cis). No class structure transfers between lines.
- The global c the other surfaces prefer (0.16–0.22 normalised) is too large for H1 (+0.014 … +0.022 vs shipped when H1
  is held out). That is the power effect of §4 again: the folds are not calibrated for magnitude.
- The one rule that "transfers" among cancer folds, direction × expression (−0.061 / −0.104 vs global on the held-out
  Orion folds), learns the Orion expression-direction artefact (§5.1). It costs H1 **+0.029** vs global (+0.043 vs
  shipped) and K562 +0.012. It is mostly the two Orion lines sharing one study's blob (highest-expression decile 0.08–0.09
  up, against 0.37–0.47 on every other surface), plus a low-expression up-excess that is partly a power asymmetry
  (H1's lowest decile goes 0.51 → 0.71 at fold power). Neither is biology that A/B/C at H1's power would carry.

**Verdict: no class-based magnitude rule beats the shipped coarse MDE on H1 out of sample, so none passes "H1 AND the
folds".** The class is closed: expression class, target strength, cis position, predicted direction, call status,
donor magnitude and donor agreement, alone or crossed. Every held-out interval spans 0 or is worse. The in-sample bounds
(≤ 0.004 single-variable, ≤ 0.010 with 24 parameters) shrink to ≤ 0.002 held out. On the board, pure magnitude levers
carry ×1 … ×2.3 of H1 (calibration rule 1), so the best point estimate (−0.0018 raw) is worth at most about +0.0005 …
+0.0012 overall, well inside the ordering noise (sd 0.0045). Not worth a slot.

## 7. What the DE channel's arithmetic means would have to become (description only, no emitter code)

The measured result is the absence of a class rule. The one construction that survives is the rider (§4.6). As a
target for the arithmetic means:

- For each (target p, gene g), the delivered per-cell arithmetic mean CPM becomes A'_pg = A_pg · 2^(0.2 · cons_pg),
  with A_pg the current one (after the coarse MDE).
- cons_pg is the mean over the panel donors covering p (K562-GWPS, Orion HCT116, Orion HEK293T) of
  log2((P_pert + 1 CPM) / (P_ctrl + 1 CPM)). P is the donor pseudobulk fraction in CPM, the exact inverse of the donor
  table's log1p-CP50k (`common480.load_donor_rows`, pc 1).
- 0 where no donor covers p or g.

On A/B/C this is an additive lfc term, like the MDE add. The MDE travels through the K-class A-channel at fixed per-gene
totals, so pds/mse stay bit-identical. Two cautions:
- The emitter "keeps every sign". Where 0.2·cons opposes a small lfc of ours, the mixed form wants a sign change that
  the A-channel may refuse; the magnitude-only form (our sign kept) is null at h1-test's power.
- Adds above the MDE's dose over-filled the classes and leaked calls (426/435: c 0.15). 0.2·|cons| exceeds 0.1 on the
  donors' strong genes, so class capacity must be checked on a real build.

Expected board value: +0.001 … +0.005 overall if priced as pure magnitude; the sign part is unpriceable on H1. Smaller
than the ordering noise, so it belongs in a DE-magnitude build as a rider. Measure it on the local H1 emitter (the 435
harness) before any slot.

## 8. Files, reproduction, memory

Scripts (all `uv run python <script>` from this directory; each writes `work/<name>.json`):

| step | script | output | peak RSS |
|---|---|---|---:|
| known answers | `ka480.py` | `ka480.json` | < 1 GB |
| H1 constructions (a)–(e) | `h1con480.py` | `h1con480.{json,npz}` | < 2 GB |
| scorer-rule truth gates | `gates480.py <name> <h5> --blk 1024` | `gate-<name>.npz` | 1.5–2.9 GB |
| K562 fold harness | `k562harness480.py` (env `K562_N_CELLS`, `K562_MIN_CELLS`, `K562_OUT`) | `k562-gwps-panel.h5` (deleted after use, 1.4 GB, 2 min to rebuild) | 0.5 GB |
| power-thinned H1 | `h1thin480.py --cells N --f F --out work/h1-X.h5` | thinned h5 (deleted after gating) | 2.2 GB |
| fold constructions | `fold480.py` | `fold480.json`, `fold480-rows.npz` | 1.2 GB |
| power dose-response | `power480.py`, `power480b.py` | `power480*.json` | < 1 GB |
| replicate anchors | `replicate480.py <name> <h5> <gate>` | `replicate480.json` | 0.8–3.0 GB |
| rider | `rider480.py` | `rider480.json` | < 1 GB |
| small quoted checks | `extras480.py {aplus45,crossstudy,riderclass,riderrobust,updown}` | `extras480-<name>.json` | < 1 GB |
| decomposition | `decomp480.py`, `decomp480b.py` | `decomp480*.json` | < 1 GB |
| class rules | `rules480.py`, `cv480.py`, `loso480.py` | `rules480.json`, `cv480.json`, `loso480.json` | < 1 GB |
| sign value | `signvalue480.py` | `signvalue480.json` | < 1 GB |
| audit, H1 (independent code, 19 rows) | `audit480.py` | `audit480.json` | 0.22 GB |
| audit, folds (independent code, 15 rows) | `audit480b.py` | `audit480b.json` | 0.41 GB |
| board members, nmae-only steps | `board480.py` | `board480.json` | 0.40 GB |
| fid known answer and the fid link | `fid480.py` | `fid480.json` | 0.12 GB |
| pots per surface | `pots480.py` (reads `decomp480*.json`) | `pots480.json` | < 0.1 GB |
| truth expression-direction profile | `artefact480.py` | `artefact480.json` | 0.15 GB |

Inputs read (never modified): `323/work/reference.npz`, `426/work/h1/stats-base_s0_bps.npz`,
`435/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz`, `434/work/allcells434.npz`,
`467/work/strength467.json`, `425/work/tss_v44.tsv`, `098/work/refde_orion.npz`, `data/donors/{k562_gwps,
orion_hct116,orion_hek293t}`, `data/harness/{h1-test,orion-hct116-thinned,orion-hek293t-thinned}.h5ad`,
`data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad`, and for resume 2 `data/leaderboard/*.json` and
`435/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json` (the fid known answer). `data/donors/h1_2025*` was never used.

Rebuild the thinned H1 gates with, for example:
`uv run python h1thin480.py --cells 250 --f 1.0 --out work/h1-c250.h5 && uv run python gates480.py h1_c250 work/h1-c250.h5 --blk 1024`
(levels used: c250/c100/c150 at f 1.0; c400_d56/c250_d56/c150_d56 at f 0.56; h1_thin = 150 cells at f 0.6; seed 0).
