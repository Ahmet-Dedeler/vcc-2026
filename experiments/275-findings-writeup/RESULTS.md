# 275 — Findings note written (`docs/FINDINGS-2026-09.md`, 2,163 words of prose plus 13 tables); every number it quotes re-read from the saved JSON/CSV it came from — **447 numbers checked, 447 match** (`work/verify275.log`) — and three discrepancies between RESULTS.md headlines and their own tables found and reported below, not fixed. Writing and re-verification only.

**Date 2026-09-06 (03:46–04:15Z; 06:46–07:15 local).** `uv run python` only, to read saved JSON/CSV. **No build, no `vcc`, no
`cell-eval2`, no fold run, no model name, no submission, nothing under `submissions/`, `data/` or `scripts/`.** Files written:
`docs/FINDINGS-2026-09.md` (the deliverable) and, inside this directory, `work/explore275.py` (a read-only JSON path explorer),
`work/verify275.py` (the checklist generator) and `work/verify275.log` (its output). Peak RSS: the verifier reads ≈ 1.7 MB of JSON
and 0.1 MB of CSV; `/usr/bin/time -l` reports 30 MB maximum resident set size for the run (0.04 s wall). No `__pycache__` was
written anywhere. `git status` shows no other change of mine; the files newer than the mandate outside this directory belong to the
concurrently running 274. Wall clock ≈ 30 min of the 3-hour box.

---

## 0. What the note contains, in one table

| section | content | numbers (all verified) |
|---|---|---|
| §1 Setting | the task, the four instruments, the two tracks, what "board-validated" vs "fitted calibration" means | 5.6e-17 (KA-2), 16 s.f., 85 × 400, 82 × 12 |
| F1 | reliability of a 150–200-cell response on the panel vs the Flex-400 reference; the guide construct; reliability as a target property | 267 §3 table R (medians +0.047 / +0.028 / +0.028; 5.6 / 32 / 34 % ≤ 0; SB 0.09 / 0.05 / 0.05), 271 §1.2 D1, 273 §2 (+0.243, 0/85; +0.114 at 80/80; guide +0.014), 159's 0.330 and 0.478 |
| F2 | the donor-selection oracle is reference noise; what survives | 267 §4 (+0.056 / −0.018 / −0.022; retained −0.32 [−0.49, −0.16]; T3 +0.021 / +0.027), 263 §4–5 (78 cells, 44–61 %), 265 (p 0.011, z +2.3, +0.0121 [−0.0035, +0.0278]; cancer p 0.68 / 0.12 / 0.50) |
| F3 | the operator's in-sample noise learning; the spectrum; un-learning on PLOLO and on Flex; reliability weighting | 272 §1.2 S1 (76.9–78.7 %, 87.5–91.5 %, Σn/Σλ 0.88–1.00), 271 §1.2 / 272 §2.2 (−0.13…−0.24 vs +0.007…+0.022), 272 §1.3 / §2.3 arms, 273 §4.1–4.2 |
| F4 | the 0.02 knee on both references; regimes; deficit decomposition; fold consistency | 262 §4.2 / §5 (T5c, T5d), 273 §3 (slope 6.0, 12.9 / 48 / 39 %) |
| F5 | Broad-16 variance decomposition, rank, batch, LOLO | 268 §2–4 (0.091 / 0.115 / 0.794; null 0.028; rank 27; c1 +0.018 vs +0.001 / −0.007; full-response gain 0) |
| F6 | K562-GWPS DE sparsity and the matched-control check | 260 §2c (median 5, 85/272, 58 zero; 400 cells 7 / 9 / 6; KA 1,205 / 425 / 342), 266 §3–4 (6.0 / 5.7 / 6.3; 17 / 17 / 18) |
| F7 | the two board points; class factors; PLOLO's calibration; the σ correction | BR237 (−0.004783 → −0.00716; k_cos 0.962), BR249 (+0.01492 raw = 1.22×; k_cos 3.76; +0.0042634), 258 §3, 257 R1 (+0.00472), 162 (0.001494, 6/6, p 0.0028), 265 |
| §3 | established / not established: unit resampled, n, replication, the preregistered test on the reserved data | 15 rows |
| §4 | implications, ≤ 10 lines per track | — |
| §5 | methods: the statistic, split units, bootstrap and permutation conventions, the pessimistic pricing rule, the leakage rules | — |

## 1. The checklist (`work/verify275.py` → `work/verify275.log`)

One row per number: id, source file and JSON path (or the CSV column and the computation), the stored value, the value as quoted
in the note, OK if |stored − quoted| ≤ the row's tolerance (the rounding the note uses), and the RESULTS section the number comes
from. **447 rows, 447 OK** — 80 for F1, 54 for F2, 114 for F3, 65 for F4, 44 for F5, 40 for F6, 50 for F7. The sources:

| file | what was re-read |
|---|---|
| `experiments/267-…/work/analysis267.json` | reliability by fold × split (medians, p10/p90, share ≤ 0, cells per half, by 262 group, ρ with `pds` and with cells), the four oracle tables (in-sample, cross-validated, random, cross − random, retained fraction, abstaining), strata (group (e), T3 on both splits), pick agreement, sanity |
| `experiments/271-…/work/{analysis271,diagwht271,price271}.json` | donor-row reliability D1 (raw, SB, whitened), the estimator arms and their controls, retention, the price of `ebrow`, the BR249 retrodiction |
| `experiments/272-…/work/{analysis272,prod272}.json` | the spectra of the four operators, in-sample and cross-fitted anti-correlation, the Wiener / floor / hard arms and controls, construction B on two folds, the production families' energy kept |
| `experiments/273-…/work/{analysis273,reliability273,instrument273}.json` | Flex reliability by unit, KA-1 / KA-2 / KA-2a, the emission gap, the conversion curve and regimes, the (a) / (b) / (c) re-reads and their controls |
| `experiments/262-…/work/{analysis262,extra262}.json` | the conversion curve and deficit shares, the partition by mechanism, covariate ranking, bottom-decile overlap and cross-fold Spearman, per-fold `pds` |
| `experiments/265-…/work/boot265.json` | the iPSC-fold bootstrap and permutations, the seed diagnostic, the cancer-fold permutations and pooled bootstraps, 263's two-donor cells |
| `experiments/268-…/work/{decomp268_cpm5,lolo268_K2}.json` | the variance decomposition (cell halves, guide halves, OR null), batch test, SVD rank, components, the LOLO arms and controls |
| `experiments/260-…/work/*.csv`, `de_power_summary.json`; `experiments/266-…/work/analysis266.json` | medians and gate counts recomputed from the per-target CSVs; the 400-cell curve; the strong-target KA; 266's nine scheme × seed cells and diagnostics |
| `data/board-entries/{iUziqE1gewfcPKzYUDaB,ydAemVdWmsTx9xhDWcou,Xv0ayvuoHpsKYWv2iu42}.json` | the six raw and scaled members of `h1union`, ARM6 and `cd4svd5`; the member deltas, k_pds and k_cos recomputed from them |
| `experiments/257-…/work/price257.json`, `experiments/162-…/work/retrodict.json`, `experiments/159-…/work/{reliability,crossline}.json` | the fitted-factor retrodiction of BR249 and PLOLO's Δ`pds`; PLOLO's LOO RMS, sign count, permutation p; 159's measured 200/200 value and its fitted 400-cell extrapolations |

**Known-answer rows of this write-up** (the checks that would have caught a mis-transcription): the three board entries' member
deltas recompute BR237's overall −0.0036403 and its `nmae` leg −0.00716 (ratio 1.497 to the identity's −0.004783) and BR249's
+0.0042634 with `pds` +0.00545 and `mse` −0.00080 at 5e-8 / 5e-6; k_pds 1.22 (BR249) and 2.48 / 3.79 (ARM6) and k_cos 3.76 /
0.962 recompute from the entries and the saved PLOLO increments; 257's, 271's and 272's price recipes agree on +0.004722; 266's
pooled seed-0 counts equal 260's on 50/50; 273's KA-2 reproduces cell-eval2's per-target `pds` at 5.55e-17 with the mean identical
to 16 digits.

## 2. Discrepancies found between a RESULTS.md headline and its own tables or saved files (reported, not fixed)

1. **R159 §1.3 calls h1-test's 400-cell residual reliability "0.481 (measured, not fitted)"; it is fitted.** 159's
   `work/reliability.json` holds no 400-cell measurement for h1-test — the largest half is 200 cells per side, where the measured
   residual reliability is 0.3300 (`h1test/resid_disjoint/r_mean[9]`); its linear-law extrapolation to 400 is 0.478
   (`h1test/resid_disjoint/linear/rho_400`) and its power-law one 0.458. Spearman–Brown from the measured 0.330 gives 0.496 (mean) /
   0.391 (median) (R273 §2 v). The number is repeated as "MEASURED at 400" in R271 §5 (table Q and the text) and in READ-FIRST §6
   ("an H1-like Flex reference 0.48 (√ρ 0.69)"). R273 §2 (v) flagged it, but says the 0.481 "is not derivable from its h1-test rows —
   the only 0.480 in 159's json is an HCT116 full-space extrapolation" (`crossline.json` `rho150/A_hct116/full/power/rho_400` =
   0.4802); the nearer source is h1-test's own linear-law 0.478. The note quotes 0.330 measured and 0.478 fitted (F1 limit).
2. **R273 §4.1 mislabels the production energy kept by the floor family.** Its parenthetical reads "floor c 0.25 / 1 / 4 = 95 / 83 /
   68 % of the production rows' (272 §3)"; `prod272.json` `families/floor:*/whitened_energy_kept` gives c 0.25 = 94.9 %, c 0.5 =
   83.1 %, c 1 = 67.7 %, c 2 = 49.9 % (R272 §3's own "94.9 / 83.1 / 67.7 / 49.9 / 32.4 % at c 0.25 … 4"), so the three values belong
   to c 0.25 / 0.5 / 1, not 0.25 / 1 / 4. The note uses 68 % for c = 1.
3. **R272 §2.3 (table B1 and the headline) prints the proper cross-fit interval as [−0.0050, +0.0057]; the stored upper bound is
   0.0056477**, which rounds to +0.0056 (`analysis272.json` `refit/orion_hct116/arms/half:xfit|half:inc/boot/hi`). A rounding slip of
   one unit in the fourth decimal; the note quotes +0.0056.

Smaller notes, not discrepancies: (i) R271's headline gives the K562-GWPS donor rows' median reliability as +0.047 (267's reference
value); its own D1 rows read +0.048 (0.0476 on both folds where K562 is a donor) — the note uses +0.048 for the donor row and +0.047
for the reference. (ii) R262 §0 states two different anti-aligned shares that a reader could conflate: 17.5 % of targets with
base·truth ≤ 0 carrying 53.1 % of the deficit (T5d) and 16 % with delivered own cosine ≤ 0 carrying 57 % (the conversion curve);
both recompute from `extra262.json` (139/853 = 16.3 %, 57.1 %; 149/853 = 17.5 %, 53.1 %). The note uses the base·truth pair and
says which. (iii) R263's "8 σ" and R264's "σ" columns are already corrected by R265 and by R264's addendum; the note carries the
corrected statement (2.3 σ on the permutation null). (iv) The mandate's own summary line for F3 ("cross-fit +0.02 vs −0.14/−0.24")
matches the HCT116 fold's numbers (+0.022 / +0.007 vs −0.142 / −0.243).

## 3. What the note leaves out, and why

- 269 and 270 enter only through their headline numbers (269: +0.0001 [−0.0024, +0.0028] on held context and family; the leakage
  rule of R269 §5 with +0.033 / −0.014 [−0.023, −0.005]; 270: AUC 0.89 [0.77, 0.97], ρ +0.36 [+0.00, +0.63], 40 targets, 4 anti-aligned).
  Their RESULTS were read for those sections only; they are not among the seven findings the mandate names, and the note's §3 table
  records 270's feature as a weak positive to confirm on held targets, not as a finding.
- The Flex-400 statements (F1, F3, F4) rest on one embryonic line and 85 targets; the note says so in every limit line and in §4.
- The board calibrations (F7) are one point per class; the note labels every retrodiction as a fitted calibration.
- Nothing was re-measured. Where a saved number and a RESULTS.md sentence disagreed, the saved number is quoted and the disagreement
  is listed in §2; no RESULTS.md was edited.
- Not verified against a saved file, quoted from RESULTS/BOARD text only: BR237's class band +0.0022 … +0.0038 and its direction-leg
  sum +0.00352 (overall minus `nmae`; `pds` + `mse` recompute from the entries at +0.003356), P249's predicted bands, 243's 0.107332 cosine, 258's DE-three rates
  0.0425 / 0.185, the 0.24 recipe, 271 §5's 400-cell law table, 272 §1.1's half-difference signal shares (0.105 / 0.054 / 0.064),
  R263's 78-cell and 44–61 % summary, R188's median-89 and R221's 58/85 comparators, AGENTS §8.1b's sd 0.004502, and the reserved
  Broad-16 line names. Each is cited to its section.

## 4. Files

| file | what |
|---|---|
| `docs/FINDINGS-2026-09.md` | the deliverable: setting, F1–F7, the established / not-established table, implications, methods appendix |
| `work/explore275.py` | read-only JSON explorer used to find the paths (`keys <file> [depth]`, `grep <file> <regex>…`) |
| `work/verify275.py` | the checklist generator: 447 (label, file, path, quoted value, tolerance, section) rows plus the computed board deltas and class factors |
| `work/verify275.log` | its output, the markdown table of all 447 rows with stored vs quoted values |

---

## Verdict

**DONE — writing only.** `docs/FINDINGS-2026-09.md` states seven findings with the unit, n, statistic, interval and limit of each,
a table of what is and is not established with the preregistered test each could get on the reserved data, ten-line implications
per track, and a methods appendix; all 447 numbers it quotes match the saved files they came from, and the three places where a
RESULTS.md headline disagrees with its own tables (159's "measured" 0.481, 273's floor-family labels, 272's rounded upper bound) are
reported above and left as they are. No build, no submission.
