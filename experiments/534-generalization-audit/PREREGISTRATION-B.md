# 534 Part B: the sealed mock final (bulk channel). Preregistration

Written 2026-09-29 13:53Z by the session "Public data scope [2e83f3]" at Ahmet's request, before any fresh panel exists.
Coordinator OK and heavy window: its message of 2026-09-29 ~13:55Z (free until 2026-09-30 00:30Z).

## 1. Question

Does the October bulk recipe's gain hold on Arc-like panels that nobody has ever used for a decision?

The recipe was built step by step on A/B/C's board and re-priced on three October-like panels (`real`, `efp`, `efr`;
438, 522, 540–546). Those three panels have since informed dozens of choices, so they are tuned-on: their readings are
optimistic by an unknown amount. This part builds four fresh panels, seals them before drawing, and prices the recipe's
milestones on them once.

Scope: the bulk channel only (pds, mse). That is about 73 % of the current overall (446bd2: pds 0.137 + mse 0.061 of
0.270). The DE members are out of scope: there is no validated cancer-line DE surface (AGENTS 6.5; FINDINGS 8.1h).

## 2. Reserve R6, declared here before drawing

Four panels of 300 targets, drawn with 238's `bundle238.py` rules: the same gene axis, the same relabelled A/B/C control
files, the essential list `data/donors/k562_essential/perts.txt` subtracted. There is one change. The eligible pool also
excludes every target of every panel file on disk at declaration time, plus each earlier R6 draw. So the four panels
are disjoint from every tuned-on panel and from each other.

Panel files excluded (sha256 prefix of the sorted target list; 7 distinct panels in 9 files):
`data/challenge/pert_counts.csv`; `experiments/180-def-readiness/work/def-bundle/`;
`experiments/199-def-runbook/work/{def-bundle, def-bundle-r, def-bundle-badaxis}/`;
`experiments/238-essential-free-panel/work/{ef-bundle, ef-bundle-r}/`;
`experiments/247-oct22-guards/work/census/collapsed-bundle/`; `experiments/323-current-family-validation/work/bundle/`.
The drawing script re-runs `find experiments data -name pert_counts.csv` at draw time and excludes all of them.

| panel | rule (238) | pool | seed |
|---|---|---|---|
| S1 | efr rule | K562 ∩ HCT116 ∩ HEK293T ∩ CD4 atlas | 534101 |
| S2 | efp rule | K562 ∪ HCT116 ∪ HEK293T | 534201 |
| S3 | efr rule | as S1 | 534102 |
| S4 | efp rule | as S2 | 534202 |

Drawn in the order S1, S2, S3, S4. Each draw's target list and sha256 go into `work/partB/draws.json` at draw time.

## 3. Access rule

R6 serves exactly two claims:
- **C10a** (now): the milestone table in §4.
- **C10b** (once, after the recipe's last A/B/C read on 10-05): the final October recipe's exact bulk configuration.

No fitting, selection, dose search or debugging may read an R6 panel's scores. The known answers (§5) run on `efp`,
never on R6. Every access is appended to `docs/CONFIRMATION-LEDGER.md` under R6 / C10.

## 4. Designs (lib438 names; 522's route for the cap)

Cumulative milestones, each priced against the one before it:

| id | design | what it adds | October switch |
|---|---|---|---|
| M0 | `PP\|K3` (= `h1c1\|K3`, q 1) | the 410b-era bulk: per-target 410 key, no sharpening, no stack | baseline |
| M1 | `h1c1\|K3+sh1.25` | sharpening q 1.25 (417a) | q |
| M2 | `h2c1\|K3+sh1.25` | H1 ×2 stack (417e's H1 half) | cd4_stack (H1 half) |
| M3 | `h2c1\|K3+sh1.5` | q 1.5 (433k) | q |
| M4 | `h2c1\|K3+sh1.75` | q 1.75 (446b) | q = 1.75 |
| M5 | `h2c1\|V1_5_3+sh1.75` | LODO key, clip (0.5, 3) (446bj2) | lodo = 1 |
| M6 | M5 with the percentile cap P = 96.64697522769474 (522's `pctPROD`) | the high cap (495/522) | lodo_hi_pct (PENDING) |

Controls:
- **C-a:** `h2c1.5|K3+sh1.75` vs M4 (CD4 ×1.5, the half that 438 found A/B/C-only). Predicted ≤ 0.
- **C-b:** the LODO key permuted over targets (438's CTRL construction) vs M5. Predicted < 0.

## 5. Known answers (run before any R6 number exists; abort on failure)

- A copy of `lib438.py` with `WORK` and every cache path pointed at `experiments/534-generalization-audit/work/partB/`,
  so nothing is written into 438's or 522's directories.
- On `efp`: M0–M5 and C-a reproduce 438's stored values, and M6 reproduces 522's `pctPROD` row, to |d| ≤ 1e-12
  (438 and 522 report |d| = 0 on their own rows).
- Every R6 fold operator passes `plolo.Operator._check_axis` for its own panel. Its `axis_sha256` must differ from
  A/B/C's (`edf6ba88…`) and from every tuned-on panel's.

## 6. Units and statistics

- Primary: 438's LEVELLED route at k = 1, Δoverall in board units.
- Projection: the calibrated transfer from the record: shape steps k_sh 0.65 (mse) / 0.42 (pds); composition k 1.3;
  key steps k 1.31.
- Per panel: the mean over its three cancer folds (K562, HCT116, HEK293T held out in turn).
- Pooled over R6: a paired target bootstrap, 4,000 draws, seed 534. Each step reports Δ, n+ folds and the 95 % interval.
- Reported beside the same step's Δ on A/B/C, `real`, `efp` and `efr`, taken from the 438/522 records (not recomputed).

## 7. Reading rules (frozen now)

1. **Step verdict** (R6 pooled):
   - TRANSFERS if Δ > 0 and the 95 % lower bound > 0;
   - FRAGILE if Δ > 0 and the lower bound ≤ 0;
   - FAILS if Δ ≤ 0.
2. **Retention** = R6 pooled Δ ÷ A/B/C Δ, per step and for the cumulative M0 → M6.
3. **Tuned-on check.** `real`/`efp`/`efr` count as OVER-READ if, for ≥ 2 of the three October switches in scope (M3→M4,
   M4→M5, M5→M6), the R6 pooled Δ is below half the mean of their `real`/`efp`/`efr` Δ.
4. **Controls.** C-a and C-b must read ≤ 0 pooled. If C-a reads > 0 with its interval above 0, 438's CD4 ×1 decision
   is flagged to the coordinator. Nothing is changed here.
5. **Consequence.**
   - An October switch that FAILS on R6 goes to the coordinator as a candidate to switch off.
   - One that is FRAGILE is flagged.
   - This experiment changes no production file, recipe value or slot.

## 8. Prior expectation (from 438 / 522, written before the draw)

| step | expected on R6 |
|---|---|
| M1 | > 0, retention 0.4–0.65 |
| M2 | small > 0, likely FRAGILE (438: +0.0003 … +0.0051 by panel) |
| M3, M4 | > 0, retention ≈ 0.5 |
| M5 | > 0, retention ≈ 1 |
| M6 | > 0 (522: +5 … +6 cosine points over cap 6 on October-like panels) |
| C-a, C-b | < 0 |

A surprise is any sign flip on M3–M6, or a cumulative retention below 0.3.

## 9. Memory and time

- A heavy process starts only when ALL of these hold: `memory_pressure` ≥ 40 % free; `vm_stat` free+inactive ≥ 33 GiB;
  `def_pipeline.other_builders()` is empty; no bash process runs `run_arm_generic.sh`, `run_ka_*`, `def_day.sh`,
  `chain4NN`/`chain5NN*.sh`, `relevel_sp2.sh`, `run_dump4*`, `build527.sh` or `fix526.sh`.
- Heavy processes are operator fits and scoring. Operator fits run `scripts/build_union_whitener.py` at
  `--max-rss-gib 14` (~12–13 GiB `ru_maxrss`, as 238 measured); scoring stays ≤ 9 GB.
- Strictly serial, each under `/usr/bin/time -l`.
- No heavy process starts after 2026-09-29 23:45Z, or during 2026-09-30 00:30–06:00Z. If the work is unfinished, write
  state and resume after 06:00Z.
