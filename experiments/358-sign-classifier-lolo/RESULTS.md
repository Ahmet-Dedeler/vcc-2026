# 358 — Label-free sign classification adds nothing over donor consensus: DEAD

**Status: NEGATIVE (COMPLETE).** Both classifiers lose pooled weighted sign vs
the cross-line diagonal (−0.0013 LR, −0.0011 HGB); shuffled-recipient basal
matches or beats own-recipient basal; the regulon prior contributes
+0.000004. All three PLAN.md kill clauses fire. No follow-up, no production
change, no submission.

## Known-answer reproduction (own code, `ka358.py`)

Independent reimplementation (own loader, split, diagonal, endpoints — no
import of core334/core341/run334) of 341's pooled means over 6 lines ×
1,024 targets × 512 coordinates:

| arm / endpoint | reproduced | recorded (341) | diff |
|---|---|---|---|
| source-mean weighted sign | 0.589398 | .589398 | 2.4e-07 |
| source-mean top-64 | 0.642517 | .642517 | 3.0e-07 |
| source-mean own cosine | 0.195532 | .195532 | 6.9e-08 |
| source-mean retrieval | 0.790095 | .790095 | 6.6e-08 |
| cross-line diagonal weighted sign | 0.589398 | .589398 | 2.4e-07 |
| cross-line diagonal top-64 | 0.642517 | .642517 | 3.0e-07 |
| cross-line diagonal own cosine | 0.198735 | .198735 | 9.2e-08 |
| cross-line diagonal retrieval | 0.794614 | .794614 | 1.5e-07 |

**8/8 MATCH** (tolerance 5e-06; residual is rounding of the recorded values).
The report phase additionally recomputes the diagonal from saved
consensus/gain and reprints 0.589398 / 0.794614 exactly, validating the
endpoint path end to end.

## What was measured (inputs, estimator, controls)

Six-line source frame (`selection.json`: 1,024 targets, 512 coordinate genes),
unit-RMS response directions, LOLO over lines × the frame's two 512-target
blocks (train = other half, donor test-half responses never touched — 341's
leakage regime). Per (outer, block): 1,310,720 train rows (5 inner recipients
× 512 targets × 512 genes; consensus/agreement over the other 4 donors,
leakage-matched), 15 label-free features (donor consensus, its sign, |·|,
agreement fraction; recipient `b_ctrl` log1p-scale basal of g and of t with
line-median imputation + missing indicator for 109–138/1024 absent targets in
4 banks; donor-mean basal and recipient−donor diff; DoRothEA ABC prior:
`is_tf(t)` (19/1024), t→g presence+sign (50 pairs, 17 signed), log1p shared
regulators (4.1% nonzero), g→t presence+sign (82 pairs, 20 signed)). Label =
sign of observed response, weight = |response| (all 512 coords; 334 defines no
noise floor). Labels 49.3% positive overall (45.5–51.8% per line), so no class
weighting. Delivered object = predicted sign × |5-donor consensus|
(magnitude-preserving, §8.1e-style: retrieval guard differs from diagonal only
on flipped signs). Fits: (1) cross-line diagonal; (2) L2 logistic (standardised,
lbfgs); (3) HGB classifier (150/0.05/15-leaf/4-deep, mirror of core334.model);
(A) LR+HGB with gene-indexed recipient basal columns under one fixed
derangement (donor-mean untouched; target-indexed `basal_t` untouched —
disclosed limit); (B) LR on donor-consensus cols only; (B2) LR donor+prior;
(C) LR on full features with train labels under one fixed target derangement
per (outer, block). Endpoints per held line over all 1,024 targets:
|response|-weighted sign (primary), top-64 weighted sign, retrieval guard;
pooled mean + per-line means + 10,000-draw paired target bootstrap and crossed
line/target intervals (core341.interval reimplemented, seed 3580909).

## Results tables (units: sign accuracy / retrieval fractions; deltas pooled)

| arm | weighted sign | top-64 | retrieval |
|---|---:|---:|---:|
| diag (baseline) | 0.589398 | 0.642517 | 0.794614 |
| lr_full (fit 2) | 0.588109 | 0.647172 | 0.783542 |
| hgb_full (fit 3) | 0.588313 | 0.645591 | 0.785501 |
| lr_shufrecip (A) | 0.588280 | 0.647324 | 0.783202 |
| hgb_shufrecip (A) | 0.588774 | 0.645957 | 0.784385 |
| lr_donoronly (B) | 0.588646 | 0.646705 | 0.790249 |
| lr_donorprior (B2) | 0.588650 | 0.646663 | 0.790323 |
| lr_shuflabel (C) | 0.567261 | 0.618409 | 0.743884 |

| comparison / endpoint | mean | target 95% | crossed 95% | pos lines |
|---|---|---|---|---|
| lr_full − diag / w | **−0.001288** | [−0.002165,−0.000421] | [−0.006342,+0.003522] | 2/6 |
| hgb_full − diag / w | **−0.001085** | [−0.001903,−0.000261] | [−0.005809,+0.003355] | 2/6 |
| lr_full − diag / top64 | +0.004655 | [+0.003247,+0.006041] | [−0.004217,+0.013652] | 5/6 |
| hgb_full − diag / top64 | +0.003075 | [+0.001726,+0.004398] | [−0.004584,+0.010347] | 5/6 |
| lr_full − diag / retr | −0.011072 | [−0.012505,−0.009642] | [−0.015730,−0.006586] | 0/6 |
| hgb_full − diag / retr | −0.009113 | [−0.010823,−0.007447] | [−0.014289,−0.004366] | 0/6 |
| lr_shufrecip − lr_full / w | +0.000171 | [−0.000009,+0.000347] | [−0.000246,+0.000586] | 3/6 |
| hgb_shufrecip − hgb_full / w | +0.000461 | [+0.000068,+0.000849] | [−0.000166,+0.001101] | 6/6 |
| lr_donoronly − diag / w | −0.000752 | [−0.001650,+0.000142] | [−0.006402,+0.004823] | 2/6 |
| lr_donorprior − lr_donoronly / w | +0.000004 | [−0.000072,+0.000082] | [−0.000128,+0.000129] | 3/6 |
| lr_shuflabel − diag / w | −0.022137 | [−0.023798,−0.020467] | [−0.030609,−0.014478] | 0/6 |

Per held line, weighted sign (diag / lr_full / hgb_full / Δlr / Δhgb):

| line | diag | lr_full | hgb_full | Δlr | Δhgb |
|---|---:|---:|---:|---:|---:|
| k562_gwps | 0.588676 | 0.576782 | 0.577706 | −0.011894 | −0.010969 |
| orion_hct116 | 0.579402 | 0.577183 | 0.578203 | −0.002219 | −0.001199 |
| orion_hek293t | 0.563000 | 0.561180 | 0.561754 | −0.001819 | −0.001246 |
| rpe1 | 0.609090 | 0.617384 | 0.616843 | +0.008295 | +0.007753 |
| nadig_jurkat | 0.591160 | 0.589174 | 0.588369 | −0.001987 | −0.002792 |
| nadig_hepg2 | 0.605059 | 0.606952 | 0.607003 | +0.001894 | +0.001944 |

## Verdict against the PLAN.md kill criterion

DEAD on all three clauses: (1) fit 2/3 minus diagonal pooled weighted sign =
−0.001288/−0.001085 ≤ +0.005 (negative, not just small); (2) control A is
within — indeed slightly above — its fit (shuffled ≥ own on both learners, so
recipient basal columns carry zero resolvable signal); (3) there is no positive
gain for control C to explain, and C's −0.022137 proves the only learnable
structure is the donor-label consensus the diagonal already captures. The ALIVE
bar (+0.015 with A resolved below) is missed by ~0.016. Secondary notes, all
consistent with dead: B2−B = +0.000004 (regulon prior inert); the top-64 lift
(+0.0047, crossed CI unresolved) is already 90% present in donor-only refit
(+0.0042), i.e. not recipient-driven; retrieval guard loses ~0.01 on 0/6 lines
(sign flips cost discrimination); the lone positive line (rpe1, lowest depth)
is matched by shuffled-recipient (+0.009) and donor-only (+0.010), so it is a
donor-refit artefact, not basal conditioning. Mechanism in one line: with
balanced labels and all-positive diagonal gains, the consensus sign is the
sufficient statistic, and 15 label-free features cannot improve it.

## Limits

Six-line development frame only — not an A/B/C forecast; no score gain is
claimed. Agreement-as-fraction reconciles train (/4) vs test (/5) donor counts
by construction. `basal_t` untouched by control A's gene derangement
(target-indexed; disclosed). b_ctrl is log1p(50000·p) scale (build scripts),
used as the "log CPM" feature. Prior restricted to DoRothEA ABC static
indicators; STRING/CollecTRI not tried (same sparsity regime expected).
Intervals use seed 3580909, not 341's seed — point estimates are seed-free.
HGB/LR hyperparameters fixed a priori (core334 mirror); no tuning was run, and
tuning is not motivated by a −0.001 deficit.

## Exact commands to reproduce

```bash
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/358-sign-classifier-lolo/ka358.py       # KA: expect KA_OVERALL: PASS
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/358-sign-classifier-lolo/run358.py cache # ~1 s, writes work/cache358.npz
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/358-sign-classifier-lolo/run358.py run   # 132 s, 12 resumable pred files
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/358-sign-classifier-lolo/run358.py report# prints tables, writes results358.json + per-target358.csv
```

## Peak memory and wall time

Cache 0.7 s; run 132.5 s wall (12 blocks × ~11 s: 5 LR + 2 HGB fits each);
report < 60 s. `memory_pressure` free 70–72% throughout (start 71%, end 72%);
largest live objects two 68 MB float32 feature blocks per block, freed per
block — well under the 12 GiB cap; no process approached the 4 GiB gate.
Scratch retained in `work/`: `cache358.npz` (24 MB), 12 prediction files
(~12 MB each, keep: they are the result), `results358.json`,
`per-target358.csv` (6 lines × 1024 targets × 24 metrics).
