# 547 — Independent check: how much do the DE members move on the board when only the bulk changes?

Coordinator, 2026-09-30 11:38Z. Checks worker 537's §0 A claim (v0, 2026-09-30) before anyone acts on it: the DE
channel's board noise for arms that keep the DE design is far below AGENTS 11.1 rule 13's "sd ≈ 0.0045 overall between
digest-identical arms". `pairs547.py` → `pairs547.json`. Inputs: our 34 captured entries with package metas.

## Method
Pairs of our own board entries whose package metas differ ONLY in bulk keys (`bulk_ctx_scale`, `bulk_load(_sha256)`,
`bulk_sharpen_q/_pt`, and `sameline` when neither side carries a programme). The 458 family is excluded: it puts Arc's
programme, a DE block prior, under `sameline` (a first pass that kept it found reach −0.027, fid −0.0055: that is 458d's
real DE effect, not noise). Overall units = raw / |r − b| / 6 with AGENTS §4's mid-range spans.

## Result (13 pairs)
| member | rms raw | max raw | rms overall |
|---|---:|---:|---:|
| reach | 0.001235 | 0.002314 | 0.000230 |
| nmae | 0.001285 | 0.002282 | 0.000356 |
| fid | 0.000022 | 0.000050 | 0.000012 |
| jac | 0.000016 | 0.000025 | 0.000007 |
| **DE total (quadrature)** | | | **0.00042** |

- The nmae part is mostly SYSTEMATIC, not noise: the four pairs that change `bulk_sharpen_q` move nmae by 0.0013–0.0023 raw;
  the pairs that do not change it move nmae by ≤ 0.00015 raw. Without sharpening changes the DE total is ~0.00025 overall.
- So for arms that keep the emission law, the DE members vary by ~0.0002–0.0004 overall, **10–20× below 0.0045**. This agrees
  with 537's figure (reach rms 0.0013 raw, 0.00024 overall) and with FINDINGS' own notes (468: prior-only arms share emission
  draws and "carry none of the 0.0045 re-roll sd"; H1 seed sd of a paired effect 0.0002).

## What it does not settle
- It measures bulk-only pairs, which keep every DE draw. A DE-knob arm that re-draws the emission carries the re-roll noise
  of its own design; 537 estimates that from H1's three seeds at board scale (~0.0001–0.0002 overall). Its known answer
  (230 §6's 0.0045 reproduced from 230's own pairs as a spread of emission DESIGNS) is due in 537's RESULTS by 10-02.
- Until 537's version is in, AGENTS rule 13 and FINAL-PHASE-FEEDBACK-PLAN §1.4 stay as written; preregistered bins are not
  changed retroactively.

## Correction (2026-09-30T21:09Z, after 537 v4 §2.6)
"It measures bulk-only pairs, which keep every DE draw" above is **wrong**. The builder makes ONE emission stream per file and
each block's ballast count K sets how many carriers that block draws from it, so a bulk change that moves K on any block
shifts the stream for every later block and re-draws pass 1 (537 §2.6, code + H1 test: a bulk-only step re-draws pass 1 on
83–84 of 85 targets; same-seed per-target reach correlation with its base ≈ 0; paired sd 0.00104 vs 0.00094 predicted for two
independent draws). So these 13 bulk-only board pairs ARE re-rolls of the DE emission, and their ~0.0003 overall is a board
measurement of the re-roll (537: reach rms 0.00130 raw on the board vs H1's re-roll scaled by √(85/900), 0.00123 raw). DE
knobs that leave K alone keep pass 1 and read below it. 548's refusal to scale H1 to the board is answered by this; 230's
0.0045 remains the prior spread for a change of emission LAW. Rule 13 is rewritten after 537 §2.5's preregistered check on
the 10-01 … 10-04 reads.
