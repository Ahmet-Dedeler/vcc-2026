# 656 — re-price the same-line dose w under production's zero repair (PREREGISTRATION)

Council session [7772dc], assigned by the coordinator after 654's sweep. Written and committed BEFORE any 656 run.
Class: transfer (same-line rows, 450's mechanism). This is a re-price of a DECIDED October value with one changed construction
detail. Low priority: it runs in gate602 behind 655 and 618's KA.

## Why
- `sameline_w = 0.5` (DECIDED 09-27; live on Oct 22 whenever a D/E/F line matches a bank line) rests on 450 §10's
  weaker-bank table. That table is seed 0 only, run without zero repair, in board units vs 435's base:

  | bank (flip p) | w 0.5 | w 1 | w 2 |
  |---|---:|---:|---:|
  | p 0.2 (gate signs ≈ 0.66) | +0.0323 | +0.0316 | +0.0268 |
  | p 0.3 (≈ 0.61) | +0.0166 | +0.0086 | +0.0007 |

- 654's sweep found scored zero-total cells in those arms, rising with w:
  - p 0.3: w 0.5 / 1 / 2 = 1 / 2 / 4 in the f 0.35 pools;
  - p 0.2: w 0.5 / 1 / 2 = 1 / 1 / 2.
- Production floors them (zero repair, −5). Without the repair the harness charges ~+30 error each, which is ~0.0014
  overall-eq per cell by 654's scale. So the higher w were penalised more, and the w choice may have been biased toward 0.5.

## Arms (one change vs 450 §10: `--zero-repair -5`)
- `run656c.py` is a copy of 450's `run450c.py`. It adds 450's dir to `sys.path` and a local numba cache; outputs go to
  656's `work/`. Diff recorded.
- w ∈ {0.5, 1, 2} × flip p ∈ {0.2, 0.3} (`RUN450_FLIP`, `RUN450_NOISE_SEED=0`: the same weakened bank as 450) × emission
  seeds {0, 1, 2}. CLEAN4 default prior, full mode (rows + per-target prior). 450's flags `--bulk ps --mde-c 0.1
  --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse`, plus `--zero-repair -5`. 18 runs.
- zs is NOT added. 450 §10 had none, so the one change is the repair. Production also has zs; that stays a recorded
  difference.

## Known answer (structural, per seed-0 arm)
- Each zero-repair seed-0 run equals 450's stored non-repaired seed-0 run (q, p, lfc, sums) on every target without a
  non-own zero total. They may differ only inside those targets.
- A failure stops the read.

## Measure
- 450's own board units (`read450.py`'s SLOPE on the harness's raw six members at f 1, vs 435's base, paired by seed), for
  every arm.
- The decision quantity per bank: Δ(w) = overall(w) − overall(0.5), paired by seed. Report its mean and per-seed sd.
- Per-arm zero-total counts (non-own; inside the f 1 / f 0.35 pools).

## Decision rule (fixed now)
- **The October value moves** from 0.5 to w* ∈ {1, 2} only if, at BOTH p 0.2 and p 0.3, the mean Δ(w*) > 0 AND
  > 2 × its per-seed sd. If both w qualify, take the larger mean at p 0.3.
- If that holds at p 0.2 only: record it. A per-context w needs builder support, so there is no October change.
- **Otherwise 0.5 stands**, now with zero repair and 3 seeds behind it.

## Prediction
- 0.5 stands. At p 0.3 the seed-0 gap to w 2 (0.016) is ~4× the expected artifact (~0.004).
- At p 0.2 w 1 ties 0.5 (gap 0.0007 at seed 0, inside the artifact), and w 2 stays below.
