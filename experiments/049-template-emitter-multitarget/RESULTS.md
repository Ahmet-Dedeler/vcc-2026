# 049 — Fixed-margin template emitter across three targets

Date: 2026-08-22

## Result

The fixed-margin template emitter reproduced experiment 043's favorable MED25 result and
extended it consistently to MED15 and PLAGL2.

Across all three targets it:

- preserved every gene's column total and scorer pseudobulk exactly;
- recovered HCT116 control library-shape and detection behavior;
- improved the simple DE-set Jaccard, significant-intersection direction accuracy, and top-500
  overlap on **3/3 targets**;
- worsened the LFC NMAE proxy on **3/3 targets**.

This is stronger evidence that realistic cell layout is useful for the DE metrics, but it is
still not sufficient evidence to replace quasi emission in production. The consistent LFC
tradeoff needs a native same-context prediction and official count-level scoring.

## Critical limitation: H1 prediction labels versus HCT116 truth

The existing experiment-040 prediction file is:

```text
experiments/040-union-gamma/out/runs/g0.75_q1.0_a1.08_spread0.25/pred.h5ad
```

Despite old HCT116 labeling in experiment-040 logs, this file has the H1-test-style layout:
52,400 cells, 18,080 genes, 85 perturbations, and 18,400 controls. It also inherits the
quarantined H1 target-mask mismatch documented in experiment 040.

The real reference used here is:

```text
data/harness/orion-hct116.h5ad
```

It has 38,200 cells, 18,533 genes, 132 perturbations, and 18,400 controls. The two files share
18,077 genes and exactly three perturbations: MED15, MED25, and PLAGL2. Therefore this is an
emitter-layout diagnostic using identical model pseudobulk in both arms, **not** a clean model
or cross-context transfer validation.

The existing arm is specifically experiment 040's `quasi` output with global spread 0.25. It is
the current quasi/even fixed-column family comparator; this experiment does not synthesize a
second separate `even` arm.

## Method

For each target independently:

1. Read only its 400 predicted rows from the backed experiment-040 H5AD.
2. Preserve that block's integer column totals exactly.
3. Select the same 400 HCT116 control templates by library-size quantiles.
4. Reallocate counts with experiment 043's fixed-margin template emitter.
5. Compare quasi and template cells against the target's real HCT116 cells.

DE diagnostics use two-sided SciPy Mann-Whitney U tests against the same 2,000
library-quantile-stratified HCT116 controls, BH correction at 0.05, genes above 5 CPM in
controls, and exclusion of the target gene. This is deliberately a bounded proxy, not
`cell-eval2`.

MED25 uses seed 4300, exactly matching experiment 043. Its reported metrics reproduce experiment
043 bit-for-bit. MED15 and PLAGL2 use deterministic seeds 4400 and 4500.

## Exact pseudobulk

For MED15, MED25, and PLAGL2:

| check | result |
|---|---:|
| Maximum column-total difference | **0 counts** |
| Maximum scorer bulk-log difference | **0.0** |
| Full scorer bulk vector array-equal | **yes** |

The two arms therefore present identical plug-in pseudobulk profiles to PDS and MSE. Differences
below come only from rearranging those counts across cells.

## Per-target results

Lower is better for library-shape KS, detection-rate MAE, and LFC NMAE. Higher is better for the
other three diagnostics.

### Library-shape KS versus HCT116 controls

| target | quasi | template |
|---|---:|---:|
| MED15 | 0.490163 | **0.001630** |
| MED25 | 0.492609 | **0.001630** |
| PLAGL2 | 0.492663 | **0.001576** |
| **macro** | 0.491812 | **0.001612** |

### Gene detection-rate MAE versus the real target

| target | quasi | template |
|---|---:|---:|
| MED15 | 0.223504 | **0.082485** |
| MED25 | 0.223246 | **0.082965** |
| PLAGL2 | 0.227875 | **0.081438** |
| **macro** | 0.224875 | **0.082296** |

### BH<0.05 DE-set Jaccard versus the real target

| target | quasi | template |
|---|---:|---:|
| MED15 | 0.042409 | **0.049381** |
| MED25 | 0.054588 | **0.063025** |
| PLAGL2 | 0.012769 | **0.014795** |
| **macro** | 0.036589 | **0.042400** |

### Direction accuracy in the BH-significant intersection

| target | quasi | template |
|---|---:|---:|
| MED15 | 0.404580 | **0.438746** |
| MED25 | 0.557087 | **0.624444** |
| PLAGL2 | 0.474576 | **0.523810** |
| **macro** | 0.478748 | **0.529000** |

### Top-500 DE set overlap

| target | quasi | template |
|---|---:|---:|
| MED15 | 0.048 | **0.160** |
| MED25 | 0.028 | **0.174** |
| PLAGL2 | 0.040 | **0.120** |
| **macro** | 0.038667 | **0.151333** |

### LFC NMAE proxy

| target | quasi | template |
|---|---:|---:|
| MED15 | **9.48668** | 9.65042 |
| MED25 | **8.82827** | 8.98818 |
| PLAGL2 | **10.20828** | 10.38561 |
| **macro** | **9.50774** | 9.67474 |

## Macro interpretation

Relative to quasi, the template emitter:

- reduced library-shape KS by **99.67%**;
- reduced detection-rate MAE by **63.40%**;
- increased DE-set Jaccard by **0.00581 absolute / 15.88% relative**;
- increased direction accuracy by **5.03 percentage points**;
- increased top-500 overlap from 0.0387 to 0.1513, a **3.91x** ratio;
- worsened LFC NMAE by **0.1670 / 1.76%**.

The template arm reduced the mean predicted BH<0.05 set from 9,248.7 to 7,054.7 genes. It still
over-calls heavily relative to the real targets, whose significant sets contain 416, 551, and
120 genes. The better Jaccard and sign accuracy are therefore meaningful but do not solve DE
calibration.

## Decision

**Promising emitter, not production-approved.**

The multi-target consistency upgrades experiment 043 from a one-target curiosity to a credible
DE-layout improvement. However:

- the H1-versus-HCT116 mismatch prevents a clean estimate of leaderboard value;
- LFC NMAE consistently moves in the wrong direction;
- no official six-metric count-level score was run.

The next scientifically valid test is the same quasi-versus-template comparison on a native
HCT116 prediction artifact, followed by one serial guarded `cell-eval2` run in cloud or another
isolated high-memory environment. This experiment itself makes no production-code or submission
change.

## Resource and validation record

- Final diagnostic wall time: **11.58 seconds**
- Script self-reported peak RSS: **0.973 GiB**
- Outer process-tree guard peak RSS: **1.00 GiB**
- Guard threshold: **3.0 GiB**
- System free memory during run: **87–88%**
- Tests: **4 passed**
- Ruff: **passed**
- `cell-eval2`: **not invoked**
- Submission: **not built or sent**

Command:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 3.0 --min-free-percent 20 --poll-seconds 0.25 -- \
  uv run python experiments/049-template-emitter-multitarget/run.py \
  --max-rss-gib 2.8
```

Machine-readable results are in `results.json`.
