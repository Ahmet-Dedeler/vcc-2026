# 067 — Context-control variance transport

Date: 2026-08-23

## Decision

**Reject.** Control-variance transport is not a reliable addition to the
source-only union-whitened predictor.

The leakage-safe selector chose three different candidates across the three
outer folds, including identity for HCT116. Applying each replicate-A-selected
candidate to its held donor gave macro PDS deltas of **-0.000249419** with
control replicate A and **-0.000394082** with the frozen replicate-B control
statistics. No selected nonidentity result exceeded the 95th percentile of its
gene-permuted ratio null.

## Protocol

- Outer held donors: `k562_gwps`, `orion_hct116`, `orion_hek293t`.
- Gallery: 259 VCC targets measured in all three donors.
- Prediction: fixed equal average of the two source-donor effects.
- Geometry: experiment 062's fixed source-only union whitener, gamma 0.65.
- Removed coordinates: all 300 VCC target genes.
- Control features: two disjoint deterministic samples of 512 cells per donor.
- Scale estimators:
  - winsorized SD of per-cell `log1p(CPM)`, clipping each gene at 1%/99%;
  - scaled MAD, `1.4826 * median(abs(x - median(x)))`.
- Transport ratio:

  ```text
  r[g] = clip(scale_target[g] /
              geometric_mean(scale_source_donors[g]), 0.5, 2.0)
  ```

- Grid: eta `{0, 0.25, 0.5, 1.0}` and placement immediately before or after
  the fixed whitener. Eta zero is represented once by exact identity.
- Selection: for each outer fold, predict each source donor from the other on
  the same deterministic 96-target subset. Maximize the worse directional
  delta, then the mean delta. Held-donor perturbation effects are not available
  to this selection.
- Stability: replicate A selects; replicate B evaluates the frozen choice.
- Null: 100 outer-fold permutations of the selected gene-ratio coordinates.

This is distinct from basal-mean weighting. Donor weights do not change and no
control mean is multiplied into the prediction; only the coordinates of the
already-built effect vector are rescaled by control-cell variability.

## Primary fold results

All values below are raw effect-level rank PDS over the 259-target gallery.

| held donor | source-only selection A | selection B | baseline | selected PDS A | delta A | selected PDS B | delta B | null p95 | empirical p |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| K562-GWPS | `mad:post:eta=0.5` | identity | 0.799796474 | 0.801726976 | +0.001930502 | 0.802011314 | +0.002214839 | +0.005030529 | 0.990099 |
| Orion HCT116 | identity | identity | 0.800574661 | 0.800574661 | +0.000000000 | 0.800574661 | +0.000000000 | +0.000000000 | 1.000000 |
| Orion HEK293T | `sd:pre:eta=0.25` | `sd:pre:eta=0.25` | 0.782332166 | 0.779653408 | -0.002678758 | 0.778935081 | -0.003397085 | +0.001379037 | 0.990099 |
| **macro** | non-unanimous | non-unanimous | **0.794234434** | **0.793985015** | **-0.000249419** | **0.793840352** | **-0.000394082** | — | — |

For K562, the observed positive delta was actually below almost the entire
coordinate-permutation null: null mean **+0.003825536**, p95
**+0.005030529**, empirical upper-tail p **0.990099**. The biological
gene-to-gene alignment therefore did not explain that gain. HEK293T was
negative before the null test.

## Complete fixed-grid comparison

This table evaluates every fixed candidate on the outer folds and is therefore
exploratory, not eligible for selection. `positive` is the number of positive
outer folds out of three.

| candidate | replicate-A macro delta | A positive | replicate-B macro delta | B positive |
| --- | ---: | ---: | ---: | ---: |
| `mad:post:eta=0.5` | +0.001007652 | 3 | +0.000668443 | 2 |
| `mad:post:eta=0.25` | +0.000658466 | 2 | +0.000419024 | 2 |
| `mad:post:eta=1` | +0.000274361 | 2 | -0.000174593 | 2 |
| `mad:pre:eta=0.25` | +0.000064849 | 2 | -0.000104756 | 1 |
| identity | +0.000000000 | 0 | +0.000000000 | 0 |
| `mad:pre:eta=0.5` | -0.000314268 | 1 | -0.000493849 | 1 |
| `sd:post:eta=0.25` | -0.000648489 | 1 | -0.000528768 | 1 |
| `sd:post:eta=0.5` | -0.000997675 | 1 | -0.000818094 | 1 |
| `sd:pre:eta=0.25` | -0.001606257 | 1 | -0.001481548 | 1 |
| `sd:post:eta=1` | -0.003277364 | 1 | -0.002728642 | 1 |
| `sd:pre:eta=0.5` | -0.003446969 | 1 | -0.003092794 | 1 |
| `mad:pre:eta=1` | -0.003935830 | 0 | -0.003995690 | 1 |
| `sd:pre:eta=1` | -0.008784532 | 1 | -0.008091148 | 1 |

The tempting post-hoc result is `mad:post:eta=0.5`, but it is still below the
predeclared +0.002 macro gate, loses HCT116 on replicate B
(-0.000778187), and was not selected by all source-only folds. It is not a
shipping result.

## Control-feature stability and coverage

Replicate stability is measured across genes with positive finite scales in
both disjoint 512-cell samples.

| donor | estimator | valid genes | Pearson of log scales | median absolute log ratio |
| --- | --- | ---: | ---: | ---: |
| K562-GWPS | SD | 7,679 | 0.973010438 | 0.020932287 |
| K562-GWPS | MAD | 2,043 | 0.963936713 | 0.077583790 |
| Orion HCT116 | SD | 11,960 | 0.968451919 | 0.029588476 |
| Orion HCT116 | MAD | 4,574 | 0.961313286 | 0.089499854 |
| Orion HEK293T | SD | 12,960 | 0.974676966 | 0.023114927 |
| Orion HEK293T | MAD | 5,476 | 0.956152714 | 0.077471979 |

The K562-selected MAD ratio was valid on only **1,942 / 17,806** outer-axis
genes in replicate A; 15,864 coordinates therefore remained identity and 468
of the valid ratios hit the upper clip. The HEK293T-selected SD ratio covered
11,438 / 17,942 genes, with 413 low and 317 high clips.

## Leakage audit

- The outer held donor's perturbation matrix is read only after the
  source-to-source candidate has been fixed.
- Inner selection uses only the two source donors' effects and controls.
- Held controls are used only to instantiate the already-selected inference
  transform, matching the challenge setting.
- Every fold-specific whitener declares exactly the two source donors; its
  donor list, axis length, and axis SHA-256 are checked before use.
- All 300 panel-target coordinates are excluded from whitening and retrieval.
- The outer all-candidate grid is labeled optimistic and is not used for the
  ship decision.
- No `cell-eval2`, cloud job, submission, network download, or challenge
  context A/B/C data was used.

## Relation to earlier experiments

- Experiment 010's failed basal-profile donor weighting is not repeated:
  source weights stay equal.
- Experiment 012 weights genes by cross-donor perturbation agreement; this
  experiment instead uses inference-time control-cell variability.
- Experiment 031's variance channel changes emitted single-cell variance for
  DE metrics; this experiment changes effect-vector coordinates and evaluates
  PDS only.
- Experiment 038 supplies the union-axis geometry and fixed gamma 0.65.
- Experiment 041 identifies controls through cross-platform basal means but
  does not transform perturbation vectors.
- Experiment 046 learned additive/PCA context maps and failed its specificity
  control. This experiment tests only a diagonal, control-derived ratio and
  includes a coordinate-permutation null.

## Resource and implementation checks

- Successful Orion extraction: **0.59 GiB** peak process-tree RSS; maximum
  worker self peak **0.529053 GiB**.
- Successful control-stat build: **0.35 GiB** process-tree peak,
  **0.321182 GiB** self peak.
- Successful three-fold evaluation: **0.58 GiB** process-tree peak,
  **0.548752 GiB** self peak.
- Direct bounded reads reproduced prior donor median centered norms within
  0.000024 absolute:
  - K562 6.816785 vs prior 6.816794;
  - HCT116 7.129618 vs prior 7.129641;
  - HEK293T 6.221559 vs prior 6.221577.
- Operator symmetry residuals were at most **2.384186e-05**.
- Unit tests, Ruff, and Python compilation pass.

Resource-development note: an abandoned one-process Parquet prototype was
killed by the guard after a transient sampled process-tree RSS of 0.81 GiB.
No artifact from that attempt was used. The final implementation separates
files into fresh workers and streams both donor arrays and whiteners directly;
all successful final runs stayed at or below 0.59 GiB.

## Artifacts

- `run.py` — extraction, strict holdout selection, evaluation, null, and
  resource-bounded readers.
- `test_variance_transport.py` — synthetic tests for direct NPZ/NPY reads,
  blockwise operators, transport algebra, clipping, selection, and robust
  scales.
- `control_stats.npz` / `control_stats.json` — control variability and
  extraction audit.
- `results.json` — complete folds, candidate grids, ratio metadata, hashes,
  nulls, and decision.
- `per_target_pds.csv` — identity and selected per-target PDS for both control
  replicates.
