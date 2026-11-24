# 048 — Source-only low-rank geometry passes

**Date:** 2026-08-22

## Verdict

The salvageable geometry signal from experiment 046 is real and does not
require held-out-target effects.

A strictly source-only, nested LODO selection improved exact macro PDS from
**0.71995 to 0.77695 (+0.05700)**. The gain was positive in all three outer
donors. An independent perturbation-fold seed reproduced it at
**0.77887 (+0.05892)**.

This nearly reproduces experiment 046's target-derived permutation/wrong-map
gains (+0.05656 / +0.05786) without using any outer-donor effect at fit or
hyperparameter-selection time. The mechanism is generic source covariance
geometry, not target-context adaptation.

## Strict protocol

- Outer donors: K562-gwps, Orion-HCT116, Orion-HEK293T.
- Outer evaluation: 259 shared VCC-panel perturbations.
- Fit/selection panel: 600 shared perturbations with zero VCC-panel overlap.
- Score axis: 7,312 genes after excluding all 300 VCC target genes.
- Shaping axis: 5,582 genes after additionally excluding all 600 selection
  target genes.
- For each outer donor, only the other two donors are available.
- Inner validation uses both ordered source-to-source directions and five
  perturbation folds per direction.
- Every PCA basis is fit on the predictor source donor's training folds only.
- Centering and gene scaling are estimated from separate background
  perturbations, excluding both the 600 selection perturbations and every VCC
  target.
- Final source PCs are exact through a sample-space Gram eigendecomposition;
  no gene-by-gene covariance or full eigendecomposition is formed.
- No `cell-eval2`, count emission, submission build, or production artifact.

## Primary nested result

| outer donor | source-only selection | identity | selected | delta | fixed current-style comparator | selected random null |
|---|---|---:|---:|---:|---:|---:|
| K562-gwps | rank 16, gamma 4, ridge 0.1 | 0.74654 | 0.78836 | **+0.04183** | 0.77994 | 0.74598 |
| Orion-HCT116 | rank 128, gamma 1, ridge 0 | 0.70296 | 0.76919 | **+0.06624** | 0.75985 | 0.70326 |
| Orion-HEK293T | rank 64, gamma 2, ridge 0 | 0.71036 | 0.77329 | **+0.06293** | 0.75411 | 0.71102 |
| **macro** | nested | **0.71995** | **0.77695** | **+0.05700** | **0.76463** | **0.72009** |

The fixed current-style comparator is the bounded rank-256, gamma-0.65
source-covariance transform. It improved macro PDS by **+0.04468**. Nested
selection added another **+0.01232** over that comparator.

The selected source transform improved 29/30 inner folds. Its per-validation
donor gains, pooled without using VCC effects, were:

| validation donor | inner PDS delta |
|---|---:|
| K562-gwps | +0.05125 |
| Orion-HCT116 | +0.07165 |
| Orion-HEK293T | +0.04845 |

## Independent split

Seed 148 changed only nearby hyperparameters:

- K562: rank 16, gamma 4, ridge 0 instead of ridge 0.1.
- HCT116: rank 128, gamma 1.25 instead of gamma 1.
- HEK293T: the same rank 64, gamma 2 transform.

Macro PDS was **0.77887 (+0.05892)**, all three outer deltas were positive, and
the pooled source-only selection was positive in **30/30** inner folds.

## Deployable source-only choice

Pooling all six ordered source-to-source directions, while still selecting on
non-VCC perturbations only, picks:

```text
source spectral suppression: rank 64, gamma 2, ridge 0
```

This choice was identical for both split seeds. Its source-only selection gain
was +0.05711 on seed 48 and +0.05615 on seed 148. As a diagnostic, it scores
0.78068 on the three donors' VCC panel (+0.06073), but that number is not the
strict nested estimate because all three donors contributed non-VCC data to
the global hyperparameter choice.

## Controls

- Selected random-subspace macro delta: **+0.00013** on seed 48 and
  **-0.00031** on seed 148.
- Perturbation-label permutation changed PDS by exactly **0.0** in every outer
  fold, as a covariance-only method should.
- Maximum score-matrix difference after row permutation:
  **1.56e-8** on seed 48 and **1.77e-8** on seed 148.
- The result therefore depends on the learned source covariance directions,
  not arbitrary dimensionality reduction or perturbation labels.

## Interpretation

The strongest candidates suppress a relatively small set of high-variance
source-effect PCs. That is consistent with the hubness explanation behind
spectral whitening: shared growth/stress/ribosomal programs dominate cosine
geometry and obscure perturbation identity.

This is a genuine source-only lead, but it is **not yet proof of a production
improvement over the existing full-rank/union-axis whitener**. Experiment 048's
fixed comparator is deliberately low-rank and train-panel-safe, whereas the
production whitener uses a much larger basis and axis. The next valid test is
to compose or replace the full whitener with the rank-64 source suppression on
a clean held-out context. No production file should change from this result
alone.

## Resource and validation record

- Internal peak RSS: **2.919 GiB**
- External process-tree peak: **2.90 GiB**
- Hard external cap: **4 GiB**
- Soft in-process stop: **3.5 GiB**
- System free memory remained **88%**
- Tests: **5 passed**
- Ruff: passed
- Python compilation: passed
- No heavy process remained after either run

## Files

- `README.md` — protocol and commands
- `run.py` — memory-mapped preprocessing, exact source subspace fit, nested
  LODO, comparators, and controls
- `test_source_only_subspace.py` — synthetic correctness/leakage tests
- `results.json` — primary seed 48
- `results_seed148.json` — independent perturbation-fold split
