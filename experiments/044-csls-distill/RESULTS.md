# 044 — Donor-trained CSLS distillation does not transfer

**Date:** 2026-08-22

## Verdict

**Do not ship this transform.** A low-rank donor-trained approximation to the
CSLS candidate-density term improved the inner validation task, but failed on
the held-out VCC panel in every outer fold.

| outer held donor | baseline PDS | nested selected | delta | exact-CSLS oracle |
|---|---:|---:|---:|---:|
| K562-gwps | 0.74645 | 0.74136 | **-0.00509** | 0.75064 |
| Orion-HCT116 | 0.70329 | 0.70157 | **-0.00172** | 0.69522 |
| Orion-HEK293T | 0.70983 | 0.69321 | **-0.01661** | 0.68653 |
| **macro** | **0.71985** | **0.71205** | **-0.00781** | **0.71079** |

The result is not a small miss hidden inside selection noise. The inner
selection gains were +0.00756, +0.00545, and +0.00807 for the three outer
folds, respectively; all three reversed sign on the outer VCC perturbations.

## What was tested

Exact CSLS ranks candidate `j` for query `i` by

```text
cos(q_i, t_j) - 0.5 * density(t_j)
```

after dropping terms that are constant within a query row. Exact CSLS is not
submission-valid because `density(t_j)` requires the hidden real-effect
gallery.

The experiment distilled that term into transformations which require only a
submitted query vector:

```text
affine:   q' = q + eta * b
spectral: q' = q + eta * V diag(a) V^T q
combined: both terms
```

The coefficients were ridge-fitted to exact-CSLS density penalties from public
donor pseudo-transfer folds. `V` was limited to rank 8/16/32 and obtained by
randomized range finding, not a full covariance eigendecomposition.

## Leakage controls

- Outer test donors: K562-gwps, Orion-HCT116, and Orion-HEK293T.
- Transform-training donors additionally included RPE1 and Nadig-Jurkat.
- Training used 600 common essential-gene perturbations.
- Outer evaluation used 259 shared 2026 VCC-panel perturbations.
- **Training/evaluation perturbation overlap was exactly zero.**
- Each outer donor was excluded from both fitting and hyperparameter selection.
- Hyperparameters were selected by another leave-one-donor-out loop inside the
  four outer-training donors.
- All target genes from both panels were excluded from the correction basis.

This is stricter than ordinary LODO: it tests transfer across both donor context
and perturbation class.

## Why it failed

The inner essential-gene galleries consistently had positive exact-CSLS
headroom. The unwhitened VCC-panel outer galleries did not:

- K562: exact CSLS +0.00419
- HCT116: exact CSLS -0.00807
- HEK293T: exact CSLS -0.02330

So the density pattern learned from common-essential perturbations is not a
stable property of a gene-effect vector. It changes with the perturbation
gallery and with the upstream geometry. The model learned a real inner-fold
signal and still transferred it in the wrong direction.

The best fixed configuration after inspecting the outer folds was a rank-16
spectral-band transform at only **+0.00178 macro PDS**. That number is
optimistically selected on the same folds and is below the project's practical
resolution threshold; it is not evidence for a usable lever.

## Scope of the negative result

This does **not** refute the earlier exact-CSLS diagnostic after full spectral
whitening. That run showed roughly +0.013 to +0.016 macro PDS of scorer-side
headroom. It does show that:

1. prediction-gallery density is not enough;
2. donor-trained density on disjoint essential perturbations is not enough;
3. a low-rank affine or diagonal spectral correction does not transport that
   headroom to the VCC panel before full whitening.

The next valid test, if revisited, is to serialize source-only **post-whitening
query matrices per donor fold when those operators already exist**, then run
this same nested protocol on those cached matrices. Recomputing full fold
eigensystems locally is intentionally out of scope after the memory incident.

## Resource use

- No `cell-eval2`
- No count generation
- No full eigendecomposition
- One process
- Peak live RSS observed by the guard: **0.89 GiB**
- Experiment-local cache: **211 MiB**

Artifacts:

- `run.py` — reproducible nested protocol
- `results.json` — full fold/config output
- `run.log` — execution log and memory-guard trace
- `cache/` — compact memory-mapped effect matrices
