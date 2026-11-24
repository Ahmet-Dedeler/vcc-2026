# 075 — Post-whitening reciprocal local query centering

**Date:** 2026-08-23

## Status

**Completed successfully; statistically positive but narrowly below the
predeclared practical promotion threshold. Advance only to independent frozen
confirmation.**

The selected `k=64, beta=0.5` transform improved all three outer donor folds
and its target-cluster bootstrap interval excluded zero. Its macro gain was
`+0.002943142`, only `0.000056858` below the frozen `+0.003` gate. The
machine-readable `promotion_gate.passed` is therefore `false`; this result
does not authorize a production change.

## Predeclared method

Subtract a fraction of each post-whitened query's reciprocal-kNN centroid,
renormalize, and restore the exact original row norm. The transform is
query-only and does not use the hidden reference gallery.

Candidate grid:

- `k`: 8, 16, 32, 64
- `beta`: 0.10, 0.25, 0.50, 0.75
- minimum reciprocal degree before deterministic top-k fallback: 2

Selection and the six-part promotion gate are frozen in `README.md` and
implemented in `run.py` before any outer truth is opened.

## Prior-work audit

The correction is distinct from:

- global component removal and donor-covariance whitening (017–025);
- row/column magnitude sharpening (030);
- hidden-gallery CSLS-density distillation (044 and 073);
- source-covariance subspaces (048 and 052);
- context-control variance transport (067);
- covariance realignment to PerturbMatch (069).

## Result

The deployment candidate was selected entirely from the auxiliary teacher
roles and then frozen before outer VCC-panel truth was opened.

| held donor | baseline PDS | centered PDS | delta |
|---|---:|---:|---:|
| K562-GWPS | 0.795412 | 0.796714 | +0.001302 |
| Orion HCT116 | 0.792493 | 0.796549 | +0.004056 |
| Orion HEK293T | 0.783155 | 0.786627 | +0.003472 |
| **macro** | **0.790353** | **0.793297** | **+0.002943** |

Target-cluster bootstrap:

- 20,000 replicates, seed `7502026`;
- 95% interval `[+0.001048, +0.004899]`;
- probability of positive mean `0.99865`.

Five of six promotion conditions passed:

- nonidentity candidate selected;
- every source-role teacher cell nonnegative;
- source-role nested macro gain at least `+0.001`;
- all three outer donor deltas nonnegative;
- bootstrap lower bound positive;
- **failed only:** outer macro gain `+0.002943 < +0.003`.

The miss was too small to treat as a scientific rejection by itself, but the
threshold was fixed before evaluation and could not be rewritten afterward.
Experiment 079 therefore evaluated the frozen rule once on untouched
KOLF2.1J truth.

That confirmation returned only `+0.0002962724` raw PDS, with bootstrap 95%
interval `[-0.0022893773, +0.0028549881]`. It failed both the positive
lower-bound and `+0.0015` practical-effect gates. **Local centering is now
rejected for production.**

## Execution

- successful run id: `20260823T062845Z`;
- worker: Spot `i-0e05556f09944111f`, `r7i.2xlarge`;
- peak evaluator RSS: `0.302 GiB`;
- exit code: `0`;
- all uploaded artifact hashes verified locally.

Two pre-compute integration failures preceded the successful run:

1. the worker role could not read the shared `output/*` operator prefix;
   the launcher now stages the exact experiment-073 verified operators under
   the run's readable `input/*` prefix;
2. the evaluator expected the obsolete metadata key
   `full_eigendecompositions`; it now requires the current
   `whitening_or_full_covariance_eigendecompositions=0` field.

Neither failed attempt opened held truth or produced a biological result.
No count matrix or challenge submission was created.
