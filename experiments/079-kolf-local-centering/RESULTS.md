# 079 — KOLF confirmation of frozen local centering

**Date:** 2026-08-23

## Status

**Completed successfully and rejected.**

The one frozen experiment-075 candidate produced only a negligible positive
delta on untouched KOLF truth, with a confidence interval crossing zero. It
failed two of the three predeclared promotion conditions and must not be added
to production.

Ten synthetic checks pass under a 1 GiB process-tree ceiling. Measured peak
RSS was 0.061 GiB. Ruff, Python compilation, shell syntax, dependency hashes,
and the launcher's fail-closed approval gate also pass.

## Frozen comparison

- Baseline: exact three-donor production-v9 effect query after the exact
  union-axis `gamma=0.65` whitener.
- Candidate: experiment 075's frozen reciprocal local centering at
  `k=64`, `beta=0.5`, minimum reciprocal degree `2`.
- Held truth: KOLF2.1J only.
- Gallery: all 273 KOLF/VCC shared targets.
- Score axis: all 17,806 KOLF/VCC shared non-panel genes.
- Candidate count exposed to KOLF: exactly one.

## Predeclared gate

Promotion is allowed only if:

```text
raw PDS delta > 0
bootstrap lower 95% > 0
raw PDS delta >= 0.0015
```

## Result

| candidate | raw PDS |
|---|---:|
| exact v9 identity | `0.6064291101` |
| reciprocal centering | `0.6067253825` |
| **delta** | **`+0.0002962724`** |

The 20,000-replicate target-cluster bootstrap interval was
`[-0.0022893773, +0.0028549881]`, with probability of a positive mean
`0.5805`.

The result passes only `delta > 0`. It fails both the positive lower-bound
condition and the minimum `+0.0015` practical effect condition. This
independent confirmation closes the local-centering candidate despite the
near-threshold experiment-075 result.

## Execution

- run id: `20260823T064635Z`;
- worker: Spot `i-0e05556f09944111f`, `r7i.2xlarge`;
- peak process-tree RSS: `14.646 GiB`;
- exit code: `0`;
- all uploaded artifact checksums verified locally;
- no count emission, `cell-eval2`, or challenge submission.
