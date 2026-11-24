# 109 — Gauge-invariant global guide nuisance

Date: **2026-08-30**

## Decision

**NO-GO at the source-only gate.**

The corrected guarded run selected exact identity:

```text
response_rank = 0
ridge = 0
eta = 0
weight_power = 0
huber_c = 2
```

On 102 untouched test targets and 204 guide holdouts:

| quantity | value |
| --- | ---: |
| cosine delta | `0.0` |
| normalized-MSE delta | `0.0` |
| target-cluster bootstrap lower 95% | `0.0` |
| within-target permutation p95 | `0.0` |
| cross-target permutation p95 | `0.0` |
| direction-preserving shrink delta | `0.0` |

Only the nonpositive normalized-MSE gate passed. Every gate requiring positive
source signal failed. No deployment model or singleton correction artifact was
created, held-context evaluation was not authorized, and no count generation
or challenge submission occurred.

## Corrected source extraction

The repaired extractor keeps two estimands separate:

- guide effects use only guide categories with at least 20 cells;
- compact-target reconstruction pools every cell assigned to each target with
  at least 20 total cells.

The complete 1,989,578-cell source pass produced:

| quantity | value |
| --- | ---: |
| genes | 8,248 |
| control cells | 75,328 |
| guide categories | 11,187 |
| retained guide pairs | 10,581 |
| retained targets | 9,787 |
| multi-guide targets | 775 |
| panel-covered targets | 269 |
| panel multi-guide targets | 8 |
| maximum compact reconstruction error | `4.76837158203125e-07` |
| frozen tolerance | `2e-5` |
| peak guarded RSS | `3.904 GiB` |

The source H5AD SHA-256 was
`d3269d9c863d96555eba768b42ab330cc1b57493cd33b4ac8fe7fbab8278f104`.
The immutable rejected bundle is under `work/frozen-source/`; its
`source-run.json` SHA-256 is
`97de09464eb03965c0e147c1e0c46385322457e46e0ede2b72e876df63977aff`.

## Interpretation

Guide-sequence features do explain a tune-set component under some
nonidentity configurations, but none survives the frozen joint
cosine/NMSE selection into the untouched source test. The deployable
gauge-invariant correction therefore has no validated signal. Do not open
held-context truth or revisit this family without a materially different
identifiable mechanism.
