# Experiment 135 — TeloHAEC atlas confirmation

Status: **terminal NO-GO**

The exact Experiment 122 `equal_mean_agreement` panel atlas was evaluated
once on 65 previously unaggregated TeloHAEC panel targets. The comparator was
the exact staged-a15 effect in the same pooled TeloHAEC control frame, with
one production-like `default_rng(0)` library replay across all 300 targets.
The candidate differed only by the exact Experiment 128 integer-budget
transport.

## Result

| endpoint | staged a15 | candidate | delta | frozen gate |
| --- | ---: | ---: | ---: | ---: |
| signed credit | 0.010202 | 0.014268 | **+0.004066** | >= +0.005 |
| top-20% magnitude-ranked sign precision | 0.515000 | 0.509755 | **-0.005245** | >= +0.010 |

The signed-credit delta exceeded both null p95 values:

- target-label p95 `+0.003255`, add-one p `0.025`;
- response-label p95 `+0.002738`, add-one p `0.020`.

But its target-bootstrap 95% interval was
`[-0.000456, +0.008541]`, so the lower-bound gate failed. More importantly,
precision moved in the wrong direction. Its bootstrap interval was
`[-0.020571, +0.010255]`, and null p-values were `0.355` and `0.260`.

Median integer-transport sign retention was `0.890`, below the frozen `0.900`
gate. Every integer/scorer-space invariant passed, including all observed and
null transports; there were no fallback rows.

Guide replication was noisy even after strict QC: median held-out-guide sign
agreement was `0.555`, median per-target minimum agreement was `0.508`, and
median worst held-out-guide normalized error was `2.490`. The 65 targets had a
median of 8 guides and 67 eligible cells.

## Decision

`promotion_passes=false`. Do not integrate, rescale, subset, or otherwise tune
this atlas against the consumed TeloHAEC outcome. Do not build or submit an
a15-plus-panel candidate. Experiment 122's broad signed signal is real, but
the fixed atlas does not improve the directional-precision endpoint and does
not clear its preregistered incremental gain or stability gates.

Only aggregate results and hashes were released. No truth matrix, per-target
outcome, null sample, or bootstrap sample was persisted.

## Integrity

- implementation freeze SHA-256:
  `0f6a218d590853258f1456163089c60581b773fa2ab0a4064eead4ad1426a743`
- result SHA-256:
  `b6fc074904eaf11b6f954a9ba21d9425e1b67c096338d4c8ca9b8412dc0e99b5`
- artifact-manifest SHA-256:
  `0aaf2a96bfd6345c1fb5cae61fb607e8767b6ce012221cd561d28c9c05d33d81`
- all 10 manifested output files verified byte-for-byte;
- output directory and files are read-only;
- no challenge build, submission, `cell-eval2`, commit, or push occurred.
