# Experiment 152 results

Decision: **CLOSE_EXACT_GENE_WISE_CONSENSUS_RAMP_FAMILY** — terminal NO-GO.

The gene-wise consensus ramp gates the dual-moment relaxation *within* each
target instead of between targets. It was preregistered and hash-sealed before
Experiment 150's result existed. It is rejected on its own frozen gates, and
its preregistration closes this family on failure.

## Comparison A — gene-wise candidate minus Experiment 150 incumbent

120 HEK-HCT common targets, frozen HCT order.

| Metric | Scaled delta |
| --- | ---: |
| PDS | `0.000000` |
| NMAE | `+0.007151` |
| Jaccard | `+0.014893` |
| MSE | `+0.000321` |
| Fidelity | `-0.001295` |
| Reach | `-0.040483` |
| **Projected average** | **`-0.003235`** |

Projected overall: incumbent `0.154002`, gene-wise candidate `0.150767`.

20,000-draw paired target bootstrap (seed `15200`): estimate `-0.003235`,
95% interval `[-0.008406, +0.001811]`, `P(delta > 0) = 0.10510`.

| Gate | Required | Observed | |
| --- | ---: | ---: | --- |
| per-target PDS exact | true | true | PASS |
| projected avg gain | `>= +0.006` | `-0.003235` | **FAIL** |
| bootstrap 95% lower | `>= 0.0` | `-0.008406` | **FAIL** |
| improved among nmae/reach/fidelity | `>= 2` | `1` | **FAIL** |
| worst individual DE delta | `>= -0.015` | `-0.040483` | **FAIL** |
| scaled MSE delta | `>= 0.0` | `+0.000321` | PASS |

## Comparison B — versus the best frozen comparator

The preregistered comparator is the higher projected-average of Experiment 150's
uniform lambda-1 candidate and Experiment 151's target-wise selector. On these
120 targets those score `0.150051` and `0.149658`, so the comparator is
**Experiment 150's uniform candidate**.

Estimate `+0.000716`, 95% interval `[-0.000112, +0.001599]`,
`P(delta > 0) = 0.95370`.

| Gate | Required | Observed | |
| --- | ---: | ---: | --- |
| projected avg gain | `>= +0.001` | `+0.000716` | **FAIL** |
| bootstrap 95% lower | `>= -0.001` | `-0.000112` | PASS |

Gene-wise placement genuinely beat uniform placement — `+0.000716` with
`P(delta > 0) = 0.954` — but by less than the preregistered margin, and
Comparison A is the binding constraint regardless. On this target set the
ordering is uniform `0.150051` < gene-wise `0.150767`, with target-wise last at
`0.149658`.

## The decisive result: harm is proportional to dose

Three arms of the same emitter now differ only in how much relaxation they
apply. Mean lambda is computed over all 132 targets.

| exp | placement | mean lambda | delta vs incumbent |
| --- | --- | ---: | ---: |
| 150 | uniform | `1.0000` | `-0.005045` |
| 151 | target-wise, two-level | `0.7727` | `-0.004344` |
| 152 | gene-wise ramp | `0.6170` | `-0.003235` |

Least squares over the three points:

```
delta = -0.004606 * mean_lambda - 0.000539     R^2 = 0.9449
```

Residuals are `[+0.000100, -0.000246, +0.000146]` — an order of magnitude below
the bootstrap intervals. **At zero dose the fit predicts `-0.000539`, which is
indistinguishable from doing nothing.** The relaxation's damage on HCT116 is
a linear function of how much of it you apply, and no forward dose reaches the
`+0.006` gate.

Two independent placement schemes — between targets (151) and within targets
(152) — both failed to beat the dose trend by a meaningful margin. Cross-context
relaxation consensus is not a useful localiser at either granularity.

Reach remains the single failing member (`-0.0405` here, `-0.0436` in 151,
`-0.0451` in 150) while NMAE and Jaccard improve in every arm. The emitter
trades directional depth for lfc calibration; that trade wins on H1 and HEK293T
and loses on HCT116.

**Do not** reopen this family by tuning percentiles, cutoffs, or coefficients on
HCT outcomes — the preregistration forbids exactly that, and the dose fit
explains the whole observed range without needing a new parameter.

## Execution integrity

- Preregistration SHA-256:
  `e560d0acca01d318490f6fc70e0d362222e6f7555f3174924e06db75c21b6f25`
- Scoring authorization SHA-256:
  `b83a6b7ea00d2e3821dd4263b07316faa415c8f17be38125b2657dd9d6580a9e`
- Candidate SHA-256:
  `797f85202737646f9989fdab27b6971cbaf52cef8358b98d3d28941f5acee8c5`
- Candidate results SHA-256:
  `3547d781e6eb743d7ab5569980e8571ae60ebe63917fdc74a2527cb127167d14`
- Only the new candidate arm was scored; Experiment 150's and 151's frozen
  frames were reused. Parallel scorers used: no.
- Peak scorer RSS `16.615 GiB`; minimum system-free memory `61%`;
  elapsed `4593.37 s`; final scorer process count `0`.
- Challenge build created: no. Challenge submission created: no.

## Provenance note

The preregistration and the source build were both frozen before any HCT116
outcome existed. Only the scoring step was authorized afterwards, in
`SCORING_AUTHORIZATION.v1.json`, which binds 17 artifact hashes plus the
scorer's own source. No gate, lambda, seed, or comparator was changed after
Experiment 150 or 151 was read.
