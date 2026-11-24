# Experiment 150 — score-level dual-moment arbitration on held HCT116

**Decision: `TERMINAL_NO_GO`.**

The fixed dual-moment transform did not transfer to held HCT116. It preserved
every target's PDS exactly and left scaled MSE unchanged, but the projected
six-metric average fell by `0.00504491`.

## Primary result

| Metric | Incumbent scaled | Candidate scaled | Delta |
| --- | ---: | ---: | ---: |
| PDS | 0.638241 | 0.638241 | 0.000000 |
| MSE | 0.000785 | 0.000785 | 0.000000 |
| NMAE | 0.110370 | 0.117299 | +0.006929 |
| Fidelity | 0.201337 | 0.194572 | -0.006765 |
| Reach | 0.032192 | -0.012936 | **-0.045128** |
| Jaccard | -0.040861 | -0.026167 | +0.014694 |
| **Average** | **0.157011** | **0.151966** | **-0.005045** |

Paired target bootstrap, 20,000 draws, seed `15000`:

```text
estimate              -0.0050449084
95% interval          [-0.0099828940, -0.0002777040]
positive draw fraction 0.01935
```

The interval is entirely negative. This is not a borderline gate miss.

## Frozen gates

| Check | Result |
| --- | --- |
| Per-target PDS exactly equal | PASS |
| Projected average delta at least +0.004 | **FAIL** |
| Bootstrap lower bound above zero | **FAIL** |
| At least two nonnegative DE metrics | PASS |
| No DE metric below -0.015 scaled | **FAIL** |
| Scaled MSE does not decrease | PASS |

The failed DE floor is reach at `-0.045128` scaled.

The recorded MSE gate passes by `+2.33e-16`; an independent ratio-of-sums
reduction gives `-5.83e-16`. This is a floating-point-ordering tie at machine
precision and cannot affect the no-go decision.

## Interpretation

Experiments 148 and 149 showed positive total-score movement from the same
emitter, but the carrying metrics differed by context. HCT116 resolves that
heterogeneity against a universal transform: the NMAE and Jaccard gains do not
pay for the reach and fidelity losses.

This result rejects uniform lambda `1.0` as a production change. It does not
invalidate the exact-margin mechanism itself. Experiments 151 and 152 were
preregistered and source-build locked before this outcome existed to test
source-only target-wise and gene-wise confidence restrictions.

## Integrity and process

- Both 132-target arms completed serially under the 24 GiB guard.
- Incumbent peak RSS: `15.225 GiB`; candidate peak RSS: `15.931 GiB`.
- Both process groups reached zero members with return code 0.
- All 132 pseudobulks were exact at construction.
- No challenge prediction or submission was created.

Authoritative result:
`work/scoring-v1/results.json`.
