# Experiment 118 results

Status: **source-only signal confirmed in aggregate, production candidate
NO-GO**

## Source-only stage

Three-condition CD4 leave-one-condition-out validation selected
`unweighted_eta1p00`.

- worst-fold signed credit: `+0.3534167282`;
- worst-fold top-20% weighted sign precision: `0.7250566922`;
- target-label complete-selection null p95: `0`;
- response-label complete-selection null p95: `0`;
- all 300 panel rows have exact zero sum;
- every target diagonal and all 18,233 non-panel coordinates are exact zero.

The source freeze is bound by:

```text
40a0ef4ec5dbbd919f67877c2818f6fc7534ae440547031b149839d97d29252f
```

## One-shot Pan confirmation

The evaluator was repaired and re-frozen before Pan values were loaded. Its
final implementation SHA-256 was:

```text
c3ca0ac8ff7c52e6460fbc28878729d81430eb938e6d1bd65ddd3b5b8afe3664
```

The one-shot run covered all 64 exact Pan/VCC target intersections.

| Statistic | Observed |
| --- | ---: |
| Target-equal signed credit | **+0.0237657838** |
| Signed-credit bootstrap 95% | **[+0.0114574525, +0.0371985758]** |
| Target-label credit null p95 | +0.0071981548 |
| Response-label credit null p95 | +0.0074868151 |
| Top-20% weighted sign precision | **0.3547940640** |
| Precision bootstrap 95% | [0.2772146287, 0.4303335033] |
| Target-label precision null p95 | 0.4037518492 |
| Response-label precision null p95 | 0.3265178210 |

The broad signed-credit signal is real under this test: it is positive with a
positive bootstrap lower bound and exceeds both frozen credit nulls. The
confidence ordering did not transfer: precision missed the `0.55` gate and
did not beat the target-label null.

The direct additive linear-mean application preserved non-panel identity,
target-diagonal identity, panel and total mass, and both normalized PDS-visible
profiles exactly. It failed the required nonnegativity invariant, so the
candidate is not directly emit-able.

## Decision

**NO-GO for the frozen candidate.** Do not tune its dose, clip negative means,
or change its ranking using the revealed Pan result. The result supports only
a new, independently validated family that derives positivity and confidence
rules without Pan-outcome tuning.
