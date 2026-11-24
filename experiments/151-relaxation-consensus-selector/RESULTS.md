# Experiment 151 results

Decision: **RETAIN_EXPERIMENT_150_UNIFORM_DECISION** — terminal NO-GO.

The target-wise relaxation-consensus selector did not rescue the dual-moment
emitter on held HCT116. It was preregistered and hash-sealed before Experiment
150's result existed, and it is rejected on its own frozen gates.

## Comparison A — gated candidate minus Experiment 150 incumbent

120 HEK-HCT common targets, frozen HCT order.

| Metric | Scaled delta |
| --- | ---: |
| PDS | `0.000000` |
| NMAE | `+0.006315` |
| Jaccard | `+0.014645` |
| MSE | `+0.000341` |
| Fidelity | `-0.003736` |
| Reach | `-0.043628` |
| **Projected average** | **`-0.004344`** |

Projected overall: incumbent `0.154002`, gated candidate `0.149658`.

A 20,000-draw paired target bootstrap (seed `15100`) gave estimate
`-0.004344`, 95% interval `[-0.009484, +0.000669]`, `P(delta > 0) = 0.04565`.

| Gate | Required | Observed | |
| --- | ---: | ---: | --- |
| per-target PDS exact | true | true | PASS |
| projected avg gain | `>= +0.004` | `-0.004344` | **FAIL** |
| bootstrap 95% lower | `>= 0.0` | `-0.009484` | **FAIL** |
| worst individual DE delta | `>= -0.015` | `-0.043628` | **FAIL** |
| nonnegative DE metrics | `>= 2` | `2` | PASS |
| scaled MSE delta | `>= 0.0` | `+0.000341` | PASS |

## Comparison B — gated candidate minus Experiment 150 uniform lambda-1

Estimate `-0.000393`, 95% interval `[-0.001291, +0.000428]`,
`P(delta > 0) = 0.18285`. Worst individual DE scaled delta `-0.001520`.

| Gate | Required | Observed | |
| --- | ---: | ---: | --- |
| per-target PDS exact | true | true | PASS |
| projected avg gain | `>= +0.001` | `-0.000393` | **FAIL** |
| bootstrap 95% lower | `>= -0.001` | `-0.001291` | **FAIL** |
| worst individual DE delta | `>= -0.01` | `-0.001520` | PASS |

## What this establishes

**Reach carries the damage, and it is diffuse across targets.** Experiment 150's
uniform lambda-1 transform lost `-0.0451` scaled reach on HCT116. Halving lambda
on 60 of 132 targets — the entire low-consensus half — recovered only `+0.0015`
of that, ending at `-0.0436`. If the reach loss were concentrated in targets
where HEK and HCT disagree about the relaxation direction, the gated arm would
have recovered far more.

So cross-context relaxation agreement, measured at the **target** level, does not
localise where the dual-moment emitter hurts. The whole-transform gain seen on
HEK293T (Experiment 149, `+0.012919`) does not survive onto HCT116 in either
uniform or target-gated form.

NMAE and Jaccard again moved the right way (`+0.0063`, `+0.0146`) while reach
and fidelity moved the wrong way — the same signature Experiment 150 produced.
The emitter is trading directional depth for lfc calibration, and on HCT116 that
trade is a net loss.

This does **not** close the gene-wise family: Experiment 152 gates the same
relaxation *within* each target rather than between targets, which is a
different axis and is unaddressed by this result.

## Execution integrity

- Preregistration SHA-256:
  `bc27f9461b8bf7b36bb301d31eb0001e1bb9ecb5e56bdff31848881328a5e2b1`
- Authorization SHA-256:
  `754b71baa772631ec8c43af45a9d7e3edc36b3f7f047dcb49bd3f9bf55ac687b`
- Candidate SHA-256:
  `24f355a7541e1dfa79c8ac0ee8003fe0f6d6eaa9e0376078f6f88c6b11f5c23a`
- Candidate results SHA-256:
  `2f6bfa7d8ebe6b3dca8cc6052daae4b97bc6a76d4135d4a8826c1ed9c780be8c`
- Primary target set: 120 targets, hash
  `ecbfb06f75a6c275cd1f06f4ea10a13cf4512ca279ef507b7bc0aca0fbbfa90a`
- Only the new candidate was scored; Experiment 150's frozen incumbent and
  uniform frames were reused. Parallel scorers used: no.
- Peak scorer RSS: `16.781 GiB`; minimum system-free memory: `66%`;
  elapsed `4508.66 s`; final scorer process count: `0`.
- Challenge build created: no. Challenge submission created: no.
