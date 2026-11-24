# Experiment 149 results

Decision: **TERMINAL NO-GO under the preregistered mechanism gates**

The fixed dual-moment emitter produced a large and statistically resolved
overall gain on held Orion HEK293T, with exact PDS preservation. It did not
produce the specific NMAE/reach pattern required by the frozen protocol, and
raw normalized MSE exceeded its tolerance. The experiment therefore cannot be
retroactively promoted.

## Held-context result

| Metric | Raw candidate minus incumbent | Scaled delta |
| --- | ---: | ---: |
| PDS | `0.000000` | `0.000000` |
| MSE | `+0.024955` | `0.000000` |
| NMAE | `-0.001114` | `+0.001853` |
| Reach | `+0.003694` | `+0.004122` |
| Fidelity | `+0.005881` | `+0.019605` |
| Jaccard | `+0.019215` | `+0.051932` |
| Average | — | **`+0.012919`** |

The 20,000-draw paired target bootstrap was strongly positive:

- estimate: `+0.0129187`;
- 95% interval: **`[+0.0099553, +0.0158571]`**;
- positive-draw fraction: `1.000`.

## Frozen gate failures

- NMAE required raw delta `<= -0.004`; observed `-0.001114`.
- Reach required raw delta `>= +0.010`; observed `+0.003694`.
- MSE required raw delta `<= +0.010`; observed `+0.024955`.

PDS equality, total projected gain, fidelity, Jaccard, and bootstrap gates
passed.

## Interpretation

The original mechanistic hypothesis did not transfer cleanly: H1 gained mostly
through NMAE and reach, while HEK293T gained mostly through fidelity and
Jaccard. The broader score-level hypothesis is nevertheless supported in two
contexts: the same fixed transformation improved projected overall score on
H1-test and held HEK293T while preserving PDS exactly.

That broader hypothesis requires a new, separately preregistered confirmation;
it cannot change this experiment's terminal decision.

## Execution integrity

- Incumbent prediction SHA-256:
  `e837f122edc831a7012a2f8558516b9fa3d786b519ff49c6e8dce8b7e8dcb9b6`
- Candidate prediction SHA-256:
  `fadd11121718c73449978a9c48d283a5577ee180e54007ea6bf921e15e68dd77`
- Result SHA-256:
  `4d70ad27841e4300ae7bd5c890eeb979fecc5a93db6e1498c0a116abf8283203`
- Incumbent peak RSS: `14.689 GiB`
- Candidate peak RSS: `15.286 GiB`
- Both scorer process groups ended with zero processes.
- Challenge build created: no
- Challenge submission created: no
