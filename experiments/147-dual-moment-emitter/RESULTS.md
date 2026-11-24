# Experiment 147 — dual-moment directional emitter

Status: **terminal NO-GO under the preregistered H1-fast gate; mechanism validated**

Date: **2026-08-31**

## Result

The emitter did exactly what it was designed to do:

- every incumbent per-gene column total was preserved exactly;
- scorer pseudobulk and raw PDS were bit-identical in all four arms;
- arithmetic per-cell CPM moved toward the fixed gamma-0.45 directional
  profile;
- all runs completed serially under the memory guard.

The strongest frozen arm, lambda `1.00`, improved the H1-fast six-member score
from `0.2520684612` to `0.2630367491`, a gain of **`+0.0109682879`**.
However, it missed the preregistered fidelity non-regression threshold by
`0.0002699742` raw, so Experiment 147 itself does not promote a candidate.

| arm | scaled avg delta | raw PDS delta | raw NMAE delta | raw reach delta | raw fidelity delta | raw Jaccard delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| lambda 0.25 | `+0.00113721` | `0` | `-0.00456608` | `+0.01368054` | `-0.00677375` | `-0.00392657` |
| lambda 0.50 | `+0.00417508` | `0` | `-0.00794168` | `+0.02207927` | `-0.00568830` | `-0.00374727` |
| lambda 1.00 | **`+0.01096829`** | **`0`** | **`-0.01184725`** | **`+0.04439623`** | `-0.00126997` | `-0.00367811` |

Lower NMAE is favorable. The MSE plug-in profile remained fixed, while the
prediction-side jackknife moved raw normalized MSE by `+0.00518` at lambda
`1.00`; it remained clamped to zero on the H1-fast scale.

## Interpretation

This is the first count-level experiment in the project to produce a large
directional gain while guaranteeing no PDS loss. It validates the underlying
two-summary mechanism:

- pooled totals can remain fixed for retrieval;
- depth-correlated within-block allocation can improve arithmetic-CPM
  direction and reach.

The exact Experiment 147 promotion rule must remain failed. The lambda-1.00
arm may be carried unchanged into a new, independently preregistered
confirmation because that is a new test of transfer, not a post-outcome
retuning of Experiment 147.

## Integrity

- Source-only gate: passed on 40 targets x 3 lambdas.
- Independent source-only audit: passed after the immutable v2 repair.
- Development arms: four, strictly serial.
- Peak scorer process-tree RSS: `10.19 GiB`.
- Minimum observed system free memory: `68%`.
- Challenge build: none.
- Submission: none.

Primary machine-readable result:
`work/h1-fast-v2/scoring-v1/results.json`.
