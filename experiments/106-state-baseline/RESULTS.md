# 106 — official STATE baseline reproduction

**Verdict: NO-GO for `state_sm` at the fixed step-20,000 checkpoint.**

The official support bundle, pinned source, isolated environment, MPS runtime
patch, training run, checkpoint freeze, control-only inference, held
operators, and scoring contracts all completed. The checkpoint hash is
`3f9a1f26138ff1ef1abb04d82640201efdfafb284b92cac88dfef5574505c703`.

## Held-context result

The held screen used the exact 120-target intersection of HCT116 and HEK293T
and removed all 300 VCC panel coordinates from PDS.

| candidate | HCT116 PDS | HEK293T PDS | macro PDS delta | HEK293T credit delta |
| --- | ---: | ---: | ---: | ---: |
| incumbent | 0.80917 | 0.79454 | 0 | 0 |
| STATE fallback | 0.50252 | 0.50028 | **-0.30046** | **-1.90811** |
| direction 10% | 0.79405 | 0.75819 | **-0.02574** | **-0.00603** |
| direction 20% | 0.75350 | 0.70359 | **-0.07349** | **-0.01365** |

All larger direction weights and every effect/composition blend were worse.
The pure STATE effect had median row norm 14.3x the incumbent on HCT116 and
20.2x on HEK293T before direction matching.

## In-domain falsification

Experiment 113 evaluated the same frozen checkpoint on the 85-target H1 test
panel, using only H1 control rows during inference and opening truth only after
the four-bank output was sealed.

| candidate | PDS | cosine | effect NMSE | directional credit |
| --- | ---: | ---: | ---: | ---: |
| incumbent | 0.86345 | 0.03554 | 2.79795 | 0.10528 |
| STATE pooled | **0.50448** | 0.04601 | 19.27952 | 0.12912 |
| direction 1% | 0.86078 | 0.03600 | 2.79713 | 0.10423 |
| direction 10% | 0.81583 | 0.04037 | 2.78940 | 0.10790 |

STATE raises mean cosine but carries almost no perturbation identity: all four
control-bank PDS values lie between 0.4940 and 0.5133. This reproduces the
project's earlier lesson that mean cosine is not a retrieval objective.

## Decision

- Do not build or submit this checkpoint.
- Do not tune a global blend weight against these held outcomes.
- Do not train experiment 110's two additional `state_sm` checkpoints. Its
  HEK293T disagreement-routing oracle is only `+0.00686` directional credit,
  below the preregistered `+0.015` minimum.
- Treat the larger current-source `state` profile as a separate hypothesis
  requiring a new cost/benefit decision; it is not rescued by these results.

No challenge submission was made.
