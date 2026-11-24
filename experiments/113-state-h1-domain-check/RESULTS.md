# 113 — STATE H1 domain check

**Verdict: the experiment-106 `state_sm` checkpoint is not a viable
perturbation predictor.**

The truth-free inference freeze used four deterministic 128-cell control
banks and the immutable step-20,000 checkpoint. H1 truth was opened only by
the separate score stage. The strict PDS mask removed the union of the 300
VCC-2026 panel genes and all 85 H1 test perturbation genes.

## Result

| candidate | PDS | cosine | effect NMSE | directional credit |
| --- | ---: | ---: | ---: | ---: |
| frozen incumbent | 0.86345 | 0.03554 | 2.79795 | 0.10528 |
| STATE pooled | **0.50448** | 0.04601 | 19.27952 | 0.12912 |
| STATE bank 0 | 0.49398 | 0.04243 | 23.42708 | 0.13137 |
| STATE bank 1 | 0.50322 | 0.04249 | 17.60136 | 0.13195 |
| STATE bank 2 | 0.50476 | 0.04043 | 18.32726 | 0.09461 |
| STATE bank 3 | 0.51331 | 0.04751 | 21.18527 | 0.11492 |

Pure STATE fails the preregistered `PDS > 0.55` learned-signal threshold.
Its slightly better mean cosine does not encode perturbation identity: PDS is
at chance in every bank and after pooling.

Small row-norm-matched direction blends also fail:

| STATE direction weight | PDS delta | cosine delta | NMSE delta | credit delta |
| ---: | ---: | ---: | ---: | ---: |
| 0.01 | `-0.00266` | `+0.00046` | `-0.00082` | `-0.00106` |
| 0.02 | `-0.00644` | `+0.00093` | `-0.00164` | `-0.00061` |
| 0.05 | `-0.01723` | `+0.00235` | `-0.00417` | `+0.00081` |
| 0.10 | `-0.04762` | `+0.00482` | `-0.00855` | `+0.00262` |

The directional-credit intervals cross zero at every weight, while PDS is
already significantly negative from 2% onward.

## Arbitration implication

A perfect-information selector over row-norm-matched STATE/a14 sign
disagreements has `+0.14504` directional-credit headroom on H1
(`[+0.11098,+0.19340]`). The analogous frozen experiment-106 HEK293T ceiling
is only `+0.00686`, below experiment 110's required `+0.015` on every primary
surface. Experiment 110 therefore stops before training two additional
20,000-step checkpoints.

`work/results-v2.json` is authoritative. `work/results.json` is retained as a
superseded score report: its means were correct, but its bootstrap serializer
did not drop the three targets with no eligible DE coordinates and emitted
null confidence intervals. Inference was not changed.

No challenge artifact or submission was produced.
