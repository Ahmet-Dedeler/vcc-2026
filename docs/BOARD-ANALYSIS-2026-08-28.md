# Board analysis, 2026-08-28 — we are winning the wrong metric

Snapshot `data/leaderboard/lb_20260828T1325.json`, 385 teams. We are **rank 16** at
`0.1543607` (down from 7 overnight — every team above us submitted on 08-28; we last
submitted 08-27).

## The board publishes the scaled members directly

This was not being used. Each entry carries `score_pds`, `score_mse`, `score_nmae`,
`score_fid`, `score_reach`, `score_jac` alongside the raw values, and

```
score_avg == mean(the six score_* fields)      max error 3.3e-16 over all 385 entries
```

So there is **no need to re-derive the scaling from anchors** — which matters, because the
anchors in `experiments/089-production-v11/validate_authorization.py` do *not* reproduce
`score_avg` under any clipping rule (max error 1.01 clipped, 358 unclipped). Use the
published `score_*` fields for every board comparison from now on.

## Where we actually stand

| team | overall | pds | mse | nmae | fid | reach | jac |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Aginglab.com *(#1)* | 0.211224 | 0.7089 | **0.1693** | 0.1838 | 0.0030 | 0.2066 | −0.0042 |
| AgingLab_BI *(#2)* | 0.176372 | 0.7038 | 0.0000 | 0.1568 | 0.0001 | 0.1930 | 0.0046 |
| Dewei Hu *(#4)* | 0.170168 | 0.6224 | 0.0000 | 0.1367 | 0.0291 | 0.2317 | 0.0010 |
| **Ahmet Dedeler *(#16)*** | **0.154361** | **0.7366** | **0.0000** | 0.0943 | −0.0287 | 0.1243 | −0.0004 |

**We have the best `pds` on the entire board — 0.7366 against the leader's 0.7089.** We lose
on everything else.

Gap to the leader, per member:

| member | gap (scaled) | gap (overall) | share of the 0.0569 total |
| --- | ---: | ---: | ---: |
| **mse** | **+0.1693** | **+0.02822** | **50%** |
| nmae | +0.0895 | +0.01491 | 26% |
| reach | +0.0822 | +0.01371 | 24% |
| fid | +0.0317 | +0.00528 | 9% |
| jac | −0.0039 | −0.00065 | −1% |
| pds | −0.0276 | −0.00461 | −8% |

## The finding: our expression predictions are worse than a constant

`expr_mse_unbiased_capped_norm` is the predicted profile's squared distance from the measured
one, as a fraction of the measured profile's distance from the control, sampling-noise
corrected. The baseline `b ≈ 0.993` is *the context's mean perturbation response assigned
identically to every perturbation*.

**Our raw is 1.3627. We are 37% worse than predicting one constant vector for all 300
perturbations.** The member is clamped to `[0, 1]`, so we score exactly 0 — as do **376 of
385 teams**. Only 9 teams are above the floor.

This is a real diagnosis, not a scaling artefact: `pds_cosine` is a *cosine* metric, so it
rewards getting each perturbation's **direction** distinguishable while being blind to its
**magnitude**. We have optimised hard against pds and the rank-based DE members, all of which
are magnitude-insensitive, and the one member that reads magnitude says our magnitudes are
wrong.

## We had it, and traded it away

From our own board history:

| date | model | overall | mse raw | mse scaled |
| --- | --- | ---: | ---: | ---: |
| 08-25 | `Alpha 0.40: probing the mse floor` | 0.139963 | **0.9640** | **+0.0251** |
| 08-27 | `Half a step back` | 0.149436 | 1.3236 | 0.0000 |
| 08-27 | `Wildflower` *(current)* | 0.154361 | 1.3627 | 0.0000 |

`Alpha 0.40` cleared the floor. Its `pds` was **higher** than today's (0.7418 vs 0.7366) —
consistent with pds being scale-insensitive — and it paid for the mse in `reach` (+0.0641 vs
+0.1243) and `nmae` (+0.0517 vs +0.0943). But that comparison is confounded: 08-25 predates
the over-dispersion work, so it is *not* a clean amplitude contrast.

## Two structural facts about the mse member worth exploiting

**1. It is clamped at 0, so we are currently paying nothing for being at 1.3627 rather than
0.9940.** Every unit of improvement inside the dead zone is worth exactly zero, which is
precisely why earlier "mse routes" measured as null and were retired. The lever is not
continuous — it is a cliff at `b ≈ 0.993`, and only a jump that clears it pays.

**2. It appears to be clamped per context, then averaged.** `We527` scores `+0.0532` on a raw
of **1.3464** — worse raw than ours, positive score. A single global monotone map cannot do
that; per-context clamping can, and the other eight positive-mse teams fit a global map to
within the residual such averaging would produce. If that reading is right, **a submission
only has to beat the baseline in ONE context to collect a third of the member, with no
penalty in the other two.**

## What this changes

The `reach` work (experiments 094/095, and the built `a10`) is sound and stays. But it is
being pursued in the member worth 24% of our gap while the member worth 50% sits untouched at
a hard zero, and the leader is the only team in the top four collecting it.

**Next: decompose our 1.3627 and find whether the mse cliff can be cleared without giving
back `reach`.** That is experiment 096.
