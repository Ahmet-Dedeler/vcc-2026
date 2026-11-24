# 077 — low-alpha union-whitened MSE crossing

Status: prepared for an isolated cloud run; no result yet.

## Question

The existing HCT116 alpha sweep tested:

- alpha `1.08`, low-noise emission: raw MSE `2.259`;
- alpha `0.70`, low-noise emission: raw MSE `1.362`;
- alpha `0.50`, low-noise emission: raw MSE `1.052`;
- alpha `0.35`, but only with high spread `0.70`: raw MSE `1.030`.

The last comparison does not answer whether low alpha can cross the official
MSE baseline, because experiment 033 later established that extra emission
noise worsens MSE. A quadratic through the low-noise points predicts that
alpha near `0.35–0.40` may reach raw MSE below the official `0.989` baseline.

This experiment tests that missing combination using the production-shaped
union whitener:

- held truth: Orion HCT116, 132 VCC-panel targets;
- donors: K562-GWPS and Orion HEK293T;
- frozen union-axis whitener from experiment 040;
- gamma `0.65`, q `1`, CPM gate `5`;
- exact-total quasi emission, spread `0.25`;
- candidate alphas: `0.35` and `0.40`;
- exact experiment-040 alpha-`1.08` aggregate as comparator.

The two candidates are scored concurrently on an isolated 128 GiB worker.
No local scorer and no challenge submission are used.

## Promotion gate

This is exploratory because the alpha neighborhood was motivated by previous
HCT116 results. A candidate may advance to independent confirmation only if:

1. raw MSE is below `0.989`;
2. raw PDS is no more than `0.003` below the fixed baseline;
3. official-anchor six-metric average improves by at least `+0.003`;
4. both adjacent alpha candidates improve the six-metric average.

Passing does not directly authorize a challenge submission.

## Launch

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
INSTANCE_ID=<running r7i.4xlarge or larger> \
bash experiments/077-low-alpha-union/aws/launch.sh
```

