# 047 — Direct primary-CD4 donor transfer results

Date: 2026-08-22

## Outcome

Primary-CD4 effects pass the predeclared transfer gate as an additional donor.
Selection used only the held HCT116 and HEK293T folds. The frozen winner was:

- prediction mode: centered mean-absolute gene normalization;
- CD4 source: equal mean of Rest, Stim8hr, and Stim48hr;
- CD4 fraction: `0.50`;
- geometry: no whitening in this experiment.

| evaluation | baseline PDS | selected PDS | delta |
|---|---:|---:|---:|
| HCT116 / HEK293T development mean | — | — | **+0.03085** |
| worst development fold | — | — | **+0.01867** |
| K562 target-role holdout | 0.74643 | 0.76613 | **+0.01970** |
| three-fold macro | — | — | **+0.02714** |

The primary gate passed. CD4 should therefore be tested in the production
union-whitened pipeline rather than rejected as a lineage-specific donor.

## Data and coverage

- Source: Marson 2025 genome-scale primary-CD4 Perturb-seq DE matrix.
- Extracted panel coverage: 291 of 300 VCC targets across the condition union.
- LODO evaluation intersection: 242 targets.
- Evaluation axis: 6,406 non-panel genes shared by all six banks.
- Compact mapped axis: 10,113 official VCC genes.
- Extraction peak RSS: 0.34 GiB.

The compact artifact is preserved under:

```text
data/external/cd4t/cloud-compact/cd4_compact/
```

The downloaded archive SHA-256 is:

```text
ed694e2d6c978d88a2a5065bc4ce2ca075604f5421c165e686702133c4354e61
```

## Important interpretation limit

K562 is a target-role holdout, not a fully independent donor holdout: K562
effects are source inputs while HCT116 and HEK293T are used for development
selection. The positive K562 result is still useful confirmation that the
frozen rule transfers to K562 as a target, but it is weaker evidence than an
outer split that excludes K562 from both fitting and selection.

## Decision

Keep the fixed `0.50` three-condition CD4 blend as the candidate. Do not yet
change the production model: first test it with v9's exact union-whitening
geometry and preserve non-CD4 fallback behavior for the nine missing targets.
