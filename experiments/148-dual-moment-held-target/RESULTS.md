# Experiment 148 results

Decision: **ADVANCE_TO_CROSS_CONTEXT**

The fixed lambda-1 dual-moment emitter transferred from the 40 H1-fast
development targets to the other 45 H1-test targets without retuning.

## Held 45-target result

| Metric | Raw candidate minus incumbent | Scaled delta |
| --- | ---: | ---: |
| PDS | `0.000000` | `0.000000` |
| NMAE | `-0.019238` | `+0.031994` |
| Reach | `+0.033646` | `+0.037552` |
| Fidelity | `-0.001413` | `-0.004711` |
| Jaccard | `-0.002440` | `-0.006595` |
| MSE | `+0.009678` | `0.000000` |
| Average | — | **`+0.009707`** |

Every frozen held-target gate passed. Per-target PDS was exactly unchanged.

An independent 100,000-draw paired target bootstrap reproduced the point
estimate and gave a 95% interval of **`[+0.001616, +0.017761]`**, with
`P(delta > 0) = 0.99029`.

## Full 85-target result

The projected average improved **`+0.007213`**. PDS remained exactly unchanged.

## Execution integrity

- Candidate SHA-256:
  `a2b8db312bde06752d8a3f5fb6aeaf0037a2cea47c5f9b745f0064aff499b37e`
- Result SHA-256:
  `90863f98958c0498d0b480f664b5ac509fd6829b546d0eab2ee48961b6e757a3`
- Peak scorer RSS: `15.322 GiB`
- Minimum system-free memory: `53%`
- Final scorer process count: `0`
- Challenge build created: no
- Challenge submission created: no

This experiment authorized a separate cross-context confirmation only.
