# Experiment 127 results — panel-only metric ceiling

Status: **complete, oracle/reference-composition analysis only**

All surfaces are reused development evidence. None is an independent
confirmation.

## Critical correction and execution state

Experiment 124 was rejected before Experiment 127 was frozen or executed.
No Experiment 127 biological NumPy value array had been loaded, no
`IMPLEMENTATION_FREEZE.json` existed, and no analysis output existed at that
point. The protocol and implementation were amended before the first freeze.

Experiment 127 does not load the Experiment 122 residual, does not use any
Experiment 124 transform, and makes no buildable-transport claim. That work is
deferred to Experiment 128.

## Main result

Across the three reused DE surfaces, only about **2.28%** of reference DE genes
and **2.28%** of absolute reference LFC weight lie on mutable off-diagonal VCC
panel coordinates.

| Reused development surface | VCC target rows | NMAE rows | Reference DE genes on panel | Absolute LFC weight on panel | Perfect-panel incumbent raw NMAE gain |
| --- | ---: | ---: | ---: | ---: | ---: |
| H1 | 5 | 4 | 0.022674 | 0.020528 | 0.025496 |
| Repaired Orion HCT116 | 132 | 127 | 0.023586 | 0.026557 | 0.026890 |
| Repaired Orion HEK293T | 132 | 125 | 0.022240 | 0.021447 | 0.021608 |
| Equal-surface macro | — | — | **0.022833** | **0.022844** | **0.024665** |

The H1 raw-to-scaled NMAE mapping was exactly recoverable from four
already-recorded same-bundle points: slope `-1.6630633627`, absolute raw span
`0.6013000000`, maximum affine-fit residual `2.36e-16`. On its four scored VCC
target rows, perfect independent correction of every mutable panel coordinate
would improve the scaled NMAE member by **0.042402**, equivalent to **0.007067**
in a six-member average if all other members stayed fixed. H1 is a very small
five-target overlap, so this is not a stable cross-context estimate.

The distinction between the two NMAE columns matters:

- `absolute LFC weight on panel` is the zero-prediction-to-perfect-panel
  contribution;
- `perfect-panel incumbent raw NMAE gain` corrects the frozen incumbent and can
  be larger when its panel prediction overshoots or has the wrong sign.

## Exact DE-table ceiling on repaired HCT116

This is the only reused surface with a pre-existing exact prediction DE table.
It belongs to the legacy repaired-HCT116 production prediction
`g0.65_q1.0_a1.08_spread0.25`, not current a14.

| Metric | Incumbent raw | Oracle raw | Maximum raw gain |
| --- | ---: | ---: | ---: |
| Fidelity | 0.540836 | 0.571784 | **0.030948** |
| Jaccard | 0.012031 | 0.012609 | **0.000579** |
| Reach | 0.077017 | 0.102925 | **0.025908** |

These are separate metric-specific optima and must not be summed as one
simultaneously attainable intervention.

- Fidelity allows the oracle to choose any correctly signed, adjudicable panel
  calls while all non-panel calls remain fixed. This can favor calls that are
  not reference-significant and therefore conflicts with Jaccard.
- Jaccard calls exactly the reference-significant mutable panel genes and no
  other panel genes.
- Reach places all correctly signed reference-DE panel genes before the
  incumbent non-panel sequence, whose relative order and signs remain fixed,
  then applies the exact deepest-prefix 90% purity rule. This is an
  ordering-respecting bound, not the invalid shortcut of adding the panel
  fraction directly.

## Interpretation

Panel-only work cannot be the main route to closing the remaining score gap:
roughly 97.7% of reference DE membership and LFC mass is outside the panel on
all three development surfaces.

The panel can still move fidelity and reach by a few raw points because:

- fidelity can benefit from correctly signed panel calls even outside the
  reference-significant set, at a Jaccard cost; and
- early correct panel entries can buffer later non-panel sign errors in
  reach's deepest-prefix calculation.

Jaccard is effectively closed for a panel-only intervention on the repaired
HCT116 table: even perfect panel membership changes it by only `0.000579` raw.

## Pan

Pan was not loaded. Its outcome is already consumed, but it lacks a
scorer-equivalent genome-wide reference significance table and matching
prediction DE table, so it cannot support these metric ceilings without
inventing a gate.

## Verification

- Protocol frozen before Experiment 127 biological value-array access.
- Synthetic tests: **5/5 passed**.
- Ruff: passed.
- Python compilation: passed.
- `cell-eval2`: not run.
- Counts/submission: not generated.
- No new held outcome, DLD1, HCASMC, or challenge truth accessed.
- No commit or push.

Key SHA-256:

| Artifact | SHA-256 |
| --- | --- |
| `preregistered_protocol.json` | `b735425e124c76b0b9f5cf388869af50446f109e89af90034f49df3d0be27616` |
| `IMPLEMENTATION_FREEZE.json` | `b95b935dd40535404e122ad5e0ade632f5d5a2202d7eafa33113de9dbb902e32` |
| `summary.json` | `8a2b6e71297512b1d2edd515d989d3ef28370a972d69d104eea3ae30f077aca1` |
| `per-target.csv` | `c1f7729762dce575dbb47e4609e4031a8cfd0eaf0f6349a426f151f4e424733d` |
| `artifact-manifest.json` | `cb4c0f006d7ea6bdb2d09cf773b3bdceb6eb9db1ee8d8a7c59c739adbde22871` |
