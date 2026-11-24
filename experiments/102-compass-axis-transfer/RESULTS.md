# 102 — donor-estimated COMPASS axis transfer

**Primary verdict: NO-GO.**

The preregistered centered reconstruction arm D improved H1 and HCT116, but
lost slightly on HEK293T and did not beat the beta-permutation null. This
experiment does not justify changing the production model.

Before any donor array or held perturbation value was opened, the implementation
was corrected so each surface's held-gallery target rows were excluded from the
global source-axis fit, along with the 300 VCC panel rows. The held target rows
remained available only as per-target source inputs after the geometry was
frozen. Synthetic invariants, Ruff, compilation, and the 24 GiB guarded
three-surface run passed.

## Preregistered primary result

Arm A is the centered donor average. Arm D replaces each donor's axial loading
with the mean centered loading on the donor-estimated target axis.

| surface | PDS A | PDS D | D - A | 95% target bootstrap |
| --- | ---: | ---: | ---: | ---: |
| H1 | 0.815826 | 0.826050 | +0.010224 | [+0.000280, +0.020728] |
| HCT116 | 0.681066 | 0.694772 | +0.013706 | [-0.001561, +0.027644] |
| HEK293T | 0.687890 | 0.687370 | **-0.000520** | [-0.004569, +0.003588] |
| **macro** | — | — | **+0.007803** | [+0.001611, +0.014104] |

The macro beta-permutation-null 95th percentile was **+0.017590**, above the
observed +0.007803. Two preregistered conditions therefore failed:

- D did not beat A on every surface;
- the observed macro gain did not beat the beta-permutation null.

Coverage was identical across arms and every forbidden coordinate remained
exactly zero.

## Post-hoc diagnostic: residual only

Arm C removes each source donor's target-specific loading on that donor's
shared-response axis and retains only the orthogonal residual. It was a
diagnostic arm and was not eligible to rescue the primary result.

| surface | PDS A | PDS C | C - A | post-hoc 95% target bootstrap |
| --- | ---: | ---: | ---: | ---: |
| H1 | 0.815826 | 0.836275 | +0.020448 | [-0.000140, +0.040196] |
| HCT116 | 0.681066 | 0.716111 | +0.035045 | [+0.012144, +0.058294] |
| HEK293T | 0.687890 | 0.693905 | +0.006014 | [-0.004742, +0.016077] |
| **macro** | — | — | **+0.020503** | [+0.009554, +0.031301] |

This is a discovery, not confirmation: C was selected after seeing diagnostic
outcomes, has no preregistered arm-specific null here, and H1/HEK293T intervals
cross zero. It motivates a separate fixed production-shaped experiment on an
independent held surface. Experiment 103 freezes that test before KOLF truth is
opened.

## Interpretation boundary

The official COMPASS implementation estimates the target cell's shared axis
from perturbations measured in that target cell. VCC 2026 supplies no such
target perturbations. Experiment 102 instead averaged source-donor axes, so its
NO-GO does not refute COMPASS and the residual-only signal concerns this
donor-axis surrogate only.

Machine-readable primary results are in `work/summary.json`, with per-surface
metrics and target-level scores under `work/<surface>/`.
