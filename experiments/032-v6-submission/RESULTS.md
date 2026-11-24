# v6 — the submission, and the joke → config map

**Built:** 2026-08-21 23:47 local · **submitted:** 2026-08-22 00:00 UTC
**Model name on the board:** `Attention is not all you needed`

Model names are public and permanent (AGENTS.md §14). The name is a fake paper title with zero
method content. This file is where the mapping lives.

## Config

| knob | value | why |
|---|---|---|
| donors | K562-gwps : Orion-HCT116 : Orion-HEK293T = **1:1:1** | exp 019 — equal weighting is best *once whitening is on*; it was the worst configuration without it |
| `gene_norm` | off (multi-donor default) | exp 019 — worth +0.0017 under whitening, inside the resolution floor, and the production gene axis differs from LODO's |
| EB shrinkage | off | exp 010/022 — loses at every donor weight, with or without whitening |
| **whitening** | **γ = 0.85**, λ = 0, full 45,531-perturbation donor basis | exp 017/019/024 — the operator, its exponent, and its basis size |
| whitening renorm | **preserve** | exp 025 — global rescaling silently down-weights the 23% of tested genes off the operator axis |
| CPM gate | **5** | exp 023 — beats 20 in six configurations out of six; the scorer tests everything above 5 |
| α | **1.0796** | emitted median effect norm 4.801, identical to v4b, and the count-level optimum (exp 023) |
| emitter | even (minimum-variance) | exp 007; the `quasi` variance channel is measured but not yet shipped (exp 027) |
| β | 0 | exp 005 — the donor mean response is the scored baseline and earns nothing |

Artefacts: `submissions/v6.h5ad` (21.4 GB, 2,674,956,615 stored entries, 56.3% of the density
cap), `submissions/v6.vcc` (936 MB). `vcc prep` peaked at **88.5 GB** — see the log; do not run
it alongside a scoring sweep.

## What it is betting on

One idea, three increments. Cumulative LODO macro PDS **0.7378 → 0.8078, +0.0700** — four
times the entire multi-donor gain and larger than everything else in the notebook combined,
replicated on an independent 4-line harness, with the DE-metric cost measured at count level
(≈ −0.006 overall against a `pds` gain worth +0.017 to +0.023).

Expected **0.115–0.125** against v3's 0.0825 and Vivai's 0.1341. The uncertainty is almost
entirely in how the LODO→board ratio extrapolates: the two calibration points available
(v1 and v3 families) are only 0.010 apart and imply slopes between 0.33 and 0.86, and this
change is six times larger than the gap they span.

## What would falsify it

If v6 lands below ~0.10, whitening did not transfer, and the most likely reason is that the
2026 contexts' effect covariance differs from the donors' more than the held-out donor lines'
do. `submissions/v4b.vcc` (4:1:1, no whitening, expected ~0.106) is built and packaged as the
restore path.
