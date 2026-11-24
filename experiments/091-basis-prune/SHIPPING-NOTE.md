# 091 — shipping note, 2026-08-27

Experiment 091 measured the essential-pruned whitening basis at **+0.002409 avg6
(t = 19.85)** against a size-matched random prune, and **+0.000876 (t = 5.39)** against the
unpruned baseline. It was never shipped because production consumes a *pre-fitted*
operator and the frozen-input guard in `make_submission_pm_cd4.py` pinned exactly one
whitener pair. It ships today.

## What was blocking it

`FROZEN_INPUT_SHA256` pinned `whitener/whitener_union.{npz,json}` to v11's hashes, so any
other operator aborted with `frozen model input hash drift` before a single input was
read. Replaced with `WHITENER_REGISTRY` + `--whitener-id`, which swaps the pair **by name**
and still fails closed on an unregistered id or a one-byte drift. See commit `b24e5e0`.

Operator: `data/donors/whitener_union_pruned.npz`, sha256
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, fitted 2026-08-25 by
`scripts/build_union_whitener.py --basis-prune essential`. 4,393 of 33,787 basis rows
dropped (k562_gwps 2,003 / hct116 1,193 / hek293t 1,197), `dropped_in_panel: []`, and the
2,021-gene common-essential list has **zero overlap** with the 300-gene challenge panel, so
the prune cannot remove signal the submission needs.

## How different is the operator, actually

Not marginally. Eigenspectra of the two covariance operators, both 17,942 × 17,942:

| | λ₁ | **λ₂** | λ₃ | λ₄ | λ₅ | trace |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| shipped (v11) | 6.2567 | **1.1381** | 0.8879 | 0.5054 | 0.4166 | 65.4 |
| pruned | 5.8946 | **0.4966** | 0.3961 | 0.3317 | 0.2988 | 59.4 |

Trace ratio `0.9073`. Median `|Δλ|/λ` across the whole spectrum **0.109**, and **0.197**
over the top 100 components.

**λ₁ barely moves while λ₂ falls 56%.** The common-essential perturbations were
contributing a strong *secondary* shared axis, and `W = (C + λI)^(−γ/2)` was spending its
budget suppressing it. Removing those rows is therefore a real change in what the operator
down-weights, not a rescaling — which is consistent with exp090 losing at t = −13.9 when
unrepresentative rows were *added*.

Expect a board delta that is small but resolvable, in either direction. This is not a
0.1% tweak.

## Caveat that must travel with the number

exp091's `+0.002409` was measured on `data/harness/orion-hct116.h5ad`, which we proved
today is **depth-confounded** (`docs/HARNESS-DEPTH-CONFOUND.md`). Both arms faced the
identical confound and the row counts were matched, so the *ranking* survives; the
magnitude and the per-metric attribution (`nmae +0.0032, reach +0.0017, fid +0.0012`) do
not, because the harness's reference direction distribution is 90.8% negative against ~0.50
on clean data. **This submission is the first clean measurement of the effect.**
