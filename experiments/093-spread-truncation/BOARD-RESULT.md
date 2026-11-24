# 093 board result — v14, and where the `fid` peak actually is

Entry `gwy86hzx9EKaK8blS8x4`, model name `Stop calling every gene`, submitted
2026-08-26 00:00:44 UTC. **Overall `0.1472728`, against v11's `0.1501827`
(`-0.00291`).** Config: v11 exactly, plus `--tail-spread 1.0 --keep-top 850`.

## Every construction guarantee held, on live data

| metric | v11 | v14 | Δ | predicted |
| --- | ---: | ---: | ---: | --- |
| `pds` | `0.736281` | `0.736281` | **`+0.000000`** | bit-identical |
| `mse` | `0.000000` | `0.000000` | `+0.000000` | clamped either way |
| `nmae` | `0.093404` | `0.093414` | **`+0.000010`** | invariant |
| `reach` | `0.104056` | `0.117787` | `+0.013731` | up |
| `jac` | `0.003053` | `0.007324` | `+0.004271` | `+0.0043` |
| `fid` | `-0.035698` | `-0.071168` | **`-0.035470`** | **`+0.0346` — wrong sign** |

`pds` moved by **exactly zero** and `nmae` by one part in 10^5, on an
independent build against a reference we have never seen. The claim that the
`quasi` emitter's exact-total repair leaves the pseudobulk untouched is now
**board-verified**, not merely asserted. `jac` landed within `0.00003` of its
prediction.

Raw `mse` also improved `1.5801 -> 1.2825` — the jackknife term `C_p` reads
cell spread, so more spread claims more of the noise correction. It stays
clamped at 0 either way, but it confirms the mechanism.

## The one wrong sign is diagnostic, not random

`fid = (k / n_pred) * min(1, n_pred / n_real)`. On the HCT116 harness
`n_pred` stayed ~8x above `n_real`, so the coverage term was pinned at 1 and
`fid` was pure precision, which rose `0.5313 -> 0.5659`. On the challenge
contexts `n_pred` fell **below** `n_real`, coverage engaged, and we took the
undershoot penalty — the asymmetry experiment 092 measured at ~35x.

Inverting: raw `fid = 0.49126` at a truncated-set precision of ~`0.566` implies
**coverage ~= 0.868**, i.e. we call about **13% too few** genes.

**This quantifies the transfer caveat.** The HCT116 reference has 150 cells per
perturbation; the challenge has 400. More cells is more power is a larger
`n_real`, so the harness's median `n_real` of 304 badly understates the
challenge's. Any future harness-tuned call-set parameter must be corrected for
this, and the correction is large.

## Two board points now bracket the peak

| | relative `n_pred` | coverage | precision | raw `fid` |
| --- | ---: | ---: | ---: | ---: |
| v11 | ~5340-equivalent | 1 (saturated) | `0.5019` | `0.5019` |
| **peak** | **= `n_real`** | **1** | **?** | **max** |
| v14 | ~2300-equivalent | `0.868` | `0.566` | `0.49126` |

The optimum is a call set roughly **15-25% larger than v14's**. Because
overshoot is nearly free and undershoot is linear, the target should sit
*above* `n_real`, not on it.

Harness rank-band precision (cumulative, by our own p_adj order) says what is
available there: `0-0.25x n_real` 0.6371, `0.25-0.5x` 0.6131, `0.5-1x` 0.6056,
`1-2x` 0.5915, `2-5x` 0.5559, `5-10x` 0.5216, `>10x` 0.5076. Precision at
`n_pred = n_real` is ~`0.615` on that harness.

**v15 tests `--keep-top 2000` at the same `--tail-spread 1.0`** — a step back
toward v11 that adds ~1,150 protected genes in rank order. Raising `keep_top`
is the right lever rather than lowering `tail_spread`: it adds calls in
descending order of predicted effect, whereas a milder tail lets *random* noise
calls survive.

## What this retires

`--keep-top 850` was experiment 092's counterfactual optimum and experiment
093's harness optimum, and it is **wrong on the board by roughly 30%**. Neither
proxy could see it, because both were measured where the coverage term never
binds. The lesson generalizes: a parameter tuned on a harness whose reference
has one third of the challenge's cells per perturbation is tuned on the wrong
side of a saturating term.
