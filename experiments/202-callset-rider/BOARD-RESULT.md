# `rc_callset` board result — **+0.0084448**, `jac` goes POSITIVE, and the identity holds a third time

Entry `31oB9yWy8CSNBDEkY4iW`, model name `Peer Review Rejected This`, **empty description**,
published 2026-09-03 ~18:05 UTC. **Overall `0.1918638367` against `rc_hybrid`'s `0.1834190095` —
`+0.0084448`, rank 11 of 643.** Preregistration: `PREREGISTRATION.md`, committed before the build.

| member | `rc_hybrid` scaled | `rc_callset` scaled | Δ scaled | Δ overall | rc raw | cs raw | Δ raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| `pds` | 0.7458708290005253 | 0.7458708290005253 | **0.000000** | 0.000000 | 0.83712 | 0.83712 | **0.00000** |
| `mse` | 0.12005661519665856 | 0.12005661519665856 | **0.000000** | 0.000000 | 0.87367 | 0.87367 | **0.00000** |
| `nmae` | 0.09429017536041247 | 0.09429042916917090 | +0.0000003 | +0.0000000 | 0.94327 | 0.94326 | −0.00000 |
| `fid` | −0.029758 | −0.029268 | +0.000490 | +0.000082 | 0.50356 | 0.50370 | +0.00015 |
| **`reach`** | 0.170308 | **0.192898** | **+0.022590** | **+0.003765** | 0.22982 | **0.24983** | +0.02001 |
| **`jac`** | −0.000253 | **+0.027335** | **+0.027589** | **+0.004598** | 0.03031 | **0.04043** | +0.01012 |
| | | | **+0.050669** | **+0.008445** | | | |

## Predictions and falsifiers

**Predicted +0.0091 to +0.0111** (192 §5.3, retrodicted level); 202's own emitted h1-test range was
**+0.0050 to +0.0124**. **Realised +0.0084448** — just under the point prediction, inside the
build's own interval.

| falsifier | result |
|---|---|
| F3 — board `pds` and `mse` reproduce `rc_hybrid`'s to ≥12 significant figures | **PASS**, bit-identical at 16 |
| F4 — board `jac` **rises** (88% of the predicted gain, no prior board point in its favour) | **PASS**, 0.03031 → **0.04043** |
| F5 — board `reach` does not fall by more than 0.005 scaled | **PASS** — it **rose** +0.02259 |
| F1/F2/F6 (pre-submit) — panel-total digests 900/900, `n_pred` 2,599 vs 2,600, config diff 4 keys | **PASS** |

## What is new here

**1. `jac` is positive for the first time in the project.** Raw 0.03031 → **0.04043**, scaled
−0.000253 → **+0.027335**. It had been pinned near zero across every arm since v1.

**2. The construction's identity is now board-proven three times.** `pds` and `mse` reproduced to
**sixteen significant figures** for the second consecutive submission, this time while the call set
changed by 79% (jaccard 0.2072 with `rc_hybrid`'s) and `n_pred` fell 4,480 → 2,599. The claim that
the emitted call set and order are free variables at fixed pseudobulk is no longer a construction
argument; it is a measurement, repeated.

**3. Both halves of the union did their job, and the design reason holds.** 192 predicted the
detectability half would buy `jac` and the protected `|H|` head would keep `reach`. `jac` +0.01012
raw and `reach` **+0.02001 raw** — `reach` was predicted to be roughly held and instead improved,
which is the one number that came in above forecast.

**4. Direction is untouched, by construction.** `pds`, `mse` and `nmae` moved by exactly zero. This
arm bought its entire +0.0084 inside the channel `AGENTS.md` §8.2 calls ORDERING.

## Position

**Rank 11 of 643, overall 0.1918638.** Gap to `cqawesome` **+0.02895** (was +0.04383 at 00:30 UTC).
**Both slots today paid: 0.1740865 → 0.1918638, +0.0177773.**
