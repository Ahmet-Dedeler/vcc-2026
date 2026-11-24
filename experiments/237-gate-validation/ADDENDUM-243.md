# ADDENDUM to the 237 preregistration — written 2026-09-04 12:35 UTC, ~11.5 h BEFORE the slot, from 243

The armed package (`submissions/229-stack6-g035.vcc`, `7a3c2d3e…`, pid in `submit-armed.log.pid`) and the
three pre-committed readings at **±0.0015 overall** (§4 of the preregistration) are unchanged. This
addendum records what 243 (`experiments/243-transfer-factor/RESULTS.md`) pre-computed so the result
is read on the instrument that was fixed before it existed, not on one chosen after.

## 1. Two predictions, one board point, and they are ~0.0045 apart
| instrument | ARM6 vs `h1union` (0.2118310) | the leg that separates them |
|---|---:|---|
| 243's class rule, direction legs only (`mse`+`pds`+DE-four), ARM6's own amplitude (−1.24%) charged | **+0.0022 … +0.0038** (point +0.0024; coverage model +0.0029; dose class +0.0044) | — |
| 229's six-member composite (same direction legs, plus 216/188's exact `nmae` identity) | **−0.00219** | `nmae` identity leg **−0.004783** |

So Δ ≥ +0.0015 means the identity mis-signed `nmae` a third time; Δ ≤ −0.0015 means the identity is
validated at scale (the 237 readings, verbatim). Between: the magnitude of the `nmae` leg is refuted
but not its sign.

## 2. Read the cosine MEASURED, never the floor (243 §3, §5.2)
The board cosine of a package on disk is a 13-second measurement (`totals243.py` + 197's identity with
the arm's own `A_c`/`B_c`), and the one-sided floor mis-reads ARM6's −1.24% amplitude as 0.30 units of
the transfer factor (floor k 0.64 for a true 0.94). Use `experiments/243-transfer-factor/work/arm6_243.json`:

| board raw `mse` | k realised | reading |
|---:|---:|---|
| ≤ 0.80768 | ≥ 1.55 | above every model |
| 0.80768 – 0.81002 | 1.25 – 1.55 | dose class (the union point was the outlier) |
| 0.81002 – 0.81135 | 1.08 – 1.25 | coverage-linear `k = 0.78 + 0.63·cov` |
| **0.81135 – 0.81392** | **0.75 – 1.08** | **new-source class (point 0.81244; predicted 0.8119–0.8124)** |
| ≥ 0.81392 | < 0.75 | below every model — the union factor itself was optimistic |

The `mse` leg is deterministic to 1e−7 on this board, so the raw `mse` bin is the clean part of the
reading whatever the six-member total does.

## 3. What each outcome does next (pre-committed)
- **Δ ≥ +0.0015 AND raw `mse` in the new-source bin**: the class rule is confirmed on a fourth point and
  the `nmae` identity is retired as a veto. Re-derive the `nmae` leg against four board points before
  any use; re-open 228's per-gene arm (+0.0028…+0.0042 on direction legs, dead only on `nmae`) as the
  next candidate, re-priced on the class rule.
- **Δ ≤ −0.0015**: the identity is validated; 228 and ARM6 stay dead a fortiori; the class rule's
  direction legs are an upper bound, not a price. Tomorrow's second slot is held.
- **Between**: downgrade the identity from a veto to a discount; re-price 228 with the realised
  `nmae` leg scaled by 228's 2.2× larger raw Δ`nmae`.
- In every branch: record the realised `k` from the measured column, and whether it falls in the
  class bin or the coverage bin — that is the fourth point the two-parameter model needs.
