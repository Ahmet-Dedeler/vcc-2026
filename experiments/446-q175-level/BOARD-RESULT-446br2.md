# 446br2 BOARD RESULT — 0.2704006 (+0.0003578 vs 446b); the repair floor −2 is partial, so −5 stays

Entry `wE7usCaR6TAlPp16TR0V`, submitted 2026-09-27 00:28Z by `arm0927.sh` slot 2 ("the pipette rack learned origami",
description empty). Captured 01:04Z, rank 6 at capture. Package `446br2-…-zr2.vcc` sha256 `9bbc8d2a…`. Reference 446b
`UDgjOKL9gI1udAexq6ue`. Bins: `449-repair-floor/PREREGISTRATION.md`. Reader: `446-q175-level/read0927.py floor`.

446br2 = 446b with the zero-gene repair floor at −2 instead of −5 (`VCC_DE_ZERO_REPAIR=-2`).

| member | Δ raw vs 446b | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.000007 | +0.000018 | +0.000003 |
| mse | +0.000070 | −0.000073 | −0.000012 |
| nmae | **−0.001413** | +0.002244 | **+0.000374** |
| fid | 0 | 0 | 0 |
| reach | −0.000039 | −0.000042 | −0.000007 |
| jac | 0 | 0 | 0 |
| **overall** | | | **+0.0003578** |

## Against the bins
- **nmae −0.0014 raw is "between bins" (−0.002 < d ≤ −0.001): partial.** The rule keeps −5 unless overall ≥ +0.0006.
  Overall is +0.00036, so **−5 stays**, and floor_ok = False: 463's floor adoption is not launched.
- Calls unmoved (fid, jac exactly 0; reach −0.00004). The bulk is pinned (|d| ≤ 0.0005).
- Predicted overall +0.0006 … +0.0026 and nmae −0.002 … −0.008: the effect is in the right direction but below range.
  The shallower floor helps the gated zero genes a little. Most of them sit at the deep end of the reference too.

## Consequences
- October keeps `VCC_DE_ZERO_REPAIR=-5` (SHIPPED_ENV unchanged).
- The 09-28 rule reads key_ok True and floor_ok False, so its zs rider goes on the LODO base (446bj2z).
