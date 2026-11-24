# 446bj2 BOARD RESULT — 0.2722766, NEW BEST (+0.0022338 vs 446b); the LODO key transfers

Entry `YvVru3DAF89Jiyr34olR`, submitted 2026-09-27 00:01Z by `arm0927.sh` slot 1 ("the petri dish hired an agent",
description empty). Scored and captured 00:27Z, rank 6. Package `446bj2-…-stacklodo05load-ctxlevelm-mde010coarse-zr5.vcc`
sha256 `c7206c46…`. Reference 446b `UDgjOKL9gI1udAexq6ue` 0.2700428. Bins: `433-coarse-mde-arm/PREREGISTRATION-0926.md`
(the 433j LODO-key section). Reader: `446-q175-level/read0927.py key`.

446bj2 = 446b with the stack dumped under 420/415's leave-one-donor-out bulk key at clip (0.5, 3) (`lodo415_lo0.5_hi3`,
`8c82bef3`) instead of 410's `wshare_lo0.3_hi3`, re-levelled.

| member | Δ raw vs 446b | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.000208 | +0.000511 | +0.000085 |
| mse | **−0.014473** | +0.015221 | **+0.002537** |
| nmae | +0.000060 | −0.000097 | −0.000016 |
| fid | −0.000022 | −0.000077 | −0.000013 |
| reach | −0.001963 | −0.002198 | −0.000366 |
| jac | +0.000016 | +0.000042 | +0.000007 |
| **overall** | | | **+0.0022338** |

Scaled members: pds 0.81523, mse 0.38042, nmae 0.17122, fid 0.00550, reach 0.23287, jac 0.02842.

## Against the bins
- **Overall +0.0022 with pds ≥ 0: "LODO key TRANSFERS"** (bin ≥ +0.001 with pds ≥ 0). **key_ok = True.**
- The prediction was pds +0.001 … +0.004 raw and mse −0.004 … −0.012 raw. The key's gain landed on mse (−0.0145, past
  the top of the range) rather than on pds (+0.0002).
  - The LODO key down-weights targets whose donors disagree, and the floor rises 0.3 → 0.5. Both move the bulk's
    amplitude profile (mse) more than its direction (pds).
  - This fits 492c's finding that a floor of 0.5 is an mse lever with literal pds ~0.
- **DE members pinned:** nmae +0.00006, fid −0.00002, jac +0.00002, reach −0.002 raw, inside the ordering noise.

## Consequences
- **October: `DEF_DAY_LODO=1`** (DEF-RUNBOOK; the day-of fit `scripts/def_lodo_key.py` is ready, rehearsed by 460dr).
- 460dr (the LODO October path, vs 453dr) becomes a 10-01 candidate. `rule1001b.py` appends it when key_ok holds.
- 446bj2 is the new incumbent for A/B/C bulk arms.
  - Later bulk arms (471c, 479's per-target sharpening pair) were designed on 446b's 410-key stack. Their reads stay
    one-variable vs 446b, since bulk members reproduce exactly.
  - Stacking any of them with the LODO key needs its own price. 479's tercile will come from whichever key is shipped.
