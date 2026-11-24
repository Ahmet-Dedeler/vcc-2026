# 498p2g BOARD RESULT: 0.2749190, a NEW BEST (+0.0010083 vs 495, +0.0026424 vs 446bj2). P2g HOLDS on the cap-6 base

Entry `JWSghc0o7MN5D6o8Iadc`, rank 15. Package
`submissions/498p2g-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc`
(sha256 `b8039f06…`). 498p2g = the cap-6 dump (495) + 479's per-target key-tercile sharpening (P2g), re-levelled (ρ A/B/C
0.97930 / 0.97930 / 0.98012). Chosen by `rule1004.py` (fixed 09-27): P = P2g, C = TRANSFERS, O = OK.

**Delivery.** `arm1004.sh` slot 1 started at 10-04 00:01:42Z. It reached 90 % (2.73 / 2.77 GB) by 00:59Z, then the hotspot lost
DNS; six resume attempts failed by 01:20:26Z (rc 1, pending upload kept). The coordinator resumed it at 05:21:47Z
(`resume1004s1.sh`, FORCE=1): SUBMIT rc 0 at 05:37:38Z, captured 05:38:10Z. After this, arm1004–1007's retry limit went
6 → 120 (a4a5cc3f + this commit).

**Reading rule** (no bins were preregistered for this pair; the rule was fixed in NEXT-SLOT-DECISION at 05:25Z, before the
score existed): P2g HOLDS unless Δ overall ≤ −0.0005 AND the bulk part < 0.

| member | Δ raw vs 495 | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.001330 | +0.002891 | +0.000482 |
| mse | −0.002697 | +0.002836 | +0.000473 |
| nmae | +0.000056 | −0.000106 | −0.000018 |
| fid | −0.000013 | −0.000041 | −0.000007 |
| reach | +0.000427 | +0.000486 | +0.000081 |
| jac | −0.000006 | −0.000015 | −0.000002 |
| **overall** | | | **+0.0010083** |

- **BIN: P2g HOLDS** → C10b's F = M6+P2g (sent to [10eca4] at 05:3xZ). sharpen_pt stays DECIDED.
- The per-target increment on cap 6 (+0.0010; pds +0.0013, mse −0.0027 raw) matches 479p2g's on cap 3 (+0.0010; pds +0.0016,
  mse −0.0037 raw). The two bulk levers compose additively: 498p2g − 446bj2 = +0.0026 ≈ 495's +0.0016 + P2g's +0.0010.
- DE: all inside 537's 3σ re-roll bands (fid −1.3e-5, jac −6e-6, reach +0.0004 raw).
