# 433k BOARD RESULT — 0.2677594, NEW BEST (+0.0046475 vs 433h); the levelled q 1.25 → 1.5 step transfers at 438's calibrated rate

Entry `T76JWD3tOajk5IYvzwUW`, submitted 2026-09-25 00:35Z by `arm0925_s2.sh` ("the incubator adopted a goldfish") after
433h's capture. The slot-2 rule passed: 433h vs 433a nmae −0.0247 ≤ +0.001 and |Δpds| 0.0005 ≤ 0.003. Scored and captured
01:02Z, rank 7. Package `433k-…-sharpen150-stackload-ctxlevel-mde010coarse-zr5.vcc` sha256 `f1dc369a…`.
Reference 433h `wMoDAi1mca08FA2Bfqkd` 0.2631120. Bins: `PREREGISTRATION-0925.md`. Reader: `read0925.py`.

| member | Δ raw vs 433h | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | **+0.004623** | +0.010343 | +0.001724 |
| mse | **−0.017898** | +0.018818 | +0.003136 |
| nmae | +0.000815 | −0.001444 | −0.000241 |
| fid | +0.000011 | +0.000028 | +0.000005 |
| reach | +0.000230 | +0.000196 | +0.000033 |
| jac | −0.000021 | −0.000056 | −0.000009 |
| **overall** | | | **+0.0046475** (vs 433a +0.0141072) |

## Against the bins
- **pds +0.0046 raw: "dose step transfers" (≥ +0.004). mse −0.0179 raw: "transfers" (≤ −0.010).**
- DE members pinned inside the ballast's realisation noise. nmae moved +0.0008 (bin ≤ 0.001), so the render change did
  not leak into the carriers.
- Overall +0.0046, inside the predicted +0.002 … +0.006.
- **Transfer calibration:** 438's grid priced the levelled stack q 1.25 → 1.5 at +0.00847 (k = 1) and **+0.00436 at 417a's
  calibrated k 0.65 / k_pds 0.42**. The board gave +0.0046, which is the calibrated number. The same calibration prices
  the next step (q 1.5 → 1.75 levelled) at +0.00648 − 0.00436 = **+0.0021** on A/B/C, and +0.0007 on the `real` and
  `efp` October-like panels.

## Consequences
- 433k is the incumbent. **The 446b gate fired GO** (pds ≥ +0.001 or mse ≤ −0.004 vs 433h).
- 09-26 is branch A:
  - slot 1 = 446b (433k + q 1.75 at matched ρ), `arm0926_s1.sh`;
  - slot 2 = the LODO-key rider on the better bulk (446bj2 if 446b ≥ 433k + 0.001 with pds ≥ −0.0005 and DE pinned,
    else 433kj), `arm0926_s2.sh`.
- October: SHIPPED q moves off 1.25. It is 1.5 now; 1.75 if 446b reads in its transfer bins. The re-level bound issue
  (DEF-RUNBOOK OPEN, 09-24 22:20Z) must be resolved first.
