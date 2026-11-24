# 446bd2 BOARD RESULT — 0.2704629 (+0.0004200 vs 446b): FLAT, DLD-1 stays out of the bulk and the October stack

Entry `R5uE5Rve3R9nWXGbeSLM`, submitted 2026-09-29 by `457/arm0929.sh` slot 2 ("The task is not so much to see what no one
has yet seen, but to think what nobody has yet thought about that which everybody sees.", description empty). The upload ran
03:08Z → ~10:40Z (a network crawl: 148 MB of 2.77 GB after 3.6 h, then ~3.9 MB/s), published 10:52Z, captured 10:52Z,
rank 7. Package `submissions/446bd2-…-stackdld1u2load-….vcc` sha256 `3ea0fdb6…`. Counterpart 446b `UDgjOKL9gI1udAexq6ue`
0.2700428. Bins: `PREREGISTRATION.md` (2026-09-24 21:40Z). Reader: `457-cd4-half-def/read0929.py dld1`.

| member | Δ raw vs 446b | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | +0.002148 | +0.004910 | +0.000818 |
| mse | −0.000471 | +0.000496 | +0.000083 |
| nmae | +0.000121 | −0.000203 | −0.000034 |
| fid | −0.000020 | −0.000071 | −0.000012 |
| reach | −0.002314 | −0.002614 | −0.000436 |
| jac | +0.000000 | +0.000002 | +0.000000 |
| **overall** | | | **+0.0004200** |

## Against the bins
- **|Δoverall| < 0.0005: FLAT → keep DLD-1 out** (simpler October recipe). The transfer bin also needed mse ≤ −0.001 raw;
  mse moved −0.00047.
- The bulk half landed inside the prediction on pds (+0.0021 raw vs −0.0005 … +0.002; the top of the band) and below it on
  mse (−0.0005 vs −0.001 … −0.004): the DLD-1 residual sharpened direction a little without buying magnitude. Bulk alone
  = +0.0009 overall, at the prediction's centre; reach −0.0023 raw (inside the ballast's 0.01 realisation band: no leak)
  took half of it back.
- DE members inside the leak band (nmae 0.0001, fid 0.00002, jac 0.0000, reach 0.0023 ≤ 0.01): no bulk-to-carrier leak. The
  head's |H| ordering follows the bulk rows, so a composition change re-draws a few head genes; that, not a leak, is the reach
  move.

## Decision
DLD-1 as the second union residual is closed for A/B/C and October (its coverage argument for a disjoint panel, 67 → 73 on
225's pessimistic draw, does not outweigh an A/B/C read at +0.0004). No pass-through is written.
