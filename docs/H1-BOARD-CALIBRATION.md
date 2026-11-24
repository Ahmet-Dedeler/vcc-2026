# H1 → board transfer, per lever class and member (compiled 2026-09-26 10:00Z from the BOARD-RESULT records)

Use this when preregistering an arm priced on the h1-test harness. Raw member deltas; "ratio" = board Δ / H1 Δ (same
sign unless marked). Sources: each lever's `BOARD-RESULT*.md` and its H1 experiment.

| lever (board entry) | class | member | H1 Δ raw | board Δ raw | ratio |
|---|---|---|---:|---:|---|
| 433a coarse MDE (vs 417e) | DE magnitude (pure nmae) | nmae | −0.00953 | −0.0168 | **1.76×** |
| | | reach | −0.004 ± 0.002 | +0.0024 | inverted (noise-sized) |
| 433h zero repair −5 (vs 433a) | DE delivery defect | nmae | ≈ −0.0005 net at H1's gate hit rate | −0.0247 | ≫ 1 (A/B/C gates held more zeros) |
| 410 welded bulk key (vs 395) | bulk magnitude, welded | mse | −0.049 | −0.080 | 1.6× |
| | | nmae (weld cost) | +0.024 | +0.054 | 2.3× |
| | | reach | −0.016 | −0.011 | 0.7× |
| | | pds | +0.0046 | −0.0038 | **inverted** (the sum-to-one render) |
| 410b ballast, same bulk (vs 395) | bulk magnitude, unwelded | mse | −0.0496 | −0.0789 | 1.6× |
| | | pds | +0.0056 | −0.0036 | **inverted** |
| 387 consensus signs through the bulk | bulk signs | pds | ≈ −0.0043 | −0.0098 | 2.3× (harm) |
| 392 CLEAN4 sign prior (vs 319p2) | DE signs (block prior) | fid | +0.0122 | +0.0082 | 0.67× |
| | | nmae | +0.0024 | +0.0041 | 1.7× (worse, same sign) |
| | | reach | +0.0175 | −0.0039 | **inverted** |
| 395 q80 protection (vs 392) | DE ordering | reach | +0.025 (394 contrast) | +0.0100 | 0.4× |
| | | fid / nmae | −0.0004 / +0.0004 (predicted) | −0.0019 / −0.0017 | small, noisy |
| 417c unwhitened carriers (vs 417a) | DE carrier design | fid | +0.0025 | +0.0012 | 0.47× |
| | | reach | +0.0014 | +0.0024 | ~1.7× (noise-sized) |
| | | nmae | −0.0068 | +0.0089 | **inverted** (2nd carrier change mis-signed) |
| 446bja Jurkat programme on A (vs 446b) | DE signs, one context | fid | +0.0032 (a third of H1's full-panel +0.0096) | +0.0002 | **0.07×** |
| | | nmae | −0.0029 | −0.0010 | 0.35× |
| | | reach | −0.014 | −0.0079 | 0.6× |
| 417a / 433k / 446b sharpening q steps | bulk shape | pds, mse | PLOLO folds (cancer) | k_sh 0.42–0.65 of the cancer folds | kolf2's fold predicted 417a |

**Rules that fall out:**
1. **Pure magnitude levers transfer larger than H1, same sign:** MDE 1.76×; mse from a bulk key 1.6×; the weld's
   nmae cost 2.3×. Use a ×1 … ×2.3 band.
2. **Delivery defects scale with the defect's prevalence on A/B/C,** not with H1's. Measure the prevalence on the real
   build (428's audit) before pricing.
3. **DE sign priors transfer weakly on fid: 0.07–0.67×.** Their reach and nmae effects are not priceable on H1: one
   inversion, one attenuation.
   - Mechanism of a prior's reach cost (2026-09-26, 462 §11): small-gate targets whose top-ranked reference gene the
     prior flips from a correct design sign. On H1, 4 targets (n_real 2–17) carry ~100 % of Pan's −0.041.
   - The count of such targets on A/B/C is unknown, which is why H1 cannot price it. Read reach per arm on the board,
     against `458/PREREGISTRATION.md`'s reach contingency.
4. **Carrier-design changes: H1 cannot price nmae** (two of two mis-signed). Price fid/reach only, and treat nmae as
   unknown with a board-sized downside.
5. **Bulk renders: H1 mis-signs pds for per-target scaling** (410, 410b). mse transfers at 1.6×.
6. **Ordering changes (reach) transfer at ~0.4×** (395).
