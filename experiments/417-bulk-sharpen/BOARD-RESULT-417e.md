# 417e board result — 0.2485853, NEW BEST (+0.0050519 vs 417a)

Entry `wPP8U70kj4mOIvtA83CO`, submitted 2026-09-24 00:01:26Z (scored 00:48Z; upload and scoring took 47 min, longer
than usual). Package `c9a29cdd…`, model name "the cell line got a lawyer", description empty. Captured to
`data/board-entries/wPP8U70kj4mOIvtA83CO.json`.

**Arm.** 417e = 417a + the bulk stack H1 ×2 + CD4 ×1.5, loaded through `VCC_BULK_LOAD` (dump `bb19f0e5…`). Bulk only;
the DE channel is pinned by ballast.

| member | raw 417e | raw 417a | Δ raw | Δ scaled | → overall |
|---|---:|---:|---:|---:|---:|
| pds | 0.858952 | 0.857544 | +0.001408 | +0.002884 | +0.000481 |
| **mse** | **0.679839** | 0.706077 | **−0.026238** | +0.027611 | **+0.004602** |
| nmae | 0.936405 | 0.936552 | −0.000147 | +0.000246 | +0.000041 |
| fid | 0.514323 | 0.514283 | +0.000040 | +0.000135 | +0.000022 |
| reach | 0.284445 | 0.284949 | −0.000504 | −0.000573 | −0.000095 |
| jac | 0.040829 | 0.040826 | +0.000003 | +0.000008 | +0.000001 |
| **overall** | **0.2485853** | 0.2435334 | | | **+0.0050519** |

**Reading (bins in NEXT-SLOT-DECISION's 09-24 rule; 420's price).**
- **The stack transfers, near the top of its range.** 420 priced it at +0.0037 (k = 1) / +0.0056 (k = 1.407) on the
  PLOLO cancer folds. The board gives +0.0051, i.e. k ≈ 1.3.
  - Bulk *composition* steps transfer at about full PLOLO rate. Sharpening, a bulk *shape* step, transferred at
    k 0.42–0.65 (417a).
  - Almost all of it is mse (−0.026 raw), as 418 predicted: energy-neutral composition, better direction on the strong
    targets that dominate mse's pooled ratio.
- **The DE channel stayed pinned** (|Δ| ≤ 0.0005 raw on every DE member): ballast unwelding holds for LOAD arms too.
  Within the band, so slot 2 went to 433a automatically (00:48:53Z).
- **For October** (438): on October-like panels the CD4 ×1.5 half turns negative and H1 ×2 holds. So the October stack
  default is H1 ×2 / CD4 ×1 with the day-of re-pricer, even though the A/B/C board just paid for both halves.
