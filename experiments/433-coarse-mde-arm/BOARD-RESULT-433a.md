# 433a board result — 0.2536523, NEW BEST (+0.0050669 vs 417e): the MDE law transfers at 1.76× H1

Entry `FCbBcvGxF1HOFPxPISUC`, submitted 2026-09-24 00:49:17Z (scored 01:27Z). Package `d879cf9f…`, model name "the
centrifuge wrote a memoir", description empty. Captured to `data/board-entries/FCbBcvGxF1HOFPxPISUC.json`.

**Arm.** 433a = 417e + `VCC_DE_MDE=0.1,0.1,0.3,coarse`, one variable. The count-balanced MDE delivers A 0.954 /
B 0.975 / C 0.991 of its dose.

| member | raw 433a | raw 417e | Δ raw | Δ scaled | → overall |
|---|---:|---:|---:|---:|---:|
| pds | 0.8589520624303233 | 0.8589520624303233 | **0 (bit-identical)** | 0 | 0 |
| mse | 0.6798392234982936 | 0.6798392234982936 | **0 (bit-identical)** | 0 | 0 |
| **nmae** | **0.919618** | 0.936405 | **−0.016786** | +0.027789 | **+0.004631** |
| fid | 0.514327 | 0.514323 | +0.000004 | +0.000026 | +0.000004 |
| reach | 0.286842 | 0.284445 | +0.002397 | +0.002648 | +0.000441 |
| jac | 0.040808 | 0.040829 | −0.000022 | −0.000061 | −0.000010 |
| **overall** | **0.2536523** | 0.2485853 | | | **+0.0050669** |

**Reading against the preregistration** (`PREREGISTRATION.md`):
- **pds identical: yes, to all 16 digits, and mse too.** The A-channel design (class splits at fixed per-gene totals,
  ballast unwelding) is exact on the server's scorer.
- **nmae −0.0168: the "transfer" bin (≤ −0.008).** The predicted range was −0.009 … −0.022, central −0.012. The
  board read is **1.76× H1's −0.00953**, inside the ×1 … ×2.3 window that pure-magnitude arms have shown (410's weld:
  2.3×).
- fid and jac are flat as predicted.
- reach **+0.0024**, where H1 said −0.004 ± 0.002. It is inside board noise for a single arm, but not the predicted
  sign.
- **Overall +0.0051**, at the top of the predicted +0.0016 … +0.0055.

**What it means.**
- nmae's L1 optimum is conditioned on significance, and a calibrated magnitude floor at the minimum detectable effect,
  sign × min(0.1/√count, 0.1), buys a pure-nmae gain. The same principle holds on the benchmark's own reference.
- The `coarse` delivery was needed to get it into context A (431, 426h).
- **The MDE rides every later arm and is in the October `SHIPPED_ENV` (437).** The 09-25 rule's "MDE transfers"
  branch applies: 433h + 433k.
- zs (433z, H1 −0.0018 on top of coarse) is now read-worthy. It is also a pure-nmae arm.
- Manuscript: `experiments/294-writeup/DRAFT-3i-significance-L1.md` gets its board line.
