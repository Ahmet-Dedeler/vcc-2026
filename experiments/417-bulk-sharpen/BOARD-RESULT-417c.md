# 417c on the board: unwhitened carriers gain signs (fid, reach) and lose more in nmae. −0.0017 vs 417a

Entry `0xe3j7aqD2bypXKkCTBA`, "the mitochondria unionized" (description empty). Submitted 2026-09-23 00:28:38Z, scored
00:54Z. Package `01291725…`. Builder fc30e8b0.

**One variable vs 417a:** `VCC_CARRIER_GAMMA=0`. The DE carriers come from an unwhitened second pass through the same
frozen CD4/H1 residuals, rows norm-matched to the shipped rows. The bulk is 417a's: identical render energy,
0.960/0.957/0.960 vs 410b.

| member | 417a raw | 417c raw | Δ raw | → overall | H1 prediction (416, 5 seeds) |
|---|---:|---:|---:|---:|---|
| pds | 0.857544 | 0.857878 | +0.000334 | +0.000110 | within noise |
| mse | 0.706077 | 0.706124 | +0.000047 | −0.000008 | within noise |
| **nmae** | 0.936552 | 0.945408 | **+0.008856** | **−0.002896** | **−0.0068 (SIGN INVERTED on the board)** |
| fid | 0.514283 | 0.515457 | +0.001174 | +0.000667 | +0.0025 (right sign, 0.47×) |
| reach | 0.284949 | 0.287319 | +0.002370 | +0.000455 | +0.0014 ± noise (right sign) |
| jac | 0.040826 | 0.040856 | +0.000030 | +0.000013 | +0.00025 |
| **overall** | **0.2435334** | **0.2418746** | | **−0.0016588** | +0.0034 |

## Reading
- **The unwhitened carriers do carry better SIGNS on A/B/C, as on H1.** fid and reach both rose, the two members that
  read signs and ordering.
- **nmae moved the other way, and by more.**
  - nmae reads the arithmetic-mean LFC magnitudes on the reference's DE genes. The unwhitened design moves magnitude
    onto the dominant shared-programme genes.
  - On H1 that helped, because the H1 reference DE sets sit on those programmes. On the three cancer lines it
    overshoots them and under-predicts the line-specific genes.
  - This is the second time H1 has mis-signed nmae for a carrier-design change. 392 was the first.
  - **Rule:** H1 cannot price nmae for carrier-direction or carrier-magnitude changes. Price such arms on fid/reach
    from H1 only, and treat nmae as unknown, with a board-sized downside (here −0.0029).
- **The obvious follow-up keeps the signs and drops the magnitude shift.** Carriers would be sign(γ 0 design) ×
  |γ 0.65 design|, a sign-only arm. The magnitude profile the board already prices stays unchanged, and only the signs
  that fid/reach just paid for move. It needs one builder line in the carrier pass.
  - nmae then moves only through the sign channel, per 232's e = 29–40 for pure sign arms.
  - H1's sign accuracy on the reference DE genes was 57.2 % vs 56.1 %.
- **417a stays the best entry.** 417d (417b + γ 0 carriers) will not be sent: it carries the same nmae loss.
