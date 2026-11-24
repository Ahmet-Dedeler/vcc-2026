# 319d3 — published result: cap −3.0 is flat against cap −1.5; the family closes

Entry `5LOjb2VLiIBh0JG6aAgg`, "a moth in the lint trap", submitted **2026-09-09 22:36:26Z** (upload 12m08s at
3.7 MiB/s), published by 23:02Z, rank 13 of ~840. Package `submissions/319d3-cap30.vcc`, SHA-256
`6ddbf76d2c7b914fae09f606f5da44262f25e415870552f211d606fa08c33d74`. Full capture
`data/board-entries/5LOjb2VLiIBh0JG6aAgg.json`. Reference for Δ: 319, `xuDelEik7zGWMXdNu0R4`, 0.21865973944.

| member | 319 (cap −1.5) | 319d3 (cap −3.0) | Δ scaled | plain's Δ at the same dose step (312→312b) |
|---|---:|---:|---:|---:|
| pds | 0.788313923 | 0.793795 | **+0.005481** | +0.00202 |
| mse | 0.181012812 | 0.181321 | +0.000308 | fell |
| nmae | 0.107514506 | 0.101252 | **−0.006263** | −0.00296 |
| fid | −0.015759119 | −0.015686 | +0.000073 | — |
| reach | 0.222191666 | 0.221041 | −0.001151 | fell |
| jac | 0.028684649 | 0.028624 | −0.000060 | — |
| **overall** | **0.218659739** | **0.218391141** | **−0.000269** | −0.00127 |

**Bin fired: 2, "−0.0005 … +0.0005, flat"** (`PREREGISTRATION-dose2.md`, frozen 2026-09-08 before the entry).
Reading, verbatim from the registration: *the cap removed plain's turnover but adds nothing beyond −1.5. The
family closes at +0.0026; no further dose earns a slot.*

**Secondaries, declared without verdict.** Both answered: cap keeps `pds` rising with dose (more than plain
did), and cap does **not** protect `nmae` at −3.0 (−0.0063, twice plain's loss) although it did at −1.5. The
dose knob trades PDS for NMAE; the cap form neutralised plain's MSE and reach losses but not that trade.
343 (same night) explains why NMAE cannot be bought by amplitude on this composition: sign accuracy on the
reference gate is 0.59.

**The 2×2 is now measured on the board.** plain: −1.5 → −3.0 = −0.00127. cap: −1.5 → −3.0 = −0.00027.
Interaction +0.00100, same sign as 326's local +0.0062 but smaller. Construction-specific response, not a
noise mechanism (COORDINATOR-REVIEW-2026-09-09.md point 4).

**Consequences.** Slot 2 of 2026-09-09 expires unused by the frozen reading; the speculative −4.0 build
(`run319d4.sh`, `PREREGISTRATION-dose3.md` §2) was stopped before its artifact and nothing was packaged.
319 remains the best-scored package. The displayed row is 319d3 (0.2183911); not to be corrected by
re-submission.
