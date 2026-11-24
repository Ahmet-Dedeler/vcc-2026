# 319p — published result: the partial-condition CD4 rows compose with the cap composition; NEW BEST

Entry `SvzrTOQJAAaNjawzSyrC`, "the stapler has tenure", submitted 2026-09-10 00:16:37Z (upload 16m32s at 2.7 MiB/s),
published 00:46Z. Package `submissions/319p-cap15-partial.vcc`, SHA-256
`e8bd9cb5796f2bcd4971b413020338a7a8eb55a64ca08824562d5f5dff29ab4e`. Capture
`data/board-entries/SvzrTOQJAAaNjawzSyrC.json`. Reference 319 `xuDelEik7zGWMXdNu0R4` 0.21865973944.
Bins frozen in `PREREGISTRATION-partial.md` (commit 1af351c1) before the build.

| member | 319 | 319p | Δ scaled | contribution |
|---|---:|---:|---:|---:|
| pds | +0.788314 | +0.788073 | -0.000241 | -0.000040 |
| mse | +0.181013 | +0.185267 | +0.004255 | +0.000709 |
| nmae | +0.107515 | +0.108123 | +0.000609 | +0.000101 |
| fid | -0.015759 | -0.015170 | +0.000589 | +0.000098 |
| reach | +0.222192 | +0.223912 | +0.001720 | +0.000287 |
| jac | +0.028685 | +0.028709 | +0.000024 | +0.000004 |
| **overall** | **0.2186597** | **0.2198190** | | **+0.0011593** |

**Bin fired: 1, Δ ≥ +0.0005 — the rider composes with cap.** Reading, verbatim: *319p becomes the best-scored
package and the October composition carries `VCC_CD4_PARTIAL_G=0.34`.* The increment is the same size and
shape as 254's on the cd4svd5 base (+0.00133, mse-led, pds flat): the mechanism is additive across the two bases.

Secondaries: mse carried it (+0.00425 scaled), reach +0.0017, nmae/fid +0.0006 each, pds −0.0002, jac flat.

**Consequences.** 319p is the incumbent. Production default: `VCC_CD4_PARTIAL_G=0.34` should be switched on in
`def_pipeline.py`'s A/B/C environment when the October driver is next touched (not done tonight; no production
edit). Next information on this axis: the dose (0.34 → 0.5), preregistered separately.
