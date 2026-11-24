# 319pA — published result: the same-line calibration rider pays pds on the board; bin 3, closed at this construction

Entry `sEa0LAsWFnXOoXSPJZ59`, "the fern unionized", submitted 2026-09-10 02:14Z (upload 12m49s), published 02:40Z.
Package `submissions/319pA-cap15-partial-ctxA035.vcc`, SHA-256 `fe325f45abfea30666fc63f81ad1739b695b6cc1e6a85e83c8b267aa972d48d9`. Capture
`data/board-entries/sEa0LAsWFnXOoXSPJZ59.json`. Reference 319p `SvzrTOQJAAaNjawzSyrC` 0.2198189778. Bins frozen in
`PREREGISTRATION.md` (commit 4ba9fff3) before the build.

| member | 319p | 319pA | Δ scaled | contribution |
|---|---:|---:|---:|---:|
| pds | +0.788073 | +0.779087 | -0.008986 | -0.001498 |
| mse | +0.185267 | +0.187842 | +0.002574 | +0.000429 |
| nmae | +0.108123 | +0.108049 | -0.000075 | -0.000012 |
| fid | -0.015170 | -0.014104 | +0.001065 | +0.000178 |
| reach | +0.223912 | +0.221647 | -0.002264 | -0.000377 |
| jac | +0.028709 | +0.028718 | +0.000009 | +0.000002 |
| **overall** | **0.2198190** | **0.2185396** | | **−0.0012794** |

**Bin fired: 3, Δ ≤ −0.0005.** Reading, verbatim: *the board's pds price exceeds the frame's prediction; 319p stays;
closed at this construction.*

**What transferred and what did not.** Only context A changed (259 of 300 targets, 4,951 genes, β=0.35). pds fell
−0.0090 scaled panel-wide, i.e. roughly −0.012 raw on A alone — the same order the H1 harness priced (346: −0.009 raw at
β=0.25, −0.020 at 0.5) and NOT what 374's Jurkat source frame predicted (no retrieval loss to β=0.5). mse gained
+0.0026 scaled and fid +0.0011; nmae was flat (H1 predicted −0.015 raw) and reach fell −0.0023 (H1 predicted +0.024).
So the H1 harness priced the discrimination cost correctly and over-predicted the DE-side gain for the real Jurkat
context; the source-frame retrieval proxy under-predicted the cost. Neither instrument is a substitute for the board.

**Consequences.** 319p remains the best package and the incumbent. The same-line calibration rider is closed at this
construction (norm-preserving blend of a shared linear map's output). A target-shared-deflated variant is a different
construction and is not authorised by this reading; it would need its own registration and H1 pricing. The displayed
row is 0.2185396; not to be corrected by re-submission. Production: `VCC_CTX_A_PRIOR` stays unset (byte-inert).
