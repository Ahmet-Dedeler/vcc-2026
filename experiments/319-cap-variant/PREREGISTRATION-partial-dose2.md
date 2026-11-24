# 319p2 — Preregistration: the partial-condition CD4 rows' SECOND dose (G 0.34 → 0.5) on the 319p composition. Frozen 2026-09-10 00:50Z before the build.

**The arm.** 319p's exact configuration plus `VCC_CD4_PARTIAL_G=0.5` (one key: `cd4_partial_g` 0.34 → 0.5). Builder
`1318cd3c…`, driver `27edc486…` (hash-gated). Package `submissions/319p2-cap15-partial50.vcc`. Reference for Δ:
**319p, `SvzrTOQJAAaNjawzSyrC`, 0.2198189778** (same panel, anchors r4).

**Why.** The partial rows composed on two bases at the same size (+0.00133 on cd4svd5, +0.00116 on cap), both
mse-led. The dose axis for this residual has one board point (0.34, mirrored from `VCC_CD4_TARGET_F`). A 1.5× step
tells whether the residual is under-dosed (mse keeps improving), saturated, or already past its optimum (246's
"concentrated" note suggests the 15 rows' information is finite). This is the cheapest information on the table
while the same-line calibration arm is engineered; it is not expected to be large.

## Bins (Δ overall vs 0.2198190)
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | under-dosed at 0.34; 319p2 becomes the incumbent; a third point (0.75) earns a slot only if mse still carries it. |
| **−0.0005 … +0.0005** | saturated; the family closes at 0.34–0.5 with the incumbent staying 319p (ship the lower dose). |
| **Δ ≤ −0.0005** | past the optimum; 0.34 stays; no further dose. |
Secondaries, no verdict: mse (expected carrier), pds (expected flat), reach.

**Stated in advance.** Displayed row becomes this entry's score; not corrected by re-submission. Upload 13–17 min.
This arm may yield its slot to the context-A calibration arm if that one is ready first; whichever is uploaded
second is read against the then-current incumbent.
