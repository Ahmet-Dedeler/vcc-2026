# 319p — Preregistration: the partial-condition CD4 rows as a rider on the 319 composition. Frozen 2026-09-09 23:46Z before the build.

**The arm.** 319's exact configuration (`VCC_NORMREL_EXP=-1.5 VCC_NORMREL_MODE=cap`, PM unset, cd4svd5 base) plus
**`VCC_CD4_PARTIAL_G=0.34`** — the flag 254 built and 254's board point validated on the cd4svd5 base
(+0.0013325, entry `PVHrJHweHD5gvOf43ftL`, all through `mse`). One config key away from 319
(`cd4_partial_g` absent → 0.34); the run script stops before prep unless the built meta shows exactly that.
Builder and driver must hash to the 319d3-run builder `1318cd3c…` and the comment-only-delta driver `27edc486…`.
Package `submissions/319p-cap15-partial.vcc`. Intended as the first 2026-09-10 UTC slot.

**Reference for Δ:** 319, `xuDelEik7zGWMXdNu0R4`, **0.21865973944**, same panel and anchor set (r4).

**Why it may not compose.** On 319 every CD4 row's donor share is re-weighted by the cap normrel rule before the
residual is sized; the partial rows enter as their own residual after that. 254's gain was expression accuracy on
the 15 partially carried targets. If the cap rule already suppresses or rescales those rows, the increment shrinks.

**323 pricing.** 323's H1 composition has 3 partial-condition targets of 85 (`work/compose.json`
`target_condition_count_histogram`), so an H1 point estimate rests on three targets and does not gate this arm.
If run, it is recorded as a sign only.

## Bins, fixed now (Δ overall vs 0.21865974)
| bin | reading |
|---|---|
| **Δ ≥ +0.0005** | the rider composes with cap. **319p becomes the best-scored package** and the October composition carries `VCC_CD4_PARTIAL_G=0.34`. Above 319d3's flat band and an order of magnitude above same-call-set reproducibility (4e-5 raw Jaccard). |
| **−0.0005 … +0.0005** | the rider does not compose with the cap rule; 254's gain was specific to the plain-share base. 319 stays; the flag stays off on the cap composition. |
| **Δ ≤ −0.0005** | interference: the re-weighted partial rows hurt. Flag off for good on this family. |

**Secondaries, no verdict:** `mse` (expected to carry it, raw −0.003 … −0.005 as in 254), `pds` (expected flat),
`reach`/`fid`/`nmae` (small positive in 254).

**Stated in advance.** The displayed row becomes this entry's score. If it wins it is the new best; if it loses
it is not corrected by re-submission. Upload takes 12–19 minutes at the observed 2.3–3.7 MiB/s.
