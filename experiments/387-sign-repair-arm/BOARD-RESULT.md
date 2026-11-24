# 387 BOARD RESULT — sign repair to the un-whitened consensus on 319p2: 0.2186049, Δ −0.00151 vs 319p2 (lower bin)

Entry `LRbQ4eOEGLjW7yK82mVa`, submitted 2026-09-13 07:06:50Z, published 07:32Z, rank 19. Package
`submissions/387-cap15-partial50-signcons.vcc` sha256 `986f7665…`, model name "the lint roller has opinions".
Full record `data/board-entries/LRbQ4eOEGLjW7yK82mVa.json`. Reference 319p2 `lqWNXX5YaYArdglu6Stn` 0.2201159.

| member | raw 319p2 → 387 | Δ raw | Δ overall points |
|---|---|---|---|
| pds | 0.854738 → 0.844894 | −0.009844 | **−0.003752** |
| mse | 0.809898 → 0.804997 | −0.004901 | +0.000859 |
| nmae | 0.934784 → 0.934080 | −0.000704 | +0.000188 |
| fid | 0.507933 → 0.508872 | +0.000938 | +0.000522 |
| reach | 0.277798 → 0.281416 | +0.003618 | +0.000667 |
| jac | 0.040903 → 0.040916 | +0.000013 | +0.000006 |

## Reading against PREREGISTRATION.md (bins frozen before the entry; addendum predicted the lower bin)
- **Overall Δ −0.00151 → bin "≤ −0.0005":** the pds cost dominates, as PLOLO predicted; `cons` does not join the bulk.
- **DE gain (nmae + fid + reach) = +0.00138 overall → "0 ± 0.002":** the consensus signs' DE gain through the bulk is
  inside the board's DE noise band — positive on all three members but not a measured gain.
- **pds exchange rate on real lines:** −0.0098 raw, 2.3× the H1 prediction (−0.0042) and well under PLOLO's 6×.
- mse +0.0049 raw in the predicted direction (H1 via the cosine mechanism predicted ≈ +0.002 overall; board +0.0009).

## What it proved
Partial un-whitening through the bulk trades pds for DE at about 2.7 : 1 against on A/B/C. The same class of DE-side
sign information must be delivered without touching the bulk — which is what 392 does (per-cell library channel,
pds/mse pinned by construction). Displayed row fell to 0.2186; not corrected.
