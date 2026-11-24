# 410 BOARD RESULT — entry `EETok0iqIbzQFg4CTYMQ` ("the houseplant filed its taxes"), published 2026-09-22 18:51Z

**0.2195152 = −0.0044212 vs 395 (0.2239363).** Package sha 433fdd6c…, prediction a8d97917…, builder 5ff6ef80 (commit ee936a44).
One variable vs 395: `VCC_PERT_SCALE=wshare_lo0.3_hi3` (WELDED per-target bulk magnitude), plus the emitter guard re-pin.

| member | raw 395 → 410 | Δ raw | scaled Δ |
|---|---|---:|---:|
| pds | 0.854738 → 0.850914 | −0.003824 | −0.00891 |
| **mse** | **0.809898 → 0.730213** | **−0.079684** | **+0.08392** |
| nmae | 0.937101 → 0.991476 | +0.054375 | −0.08736 |
| fid | 0.514228 → 0.513441 | −0.000788 | −0.00242 |
| reach | 0.283909 → 0.273245 | −0.010664 | −0.01182 |
| jac | 0.040891 → 0.040911 | +0.000020 | +0.00005 |

## Reading against the preregistration
- **Primary (Δ raw mse ≤ −0.030): TRANSFERS.** −0.0797 raw, above the predicted −0.03…−0.07 range. The per-target key
  carries A/B/C information: at 197's slope (−2.0745 per unit cosine) this is +0.038 of aggregate transfer cosine, i.e.
  0.107 → ~0.146 (+36 %), matching PLOLO's +27…+39 %.
- **The weld costs what H1 said, and more:** nmae +0.054 raw (H1 +0.024), reach −0.011 (H1 −0.016, noise-level there),
  fid ≈ 0, jac 0. The DE channel lost 0.087 + 0.012 + 0.002 scaled = the mse gain.
- **pds −0.0038 raw is NOT what H1 predicted (+0.0046).** Candidate mechanism: the floor targets (final scale 0.19×)
  render a small effect on top of an unscaled rounding/sampling pattern, degrading their retrieval rank on 300 decoys.
  PLOLO effect-level cosine is insensitive to the floor (lo 0.3 vs 0.5: +0.299 vs +0.301), so a higher floor may be free.
## Consequence
410b (ballast, same bulk) is expected at ≈ 395 + 0.0839 − 0.0089 (+ DE ≈ 0) → **≈ +0.0125 overall (~0.236)** if the
unwelding holds on the board as it did on H1.
