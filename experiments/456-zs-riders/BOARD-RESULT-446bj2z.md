# 446bj2z BOARD RESULT — 0.2728260, NEW BEST (+0.0005493 vs 446bj2); the zs rider TRANSFERS, a pure-nmae read

Entry `be90smZUDcIV7j5MgQ03`, submitted 2026-09-28 12:38Z by `456/arm0928.sh` slot 2 ("the fume hood took up yoga",
description empty), scored and captured 13:05Z. Package `446bj2z-…-zs250.vcc` sha256 `f252b815…`. Base 446bj2
`YvVru3DAF89Jiyr34olR` 0.2722766. Bins: `433-coarse-mde-arm/PREREGISTRATION-0926.md` (the 433z section; chain0928's rule
put the rider on 446bj2 because key_ok held). Reader: `456-zs-riders/read0928.py zs`.

446bj2z = 446bj2 + `VCC_DE_MAG=dn_k250_zs` (429's exact-U down-gene pass on the carriers: every per-gene total, library,
ballast row and rank sum U unchanged; only the arithmetic means of down "mass-holder" genes move).

| member | Δ raw vs 446bj2 | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | 0.000000 (bit) | 0 | 0 |
| mse | 0.000000 (bit) | 0 | 0 |
| **nmae** | **−0.002179** | +0.003296 | **+0.000549** |
| fid | 0.000000 (bit) | 0 | 0 |
| reach | 0.000000 (bit) | 0 | 0 |
| jac | 0.000000 (bit) | 0 | 0 |
| **overall** | | | **+0.0005493** |

Scaled members: pds 0.81523, mse 0.38042, nmae 0.17452, fid 0.00550, reach 0.23287, jac 0.02842.

## Against the bins
- **Δnmae −0.00218 ≤ −0.001: zs TRANSFERS; it rides the October stack.** Prediction −0.0013 … −0.0029 (overall +0.0003 …
  +0.0008): the read sits mid-range.
- **Every other member bit-identical** (the preregistered bug check: none moved by more than 1e-6). The exact-U design is
  confirmed on the board: fid, reach and jac unchanged to the last digit, so the call set, signs and ranks were held.
- Transfer: 1.19× the H1 price (−0.00183), i.e. 1.7× after A/B/C's 0.70× capacity factor. Inside the ×1 … ×2.3 window of
  pure-magnitude levers (`docs/H1-BOARD-CALIBRATION.md`).

## Consequences
- **October: `de_mag = dn_k250_zs`** (`scripts/october_recipe.py decide de_mag dn_k250_zs`). It needs 502's
  `DEF_DAY_DE_MAG` pass-through in `def_day.sh` (prepared, not applied): applied after the current build queue (chain497e
  and chain510 pin def_day e29e793a), with the inert known answer (efp relaxed == ac49d319, rider unset) and a rider-on
  rehearsal, before the 10-06 recipe dr.
- 446bj2z is the new best (0.2728260). For bulk arms the incumbent is unchanged (446bj2's bulk is bit-identical). The
  queued DE-channel arm 496a is one variable vs 446bj2 and stacks exactly on zs on H1 (496z); its 10-05 read stays vs 446bj2.
