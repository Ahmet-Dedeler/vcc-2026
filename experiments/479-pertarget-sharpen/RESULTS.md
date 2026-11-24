# 479 — per-target bulk sharpening (P2g / S_P2) on 446bj2: builder, known answers, arms (running record)

Preregistration: `PREREGISTRATION.md` (b692c2c5). Builder patch: `patch_builder479.py` (efb8c9dc → **d09a986c**).

## Known answers (chain479, 2026-09-27)
| row | request | result |
|---|---|---|
| tests | builder suite; 479 / 468 / 472 / 451 / 417 unit tests | 56/0/1; all pass |
| KA-S | scalar `VCC_BULK_SHARPEN_Q=1.75`, 446bj2's env | output `711ee43b` = 446bj2, **driver rc 0** |
| KA-T | `tercile:1.75,1.75,1.75` | output `711ee43b` = 446bj2 (terciles 100/100/100, floor_to_low 97), driver rc 1 (below) |
| KA-ST | `stercile:1.75,1.75,1.75` | output `711ee43b` = 446bj2 (terciles 100/100/100, no floor or tie moves), driver rc 1 (below) |

**The KA-T/KA-ST driver exit (found by the cloud review 499, F1).** The builder wrote the right artifact (hash-identical to
446bj2). Then the driver's arm assert (`def_pipeline.py september_rows`) did `float()` on the per-target string and raised
`ValueError` (driver rc 1). The KA runner compares hashes and did not check the rc, so it reported PASS. Any
`run_arm_generic` build with the per-target form would have stopped before its postbuild. The fixed driver is 44714b9a
(cb9e432a), and the 495/496 KA runners now fail on a nonzero rc. The arms below run the full driver path, so they are the
end-to-end proof.

## Arms (chain479b, re-armed behind chain495 at 06:40Z)
P2g (`tercile:2.5,2.0,1.75`) and S_P2 (`stercile:2.5,2.0,1.75`), each built at 446bj2's factors, ntclevel'd, and rebuilt at
factors × the EV-scale (ρ within 0.5 % of 0.9793). -m names from the MODEL-NAMES queue (#1 / #2). Slot: 10-02 (chain1002).
- **P2g round 1** (446bj2's factors), with driver 44714b9a: the per-target form went through build and postbuild end to
  end. Config diff exactly `bulk_sharpen_pt`, `bulk_sharpen_q`; terciles 100/100/100, floor_to_low 97, tie_to_mid 0;
  CLIFF OK 0.95 / 1.02 / 1.07; novelty rc 0. The driver's rc 1 is the expected single postbuild FAIL, the config-drift line.
  ρ A 0.9804, B 0.9720, C 0.9585; EV-scale gives round-2 factors A 1.0957, B 1.1505, C 1.1807.
- **P2g round 2** (A 1.0957, B 1.1505, C 1.1807): config diff exactly `bulk_ctx_scale`, `bulk_sharpen_pt`,
  `bulk_sharpen_q`; CLIFF OK 0.95 / 1.02 / 1.09; novelty rc 0; **package `5ede8f41…`**. ρ **A 0.9795, B 0.9789, C 0.9779**,
  within 0.2 % of 0.9793. READY for 10-02 s1 (-m queue #1).
- **S_P2 round 1** (446bj2's factors): config diff exactly `bulk_sharpen_pt`, `bulk_sharpen_q`; terciles 100/100/100;
  CLIFF OK 0.91 / 0.97 / 1.02. ρ **0.9368 / 0.9302 / 0.9170**: the strength-tercile q moves energy out of the render, as
  492g priced it (energy ×0.94, refunded by the re-level). EV gives round-2 factors A 1.1208, B 1.1760, C 1.2071.
- **S_P2 round 2**: config diff exactly `bulk_ctx_scale`, `bulk_sharpen_pt`, `bulk_sharpen_q`; CLIFF OK 0.95 / 1.02 / 1.09;
  novelty rc 0; **package `7287f6b9…`**. ρ **A 0.9759, B 0.9760, C 0.9757** (−0.35 %, inside the 0.5 % band). READY for
  10-02 s2 (-m queue #2).
- The pair is complete: 10-02 s1 P2g, s2 S_P2 (chain1002 applies at 10-01 20:00Z; reader `read1002.py`).
