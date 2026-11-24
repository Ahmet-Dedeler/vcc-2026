# 537 — a preregistered, slot-free test of the DE re-roll noise claim (RESULTS §2) on the board reads of 10-01 … 10-04

Written 2026-09-30 13:45Z, before any 10-01 entry exists. The claim (RESULTS §2.4): on the K-class ballast emission law, two
builds that keep the DE design and differ only in the bulk read DE members that differ by the re-roll alone, with paired sd
**fid 0.000024, jac 0.000013, reach 0.0013 raw** (B1: seven board pairs; H1 scaled to 900 target-rows agrees), i.e.
σ_DE ≈ 0.00025 overall — not 230's 0.0045.

The scheduled reads include several bulk-only arms whose DE design equals their reference's. Their DE deltas are therefore
fresh draws of exactly the quantity §2 estimates, and they cost no slot. Predictions are fixed here.

## The pairs (arm − reference; each is bulk-only with the same DE design and builder seed)
| day | pair | why it is bulk-only |
|---|---|---|
| 10-01 s1 | 460dr − 453dr | the LODO key on the October path (same DE construction) |
| 10-02 s1 | 479p2g − 446bj2 | per-target sharpening (bulk rows only) |
| 10-02 s2 | 479sp2 − 446bj2 | per-target sharpening |
| 10-03 s1 | 495 − 446bj2 | the LODO key at cap 6 |
| 10-03 s2 | 471c − 446b | the H1 union bank swap |
| 10-04 s1 | 498⟨P⟩ − 495 | the per-target form on the cap-6 dump |
| 10-04 s2 | 497dr − 460dr | the October path's percentile cap + S_P2 |
(Any other bulk-only pair of a DE-identical design read by 10-05 is added under the same bands; DE arms such as 446bj2zcp and
496a are NOT noise pairs and are excluded.)

## Predicted bands per pair (3σ of the B1 paired sd; raw members)
- |Δfid| ≤ **0.0001**; |Δjac| ≤ **0.00005**; |Δreach| ≤ **0.0040**.
- nmae carries the ballast term of each bulk change (sharpening raises the render's energy → fewer carrier cells: 433k/446b
  moved nmae +0.0008 / +0.0013 raw); its band is the arm's own preregistered nmae range widened by ±0.0005 and it is reported,
  not scored.

## Reading rules (fixed now)
- **CONFIRMED** if at least 5 of the scored pairs (or all but one, if fewer than 6 are read) are inside all three bands AND
  the rms of the new reach deltas is ≤ 0.0026 raw (B1's 95 % upper bound).
- **REFUTED** if the rms of the new reach deltas is ≥ 0.0050 raw, or at least 3 pairs fall outside the fid or jac band.
- Otherwise **UNRESOLVED**: report the pooled estimate (B1 + the new pairs) with its χ² interval.
- The checker is `noisecheck537.py` (written with this file; it takes the entry ids as arguments and reads only
  `data/board-entries/<id>.json`). Its output goes to RESULTS §2.5 as each day's reads land.
