# 479 — per-target bulk sharpening on 446bj2: the P2g / S_P2 pair (preregistration, 05:20Z 09-27, before any 479 arm is built)

## The arms (one variable each vs 446bj2, the incumbent: LODO key, q 1.75 scalar)
- **P2g** = 446bj2 with `VCC_BULK_SHARPEN_Q=tercile:2.5,2.0,1.75`: the key tercile, from the loaded LODO dump's
  r_t = ‖bulk_t‖/‖carriers_t‖, with the floor guard. LOW (unreliable rows) is sharpened at q 2.5, MID at 2.0, HIGH keeps
  1.75.
- **S_P2** = 446bj2 with `VCC_BULK_SHARPEN_Q=stercile:2.5,2.0,1.75`: the strength tercile, ‖carriers_t‖. The weakest
  rows get 2.5.
- Each is built, ntclevel'd, and re-levelled to ρ 0.9793 (within 0.5 %). Config diff vs 446bj2: exactly
  `bulk_sharpen_q` (removed), `bulk_sharpen_pt` (added), `bulk_ctx_scale`.
- On 446bj2's dump both assignments are 100/100/100. The key tercile's floor guard puts the 97 LODO-floor targets in LOW.
  The key and strength terciles share 44 of 100 LOW targets.
- Builder `d09a986c` (`patch_builder479.py`). KAs before any arm: scalar q 1.75, tercile at all-1.75 and stercile at
  all-1.75 each rebuild 446bj2 byte for byte (`chain479.sh`).

## Why (third session: 492, 492g, 493, 540c, 546)
- Sharpening extracts the peaks of noisy rows. Low-key / weak rows gain most, while high-key rows are already saturated
  at q 1.75.
- PLOLO on the LODO stack: P2g +0.0035 (3/3 cancer folds, CI [+0.0023, +0.0049]), kolf2 +0.0020; S_P2 +0.0038 (3/3),
  kolf2 +0.0032.
- The permuted control nets +0.0000 (direction −0.5 %), so the targeting is real.
- It ports to October-like panels (9/9 folds) and stacks with the cap-6 key (540c).

## Prediction (vs 446bj2, raw; pds 0.34 × literal PLOLO, mse at k_X 0.61 / 0.12 re-levelled, net of the ballast's nmae)
| arm | Δpds | Δmse | Δnmae | overall |
|---|---:|---:|---:|---:|
| P2g | +0.0026 | −0.0031 / −0.0006 | +0.0010 | **+0.0012 / +0.0008** |
| S_P2 | +0.0025 | −0.0035 / −0.0006 | +0.0006 … +0.0009 | **+0.0013 / +0.0008** |

fid, reach and jac pinned (the carriers are untouched).

## Bins (per arm, vs 446bj2)
- **Transfers:** overall ≥ +0.0008 and Δpds ≥ +0.0015 raw.
- **Partial:** overall +0.0003 … +0.0008, or Δpds +0.0005 … +0.0015.
- **Fails:** below that.
- **Ballast check:** Δnmae +0.0006 … +0.0010 raw expected. Above +0.0015 means the low tercile's ballast is under-priced;
  record it.
- **Leak check:** |Δfid|, |Δjac| > 0.001 or |Δreach| > 0.01 raw means diagnose.
- **The pair:** if both transfer, October takes the one with the larger Δpds. On a tie within 0.0005 raw it takes S_P2
  (tie-free, and the most robust on October-like panels, 493). If only one transfers, it goes. If neither, per-target
  sharpening is closed on the board.

## Slot
10-02 (s1 P2g, s2 S_P2), unless the 10-02 rule, fixed before 10-01's entries exist, says otherwise. Names: fresh jokes.
Descriptions: empty.
