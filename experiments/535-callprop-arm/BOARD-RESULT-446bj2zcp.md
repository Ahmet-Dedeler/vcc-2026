# 446bj2zcp BOARD RESULT — 0.2720318 (−0.0007942 vs 446bj2z): NULL on the primary; jac up, fid down, as the unsteered block predicts

Entry `K1Q0y2OFYjTNRqh9jKfk`, submitted 2026-10-01 00:28Z by `464/arm1001.sh` slot 2 ("I see problems and I solve them.",
queue #23, description empty), published ~01:01Z, captured 01:02:05Z, rank 9. Package `fa556a05…` (rebuild == first build
`0b3733d5` byte for byte). Counterpart 446bj2z `be90smZUDcIV7j5MgQ03` 0.2728260. Bins: `PREREGISTRATION.md` (528's, 09-29
16:45Z). Reader: `464/read1001.py cp` (09-30 11:11Z). Read 2026-10-01T01:03Z.

| member | Δ raw vs 446bj2z | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | 0 (identical) | 0 | 0 |
| mse | −2.3e-08 (0.6264081539 → 0.6264081307) | ~0 | ~0 |
| nmae | +0.000099 | +0.000028 | +0.000005 |
| fid | **−0.003580** | −0.011898 | −0.001983 |
| reach | −0.003353 | −0.003734 | −0.000622 |
| jac | **+0.003986** | +0.010839 | +0.001807 |
| **overall** | | | **−0.0007942** |

## Against the bins
- **Primary Δ(fid + jac) = −0.0001765 overall: NULL** (−0.0004 … +0.0004). No October rider for the block alone.
- Falsifiers: none (pds identical; mse −2.3e-8 raw, far below the 1e-6 falsifier: the knob is DE-only on the board to that
  precision; |Δnmae| ≪ 0.003; jac rose). (Wording corrected 2026-10-01T13:45Z, 592: mse is not bit-identical.)
- Predicted Δ(fid + jac) +0.0020 … +0.0025 (H1 at 0.72×) did not hold. H1 at f 1 had fid −0.0013, jac +0.0091, reach
  +0.0096 raw; the board gave fid −0.0036, jac +0.0040, reach −0.0034: jac at 0.44×, fid 2.8× worse, reach of opposite sign.

## Reading
- The mechanism 537 diagnosed (v3 §0 B) is visible: the block's ~550–680 swapped-in genes per context are CALLED (jac up)
  but carry no sign steering (~0.508 sign accuracy on H1), so the extra calls lower fid's directional precision.
- **537 §0 C's rescue condition is met** (NULL with Δjac raw ≥ +0.002 and Δfid < 0): the block + CLEAN4 hook arm
  **446bj2zcph** (549, READY, package `de1e4ab9…`) is the 10-07 candidate: the hook gives the CLEAN4-covered swapped-in genes
  CLEAN4's sign. Its board price must be re-derived from this read (H1 over-predicted the block 1/0.44 on jac).
- Class (AGENTS §11.1c): scorer mechanics (the block is computed from each context's own controls); tries on this switch:
  1 board arm; H1 528 + 537 (forms kept/coupled/decoupled, 3 seeds).
