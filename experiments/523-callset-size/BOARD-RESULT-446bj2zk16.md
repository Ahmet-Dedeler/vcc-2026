# 446bj2zk16 BOARD RESULT: 0.2707192 (−0.0021068 vs 446bj2z). K 1,600 HURTS; October K stays 2,200 (rc_callset_cpm off)

Entry `Imj3EXViNco8rT4eLLlE`, rank 26. Submitted 2026-10-06 05:31:20Z by `arm1006.sh` slot 2 (chain1006). It went first because
the slot-1 entry was stuck: 524dr's entry VqBgraiBxlaK0WP5efOF was cancelled while still queued (667). The scoring job waited about
12.5 h for a machine; captured 18:15:00Z. Package `submissions/446bj2zk16-…-zs250-cpm1200.vcc` (sha256 `50d07350…`). Base 446bj2z
`be90smZUDcIV7j5MgQ03` 0.2728260. Arm: 446bj2z + `VCC_RC_CALLSET_CPM=1200` (K 1,600). Bins: `PREREGISTRATION.md` with its 09-30 dose
amendment. Prediction: `experiments/662-callset-regime/PRED-662-K16.md`. Read 2026-10-06 ~18:35Z.

| member | Δ raw vs 446bj2z | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | 0 (bit-identical) | 0 | 0 |
| mse | −4.1e-08 (DE-only jackknife move) | +0.000000 | +0.0000000 |
| nmae | −0.000478 | +0.000927 | +0.0001545 |
| fid | **−0.004065** | −0.013348 | −0.0022247 |
| reach | +0.000353 | +0.000471 | +0.0000786 |
| jac | −0.000174 | −0.000690 | −0.0001151 |
| **overall** | | | **−0.0021068** |

## Against the bins
- **Falsifiers: none.** pds is bit-identical and mse −4.1e-8 raw (≤ 1e-6). jac fell while fid fell, which is not a construction
  error. |Δnmae| 0.00048 ≤ 0.002.
- **Primary Δ(fid + jac) overall = −0.0023398 → bin "≤ −0.0005": the board is closer to H1's regime than either estimate said.
  Keep 2,200.**
- **Dose rule (09-30 amendment).** Both arms are now scored: K20 (10-04) −0.00086 and K16 −0.00234. Neither reaches ≥ +0.0012, so
  nothing is adopted. **October `rc_callset_cpm` → off (K 2,200).** Tries: 2 board arms (K20, K16); H1 2 sizes × 2–3 seeds (537).

## Against the registered predictions (PRED-662-K16)
| instrument | Δfid raw | Δjac raw | Δ(fid + jac) overall | graded |
|---|---:|---:|---:|---|
| I-662 point (band −0.0130 … −0.0107) | −0.0116 | −0.0007 | −0.0067 | **FAILS** (read −0.0041, window [−0.0140, −0.0098]) |
| f 1 reading | −0.0031 | −0.0029 | −0.0030 | just outside on fid (window [−0.0040, −0.0022]; read −0.00406, 0.00006 past the edge, inside the 0.00045 raw read noise); jac far smaller |
| 611 card, f 0.35 (failed on K20 already) | +0.0069 | +0.0010 | +0.0043 | wrong sign again |

- **Grading branch "neither".** The dial puts K16's level between f 0.70 (−0.0053) and f 1.0 (−0.0031), at **f ≈ 0.87** by linear
  interpolation. For call-set SIZE the board behaves close to H1's f 1 regime, as K20 (10-04) indicated.
- I-662's single f\* ≈ 0.61 does not price both size pairs. Its call-set composition reading (cp at f\* 0.61) stays as registered
  until CPH reads (10-07).
- Follow-up from the grading rule: re-fit the level on the three pairs (cp, K20, K16) with `fit662.py` extended by one pair.
- Expected entry ≈ 0.266 (I-662); read 0.2707. The loss was a third of I-662's prediction.

## Consequences
- October K stays 2,200 (`rc_callset_cpm` off), by 523's bins on both board reads.
- 10-07's rule: K reads NEG (k16 < MID), so 10-07 s2 = 459dr. K20 is already scored, so chain1007's dry run (CPH + 459dr) holds.
- CPH (10-07 s1) is the composition read. Its SHIP bar stays as registered.
