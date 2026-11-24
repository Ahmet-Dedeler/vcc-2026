# PRED-662-CPH: 446bj2zcph (10-07 s1) under the corrected call-set instrument, at I-662 and at f 1

[21dd4c], 2026-10-04 ~17:55Z. Registered BEFORE 496a's read (10-05), K16's read (10-06) and CPH's own read (10-07).
- **Coordinator's ask:** a corrected CPH prediction at f 1, so that the read grades the corrected instrument.
- **Unchanged:** the arm, the statistics and the bins. They are 537's `PREREGISTRATION-draft-CPH.md` as frozen in 653's
  rule v4:
  - primary = Δ(fid + nmae) overall vs 446bj2zcp, with POSITIVE ≥ +0.0004 / PARTIAL +0.0001 … +0.0004 / NULL < +0.0001;
  - SHIP = Δ(fid + jac) vs 446bj2z ≥ +0.0012 AND overall vs 446bj2z > 0.
- **Added:** the corrected instrument's prediction, beside 537's f 0.35 one. Numbers come from `fit662.py` →
  `work/fit662.json`.

## Construction
- **CPH − 446bj2zcp on the board** = the hook-on-block increment. It is priced from 537's H1 runs on the cp base (3 seeds)
  from three sources:
  - **P8** (CLEAN4 attenuated to A/B/C-like quality, r 0.17): the central estimate;
  - **P1** (H1's full CLEAN4 information): the high case;
  - **P9** (CLEAN4 signs shuffled, zero information): the floor.
- **CPH − 446bj2z** = cp's MEASURED board read plus that increment. The cp read is K1Q0y2OFYjTNRqh9jKfk − be90smZUDcIV7j5MgQ03:
  - fid −0.00358, jac +0.00399, reach −0.00335, nmae +0.00010 raw;
  - −0.00079 overall, of which Δ(fid + jac) is −0.00018.
- **Multipliers:**

| | level | fid | jac | nmae: T before 496a reads | T after 496a reads (10-05) |
|---|---|---|---|---|---|
| **I-662 (instrument of record)** | f\* 0.61 | × 2.69 (one multiplier for both call-set pairs) | × 0.43 | 0.92 … 1.37, mid 1.12 (delivery pairs at f\*, h1-test + 618) | 496a's board Δnmae raw / **(−0.001979)** = 496a's own H1 pair at f\* (`t496_662.json`) |
| f 1 (as asked) | f 1 | × 2.83 (cp's; K20's at f 1 is 0.17, so f 1 has no single call-set multiplier) | × 0.44 | 1.14 … 1.74, mid 1.41 | 496a Δnmae / **(−0.001458)** (537's own denominator) |

## Registered prediction (raw Δ; overall units where stated)
| instrument / source | **primary: Δ(fid + nmae) overall vs 446bj2zcp** | Δfid vs cp | Δnmae vs cp | Δjac vs cp | **SHIP stat: Δ(fid + jac) overall vs 446bj2z** | overall vs 446bj2z (increment's reach excluded) |
|---|---:|---:|---:|---:|---:|---:|
| **I-662, P8 (central)** | **+0.0009** | +0.0010 | −0.0011 | +0.0001 | **+0.0004** | +0.0001 |
| I-662, P1 | +0.0019 | +0.0021 | −0.0027 | +0.0001 | +0.0010 | +0.0012 |
| I-662, P9 | −0.0004 | +0.0000 | +0.0015 | +0.0001 | −0.0001 | −0.0011 |
| I-662, full range (f 0.60–0.62 × T band × P9…P1) | −0.0007 … +0.0021 (P8 alone: +0.0007 … +0.0010) | | | | −0.0003 … +0.0010 | −0.0014 … +0.0013 |
| f 1, P8 | +0.0008 | +0.0006 | −0.0015 | +0.0001 | +0.0002 | +0.0000 |
| f 1, P1 / P9 | +0.0017 / −0.0007 | | | | +0.0008 / −0.0002 | +0.0009 / −0.0014 |
| 537's 10-01 re-price (f 0.35) | +0.0004 … +0.0023 (P8 × 1.95, T 1: +0.0010) | | | | −0.0001 … +0.0009 | |

- **Bins:**
  - POSITIVE on the primary at the central source under both instruments.
  - NOT SHIP under every source and both instruments: the best case is +0.0010 < +0.0012. The block does not ship, with or
    without the hook, as 537 already predicted.
  - Overall vs 446bj2z ≈ 0 at P8.
- **Falsifiers (537's, unchanged):**
  - pds / mse ≠ 446bj2zcp's by more than 1e-6 raw;
  - |Δjac| vs 446bj2zcp > 0.001 (predicted +0.0001);
  - Δfid vs 446bj2zcp < −0.0003 raw.
- **reach:** unpriced. The H1 increment at f\* is +0.0012 … +0.0022 raw, but cp's own reach was mis-signed on H1 at every
  level.
- **The formula after 496a reads** (fixed now, so 10-05's read only sets T):
  - I-662: primary = **+0.000575 + 0.000273 × T** (P8); +0.001173 + 0.000650 × T (P1); +0.000023 − 0.000359 × T (P9).
  - f 1: primary = +0.000349 + 0.000299 × T (P8); +0.000980 + 0.000502 × T (P1); −0.000099 − 0.000419 × T (P9).
  - Bars: if 496a FAILS with Δnmae between −0.0003 and 0 (T 0 … 0.15), the I-662 P8 primary is still +0.0006, because the
    fid term carries it (it stays ≥ +0.0004 down to T −0.6). That is the block-specific repair; 496a cannot test it,
    because the plain hook moves no calls.

## What this read can and cannot grade
- **It cannot separate I-662 from the f 1 reading.** Their CPH predictions differ by ≤ 0.0002 overall, inside the read
  noise (~0.00025).
  - The reason: the block's members keep their signs at every level (fid −, jac + from f 0.35 to f 1), so the level barely
    moves CPH's price.
  - The level test is K16 (10-06, `PRED-662-K16.md`): I-662 says −0.0067, the f 1 reading −0.0030 (≈ 15 read-noise sd
    apart).
- **It does grade the hook's CLEAN4 information on A/B/C under the corrected scaling.** The fid term is the block-specific
  part.

| primary vs 446bj2zcp | reading |
|---|---|
| ≥ +0.0007 | P8 or better: CLEAN4 carries A/B/C-like sign information into the swapped-in calls |
| +0.0001 … +0.0007 | between P9 and P8 |
| ≤ −0.0002 | the P9 floor: no usable sign information on A/B/C, or worse |
| outside [−0.0007, +0.0021] | the corrected instrument's CPH pricing FAILS; record which member |

- **Member checks:**
  - Δfid vs 446bj2zcp: I-662 P8 +0.0010 (P9 0.0000 … P1 +0.0021). Fid RISES vs 446bj2zcp.
  - vs 446bj2z, fid still falls: about −0.0026 = cp's −0.0036 + the repair.
- **For October (652 F2):** a POSITIVE does not make the block ship (the SHIP stat stays < +0.0012). It only says the hook is
  worth keeping IF the block is ever adopted.
