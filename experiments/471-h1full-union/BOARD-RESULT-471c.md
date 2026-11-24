# 471c BOARD RESULT: 0.2696638 (−0.0003790 vs 446b). Bin HARM (on DE noise; the bulk is flat); keep the default union

Entry `yzNoMVih0pKqoDjdqrXm`, rank 15, `495/arm1003.sh` slot 2 on 10-03.
- It started at 10:00:45Z (after slot 1's late capture). The coordinator paused it at 10:53Z over an iPhone hotspot.
- It resumed at 12:08Z after Ahmet OK'd hotspot data: 1.3 GiB in 57m51s at 403 KiB/s. SUBMIT rc 0 at 13:19:48Z, captured 13:20:07Z.
Bins: `PREREGISTRATION.md` §Bins. Reader: `495/read1003.py` (written 09-27).

| member | Δ raw vs 446b | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | −0.000045 | −0.000054 | −0.000009 |
| mse | −0.000236 | +0.000252 | +0.000042 |
| nmae | −0.000142 | +0.000237 | +0.000040 |
| fid | −0.000015 | −0.000053 | −0.000009 |
| reach | −0.002407 | −0.002673 | −0.000445 |
| jac | +0.000006 | +0.000017 | +0.000003 |
| **overall** | | | **−0.0003790** |

Scaled: pds 0.81467, mse 0.36545, nmae 0.17155, fid 0.00552, reach 0.23240, jac 0.02839.

## Against the bins
- **BIN: HARM** (overall ≤ −0.0002) → keep the default union. There is no October switch to change.
- **Reading:** the BULK effect was ~0. mse −0.00024 raw against a predicted −0.0005 … −0.0030; pds −0.00005 against −0.0005 …
  +0.0015. The overall loss is reach −0.0024 raw, which is inside 537's 3σ re-roll band (|Δreach| ≤ 0.0040). A bulk change that
  moves ballast K re-draws pass 1 (537 §2.6).
- So the honest summary is "the h1_2025_full swap does nothing measurable for the bulk": swapping a noisier union member for its
  full-data version does not pay. The prediction (+0.0001 … +0.0006 overall) missed in size, not in a way that changes any
  October decision.
