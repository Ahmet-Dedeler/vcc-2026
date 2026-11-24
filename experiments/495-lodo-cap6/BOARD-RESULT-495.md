# 495 BOARD RESULT: 0.2739107, a NEW BEST (+0.0016340 vs 446bj2). LODO key at cap 6 TRANSFERS; October lodo_hi_pct DECIDED

Entry `9Hxef2F6Clf4G2hIkNPq`, rank 14. Package
`submissions/495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5.vcc`
(sha256 `137893bd…`). Bins: `PREREGISTRATION.md` (706ad5f0). Reader: `read1003.py`, written 2026-09-27 before any 495 entry.

**Delivery incident.** `495/arm1003.sh` slot 1 started submitting at 10-03 00:01Z. The upload reached 60 %, then failed at
01:37Z ("nodename nor servname provided": DNS lost). Two resume attempts followed, then the Mac slept on battery with the lid
closed from ~01:45Z to 09:59Z (pmset: lid wake at 09:59:27Z). The 536 resume block finished the upload on wake: 799.5 MiB in
5m12s, SUBMIT rc 0 at 09:59:38Z, captured 10:00:03Z. Slot 2 (471c) started at 10:00:45Z, inside its 12:00Z cap. AGENTS §6.2's
rule (AC + lid open through slot windows) was broken by the environment, not the scripts. The retry/resume block worked as
designed.

| member | Δ raw vs 446bj2 | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | −0.001289 | −0.002845 | −0.000474 |
| mse | −0.011991 | +0.012628 | +0.002105 |
| nmae | +0.000615 | −0.001043 | −0.000174 |
| fid | +0.000058 | +0.000193 | +0.000032 |
| reach | +0.000699 | +0.000813 | +0.000135 |
| jac | +0.000021 | +0.000058 | +0.000010 |
| **overall** | | | **+0.0016340** |

Scaled members: pds 0.81239, mse 0.39305, nmae 0.17018, fid 0.00569, reach 0.23369, jac 0.02847.

## Against the bins
- **BIN: TRANSFERS** → October `DEF_DAY_LODO_HI_PCT=96.64697522769474` (the percentile form, 543; == cap 6 on A/B/C).
- Prediction (PLOLO instrument): overall +0.0019 … +0.0029, mse −0.0140 raw. Observed +0.0016 and mse −0.0120 raw: just BELOW the
  predicted band in size, same sign. The 09-29 coordinator note made a FLAT read the instrument's failure; this read is not flat.
  It reads as "the instrument over-predicts by ~15–40 %", and that is recorded as such.
- pds −0.0013 raw (a small cost the prediction did not carry), mse −0.0120 raw (the gain).
- DE leak check (reader): pinned. In 537's preregistered noise check (3σ bands |Δfid| ≤ 0.0001, |Δjac| ≤ 0.00005, |Δreach| ≤
  0.0040 raw) this bulk-only pair sits inside every band.
- Sealed R6 (534 Part B, C10a) had priced the percentile cap at +0.00694 pooled; C10b (10-04) is its next sealed read.
