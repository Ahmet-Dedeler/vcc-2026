# ALLSW board result (10-14; read by the coordinator [400b6e] at 23:06Z 10-14)

Entry **4WNmsNfdQhYx6bfib0dR**, package c86515cc…, submitted 22:16Z (10-14's slot), published 23:04Z: **overall 0.2760851, rank 59 — NEW BEST**
(+0.0011661 vs 498p2g 0.2749190). Members: pds 0.81527, mse 0.39588, nmae 0.17599, fid 0.00558, reach 0.23536, jac 0.02843.
Reader: `read658.py allsw 4WNmsNfdQhYx6bfib0dR --hook HoQY4xrzYer5y232oh85` → work/read_allsw.json.

## Falsifier F1 fired, diagnosed benign (coordinator override, written before any decision uses this read)
- read658 stops on F1: Δpds raw (ALLSW − 498p2g) = **−3.716e-6** > 1e-6. That value is exactly **1 / (3 × 300 × 299) = 3.7161e-6**:
  ONE nearest-neighbour rank flip among the 300 × 299 cosine comparisons of one context.
- Cause: the 28 paired ±1 count swaps (11 of 900 blocks, net 0 per block) that bulkid658 found before the slot
  (DIAGNOSIS-ALLSW-bulkid.md), rounding ties under emitter 0f56209d vs 498p2g's edb5b8dc. A ±1 count in a near-tied pair can flip one
  rank. mse raw moved +0.000000, the bulk term is −0.0000013 overall: the bulk is 498p2g's to the scorer's resolution.
- So the read stands. F1's 1e-6 rule is too tight for builds whose emitter pin moved: one pds rank flip is 3.7e-6.

## Reads (hand-computed from the reader's printed deltas)
| quantity | value | bin |
|---|---:|---|
| X = 498p2g − 446bj2 (DE leak) | +0.0000573 | — |
| Y = Σ single DE reads (zs + hook) vs 446bj2 | +0.0011045 | — |
| ΔDE(ALLSW − 498p2g) | +0.0011674 | — |
| **R2: I = ΔDE(ALLSW − 498p2g) − Y** | **+0.0000629** | **ADDITIVE** (\|I\| ≤ 0.00025); I_nmae −0.000126 raw, I_fid +0.000004, I_reach +0.000185 |
| **R1 (v2 band): ΔDE(F) − ΔDE(Y) + 0.0002** | **+0.0003202** | **PASS** (≤ 0.0006; window I ∈ [−0.0008573, +0.0003427]) |
| prediction (additive) | 0.2760235 | read 0.2760851 (+0.0000616) |

## What follows (fixed by the prereg §6)
- R1 PASS: a D/E/F entry may pair one bulk switch and one DE switch (637 §1), if the members come back.
- R2 ADDITIVE: Σ-based predictions for the October recipe stay additive (zs + hook on the cap6+P2g bulk).
- 10-15's DN read (on 630arm's base) and the PIE bulk arm 702 (on ALLSW's base) can be composed additively within this noise.
