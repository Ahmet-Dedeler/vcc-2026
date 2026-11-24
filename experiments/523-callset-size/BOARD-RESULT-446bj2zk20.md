# 446bj2zk20 (K 2,000) BOARD RESULT: 0.2721097 (−0.0007163 vs 446bj2z). The call-set enlargement HURTS fid; f 0.35 call-set pricing is wrong-signed

Entry `kJfe2dP9LEJawDEf0pMm`, 10-04 slot 2. It was the coordinator's fill, since rule1004's fill list was exhausted.
- **Upload:** 05:39Z start over an iPhone hotspot. The session's offset query kept failing ([Errno 8]); arm1004 made 120 tries, then
  resume1004s2 took 23 more.
- **Finish:** SUBMIT rc 0 at 16:07:21Z, captured 16:07:42Z.
- **Package:** `submissions/446bj2zk20-…-cpm1600.vcc` (a8a593dd…).
- **Bins:** `523-callset-size/PREREGISTRATION.md` (+ the 09-30 dose amendment); instrument card
  `610-board-forensics/PREDICTIONS-611-K16-K20.md` (a995f0e2).

| member | Δ raw vs 446bj2z | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | 0 | 0 | 0 |
| mse | 0 | 0 | 0 |
| nmae | −0.000059 | +0.000150 | +0.000025 |
| fid | −0.001689 | −0.005578 | −0.000930 |
| reach | +0.000600 | +0.000712 | +0.000119 |
| jac | +0.000185 | +0.000419 | +0.000070 |
| **overall** | | | **−0.0007163** |

- **Falsifiers clean:** pds / mse exactly 0, so there is no bulk leak; nmae moves −0.00006 raw, well within ±0.002.
- **Primary Δ(fid + jac) = −0.00086 overall.** K 2,000 cannot be adopted: the dose rule needs ≥ +0.0012. The decision waits for
  K16 on 10-06, but with H1's shape (K 2,000 ≥ K 1,600 everywhere) K16 is expected no better, so October keeps 2,200.
- **Instrument card graded: FAIL.** It predicted Δfid +0.0070 raw (+0.0034 … +0.0097), Δjac +0.0012, Δ(fid + jac) +0.0044 overall.
  The board gives Δfid −0.0017 raw, the WRONG SIGN. By the card's own rule, every f 0.35 call-set price is wrong-signed: for
  call-set size changes the board behaves like H1's f 1 regime (fid −0.010 … −0.019 there). This agrees with 446bj2zcp's fid
  −0.0036 (10-01).
- **Consequence:** CPH (10-07 s1) was priced only in the f 0.35 call-set mode; expect its fid to fall as well. Its October case
  already needed a written reason (652 F2).
