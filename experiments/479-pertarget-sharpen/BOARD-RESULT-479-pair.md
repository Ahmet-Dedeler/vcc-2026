# 479 pair BOARD RESULT (10-02): P2g 0.2733131 TRANSFERS (new best); S_P2 0.2715625 FAILS → October sharpen_pt = P2g

Entries: **479p2g** `ioNMsGwrG2EvYuoB5Zu7` (slot 1, published ~00:27Z, uploaded 2.6 GiB in 11 m 43 s) and **479sp2**
`NPiJoPH47VYhiWQHhoOG` (slot 2, published ~00:55Z). Both vs 446bj2 `YvVru3DAF89Jiyr34olR` (0.2722766). Bins:
`PREREGISTRATION.md` (b692c2c5). Reader: `read1002.py` (written before any 10-02 entry). Read 2026-10-02T00:58Z.

| arm | overall | Δ vs 446bj2 | pds raw | mse raw | nmae raw | fid raw | reach raw | jac raw | bin |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| P2g = tercile:2.5,2.0,1.75 (key-ratio terciles) | 0.2733131 | **+0.0010364** | +0.001587 | −0.003739 | +0.000120 | +0.000004 | −0.000836 | −0.000011 | **TRANSFERS** |
| S_P2 = stercile:2.5,2.0,1.75 (strength terciles) | 0.2715625 | −0.0007141 | −0.001033 | −0.000287 | +0.001798 | −0.000000 | +0.000777 | −0.000017 | **FAILS** |

- P2g: bulk-only gain (pds +0.000578, mse +0.000657 overall), DE in noise, leak check clean; ballast nmae +0.00012 (below
  the expected +0.0006 … +0.0010). **New best entry: 0.2733131** (446bj2z was 0.2728260; P2g carries no zs rider).
- S_P2: loses pds, and its low tercile's ballast is under-priced (nmae +0.0018 raw, above the +0.0015 check).
- **October (pair rule): `sharpen_pt` DECIDED `tercile:2.5,2.0,1.75`** (october_recipe decide, 00:58:13Z). Held-out record:
  on 438's October-like panels S_P2 was best (+1.32 %), plain P2g +0.48 % (2/3; MID empty on 2 efp folds), carried-only P2g
  +0.84 % (3/3): the board and the off-board instruments disagree on S_P2; the board's FAIL is a ballast/nmae and pds cost
  on A/B/C that the bulk-only instruments do not see.
- **Upgrade queued:** 493h/492's carried-only cut points (non-carried targets → MID; read `carried` from the key artifact)
  as a builder knob with known answers before Oct 22 (on A/B/C it moves 4 targets, 98/104/98).
- Both are bulk-only pairs inside 537's preregistered 3σ DE noise bands (reach −0.0008 / +0.0008, fid ~0, jac ≤ 2e-5 raw).
- Class (AGENTS §11.1c): transfer (shape). Tries on this switch: 2 board arms; PLOLO + 9 October-like folds (492/493).
