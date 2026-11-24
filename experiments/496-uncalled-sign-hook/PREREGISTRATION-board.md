# 496a — the uncalled-gene sign hook on the board: preregistration (2026-09-27T06:57Z, before 496a is built and before any 09-28 entry exists)

**Arm.** 496a = 446bj2 + `VCC_DE_MDE_SIGN=clean4:0.05` (chain496p: emitter 1658ae2c, builder 24b95dec; KA-H1 and KA-1 must pass
first). A DE-only change: the carriers' class split moves; per-gene totals, the ballast and the bulk do not. Config diff vs
446bj2: exactly `de_mde_sign` + `de_sign_emitter_sha256`. Slot: 10-05 s1 (`rule1005.py`). -m: queue #8.

**Prediction (vs 446bj2, raw).** H1 (496, 3 seeds): nmae −0.00146 (seed sd 0.00002); the MDE family transferred at 1.76×
(433a), so the board range is **nmae −0.0015 … −0.0035**. fid and jac pinned. reach is not priced: on H1 it is a lottery on
two small-gate targets (496 addendum). Overall about +0.0004 … +0.0010 before reach.

**Falsifier (read first).** pds and mse raw must equal 446bj2's exactly (they reproduce to 16 significant figures on the
board). Any difference is a bug: diagnose before reading anything else.

**Bins (Δnmae raw vs 446bj2):**
- **TRANSFERS:** ≤ −0.0010. October adds the hook (`VCC_DE_MDE_SIGN=clean4:0.05` through a `def_day.sh` pass-through, with
  its own rehearsal), stacked with zs if zs transfers (496z: additive on H1).
- **PARTIAL:** −0.0010 … −0.0003. Recorded; not shipped without a second read.
- **FAILS:** above −0.0003. Closed on the board.
- **Leak check:** |Δfid| or |Δjac| > 0.001, or |Δreach| > 0.01, means diagnose (report reach in every case).
