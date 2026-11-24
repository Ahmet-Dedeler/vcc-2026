# PRED-661-dnmag: DNMAG (d4cd5f's DN module, `d_m200`) on 630's entry — 618-registered board prediction

[21dd4c], 2026-10-03 ~19:05Z, committed BEFORE DNMAG's board read. Coordinator's 661 item (2). The arm is DN added on top
of the 630 entry, so its effect is (630 + DN) − 630.

## Source (618; `618/work/stacks618b.json`: "(b) 630 + DN" minus "630 alone (×1)", the same 3 seeds and the same plain base)
| 618 raw Δ (DN increment over 630) | f 1 | f 0.5 | **f 0.35** | cmp |
|---|---|---|---|---|
| nmae | −0.00104 | −0.00319 | **−0.00491** | −0.00461 |
| fid / jac / reach | 0 / 0 / 0 | 0 / 0 / 0 | **0 / 0 / 0** | 0 / 0 / 0 |
| overall (nmae only) | +0.00029 | +0.00087 | **+0.00135** | +0.00127 |

DN's increment is strongly regime-dependent: it acts on DOWN magnitudes of the strongest responders, which f 0.35 / cmp
weight most.

## Board prediction (619: board nmae ≈ 618 f 0.35 × 1.0–1.2)
- **Point:** board Δnmae (DNMAG entry − 630 entry) = **−0.0054 raw** (×1.1), i.e. **+0.0015 overall**.
- **Band:** **[−0.0059, −0.0049] raw**, i.e. **+0.0013 … +0.0016 overall**.
  - Sensitivity: if the board reads like f 0.5, it is +0.0009 … +0.0011. A board value between +0.0009 and +0.0013 is
    read as "the board's regime sits between f 0.5 and f 0.35", not as a transfer failure.
- **Members:**
  - **fixed:** pds and mse (sums pinned);
  - **fixed:** fid, jac and reach (rank-neutral on 618: 0 Δ in every seed);
  - **moving:** nmae only.
- **Reading:**
  - inside the band confirms 618 for DOWN-magnitude levers;
  - below +0.0009 says DN's 618 gain does not reach the board;
  - a pds / mse move or |Δjac| > 0.0005 falsifies the rank-neutral premise.

## I-662 ALTERNATIVE (added 2026-10-04 ~17:55Z, BEFORE DNMAG's read; the f 0.35 prediction above stands as registered)
- **Why:** K20's read put the board's level for call-set changes at f ≈ 0.61 (662). DN is the most level-dependent lever
  priced so far on 618:
  - its increment is −0.00491 at f 0.35, **−0.00204 at f 0.61** and −0.00104 at f 1 (`alt662.json`; known answer vs
    stacks618b exact);
  - so this read tells which level the board's nmae reads at for magnitude levers.
- **I-662 board prediction:**
  - 618 f\* × 1.05–1.30 (mid 1.17): Δnmae **−0.0024 raw** (band −0.0027 … −0.0020), i.e. **+0.00065 overall
    (+0.00056 … +0.00075)**.
  - The f 0.35 line gives +0.0015 (+0.0013 … +0.0016).
- **Reading (replaces "below +0.0009 says DN's 618 gain does not reach the board"):**

| board value | reading |
|---|---|
| +0.0012 … +0.0017 | nmae reads at f 0.35 for magnitude levers |
| +0.0005 … +0.0008 | nmae reads at f ≈ 0.61, one level for all members (I-662) |
| between +0.0008 and +0.0012 | between the two levels |
| < +0.0004 | DN's gain does not reach the board at either level |
