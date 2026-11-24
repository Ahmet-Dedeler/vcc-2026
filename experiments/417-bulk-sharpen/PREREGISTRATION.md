# 417 preregistration: bulk-only magnitude sharpening through ballast

Written 2026-09-22 20:35Z, before any 417 artifact existed. Author: PI session.

## What changed in the plan, and why

The 09-23 first slot was going to 410c (the uncapped key through ballast; waiter `arm410c.sh` was armed). Worker 414
(`experiments/414-bulk-reprice/RESULTS.md`) then priced 410c at about flat against 410b (−0.001 to +0.001 overall: a
direction gain of −0.009 raw mse, cancelled by +0.0033 raw of amplitude past context A's cliff and −0.0026 board-emulated
pds). It priced 169's bulk sharpening at +0.006 to +0.016 overall. The sharpening slot is worth more and says more, so it
takes the first slot. 410c is not dropped: `arm417a.sh` falls back to the untouched `arm410c.sh` if 417a fails a gate.

## The arm (417a)

- **One variable vs 410b:** `VCC_BULK_SHARPEN_Q=1.25` (builder 400ef6e3, commit 907d6eb0). Everything else is 410b's
  environment, including the capped key `wshare_lo0.3_hi3`.
- **Transform:** `sign(h)|h|^1.25` on each context's 10-CPM render gate, over non-panel coordinates. Each row's gated
  non-panel L2 norm is restored. Panel columns, including the row's knockdown, are untouched. It is applied to the
  ballast bulk render only. The carriers keep 395's DE design, so fid, reach, jac and nmae stay pinned by construction.
- **Why protect the panel columns:** 414 priced the arm with every panel coordinate zeroed. The knockdown is the largest
  |h| in most rows. Left inside the norm, it would absorb energy under the power and cut the scored coordinates: a
  second, unpriced variable.

## Build gates (in `run417a_build.sh`; the waiter falls back to 410c if any fails)

1. The config diff vs 410b's meta is exactly `bulk_sharpen_q: <absent> -> 1.25`.
2. `cliffcheck414.py`: render-space energy A_c/A_c(410b) in [0.93, 0.995] in every context. 414 predicts about 0.96.
   A value above 1 means the implementation is not what was priced.
3. novelty rc=0 and the package sha matches.

## Predictions vs 410b (414 §3/§7, cancer folds to the kolf2-like transfer)

| member | predicted |
|---|---|
| pds raw | +0.006 … +0.019 |
| mse raw | −0.020 … −0.050 |
| nmae / reach / fid / jac raw | within the GO-410c noise band: \|Δnmae\| ≤ 0.010, \|Δreach\| ≤ 0.012, \|Δfid\| ≤ 0.003, \|Δjac\| ≤ 0.001 |
| **overall** | **+0.006 … +0.016** |

K3+sh@1.25 was not in 414's grid. K3 and Kinf take the same sharpening increment on every fold (K3+sh@1.15 − K3 =
Kinf+sh@1.15 − Kinf to 1e-4 on cancer, and to 0.001 rel on kolf2). So the Kinf+sh@1.25 prediction carries over.

## Bins (Δoverall vs 410b's 0.2365722)

- **≥ +0.006:** transfers at least at kolf2's rate. The next arm is the dose, q 1.5 (417b, built tonight).
- **0 … +0.006:** partial transfer. Hold the dose and read which member moved. The second slot goes to 410c.
- **< 0:** PLOLO's magnitude structure does not reach the board. Close the class and send 410c in the second slot. A pds
  loss with an mse gain would repeat a13.
- **Any DE member outside the noise band:** ballast did not hold the DE channel under a sharper bulk. Read the ballast
  diagnostics (K = 16 share, clip mass) before any further bulk arm.

## Addendum, 2026-09-22 21:00Z (before any 417a read): slot 2 gets a DE arm instead of 410c

416 (worker, H1 six-member, seed 0) shows a clean dose response in the carriers' whitening gamma:

| gamma | Δnmae | Δfid | Δreach |
|---|---:|---:|---:|
| 0 | −0.0069 | +0.0025 | +0.0024 |
| 0.3 | −0.0052 | +0.0015 | −0.0044 |
| 0.65 (base) | 0 | 0 | 0 |
| 0.9 | +0.0077 | −0.0009 | −0.0080 |

Realisation noise on nmae/fid is ≤ 3e-4. Reach noise is about ±0.013.

This fits the old note that whitened effects serve pds/mse while the DE ordering wants the unwhitened consensus. It
is only buildable now that ballast separates the carriers from the bulk.

**417c = 417a + `VCC_CARRIER_GAMMA=0`.** The carriers come from an unwhitened second pass through the same frozen
residuals, with rows norm-matched to the shipped rows. The bulk is 417a's. One variable vs 417a (config diff = exactly
`carrier_gamma`). Predicted DE-only change at board slopes: +0.002…+0.004 overall. pds/mse move only by ratio-mode
realisation.

Known answer, run before the arm: a dump-exit run with `VCC_CARRIER_GAMMA=0.65` must print `max|carrier-bulk| = 0`.

**Slot-2 rule, replacing the 20:35Z rule:**
- 417a DE member outside the band → 410c.
- 417a − 410b ≥ +0.006 → 417b (q 1.5).
- otherwise → 417c if its build passed its gates, else 410c.

Rationale: 414 prices 410c at about flat, so it buys less than the DE arm in expected value and in information.

## Addendum 2, 2026-09-22 22:25Z (before any 417a read): the 417b threshold rises to +0.009

420 (worker) mapped 417a's board row to an implied sharpening transfer k_sh, in its RESULTS §6.3 T-read table. At
k_sh 0.55 (417a − 410b about +0.0069), the unlevelled q 1.5 step is worth only +0.0021. That is no better than 417c
(+0.002…+0.004, an orthogonal DE axis) or the stack (+0.0037). q 1.5 pulls clearly ahead only at k_sh ≳ 0.75, which is
417a − 410b ≳ +0.009 (raw mse ≤ 0.703 or raw pds ≥ 0.8627).

**Slot-2 rule, replacing the 20:57Z rule:**
- DE leak → 410c.
- 417a − 410b ≥ **+0.009** → 417b.
- otherwise → 417c.
- 410c remains the fallback for any failed gate.

**The next-day ladder (420):** stack (H1 ×2 + CD4 ×1.5 via LOAD) → + q 1.5 → + per-context re-level. Cumulative
+0.0148 / +0.0182 (k 1 / 1.407), +0.0091 at kolf2-like transfer.

Reading check: if k_sh from pds, (raw pds − 0.851156)/0.01713/0.895, and k_sh from mse (via 187's claims) disagree by
more than 0.3, the cliff model failed. Re-price the re-level rows before any CTX arm.
