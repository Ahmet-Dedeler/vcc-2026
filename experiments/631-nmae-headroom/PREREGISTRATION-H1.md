# 631 — registered H1 predictions for the DOWN pass, written BEFORE any delivered-level number exists (2026-10-02T16:18:30Z)

Runs in flight (`work/chain_dn.sh`, started ~16:15Z; no output read at the time of writing; chain_dn.log is empty):
- **DN** = `b_k000_q20r3p10_D200`: the DOWN pass alone (dn631.dn_pass, m = 2, net |ΔAUC| ≤ 1e-4, ballast sub-gate
  payback) on 446bj2z's zs base;
- **UPDN** = `b_k100_q20r3p10_D200`: 630 v3, then the DOWN pass.
H1 real emitter, seeds 0–2, priced with price630.py at 521's regimes (raw Δnmae and overall-equivalent).

## Basis
- proxy631: ×2 on all called DOWN genes at FULL delivery is worth −0.0227 (f 1) / −0.0251 (f 0.35) raw.
- The timing run on 2 targets delivered 7 % and 11 % of the need (`work/time_dn.log`), bound by the ballast's
  sub-gate pool.
- The pass serves the most negative genes first (the head-like, highest-value ones), so the gain per delivered
  fraction should run at or above linear.

## Predictions (fixed now)
| quantity | DN | UPDN |
|---|---|---|
| mean delivery fraction (dn_frac) | 0.05 … 0.15 | same as DN ± 0.02 |
| Δnmae raw, f 1 | −0.002 … −0.005 | v3 (−0.0084) + DN, ± 0.0005 |
| Δnmae raw, f 0.35 | −0.002 … −0.006 | v3 (−0.0106) + DN, ± 0.0005 |
| overall-eq, f 0.35 | +0.0006 … +0.0017 | v3 (+0.0029) + DN, ± 0.0002 |
| fid / jac / reach | exactly 0 (0 call flips) | exactly 0 |
| pds | exact; mse within 1e-4 raw | same |

**Reading rule:**
- If DN's delivery stays ≤ 0.25, the coordinator and I decide together whether an emitter redesign is worth it (their
  condition).
- If UPDN departs from additivity by more than 0.0005 raw, the two passes interact and need separate study.
- v7 (`up7_630.py`: DOWN payments for UP receipts under the payer gene's net-U budget) is priced after this, with its own
  prediction written first.

## v7 predictions (written 2026-10-02T16:19:18Z, before any v7 or DN result is read)
Arms (H1, seeds 0–2, price630):
- **C5** = `b_k100_q20r5p10`: v3 with 5 rounds (the control for the extra rounds);
- **V7r3** = `b_k100_q20r3p10_P200`;
- **V7r5** = `b_k100_q20r5p10_P200`.
v7 pays 630's UP receipts with counts of called DOWN genes (×2 target cap) under each payer gene's net-U budget; any
remainder is paid in sub-gate counts. Move size is bounded by sub-gate as in v3, so it is never short.

| quantity | prediction |
|---|---|
| DOWN-paid share of payments (dn_paid / (dn_paid + sub_paid)) | 0.05 … 0.40 |
| V7r3 vs v3, overall-eq f 0.35 | +0.0000 … +0.0008 (DOWN Δ on payer genes; UP unchanged: same rounds) |
| V7r5 vs C5, overall-eq f 0.35 | +0.0002 … +0.0015 (saved sub-gate buys extra UP moves, plus the DOWN Δ) |
| C5 vs v3 | within ±0.0003 (more rounds without new currency do little) |
| fid / jac / reach | exactly 0 for every arm |
