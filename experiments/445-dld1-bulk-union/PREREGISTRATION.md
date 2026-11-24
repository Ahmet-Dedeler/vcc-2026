# 445r — DLD-1 as the bulk's second union residual (g 0.2). Preregistration, 2026-09-24 21:40Z
Written before any 09-25/26/27 entry exists. The arm is one variable on whichever bulk is the incumbent when it is built
(433h, 433k or 446b): that base's exact config with `VCC_BULK_LOAD` swapped to the DLD-1 stack dump. Its config diff vs
the base is exactly `bulk_load` + `bulk_load_sha256`, plus `bulk_ctx_scale` only if the re-level has to follow.

## Construction
- The dump is `work/dump_h1x2_cd4x15_dld1u2g020.npz` (sha256 `575b4fda…`): the stack dump's env (H1 ×2, CD4 ×1.5, key
  wshare_lo0.3_hi3) plus `VCC_UNION2_SOURCES=dld1_lowmoi VCC_UNION2_G=0.2`. That is 248's second-union slot, sized on
  DLD-1's own carrying rows at its own dose.
- **KA:** the same dump without union 2, remade tonight with builder `63dfd7d0` / driver `5463fe28`, equals 417's
  `dump_h1x2_cd4x15.npz` array for array.
- **Re-level:** ntclevel ρ is re-measured on the built arm. If any context moves more than 0.5 % in scale from its base,
  the factors follow at matched ρ (ρ ∝ f^1.86).

## Price (PLOLO on the q 1.5 stack, 445 table; realised KMAX 16)
- Cancer folds: rel cos +0.87 / +0.83 / +0.72 %, β-emulated pds +0.0016 / +0.0003 / +0.0011.
- kolf2: +0.11 % / −0.0014.
- The permuted-row control loses about 1 % cos on every cancer fold.
- Composition steps transferred at k ≈ 1–1.3 (417e).
- **Prediction: overall +0.0005 … +0.0015; mse −0.001 … −0.004 raw; pds −0.0005 … +0.002 raw; DE members pinned.**

## Bins (vs the base's board row)
- Overall ≥ +0.0005 with mse Δ ≤ −0.001 means DLD-1 transfers. It joins the bulk and the October stack: DLD-1 gains
  coverage on a disjoint panel, 67 → 73 on 225's pessimistic draw.
- |Δoverall| < 0.0005 means flat; keep it out (simpler October recipe).
- Δoverall ≤ −0.0005 means close.
- A DE member outside the ballast's realisation noise (nmae ≤ 0.001, fid ≤ 0.001, jac ≤ 0.001, reach ≤ 0.01) means a
  bulk-to-carrier leak: stop.
