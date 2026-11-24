# 471 — the denoised H1 bank (h1_2025_full) in the stack's H1 union: 471c is built and gated (05:20Z 09-27)

## What was built
- **b_ctrl for h1_2025_full** (`bctrl471.py`): the recomputed delta equals 458d's bank bit for bit, and the cell counts
  equal the official CSVs. `b_ctrl.npy` sha `2e414a9e…`.
- **Builder `4b8aa34e`** (committed 7cde34f3): `h1_2025_full` registered in `UNION_SOURCE_REGISTRY`.
  - KA-1: unset rebuilds 446b byte for byte.
  - Dump KA: the default union re-dumps 446b's stack array for array.
- **The new dump** (`dump_h1fx2_cd4x15.npz`, `30086b5b…`): 3 panel targets newly covered (HSBP1, SHPRH, ZFP62) and 22
  cleaner rows. The energy-preserving key rescales the other 220 rows by 0.988146.
- **Arms:**
  - 471 (446b's factors): built, not submitted (confounded by the rescale);
  - 471b (factors ÷ 0.988146): ρ A 1.005, B 0.996, C 1.026, past the cliff, **never ship** (`37d1919d`). The key had
    already kept the total energy;
  - **471c** (factors from 471b's ntclevel EV-scale: A 1.0993, B 1.1557, C 1.1587): config diff vs 446b exactly
    `bulk_ctx_scale`, `bulk_load`, `bulk_load_sha256`; CLIFF OK; **ρ A 0.9812, B 0.9804, C 0.9745**, within 0.5 % of
    0.9793 in every context (C at −0.49 %); novelty rc 0; package **`137893bd…`**.

## Status
- 471c is ready for a slot from 10-02. Bins: `PREREGISTRATION.md`, vs 446b.
- It is built on 446b's 410-key stack. The incumbent is now 446bj2 (the LODO key), but the read stays one variable vs
  446b, since bulk members reproduce exactly. October would stack the union swap on the LODO stack.
