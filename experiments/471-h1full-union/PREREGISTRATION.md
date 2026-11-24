# 471 — the denoised H1 bank in the stack's H1 union (preregistration, 2026-09-26 17:22Z, before 471 is built)

## The arm
471 = 446b with its stack dump rebuilt with `VCC_UNION_SOURCES=gse295214_h1,h1_2025_full` in place of the default
`gse295214_h1,h1_2025`, loaded through `VCC_BULK_LOAD`. Everything else is 446b's, including its re-level factors
A=1.1004, B=1.1517, C=1.1720.
- The DE channel is untouched: LOAD swaps only the bulk, and the carriers keep the build's own design.
- Expected config diff vs 446b: exactly `bulk_load`, `bulk_load_sha256`.
- `chain471.sh` runs the known answers before it builds:
  - KA-1: the registry addition leaves the unset build = 446b byte for byte;
  - the dump KA: the default union re-dumps 446b's stack array for array.
- Name: "the hemocytometer lost count". Description: empty.

## Why
- `h1_2025` is the 400-cell harness bank: 255 perturbations and 18,399 control rows, about a third of which is one
  control sample's noise (462c).
- `h1_2025_full` sums all public VCC 2025 H1 cells: ~1,071 per perturbation, 38,176 controls. Disjoint-control
  reliability is 0.69 vs 0.39.
- The union gains 3 panel targets (77 → 80: HSBP1, SHPRH, ZFP62). 22 rows (20 from h1_2025 plus ADNP and BRPF1) get
  cleaner.
- The H1 ×2 union is part of 417e's stack (mse −0.0262 raw for the whole stack). This is a small refinement of one of
  its two H1 sources.
- It carries to October to the extent D/E/F's panel overlaps the H1 screen.

## Prediction (471 − 446b, raw)
- mse −0.0005 … −0.0030; pds −0.0005 … +0.0015. Overall +0.0001 … +0.0006.
- The DE members stay pinned up to the ballast's drift: nmae ≤ +0.0015, |fid| ≤ 0.001, |jac| ≤ 0.001, |reach| ≤ 0.01.
- The re-level: with 446b's factors kept, the cliff check (A_vs_410b within [0.85, 1.15]) must pass. A small bulk change
  moves ρ by < 1 %.

## Bins
- **Adopt** (use `h1_2025_full` in the October stack's union): overall ≥ +0.0002 with mse ≤ 0 raw.
- **Flat:** |overall| < 0.0002. Keep the default union.
- **Harm:** overall ≤ −0.0002. Keep the default, and note that noise reduction in a union member does not pay.
- A DE member outside the pinned bounds means a leak. Diagnose it.

## Slot
Not before 10-02. The rule for that day is fixed before the previous day's entries exist. It needs `build_471`'s gate
lines (PASS, CLIFF OK, novelty rc=0).

## Addendum 21:36Z 09-26 (before any 471 entry exists): the arm to submit is 471b, re-levelled
The new dump differs from 446b's in more than the union rows (`chain471.log`):
- On the 220 rows outside every H1 union, the bulk is exactly a global scale of 446b's: × 0.988146 (residual 0). The
  carriers are × 0.99984.
- The per-target key is energy-preserving, so the three newly covered targets (HSBP1, SHPRH, ZFP62) take energy
  (row change ~100 %) and every other row shrinks 1.2 %.
- With 446b's re-level factors, the arm therefore mixes the bank swap with a −1.2 % bulk rescale below the cliff. 414's
  slope puts that at roughly +0.004 raw mse, against the swap's predicted −0.0005 … −0.003.

**471b** is the same dump with 446b's factors divided by the scale: A 1.1136, B 1.165516, C 1.186059.
- The 220 non-union rows render at 446b's amplitude (to ~1e-6), so the only variable is the 80 union rows' content.
- Expected config diff vs 446b: exactly `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale`.
- `ntclevel417` on the built artifact must give ρ within 0.5 % of 446b's matched level (0.9793) in every context.
- Name: "the pH meter joined a union". Description: empty.

**Prediction and bins: unchanged** (vs 446b: mse −0.0005 … −0.0030 raw; overall +0.0001 … +0.0006; adopt at
≥ +0.0002 with mse ≤ 0). The unlevelled 471 package is built but will not be submitted: it is the confounded arm.

## Addendum 23:54Z 09-26 (before any 471 entry exists): 471b fails its ρ check; the arm to submit is 471c
471b (factors ÷ 0.988146) built cleanly, with config diff exactly `bulk_ctx_scale`, `bulk_load`, `bulk_load_sha256`
and CLIFF OK. But ntclevel reads **ρ A 1.0050, B 0.9960, C 1.0261**, against 446b's matched 0.9793. That fails the
preregistered ±0.5 %, and A and C are past the cliff, so **471b is never shipped** (package `37d1919d`).
- Why: the energy-preserving key had already kept the total energy. The union rows gained what the 220 other rows lost
  (1.2 %), so restoring those rows' amplitude over-levelled every context. The correction in the previous addendum was
  wrong.
- **471c** takes 471b's factors × ntclevel's EV-scale (0.9872 / 0.9916 / 0.9769): A 1.0993, B 1.1557, C 1.1587.
  - ρ then reads 1.0050 × 0.9872² = 0.9795, 0.9960 × 0.9916² = 0.9793, 1.0261 × 0.9769² = 0.9792, matched to 446b's level.
  - This is how 446b itself was levelled.
- Expected config diff vs 446b: exactly `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale`. ntclevel on 471c must give ρ
  within 0.5 % of 0.9793 in every context.
- Name: "the pH meter joined a union" moves to 471c. Bins unchanged.
