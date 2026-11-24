# 433z / 433j — the 09-26 candidates. Preregistration, 2026-09-23 16:55Z (no 09-24, 09-25 or 09-26 entry exists)

Both are one variable against 433h (`e7495644…`, a 09-25 candidate). Which of them goes, and when, is decided after
the 09-24/25 reads by `docs/NEXT-SLOT-DECISION.md`. Their readings and bins are fixed now.

## 433z (`b4598c84…`) = 433h + `VCC_DE_MAG=dn_k250_zs` (config diff exactly `de_mag`)
429's exact-U down-gene pass on the carriers. Every per-gene total, library, ballast row and rank sum U is unchanged;
only the arithmetic means of down-regulated "mass-holder" genes move.
- H1 prices it on top of the coarse MDE at nmae **−0.00183 ± 0.00007** (435, 3 seeds), with the other five members
  bit-identical.
- A/B/C capacity is 0.70× H1's (429: A 0.91×, B 0.14×, C 1.04×). The pure-magnitude transfer window is ×1 … ×2.3.
- **Predictions (433z − 433h, raw):**
  - pds and mse identical (the bulk is untouched; CLIFF is identical to 4 decimals);
  - fid, reach and jac identical to within float noise (U exact);
  - **nmae −0.0013 … −0.0029**;
  - overall **+0.0003 … +0.0008**.
- **Bins:**
  - Δnmae ≤ −0.001 means zs transfers, and it rides the October stack.
  - −0.001 < Δnmae ≤ 0 means below readability (the ordering channel's sd does not apply to a pure-nmae arm, but
    0.001 is the smallest H1-grounded effect). Keep it only if the board reads ≤ 0.
  - Δnmae > 0 means drop zs.
  - Any other member moving by more than 1e-6 means a bug: stop.

## 433j (`e1a34076…`) = 433h with the stack dumped under 420's LODO key at clip(0.5, 3) (config diff exactly `bulk_load`, `bulk_load_sha256`)
The per-target bulk scale uses 420's leave-one-donor-out key instead of 410's, clipped at (0.5, 3); the stack
composition is the same. The DE channel is pinned (same carriers; ballast).
- 420's PLOLO price for the key step is about **+0.003 overall** at the cancer folds' rate. The measured transfer for
  bulk-shape steps is k_sh ≈ 0.42–0.65 of those folds (417a), so +0.0013 … +0.002.
- CLIFF: A 0.9637 / B 1.0346 / C 1.0846, against 433h's 0.9638 / 1.0337 / 1.0837. The energies are unchanged to 0.1 %,
  so the re-level still holds.
- **Predictions (433j − 433h, raw):**
  - **pds +0.001 … +0.004**;
  - **mse −0.004 … −0.012**;
  - DE members within the ballast's realisation noise (nmae ≤ 0.001, reach ≤ 0.01);
  - overall **+0.001 … +0.003**.
- **Bins:**
  - overall ≥ +0.001 with pds ≥ 0 means the LODO key transfers. It replaces 410's key in the October stack; the key
    script is `scripts/def_pert_scale_key.py`, and the LODO variant needs its own day-of fit (write it then).
  - |Δoverall| < 0.001 means flat; keep 410's key.
  - overall ≤ −0.001 means close the LODO key.
