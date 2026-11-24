# 710: a board point for H1 above ×2 (H1 ×2.5 on ALLSW), before 708 could raise it in the dark (preregistration)

Written 2026-10-15 ~22:05Z by the science session [d80a10], on the coordinator [2d6f46]'s YES (10-15 ~22:00Z). No dump,
build, package or entry exists at this time.

## Why
708 Part A (day-of re-pricing of the stack weights on D/E/F's panel) will most plausibly propose H1 above ×2. On A/B/C's
panel, 438's instrument picks H1 ×2.5 / CD4 ×1.5 under the D/E/F rules (`abcdef`, +0.00044 over H1 ×2 at CD4 ×1.5) and
stays at ×2 under the A/B/C rules (`abc`). So the choice sits on a +0.0004 margin. H1 has board points only at ×1 (h1union,
09-04) and ×2 (417e, 09-24). This entry gives the instrument's H1 step one A/B/C board point, on Arc's own lab, panel and
truth, before 708 could act on it.

## The arm (one variable vs ALLSW, entry `4WNmsNfdQhYx6bfib0dR`, 0.2760851)
- **710h125** = ALLSW's exact env, with `VCC_BULK_LOAD` = a re-dump of ALLSW's stack at **`VCC_H1_UNION_G=0.70`**
  (H1 ×2.5). ALLSW's stack is 495's cap-6 dump, d58750b8, at 0.56.
  - Every other dump knob is 495's: NORMREL_EXP −1.5, CD4_PARTIAL_G 0.5, NORMREL_MODE cap, DE_SIGN_PRIOR clean4q80,
    PERT_SCALE lodo415_lo0.5_hi6 (dc374281), PERT_SCALE_BALLAST 1, CD4_TARGET_F 0.51 (CD4 ×1.5).
  - Then the usual re-level: r1 at ALLSW's factors → ntclevel417 → r2 (rho within 0.5 % of 0.9793; at most one r3).
- **Not bulk-only.** The dump's `carriers` are the stack composition itself (705's KA-cap6 showed bulk = carriers × w·c
  exactly). So the H1 dose moves the DE design too. That is exactly what October's `DEF_DAY_H1_STACK=x2.5` would do (690i
  acts in the dump subshell), so this entry tests the real October implementation, both channels.
- **Known answer (before the arm's dump is used):** a re-dump at `VCC_H1_UNION_G=0.56` on today's builder must reproduce
  d58750b8's bulk, carriers, perts and genes arrays (np.array_equal). Otherwise STOP, because the dump path has drifted.
- **Expected config diff vs ALLSW:**
  - r1: `bulk_load`, `bulk_load_sha256`, `bulk_sharpen_pt`. The tercile edges scale with the key's energy constant when
    the composition changes; the tiers do not.
  - r2 adds `bulk_ctx_scale`, plus `de_sign_emitter_sha256` if the pin moved.

## Class
**Transfer**: the composition weights, a choice made on held-out folds. Its October use goes through 708's rule on D/E/F's
own panel. This A/B/C read can only block or calibrate.

## Primary statistic, prediction, bins (board, vs ALLSW's own entry)
- **Primary: Δoverall (all six members).** Both channels move. Reported beside it: Δ(pds + mse) (exact) and the DE members
  (true effect plus re-roll, σ ≈ 0.00025).
- **Prediction:**
  - Bulk part ≈ +0.0001 … +0.0004. 545, under the production cap-6 key: H1 0.70 at CD4 ×1.5 gives +0.34 % direction
    (2/3 folds) and Δpds_lit −0.0010, net ≈ +0.0002. 438's K3 table gives +0.00046.
  - DE part: no instrument prices it. Prior 0 ± 0.0005 (more H1-union rows in the carriers).
  - Total ≈ +0.0002, in −0.0005 … +0.0008.
- **Bins (the coordinator's):**
  - **HARM, Δoverall ≤ −0.0003:** BLOCKS any 708 move that raises H1 above ×2 (708 AMENDMENT 2; decide708's R710=HARM).
  - **NULL, −0.0003 … +0.0003:** no block. 708's rule decides on D/E/F's panel alone.
  - **POSITIVE, ≥ +0.0003:** the instrument's H1 step is calibrated on Arc's panel. No block.
- **Falsifiers (STOP and diagnose before reading):** the dump KA fails; a config diff beyond the expected keys; a final rho
  outside the band; a zero-gene audit > 0.

## Tries
- Board, on the H1 dose: ×1 (09-04) and ×2 (09-24) exist; this is the first point above ×2 (one try).
- Local: 545's grid (4 H1 × 6 CD4, cap-6 key) and 438's grid (K3), both chose from many points. The board read is one
  fresh draw, so its bins stand uncorrected.
