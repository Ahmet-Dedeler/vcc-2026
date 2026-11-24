# 495 — the LODO key's cap raised 3 → 6 (the third session's 540), built on 446bj2 (plan, 04:42Z 09-27)

**Why.** The third session's 540 (preregistered): the LODO key's cap of 3 binds on 49 of 294 carried targets and costs.
- Raising it is better on every fold: direction +1.7 % at cap 4 and +2.6 % at cap 6 (3/3 cancer folds; kolf2 +3.8 % /
  +8.3 %).
- Calibrated at the key class the board just measured (446b → 446bj2: k_X 1.31, literal pds): cap 6 prices **+0.0025
  overall net vs 446bj2** (+0.0019 at k 1). All four preregistered rules pass.
- A/B/C energy moves ≤ 1.6 %; max row ×2.04. It is orthogonal to 479: capped targets are all in the HIGH tercile.

**Artifact** (`work/pert_scale_415pp_lo0.5_hi6.npz`, sha256 `dc374281…`):
- rule: w = clip(key / median(key[carried]), 0.5, 6) on carried targets, 1 elsewhere. That is `def_lodo_key.py`'s rule
  at `--hi 6`, applied to the pinned `8c82bef3` artifact's own key array;
- KA: the same rule at (0.5, 3) reproduces `8c82bef3`'s w exactly (max |d| 0);
- 294 carried targets; 49 now above 3; 10 at 6 (max key/median 9.99).

**Plan** (after 479's KAs and arms):
1. Register `lodo415_lo0.5_hi6` in the builder's `PERT_SCALE_ARTIFACTS` (a pure addition), with KA-1: unset == 446b.
2. Dump the stack with `VCC_PERT_SCALE=lodo415_lo0.5_hi6` (`dump_h1x2_cd4x15_lodo05hi6.npz`). The dump KA is the same
   env with `lodo415_lo0.5_hi3`, which must equal 446bj2's dump array for array.
3. Arm 446bj2 + the cap-6 dump, then ntclevel, then a rebuild at the EV-scale factors (ρ within 0.5 % of 0.9793). Config
   diff vs 446bj2: `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale`.
4. October: `def_day.sh` stage 2b passes `--hi 6` to `def_lodo_key.py` (a `DEF_DAY_LODO_HI` switch) if the board
   confirms.

## Build record (chain495, 2026-09-27)
- Builder d09a986c → **b4975024** (the registry entry; tests 56/0 and 479/468/472/451/417 pass).
- **KA-1 PASS:** 446bj2's env rebuilds 446bj2's output `711ee43b` byte for byte, driver rc 0.
- **Dump KA PASS:** 446bj2's dump env with `lodo415_lo0.5_hi3` re-dumps `dump_h1x2_cd4x15_lodo05.npz` array for array
  (bulk, carriers, perts, genes).
- **Cap-6 dump** `dump_h1x2_cd4x15_lodo05hi6.npz`, sha `d58750b8…`: the carriers are identical (the DE channel is untouched),
  and all 300 bulk rows change.
  - The energy-preserving constant falls to 0.7189. The capped targets' weight is up to 2× what it was, so their rows grow
    by up to ×1.4377, and every uncapped row shrinks to 0.7189 of its 446bj2 norm (the median).
  - Total non-panel energy ratio 0.9978. The redistribution toward the strong-key targets is the arm; the overall level
    goes back to ρ 0.9793 through the two-round re-level.
- **Round 1** (446bj2's factors A 1.0962 B 1.1461 C 1.1681): config diff exactly `bulk_load`, `bulk_load_sha256` (PASS).
  ntclevel ρ **A 1.0316, B 1.0229, C 1.0741**, all past the cliff. The cap-6 key moves energy onto the strong-key targets,
  and ntclevel's between-construct statistic rises (B 0.241 vs claim 0.234 on A), although the dump's total non-panel
  energy is 0.998 of 446bj2's. Cliff check vs 410b A 1.00 / B 1.07 / C 1.20 ("CLIFF FAIL" at the 1.15 bound, expected
  for an over-levelled round 1). EV-scale 0.9743 / 0.9785 / 0.9549 gives round-2 factors **A 1.0680, B 1.1215, C 1.1154**.
- **Round 2** (A 1.0680, B 1.1215, C 1.1154): config diff exactly `bulk_ctx_scale`, `bulk_load`, `bulk_load_sha256`;
  CLIFF OK (0.95 / 1.03 / 1.07); novelty rc 0; package `2d931b0c…`. ρ **A 0.9805, B 0.9812, C 0.9627**. A and B are within
  0.5 % of 0.9793; **C is −1.7 %**, outside the band chain1003 requires. relevel_iter's secant (C's local exponent 2.373)
  gives **C 1.1235**, with A and B unchanged.
- **Round 3** (`relevel3/chain495.sh`, after chain479b): the same tag, name and stem at A 1.0680, B 1.1215, C 1.1235.
  The round-2 build dir and package are kept as `build_495.r2` and `*.r2.vcc`. Its gate regex name makes chain476 wait for it.
  - **READY (2026-09-27 11:26Z).** Config diff exactly `bulk_ctx_scale`, `bulk_load`, `bulk_load_sha256` (PASS); CLIFF OK
    (0.95 / 1.03 / 1.09); novelty rc 0; package **`ce0a41af…`**. ρ **A 0.9805, B 0.9812, C 0.9788**: all within 0.5 % of
    0.9793 (C's secant step landed at −0.05 %). EV-scale 0.9994 / 0.9991 / 1.0003, so no fourth round.
  - chain1003's own readiness test (logs, package, novelty re-check, ρ from `build_495/*.ntclevel.json`) was dry-run by
    hand at 11:28Z: all pass, so 10-03's s1 is 495 unless something changes the package before 10-02 20:00Z.
- Builder b4975024 committed by its committer (643d11e8) with KA-1 and the dump KA.
