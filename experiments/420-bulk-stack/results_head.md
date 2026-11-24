# 420: the bulk STACK, priced as a stack. Composition, key and dose add up; re-levelling refunds sharpening's amplitude cost; the next arm is the stack via VCC_BULK_LOAD unless 417a reads a strong sharpening transfer

Status: **COMPLETE** (2026-09-22 22:15Z). PLOLO on four folds (cancer k562_gwps / orion_hct116 / orion_hek293t decide, kolf2
reported without a vote), the real A/B/C panel, 418's machinery. No vcc, no cell-eval2, no builder edit, no build, no board
access. Nothing was written outside `experiments/420-bulk-stack/`. Every process ran on one fold at a time and peaked at
≤ 6.56 GiB RSS (`/usr/bin/time -l`: max 7.04e9 bytes).

## 0. The answer in twelve lines

1. **Every known answer reproduces with |d| = 0 before any arm number is read** (§1). That covers 414's `Kinf+sh@1.25`
   and `Kinf+sh@1.5` rows, and 418's `PP|K3`, `PP|K3+sh`, `b:h1x2|K3+sh` and `a:cd4x1.5|K3+sh` rows on all four folds.
   It also covers 418's V1 and K3 weight vectors, 418's price table through my pricer, 418's A/B/C energies, and 414's
   ballast emulation (bit-identical on 900 blocks).
2. **An instrument fix comes first.** The render energy of the *ballast-realised* bulk predicts the shipped artifacts
   2–4× better than 418's *design* route (§2). The test covers five measured ratios, including 417b (q 1.5), which was
   built tonight and never fitted. The design route misses 417b/417a by +1.1 to +1.3 %, the realised route by −0.4 to −0.7 %.
   The primary price below is therefore anchored at 417a's **measured** energy (ρ = 0.991 / 0.922 / 0.880, all three
   contexts below the mse cliff). 418's route is kept as a column; on its own rows it reproduces 418 exactly.
3. **(a) The stack H1 ×2 + CD4 ×1.5 is worth +0.0037 (k = 1) / +0.0056 (k = 1.407) vs 417a**, 3/3 cancer folds.
   - History-calibrated along H1 then CD4: **+0.0044**.
   - That is 87 % of its two halves summed: H1 ×2 +0.0025 / +0.0038, CD4 ×1.5 +0.0017 / +0.0025.
   - Energy-neutral: 0.995 / 0.991 / 0.992 of 417a.
   - Weaknesses: hek293t adds only +0.0009, and board-emulated pds falls −0.0025 (P(Δpds > 0) = 0.81 / 0.17 / 0.01).
   - Neither stack neighbour (H1 ×1.5, CD4 ×1.25) does better.
4. **(b) Sharpening dose, vs 417a, unlevelled.**
   - q 1.35: +0.0031 / +0.0038. q 1.5: +0.0063 / +0.0080. q 1.75: +0.0087 / +0.0113. q 2: +0.0084 / +0.0112.
   - All are 3/3 with bootstrap P = 1.00. q 1.75 and q 2 are flagged **LO**: energy 0.82 and 0.73 of 417a.
   - This is the one lever whose transfer is unmeasured. At kolf2-like transfer (rel ×0.371, pds ×0.332), q 1.5 is worth
     only +0.0006 and q ≥ 1.75 is negative, because each dose step lowers render energy further below the cliff.
5. **Re-levelling refunds that cost (d).** Each design's EV-optimal `VCC_BULK_CTX_SCALE` (187's cliff model, factor read
   off the realised energy map, pds measured on PLOLO):
   - worth +0.0019 on 417a itself (`A=1.002,B=1.042,C=1.067`);
   - +0.0043 on q 1.5 (`A=1.048,B=1.095,C=1.119`);
   - +0.0072 on q 1.75.
   A re-levelled q step vs re-levelled 417a is therefore ≈ 0 at zero transfer (q ≤ 1.75), +0.0031 (q 1.5) and +0.0049
   (q 1.75) at kolf2-like transfer, and +0.008 / +0.013 at cancer-like transfer. **Re-levelled sharpening does not lose on
   this model.**
6. **(b) KMAX 24 / 32 is a small rider.**
   - Bulk: +0.0005 / +0.0008 at q 1.25, +0.0009 / +0.0014 at q 1.5.
   - The capped share barely falls (20.1 → 17.2 → 15.8 % at q 1.25): most capped blocks have r_min = 0 (the sharpened
     design zeroes a gene the carriers express), and no KMAX realises those.
   - DE channel: the ballast's weight in the per-cell means rises 0.28 / 0.53 pp of 410's full weld. That is ≤ 0.0001
     overall at the weld's measured −0.0171. The Wilcoxon-rank channel is not bounded: flag it for H1.
7. **(c) Keys on the sharpened bulk, vs 417a.**
   - LODO V1 (`lodo415_lo0.3_hi3`): +0.0019 / +0.0027, 3/3.
   - Floor 0.5 (`wshare_lo0.5_hi3`): +0.0012 / +0.0013.
   - LODO at floor 0.5 (new artifact): +0.0031 / +0.0038, 3/3.
   - The floor's whole value is the β-emulated pds, which is fitted to one board point. Literal pds reads +0.0001, and on
     literal pds LODO at 0.5 is +0.0020, no better than V1.
   - The key gain holds or grows at every q (V1 +0.0019 → +0.0024 from q 1.25 to 2).
   - Control (V1 permuted): −0.022, 0/3.
8. **(d) The best combination found** is the stack + LODO V1 at floor 0.5 + q 1.5 + re-levelling
   (`A=1.045,B=1.092,C=1.117`): **+0.0178 (k = 1) / +0.0222 (k = 1.407) vs 417a**, 3/3.
   - Without re-levelling: +0.0132 / +0.0175.
   - Taking q to 1.75 reaches +0.0231, but the factors reach 1.10–1.18, realisation cos q1 is 0.97, and q 1.75's transfer
     is the least known.
9. **The single recommended next bulk arm (one variable from 417a) is the stack via `VCC_BULK_LOAD`.** The load is 417's
   overnight dump `dump_h1x2_cd4x15.npz`, made with `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51` and `wshare_lo0.3_hi3`.
   - Predicted **+0.0037 … +0.0056 (history-calibrated +0.0044)**, 3/3.
   - It does not depend on the sharpening transfer that 417a is about to measure, and it needs no builder edit.
   - **Exception:** if 417a's board row implies k_sh ≥ 0.75 (raw mse ≤ 0.703 or raw pds ≥ 0.863), 417b (q 1.5, already
     built) goes first. It is then worth +0.0037 … +0.0091 (§6).
10. **Next three bulk arms, each one variable from the previous.** Default: stack → + q 1.5 → + re-levelling. Steps
    +0.0037 / +0.0063 / +0.0049 at k = 1, **cumulative +0.0148 (k = 1) / +0.0182 (k = 1.407) / +0.0091 at kolf2-like
    sharpening**. The strong branch (417b → re-level → stack) ends at the same bulk and the same totals.
11. **The 09-23 slot-2 threshold should rise.** The PI's rule sends 417b if 417a − 410b ≥ +0.006. That line sits at
    k_sh ≈ 0.47, where unlevelled q 1.5 is worth ≈ +0.0015, less than the stack. Read the noise-free bulk members instead
    (T-read): 417b first only when raw mse ≤ 0.703 or raw pds ≥ 0.863.
12. **For October**, re-levelling is an A/B/C-only lever: its cliffs are calibrated from A/B/C board history, and D/E/F
    have no board. Composition, key and sharpening port; the cliff optimum does not, unless it can be predicted from
    controls.

---

## 1. Known answers (every script aborts on failure; all run before any arm number)

| # | row | expected | measured |
|---|---|---|---|
| KA-414 | 228 composition, Kinf + sharpening q 1.25 (arms414 + board414) and q 1.5 (extra414): cos_agg, pds, board-emulated pds lam 0.5 / 1, render A/B per frame | 1e-12 | **\|d\| = 0**, 4/4 folds (`comp420.py`) |
| KA-PP | `compose_arm` at the defaults == `sens415.compose_pp` (final rows and merged residual) | bit-identical | pass, 4/4 folds + A/B/C |
| KA-V1 | my LODO-V1 weights on PP (clip 0.3, 3) and the fold K3 weights == 418's saved `V1__PP` / `W3` | bit-identical | pass, 4/4 folds |
| KA-418 | `PP\|K3`, `PP\|K3+sh`, `b:h1x2\|K3+sh`, `a:cd4x1.5\|K3+sh` == `comp418-<fold>.json` (cos_agg, pds, pds lam 0.5 / 1, A/B per frame) | \|d\| = 0 | **\|d\| = 0**, 4/4 folds |
| KA-REAL | my row-scorer on the designed rows == `lib418.score` | 1e-12 | \|d\| = 0 on 26 rows per fold (every KA row and every 7th grid row) |
| KA-A/B/C | my design energies of PP\|none, PP\|K3, PP\|K3+sh, h1x2\|K3+sh, cd4x1.5\|K3+sh == `abc418.json`; PP constants vs the 410b meta | 1e-12; < 1e-5 | \|d\| = 0; 8.87e-6 (418's value) |
| KA-414b | `realise_k(KMAX 16)` == `ballast414.realise` on every block of PP\|K3+sh; PP\|K3 / PP\|K3+sh block stats == `ballast418.json` | bit-identical; 1e-12 | pass on 900 blocks; ≤ 9.8e-15 |
| KA-A (analysis) | 418's A/B/C-energy prices of H1 ×2 and CD4 ×1.5 vs 417a (k = 1 and 1.407, four folds) through my pricer on the d418 route | 1e-12 | **\|d\| = 0** |
| KA-LEVEL | the unscaled design re-scored inside each `--level` run == the main grid (realised K 16/24/32) | \|d\| = 0 | pass, 4/4 folds |
| KA-KEY | key artifact shas pinned. The new LODO (0.5, 3) artifact re-clips 415's key vector; re-clipping at (0.3, 3) reproduces the shipped w | bit-identical | pass |
| CTRL | the V1 key permuted over targets, on the sharpened PP bulk | negative | **−0.022 (k = 1), 0/3**; rel cos −0.36 |

## 2. Instruments

**Reused unchanged from 418:**
- the PP composition, with the composition constants as knobs;
- the 410-formula key per fold (composition-free);
- 415's LODO-V1 key, computed on PP and then **held fixed** across compositions (that is what a pinned builder artifact does). "V1own", recomputed per composition, is a secondary column and would need a new artifact;
- 417a's sharpening (`sharpen(mode="row")` on each frame's gate);
- the a23 read-back with the board-emulated pds (414's β at λ 0.5);
- the cliff-aware mse model (187's per-context identity, X solved from 410b's board raw);
- overall = (2.2138 Δpds_bd − 1.051 Δmse_raw)/6;
- k = 1 and k = 1.407 as the direction-transfer range.

**New 1: ballast-realised energies (the primary route, "r417a").**
- Every design is realised through 414's ballast emulation, on PLOLO for pds and cos and on the real A/B/C panel for
  render energy.
  - Carriers: PP unkeyed and unsharpened, as in 417a.
  - KMAX 16, plus 24 and 32 as an emulation column.
- The mse reference is 417a at its **measured** artifact energy (`cliffcheck414.json`, "the pigeon"). An arm's energy is
  417a's measured energy times (arm realised / 417a realised) on A/B/C.
- Validation against every measured artifact ratio (render-space A_c; err = predicted − measured):

{{T-val}}

**New 2: re-levelling as priced arms.**
- Per design: the EV-optimal per-context factor under 187's cliff-position sd (ctxlevel417's model; SD A/B/C 0.0025 /
  0.0004 / 0.0013), on the arm's own A_c/B_c and X.
- The builder factor is read off the realised energy map A(f)/A(1), B(f)/B(1). It was measured on A/B/C for each design
  at f = 0.94 … 1.20. The builder bound is [0.8, 1.2].
- The scaled design is then re-scored on PLOLO (`comp420.py --level`, realised K 16), so its pds and cos changes are
  measured, not assumed.

**New 3: DE leak.**
- The ballast enters a block's per-cell arithmetic mean CPM (what nmae/fid read) as
  `mean = p_de + L_K (p_realised − p_de)`, with `L_K = (K + t_c/LB)/400`.
  - L_K is 2.6 % at K 2, 6.0 % at K 16, 8.0 % at K 24 and 9.9 % at K 32.
- `abc420.py --leak` reports mean_blocks L_K‖log2 mean − log2 p_de‖₁ over genes above 5 CPM, in units of 410's full weld
  (L = 1).
- **410b sits at 0.0333, and its board DE members moved +0.0002 overall vs 395.** So a linear 410-weld scaling (the weld
  cost −0.0171) is a conservative bound.

**Transfer classes for the "class-calibrated" column:**
- composition: 418's history constants (CD4 k_mse 1.41 / k_pds 0.88; H1 0.92 / 0.70);
- key: k 1.407 (410's measured key transfer);
- sharpening: **kolf2-like**, rel ×0.371 and pds ×0.332. These are 414's kolf2 / cancer ratios at q 1.25, a pessimistic
  reading until 417a's board row measures it (T-read).

