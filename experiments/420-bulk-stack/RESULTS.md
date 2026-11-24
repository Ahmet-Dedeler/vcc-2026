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

### T-val. Energy-route validation on the shipped artifacts (A/B/C)

| ratio | measured | 418 design route (err) | realised K16 route (err) |
|---|---|---|---|
| 410b/395 | 1.0234 / 1.0222 / 1.0225 | 1.0271 / 1.0271 / 1.0270 (+0.0036 / +0.0049 / +0.0045) | 1.0247 / 1.0239 / 1.0239 (+0.0013 / +0.0017 / +0.0014) |
| 410c/410b | 1.0099 / 1.0084 / 1.0149 | 1.0173 / 1.0158 / 1.0190 (+0.0073 / +0.0074 / +0.0041) | 1.0140 / 1.0115 / 1.0148 (+0.0041 / +0.0031 / -0.0000) |
| 417a/410b | 0.9609 / 0.9578 / 0.9612 | 0.9731 / 0.9672 / 0.9706 (+0.0122 / +0.0094 / +0.0094) | 0.9587 / 0.9553 / 0.9582 (-0.0022 / -0.0025 / -0.0030) |
| 417b/417a | 0.9238 / 0.9185 / 0.9242 | 0.9373 / 0.9310 / 0.9354 (+0.0134 / +0.0125 / +0.0112) | 0.9194 / 0.9134 / 0.9177 (-0.0044 / -0.0050 / -0.0066) |
| 417b/410b | 0.8877 / 0.8797 / 0.8883 | 0.9120 / 0.9004 / 0.9079 (+0.0243 / +0.0207 / +0.0196) | 0.8815 / 0.8726 / 0.8793 (-0.0062 / -0.0071 / -0.0091) |

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


## 3. Ranked tables (all vs 417a's bulk = PP + 410 key clip(0.3, 3) + sharpening q 1.25)

**Conventions**
- Every Δoverall is (2.2138 Δpds_bd − 1.051 Δraw mse)/6 on the primary route (r417a: ballast-realised, KMAX 16, anchored
  at 417a's measured energy). Values are per fold in the order k562_gwps / orion_hct116 / orion_hek293t, with the
  cancer-fold mean and n+ = folds positive.
- "418 route" is 418's design-energy route. On 418's own rows it reproduces 418's table exactly (KA-A).
- "A/B/C energy vs 417a" is the realised render energy per context. Flags: **HI** > 1.03 and **LO** < 0.90 in any
  context.
- The builder column lists only what changes vs 417a's env: `VCC_BULK_SHARPEN_Q` 1.25, `VCC_PERT_SCALE=wshare_lo0.3_hi3`,
  PP composition.
  - `LOAD[…; PS …]` means `VCC_BULK_LOAD=<dump.npz>:<sha256>`, with a dump made by `VCC_BULK_DUMP` under 417a's env plus
    the listed composition knobs (`VCC_H1_UNION_G`, `VCC_CD4_TARGET_F`) and that `VCC_PERT_SCALE`. The key lives inside
    the dump, and sharpening applies at render time, so dumps are q-agnostic.
  - `CTX=` means `VCC_BULK_CTX_SCALE=`.
  - `*` means an artifact that is not wired into the builder yet. `lodo415_lo0.5_hi3` = `work/pert_scale_415pp_lo0.5_hi3.npz`,
    sha256 `8c82bef37d0a…` (full sha in `work/pert_scale_415pp_lo0.5_hi3.json` and `analyse420.json`). It is 415's key
    vector re-clipped at (0.5, 3), i.e. w = max(w_lo0.3, 0.5).

The full 150-design × 4-fold grid (6 compositions × 4 keys × 6 q, plus own-key and control rows) is in
`work/analyse420.log` ("ALL DESIGNS") and `work/analyse420.json`.

### 3.1 The one-variable steps from 417a (the ranking that decides the next arm)

### T-one. One-variable steps from 417a with class-calibrated transfer and bootstrap

| step | k=1 per fold | k=1 mean n+ | k=1.407 mean | class-cal. | Δpds_bd (emulated) / literal | ov with literal pds | bootstrap P(Δrel>0)/P(Δpds_bd>0) k562 · hct116 · hek293t | kolf2 k=1 |
|---|---|---|---:|---:|---|---:|---|---:|
| q 1.25 -> 1.35 | +0.0031 / +0.0028 / +0.0033 | **+0.00306** 3/3 | +0.00384 | +0.00051 | +0.00549 / +0.00558 | +0.00310 | 1.00/1.00 · 1.00/1.00 · 1.00/1.00 | +0.0005 |
| q 1.25 -> 1.5 (= 417b) | +0.0062 / +0.0058 / +0.0068 | **+0.00626** 3/3 | +0.00796 | +0.00060 | +0.01227 / +0.01274 | +0.00644 | 1.00/1.00 · 1.00/1.00 · 1.00/1.00 | +0.0003 |
| q 1.25 -> 1.75 | +0.0079 / +0.0080 / +0.0101 | **+0.00868** 3/3 | +0.01130 | -0.00051 | +0.02086 / +0.02172 | +0.00900 | 1.00/1.00 · 1.00/1.00 · 1.00/1.00 | -0.0013 |
| q 1.25 -> 2.0 | +0.0070 / +0.0079 / +0.0102 | **+0.00838** 3/3 | +0.01124 | -0.00261 | +0.02662 / +0.02742 | +0.00867 | 1.00/1.00 · 1.00/1.00 · 1.00/1.00 | -0.0043 |
| key K3 -> V1 (lodo415_lo0.3_hi3) | +0.0029 / +0.0005 / +0.0024 | **+0.00192** 3/3 | +0.00265 | +0.00265 | -0.00007 / +0.00002 | +0.00195 | 1.00/0.46 · 0.95/0.07 · 1.00/0.83 | +0.0042 |
| key K3 -> K5_3 (wshare_lo0.5_hi3) | +0.0009 / +0.0010 / +0.0018 | **+0.00122** 3/3 | +0.00129 | +0.00129 | +0.00280 / +0.00012 | +0.00023 | 1.00/0.97 · 0.79/0.97 · 0.85/0.97 | -0.0007 |
| key K3 -> V1_5_3 (new artifact) | +0.0033 / +0.0019 / +0.0039 | **+0.00305** 3/3 | +0.00379 | +0.00379 | +0.00303 / +0.00015 | +0.00199 | 1.00/0.88 · 0.96/0.95 · 1.00/0.98 | +0.0032 |
| comp H1 x2 | +0.0037 / +0.0030 / +0.0009 | **+0.00254** 3/3 | +0.00379 | +0.00252 | -0.00218 / -0.00207 | +0.00258 | 1.00/0.59 · 1.00/0.06 · 1.00/0.02 | +0.0165 |
| comp CD4 x1.5 | +0.0030 / +0.0018 / +0.0004 | **+0.00172** 3/3 | +0.00252 | +0.00255 | -0.00070 / -0.00166 | +0.00137 | 1.00/0.83 · 1.00/0.35 · 1.00/0.12 | +0.0011 |
| comp stack H1x2+CD4x1.5 | +0.0059 / +0.0044 / +0.0009 | **+0.00369** 3/3 | +0.00563 | +0.00369 (H1-then-CD4 class path +0.00441) | -0.00248 / -0.00329 | +0.00340 | 1.00/0.81 · 1.00/0.17 · 1.00/0.01 | +0.0160 |
| comp stack H1x1.5+CD4x1.5 | +0.0047 / +0.0034 / +0.0008 | **+0.00299** 3/3 | +0.00443 | +0.00299 | -0.00135 / -0.00205 | +0.00274 | 1.00/0.86 · 1.00/0.26 · 1.00/0.04 | +0.0091 |
| comp stack H1x2+CD4x1.25 | +0.0044 / +0.0033 / +0.0006 | **+0.00277** 3/3 | +0.00443 | +0.00277 | -0.00194 / -0.00246 | +0.00258 | 1.00/0.81 · 1.00/0.13 · 1.00/0.01 | +0.0156 |
| CONTROL V1 permuted | -0.0244 / -0.0196 / -0.0227 | **-0.02225** 0/3 | -0.03124 | -0.03124 | -0.00074 / -0.00055 | -0.02218 | 0.00/0.87 · 0.00/0.35 · 0.00/0.26 | -0.0276 |

Reading T-one:
- **Class-calibrated.** Composition steps use 418's history constants. The stack uses H1 then CD4, each with its own
  family constants. Key steps use k 1.407. Sharpening uses the kolf2-like transfer, which is pessimistic.
- **The composition steps carry a real board-emulated pds cost.** H1 ×2: P(Δpds_bd > 0) = 0.06 / 0.02 on the Orion
  folds. It is paid for by mse. kolf2 reads the H1 arms 5–10× higher; that fold is biased (iPSC truth vs the hESC source,
  418 §3).
- **Keys.** The floor-0.5 variants gain only through the emulated pds (literal +0.0001 / +0.0002).

### T-a. The stack and its halves (vs 417a; q 1.25, key K3)

| design | Δrel cos | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | **k=1 mean** n+ | Δov k=1.407 per fold | **k=1.407 mean** n+ | 418 route k=1 / 1.407 | A/B/C energy vs 417a | kolf2 k=1 | builder change vs 417a |
|---|---:|---:|---|---|---|---|---|---|---|---:|---|
| stack\|K3+sh1.25 | +0.0787 | -0.00248 | -0.0263 / -0.0373 | +0.0059 / +0.0044 / +0.0009 | **+0.00369** 3/3 | +0.0080 / +0.0066 / +0.0023 | **+0.00563** 3/3 | +0.00379 / +0.00568 | 0.995/0.991/0.992  | +0.0160 | `LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3]` |
| h1x2\|K3+sh1.25 | +0.0504 | -0.00218 | -0.0191 / -0.0262 | +0.0037 / +0.0030 / +0.0009 | **+0.00254** 3/3 | +0.0050 / +0.0046 / +0.0017 | **+0.00379** 3/3 | +0.00198 / +0.00321 | 1.009/1.012/1.011  | +0.0165 | `LOAD[H1_G .56; PS wshare_lo0.3_hi3]` |
| cd4x1.5\|K3+sh1.25 | +0.0321 | -0.00070 | -0.0113 / -0.0158 | +0.0030 / +0.0018 / +0.0004 | **+0.00172** 3/3 | +0.0039 / +0.0026 / +0.0010 | **+0.00252** 3/3 | +0.00148 / +0.00224 | 1.002/0.996/0.997  | +0.0011 | `LOAD[CD4_F .51; PS wshare_lo0.3_hi3]` |
| stack_h1x1.5\|K3+sh1.25 | +0.0586 | -0.00135 | -0.0199 / -0.0282 | +0.0047 / +0.0034 / +0.0008 | **+0.00299** 3/3 | +0.0063 / +0.0051 / +0.0019 | **+0.00443** 3/3 | +0.00297 / +0.00438 | 0.999/0.994/0.994  | +0.0091 | `LOAD[H1_G .42, CD4_F .51; PS wshare_lo0.3_hi3]` |
| stack_cd4x1.25\|K3+sh1.25 | +0.0680 | -0.00194 | -0.0199 / -0.0294 | +0.0044 / +0.0033 / +0.0006 | **+0.00277** 3/3 | +0.0062 / +0.0053 / +0.0018 | **+0.00443** 3/3 | +0.00292 / +0.00455 | 0.979/0.978/0.978  | +0.0156 | `LOAD[H1_G .56, CD4_F .425; PS wshare_lo0.3_hi3]` |

### T-b. Sharpening dose on PP|K3 (and PP|V1), vs 417a

| design | Δrel cos | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | **k=1 mean** n+ | Δov k=1.407 per fold | **k=1.407 mean** n+ | 418 route k=1 / 1.407 | A/B/C energy vs 417a | kolf2 k=1 | builder change vs 417a |
|---|---:|---:|---|---|---|---|---|---|---|---:|---|
| PP\|K3+sh1.35 | +0.0318 | +0.00549 | -0.0059 / -0.0103 | +0.0031 / +0.0028 / +0.0033 | **+0.00306** 3/3 | +0.0039 / +0.0035 / +0.0040 | **+0.00384** 3/3 | +0.00382 / +0.00470 | 0.971/0.969/0.971  | +0.0005 | `Q=1.35` |
| PP\|K3+sh1.5 | +0.0721 | +0.01227 | -0.0099 / -0.0196 | +0.0062 / +0.0058 / +0.0068 | **+0.00626** 3/3 | +0.0081 / +0.0073 / +0.0085 | **+0.00796** 3/3 | +0.00782 / +0.00978 | 0.919/0.913/0.918  | +0.0003 | `Q=1.5` |
| PP\|K3+sh1.75 | +0.1173 | +0.02086 | -0.0056 / -0.0206 | +0.0079 / +0.0080 / +0.0101 | **+0.00868** 3/3 | +0.0108 / +0.0104 / +0.0127 | **+0.01130** 3/3 | +0.01150 / +0.01457 | 0.825/0.815/0.821 **LO** | -0.0013 | `Q=1.75` |
| PP\|K3+sh2 | +0.1359 | +0.02662 | +0.0083 / -0.0081 | +0.0070 / +0.0079 / +0.0102 | **+0.00838** 3/3 | +0.0101 / +0.0106 / +0.0130 | **+0.01124** 3/3 | +0.01225 / +0.01566 | 0.734/0.725/0.728 **LO** | -0.0043 | `Q=2` |
| PP\|V1+sh1.5 | +0.1019 | +0.01270 | -0.0211 / -0.0349 | +0.0090 / +0.0066 / +0.0095 | **+0.00838** 3/3 | +0.0120 / +0.0086 / +0.0118 | **+0.01080** 3/3 | +0.00991 / +0.01260 | 0.926/0.921/0.924  | +0.0047 | `Q=1.5 · PS=lodo415_lo0.3_hi3` |
| PP\|V1+sh1.75 | +0.1454 | +0.02188 | -0.0166 / -0.0353 | +0.0110 / +0.0092 / +0.0127 | **+0.01098** 3/3 | +0.0150 / +0.0119 / +0.0158 | **+0.01425** 3/3 | +0.01378 / +0.01753 | 0.835/0.826/0.831 **LO** | +0.0033 | `Q=1.75 · PS=lodo415_lo0.3_hi3` |

### T-c. Key variants on the sharpened bulk (q 1.25), vs 417a

| design | Δrel cos | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | **k=1 mean** n+ | Δov k=1.407 per fold | **k=1.407 mean** n+ | 418 route k=1 / 1.407 | A/B/C energy vs 417a | kolf2 k=1 | builder change vs 417a |
|---|---:|---:|---|---|---|---|---|---|---|---:|---|
| PP\|V1+sh1.25 | +0.0298 | -0.00007 | -0.0111 / -0.0153 | +0.0029 / +0.0005 / +0.0024 | **+0.00192** 3/3 | +0.0040 / +0.0009 / +0.0031 | **+0.00265** 3/3 | +0.00163 / +0.00235 | 1.004/1.005/1.005  | +0.0042 | `PS=lodo415_lo0.3_hi3` |
| PP\|K5_3+sh1.25 | +0.0030 | +0.00280 | -0.0011 / -0.0015 | +0.0009 / +0.0010 / +0.0018 | **+0.00122** 3/3 | +0.0010 / +0.0010 / +0.0019 | **+0.00129** 3/3 | +0.00120 / +0.00126 | 1.000/1.000/1.000  | -0.0007 | `PS=wshare_lo0.5_hi3` |
| PP\|V1_5_3+sh1.25 | +0.0296 | +0.00303 | -0.0110 / -0.0152 | +0.0033 / +0.0019 / +0.0039 | **+0.00305** 3/3 | +0.0044 / +0.0024 / +0.0046 | **+0.00379** 3/3 | +0.00275 / +0.00346 | 1.005/1.005/1.005  | +0.0032 | `PS=lodo415_lo0.5_hi3*` |
| CTRL:PP\|V1shuf+sh1.25 | -0.3637 | -0.00074 | +0.1254 / +0.1768 | -0.0244 / -0.0196 / -0.0227 | **-0.02225** 0/3 | -0.0347 / -0.0275 / -0.0316 | **-0.03124** 0/3 | -0.02296 / -0.03219 | 1.004/1.005/1.005  | -0.0276 | - |

### 3.2 (d) Best combinations

Top of the whole grid (no re-levelling). The q ≥ 1.75 rows are all **LO**, so read them with §4 and §6.

### T-d. Best combinations found (top 15 of the 150-design grid by k=1 mean), vs 417a (no re-levelling)

| design | Δrel cos | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | **k=1 mean** n+ | Δov k=1.407 per fold | **k=1.407 mean** n+ | 418 route k=1 / 1.407 | A/B/C energy vs 417a | kolf2 k=1 | builder change vs 417a |
|---|---:|---:|---|---|---|---|---|---|---|---:|---|
| stack\|V1_5_3+sh1.75 | +0.2195 | +0.02313 | -0.0415 / -0.0698 | +0.0163 / +0.0150 / +0.0161 | **+0.01580** 3/3 | +0.0221 / +0.0198 / +0.0204 | **+0.02075** 3/3 | +0.01823 / +0.02357 | 0.842/0.832/0.835 **LO** | +0.0185 | `Q=1.75 · LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack\|V1_5_3+sh2 | +0.2298 | +0.02946 | -0.0261 / -0.0543 | +0.0149 / +0.0147 / +0.0167 | **+0.01544** 3/3 | +0.0206 / +0.0195 / +0.0210 | **+0.02037** 3/3 | +0.01877 / +0.02419 | 0.762/0.754/0.752 **LO** | +0.0149 | `Q=2 · LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack_h1x1.5\|V1_5_3+sh1.75 | +0.2031 | +0.02372 | -0.0363 / -0.0625 | +0.0154 / +0.0140 / +0.0159 | **+0.01511** 3/3 | +0.0208 / +0.0183 / +0.0199 | **+0.01969** 3/3 | +0.01756 / +0.02254 | 0.843/0.831/0.834 **LO** | +0.0116 | `Q=1.75 · LOAD[H1_G .42, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack_cd4x1.25\|V1_5_3+sh1.75 | +0.2095 | +0.02408 | -0.0354 / -0.0622 | +0.0153 / +0.0140 / +0.0159 | **+0.01509** 3/3 | +0.0208 / +0.0185 / +0.0200 | **+0.01978** 3/3 | +0.01758 / +0.02270 | 0.828/0.820/0.823 **LO** | +0.0181 | `Q=1.75 · LOAD[H1_G .56, CD4_F .425; PS lodo415_lo0.5_hi3*]` |
| h1x2\|V1_5_3+sh1.75 | +0.1911 | +0.02420 | -0.0349 / -0.0597 | +0.0151 / +0.0138 / +0.0162 | **+0.01504** 3/3 | +0.0202 / +0.0179 / +0.0200 | **+0.01938** 3/3 | +0.01774 / +0.02255 | 0.851/0.846/0.849 **LO** | +0.0192 | `Q=1.75 · LOAD[H1_G .56; PS lodo415_lo0.5_hi3*]` |
| stack_cd4x1.25\|V1_5_3+sh2 | +0.2211 | +0.03051 | -0.0207 / -0.0476 | +0.0142 / +0.0139 / +0.0165 | **+0.01488** 3/3 | +0.0197 / +0.0185 / +0.0206 | **+0.01959** 3/3 | +0.01836 / +0.02359 | 0.749/0.742/0.742 **LO** | +0.0146 | `Q=2 · LOAD[H1_G .56, CD4_F .425; PS lodo415_lo0.5_hi3*]` |
| h1x2\|V1_5_3+sh2 | +0.2041 | +0.03054 | -0.0203 / -0.0455 | +0.0143 / +0.0138 / +0.0164 | **+0.01483** 3/3 | +0.0195 / +0.0180 / +0.0203 | **+0.01924** 3/3 | +0.01854 / +0.02351 | 0.768/0.763/0.765 **LO** | +0.0154 | `Q=2 · LOAD[H1_G .56; PS lodo415_lo0.5_hi3*]` |
| stack_h1x1.5\|V1_5_3+sh2 | +0.2154 | +0.02969 | -0.0213 / -0.0476 | +0.0141 / +0.0137 / +0.0163 | **+0.01468** 3/3 | +0.0195 / +0.0180 / +0.0204 | **+0.01929** 3/3 | +0.01803 / +0.02314 | 0.761/0.751/0.750 **LO** | +0.0082 | `Q=2 · LOAD[H1_G .42, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack\|V1+sh1.75 | +0.2181 | +0.01989 | -0.0404 / -0.0684 | +0.0155 / +0.0137 / +0.0140 | **+0.01441** 3/3 | +0.0212 / +0.0184 / +0.0184 | **+0.01933** 3/3 | +0.01689 / +0.02220 | 0.839/0.829/0.832 **LO** | +0.0196 | `Q=1.75 · LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.3_hi3]` |
| stack\|V1+sh2 | +0.2280 | +0.02619 | -0.0245 / -0.0524 | +0.0140 / +0.0136 / +0.0143 | **+0.01396** 3/3 | +0.0196 / +0.0183 / +0.0187 | **+0.01884** 3/3 | +0.01733 / +0.02270 | 0.758/0.749/0.748 **LO** | +0.0158 | `Q=2 · LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.3_hi3]` |
| cd4x1.5\|V1_5_3+sh1.75 | +0.1786 | +0.02372 | -0.0287 / -0.0517 | +0.0140 / +0.0122 / +0.0151 | **+0.01378** 3/3 | +0.0189 / +0.0158 / +0.0188 | **+0.01781** 3/3 | +0.01622 / +0.02065 | 0.845/0.831/0.835 **LO** | +0.0035 | `Q=1.75 · LOAD[CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack_h1x1.5\|V1+sh1.75 | +0.2018 | +0.02056 | -0.0352 / -0.0611 | +0.0147 / +0.0128 / +0.0138 | **+0.01375** 3/3 | +0.0200 / +0.0170 / +0.0179 | **+0.01830** 3/3 | +0.01625 / +0.02120 | 0.840/0.828/0.831 **LO** | +0.0126 | `Q=1.75 · LOAD[H1_G .42, CD4_F .51; PS lodo415_lo0.3_hi3]` |
| stack_cd4x1.25\|V1+sh1.75 | +0.2082 | +0.02090 | -0.0343 / -0.0609 | +0.0145 / +0.0128 / +0.0138 | **+0.01372** 3/3 | +0.0199 / +0.0173 / +0.0179 | **+0.01838** 3/3 | +0.01627 / +0.02135 | 0.825/0.817/0.820 **LO** | +0.0192 | `Q=1.75 · LOAD[H1_G .56, CD4_F .425; PS lodo415_lo0.3_hi3]` |
| stack\|K5_3+sh1.75 | +0.1919 | +0.02212 | -0.0317 / -0.0564 | +0.0133 / +0.0138 / +0.0141 | **+0.01372** 3/3 | +0.0180 / +0.0182 / +0.0179 | **+0.01804** 3/3 | +0.01626 / +0.02095 | 0.838/0.827/0.831 **LO** | +0.0136 | `Q=1.75 · LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.5_hi3]` |
| h1x2\|V1+sh1.75 | +0.1900 | +0.02083 | -0.0338 / -0.0584 | +0.0143 / +0.0125 / +0.0140 | **+0.01360** 3/3 | +0.0193 / +0.0166 / +0.0178 | **+0.01791** 3/3 | +0.01635 / +0.02113 | 0.848/0.842/0.846 **LO** | +0.0202 | `Q=1.75 · LOAD[H1_G .56; PS lodo415_lo0.3_hi3]` |

The same ranking restricted to q ≤ 1.5 with the energy inside [0.90, 1.03] of 417a:

### T-d2. Best combinations restricted to q <= 1.5 (energy within [0.90, 1.03] of 417a)

| design | Δrel cos | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | **k=1 mean** n+ | Δov k=1.407 per fold | **k=1.407 mean** n+ | 418 route k=1 / 1.407 | A/B/C energy vs 417a | kolf2 k=1 | builder change vs 417a |
|---|---:|---:|---|---|---|---|---|---|---|---:|---|
| stack\|V1_5_3+sh1.5 | +0.1807 | +0.01347 | -0.0470 / -0.0714 | +0.0149 / +0.0125 / +0.0122 | **+0.01320** 3/3 | +0.0199 / +0.0167 / +0.0158 | **+0.01748** 3/3 | +0.01449 / +0.01898 | 0.925/0.918/0.920  | +0.0205 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack_h1x1.5\|V1_5_3+sh1.5 | +0.1624 | +0.01434 | -0.0414 / -0.0634 | +0.0139 / +0.0115 / +0.0122 | **+0.01254** 3/3 | +0.0185 / +0.0151 / +0.0156 | **+0.01639** 3/3 | +0.01383 / +0.01789 | 0.928/0.919/0.922  | +0.0133 | `Q=1.5 · LOAD[H1_G .42, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| stack\|V1own+sh1.5 | +0.1840 | +0.01087 | -0.0478 / -0.0726 | +0.0142 / +0.0119 / +0.0110 | **+0.01238** 3/3 | +0.0193 / +0.0162 / +0.0148 | **+0.01673** 3/3 | +0.01373 / +0.01831 | 0.923/0.916/0.918  | +0.0215 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS <own-composition LODO>*]` |
| stack_cd4x1.25\|V1_5_3+sh1.5 | +0.1697 | +0.01414 | -0.0404 / -0.0632 | +0.0135 / +0.0115 / +0.0118 | **+0.01230** 3/3 | +0.0182 / +0.0154 / +0.0153 | **+0.01629** 3/3 | +0.01368 / +0.01790 | 0.910/0.905/0.908  | +0.0201 | `Q=1.5 · LOAD[H1_G .56, CD4_F .425; PS lodo415_lo0.5_hi3*]` |
| h1x2\|V1_5_3+sh1.5 | +0.1506 | +0.01427 | -0.0398 / -0.0603 | +0.0131 / +0.0113 / +0.0123 | **+0.01223** 3/3 | +0.0174 / +0.0147 / +0.0154 | **+0.01583** 3/3 | +0.01373 / +0.01758 | 0.937/0.935/0.938  | +0.0211 | `Q=1.5 · LOAD[H1_G .56; PS lodo415_lo0.5_hi3*]` |
| stack\|V1+sh1.5 | +0.1800 | +0.01045 | -0.0464 / -0.0708 | +0.0143 / +0.0112 / +0.0104 | **+0.01199** 3/3 | +0.0193 / +0.0153 / +0.0141 | **+0.01625** 3/3 | +0.01332 / +0.01779 | 0.923/0.916/0.918  | +0.0214 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.3_hi3]` |
| stack_h1x1.5\|V1own+sh1.5 | +0.1652 | +0.01172 | -0.0420 / -0.0644 | +0.0134 / +0.0106 / +0.0110 | **+0.01169** 3/3 | +0.0181 / +0.0143 / +0.0145 | **+0.01560** 3/3 | +0.01300 / +0.01714 | 0.927/0.918/0.920  | +0.0139 | `Q=1.5 · LOAD[H1_G .42, CD4_F .51; PS <own-composition LODO>*]` |
| stack_cd4x1.25\|V1own+sh1.5 | +0.1710 | +0.01149 | -0.0406 / -0.0635 | +0.0130 / +0.0108 / +0.0103 | **+0.01134** 3/3 | +0.0177 / +0.0148 / +0.0137 | **+0.01536** 3/3 | +0.01278 / +0.01705 | 0.908/0.903/0.906  | +0.0212 | `Q=1.5 · LOAD[H1_G .56, CD4_F .425; PS <own-composition LODO>*]` |
| stack_h1x1.5\|V1+sh1.5 | +0.1617 | +0.01135 | -0.0408 / -0.0627 | +0.0134 / +0.0103 / +0.0104 | **+0.01134** 3/3 | +0.0179 / +0.0138 / +0.0138 | **+0.01517** 3/3 | +0.01266 / +0.01671 | 0.927/0.918/0.920  | +0.0142 | `Q=1.5 · LOAD[H1_G .42, CD4_F .51; PS lodo415_lo0.3_hi3]` |
| cd4x1.5\|V1_5_3+sh1.5 | +0.1364 | +0.01466 | -0.0333 / -0.0518 | +0.0124 / +0.0098 / +0.0116 | **+0.01124** 3/3 | +0.0163 / +0.0126 / +0.0145 | **+0.01448** 3/3 | +0.01251 / +0.01596 | 0.932/0.921/0.924  | +0.0049 | `Q=1.5 · LOAD[CD4_F .51; PS lodo415_lo0.5_hi3*]` |

With re-levelling on top: each design at its own EV-optimal `VCC_BULK_CTX_SCALE`, with the pds and cos of the scaled design
measured on PLOLO (realised, K 16):

### T-lev. Per-context re-levelling on top (VCC_BULK_CTX_SCALE, EV-optimal; pds measured on PLOLO for the scaled design)

| design | factors (k=1) | render scale A/B/C | level-only Δov per fold | level-only mean | Δraw mse | Δpds_bd | design + level vs 417a, k=1 per fold | mean n+ | k=1.407 mean |
|---|---|---|---|---:|---:|---:|---|---|---:|
| PP\|K3+sh1.25 | `A=1.002,B=1.042,C=1.067` | 1.0019/1.0406/1.0648 | +0.0019 / +0.0019 / +0.0020 | +0.00194 | -0.0107 | +0.00018 | +0.0019 / +0.0019 / +0.0020 | **+0.00194** 3/3 | +0.00186 |
| PP\|K3+sh1.35 | `A=1.018,B=1.060,C=1.084` | 1.0173/1.0574/1.0804 | +0.0028 / +0.0028 / +0.0028 | +0.00277 | -0.0155 | +0.00016 | +0.0058 / +0.0056 / +0.0060 | **+0.00583** 3/3 | +0.00652 |
| PP\|K3+sh1.5 | `A=1.048,B=1.095,C=1.119` | 1.0453/1.0895/1.1119 | +0.0043 / +0.0043 / +0.0044 | +0.00433 | -0.0244 | +0.00016 | +0.0105 / +0.0101 / +0.0112 | **+0.01059** 3/3 | +0.01218 |
| PP\|K3+sh1.75 | `A=1.114,B=1.168,C=1.194` | 1.1036/1.1529/1.1754 | +0.0069 / +0.0071 / +0.0076 | +0.00721 | -0.0410 | +0.00009 | +0.0149 / +0.0152 / +0.0177 | **+0.01589** 3/3 | +0.01837 |
| PP\|K3+sh2 | `A=1.197,B=1.200,C=1.200` | 1.1700/1.1742/1.1730 | +0.0078 / +0.0080 / +0.0087 | +0.00819 | -0.0464 | +0.00017 | +0.0148 / +0.0159 / +0.0190 | **+0.01656** 3/3 | +0.02149 |
| PP\|V1_5_3+sh1.25 | `A=1.000,B=1.039,C=1.064` | 1.0000/1.0378/1.0620 | +0.0018 / +0.0018 / +0.0018 | +0.00180 | -0.0103 | -0.00000 | +0.0051 / +0.0037 / +0.0057 | **+0.00485** 3/3 | +0.00554 |
| PP\|V1_5_3+sh1.5 | `A=1.043,B=1.089,C=1.113` | 1.0406/1.0838/1.1064 | +0.0041 / +0.0041 / +0.0040 | +0.00405 | -0.0236 | -0.00023 | +0.0137 / +0.0120 / +0.0154 | **+0.01369** 3/3 | +0.01607 |
| PP\|V1_5_3+sh1.75 | `A=1.104,B=1.157,C=1.183` | 1.0947/1.1430/1.1658 | +0.0068 / +0.0067 / +0.0067 | +0.00670 | -0.0397 | -0.00069 | +0.0186 / +0.0169 / +0.0217 | **+0.01907** 3/3 | +0.02232 |
| stack\|V1_5_3+sh1.25 | `A=1.004,B=1.046,C=1.070` | 1.0039/1.0446/1.0678 | +0.0023 / +0.0023 / +0.0023 | +0.00230 | -0.0131 | +0.00003 | +0.0117 / +0.0082 / +0.0062 | **+0.00873** 3/3 | +0.01147 |
| stack\|V1_5_3+sh1.5 | `A=1.045,B=1.092,C=1.117` | 1.0426/1.0868/1.1101 | +0.0047 / +0.0047 / +0.0045 | +0.00464 | -0.0268 | -0.00018 | +0.0196 / +0.0172 / +0.0167 | **+0.01784** 3/3 | +0.02223 |
| stack\|V1_5_3+sh1.75 | `A=1.101,B=1.155,C=1.183` | 1.0922/1.1412/1.1657 | +0.0075 / +0.0074 / +0.0071 | +0.00731 | -0.0430 | -0.00061 | +0.0238 / +0.0224 / +0.0232 | **+0.02311** 3/3 | +0.02826 |
| stack\|K3+sh1.25 | `A=1.005,B=1.047,C=1.071` | 1.0049/1.0456/1.0688 | +0.0024 / +0.0023 / +0.0023 | +0.00234 | -0.0130 | +0.00016 | +0.0082 / +0.0067 / +0.0032 | **+0.00603** 3/3 | +0.00796 |
| stack\|K3+sh1.5 | `A=1.047,B=1.095,C=1.119` | 1.0444/1.0897/1.1120 | +0.0047 / +0.0047 / +0.0048 | +0.00470 | -0.0267 | +0.00010 | +0.0159 / +0.0151 / +0.0129 | **+0.01464** 3/3 | +0.01817 |

- **The level-only gain is nearly identical across folds by construction.** It is an A/B/C-energy computation; only X
  differs by fold. So its 3/3 is not independent evidence. Its risk sits in 187's cliff positions (claim_c), not in the
  folds.
- **The scaled designs cost a little direction.** Realised cos falls 0.3–0.9 % and pds_bd −0.0002 to −0.0007 at
  factors ≥ 1.09. The render is nonlinear and more blocks cap. Both costs are already inside the numbers.
- **q 2 cannot be fully levelled.** Its factors hit the builder's 1.2 cap.

**417b re-levelled from its MEASURED artifact energy** (`cliffcheck414.json`, "my-thesis-advisor-is-a-cat" =
`submissions/417b-…-sharpen150`; ρ = 0.916 / 0.847 / 0.813):
- the EV-optimal factors are `VCC_BULK_CTX_SCALE=A=1.045,B=1.092,C=1.115`, at k = 1 and at k = 1.407;
- worth −0.0264 raw mse, **+0.0046 overall** (the mse part; PLOLO puts the pds change at +0.0002);
- these, not the design-route factors, are the ones to build with.

## 4. (b) Sharpening dose: ballast realisability, KMAX 24 / 32, and the DE channel

### T-kmax. Ballast realisability and the KMAX 24 / 32 emulation (r417a route, cancer means; DE leak in units of 410's full weld)

| design | capped K16 / K24 / K32 | cos(realised, design) q1 (K16) | Δov vs 417a k=1: K16 / K24 / K32 | k=1.407: K16 / K24 / K32 | DE leak K16 / K24 / K32 | mean K (K16 / K32) |
|---|---|---:|---|---|---|---|
| PP\|K3 | 10.3 / 8.4 / 7.8 % | 0.9992 | -0.01278 / -0.01269 / -0.01265 | -0.01510 / -0.01497 / -0.01491 | 0.0333 / 0.0358 / 0.0378 | 4.00 / 5.30 |
| PP\|K3+sh1.25 | 20.1 / 17.2 / 15.8 % | 0.9912 | +0.00000 / +0.00047 / +0.00075 | +0.00000 / +0.00059 / +0.00094 | 0.0318 / 0.0346 / 0.0371 | 5.46 / 8.12 |
| PP\|K3+sh1.35 | 23.9 / 20.8 / 19.9 % | 0.9887 | +0.00306 / +0.00367 / +0.00402 | +0.00384 / +0.00459 / +0.00502 | 0.0318 / 0.0347 / 0.0373 | 6.03 / 9.30 |
| PP\|K3+sh1.5 | 28.8 / 26.1 / 24.3 % | 0.9803 | +0.00626 / +0.00712 / +0.00762 | +0.00796 / +0.00902 / +0.00962 | 0.0336 / 0.0371 / 0.0403 | 6.78 / 10.85 |
| PP\|K3+sh1.75 | 36.2 / 33.1 / 31.4 % | 0.9721 | +0.00868 / +0.00984 / +0.01053 | +0.01130 / +0.01270 / +0.01353 | 0.0405 / 0.0461 / 0.0514 | 7.91 / 13.11 |
| PP\|K3+sh2 | 42.0 / 39.2 / 37.8 % | 0.9644 | +0.00838 / +0.00978 / +0.01061 | +0.01124 / +0.01290 / +0.01388 | 0.0485 / 0.0567 / 0.0644 | 8.73 / 14.89 |
| PP\|V1+sh1.25 | 19.8 / 16.8 / 15.8 % | 0.9912 | +0.00192 / +0.00237 / +0.00248 | +0.00265 / +0.00323 / +0.00340 | 0.0321 / 0.0348 / 0.0371 | 5.64 / 8.26 |
| PP\|V1+sh1.5 | 30.8 / 27.2 / 25.2 % | 0.9821 | +0.00838 / +0.00925 / +0.00973 | +0.01080 / +0.01187 / +0.01246 | 0.0347 / 0.0385 / 0.0419 | 7.18 / 11.44 |
| stack\|K3+sh1.25 | 18.6 / 16.3 / 15.3 % | 0.9913 | +0.00369 / +0.00413 / +0.00438 | +0.00563 / +0.00617 / +0.00648 | 0.0353 / 0.0385 / 0.0414 | 5.48 / 8.03 |
| stack\|V1+sh1.5 | 30.3 / 26.4 / 25.0 % | 0.9814 | +0.01199 / +0.01278 / +0.01323 | +0.01625 / +0.01722 / +0.01777 | 0.0390 / 0.0435 / 0.0476 | 7.25 / 11.42 |

**Bulk side.**
- Raising KMAX recovers part of each dose's realisation loss: +0.0005 / +0.0008 at q 1.25 and +0.0009 / +0.0014 at q 1.5
  (K 24 / K 32).
- It cannot recover more. Among capped blocks, most have r_min = 0: the sharpened (or keyed) design sends a gene to zero
  that the carriers express at > 1e-9. No finite K realises that; the builder clips the gene at the carriers' total.
- The capped share falls only 20.1 → 17.2 → 15.8 % (q 1.25) and 28.8 → 26.1 → 24.3 % (q 1.5).

**DE side (the carrier count drops to 376 / 368 on capped blocks only).** Two channels:
1. **Arithmetic mean CPM** (nmae's log2FC, fid's direction).
   - Here the ballast leaks L_K (p_realised − p_de) into the mean. The leak for 417a's design is 0.0318 of 410's full
     weld at K 16, 0.0346 at K 24 and 0.0371 at K 32. For q 1.5 it is 0.0336 / 0.0371 / 0.0403.
   - 410b sat at 0.0333, and the board moved its DE members **+0.0002** overall vs 395.
   - At the weld's measured −0.0171 overall, a linear bound puts the extra leak of KMAX 32 at ≤ −0.0001 (q ≤ 1.5).
   - Re-levelling barely moves the leak: 417a re-levelled 0.0328, q 1.5 re-levelled 0.0346, q 1.75 re-levelled 0.0413
     (`abc420_leak.json`).
2. **Wilcoxon ranks** (fid/jac/reach calls).
   - The ballast share of the rank sample rises from mean K/400 = 1.4 % (q 1.25, K 16) to 2.0 % (K 32). On capped blocks
     it goes from 4 % to 8 % of cells.
   - The emitter sizes its call set on 400 − K carriers, so the carriers themselves adapt. What 8 % of deep,
     bulk-composition cells do to the U statistic of marginal genes is **not bounded here**.

**Verdict.** KMAX 32 is worth ≤ +0.001 on the bulk at the doses that matter, against a mean-channel cost ≤ 0.0001. It
needs a builder constant change and an H1 check of the rank channel (416/419 territory). It is a rider for a q ≥ 1.5 arm,
not an arm.

## 5. Stacking arithmetic: what adds and what doesn't

### T-int. Interactions: each lever vs the same q without it (Δov k=1 / k=1.407, cancer means; all 3/3 unless marked)

| lever | q 1.25 | q 1.35 | q 1.5 | q 1.75 | q 2 |
|---|---|---|---|---|---|
| V1 vs K3 | +0.00192 / +0.00265 | +0.00203 / +0.00277 | +0.00216 / +0.00288 | +0.00235 / +0.00301 | +0.00237 / +0.00295 |
| K5_3 vs K3 | +0.00122 / +0.00129 | +0.00127 / +0.00135 | +0.00148 / +0.00158 | +0.00175 / +0.00188 | +0.00182 / +0.00195 |
| V1_5_3 vs K3 | +0.00305 / +0.00379 | +0.00316 / +0.00390 | +0.00342 / +0.00415 | +0.00375 / +0.00443 | +0.00381 / +0.00442 |
| V1_5_3 vs V1 | +0.00114 / +0.00113 | +0.00113 / +0.00114 | +0.00126 / +0.00127 | +0.00140 / +0.00143 | +0.00144 / +0.00147 |
| stack vs PP (K3) | +0.00369 / +0.00563 | +0.00372 / +0.00562 | +0.00375 / +0.00556 | +0.00356 / +0.00514 | +0.00337 / +0.00473 |
| stack vs PP (V1_5_3) | +0.00340 / +0.00541 | +0.00353 / +0.00551 | +0.00367 / +0.00557 | +0.00357 / +0.00527 | +0.00347 / +0.00494 |
| H1x2 vs PP (K3) | +0.00254 / +0.00379 | +0.00260 / +0.00383 | +0.00268 / +0.00384 | +0.00278 / +0.00379 | +0.00287 / +0.00373 |
| CD4x1.5 vs PP (K3) | +0.00172 / +0.00252 | +0.00171 / +0.00249 | +0.00169 / +0.00243 | +0.00149 / +0.00214 | +0.00129 / +0.00186 |
| stack vs H1x2 (K3) | +0.00117 / +0.00188 | +0.00114 / +0.00183 | +0.00109 / +0.00175 | +0.00079 / +0.00137 | +0.00051 / +0.00101 |

- **The key and the composition are orthogonal to the dose.**
  - The V1 key's gain grows with q: +0.0019 → +0.0024.
  - The stack's gain is flat to q 1.5 (+0.0037 → +0.0038) and falls only at q ≥ 1.75.
  - The CD4 half fades with dose (+0.0017 → +0.0013 at q 2). The H1 half does not.
- **The stack is 87 % additive at q 1.25.** Its CD4 increment on top of H1 ×2 is +0.0012 (3/3).
- **The key composes with the stack at ≈ 95 %.** Stack + V1_5_3 is +0.0064, against +0.0037 + +0.0031.
- **Along any path, cumulative step prices match the direct end-vs-417a price to ±0.0002** (the PATHS block of
  `analyse420.log`). Stacks here price as sums, with the stated redundancy.

## 6. Recommendation

### 6.1 The single next bulk arm, one variable from 417a

**The stack via `VCC_BULK_LOAD`:**

    417a's env + VCC_BULK_LOAD=experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz:<its sha256>

- **Dump.** 417's overnight chain makes it with `VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51 VCC_PERT_SCALE=wshare_lo0.3_hi3`.
- **Gate.** 417's KA-LOAD (same-env dump reproduces 417a's X) must PASS first.
- **Config diff vs 417a:** exactly `bulk_load`.
- **Predicted vs 417a:**
  - Δoverall **+0.0037 (k = 1) / +0.0056 (k = 1.407) / +0.0044 (history-calibrated)**, 3/3 cancer folds;
  - Δraw mse −0.026 … −0.037, Δpds_bd −0.0025;
  - render energy 0.995 / 0.991 / 0.992 of 417a;
  - realisation: 18.6 % of blocks at K 16 (417a: 20.1 %), cos q1 0.991;
  - DE members unchanged by construction; leak 0.0353 of the weld against 417a's 0.0318.
- **Why this arm and not the others:**
  - it is the largest step whose value does not depend on the sharpening transfer (unmeasured until 417a reads);
  - it needs no builder edit;
  - its class has two board-calibrated steps (418 §5).
- **Where it could disappoint:** hek293t (+0.0009) and a significant pds cost on the Orion folds.

**Exception, decided by 417a's board row:** if the bulk members imply k_sh ≥ 0.75, send 417b (q 1.5, built) first. That
means raw mse ≤ 0.703 or raw pds ≥ 0.863 (T-read below), and then q 1.5 is worth +0.0037 … +0.0091 vs 417a. Between
k_sh 0.55 and 0.75 the two arms are tied (q 1.5 +0.0021 … +0.0037); below 0.55 the stack wins.

### 6.2 The next three bulk arms (each one variable from the previous)

### T-seq. Sequences (each arm one variable from the previous)

**STRONG (417a reads k_sh >= 0.75)** — cumulative vs 417a: k=1 +0.01479, k=1.407 +0.01821, class-calibrated +0.00914

| step | Δrel | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | k=1 mean n+ | k=1.407 | class-cal. | kolf2 k=1 | builder |
|---|---:|---:|---|---|---|---:|---:|---:|---|
| 417b: q 1.5 | +0.0721 | +0.01227 | -0.0099 / -0.0196 | +0.0062 / +0.0058 / +0.0068 | **+0.00626** 3/3 | +0.00796 | +0.00060 | +0.0003 | `Q=1.5` |
| + levelling | -0.0092 | +0.00016 | -0.0248 / -0.0234 | +0.0044 / +0.0044 / +0.0045 | **+0.00441** 3/3 | +0.00415 | +0.00441 | +0.0044 | `Q=1.5 · CTX=A=1.048,B=1.095,C=1.119` |
| + stack (LOAD, factors re-derived) | +0.0715 | -0.00208 | -0.0279 / -0.0392 | +0.0055 / +0.0051 / +0.0018 | **+0.00413** 3/3 | +0.00609 | +0.00413 | +0.0169 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3] · CTX=A=1.047,B=1.095,C=1.119` |

**WEAK (417a reads k_sh < 0.75)** — cumulative vs 417a: k=1 +0.00883, k=1.407 +0.01148, class-calibrated +0.00963

| step | Δrel | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | k=1 mean n+ | k=1.407 | class-cal. | kolf2 k=1 | builder |
|---|---:|---:|---|---|---|---:|---:|---:|---|
| stack (LOAD) | +0.0787 | -0.00248 | -0.0263 / -0.0373 | +0.0059 / +0.0044 / +0.0009 | **+0.00369** 3/3 | +0.00563 | +0.00369 | +0.0160 | `LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3]` |
| + key V1_5_3 (new dump) | +0.0296 | +0.00204 | -0.0115 / -0.0161 | +0.0036 / +0.0016 / +0.0032 | **+0.00277** 3/3 | +0.00357 | +0.00357 | +0.0042 | `LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.5_hi3*]` |
| + levelling | -0.0029 | +0.00003 | -0.0134 / -0.0130 | +0.0024 / +0.0024 / +0.0023 | **+0.00237** 3/3 | +0.00228 | +0.00237 | +0.0029 | `LOAD[H1_G .56, CD4_F .51; PS lodo415_lo0.5_hi3*] · CTX=A=1.004,B=1.046,C=1.070` |

**WEAK-alt (existing key only)** — cumulative vs 417a: k=1 +0.01483, k=1.407 +0.01818, class-calibrated +0.00914

| step | Δrel | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | k=1 mean n+ | k=1.407 | class-cal. | kolf2 k=1 | builder |
|---|---:|---:|---|---|---|---:|---:|---:|---|
| stack (LOAD) | +0.0787 | -0.00248 | -0.0263 / -0.0373 | +0.0059 / +0.0044 / +0.0009 | **+0.00369** 3/3 | +0.00563 | +0.00369 | +0.0160 | `LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3]` |
| + q 1.5 | +0.0637 | +0.01273 | -0.0090 / -0.0184 | +0.0054 / +0.0061 / +0.0073 | **+0.00627** 3/3 | +0.00792 | +0.00058 | +0.0001 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3]` |
| + levelling | -0.0080 | +0.00010 | -0.0276 / -0.0262 | +0.0048 / +0.0049 / +0.0049 | **+0.00487** 3/3 | +0.00463 | +0.00487 | +0.0057 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3] · CTX=A=1.047,B=1.095,C=1.119` |

**LEVEL-FIRST** — cumulative vs 417a: k=1 +0.01471, k=1.407 +0.01827, class-calibrated +0.00911

| step | Δrel | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | k=1 mean n+ | k=1.407 | class-cal. | kolf2 k=1 | builder |
|---|---:|---:|---|---|---|---:|---:|---:|---|
| levelling on 417a | -0.0030 | +0.00018 | -0.0107 / -0.0102 | +0.0019 / +0.0019 / +0.0020 | **+0.00194** 3/3 | +0.00186 | +0.00194 | +0.0019 | `CTX=A=1.002,B=1.042,C=1.067` |
| + q 1.5 (factors re-derived) | +0.0655 | +0.01224 | -0.0236 / -0.0331 | +0.0086 / +0.0081 / +0.0092 | **+0.00864** 3/3 | +0.01031 | +0.00305 | +0.0028 | `Q=1.5 · CTX=A=1.048,B=1.095,C=1.119` |
| + stack (LOAD, factors re-derived) | +0.0715 | -0.00208 | -0.0279 / -0.0392 | +0.0055 / +0.0051 / +0.0018 | **+0.00413** 3/3 | +0.00609 | +0.00413 | +0.0169 | `Q=1.5 · LOAD[H1_G .56, CD4_F .51; PS wshare_lo0.3_hi3] · CTX=A=1.047,B=1.095,C=1.119` |

- **Default (WEAK-alt): stack → + q 1.5 → + re-levelling.** One env key per step: `bulk_load`, then `bulk_sharpen_q`,
  then `bulk_ctx_scale`.
  - Predicted steps +0.0037 / +0.0063 / +0.0049 at k = 1.
  - **Cumulative +0.0148 (k = 1), +0.0182 (k = 1.407), +0.0091 at kolf2-like sharpening.**
  - At zero sharpening transfer, the path still ends ≈ +0.006: the re-levelling refunds q's amplitude cost.
  - Re-derive step 3's factors from step 2's measured artifact energy with the same EV rule. The design route gives
    `A=1.047,B=1.095,C=1.119` for this bulk; for 417b's own bulk the measured route gives `A=1.045,B=1.092,C=1.115`.
- **If 417a reads k_sh ≥ 0.75 (STRONG): 417b → + re-levelling → + stack.** It reaches the same bulk and the same totals,
  collecting the transfer-dependent step first.
- **Fourth arm, either branch: the key.**
  - LODO V1 at floor 0.5 adds +0.0033 … +0.0042 on the stacked q 1.5 bulk. It needs the new artifact wired in and a new
    dump.
  - `lodo415_lo0.3_hi3` (already wired) adds +0.0021 … +0.0029, and does not lean on the β emulation.
- **Not recommended now.**
  - q ≥ 1.75: LO energy, factors ≥ 1.10 to relevel, and the least-known transfer. Take it only after 417b's read confirms
    k_sh.
  - KMAX 24/32 as a standalone arm.
  - The own-composition LODO key: it adds ≤ +0.0004 over the fixed PP key and would need a new artifact per composition.

### 6.3 Reading 417a's board row (for the PI at ~00:30Z)

### T-read. Reading 417a's board row: implied sharpening transfer, and the q steps' worth at that transfer

PLOLO (realised, cancer mean) 417a vs 410b: rel cos +0.1013, board-emulated pds +0.01713; measured render amplitude 417a/410b 0.9798. k_pds is tied to k_sh at the kolf2 ratio (x0.895).

| k_sh | predicted 417a raw mse | predicted 417a raw pds | Δoverall 417a vs 410b (bulk members) | q 1.35 | q 1.5 | q 1.75 | q 2 | levelled q 1.35 | levelled q 1.5 | levelled q 1.75 | levelled q 2 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.000 | 0.7263 | 0.8512 | +0.0008 | -0.00086 | -0.00245 | -0.00546 | -0.00848 | +0.00001 | +0.00003 | +0.00002 | -0.00201 |
| 0.200 | 0.7201 | 0.8542 | +0.0030 | -0.00012 | -0.00081 | -0.00279 | -0.00531 | +0.00074 | +0.00166 | +0.00264 | +0.00111 |
| 0.371 | 0.7148 | 0.8568 | +0.0049 | +0.00051 | +0.00060 | -0.00051 | -0.00261 | +0.00137 | +0.00305 | +0.00489 | +0.00378 |
| 0.550 | 0.7093 | 0.8596 | +0.0069 | +0.00118 | +0.00208 | +0.00187 | +0.00022 | +0.00203 | +0.00451 | +0.00724 | +0.00657 |
| 0.750 | 0.7031 | 0.8627 | +0.0091 | +0.00192 | +0.00372 | +0.00454 | +0.00339 | +0.00276 | +0.00613 | +0.00986 | +0.00969 |
| 1.000 | 0.6953 | 0.8665 | +0.0119 | +0.00285 | +0.00578 | +0.00787 | +0.00734 | +0.00368 | +0.00817 | +0.01314 | +0.01359 |
| 1.200 | 0.6891 | 0.8695 | +0.0141 | +0.00359 | +0.00743 | +0.01054 | +0.01051 | +0.00441 | +0.00980 | +0.01577 | +0.01671 |
| 1.407 | 0.6827 | 0.8727 | +0.0164 | +0.00436 | +0.00914 | +0.01330 | +0.01379 | +0.00518 | +0.01148 | +0.01849 | +0.01994 |

- **The pds read is the cleaner one.** k_pds = (raw pds − 0.851156) / 0.01713 does not go through the cliff model, and
  k_sh ≈ k_pds / 0.895.
- **The mse read goes through 187's claims.** k_sh = (solve_X(417a measured A/B, raw mse) / (X_410b × 0.9798) − 1) / 0.1013.
- **If the two reads disagree by more than ~0.3 in k_sh,** the cliff model (and therefore every re-levelling number here)
  is what failed. Re-price the levelled rows before building a CTX arm.

## 7. Limits

- **Sharpening transfer is unmeasured.** Every q row above k = 1 / 1.407 is an assumption. 414's kolf2 reads it at 0.37
  (rel) and 0.33 (pds) of the cancer folds. 417a is the first board point.
- **Composition and key transfer rest on thin calibration.** 418's history constants come from one board step per family.
  k_key 1.407 comes from 410's single key step. FINDINGS H.5 says to treat a once-measured constant as a lower bound on
  magnitude.
- **Re-levelling is a model prediction.** It rests on 187's claim_c (retrodicted 410/410b to +0.0002 raw) and on the
  realised-energy map. It is A/B/C-specific; nothing here calibrates D/E/F cliffs.
- **The energy route has a small residual.** The realised route over-predicts sharpening's energy drop by 0.2–0.9 %
  (417a/410b, 417b/417a, 417b/410b). The error is in the conservative direction for the unlevelled q rows (it
  overstates their amplitude cost). For the levelled rows the factors are ~0.3 % high: use the measured-artifact factors
  (§3.2).
- **The β emulation (board-emulated pds) is fitted to one board point** (410b − 395). The floor-0.5 key's value is
  entirely inside it.
- **The ballast emulation is 414's.**
  - It emulates the design profile, not ratio mode on the carriers' delivered composition. It skips the jackknife boost
    and the 999k trim, and uses a constant carrier total.
  - The energy validation above is its end-to-end check.
- **kolf2 cannot vote**, and it over-reads H1 arms. The cancer folds are two-donor folds, blind to the normrel family
  (418 §6); normrel was not varied here.
- **Parallelism.** Up to three fold processes (each one fold, ≤ 6.6 GiB) and one A/B/C process ran concurrently. The
  memory gate (≥ 40 % free, ≥ 20 GiB free+inactive) was checked before each and never went below 67 % / 22.7 GiB.

## 8. Files and reproduction

```
cd experiments/420-bulk-stack            # everything writes only into ./work (PYTHONDONTWRITEBYTECODE=1 everywhere)
uv run python abc420.py                  # A/B/C: KAs, energy-route validation rows, 146 designs x {design, realised K16/24/32}
bash chain420.sh k562_gwps orion_hct116 orion_hek293t kolf2   # comp420.py <fold>: KAs + the 150-design grid (+own-key, control)
uv run python analyse420.py              # prices, KA-A, validation, tables, re-levelling factors -> work/level_specs.txt
uv run python abc420.py --level <designs...>   # realised energy maps at f = 0.94..1.20 -> work/abc420_level.json
bash chain420L.sh k562_gwps orion_hct116 orion_hek293t kolf2  # comp420.py <fold> --level <specs>: scaled designs' pds/cos
uv run python abc420.py --leak <designs...>    # DE leak in units of 410's weld -> work/abc420_leak.json
uv run python analyse420.py && uv run python table420.py > work/table420.md && uv run python assemble420.py
```

- Libraries: `lib420.py` has the KMAX-parameterised `realise_k`, the realised rows, the key artifacts and the new LODO
  (0.5, 3) artifact. 418's `lib418.py` / `analyse418.py` are imported read-only.
- Outputs:
  - `work/comp420-<fold>.json` / `arrays420-<fold>.npz` (per-target arrays for the bootstrap);
  - `work/comp420L-<fold>.json`;
  - `work/abc420.json`, `abc420_level.json`, `abc420_leak.json`;
  - `work/analyse420.{json,log}`, `work/table420.md`;
  - `work/pert_scale_415pp_lo0.5_hi3.{npz,json}` (candidate artifact).
- Peak RSS (`/usr/bin/time -l`, max resident):
  - comp420: 6.65 / 7.00 / 7.04 / 5.26 GB (k562 / hct116 / hek293t / kolf2);
  - comp420 --level: ≤ 6.56 GiB (internal peak);
  - abc420: 6.09 GB; --level / --leak ≤ 5.43 GB.
- **Writes.**
  - Nothing outside this directory. Fold caches are APFS clones of 418's in `work/plolo/`.
  - `find -newer` over every imported directory lists only other agents' files: the builder and its `.pyc` from 417's
    builds, and 414's `cliffcheck414.json`, which 417 updated with 417b.
