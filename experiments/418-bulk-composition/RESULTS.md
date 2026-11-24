# 418: re-tuning the BULK composition under ballast. The bulk wants MORE residual (H1 ×2, CD4 ×1.5), not less; every other shipped constant is at its bulk optimum

Status: **COMPLETE** (2026-09-22/23, ~21:15Z). PLOLO on four folds, the real A/B/C panel, and a retrodiction of ten
historical board steps. No vcc, no cell-eval2, no builder edit, no build. Nothing written outside
`experiments/418-bulk-composition/`. Every process peaked at ≤ 7.09 GB RSS (one fold per process).

## 0. The answer in ten lines

1. **All known answers reproduce with |d| = 0 before any arm number is read** (§1): 414's `base`, `Kinf` and
   `Kinf+sh@1.25` rows, 415's PP rows (composition constants, base aggregate cosine, K410 and V1_PP gains), and
   bit-identity of my composition code with 415's `compose_pp` on four folds and on A/B/C (KA-PP 8.9e-6 vs the 410b meta).
2. **Best one-variable bulk change vs 417a: H1-union g ×2** (`VCC_H1_UNION_G` 0.28 → 0.56).
   - Predicted Δoverall **+0.0020 (k = 1) to +0.0032 (k = 1.407)**, positive on 3/3 cancer folds at both k.
   - History-calibrated **+0.0020**; amplitude-matched +0.0022.
   - It costs pds: Δpds_bd −0.0021; the bootstrap puts P(Δpds > 0) at 0.06 and 0.02 on the two Orion folds.
   - It pays in mse: Δraw −0.016 to −0.023.
   - Real-panel render energy rises 1.0–1.2 %; the price already includes that cost.
   - The H1 dose peaks at ×2. ×2.5 and ×3 fall to 2/3 folds as energy climbs 3–5 %.
3. **Second best: CD4 fraction ×1.5** (`VCC_CD4_TARGET_F` 0.34 → 0.51).
   - Predicted **+0.0015 to +0.0022** vs 417a, 3/3 folds at both k.
   - History-calibrated **+0.0023**, the highest of any arm.
   - Δpds_bd −0.0008 (every per-fold CI spans 0). Render energy neutral on A/B/C (1.001/0.995/0.995).
   - The CD4 dose peaks at ×1.5. ×1.75 and ×2 drop to 2/3 folds.
   - Of the two, this is the more robust arm (§4).
4. **The mandate's premise inverts.** Under the current composition (SVD-5 deflation plus the key), the bulk's pds does
   not turn over near f ≈ 0.20. It peaks at ×1 to ×1.25 of the shipped CD4 dose, and mse keeps improving to ×1.5.
   Lowering either residual loses 3/3 folds: CD4 ×0.5 −0.0055, H1 ×0 −0.0057.
5. **Deflation k, partial g and the normrel exponent are at or next to their bulk optimum.**
   - k 0/3/8: −0.0038 / −0.0008 / −0.0009.
   - g 0/0.25/0.75/1: −0.0009 / −0.0003 / −0.0007 / −0.0022.
   - normrel −0.75 / off: −0.0001 / −0.0008.
   - All values are vs 417a at k = 1. No re-tune.
6. **Own idea f1 (the key scales the donor base only; the residual goes in unkeyed) is strongly negative:** −0.0077 to −0.0104,
   0/3 folds. The residual must ride the key.
7. **Own idea f2 (normrel PLAIN instead of cap) is an amplitude lever and nothing more.**
   - It gives +0.0008 vs 410b, entirely from −1.9 % render energy. Its amplitude-matched value is +0.0001.
   - Vs 417a it is −0.0001: sharpening already takes that energy out.
8. **Instrument blind spot: normrel.** Cancer folds have two source donors, so they cannot emulate three-donor per-target
   weighting.
   - They get the sign of 2 of the 4 historical normrel board steps wrong. kolf2, with three production donors, gets
     4 of 4 right on pds.
   - kolf2 reads normrel −3 at +0.0004 to +0.0005, and the board history (pre-key) read it at +0.0010. The cancer-fold
     −0.0011 for normrel −3 should not be believed (§6).
9. **Retrodiction validates the dose families** (§5): the board's overall (pds + mse) sign for CD4 dose, H1, deflation
   and the three partial steps is right 5/6 on the A/B/C-energy route (6/6 on the PLOLO-energy route). The one miss is
   partial 0.34 → 0.5: board +0.0002, predicted −0.0001.
   - Across all ten steps the overall sign is right 8/10, corr +0.87.
   - Fitted class transfer: CD4 dose k_mse 1.41 / k_pds 0.88; H1 0.92 / 0.70.
   - So the k = 1…1.407 range brackets the dose arms.
10. **Both winners need a builder change to be bulk-only.** The knobs act on the one shared matrix that the carriers and
    the ballast render both use (§4).
    - Ballast realises both bulk designs at the same level as 410b/417a: K = 16 on 8.9–10.3 % of blocks unsharpened and
      19–20 % sharpened, realised-vs-design cos q1 ≥ 0.991.
    - Stacking the two was not measured (one-variable mandate); it is the obvious next PLOLO row.

---

## 1. Known answers (all run before any arm; every script aborts on failure)

| # | row | expected | measured |
|---|---|---|---|
| KA-414 | 228's composition: `base`, `Kinf`, `Kinf+sh@1.25` through my code path (a23 read-back, board-emulated pds lam 0.5/1, render A_c/B_c) | `arms414-<fold>.json`, `board414-<fold>.json`, 1e-12 | **\|d\| = 0**, 4/4 folds (`comp418.py`) |
| KA-415a | my `pre_normrel(cap, −1.5)` == `sens415.Donors.pre` | bit-identical | `np.array_equal` True, 4 folds + A/B/C |
| KA-415b | my `compose_arm` at the defaults == `sens415.compose_pp` (final row, merged residual, post-add s) | bit-identical | **max \|d\| = 0**, 4 folds + A/B/C |
| KA-415c | composition constants (medians, operator preserve, f_cd4/f_h1/f_partial, s) | `sens415.json` constants_pp | max \|d\| = 0, 4/4 folds |
| KA-415d | effect-level base aggregate cosine; K410 and V1_PP gains at clip (0.3,3) and (0.3,∞) | `sens415.json`, 1e-12 | \|d\| ≤ 2.2e-16, 4/4 folds (e.g. k562 V1_PP clip3 +0.306111857) |
| KA-PP | A/B/C PP constants vs the 410b build meta; deflation energy and U_k sha | 415's 8.9e-6 | **8.87e-6**; 0.198660465061 and `f06e0262…` exact; key npz sha `d7f76a22…` asserted |
| KA-X | X_410b and X_395 solved from the board raws on 414's measured A_c/B_c | `analyse414.json` Xb | \|d\| < 1e-12 |
| KA-main | the dose-continuation runs' reference rows (PP\|none, PP\|K3, PP\|K3+sh, PP\|V1, PP\|V1+sh, the shuffled-key control) | `comp418-<fold>.json` | **\|d\| = 0**, 4/4 folds |
| KA-ballast | the ballast emulation's design render == abc418's | 1e-9 | pass. PP\|K3 emulates K median 2 / max 16, as 410b's meta (`ballast_K` median 2, max 16) |
| CTRL | the fold key permuted over targets | negative | rel cos −0.36 / −0.36 / −0.38 / −0.49 (k562/hct116/hek293t/kolf2) |

The 410b analog on PLOLO (PP|K3 vs PP|none) reads rel +0.235 / +0.282 / +0.240. 417a's sharpening on PP|K3 reads rel
+0.106 to +0.122, Δpds_bd +0.016 to +0.019, overall +0.015 to +0.017 (cancer, PLOLO energies, k = 1.407). That is 414's
recommendation re-read on the shipped composition.

## 2. Instruments (414's, plus two additions)

- **Composition.** PP is 415's reconstruction of the shipped bulk (normrel cap −1.5 → whitening γ 0.65 → CD4 SVD-5 @
  target f 0.34 + H1 union @ g 0.28 + CD4 partial @ g 0.5 → one post-add preserve, panel zeroed). Each arm changes one
  constant, exactly as the builder's env knob would.
- **Key.** The fold's 410-formula key at clip (0.3, 3) is held fixed and applied with `apply_pert_scale`'s arithmetic.
  The 410 key is composition-free, so only the Frobenius constant c moves. The secondary column is each arm's own LODO-V1
  key (415's construction), recomputed on its composition with builder `cols=None`.
- **Sharpening.** 417a's is sign|h|^1.25 on each frame's gate, with the row's gated norm restored. On panel-zeroed rows
  this equals `patch_sharpen417.bulk_sharpen`.
- **Read-back.** The a23 read-back in the challenge frames gives rel aggregate cosine, pds, render A_c/B_c, and the
  board-emulated pds (414's β at λ 0.5).
- **mse.** 414's cliff-aware route (187's model; the reference sits at 410b's measured A_c/B_c and board-calibrated X).
  - **k is defined on the arm's OWN relative cosine change vs its reference:** X_arm = X_ref (1 + k r) √(Ā ratio). k = 1
    and 1.407 are reported as the range, as mandated. The retro's fitted k has the same definition.
- **Overall.** Δoverall = (2.2138 Δpds_bd − 1.051 Δmse_raw)/6.
- **Addition 1: real-panel render energies (`abc418.py`, the primary route).** Every arm is built on the real A/B/C panel
  (300 targets, three production donors, production operator) with the shipped 410 weights, then rendered in each frame.
  - The resulting A_c/B_c ratios replace the PLOLO-fold ratios in the cliff route; this is 414's move for sharpening
    (`sharp410c.py`).
  - It matters. The 2-donor folds get the sign of the energy change wrong for H1 ×2 (PLOLO −1 %, A/B/C +1 %),
    normrel −3 (+2 % vs flat) and normrel plain (+0.6 % vs −1.9 %).
  - The PLOLO-energy numbers are kept as a column.
- **Addition 2: board-history retrodiction (§5).** Ten historical one-variable composition steps have board pds/mse.
  Each is re-built on PLOLO without a key, and its board Δ is predicted from the base entry's own board raw.
  - This validates the instrument for this arm class.
  - It fits class transfer constants (k_pds, k_mse) per family. They drive the "hist" column (a: the CD4-dose step;
    b: the H1 step; c: the deflation step; d: the three partial steps, pooled).
- **Amp-matched column.** The arm is rescaled per context to the reference's render energy, so the column reads its
  direction value only.

## 3. Ranked tables (cancer folds decide; kolf2 is reported, no vote)

Columns, per fold in the order k562 / hct116 / hek293t:
- The **primary** price is the A/B/C-energy route with board-emulated pds, at k = 1 and at k = 1.407.
- "hist" is the history-calibrated price; "amp-matched" is direction only.
- "PLOLO-energy k=1" is 414's original route. "A/B/C render A ratio" is the real-panel energy vs the reference.
- **kolf2 on the H1 arms is biased.** kolf2 is an iPSC line, the same kind of cell as the H1 source, so it reads H1 about
  8× higher than the cancer folds do. Ignore it there.

#### 3.1 Ranked vs 417a's bulk (PP + fixed key + sharpening q 1.25)

Δoverall = (2.2138 Δpds_bd − 1.051 Δmse)/6. A/B/C-energy cliff route; board-emulated pds. Per fold k562 / hct116 / hek293t.

| # | arm | Δrel cos | Δpds_bd | Δmse (k=1 / 1.407) | **k=1** per fold | **k=1** mean [range] +n | **k=1.407** per fold | **k=1.407** mean [range] +n | hist | amp-matched | PLOLO-energy k=1 | A/B/C render A ratio | kolf2 k=1 / 1.407 |
|---:|---|---:|---:|---|---|---|---|---|---:|---:|---:|---|---|
| 1 | H1 g ×2 (`VCC_H1_UNION_G=0.56`) | +0.0487 | -0.00214 | -0.0158 / -0.0228 | +0.0031 / +0.0025 / +0.0004 | **+0.00198** [+0.0004, +0.0031] 3/3 | +0.0044 / +0.0040 / +0.0012 | **+0.00321** [+0.0012, +0.0044] 3/3 | +0.00196 | +0.00221 | +0.00208 | 1.0096/1.0113/1.0115 | +0.0159 / +0.0215 |
| 2 | H1 g ×2.5 (`VCC_H1_UNION_G=0.7`) | +0.0624 | -0.00359 | -0.0177 / -0.0268 | +0.0034 / +0.0026 / -0.0006 | **+0.00177** [-0.0006, +0.0034] 2/3 | +0.0051 / +0.0046 / +0.0004 | **+0.00336** [+0.0004, +0.0051] 3/3 | +0.00184 | +0.00252 | +0.00235 | 1.0287/1.0321/1.0322 | +0.0221 / +0.0302 |
| 3 | CD4 f ×1.5 (`VCC_CD4_TARGET_F=0.51`) | +0.0301 | -0.00076 | -0.0101 / -0.0144 | +0.0028 / +0.0016 / +0.0001 | **+0.00148** [+0.0001, +0.0028] 3/3 | +0.0036 / +0.0024 / +0.0007 | **+0.00224** [+0.0007, +0.0036] 3/3 | +0.00227 | +0.00163 | +0.00134 | 1.0010/0.9949/0.9953 | +0.0011 / +0.0014 |
| 4 | H1 g ×1.5 (`VCC_H1_UNION_G=0.42`) | +0.0278 | -0.00074 | -0.0099 / -0.0139 | +0.0020 / +0.0018 / +0.0006 | **+0.00146** [+0.0006, +0.0020] 3/3 | +0.0027 / +0.0026 / +0.0011 | **+0.00216** [+0.0011, +0.0027] 3/3 | +0.00140 | +0.00145 | +0.00138 | 0.9999/1.0004/1.0004 | +0.0087 / +0.0116 |
| 5 | H1 g ×3 (`VCC_H1_UNION_G=0.84`) | +0.0704 | -0.00501 | -0.0182 / -0.0286 | +0.0034 / +0.0023 / -0.0017 | **+0.00135** [-0.0017, +0.0034] 2/3 | +0.0054 / +0.0047 / -0.0007 | **+0.00316** [-0.0007, +0.0054] 2/3 | +0.00152 | +0.00249 | +0.00167 | 1.0439/1.0489/1.0489 | +0.0274 / +0.0377 |
| 6 | CD4 f ×1.25 (`VCC_CD4_TARGET_F=0.425`) | +0.0186 | +0.00002 | -0.0065 / -0.0091 | +0.0018 / +0.0012 / +0.0005 | **+0.00114** [+0.0005, +0.0018] 3/3 | +0.0023 / +0.0017 / +0.0009 | **+0.00161** [+0.0009, +0.0023] 3/3 | +0.00161 | +0.00118 | +0.00105 | 0.9996/0.9966/0.9968 | +0.0007 / +0.0009 |
| 7 | CD4 f ×1.75 (`VCC_CD4_TARGET_F=0.595`) | +0.0356 | -0.00208 | -0.0108 / -0.0159 | +0.0030 / +0.0014 / -0.0010 | **+0.00112** [-0.0010, +0.0030] 2/3 | +0.0040 / +0.0024 / -0.0004 | **+0.00201** [-0.0004, +0.0040] 2/3 | +0.00211 | +0.00151 | +0.00126 | 1.0071/0.9979/0.9987 | +0.0008 / +0.0011 |
| 8 | CD4 f ×2 (`VCC_CD4_TARGET_F=0.68`) | +0.0364 | -0.00371 | -0.0112 / -0.0164 | +0.0030 / +0.0009 / -0.0021 | **+0.00059** [-0.0021, +0.0030] 2/3 | +0.0041 / +0.0020 / -0.0015 | **+0.00151** [-0.0015, +0.0041] 2/3 | +0.00167 | +0.00098 | +0.00050 | 1.0053/0.9934/0.9945 | +0.0005 / +0.0006 |
| 9 | normrel cap exp -0.75 (`VCC_NORMREL_EXP=-0.75`) | -0.0006 | -0.00028 | -0.0000 / +0.0001 | -0.0001 / +0.0001 / -0.0004 | **-0.00010** [-0.0004, +0.0001] 1/3 | -0.0000 / +0.0001 / -0.0004 | **-0.00011** [-0.0004, +0.0001] 1/3 | – | -0.00014 | -0.00008 | 0.9924/0.9930/0.9934 | -0.0005 / -0.0005 |
| 10 | f2: normrel PLAIN −1.5 (`VCC_NORMREL_MODE=plain`) | +0.0004 | +0.00041 | +0.0016 / +0.0015 | -0.0002 / +0.0001 / -0.0002 | **-0.00012** [-0.0002, +0.0001] 1/3 | -0.0002 / +0.0001 / -0.0002 | **-0.00011** [-0.0002, +0.0001] 1/3 | – | +0.00018 | +0.00016 | 0.9815/0.9816/0.9815 | +0.0010 / +0.0010 |
| 11 | partial g 0.25 (`VCC_CD4_PARTIAL_G=0.25`) | -0.0016 | -0.00028 | +0.0009 / +0.0012 | -0.0004 / +0.0003 / -0.0008 | **-0.00027** [-0.0008, +0.0003] 1/3 | -0.0004 / +0.0003 / -0.0008 | **-0.00031** [-0.0008, +0.0003] 1/3 | -0.00036 | -0.00021 | -0.00023 | 0.9892/0.9899/0.9900 | -0.0000 / -0.0000 |
| 12 | partial g 0.75 (`VCC_CD4_PARTIAL_G=0.75`) | -0.0021 | -0.00040 | +0.0033 / +0.0036 | -0.0004 / -0.0015 / -0.0002 | **-0.00073** [-0.0015, -0.0002] 0/3 | -0.0005 / -0.0016 / -0.0003 | **-0.00078** [-0.0016, -0.0003] 0/3 | -0.00086 | -0.00027 | -0.00034 | 1.0141/1.0131/1.0124 | -0.0009 / -0.0011 |
| 13 | CD4 deflation k 3 (`VCC_CD4_DEFLATE_K=3`) | +0.0045 | -0.00261 | -0.0011 / -0.0017 | -0.0011 / -0.0001 / -0.0011 | **-0.00078** [-0.0011, -0.0001] 0/3 | -0.0009 / +0.0001 / -0.0012 | **-0.00066** [-0.0012, +0.0001] 1/3 | -0.00014 | -0.00068 | -0.00071 | 1.0033/1.0032/1.0037 | -0.0009 / -0.0007 |
| 14 | normrel cap exp 0 (`VCC_NORMREL_EXP=unset`) | -0.0055 | -0.00151 | +0.0013 / +0.0021 | -0.0002 / -0.0006 / -0.0016 | **-0.00078** [-0.0016, -0.0002] 0/3 | -0.0001 / -0.0008 / -0.0018 | **-0.00092** [-0.0018, -0.0001] 0/3 | – | -0.00089 | -0.00087 | 0.9945/0.9956/0.9964 | -0.0015 / -0.0016 |
| 15 | CD4 deflation k 8 (`VCC_CD4_DEFLATE_K=8`) | -0.0054 | -0.00147 | +0.0020 / +0.0027 | -0.0007 / +0.0003 / -0.0022 | **-0.00089** [-0.0022, +0.0003] 1/3 | -0.0009 / +0.0001 / -0.0023 | **-0.00102** [-0.0023, +0.0001] 1/3 | -0.00198 | -0.00088 | -0.00086 | 1.0005/1.0009/1.0006 | -0.0003 / -0.0005 |
| 16 | partial g 0 (`VCC_CD4_PARTIAL_G=unset`) | -0.0065 | -0.00134 | +0.0023 / +0.0033 | -0.0009 / -0.0001 / -0.0018 | **-0.00091** [-0.0018, -0.0001] 0/3 | -0.0010 / -0.0002 / -0.0020 | **-0.00107** [-0.0020, -0.0002] 0/3 | -0.00129 | -0.00090 | -0.00092 | 0.9908/0.9918/0.9919 | -0.0003 / -0.0005 |
| 17 | normrel cap exp -3 (`VCC_NORMREL_EXP=-3`) | -0.0090 | -0.00133 | +0.0035 / +0.0048 | -0.0003 / -0.0018 / -0.0013 | **-0.00111** [-0.0018, -0.0003] 0/3 | -0.0004 / -0.0020 / -0.0016 | **-0.00133** [-0.0020, -0.0004] 0/3 | – | -0.00104 | -0.00154 | 1.0008/0.9994/0.9983 | +0.0004 / +0.0002 |
| 18 | partial g 1 (`VCC_CD4_PARTIAL_G=1`) | -0.0077 | -0.00081 | +0.0106 / +0.0117 | -0.0015 / -0.0035 / -0.0014 | **-0.00215** [-0.0035, -0.0014] 0/3 | -0.0017 / -0.0037 / -0.0016 | **-0.00235** [-0.0037, -0.0016] 0/3 | -0.00262 | -0.00076 | -0.00158 | 1.0429/1.0400/1.0388 | -0.0027 / -0.0030 |
| 19 | H1 g ×0.5 (`VCC_H1_UNION_G=0.14`) | -0.0358 | -0.00013 | +0.0121 / +0.0172 | -0.0029 / -0.0026 / -0.0010 | **-0.00216** [-0.0029, -0.0010] 0/3 | -0.0038 / -0.0037 / -0.0017 | **-0.00306** [-0.0038, -0.0017] 0/3 | -0.00197 | -0.00226 | -0.00240 | 0.9952/0.9947/0.9943 | -0.0101 / -0.0133 |
| 20 | CD4 f ×0.75 (`VCC_CD4_TARGET_F=0.255`) | -0.0272 | -0.00102 | +0.0104 / +0.0143 | -0.0029 / -0.0021 / -0.0015 | **-0.00220** [-0.0029, -0.0015] 0/3 | -0.0037 / -0.0029 / -0.0021 | **-0.00289** [-0.0037, -0.0021] 0/3 | -0.00284 | -0.00208 | -0.00201 | 1.0054/1.0084/1.0078 | -0.0017 / -0.0021 |
| 21 | H1 g ×4 (`VCC_H1_UNION_G=1.12`) | +0.0729 | -0.00770 | +0.0032 / -0.0078 | -0.0006 / -0.0020 / -0.0076 | **-0.00340** [-0.0076, -0.0006] 0/3 | +0.0016 / +0.0007 / -0.0067 | **-0.00147** [-0.0067, +0.0016] 2/3 | -0.00295 | +0.00166 | -0.00693 | 1.1042/1.1132/1.1128 | +0.0331 / +0.0473 |
| 22 | CD4 deflation k 0 (`VCC_CD4_DEFLATE_K=unset`) | +0.0073 | -0.01120 | -0.0018 / -0.0029 | -0.0026 / -0.0031 / -0.0058 | **-0.00381** [-0.0058, -0.0026] 0/3 | -0.0022 / -0.0029 / -0.0058 | **-0.00363** [-0.0058, -0.0022] 0/3 | -0.00330 | -0.00367 | -0.00382 | 1.0051/1.0054/1.0059 | -0.0008 / -0.0002 |
| 23 | CD4 f ×0.5 (`VCC_CD4_TARGET_F=0.17`) | -0.0642 | -0.00378 | +0.0237 / +0.0330 | -0.0070 / -0.0051 / -0.0045 | **-0.00555** [-0.0070, -0.0045] 0/3 | -0.0087 / -0.0068 / -0.0060 | **-0.00717** [-0.0087, -0.0060] 0/3 | -0.00700 | -0.00540 | -0.00556 | 1.0072/1.0123/1.0116 | -0.0040 / -0.0050 |
| 24 | H1 g ×0 (`VCC_H1_UNION_G=unset`) | -0.0793 | -0.00245 | +0.0274 / +0.0387 | -0.0076 / -0.0064 / -0.0032 | **-0.00570** [-0.0076, -0.0032] 0/3 | -0.0096 / -0.0088 / -0.0047 | **-0.00769** [-0.0096, -0.0047] 0/3 | -0.00502 | -0.00581 | -0.00600 | 0.9962/0.9953/0.9949 | -0.0238 / -0.0304 |
| 25 | f1: key on the donor base only (no knob) | -0.1097 | -0.00163 | +0.0406 / +0.0562 | -0.0074 / -0.0098 / -0.0060 | **-0.00771** [-0.0098, -0.0060] 0/3 | -0.0100 / -0.0131 / -0.0082 | **-0.01044** [-0.0131, -0.0082] 0/3 | – | -0.00741 | -0.00737 | 0.9765/0.9819/0.9811 | -0.0152 / -0.0204 |


#### 3.2 Ranked vs 410b's bulk (PP + fixed key, no sharpening)

Δoverall = (2.2138 Δpds_bd − 1.051 Δmse)/6. A/B/C-energy cliff route; board-emulated pds. Per fold k562 / hct116 / hek293t.

| # | arm | Δrel cos | Δpds_bd | Δmse (k=1 / 1.407) | **k=1** per fold | **k=1** mean [range] +n | **k=1.407** per fold | **k=1.407** mean [range] +n | hist | amp-matched | PLOLO-energy k=1 | A/B/C render A ratio | kolf2 k=1 / 1.407 |
|---:|---|---:|---:|---|---|---|---|---|---:|---:|---:|---|---|
| 1 | H1 g ×2 (`VCC_H1_UNION_G=0.56`) | +0.0549 | -0.00239 | -0.0156 / -0.0226 | +0.0031 / +0.0023 / +0.0001 | **+0.00184** [+0.0001, +0.0031] 3/3 | +0.0044 / +0.0038 / +0.0010 | **+0.00307** [+0.0010, +0.0044] 3/3 | +0.00185 | +0.00212 | +0.00245 | 1.0090/1.0107/1.0109 | +0.0154 / +0.0210 |
| 2 | CD4 f ×1.5 (`VCC_CD4_TARGET_F=0.51`) | +0.0341 | -0.00063 | -0.0102 / -0.0145 | +0.0033 / +0.0016 / -0.0002 | **+0.00155** [-0.0002, +0.0033] 2/3 | +0.0042 / +0.0024 / +0.0004 | **+0.00231** [+0.0004, +0.0042] 3/3 | +0.00233 | +0.00168 | +0.00192 | 1.0010/0.9954/0.9959 | +0.0009 / +0.0013 |
| 3 | H1 g ×2.5 (`VCC_H1_UNION_G=0.7`) | +0.0718 | -0.00398 | -0.0172 / -0.0265 | +0.0034 / +0.0022 / -0.0010 | **+0.00154** [-0.0010, +0.0034] 2/3 | +0.0052 / +0.0043 / +0.0001 | **+0.00317** [+0.0001, +0.0052] 3/3 | +0.00165 | +0.00245 | +0.00216 | 1.0284/1.0311/1.0317 | +0.0212 / +0.0293 |
| 4 | H1 g ×1.5 (`VCC_H1_UNION_G=0.42`) | +0.0310 | -0.00096 | -0.0099 / -0.0138 | +0.0020 / +0.0016 / +0.0005 | **+0.00138** [+0.0005, +0.0020] 3/3 | +0.0027 / +0.0025 / +0.0010 | **+0.00207** [+0.0010, +0.0027] 3/3 | +0.00134 | +0.00134 | +0.00145 | 0.9993/1.0001/1.0001 | +0.0083 / +0.0113 |
| 5 | CD4 f ×1.75 (`VCC_CD4_TARGET_F=0.595`) | +0.0418 | -0.00171 | -0.0113 / -0.0166 | +0.0036 / +0.0013 / -0.0009 | **+0.00134** [-0.0009, +0.0036] 2/3 | +0.0047 / +0.0024 / -0.0003 | **+0.00227** [-0.0003, +0.0047] 2/3 | +0.00235 | +0.00173 | +0.00194 | 1.0071/0.9984/0.9992 | +0.0007 / +0.0010 |
| 6 | CD4 f ×1.25 (`VCC_CD4_TARGET_F=0.425`) | +0.0206 | +0.00023 | -0.0064 / -0.0090 | +0.0021 / +0.0012 / +0.0003 | **+0.00120** [+0.0003, +0.0021] 3/3 | +0.0026 / +0.0017 / +0.0007 | **+0.00165** [+0.0007, +0.0026] 3/3 | +0.00164 | +0.00123 | +0.00137 | 0.9997/0.9970/0.9973 | +0.0005 / +0.0008 |
| 7 | CD4 f ×2 (`VCC_CD4_TARGET_F=0.68`) | +0.0451 | -0.00328 | -0.0125 / -0.0182 | +0.0038 / +0.0011 / -0.0020 | **+0.00098** [-0.0020, +0.0038] 2/3 | +0.0051 / +0.0022 / -0.0013 | **+0.00198** [-0.0013, +0.0051] 2/3 | +0.00213 | +0.00136 | +0.00148 | 1.0051/0.9935/0.9948 | +0.0005 / +0.0008 |
| 8 | f2: normrel PLAIN −1.5 (`VCC_NORMREL_MODE=plain`) | +0.0001 | +0.00032 | -0.0038 / -0.0038 | +0.0006 / +0.0011 / +0.0007 | **+0.00079** [+0.0006, +0.0011] 3/3 | +0.0006 / +0.0011 / +0.0007 | **+0.00079** [+0.0006, +0.0011] 3/3 | – | +0.00012 | -0.00008 | 0.9813/0.9812/0.9811 | +0.0022 / +0.0022 |
| 9 | H1 g ×3 (`VCC_H1_UNION_G=0.84`) | +0.0829 | -0.00549 | -0.0150 / -0.0258 | +0.0029 / +0.0015 / -0.0026 | **+0.00061** [-0.0026, +0.0029] 2/3 | +0.0050 / +0.0039 / -0.0015 | **+0.00250** [-0.0015, +0.0050] 2/3 | +0.00082 | +0.00250 | +0.00022 | 1.0436/1.0471/1.0480 | +0.0259 / +0.0361 |
| 10 | normrel cap exp -0.75 (`VCC_NORMREL_EXP=-0.75`) | -0.0009 | -0.00014 | -0.0016 / -0.0015 | +0.0003 / +0.0003 / +0.0001 | **+0.00023** [+0.0001, +0.0003] 3/3 | +0.0003 / +0.0002 / +0.0000 | **+0.00021** [+0.0000, +0.0003] 3/3 | – | -0.00010 | +0.00008 | 0.9911/0.9917/0.9922 | -0.0005 / -0.0005 |
| 11 | partial g 0.25 (`VCC_CD4_PARTIAL_G=0.25`) | -0.0019 | -0.00019 | -0.0015 / -0.0013 | +0.0001 / +0.0008 / -0.0003 | **+0.00020** [-0.0003, +0.0008] 2/3 | +0.0001 / +0.0008 / -0.0004 | **+0.00016** [-0.0004, +0.0008] 2/3 | +0.00010 | -0.00018 | -0.00014 | 0.9898/0.9906/0.9906 | +0.0003 / +0.0003 |
| 12 | normrel cap exp 0 (`VCC_NORMREL_EXP=unset`) | -0.0060 | -0.00116 | +0.0002 / +0.0010 | +0.0001 / -0.0005 / -0.0010 | **-0.00046** [-0.0010, +0.0001] 1/3 | +0.0001 / -0.0007 / -0.0012 | **-0.00060** [-0.0012, +0.0001] 1/3 | – | -0.00076 | -0.00050 | 0.9923/0.9932/0.9940 | -0.0018 / -0.0019 |
| 13 | partial g 0 (`VCC_CD4_PARTIAL_G=unset`) | -0.0069 | -0.00141 | +0.0004 / +0.0013 | -0.0005 / +0.0003 / -0.0016 | **-0.00059** [-0.0016, +0.0003] 1/3 | -0.0006 / +0.0002 / -0.0018 | **-0.00074** [-0.0018, +0.0002] 1/3 | -0.00095 | -0.00090 | -0.00084 | 0.9918/0.9930/0.9929 | -0.0002 / -0.0004 |
| 14 | partial g 0.75 (`VCC_CD4_PARTIAL_G=0.75`) | -0.0012 | -0.00027 | +0.0032 / +0.0034 | -0.0004 / -0.0014 / -0.0002 | **-0.00066** [-0.0014, -0.0002] 0/3 | -0.0004 / -0.0014 / -0.0002 | **-0.00069** [-0.0014, -0.0002] 0/3 | -0.00073 | -0.00016 | -0.00041 | 1.0132/1.0120/1.0119 | -0.0008 / -0.0009 |
| 15 | CD4 deflation k 8 (`VCC_CD4_DEFLATE_K=8`) | -0.0064 | -0.00106 | +0.0022 / +0.0030 | -0.0007 / +0.0005 / -0.0020 | **-0.00077** [-0.0020, +0.0005] 1/3 | -0.0009 / +0.0003 / -0.0021 | **-0.00091** [-0.0021, +0.0003] 1/3 | -0.00189 | -0.00074 | -0.00060 | 1.0010/1.0014/1.0012 | -0.0002 / -0.0004 |
| 16 | CD4 deflation k 3 (`VCC_CD4_DEFLATE_K=3`) | +0.0056 | -0.00279 | -0.0010 / -0.0017 | -0.0011 / +0.0001 / -0.0016 | **-0.00085** [-0.0016, +0.0001] 1/3 | -0.0009 / +0.0003 / -0.0017 | **-0.00073** [-0.0017, +0.0003] 1/3 | -0.00015 | -0.00072 | -0.00086 | 1.0035/1.0030/1.0030 | -0.0008 / -0.0006 |
| 17 | normrel cap exp -3 (`VCC_NORMREL_EXP=-3`) | -0.0086 | -0.00101 | +0.0031 / +0.0042 | -0.0002 / -0.0013 / -0.0012 | **-0.00092** [-0.0013, -0.0002] 0/3 | -0.0004 / -0.0015 / -0.0014 | **-0.00111** [-0.0015, -0.0004] 0/3 | – | -0.00084 | -0.00162 | 1.0014/1.0001/0.9991 | +0.0005 / +0.0004 |
| 18 | H1 g ×0.5 (`VCC_H1_UNION_G=0.14`) | -0.0385 | +0.00021 | +0.0112 / +0.0160 | -0.0027 / -0.0022 / -0.0008 | **-0.00187** [-0.0027, -0.0008] 0/3 | -0.0035 / -0.0033 / -0.0014 | **-0.00273** [-0.0035, -0.0014] 0/3 | -0.00172 | -0.00203 | -0.00229 | 0.9956/0.9949/0.9950 | -0.0101 / -0.0133 |
| 19 | CD4 f ×0.75 (`VCC_CD4_TARGET_F=0.255`) | -0.0289 | -0.00085 | +0.0100 / +0.0137 | -0.0029 / -0.0021 / -0.0012 | **-0.00207** [-0.0029, -0.0012] 0/3 | -0.0036 / -0.0028 / -0.0018 | **-0.00272** [-0.0036, -0.0018] 0/3 | -0.00268 | -0.00192 | -0.00199 | 1.0051/1.0075/1.0071 | -0.0016 / -0.0020 |
| 20 | partial g 1 (`VCC_CD4_PARTIAL_G=1`) | -0.0056 | -0.00067 | +0.0105 / +0.0112 | -0.0016 / -0.0033 / -0.0014 | **-0.00209** [-0.0033, -0.0014] 0/3 | -0.0017 / -0.0034 / -0.0015 | **-0.00222** [-0.0034, -0.0015] 0/3 | -0.00239 | -0.00055 | -0.00168 | 1.0413/1.0382/1.0379 | -0.0024 / -0.0026 |
| 21 | CD4 deflation k 0 (`VCC_CD4_DEFLATE_K=unset`) | +0.0093 | -0.01228 | -0.0019 / -0.0031 | -0.0032 / -0.0030 / -0.0063 | **-0.00419** [-0.0063, -0.0030] 0/3 | -0.0028 / -0.0029 / -0.0062 | **-0.00398** [-0.0062, -0.0028] 0/3 | -0.00360 | -0.00401 | -0.00436 | 1.0047/1.0041/1.0043 | -0.0011 / -0.0006 |
| 22 | CD4 f ×0.5 (`VCC_CD4_TARGET_F=0.17`) | -0.0676 | -0.00349 | +0.0223 / +0.0309 | -0.0070 / -0.0048 / -0.0038 | **-0.00520** [-0.0070, -0.0038] 0/3 | -0.0086 / -0.0064 / -0.0051 | **-0.00671** [-0.0086, -0.0051] 0/3 | -0.00655 | -0.00503 | -0.00542 | 1.0062/1.0104/1.0096 | -0.0041 / -0.0051 |
| 23 | H1 g ×0 (`VCC_H1_UNION_G=unset`) | -0.0843 | -0.00191 | +0.0257 / +0.0364 | -0.0071 / -0.0057 / -0.0028 | **-0.00521** [-0.0071, -0.0028] 0/3 | -0.0090 / -0.0080 / -0.0042 | **-0.00708** [-0.0090, -0.0042] 0/3 | -0.00461 | -0.00531 | -0.00562 | 0.9968/0.9957/0.9958 | -0.0240 / -0.0305 |
| 24 | f1: key on the donor base only (no knob) | -0.1209 | -0.00085 | +0.0320 / +0.0471 | -0.0051 / -0.0080 / -0.0046 | **-0.00591** [-0.0080, -0.0046] 0/3 | -0.0077 / -0.0113 / -0.0068 | **-0.00857** [-0.0113, -0.0068] 0/3 | – | -0.00696 | -0.00630 | 0.9743/0.9784/0.9772 | -0.0140 / -0.0192 |
| 25 | H1 g ×4 (`VCC_H1_UNION_G=1.12`) | +0.0914 | -0.00864 | +0.0211 / +0.0089 | -0.0037 / -0.0054 / -0.0116 | **-0.00689** [-0.0116, -0.0037] 0/3 | -0.0012 / -0.0025 / -0.0105 | **-0.00475** [-0.0105, -0.0012] 0/3 | -0.00638 | +0.00180 | -0.01134 | 1.1039/1.1095/1.1109 | +0.0289 / +0.0429 |


#### 3.3 Secondary: each arm with its OWN LODO-V1 key recomputed on its composition (vs PP + V1, and + sharpening)

| arm | vs PP+V1: k=1 mean [range] +n | k=1.407 mean | vs PP+V1+sh: k=1 mean [range] +n | k=1.407 mean |
|---|---|---:|---|---:|
| H1 g ×2 (`VCC_H1_UNION_G=0.56`) | +0.00191 [-0.0002, +0.0034] 2/3 | +0.00318 | +0.00211 [+0.0002, +0.0034] 3/3 | +0.00337 |
| H1 g ×2.5 (`VCC_H1_UNION_G=0.7`) | +0.00178 [-0.0010, +0.0038] 2/3 | +0.00348 | +0.00200 [-0.0007, +0.0037] 2/3 | +0.00365 |
| CD4 f ×1.5 (`VCC_CD4_TARGET_F=0.51`) | +0.00184 [+0.0000, +0.0037] 3/3 | +0.00272 | +0.00168 [+0.0001, +0.0030] 3/3 | +0.00257 |
| H1 g ×3 (`VCC_H1_UNION_G=0.84`) | +0.00092 [-0.0023, +0.0033] 2/3 | +0.00291 | +0.00164 [-0.0015, +0.0035] 2/3 | +0.00353 |
| H1 g ×1.5 (`VCC_H1_UNION_G=0.42`) | +0.00145 [+0.0004, +0.0023] 3/3 | +0.00216 | +0.00153 [+0.0005, +0.0023] 3/3 | +0.00226 |
| CD4 f ×1.75 (`VCC_CD4_TARGET_F=0.595`) | +0.00172 [-0.0005, +0.0042] 2/3 | +0.00285 | +0.00145 [-0.0005, +0.0034] 2/3 | +0.00254 |
| CD4 f ×1.25 (`VCC_CD4_TARGET_F=0.425`) | +0.00129 [+0.0004, +0.0023] 3/3 | +0.00179 | +0.00121 [+0.0005, +0.0020] 3/3 | +0.00172 |
| CD4 f ×2 (`VCC_CD4_TARGET_F=0.68`) | +0.00128 [-0.0016, +0.0045] 2/3 | +0.00247 | +0.00081 [-0.0020, +0.0034] 2/3 | +0.00192 |
| f2: normrel PLAIN −1.5 (`VCC_NORMREL_MODE=plain`) | +0.00075 [+0.0006, +0.0011] 3/3 | +0.00075 | -0.00021 [-0.0004, +0.0001] 1/3 | -0.00020 |
| normrel cap exp -0.75 (`VCC_NORMREL_EXP=-0.75`) | +0.00012 [+0.0001, +0.0002] 3/3 | +0.00006 | -0.00021 [-0.0004, +0.0001] 1/3 | -0.00026 |
| partial g 0.25 (`VCC_CD4_PARTIAL_G=0.25`) | +0.00011 [-0.0006, +0.0009] 2/3 | +0.00005 | -0.00043 [-0.0012, +0.0003] 1/3 | -0.00049 |
| CD4 deflation k 3 (`VCC_CD4_DEFLATE_K=3`) | -0.00065 [-0.0012, +0.0002] 1/3 | -0.00049 | -0.00055 [-0.0009, -0.0000] 0/3 | -0.00040 |
| partial g 0.75 (`VCC_CD4_PARTIAL_G=0.75`) | -0.00065 [-0.0014, -0.0002] 0/3 | -0.00067 | -0.00070 [-0.0015, -0.0002] 0/3 | -0.00074 |
| normrel cap exp 0 (`VCC_NORMREL_EXP=unset`) | -0.00044 [-0.0008, -0.0001] 0/3 | -0.00059 | -0.00082 [-0.0016, -0.0004] 0/3 | -0.00098 |
| CD4 deflation k 8 (`VCC_CD4_DEFLATE_K=8`) | -0.00073 [-0.0019, +0.0007] 1/3 | -0.00086 | -0.00092 [-0.0023, +0.0006] 1/3 | -0.00104 |
| normrel cap exp -3 (`VCC_NORMREL_EXP=-3`) | -0.00077 [-0.0012, -0.0000] 0/3 | -0.00090 | -0.00100 [-0.0015, -0.0000] 0/3 | -0.00117 |
| partial g 0 (`VCC_CD4_PARTIAL_G=unset`) | -0.00068 [-0.0020, +0.0005] 1/3 | -0.00083 | -0.00108 [-0.0025, +0.0000] 1/3 | -0.00124 |
| partial g 1 (`VCC_CD4_PARTIAL_G=1`) | -0.00190 [-0.0032, -0.0010] 0/3 | -0.00196 | -0.00193 [-0.0034, -0.0010] 0/3 | -0.00206 |
| CD4 f ×0.75 (`VCC_CD4_TARGET_F=0.255`) | -0.00220 [-0.0033, -0.0015] 0/3 | -0.00292 | -0.00222 [-0.0030, -0.0018] 0/3 | -0.00297 |
| H1 g ×0.5 (`VCC_H1_UNION_G=0.14`) | -0.00197 [-0.0031, -0.0008] 0/3 | -0.00288 | -0.00232 [-0.0034, -0.0011] 0/3 | -0.00328 |
| H1 g ×4 (`VCC_H1_UNION_G=1.12`) | -0.00626 [-0.0107, -0.0034] 0/3 | -0.00392 | -0.00280 [-0.0068, -0.0007] 0/3 | -0.00069 |
| CD4 deflation k 0 (`VCC_CD4_DEFLATE_K=unset`) | -0.00452 [-0.0068, -0.0034] 0/3 | -0.00432 | -0.00413 [-0.0065, -0.0027] 0/3 | -0.00395 |
| CD4 f ×0.5 (`VCC_CD4_TARGET_F=0.17`) | -0.00543 [-0.0074, -0.0044] 0/3 | -0.00711 | -0.00589 [-0.0075, -0.0049] 0/3 | -0.00770 |
| H1 g ×0 (`VCC_H1_UNION_G=unset`) | -0.00553 [-0.0079, -0.0031] 0/3 | -0.00753 | -0.00619 [-0.0086, -0.0037] 0/3 | -0.00833 |


## 4. The two best one-variable bulk changes, priced vs 417a's design (the deliverable)

| | **H1-union g ×2** | **CD4 fraction ×1.5** |
|---|---|---|
| builder knob today | `VCC_H1_UNION_G` 0.28 → **0.56** | `VCC_CD4_TARGET_F` 0.34 → **0.51** (solved fraction 0.346 → 0.519 on A/B/C) |
| Δoverall vs 417a, k = 1, per fold (k562/hct116/hek293t) | +0.0031 / +0.0025 / +0.0004 → **+0.0020**, 3/3 | +0.0028 / +0.0016 / +0.0001 → **+0.0015**, 3/3 |
| Δoverall vs 417a, k = 1.407 | +0.0044 / +0.0040 / +0.0012 → **+0.0032**, 3/3 | +0.0036 / +0.0024 / +0.0007 → **+0.0022**, 3/3 |
| history-calibrated (retro class constants) | +0.0020 (H1 step: k_mse 0.92, k_pds 0.70) | **+0.0023** (CD4-dose step: k_mse 1.41, k_pds 0.88) |
| amplitude-matched (direction only) / PLOLO-energy route | +0.0022 / +0.0021 | +0.0016 / +0.0013 |
| Δrel cos / Δpds_bd / Δraw mse (k 1…1.407) | +4.9 % / **−0.0021** / −0.016…−0.023 | +3.0 % / −0.0008 / −0.010…−0.014 |
| bootstrap P(Δrel>0), P(Δpds_bd>0) per fold | 1.00/0.58, 1.00/0.06, 1.00/0.02 | 1.00/0.81, 1.00/0.35, 1.00/0.12 |
| with its own LODO-V1 key, vs PP+V1+sh | +0.0021 (3/3) … +0.0034 | +0.0017 (3/3) … +0.0026 |
| vs 410b (no sharpening) | +0.0018 … +0.0031 (3/3) | +0.0016 (2/3) … +0.0023 (3/3) |
| dose neighbours vs 417a (k = 1) | ×1.5 +0.0015 (3/3), ×2.5 +0.0018 (2/3), ×3 +0.0014 (2/3), ×4 −0.0034 | ×1.25 +0.0011 (3/3), ×1.75 +0.0011 (2/3), ×2 +0.0006 (2/3) |
| what the knob does on the real panel | moves the direction of **26 %** of rows (the 77 H1-carrying targets; median cos 0.978 on them) and adds render energy +1.0/+1.1/+1.1 % (grows with dose: +3 % at ×2.5, +10 % at ×4) | moves the direction of **92 %** of rows (median cos 0.990), render energy 1.001/0.995/0.995 (neutral) |
| ballast realisability (A/B/C, sharpened) | K = 16 on 19.8 % of blocks (417a 20.1 %), cos(realised, design) q1 0.9918 | 19.3 %, q1 0.9915 |
| October (D/E/F) | H1 coverage on a new panel is unknown (64 + 20 of today's 300; this panel was drawn from an H1 pilot screen), so the dose may act on few rows | CD4 is genome-wide (291/300 today); the dose carries over |

**Which one first?** H1 ×2 has the higher k-range, but it carries a real pds cost on two of three folds and an energy
excursion that puts context A further over its cliff. It also depends on today's H1 coverage. CD4 ×1.5 is
energy-neutral, and its pds CIs span zero. Its family has transferred to the board at or above PLOLO twice: cd4f034's
cosine k 1.405, and the retro k_mse 1.41. On the history calibration the two are tied (+0.0020 vs +0.0023). **For one
slot I would take CD4 ×1.5 as the lower-variance arm and keep H1 ×2 as the follow-up.** Their stack was not measured, and
both add residual through one shared post-add preserve, so expect partial redundancy.

**Both knobs act on the carriers too.** The composition stage builds ONE shared effect matrix. The carriers render its
unkeyed rows (the DE design) and the ballast renders its keyed (and, in 417a, sharpened) rows. Setting either env var
therefore also moves the carriers' DE table. That welded version is priced nowhere here; it is 416's territory. For
reference, the welded CD4 dose step 0.10 → 0.34 moved fid/reach/jac +0.0061 overall on the board. A **bulk-only** arm
needs a builder change: compose a second matrix for `bulk_effects` with its own `VCC_BULK_CD4_TARGET_F` /
`VCC_BULK_H1_UNION_G`, then apply the pinned key and the 417 sharpening to it, while the carriers keep the shipped
matrix. The ballast table above shows the carriers' unchanged design does not constrain these bulk designs beyond
410b/417a.

### 4.1 Which knobs are scale-like and which change direction (measured on the real A/B/C panel, `abc418.json`)

A pure scale knob would move only the Frobenius norm and render energy, with every row's cosine at 1. **None of the five
is pure scale.**

| knob | rows whose direction moves | median cos on those rows | effect-space Frobenius ratio | render energy A/B/C (keyed) | class |
|---|---:|---:|---:|---|---|
| `VCC_CD4_TARGET_F` ×1.5 (×0.5) | 92 % | 0.990 (0.987) | 0.992 (1.011) | 1.001/0.995/0.995 (1.006/1.010/1.010) | **direction**; energy ≤ 1 % |
| `VCC_H1_UNION_G` ×2 (×0) | 26 % | 0.978 (0.969) | 1.005 (1.000) | 1.009/1.011/1.011 (0.997/0.996/0.996) | **direction on H1 rows, plus an amplitude side effect** that grows with dose |
| `VCC_CD4_DEFLATE_K` 0 / 3 / 8 | 92 % | 0.995 / 0.999 / 0.999 | 1.000 / 1.001 / 1.001 | +0.4 % / +0.3 % / +0.1 % | **direction** |
| `VCC_CD4_PARTIAL_G` 0 / 1.0 | 5 % (15 partial targets) | 0.922 / 0.956 | 0.999 / 1.014 | −0.8 % / **+4.0 %** | few rows with large norms: **amplitude-dominated** |
| `VCC_NORMREL_EXP` (cap) off / −3 | 100 % | 0.997 / 0.997 | 1.000 / 0.996 | −0.7 % / ±0.1 % | **direction** (small) |
| `VCC_NORMREL_MODE=plain` (f2) | 100 % | 0.998 | 0.992 | **−1.9 %** | **scale-like in effect**: its whole price is the energy drop (amp-matched ≈ 0) |

## 5. Retrodiction: ten historical board composition steps (cancer folds, A/B/C-panel energies)

Each step is rebuilt on PLOLO without a key, and its board raw Δpds / Δmse is predicted from the base entry's own board
raw. The board raws come from `data/board-entries/*.json`. rc_callset and cd4f034 come from the `data/leaderboard/`
snapshots, because the 210 ledger table's raw pds is on another anchor.


| step | board Δpds | PLOLO Δpds | board Δmse | pred Δmse k=1 | pred k=1.407 | k_fit (mse) | board Δ(pds+mse) overall | pred k=1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CD4 target f 0.10 -> 0.34 | -0.00371 | -0.00422 | -0.03726 | -0.02638 | -0.03726 | 1.41 | +0.00516 | +0.00306 |
| H1 union g 0 -> 0.28 | +0.00344 | +0.00491 | -0.01928 | -0.02100 | -0.02938 | 0.92 | +0.00465 | +0.00549 |
| CD4 deflation k 0 -> 5 | +0.01492 | +0.01242 | +0.00454 | +0.00180 | +0.00218 | 3.96 | +0.00471 | +0.00427 |
| normrel off -> plain -1.5 | +0.00551 | +0.00208 | -0.00509 | -0.00085 | -0.00091 | 27.66 | +0.00293 | +0.00092 |
| normrel plain -1.5 -> -3 | +0.00091 | -0.00242 | +0.00360 | +0.00435 | +0.00597 | 0.81 | -0.00029 | -0.00166 |
| normrel plain -> cap (-1.5) | -0.00091 | -0.00058 | -0.00084 | -0.00073 | -0.00081 | 1.53 | -0.00019 | -0.00008 |
| normrel cap -1.5 -> -3 | +0.00245 | -0.00090 | -0.00029 | +0.00537 | +0.00655 | -0.94 | +0.00096 | -0.00127 |
| partial g 0 -> 0.34 (no normrel) | -0.00000 | +0.00178 | -0.00449 | -0.00298 | -0.00399 | 1.61 | +0.00079 | +0.00118 |
| partial g 0 -> 0.34 (cap) | -0.00009 | +0.00169 | -0.00404 | -0.00244 | -0.00347 | 1.63 | +0.00067 | +0.00105 |
| partial g 0.34 -> 0.5 | -0.00042 | -0.00003 | -0.00194 | +0.00053 | +0.00023 | 4.38 | +0.00019 | -0.00010 |

sign agreement: overall 8/10, pds 6/10, mse 8/10; corr overall +0.871, pds +0.903, mse +0.957


- **Dose families (CD4 f, H1 g, deflation k, partial g, 6 steps).** The overall board sign is right on 5/6 (6/6 with
  PLOLO energies). mse is right on 5/6. Both misses are the partial 0.34 → 0.5 step: the board moved mse −0.0019, and
  the A/B/C route predicts +0.0005. Magnitudes are within about ±0.002 overall.
- **Transfer.** Pooled class constants: CD4 dose k_pds 0.88, k_mse 1.41; H1 0.70 / 0.92; deflation 1.20 / 3.96 (PLOLO
  under-states deflation's mse cost); partial k_mse 1.97.
- **Normrel (4 steps).** The cancer folds get the pds sign wrong on plain −1.5 → −3 and cap −1.5 → −3 (2 of 4). They
  read off → plain at 0.4× the board's pds gain, and they predict an mse cost for −3 that the board did not show.
- **kolf2 on normrel.** kolf2 is the only fold whose normrel acts on three donors, as it does on A/B/C. It gets all four
  normrel pds signs right (below). It is off on the dose families: it over-reads the H1 step about 5× and gets the CD4-dose pds sign wrong.

kolf2-only retrodiction (A/B/C energies):


| step | board Δpds | PLOLO Δpds | board Δmse | pred Δmse k=1 | pred k=1.407 | k_fit (mse) | board Δ(pds+mse) overall | pred k=1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CD4 target f 0.10 -> 0.34 | -0.00371 | +0.00446 | -0.03726 | -0.04040 | -0.05698 | 0.92 | +0.00516 | +0.00872 |
| H1 union g 0 -> 0.28 | +0.00344 | +0.02216 | -0.01928 | -0.08885 | -0.12485 | 0.21 | +0.00465 | +0.02374 |
| CD4 deflation k 0 -> 5 | +0.01492 | +0.00431 | +0.00454 | +0.00836 | +0.01140 | 0.49 | +0.00471 | +0.00013 |
| normrel off -> plain -1.5 | +0.00551 | +0.00829 | -0.00509 | -0.00358 | -0.00476 | 1.52 | +0.00293 | +0.00369 |
| normrel plain -1.5 -> -3 | +0.00091 | +0.00277 | +0.00360 | +0.00297 | +0.00402 | 1.24 | -0.00029 | +0.00050 |
| normrel plain -> cap (-1.5) | -0.00091 | -0.00392 | -0.00084 | +0.00035 | +0.00072 | -0.34 | -0.00019 | -0.00151 |
| normrel cap -1.5 -> -3 | +0.00245 | +0.00263 | -0.00029 | +0.00329 | +0.00363 | -3.28 | +0.00096 | +0.00039 |
| partial g 0 -> 0.34 (no normrel) | -0.00000 | +0.00020 | -0.00449 | -0.00288 | -0.00385 | 1.68 | +0.00079 | +0.00058 |
| partial g 0 -> 0.34 (cap) | -0.00009 | +0.00029 | -0.00404 | -0.00232 | -0.00331 | 1.71 | +0.00067 | +0.00051 |
| partial g 0.34 -> 0.5 | -0.00042 | -0.00015 | -0.00194 | +0.00055 | +0.00026 | 4.47 | +0.00019 | -0.00015 |

sign agreement: overall 8/10, pds 7/10, mse 7/10; corr overall +0.644, pds +0.277, mse +0.733


## 6. The normrel exponent: the cancer folds cannot price it; the evidence says −3 might be worth +0.0005–0.001

The cancer folds' 2 source donors turn per-target donor weighting into a two-way choice. The A/B/C panel weights
three. Evidence on `VCC_NORMREL_EXP=-3` (cap), bulk only:
- The board, pre-key, measured cap −1.5 → −3 at pds +0.00245 raw and mse −0.0003 raw, **+0.0010 overall on the bulk
  members**. The six-member entry lost −0.0003 only because nmae paid (319d3, 2026-09-09).
- kolf2 (3 donors) reads **+0.0005 vs 410b and +0.0004 vs 417a** (k = 1).
- The cancer folds read −0.0009 / −0.0011, but they got this very step's pds sign wrong in the retro.
- On the real panel it does not change render energy (1.001/1.000/0.999), unlike the 2-donor folds' +2 %.

So normrel −3 is a plausible small bulk gain (+0.0005…+0.001) that PLOLO's cancer folds are blind to. It is not
recommended on this evidence. It needs a 3-donor surface or a board read. Exponent −0.75 and "off" are flat to negative
on every surface.

## 7. Own ideas (f)

- **f1: the key on the donor base only, with the CD4/H1/partial residual added unkeyed.** The rationale was that 415
  showed the residual's projection is informative on its own, and that the 410 key measures donor agreement only.
  - Result: **−0.0059 to −0.0086 vs 410b and −0.0077 to −0.0104 vs 417a, 0/3 folds**; rel cos −11 %.
  - The residual rows' reliability co-varies with the donor base's shared energy, so the key must scale the whole row.
    That is consistent with 415, where adding the residual term to the key's numerator improved it.
  - Closed.
- **f2: normrel PLAIN mode at −1.5.** The cap was chosen on the welded board by its DE members: 319 vs 312 was pds
  −0.0021 raw, nmae and reach +. So under ballast the plain mode was a natural bulk candidate.
  - Result: +0.0008 vs 410b (3/3), but only through the −1.9 % render-energy drop (amp-matched +0.0001). Vs 417a it is
    −0.0001, because sharpening already brings context A down toward its cliff.
  - Take the per-context re-levelling knob 414 proposed instead. Closed as a composition lever.

## 8. Limits

- **Direction transfer is unmeasured for these exact arms.** The k range 1…1.407 and the history constants bracket it.
  The retro constants come from one board step per family (two or three for partial), with bases that had no key and no
  sharpening. FINDINGS H.5: treat a class constant measured once as a lower bound on magnitude.
- **The board-emulated pds** (β at λ 0.5) is 414's one-parameter emulation fitted to the key's board pds cost. For the
  winners, literal and emulated pds agree within 0.0001.
- **A/B/C render energies** come from `plolo.deliver` on the composed design, not from the production emitter. 414
  validated the PLOLO-ratio route to ≤ 0.3 % on 410c. The A/B/C route uses the same renderer on the real panel. PLOLO's
  own energy route is kept as a column, and it ranks the winners identically.
- **The LODO-V1 secondary** uses A/B/C energy ratios of the shipped key, since only that key exists on the panel. Its
  direction and pds are its own.
- **kolf2 cannot vote** by mandate. It is also structurally biased on H1 arms (iPSC truth vs the hESC source), and it is
  the only honest surface for normrel.
- **The welded (both-channel) versions of the knobs** and the stack of the two winners were not priced.

## 9. Files and reproduction

```
cd experiments/418-bulk-composition           # everything writes only into ./work
bash chain418.sh k562_gwps orion_hct116 orion_hek293t kolf2   # comp418.py <fold>: KAs + 21 arms x {K3, K3+sh, V1, V1+sh} + 11 retro rows
uv run python abc418.py                       # the real A/B/C panel: KA-PP, every arm's render A_c/B_c and row shape
bash chain418b.sh                             # abc418 + the dose continuation (comp418.py <fold> --extra: H1 x2.5/3/4, CD4 x1.75)
uv run python ballast418.py                   # ballast realisability of the top arms on A/B/C
uv run python analyse418.py                   # both mse routes, retro + class constants, bootstrap -> work/analyse418.{json,log}
uv run python table418.py > work/table418.md  # the §3/§5 tables
```

- Libraries: `lib418.py` holds the parameterised composition, `pre_normrel`, the LODO-V1 key and the read-back.
- Fold caches: APFS clones of 414's (`work/plolo/`).
- Outputs: `work/comp418[x]-<fold>.json`, `work/arrays418[x]-<fold>.npz` (per-target arrays for the bootstrap),
  `work/abc418.json`, `work/ballast418.json`, `work/analyse418.{json,log}`, `work/table418.md`.
- Peak RSS per process (`/usr/bin/time -l`), all under the 8 GB mandate:
  - comp418: 6.10 / 6.69 / 7.03 / 7.09 GB (hct116 / k562 / hek293t / kolf2).
  - comp418x: 6.07–6.54 GB.
  - abc418: 6.17 GB. ballast418: 5.34 GB.
- Every heavy step ran behind the memory gate (≥ 40 % free and ≥ 20 GiB free+inactive).
- **One write outside this directory, since removed.** Importing `ballast414` created
  `experiments/414-bulk-reprice/__pycache__/ballast414.cpython-312.pyc`, a bytecode cache that was gitignored and new.
  I deleted it. The chain scripts now set `PYTHONDONTWRITEBYTECODE=1`; run `ballast418.py` with it set too.
  Nothing else outside `418/` changed: `find -newer MANDATE.md` over the imported experiments lists only other agents'
  410 files.
