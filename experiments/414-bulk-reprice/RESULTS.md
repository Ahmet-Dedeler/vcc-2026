# 414: bulk-only arms re-priced under ballast. Sharpening (169) is the lever it reopens; 410d's built artifact is an amplitude blow-up (NOGO)

Status: **COMPLETE** (2026-09-22, 20:35Z). PLOLO on four folds, three board artifacts streamed read-only, and two new
instruments, each anchored to a board number. No vcc, no cell-eval2, no builder edit, no build. Nothing written outside
`experiments/414-bulk-reprice/`. **Read `URGENT-410d.md` first if 410d has not been uploaded yet.**

## 0. The answer in nine lines

1. **Every known answer reproduces exactly.** 410's check410 key gains are +0.299/+0.267/+0.287 (clip 3) and
   +0.386/+0.341/+0.350 (uncapped), with |d| = 0. 228's arms228 rows (BASE, gain_reg t 0.1/0.25/0.5/1, cos_agg and pds)
   match to 1e-12 on 4/4 folds. 169's sharpening rows (q 1.00/1.15) have |d| = 0. 187's model gives X23 0.017930 and
   cliffs 0.9823/1.0129/1.0386. The 395 artifact's A_c/B_c equal 393's decomposition, |d| = 0: the 395 bulk *is* 319p2's.
2. **Recommended next bulk arm, one variable on 410c: 169's magnitude sharpening on the bulk only, q = 1.25.** The
   transform is `sign(h)|h|^q`, with each row renormalised to its own L2 norm, on each context's CPM gate, applied to
   `bulk_effects` only (the carriers stay untouched). Predicted vs 410c:
   - pds +0.0063 to +0.0191 raw.
   - raw mse −0.020 to −0.050.
   - **overall +0.0055 (kolf2-like transfer) to +0.0151 (cancer folds)**, with render-space energy measured on the real
     410c bulk.
   - Positive on 4/4 folds, with paired-bootstrap P(Δ>0) = 1.00 on cos and ≥ 0.99 on pds.
   - Both controls go negative on 4/4 folds: flattening (q 0.87) −0.0087 overall, permuted multipliers −0.0016.
3. **410d as preregistered (228's gain on the final bulk row, Frobenius-preserved in effect space) is an amplitude arm.
   Measured on its built artifact, it is worse than predicted.**
   - Render-space energy vs 410c is **×1.197 / ×1.182 / ×2.127** (A/B/C). A and B match PLOLO kolf2's prediction
     (×1.19/×1.19/×1.20; kolf2's multiplier is the shipped artifact bit-for-bit).
   - C doubles because one target explodes. **SYVN1-C's row energy goes from 229 to 4,858.** m ≈ 4–5 on the UPR genes,
     and HSP90B1 and HSPA5 take 50.5 % and 12.7 % of that block's counts (2.4 % and 0.3 % in 410c).
   - Every context ends up over its cliff (ρ 1.25/1.15/1.98). The amplitude-only raw-mse cost is **+0.50 vs 410b**
     (context C +1.12), about −0.03 to −0.04 overall. The PI had already written NOGO-410d at 20:06Z from this
     worker's first note.
   - **410d2's render-gated matching still gives HSP90B1 a ×2.0 boost, against 410d's ×2.2 that blew up.** Its
     emulation cannot settle C (known-answer row fails there, §4). Clipping m to [1/2, 2] is the safe variant.
   - The fair gain arm is amplitude-matched and bounded: +0.0025 to +0.0046 overall. The 197 slope does not see
     amplitude and calls 410d +0.004.
4. **410c: measured, and predicted about flat.** Its artifact's render-space energy is +0.99/+0.84/+1.49 % over 410b.
   PLOLO predicted +1.3/+1.1/+1.2 %. Context A sits further over its cliff (rho_A 1.043 against 1.032).
   - Predicted Δraw mse −0.006 vs 410b; the amplitude-only cost is +0.0033 of that.
   - Board-emulated pds is −0.0026: 42 % of targets sit at the 0.13× floor, under a pattern that did not scale.
   - Net overall ≈ +0.000 (range −0.001 to +0.001).
5. **(c) amplitude: the optimum does not move above A's cliff. It sits on the cliff for any X we have.** Past the cliff,
   raw mse costs +1.78 per unit scale; below it, raw mse falls 0.31 per unit. The optimum moves only if X > s·A_c, and
   X/A_c is about 0.15.
   - The key moved the BULK past A's cliff instead: +2.3 % render-space energy. 410c's bulk is 2.2 % (in scale) above it.
   - Per-context bulk scales to the cliffs on 410c's bulk (A 0.978, B 1.0136, C 1.0386; α-equivalents
     0.9604/0.9954/1.0200): raw mse −0.0185, **overall +0.0032**.
   - Context A alone gives +0.0023. One global scale of 0.978 gives +0.0015.
6. **A higher floor removes the key's board pds cost for free.** lo 0.5 (uncapped) vs lo 0.3: cosine +0.2 %,
   board-emulated pds +0.0028, overall +0.0012. On kolf2 the effect is flat.
7. **The gene gain's DIRECTION is real (+3.9 % cos on cancer, +4.5 % on kolf2), stacks with the key, and survives the
   read-back.** Its pds is not robust (+0.0046 cancer, −0.0018 kolf2). The pre-residual (228) and final-row (410d)
   placements carry the same energy excursion (+10–12 %).
8. **Two new instruments, both anchored to the board:**
   - A cliff-aware mse route: 187's model on A_c/B_c measured from the shipped artifacts, with X solved per artifact
     from its board raw. Calibrated on 410b, it retrodicts the welded 410's board raw to +0.0002.
   - A board-emulated pds: the artifact's non-scaling shared pattern β added to PLOLO rows at λ = 0.5. It reproduces the
     key's pds cost (−0.0026 against the board's −0.0036). PLOLO's read-back alone reads −0.0006.
9. **Peak RSS exceeded the 8 GB mandate on two multi-fold runs** (extra414 8.58 GB, ballast414 8.73 GB). Caches
   accumulate across folds in one process. Memory pressure stayed ≥ 59 % free and no build was affected. Run one fold per
   process next time.

---

## 1. Known-answer rows (all before any new number; each script aborts on failure)

| # | row | expected | measured |
|---|---|---|---|
| KA-1 | check410 key gains, effect level (arms414.py `ka_check410`) | 0.299/0.267/0.287/0.500 (clip 3), 0.386/0.341/0.350/0.936 (∞) | identical, \|d\| = 0 |
| KA-2 | 228 `arms228-<fold>.json`: BASE, gain_reg@0.1/0.25/0.5/1, cos_agg and pds | 1e-12 | \|d\| = 0, 4/4 folds |
| KA-3 | 169 `plolo169.json` a23 stack, sharpening q 1.00 / 1.15, pds and sse per frame | 1e-9 | \|d\| = 0, 4/4 folds |
| KA-4 | 187 `joint187.Model(F_FIT)` | X23 0.017930; cliffs 0.9823/1.0129/1.0386 | identical |
| KA-5 | 395 artifact A_c, B_c vs 393's decomposition | (393 streamed 319p2, now swept) | \|d\| = 0: 395 carries 319p2's bulk exactly |
| KA-6 | block stream vs an independent anndata dense read, 3 blocks × 5 artifacts (395/410/410b/410c/410d); Σw = 300/18532 | exact | \|d\| = 0 |
| KA-7 | board414 / extra414 / ballast414: λ = 0 and design rows == arms414 | 1e-12 | \|d\| = 0 |
| CTRL | key permuted over targets (Kinf_shuf, K3_shuf) | negative | rel cos −0.31 / −0.23, 0/4 folds positive |

## 2. Instruments

**2.1 Composition and read-back.** This is 228's `price` (whitened base at α 0.982, CD4 @ 0.34 plus the H1 union @ 0.28,
panel coordinates zeroed). The 410 key is applied to the final row by the builder's `apply_pert_scale` arithmetic
(w·c, non-panel Frobenius preserved, inner products over the operator axis as in `build_weights410.py`). Delivery is
through `plolo.deliver` (a23 minnorm) in the challenge frames A/B/C. The metrics are pds, cos_agg = sqrt(1−sse), and the
render-space row energy A and centred energy B over non-panel genes. In this composition the key reads +0.244/+0.274/+0.272
(clip 3) and +0.247/+0.313/+0.309 (∞). The effect-level check410 numbers are 0.299/0.267/0.287.

**2.2 Cliff-aware mse (analyse414.py).** 187's per-context identity `MSE_c − 1 = [A_c + Kc_c − min(claim_c, B_c) − 2X]/ΣD`
uses claim_c/Kc_c/ΣD from 187's retrodicted model, and A_c/B_c streamed off the shipped artifacts (`decomp414.py`,
cell_eval2's own `_row_weights`/`_across_pert_budget`).

| artifact | A_c (A/B/C) | ρ_c = B_c/claim_c | cliff s_c* = √(claim/B) | board raw | X (solved) |
|---|---|---|---|---:|---:|
| 395 | 0.237515 / 0.225984 / 0.236052 | 1.0078 / 0.9416 / 0.8951 | 0.996 / 1.031 / 1.057 | 0.809898 | 0.025774 |
| 410 (welded) | 0.243571 / 0.231649 / 0.242008 | 1.0342 / 0.9658 / 0.9183 | 0.983 / 1.018 / 1.044 | 0.730213 | 0.035882 |
| 410b | 0.243085 / 0.230999 / 0.241362 | 1.0319 / 0.9629 / 0.9156 | 0.984 / 1.019 / 1.045 | 0.730994 | 0.035723 |
| 410c (artifact, `cliffcheck414.py`) | 0.245499 / 0.232935 / 0.244949 | 1.0426 / 0.9713 / 0.9296 | 0.979 / 1.015 / 1.037 | — | — |

- **Calibration.** The board's relative aggregate-cosine gain of 410b over 395 is X ratio / √(A ratio) = **+0.3705**
  (197's slope reads +0.3555). PLOLO's K3 reads +0.2633, so **k_key = 1.407**.
- **Retrodiction.** From 395's board raw plus PLOLO's K3 energy ratio and gain: at k_key, raw 0.73046 (board 410 0.73021,
  410b 0.73099). At k = 1, raw 0.75481, off by 0.024.
- **Render-space energy of the key.** PLOLO K3 reads 1.021/1.019/1.019. The artifacts read 1.026/1.025/1.025 (410) and
  1.023/1.022/1.023 (410b). For 410c, PLOLO predicted A_c(410c)/A_c(410b) = 1.013/1.011/1.012 and the artifact measures
  1.010/1.008/1.015.

**2.3 Board-emulated pds (shared414.py → board414.py / extra414.py).** Fitting d^410b_p = w_p s_p + β, d^395_p = s_p + β on
the artifacts (w = the shipped final scale) finds a pattern β common to all rows that the key did not scale. ‖β‖² is
0.047/0.051/0.052, about 10 % of a floor row's energy. Floor rows' cosine with β is 0.30–0.32 in 410b, against 0.06 in 395.
Floor rows' cosine with their own 395 row is 0.94–0.95, against 0.998 for the rest.

PLOLO's own render produces a far smaller common part (‖β_r‖² 0.002–0.003, cos with β ≈ +0.25), so β is emitter-side and
is not double-counted. Adding β, rescaled to PLOLO's row energy, to every delivered row gives K3 − base pds:
- λ 0: −0.0006
- λ 0.5: −0.0026
- λ 1: −0.0091
- board 410b − 395: −0.0036

λ = 0.5 is the calibrated read (one number fitted to one board point). On kolf2, β *raises* pds (its truth aligns with
β), so kolf2's board-emulated column is not interpretable.

**2.4 Ballast realisability (ballast414.py).** This emulates the builder's `ballast_count`: K = ⌈s_c(1/r_min − 1)/995k⌉
capped at 16, and genes below the carriers' totals are clipped. Every design here is realised with small loss:

| design | K capped | realised − design, rel cos | realised − design, pds |
|---|---:|---:|---:|
| Kinf | 8.4 % | −0.007 | −0.00001 |
| Kinf+sh@1.25 | 15.5 % | −0.017 | −0.00023 |
| Kinf+sh@1.5 | 21.6 % | −0.027 | −0.00058 |

Clip mass q90 ≤ 1.3e-5 of a block. The added capped blocks are the thing to watch: `ballast_jk_ratio` in the build meta.

## 3. Ranked table: bulk-only designs vs 410c (Kinf); three cancer folds unless marked

Columns:
- **rel**: relative aggregate cosine vs the 228 base.
- **Δpds literal**: PLOLO read-back; **board-emul.**: λ = 0.5.
- **Δraw mse 197**: −2.0745·0.107·k_key·Δrel. For k = 1, multiply the mse part by 1/1.407.
- **cliff**: §2.2's route.
- **Δoverall mandate**: 197 slope plus literal pds. **primary**: cliff route plus board-emulated pds (literal pds where
  emulation was not run).
- **levelled**: the same design with each context's bulk scale set to its cliff.
- **render A ratio**: render-space energy vs base (Kinf 1.032).

| arm (bulk only, on 410c) | rel cos vs base | Δrel vs 410c | Δpds literal | Δpds board-emul. | Δraw mse 197 (k_key) | Δraw mse cliff | Δoverall mandate | **Δoverall primary** | Δoverall if levelled | render A ratio | folds + (rel/pds) | kolf2 Δrel / Δpds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| `Kinf+sh@1.5` | +0.539 | +0.250 | +0.0346 | +0.0347 | -0.0781 | -0.0747 | +0.0264 | **+0.0259** | +0.0303 | 0.922 | 3/3, 3/3 | +0.063 / +0.0093 |
| `Kinf+sh@1.35` | +0.478 | +0.188 | +0.0258 | +0.0259 | -0.0588 | -0.0645 | +0.0198 | **+0.0209** | +0.0235 | 0.968 | 3/3, 3/3 | +0.051 / +0.0077 |
| `K0.5_inf+sh@1.25` | +0.433 | +0.143 | +0.0191 | +0.0221 | -0.0448 | -0.0512 | +0.0149 | **+0.0171** | +0.0196 | 0.994 | 3/3, 3/3 | +0.041 / +0.0065 |
| `Kinf+sh@1.25` | +0.429 | +0.140 | +0.0190 | +0.0191 | -0.0436 | -0.0501 | +0.0147 | **+0.0158** | +0.0182 | 0.993 | 3/3, 3/3 | +0.040 / +0.0063 |
| `K0.5_inf+sh@1.15` | +0.379 | +0.089 | +0.0117 | +0.0147 | -0.0279 | -0.0312 | +0.0092 | **+0.0109** | +0.0137 | 1.013 | 3/3, 3/3 | +0.026 / +0.0042 |
| `Kinf+sh@1.15` | +0.376 | +0.086 | +0.0117 | +0.0118 | -0.0269 | -0.0303 | +0.0090 | **+0.0096** | +0.0125 | 1.013 | 3/3, 3/3 | +0.026 / +0.0041 |
| `K3+sh@1.15` | +0.353 | +0.064 | +0.0114 | +0.0142 | -0.0199 | -0.0249 | +0.0077 | **+0.0096** | +0.0121 | 1.007 | 3/3, 3/3 | -0.144 / +0.0037 |
| `Kinf+shL@1.15` (lfc space) | +0.362 | +0.072 | +0.0109 | +0.0109 | -0.0226 | -0.0255 | +0.0080 | **+0.0085** | +0.0114 | 1.013 | 3/3, 3/3 | +0.022 / +0.0036 |
| `Kinf+sh@1.1` | +0.347 | +0.058 | +0.0079 | +0.0080 | -0.0181 | -0.0202 | +0.0061 | **+0.0065** | +0.0095 | 1.020 | 3/3, 3/3 | +0.018 / +0.0026 |
| `Kinf+sh@1.05` | +0.319 | +0.029 | +0.0040 | +0.0040 | -0.0091 | -0.0101 | +0.0031 | **+0.0033** | +0.0064 | 1.027 | 3/3, 3/3 | +0.009 / +0.0012 |
| `K0.7_inf` | +0.293 | +0.003 | +0.0001 | +0.0036 | -0.0010 | -0.0011 | +0.0002 | **+0.0015** | +0.0047 | 1.031 | 2/3, 3/3 | -0.003 / +0.0001 |
| `K0.5_inf` | +0.292 | +0.002 | +0.0001 | +0.0028 | -0.0008 | -0.0008 | +0.0002 | **+0.0012** | +0.0044 | 1.032 | 3/3, 3/3 | -0.000 / -0.0000 |
| `K0.4_inf` | +0.291 | +0.001 | +0.0000 | +0.0018 | -0.0004 | -0.0004 | +0.0001 | **+0.0007** | +0.0040 | 1.032 | 3/3, 3/3 | +0.000 / -0.0000 |
| `K0.5_3` | +0.264 | -0.025 | -0.0000 | +0.0039 | +0.0079 | +0.0055 | -0.0014 | **+0.0005** | +0.0032 | 1.019 | 1/3, 3/3 | -0.173 / -0.0002 |
| **`Kinf` (= 410c)** | +0.290 | 0 | 0 | 0 | 0 | 0 | 0 | **0** | +0.0032 | 1.032 | — | — |
| `K3` (= 410b) | +0.263 | -0.026 | -0.0001 | +0.0026 | +0.0082 | +0.0058 | -0.0015 | **-0.0001** | +0.0027 | 1.019 | 0/3, 3/3 | -0.171 / -0.0003 |
| `Kinf+gf@0.1` | +0.341 | +0.052 | +0.0046 | +0.0049 | -0.0161 | +0.0141 | +0.0045 | **-0.0006** | +0.0080 | 1.110 | 3/3, 3/3 | +0.050 / +0.0004 |
| **`Kinf+gf@0.25` (= 410d)** | +0.341 | +0.051 | +0.0044 | +0.0046 | -0.0161 | +0.0270 | +0.0044 | **-0.0030** | +0.0078 | 1.125 | 3/3, 2/3 | +0.045 / -0.0018 |
| `gi@0.1+Kinf` | +0.318 | +0.028 | +0.0042 | – | -0.0089 | +0.0178 | +0.0031 | **-0.0016** | +0.0064 | 1.101 | 3/3, – | +0.023 / +0.0010 |
| `gi@0.25+Kinf` | +0.320 | +0.030 | +0.0051 | +0.0052 | -0.0095 | +0.0316 | +0.0035 | **-0.0036** | +0.0069 | 1.122 | 3/3, 3/3 | +0.022 / +0.0012 |
| `Kinf+gf@0.5` | +0.318 | +0.029 | +0.0020 | – | -0.0090 | +0.0394 | +0.0023 | **-0.0062** | +0.0056 | 1.131 | 3/3, – | +0.018 / -0.0052 |
| `gi@0.5+Kinf` | +0.308 | +0.018 | +0.0038 | – | -0.0056 | +0.0513 | +0.0024 | **-0.0076** | +0.0056 | 1.140 | 3/3, – | +0.015 / -0.0000 |
| `Kinf+gf@1` | +0.275 | -0.014 | -0.0021 | – | +0.0044 | +0.0471 | -0.0015 | **-0.0090** | +0.0016 | 1.123 | 0/3, – | -0.013 / -0.0079 |
| `gi@1+Kinf` | +0.284 | -0.006 | +0.0005 | – | +0.0018 | +0.0645 | -0.0001 | **-0.0111** | +0.0030 | 1.146 | 0/3, – | +0.010 / -0.0003 |
| `Kinf+gf@0.25+sh@1.15` | +0.405 | +0.115 | +0.0142 | – | -0.0361 | -0.0256 | +0.0115 | **+0.0097** | +0.0150 | 1.073 | 3/3, – | +0.052 / -0.0015 |
| `Kinf@x0.97` (global scale) | +0.291 | +0.001 | +0.0000 | -0.0003 | -0.0004 | -0.0068 | +0.0001 | **+0.0011** | +0.0032 | 0.972 | 3/3, 0/3 | +0.000 / -0.0000 |
| `Kinf@x1.03` | +0.288 | -0.001 | -0.0000 | – | +0.0004 | +0.0237 | -0.0001 | **-0.0042** | +0.0032 | 1.093 | 0/3, – | -0.000 / +0.0000 |
| CONTROL `Kinf+shperm@1.15` | +0.273 | -0.016 | -0.0019 | -0.0018 | +0.0051 | +0.0057 | -0.0016 | **-0.0016** | +0.0016 | 1.031 | 0/3, 0/3 | -0.012 / -0.0022 |
| CONTROL `Kinf+sh@0.87` (flatten) | +0.211 | -0.079 | -0.0110 | -0.0108 | +0.0247 | +0.0266 | -0.0084 | **-0.0087** | -0.0053 | 1.039 | 0/3, 0/3 | -0.028 / -0.0033 |
| CONTROL `Kinf_shuf` | -0.308 | -0.597 | -0.0001 | – | +0.1865 | +0.1794 | -0.0327 | **-0.0314** | -0.0308 | 0.988 | 0/3, – | -0.558 / -0.0003 |
| CONTROL `Kinf+gf_shuf@0.25` | +0.037 | -0.253 | -0.0270 | – | +0.0790 | +0.0741 | -0.0238 | **-0.0229** | -0.0211 | 0.999 | 0/3, – | -0.135 / -0.0049 |
| `sh@1.15` (no key) | +0.085 | -0.205 | +0.0118 | +0.0163 | +0.0640 | +0.0582 | -0.0069 | **-0.0042** | -0.0025 | 0.998 | 0/3, 3/3 | -0.389 / +0.0046 |
| `base` (no key) | 0 | -0.290 | +0.0005 | +0.0052 | +0.0904 | +0.0854 | -0.0157 | **-0.0130** | -0.0115 | 1.000 | 0/3, 3/3 | -0.413 / +0.0004 |

Notes:
- `gi@t` = 228's placement: gain on the whitened base before the residuals, median-row-norm preserve.
- `gf@t` = 410d's placement: gain on the final bulk row, Frobenius preserve.
- `K<lo>_<hi>` = the key at clip(lo, hi).
- `sh@q` = h-space sharpening with row renormalisation; `shL` = lfc space.
- Full tables are in `work/analyse414.log`, fold-level rows in `work/{arms,board,extra,ballast}414-<fold>.json`.

## 4. (a) 228's gain on top of the key: the direction is real, but both placements carry an amplitude excursion

**Both placements.** On cancer folds: (i) `gi@0.25+Kinf` rel +0.030 vs 410c, pds +0.0051; (ii) `Kinf+gf@0.25` rel
+0.051, pds +0.0044. The direction gain holds at every t ≤ 0.5, and paired bootstrap on hct116/hek293t puts P(Δrel>0) at
1.00. It stacks with the key (Kinf alone +0.290 → +0.341).

**The problem is where the energy lands.** Both placements preserve a norm in EFFECT space, and the render zeroes every
gene below each context's CPM gate. m is largest on expressed, high-reliability genes. On kolf2, the energy share on
genes gated out of every challenge context falls from 18.0 % to 6.3 %, and the top-1 % m genes go from 2.6 % to 13.8 % of
the energy (`gainconst414.py`). The render-space energy (what the scorer's A_c/B_c read) therefore rises:

| fold set | render-space energy vs 410c (A/B/C) | cos ratio | Δraw mse vs 410c (cliff) | Δoverall | amplitude-matched Δoverall |
|---|---|---:|---:|---:|---:|
| cancer (fold-specific gains) | 1.084 / 1.092 / 1.096 | 1.040 | **+0.0272** | **−0.0031** | +0.0046 |
| kolf2 (= **the shipped 410d m**, max \|d\| 0; energy constant 0.545 vs the build's 0.5325) | 1.191 / 1.193 / 1.199 | 1.045 | **+0.117** | **−0.021** | +0.0025 |

With ρ_A already 1.043, context A pays the whole excursion. B and C are pushed over their cliffs too (ρ after 410d:
1.13/1.06/1.02 cancer, 1.24/1.16/1.11 kolf2). 228 saw exactly this (+7.4 % delivered row norm) and priced the gain only at
@adel α 0.889. 410d's preregistration quotes 228's @adel price but builds the unmatched arm. `price410d.py` has the numbers.

**Pre-upload check:**

    uv run python experiments/414-bulk-reprice/cliffcheck414.py <artifact.h5ad>

- 10 s, < 3 GB RSS, read-only.
- Validated on 410c: it measures render-space energy vs 410b of 1.0099/1.0084/1.0149, against PLOLO's predicted
  1.013/1.011/1.012.
- If 410d reads ≳ 1.03 over 410c in any context, the upload measures an amplitude arm, not the gain.

The honest gain arm is amplitude-matched: bulk × 0.960/0.957/0.955 (cancer) or 0.916/0.915/0.913 (kolf2) per context.
It is predicted at +0.0025 to +0.0046 overall. That is a builder edit and is not mine to make.

**Measured on the built 410d artifact** (`cliffcheck414.py`, `diag410d.py`, `mass410d.py`, all read-only streams):

| context | A_c 410d / 410c | ρ_c 410c → 410d | amplitude-only Δraw vs 410b | where the increase is |
|---|---:|---|---:|---|
| A | ×1.197 | 1.043 → 1.246 | +0.221 | spread (median target ×0.99; top-10 genes 17 %: SDF2L1, MANF, HERPUD1 …); SYVN1 row 195 → 619 |
| B | ×1.182 | 0.971 → 1.147 | +0.155 | same pattern; SYVN1 186 → 506 |
| C | **×2.127** | 0.930 → **1.978** | **+1.124** | **SYVN1 229 → 4,858**; its block is 50.5 % HSP90B1 + 12.7 % HSPA5 counts |

SYVN1 knockdown triggers the unfolded protein response, which is very reproducible across lines, so the reliability gain
m is 4–5 on exactly those genes. HSP90B1 is at 668 CPM in context C. Boosted ×2.2 it takes half the block, and the mass
balance pushes every other gene down about one log unit.

**410d2 emulation** (`emul410d2.py`): each variant is m applied to 410c's rendered bulk, then re-rendered.

| variant | constant | energy vs 410c A / B / C | largest row ratio |
|---|---:|---|---|
| 410d (KNOWN ANSWER) | 0.5325 | ×1.166 / ×1.174 / ×1.180, vs measured ×1.197/×1.182/×2.127 → **fails in C** | ×1.6 |
| 410d2 (render-gated) | 0.4846 | ×0.971 / ×0.978 / ×0.984 (C unreliable) | ×1.3 |
| m clipped to [1/3, 3] | 0.517 | ×0.987 / ×0.995 / ×0.994 | ×1.2 |
| **m clipped to [1/2, 2]** | 0.616 | ×0.994 / ×1.002 / ×0.995 | ×1.1 |

The rendered stand-in understates SYVN1's true effect, because the render's mass fix already compressed it. The
emulation is therefore trustworthy only in A and B. **The GO rule for any 410d2: cliffcheck energy ratio 0.97…1.03 in
every context AND diag max per-target ratio ≲ 1.5.**

## 5. (b) Sharpening and the other closed magnitude reshapes

**169's sharpening (h space, row renormalisation).** It was closed only on nmae: HCT116 −0.00023 and h1-test −0.01289.
Ballast removes that cost. On the 410c bulk:

- Dose response is monotone to q = 1.5 on every fold. Rel cos over Kinf: +0.029/+0.058/+0.086/+0.140/+0.188/+0.250 at
  q 1.05/1.1/1.15/1.25/1.35/1.5.
- pds rises at every dose (1.05/1.1/1.15/1.25/1.35/1.5): +0.0040/+0.0079/+0.0117/+0.0190/+0.0258/+0.0346 raw.
- Per fold at q = 1.25, vs Kinf: cosine ×1.103/×1.111/×1.111 (k562/hct116/hek293t); pds +0.018/+0.018/+0.021 literal,
  +0.019/+0.017/+0.021 board-emulated.
- Paired target bootstrap at q 1.15 and 1.25: P(Δrel>0) = 1.00 and P(Δpds>0) ≥ 0.99 on all four folds, kolf2 included.
- **kolf2 transfers at about 1/3:** q 1.25 gives +4.0 % cos and +0.0063 pds, against +10.8 % and +0.019 on cancer.
- It lowers render-space energy (0.967/0.963/0.958 at q 1.25), so it also pulls context A back toward its cliff.
- The β-emulated pds gain is the same as the literal one. The pattern that costs the key its pds does not touch
  sharpening, because row norms are unchanged.
- **Controls fail as they should, on 4/4 folds:**
  - q = 0.87 (flatten): rel −0.079, pds −0.011.
  - Within-row permuted multipliers |h|^(q−1): rel −0.016, pds −0.0019. So the gain needs the multiplier aligned with |h|.
- lfc space (shL@1.15) reads slightly below h space (+0.0085 vs +0.0096 overall).
- **Precedent to respect.** a13 was the other per-gene magnitude reweighting. It lost pds −0.029 scaled on the board, a
  bulk loss that ballast would not rescue. Its key, donor agreement, has AUC 0.498 for the reference's DE set, against
  |h|'s 0.597 (168). Its harness surfaces disagreed in sign (A13-OPERATING-POINT §0). Sharpening's four folds, including
  the leak-free 3-donor kolf2, agree in sign, and its controls go negative.

**Other reshapes closed on nmae:**
- **158 (per-target scale, donor-cosine key):** superseded by 410's key. K3s/Kinfs (scored-column inner products) ≈ K3/Kinf.
- **318 (per-target redistribution):** its own finding was that it is not separable from its shuffles. The 410 key is the
  version that separates (shuffle −0.31 rel).
- **224's diagonal:** = 228's gain, §4.
- **187 (amplitude):** §6.
- **a13's agreement shrink (096):** not re-priced. Its board pds loss was on the bulk, and its 3-donor unanimity estimator
  cannot run on the 2-donor cancer folds (A13-OPERATING-POINT §0).

## 6. (c) Amplitude on the key-weighted bulk

The model says the optimum is ON each context's cliff whatever X is. Against 395's X or 410c's X (+43 %) the argmin is the
same to 4 decimals. The optimum could move above A's cliff only if X exceeded s·A_c; here X/A_c ≈ 0.15. Slopes per unit
bulk scale on 410c are −0.31 below the cliff and +1.78/+1.74/+1.89 above it. **So the answer is no.** What moved is the
bulk: the key's render-space energy put context A 2.2 % (in scale) over its cliff.

| context | ρ_c (410c) | cliff bulk scale | α-equivalent | Δraw on context | EV-optimal (187's cliff sd) |
|---|---:|---:|---:|---:|---|
| A | 1.0454 | **0.9780** | 0.9604 | −0.0396 | 0.9755 (−0.0383) |
| B | 0.9733 | **1.0136** | 0.9954 | −0.0042 | 1.0130 (−0.0039) |
| C | 0.9270 | **1.0386** | 1.0200 | −0.0118 | 1.0370 (−0.0112) |
| mean | | | | **−0.0185 → overall +0.0032** | −0.0178 → +0.0031 |

Context A alone: +0.0023. One global scale of 0.978: +0.0015. The builder refuses `VCC_ALPHA=knee` with
`VCC_PERT_SCALE`, so a per-context bulk scale (under ballast it touches only pds/mse, and pds is scale-invariant) needs a
builder knob. Sharpening at q 1.25 already removes about 70 % of A's excess (ρ_A 1.045 → 1.011), so re-levelling after
sharpening adds only +0.0024 (table column "levelled").

## 7. Recommendation: ONE bulk arm, one variable on 410c

**`VCC_BULK_SHARPEN_Q=1.25`**

- **Definition:** `sign(h)|h|^1.25` on each context's CPM gate, each row renormalised to its own gated L2 norm (169's
  `sharp.sharpen(mode="row")`, the transform 169's unshipped `sharpen-builder.diff` defines).
- **Scope:** `bulk_effects` / the ballast profile only. The carriers keep 395's design, and nothing else changes.
- **Render safety on the REAL 410c bulk** (`sharp410c.py`, 900 blocks re-rendered with production's
  `profile_from_effect`):
  - At q 1.25, render-space energy is ×0.990/×0.988/×0.982 of 410c.
  - No row moves by more than ×1.001, and the largest single-gene mass share is unchanged (0.011/0.029/0.022).
  - At q 1.5: ×0.95/×0.95/×0.94, rows ≤ ×1.002.
  - Row renormalisation makes sharpening immune to the SYVN1-type blow-up.
  - Re-priced with these real-bulk energy ratios (`pricesharp414.py`), q 1.25 gives **+0.0055 (kolf2-like) to +0.0151
    (cancer)** overall, against +0.0058…+0.0158 on PLOLO's energy ratios.
- **Build checks:**
  - The config diff vs 410c's meta is exactly one key.
  - `cliffcheck414.py` must show render-space energy 0.97…1.00 of 410c in every context.
  - `diag410d.py <410c> <arm>` must show the max per-target energy ratio ≲ 1.1.
  - `ballast_jk_ratio` stays at 410c's level. About 15 % of blocks will hit K = 16, against 8 %.

**Predictions vs 410c:**
- pds +0.006…+0.019 raw.
- raw mse −0.020…−0.050.
- DE members within noise (|Δnmae| ≤ 0.010, |Δreach| ≤ 0.012, |Δfid| ≤ 0.003, |Δjac| ≤ 0.001, the GO-410c rule).
- **overall +0.0055…+0.0151.**

**Bins (Δoverall vs 410c):**
- ≥ +0.006: transfers at least at kolf2's rate. Next arm q 1.5 (predicted +0.009…+0.026).
- 0…+0.006: partial transfer. Hold the dose and read which member moved.
- < 0: PLOLO's magnitude structure does not reach the board. Close the class. A pds loss with an mse gain would repeat a13.

Why q 1.25 and not 1.5: 1.5 is predicted higher, but it doubles the K-capped blocks (21.6 %) and the realisation loss
(−2.7 % rel). At 1.25 a negative result is still decisive for the class.

**Secondary arms, in order:**
1. Floor lo 0.5 (`K0.5_inf`): +0.0012, free on cosine.
2. Per-context bulk re-levelling to the cliffs: +0.0032; after sharpening, +0.0024.
3. The amplitude-matched gene gain: +0.0025…+0.0046.

## 8. Limits

- **kolf2 is the only non-Replogle/Orion truth.** It reads sharpening at about 1/3 of the cancer folds, so the board could
  sit anywhere in that range. k_key (1.407) is measured for the key class only. The sharpening and gain classes carry k as
  an unknown between kolf2-like and cancer-like.
- **The board-emulated pds is one parameter (λ) fitted to one board point** (410b − 395). It is a mechanism emulation, not
  a validated instrument. It is uninterpretable on kolf2, where β aligns with the truth.
- **The cliff model's claim_c are 187's retrodicted proxies**, assumed saturated (410b's carriers' jackknife margin ≥ 1.5 on
  H1, 1.58–1.86 on A/B/C). A_c/B_c for an arm come from PLOLO ratios on 410b's measured values. The ratio route was
  validated on 410c's artifact to ≤ 0.3 %.
- **The composition is 228's** (CD4 0.34 + H1 union 0.28), not 395's full stack (CD4 deflate k5, normrel cap, sign prior
  on carriers). The bulk-side differences are small, and the 410c energy check confirms the ratio transfers.
- **Sharpening's pds on h1-test's DE-table surface (168) was about 10× smaller than on PLOLO.** 169 argued that surface is a
  different statistic, and H1 also mis-signed the key's board pds (+0.0046 vs −0.0036). H1 was not re-run here.
- **Peak RSS exceeded 8 GB twice** (8.58 GB extra414, 8.73 GB ballast414; four folds in one process). All other runs:
  arms414 3.95 GB (one fold) / 7.29 GB (three folds), board414 6.71 GB, decomp414/shared414/cliffcheck414 < 2.8 GB.

## 9. Files and reproduction

```
cd experiments/414-bulk-reprice            # all scripts write only into ./work
uv run python arms414.py                   # KA-1/2/3 + the main arm grid, 4 folds (-> arms414-<fold>.json, arrays414-*.npz)
uv run python decomp414.py                 # A_c/B_c on the 395 / 410 / 410b artifacts (read-only stream)
uv run python shared414.py                 # the non-scaling shared pattern beta (-> beta414.npz)
uv run python board414.py                  # board-emulated pds, lam 0/0.5/1/2
uv run python extra414.py                  # floors, sharpening controls and doses
uv run python ballast414.py                # ballast realisability of each design
uv run python analyse414.py                # both mse routes, bootstrap, amplitude (c) -> analyse414.{json,log}
uv run python price410d.py "Kinf+gf@0.25" "Kinf+sh@1.25"   # arm vs 410c on cancer and kolf2 energy ratios
uv run python gainconst414.py kolf2        # the gain's effect-space energy constant vs the build's
uv run python cliffcheck414.py <artifact.h5ad>             # pre-upload: render-space energy + rho_c per context (< 3 GB)
uv run python diag410d.py <ref.h5ad> <artifact.h5ad>       # pre-upload: per-gene / per-target energy vs a reference
uv run python mass410d.py                  # SYVN1 blocks' top-gene mass shares, 410c vs 410d
uv run python sharp410c.py                 # sharpening render safety on the real 410c bulk
uv run python pricesharp414.py             # sharpening arms vs 410c with real-bulk energy ratios
uv run python emul410d2.py                 # 410d / 410d2 / clipped-m emulation on the real 410c bulk (KA fails in C)
uv run python table414.py                  # the §3 table
```
Run heavy PLOLO scripts one fold per process (`arms414.py --folds <f>`, `extra414.py <f>`, `ballast414.py <f>`) to stay
under 8 GB. All artifact reads are read-only h5py streams.
