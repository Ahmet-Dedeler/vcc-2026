# 415 — a leave-one-donor-out per-target key for the production bulk row: it beats 410's key on every fold

Status: **COMPLETE (PLOLO + A/B/C artifacts). Nothing built, submitted or wired into the builder.**
Preregistration: `PREREGISTRATION.md` (decision rule fixed before any LODO number existed).

**One-line recommendation:** swap 410's key for the LODO-V1 key in the next bulk arm. Use
`work/pert_scale_415pp_lo0.3_hi3.npz` (sha `5f1ae21f…`), built on the composition the builder actually scales. It beats
410's key on all three cancer folds and on kolf2, on the shipped-composition prediction: +1.2 … +5.4 % relative weighted
cosine at clip(0.3,3), +4.3 … +6.1 % at clip(0.3,∞). My rough estimate is about +0.002 overall on the board at clip 3,
but that figure is an extrapolation, not a measurement (§5). Whether it is worth a slot is 414's call.

## 1. The key
`key_p = mean_j <P^{−j}_p, d_{j,p}> / |P_p|²` over the source donors j that carry p, where:
- carried = at least 2 donors, the same set as 410; non-carried targets get w = 1.
- `d_j` = donor j's raw bank row (410's `donor_rows`).
- `P` = the full production-composition bulk row: whitened base + CD4 residual + H1-union residual.
- `P^{−j}` = the same row with the base rebuilt from the other donors only. The global constants are held at P's
  values: on-axis preserve c, residual fractions, post-add preserve s.
- `w = clip(key / median(key over carried), lo, hi)`.

With fixed constants, E<P^{−j}, d_j> = E<P, T> under donor exchangeability, so the key estimates
`s*_p = <P_p,T_p>/|P_p|²` for the whole row, not just its donor base. The CD4 and H1 parts are separate data sources, so
leaving a donor out cannot correlate them with that donor's noise. KA-4 checks the algebra: with 2 sources, no residual
and d = the rows W acts on, the LODO key equals 410's pairwise formula to 3e-15.

## 2. Known answers (every run asserts them; all PASS)
| check | result |
|---|---|
| KA-1: 410's key on 162's original fold predictions reproduces `check410.json` | rho +0.576210 / +0.600988 / +0.627041 / kolf2 +0.313594; gains at clip(0.3,3) +0.298732 / +0.267210 / +0.286726 / +0.499810 and clip(0.3,∞) +0.385891 / +0.341318 / +0.350060 / +0.936098. All |Δ| < 1e-9 |
| KA-2: re-prepared donor rows (`prep415.py`) rebuild the cached fold `pre` | rel max |Δ| 3.0e-8 … 5.3e-8 (float32 cache) |
| KA-3: the fold production composition reproduces 228/223 `base222` (a23 cos, pds) | Δ = 0.0e+00 on 4/4 folds |
| KA-4: LODO(base only, 2 sources) = 410 formula on the same rows | rel 3.1e-15 … 4.4e-15 |
| KA-410 (A/B/C): 410's key/weights recomputed = the shipped 410 artifacts (d7f76a22…, 94c048e9…) | max |Δkey| = max |Δw| = 0 |
| KA-PP (A/B/C): my reconstruction of the SHIPPED composition reproduces the 410b build meta | whitening medians, operator preserve 0.158168, fractions CD4 0.346141 / H1 0.194939 / partial 0.167823, post-add preserve 0.904698: worst relative error **8.9e-6** (float32 model rows). CD4 SVD-5 deflation **bit-exact**: U_k sha `f06e0262…`, energy removed 0.19866046506100588 |

## 3. PLOLO on 228's production composition (the preregistered test, `lodo415.py`, `work/lodo415.json`)
The relative aggregate transfer-cosine gain is shown against the unweighted row. In brackets is a 1,000× target
bootstrap. "Δ vs 410" is the paired difference, with its own bootstrap CI. "rel" is (1+g)/(1+g410) − 1, the
improvement on top of 410's weighting.

| fold | key | rho(key,s*) | gain clip(0.3,3) | Δ vs 410 [CI] | rel | gain clip(0.3,∞) | Δ vs 410 [CI] | rel |
|---|---|---:|---:|---|---:|---:|---|---:|
| k562_gwps | K410 | +0.572 | +0.244 [+0.166,+0.333] | | | +0.280 [+0.170,+0.395] | | |
| | **LODO-V1** | **+0.635** | **+0.318** | **+0.0735 [+0.0449,+0.1097]** | +5.9 % | **+0.341** | **+0.0613 [+0.0290,+0.0961]** | +4.8 % |
| orion_hct116 | K410 | +0.601 | +0.283 [+0.182,+0.368] | | | +0.335 [+0.216,+0.445] | | |
| | **LODO-V1** | **+0.644** | **+0.302** | **+0.0185 [−0.0019,+0.0456]** | +1.4 % | **+0.389** | **+0.0540 [+0.0151,+0.0944]** | +4.1 % |
| orion_hek293t | K410 | +0.619 | +0.293 [+0.189,+0.408] | | | +0.347 [+0.204,+0.500] | | |
| | **LODO-V1** | **+0.675** | **+0.325** | **+0.0318 [+0.0009,+0.0723]** | +2.5 % | **+0.397** | **+0.0499 [+0.0132,+0.0995]** | +3.7 % |
| kolf2 (not in rule) | K410 | +0.370 | +0.427 [+0.282,+0.574] | | | +0.721 [+0.401,+0.984] | | |
| | LODO-V1 | +0.426 | +0.491 | +0.0641 [+0.0100,+0.1412] | +4.5 % | +0.828 | +0.1073 [+0.0240,+0.2000] | +6.2 % |

**Decision rule: PASS.** V1 beats K410 at clip(0.3,3) on all three cancer folds, and also at clip(0.3,∞), so both
artifacts were written. Seven of the eight intervals exclude 0. The exception is orion_hct116 at clip 3, whose lower
bound is −0.0019 and whose bootstrap P(Δ > 0) is 0.963.

Context rows, for the cancer folds in order:
- Oracle: +0.576 / +0.680 / +0.744.
- Base aggregate cosine of the production row: 0.0402 / 0.0308 / 0.0245. The original whitened-only row was 0.0270 /
  0.0212 / 0.0187, so the CD4 + H1 residual itself carries transferable signal.
- Residual share of the row norm: 0.34–0.35.
- Negative control, V1 weights permuted among carried targets (200×): −0.19 at clip 3 and −0.30 at ∞ on every fold.
  The key's information is real.

**Where the gain comes from (decomposition, same folds, clip(0.3,3) Δ vs 410):**
| row | k562 | hct116 | hek293t | kolf2 | reading |
|---|---:|---:|---:|---:|---|
| D1: 410's formula on the prepared (centred) rows, base-only | +0.0053 | +0.0011 | −0.0015 | +0.0108 | row preprocessing adds nothing |
| D2: LODO numerator base-only, full-row denominator | **−0.0390** | −0.0029 | −0.0050 | **−0.0419** | the full denominator alone is WORSE |
| V1: + the residual term in the numerator | +0.0735 | +0.0185 | +0.0318 | +0.0641 | **the gain is <R, d_j>** |

The whole improvement comes from estimating the CD4/H1 residual's projection on the truth, which is the mandate's
hypothesis. Normalising by the full row without that term is harmful. That is the same failure mode as key2_410, a
mismatch between numerator and prediction.

Secondary variants (reported, not selected per the preregistration):
- V2 (d_j = centred, prepared rows): +0.0044 / +0.0099 / +0.0175 over V1 at clip 3.
- V3 (full pipeline rebuild per j): about the same as V2.

Both are small and inside the CIs. Centred d_j is a cheap preregistered follow-up if anyone wants the last 1 %.

## 4. Does it matter that the shipped row is not 228's composition? (`sens415.py`, `work/sens415.json`)
The 410b row differs from 228's `build_final` in four ways:
- normrel cap −1.5 donor weighting;
- PM50 applied on every PerturbMatch-observed gene (061's broad policy), not only on the columns K562 and PM share;
- panel columns kept in the base, which moves every all-rows sizing median;
- CD4 SVD-5 deflation plus CD4 partial rows at g 0.5.

Consequences for 228's composition on A/B/C: frac_cd4 is 0.314 against the shipped 0.346, and H1 is 0.187 against
0.195. "PP" is my reconstruction of the shipped composition, verified by KA-PP above. The next table scores three keys
on the **PP** fold prediction (median per-target cosine of the 228 row with the PP row: 0.960–0.991).

| fold | K410 clip3 | V1(228-comp) clip3 Δ410 | **V1(PP) clip3 Δ410 [CI]** | K410 ∞ | V1(228) ∞ Δ410 | **V1(PP) ∞ Δ410 [CI]** |
|---|---:|---:|---|---:|---:|---|
| k562_gwps | +0.240 | +0.0695 | **+0.0666 [+0.0377,+0.0982]** | +0.261 | +0.0564 | **+0.0664 [+0.0362,+0.0973]** |
| orion_hct116 | +0.295 | +0.0166 | **+0.0160 [−0.0072,+0.0409]** | +0.356 | +0.0513 | **+0.0585 [+0.0091,+0.1029]** |
| orion_hek293t | +0.263 | +0.0301 | **+0.0413 [+0.0176,+0.0738]** | +0.310 | +0.0688 | **+0.0803 [+0.0363,+0.1358]** |
| kolf2 | +0.411 | +0.0634 | **+0.1027 [+0.0410,+0.1911]** | +0.689 | +0.0976 | **+0.1103 [+0.0326,+0.2027]** |

- The result holds on the shipped composition. Both LODO keys beat K410 on every fold at both clip rules.
- V1(228) − V1(PP) is within ±0.012 on the cancer folds, and every CI spans 0. At clip ∞, PP is ahead on all four folds
  (0.007 to 0.013). On kolf2 at clip 3, PP is ahead by 0.039 [0.012, 0.071].
- On A/B/C, Spearman(key) is: PP vs the 228-comp artifact 0.973; PP vs 410 0.943; 228-comp artifact vs 410 0.947.
  Weight-vector cosine on the PP row energies at clip 3: PP vs artifact 0.991, PP vs 410 0.977.
- 19 % of targets' weights move by more than 25 % between the two LODO artifacts, and 31 % between PP and 410.

The preregistered artifact (228 composition) passed its rule. The PP artifact uses the same key definition on the row
the builder actually multiplies (`apply_pert_scale` acts on the finished shared matrix, which KA-PP reproduces), and on
PLOLO it passes the same rule. **Recommend the PP artifact.**

## 5. What it could be worth on the board (my extrapolation, NOT measured; 414 prices arms)
410 at clip(0.3,3) turned PLOLO +27…+30 % into board raw mse −0.0797, which is about +0.038 aggregate cosine
(0.107 → 0.146).
- **PP key at clip(0.3,3):** it adds rel +5.4 / +1.2 / +3.3 % (mean +3.3 %) on top of 410's weighting. That is about
  +0.005 cosine, about −0.010 raw mse, about **+0.0017 overall**. pds, and the DE channel under ballast, should not move
  in first order.
- **PP key at clip(0.3,∞) against 410b's clip 3:** rel +7.0 / +9.3 / +10.1 %, about **+0.005 overall**. The caveat is
  large: uncapped weights reach 10× (TAF4), which stresses ballast K and the jackknife credit, and 410 measured the
  weighted-row pds loss at −0.0036 raw regardless of the key. This is 414's territory.

The transfer caveat is specific to this key. The residual term uses cancer donors as the proxy for the truth's
projection on the CD4 residual. On a T-cell-like context (A), the true CD4 projection may be larger than that estimate,
so this key may be conservative there. kolf2 (iPSC, non-cancer) still gains +0.064 … +0.110, which is some reassurance
that it is not cancer-specific.

## 6. Artifacts (410's exact format: `targets, w, key, carried, n_donors` npz + json sidecar)
All are on the A/B/C panel with operator union_pruned `e2f65ecf…` and 410's Oct-22 axis guard verbatim. The builder is
untouched. Wiring an artifact into `PERT_SCALE_ARTIFACTS` is a separate, reviewed step, and the pinned sha must be the
sha256 of the npz file.

| file | sha256 | composition | clip | w q0/10/25/50/75/90/100 | floor / cap |
|---|---|---|---|---|---|
| `work/pert_scale_415pp_lo0.3_hi3.npz` **(recommended)** | `5f1ae21f530734735b23bc0a114bb7b84988f0b62221c6a11c7a64d8176547d9` | shipped (PP) | (0.3,3) | .3/.3/.311/1/2.07/3/3 | 24.5 % / 16.7 % |
| `work/pert_scale_415pp_lo0.3_hi1e+09.npz` | `5a12079f83b607bb93e14ea30aa4a8443c56deb426d26a0c7f5854a019296dfd` | shipped (PP) | (0.3,∞) | … max 9.99 | 24.5 % / – |
| `work/pert_scale_415_lo0.3_hi3.npz` (preregistered) | `845e3895f81ea0a1e771e8d7bfbb7e2e406f6cff9ad2d36f482d2b5188a4d739` | 228 | (0.3,3) | .3/.3/.3/1/2.28/3/3 | 25.9 % / 17.3 % |
| `work/pert_scale_415_lo0.3_hi1e+09.npz` (preregistered) | `ef32aa8201d679b2c3d98827fca61ec208150dfd6e5921d353a7b12a23fdbcd0` | 228 | (0.3,∞) | … max 12.6 | 25.9 % / – |

Carried 294/300: 259 targets have 3 donors, 35 have 2, and 6 have 1 (w = 1). The top of the PP key is TAF4 9.99×,
SYVN1, SIN3B, ACLY, UPF3B, SMARCA5 and PQBP1; the bottom is RAB11FIP5 and several ZNFs at the floor. The sidecar
`*.json` files hold the composition constants, the relative error against the 410b meta, and the comparisons against
410 and the 228-composition artifact.

## 7. Files
- Scripts: `prep415.py` caches the prepared donor rows; `lodo415.py` runs the PLOLO comparison with KA-1…4;
  `build_weights415.py` builds the 228-composition artifacts with KA-410; `sens415.py` holds the PP composition, KA-PP
  and the PLOLO/A/B/C sensitivity; `build_weights415pp.py` builds the PP artifacts with KA-410 and KA-PP.
- Logs and JSON are in `work/`.
- Caches: `work/prep-*.npz` and `work/full-*.npz` (8 × 44 MB, gitignored, regenerable in about 30 s). Delete them when
  convenient.
- Peak RSS ≤ 7.0 GiB on every run: the first `sens415` run streamed the Orion banks; later runs stayed at ≤ 6.1 GiB.
  No builder edit, no vcc, no cell-eval2. Everything was written inside `experiments/415-lodo-key/`.
