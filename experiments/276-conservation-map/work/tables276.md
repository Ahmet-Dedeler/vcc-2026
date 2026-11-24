## KA-3 sensitivity — 262's `agree_donor_donor` recomputed on three gene axes

| axis | genes | rho(donor-donor, pds), pooled | folds pooled |
|---|---:|---:|---|
| 276's common axis (26 sources) | 2994 | +0.306 | k562_gwps, orion_hct116, orion_hek293t |
| the three donors' own shared axis | 7543 | +0.403 | k562_gwps, orion_hct116, orion_hek293t |
| the two Orion lines' shared axis | 18106 | +0.236 | k562_gwps |

262 reports +0.381 pooled (R262 §4.1, `analysis262.json`).

## Per-fold Spearman with `cd4svd5`'s per-target `pds`

| predictor | K562 fold | HCT116 fold | HEK293T fold | kolf2 fold | pooled (3 cancer) |
|---|---:|---:|---:|---:|---:|
| C_raw (all 26 sources) | +0.508 | +0.497 | +0.469 | +0.366 | **+0.488** [+0.432, +0.540] |
| C_lolo (held line removed) | +0.411 | +0.390 | +0.365 | +0.278 | **+0.387** [+0.326, +0.444] |
| C_rel (reliability-corrected, 4 sources) | +0.419 | +0.476 | +0.497 | +0.165 | **+0.461** [+0.398, +0.522] |
| C_rel_lolo (reliability-corrected, held line out) | +0.210 | +0.265 | +0.373 | +0.165 | **+0.290** [+0.208, +0.379] |
| baseline: donor-donor agreement | +0.207 | +0.322 | +0.373 | +0.179 | **+0.307** [+0.239, +0.372] |
| baseline: CD4 condition count | +0.089 | +0.044 | +0.085 | -0.037 | **+0.056** [-0.010, +0.121] |
| baseline: target expression (control log1p) | +0.052 | +0.046 | +0.053 | +0.060 | **+0.051** [-0.015, +0.118] |
| control: measured reliability of the non-held sources | +0.200 | +0.238 | +0.390 | +0.210 | **+0.283** [+0.215, +0.346] |
| control: response magnitude of the non-held sources | +0.158 | +0.236 | +0.219 | +0.203 | **+0.207** [+0.144, +0.269] |
| control: number of pairs | +0.235 | +0.184 | +0.171 | +0.322 | **+0.197** [+0.130, +0.262] |

## What a conservation score buys, read as groups (LOLO, pooled over 3 cancer folds)

| `C_lolo` band | n (target, fold) | mean `pds` | share anti-aligned | median `C_lolo` |
|---|---:|---:|---:|---:|
| bottom decile | 85 | 0.647 | 36.5 % | -0.0059 |
| 10–30 % | 170 | 0.706 | 27.1 % | +0.0026 |
| 30–50 % | 170 | 0.794 | 19.4 % | +0.0110 |
| 50–70 % | 170 | 0.822 | 12.4 % | +0.0234 |
| 70–90 % | 170 | 0.901 | 7.6 % | +0.0432 |
| top decile | 85 | 0.895 | 4.7 % | +0.0786 |
| all | 850 | 0.799 | 17.4 % | +0.0171 |

## The same on h1-test (85 Flex-400 targets, H1 sources removed)

| `C_lolo` tercile | n | mean `pds` | median `C_lolo` |
|---|---:|---:|---:|
| bottom | 28 | 0.804 | +0.0027 |
| middle | 29 | 0.900 | +0.0380 |
| top | 28 | 0.924 | +0.0781 |

## The distribution of `C_raw` over the 300 panel targets

- cross-line, raw: mean +0.0263, median +0.0199, p10 -0.0028, p90 +0.0676, share <= 0 16.0 %
- same-line pairs (reproducibility, not conservation): mean +0.3274 on 296 targets, 10 pairs
- target-shuffled null: mean -0.00034, sd over targets 0.00133, mean 95 % band [-0.0128, +0.0122]
- reliability-corrected (4 sources with measured halves): median +0.228, p10 -0.209, p90 +0.692, share <= 0 22.5 %, n 227

## Median cross-line pair cosine, by pair (the 20 pairs with the most targets)

| pair | targets | median cosine |
|---|---:|---:|
| orion_hct116 (HCT116) — orion_hek293t (HEK293T) | 289 | +0.0474 |
| orion_hek293t (HEK293T) — cd4_stim8hr (CD4T) | 282 | +0.0192 |
| orion_hek293t (HEK293T) — cd4_stim48hr (CD4T) | 279 | +0.0215 |
| orion_hct116 (HCT116) — cd4_stim8hr (CD4T) | 278 | +0.0231 |
| orion_hek293t (HEK293T) — cd4_rest (CD4T) | 278 | +0.0164 |
| orion_hct116 (HCT116) — cd4_stim48hr (CD4T) | 275 | +0.0247 |
| orion_hct116 (HCT116) — cd4_rest (CD4T) | 274 | +0.0249 |
| orion_hek293t (HEK293T) — kolf2 (KOLF2.1J) | 268 | +0.0094 |
| orion_hct116 (HCT116) — kolf2 (KOLF2.1J) | 264 | +0.0096 |
| k562_gwps (K562) — orion_hek293t (HEK293T) | 263 | +0.0334 |
| kolf2 (KOLF2.1J) — cd4_stim8hr (CD4T) | 263 | +0.0136 |
| k562_gwps (K562) — orion_hct116 (HCT116) | 260 | +0.0334 |
| kolf2 (KOLF2.1J) — cd4_rest (CD4T) | 260 | +0.0156 |
| k562_gwps (K562) — cd4_stim8hr (CD4T) | 260 | +0.0304 |
| kolf2 (KOLF2.1J) — cd4_stim48hr (CD4T) | 259 | +0.0211 |
| k562_gwps (K562) — cd4_stim48hr (CD4T) | 258 | +0.0288 |
| k562_gwps (K562) — cd4_rest (CD4T) | 257 | +0.0290 |
| k562_gwps (K562) — kolf2 (KOLF2.1J) | 249 | +0.0170 |
| perturbmatch (K562) — orion_hek293t (HEK293T) | 245 | +0.0321 |
| perturbmatch (K562) — orion_hct116 (HCT116) | 242 | +0.0308 |
