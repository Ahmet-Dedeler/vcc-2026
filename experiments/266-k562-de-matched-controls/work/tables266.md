## KA — pooled seed 0 (260's exact control draw) vs `candidates_de_power.csv`

| field | targets identical / 50 |
|---|---:|
| `n_real_tc` | 50 |
| `n_real_uc` | 50 |
| `gate_tc` | 50 |
| `gate_uc` | 50 |
| `cells` | 50 |
| `frac_neg_tc` (|Δ| < 1e-9) | 50 |
| `target_knockdown_lfc` (|Δ| < 1e-4; means in float64 here) | 50 |
| `median_umi` | 50 |

Control median UMI (pooled s0) 11406 (260: 11,406). Genes passing the CPM > 5 gate: [np.int64(8248)] in every run; min control-mean CPM over all 450 runs 5.28.

## Per-target `n_real` (tie-corrected, target gene excluded) — 50 candidates, three schemes × three seeds

| target | cells | pooled s0 (=260) | pooled s1 | pooled s2 | gem s0 | gem s1 | gem s2 | gem+depth s0 | s1 | s2 | gem mean − pooled mean | gem+depth mean − pooled mean |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DCP2 | 241 | 346 | 389 | 387 | 377 | 379 | 370 | 315 | 314 | 313 | +1.3 | -60.0 |
| ANP32B | 260 | 134 | 154 | 136 | 152 | 147 | 151 | 179 | 174 | 176 | +8.7 | +35.0 |
| DNTTIP1 | 267 | 113 | 109 | 121 | 120 | 132 | 125 | 120 | 120 | 120 | +11.3 | +5.7 |
| FCGRT | 224 | 103 | 93 | 97 | 85 | 90 | 92 | 107 | 105 | 107 | -8.7 | +8.7 |
| HMGXB3 | 342 | 85 | 92 | 85 | 95 | 91 | 96 | 101 | 102 | 102 | +6.7 | +14.3 |
| PHF23 | 261 | 80 | 73 | 72 | 71 | 70 | 79 | 88 | 88 | 90 | -1.7 | +13.7 |
| EPC2 | 259 | 77 | 76 | 80 | 83 | 82 | 77 | 79 | 81 | 79 | +3.0 | +2.0 |
| MED13 | 268 | 65 | 60 | 57 | 68 | 63 | 58 | 72 | 69 | 67 | +2.3 | +8.7 |
| PCGF1 | 260 | 64 | 64 | 63 | 68 | 75 | 69 | 68 | 69 | 66 | +7.0 | +4.0 |
| FOXN2 | 232 | 63 | 64 | 65 | 66 | 67 | 69 | 71 | 70 | 70 | +3.3 | +6.3 |
| PAN2 | 215 | 59 | 55 | 57 | 56 | 57 | 53 | 60 | 60 | 60 | -1.7 | +3.0 |
| GIGYF2 | 251 | 54 | 49 | 51 | 60 | 56 | 61 | 66 | 66 | 64 | +7.7 | +14.0 |
| ZNF462 | 300 | 50 | 63 | 55 | 61 | 61 | 57 | 63 | 66 | 64 | +3.7 | +8.3 |
| E2F1 | 226 | 31 | 26 | 26 | 23 | 24 | 25 | 30 | 29 | 30 | -3.7 | +2.0 |
| FOXO4 | 293 | 30 | 25 | 29 | 24 | 28 | 21 | 27 | 26 | 27 | -3.7 | -1.3 |
| ZNF507 | 279 | 25 | 25 | 26 | 24 | 23 | 23 | 21 | 23 | 21 | -2.0 | -3.7 |
| PPIF | 256 | 24 | 20 | 19 | 21 | 17 | 20 | 20 | 20 | 21 | -1.7 | -0.7 |
| GPATCH8 | 284 | 9 | 10 | 8 | 9 | 6 | 8 | 9 | 9 | 9 | -1.3 | +0.0 |
| NABP2 | 400 | 9 | 9 | 10 | 7 | 6 | 6 | 7 | 7 | 7 | -3.0 | -2.3 |
| NAIF1 | 264 | 8 | 8 | 8 | 8 | 8 | 8 | 8 | 9 | 8 | +0.0 | +0.3 |
| ZHX3 | 261 | 7 | 7 | 6 | 7 | 7 | 7 | 8 | 8 | 8 | +0.3 | +1.3 |
| API5 | 206 | 7 | 7 | 7 | 11 | 10 | 7 | 7 | 7 | 7 | +2.3 | +0.0 |
| MYRF | 400 | 7 | 4 | 4 | 2 | 3 | 2 | 3 | 3 | 3 | -2.7 | -2.0 |
| TRIM35 | 382 | 6 | 6 | 6 | 6 | 6 | 5 | 5 | 5 | 5 | -0.3 | -1.0 |
| THAP7 | 299 | 6 | 7 | 7 | 8 | 7 | 8 | 8 | 8 | 8 | +1.0 | +1.3 |
| ZNF821 | 336 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | +0.0 | +0.0 |
| RABEP2 | 254 | 5 | 6 | 5 | 6 | 5 | 6 | 5 | 5 | 5 | +0.3 | -0.3 |
| SNRK | 224 | 5 | 4 | 4 | 4 | 7 | 4 | 10 | 10 | 10 | +0.7 | +5.7 |
| POR | 260 | 5 | 6 | 7 | 5 | 4 | 5 | 4 | 4 | 4 | -1.3 | -2.0 |
| BRPF1 | 230 | 4 | 4 | 4 | 4 | 4 | 4 | 3 | 3 | 3 | +0.0 | -1.0 |
| MBD4 | 200 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | +0.0 | +0.0 |
| ZNF526 | 212 | 3 | 3 | 3 | 4 | 3 | 3 | 3 | 3 | 3 | +0.3 | +0.0 |
| ZNF282 | 362 | 3 | 3 | 5 | 3 | 3 | 3 | 4 | 4 | 4 | -0.7 | +0.3 |
| ZNF251 | 219 | 3 | 3 | 3 | 2 | 3 | 2 | 3 | 3 | 3 | -0.7 | +0.0 |
| UROS | 217 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | +0.0 | +0.0 |
| ZKSCAN1 | 323 | 2 | 2 | 1 | 2 | 1 | 1 | 7 | 7 | 6 | -0.3 | +5.0 |
| TFCP2 | 330 | 2 | 2 | 3 | 2 | 2 | 2 | 2 | 2 | 2 | -0.3 | -0.3 |
| GLI3 | 324 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | +0.0 | +0.0 |
| ZNF530 | 264 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | +0.0 | +0.0 |
| ANKZF1 | 252 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | +0.0 | +0.0 |
| DTNBP1 | 327 | 1 | 1 | 1 | 0 | 2 | 1 | 2 | 1 | 1 | +0.0 | +0.3 |
| IFNGR2 | 295 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | +0.0 | +0.0 |
| SMYD2 | 205 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | +0.0 | +0.0 |
| RNF141 | 314 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | +0.0 | +0.0 |
| TWF2 | 284 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | -0.7 | -0.7 |
| CYTH2 | 203 | 0 | 2 | 0 | 0 | 0 | 0 | 2 | 2 | 2 | -0.7 | +1.3 |
| ZNF514 | 300 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | +0.0 | +0.0 |
| SMARCC2 | 272 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | +0.0 | +0.0 |
| ZNF543 | 211 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | +0.0 | +0.0 |
| MED15 | 306 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | +0.0 | +0.0 |

## Medians and gate fractions per scheme and seed (50 candidates)

| scheme | seed | median `n_real` | IQR | ≥ 10 (`nmae` gate) | ≥ 20 | ≥ 30 | ≥ 89 | zero | median paired Δ vs pooled s0 | targets up / down / tie vs pooled s0 |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|
| pooled | 0 | 6.0 | [2, 45] | 17 (34 %) | 17 | 15 | 4 | 6 | +0.0 | 0 / 0 / 50 |
| pooled | 1 | 6.0 | [2, 43] | 18 (36 %) | 17 | 13 | 5 | 4 | +0.0 | 11 / 12 / 27 |
| pooled | 2 | 6.0 | [1, 46] | 18 (36 %) | 16 | 13 | 4 | 5 | +0.0 | 13 / 14 / 23 |
| gem | 0 | 6.0 | [2, 48] | 18 (36 %) | 17 | 13 | 4 | 7 | +0.0 | 14 / 12 / 24 |
| gem | 1 | 6.0 | [2, 49] | 18 (36 %) | 16 | 13 | 5 | 6 | +0.0 | 13 / 13 / 24 |
| gem | 2 | 5.5 | [1, 46] | 17 (34 %) | 17 | 13 | 5 | 6 | +0.0 | 10 / 15 / 25 |
| gemdepth | 0 | 6.5 | [2, 52] | 18 (36 %) | 17 | 14 | 5 | 5 | +0.0 | 19 / 10 / 21 |
| gemdepth | 1 | 6.5 | [2, 52] | 18 (36 %) | 17 | 13 | 5 | 5 | +0.0 | 19 / 10 / 21 |
| gemdepth | 2 | 6.0 | [2, 52] | 18 (36 %) | 17 | 14 | 6 | 5 | +0.0 | 18 / 10 / 22 |

Seed-averaged per target (mean of the three seeds), then the median over targets and the gate fraction on the seed-mean:

| scheme | median of seed-mean `n_real` | IQR | seed-mean ≥ 10 | seed-mean ≥ 20 | median paired Δ vs pooled seed-mean | up / down |
|---|---:|---|---:|---:|---:|---|
| pooled | 6.0 | [2, 46] | 17 (34 %) | 17 | +0.0 | 0 / 0 |
| gem | 5.7 | [2, 48] | 17 (34 %) | 16 | +0.0 | 16 / 18 |
| gemdepth | 6.3 | [2, 52] | 18 (36 %) | 17 | +0.0 | 21 / 12 |

Paired log2((n+1)/(n_pooled+1)) per target, seed-means: gem: median +0.000, mean -0.063; gemdepth: median +0.000, mean +0.057

Uncorrected-variance sensitivity (`n_real_uc`), median over targets per scheme (seeds 0/1/2): pooled: 4/4/4; gem: 4/4/4; gemdepth: 4/4/4

## Diagnostics per scheme (median over 50 targets × 3 seeds unless stated)

| scheme | ctrl pool UMI CV | ctrl UMI mean | depth ratio median [min, max] | outside [0.95, 1.05] / 150 | zero-fraction gap (pert − ctrl) median | gemgroups in ctrl pool | gemgroup TV distance target vs ctrl | n_ctrl |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| pooled | 0.440 | 11719 | 1.0135 [0.886, 1.128] | 23 | -0.0007 | 267 | 0.390 | 18400–18400 |
| gem | 0.441 | 11742 | 1.0102 [0.893, 1.117] | 21 | -0.0012 | 167 | 0.005 | 18400–18400 |
| gemdepth | 0.378 | 11656 | 0.9903 [0.927, 1.037] | 3 | +0.0013 | 167 | 0.005 | 18400–18400 |

Target cells' own gemgroups: median 167 per target (min 131, max 207); pooled-scheme per-target depth ratio IQR [0.996, 1.028].
Spearman of per-target log2 change (gem+depth vs pooled) with the pooled depth ratio: -0.281; with the pooled zero-fraction gap: +0.251; with pooled n_real: +0.245.

## DE-set overlap across schemes, seed 0 — five best-covered targets by cells (+ DCP2, the largest set)

| target | cells | pooled | gem | gem+depth | ∩ all three | J(pooled, gem) | J(pooled, gem+depth) | new in gem+depth (not pooled) | their pooled p-rank median [max] | their pooled q median | lost from pooled | ρ(−log p) pooled vs gem+depth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| NABP2 | 400 | 9 | 7 | 7 | 6 | 0.60 | 0.60 | 1 | 11 [11] | 0.059 | 3 | +0.829 |
| MYRF | 400 | 7 | 2 | 3 | 2 | 0.29 | 0.43 | 0 | nan [nan] | nan | 4 | +0.834 |
| TRIM35 | 382 | 6 | 6 | 5 | 5 | 1.00 | 0.83 | 0 | nan [nan] | nan | 1 | +0.899 |
| ZNF282 | 362 | 3 | 3 | 4 | 3 | 1.00 | 0.75 | 1 | 6 [6] | 0.054 | 0 | +0.914 |
| HMGXB3 | 342 | 85 | 95 | 101 | 77 | 0.82 | 0.74 | 22 | 107 [355] | 0.065 | 6 | +0.940 |
| DCP2 | 241 | 346 | 377 | 315 | 268 | 0.82 | 0.69 | 45 | 451 [1831] | 0.069 | 76 | +0.892 |

(A pooled p-rank of k means the gene was the k-th smallest pooled p-value among the 8,248 tested, target gene included; 'same genes at a looser bar' would read as ranks just beyond the pooled set size with pooled q just above 0.05.)

- NABP2 new in gem+depth: HIST1H1C
- ZNF282 new in gem+depth: PRKAR1A
- HMGXB3 new in gem+depth: SEC62, RPL34, HEXB, DCP2, HIST1H2BC, EHMT2, RPS18, GGCT, KLHDC10, MTDH, MYC, PHF19, TMEM258, LAMTOR1, CHD2, RPL19, RPL36, EMP3, SMARCA1, MT-CO1, MT-CO2, MT-ND5
- DCP2 new in gem+depth: STMN1, MARCKSL1, ARF1, NCKAP1, PTMA, ABCC5, LRCH3, RPL34, BDP1, RHOBTB3, FAM8A1, RPS18, ME1, GCC1, DCAF13, WASHC1, ANP32B, WDR34, RPL7A, HNRNPF, JMJD1C, EIF4G2, SELENOH, HINFP, LARP4 …
