# 276 — The conservation map: for every panel target, how conserved its response direction is across the cell lines and datasets on disk, replicate-validated, and whether conservation predicts transfer

**CONSERVATION REPRODUCES ACROSS HALVES: YES** (strict, four sources with independent cells, ρ **+0.273 [+0.158, +0.388]** on 295 targets; and across two DISJOINT SETS OF CELL LINES ρ **+0.435 [+0.340, +0.522]**, 50/50 random line splits positive, median +0.394). **PREDICTS TRANSFER UNDER LOLO: YES** (ρ **+0.387 [+0.326, +0.444]** with `cd4svd5`'s per-target `pds` pooled over the three cancer folds, **AUC 0.694** against 262's anti-aligned group; +0.239 [+0.02, +0.435] and AUC 0.634 on h1-test's 85 Flex-400 targets) — **but it does NOT reliably beat the two-donor agreement the pipeline already has** (paired Δρ +0.066 [+0.004, +0.133] / +0.052 [−0.011, +0.117] / −0.028 [−0.166, +0.106] on the three readings).

The map is built from **26 single-line genetic-perturbation sources spanning 19 cell lines/states**, all already on disk, in one frozen residual space; the statistic is a mean cross-line cosine with no free parameter. Its distribution is the scientific deliverable: **a panel target's response direction agrees across cell lines at a median cosine of +0.020** (mean +0.026, p10 −0.003, p90 +0.068, 16.0 % ≤ 0) against a target-shuffled null of −0.0003 ± 0.0013; **the same targets' responses agree at +0.327 within a line** (10 same-line pairs) — conservation across cell state is an order of magnitude below reproducibility within it; and **once the measurement noise of the 150–200-cell references is divided out, the shared fraction rises to a median cosine of +0.228** (22.5 % of targets still ≤ 0, n 227).

Four honest limits on the positive verdict, all measured here:

1. **The nineteen lines do not reliably beat the two the pipeline already has.** `C_lolo` over 24 sources reads ρ +0.387; 262's two-donor agreement on the same axis reads +0.307; the paired difference is +0.066 [+0.004, +0.133] — and on the two other readings of the same comparison it does not clear zero: **+0.052 [−0.011, +0.117]** on the wider per-pair gene axis (§4.6) and **−0.028 [−0.166, +0.106]** on h1-test. Buying more contexts (programme rule 7) buys **little or nothing** on this statistic; 1 of 3 readings clears, and it is the narrowest-axis one.
2. **About half of `C_lolo` is the incumbent's own components agreeing with each other.** Dropping every source the incumbent uses (the K562/HCT116/HEK293T donors, the three CD4 conditions, the two H1 banks) leaves a 16-source, 3-pair-per-target score at ρ **+0.198 [+0.127, +0.263]** — still clearly positive, so an entirely external signal exists, but it is half the size and much thinner (POST-HOC).
3. **Conservation is entangled with measurability.** The measured reliability of the non-held sources alone reads ρ +0.283 and the response magnitude +0.207. `C_lolo` survives both (paired Δρ +0.104 [+0.026, +0.183] over reliability, +0.181 [+0.101, +0.258] over magnitude; partial ρ +0.347 and +0.353), but the reliability-corrected score under LOLO falls to +0.290 [+0.208, +0.379], i.e. inside the two-donor baseline's interval.
4. **It is a diagnosis, not a lever.** `C_lolo` predicts the base·truth cosine at ρ +0.388 — essentially the same as its correlation with `pds` — so it says *where* the whitened donor consensus points away from the held line's truth (F4's mechanism), not how to move it there. 263 already closed the arbiter-driven action for every observable arbiter.

`work/conservation276.csv` is the per-target deliverable (300 rows × 26 columns: `C_raw`, its pair count, bootstrap band and shuffled null, `C_within`, `C_rel`, the two CMap scores, the two independent halves, per-fold `C_lolo` and `pds`, CCLE breadth and level, essentiality, shared-programme share).

**Date 2026-09-06, 04:35–05:45 UTC (07:35–08:45 Ankara).** Deterministic (seed 20260906). **Measurement only: no build, no `vcc`, no `cell-eval2`, no model, no submission, nothing under `submissions/` or `data/`, no production file touched, no download** — everything was already on disk. `uv run python` throughout. Peak RSS 2.02 GiB (the whitener eigenbasis load), every other stage under 1 GiB. `find data scripts src docs submissions experiments -newer work/DEFINITION-276.md` outside this directory returns nothing but sibling agents' own logs.

---

## 0. Known answers (all read before anything else)

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 262's per-fold `cd4svd5` `pds` means | 0.815384 / 0.801659 / 0.780978 / 0.640294 (R262 KA-1) | **0.8153840 / 0.8016585 / 0.7809779 / 0.6402939** | PASS |
| **KA-2** 262's per-fold ρ(base·truth cosine, `pds`) | +0.913 / +0.885 / +0.855; kolf2 +0.795 (R262 §4.1) | **+0.9131 / +0.8846 / +0.8554; +0.7948** | PASS |
| **KA-3** 262's `agree_donor_donor` pooled ρ = +0.381, recomputed in this residual space | ±0.02 | **+0.403** on the three donors' own 7,543-gene shared axis, **+0.396** on the per-pair wide axis (§4.6); **+0.306** on 276's 2,994-gene common axis | PASS on both wide axes (+0.022 / +0.015 from 262's own convention); the 2,994-gene axis costs 0.097 of ρ — see §2.3 |
| **KA-4** 271's per-target split-half reliability medians | +0.047 / +0.028 / +0.028 for K562 / HCT116 / HEK293T (F1, R267 §3) | **+0.04706 / +0.02756 / +0.02825** | PASS |
| **KA-5** 262's anti-aligned share of (target, fold) rows | 17.5 % (R262 §0.4) | **17.4 %** (148/850) | PASS |
| **KA-6** F4's anti-aligned share on h1-test | 12.9 % (R273 §3) | **12.9 %** (11/85) | PASS |

The statistic was frozen in `work/DEFINITION-276.md`, sha256 `4e90f8ebc006af1ce67b90060e00d795648e9f729fe49063039aab9a323bd40c`, at **2026-09-06 05:04 UTC**, after the inventory and before `build276.py` computed a single cosine. Everything outside that file is labelled POST-HOC below.

---

## 1. Inventory (item 1)

43 candidate banks on disk were scanned (`inv276.py` → `work/inventory276_sources.csv`, `work/inventory276.json`). The unit everywhere is the **compact-donor table** `scripts/build_donor_table.py` writes — `delta.npy` (log-pseudobulk effect vs the source's own control, `b − b_ctrl`, `b = log1p(5e4·P/ΣP)`), `b_ctrl.npy`, `genes.txt`, `perts.txt` — which is the object production itself uses; the two LINCS/CMap banks are 270's per-(target, line) MODZ z-vectors.

### 1.1 The 26 sources in the map (19 distinct cell lines / states)

| line | sources (dataset) | panel targets | genes on the challenge axis | perturbation | unit |
|---|---|---:|---:|---|---|
| K562 | `k562_gwps` (Replogle 2022 GWPS) | 269 | 7,679 | CRISPRi | cells, 267 gem groups |
| K562 | `perturbmatch` (K562 reanalysis) | 251 | 9,373 | CRISPRi | cells |
| K562 | `jiang_K562` (Jiang 2025) | 9 | 9,517 | CRISPR-KO | bulk log2FC |
| HCT116 | `orion_hct116` (X-Atlas/Orion) | 290 | 18,106 | CRISPRi | cells, batches |
| HEK293T | `orion_hek293t` (X-Atlas/Orion) | 294 | 18,106 | CRISPRi | cells, batches |
| HEK293 | `xucao_hek293` (Xu 2023) | 5 | 18,277 | CRISPRi | cells |
| KOLF2.1J | `kolf2` | 273 | 18,106 | CRISPRi | z-scores |
| iPSC (HipSci pool) | `feng_ipsc` (Feng 2025) | 182 | 18,533 | CRISPRi | cells, median 22/target |
| iPSC-neuron | `tian_ipsc_neuron` (Tian 2021) | 6 | 18,333 | CRISPRi | cells |
| H1 hESC | `gse295214_h1` | 64 | 18,533 | CRISPRi | cells |
| H1 hESC | `h1_2025` (VCC 2025) | 20 | 18,077 | CRISPRi | cells, 10x Flex |
| melanoma (patient) | `frangieh_melanoma` (Frangieh 2021) | 5 | 15,924 | CRISPR-KO | cells |
| TeloHAEC | `igvf_telohaec` | 65 | 15,074 | CRISPRi | DE table |
| cardiomyocyte | `igvf_cardiomyocyte` | 87 | 7,731 | CRISPRi | DE table |
| A549 | `myllia_a549`, `altair_a549`, `jiang_A549` | 14 / 20 / 9 | 18,533 / 18,533 / 9,517 | CRISPRi / CRISPRi / KO | cells (2 channels) / cells / log2FC |
| DLD-1 | `dld1_lowmoi` (two low-MOI replicates) | 67 | 18,533 | CRISPRi | lfc |
| Calu-3 | `sunshine_calu3` (Sunshine 2023) | 6 | 18,333 | CRISPR | cells |
| BxPC-3, HAP1, HT-29, MCF7 | `jiang_*` (Jiang 2025) | 9 each | 9,517 | CRISPR-KO | log2FC |
| CD4 T (rest / stim 8 h / stim 48 h) | `cd4_*` (GWCD4i) | 283 / 287 / 284 | 10,113 | CRISPRi | cells, pooled donors |

Excluded by the frozen rule: line aggregates (`jiang_6line`, `jiang_5line_noK562`, the three line-collapsed CMap tables), mouse (`hspc_crispri500`, absent anyway), CRISPRa (`norman_k562`, `tian_ipsc_neuron_crispra` — the opposite perturbation direction), and every bank carrying < 5 panel targets (`k562_essential`, `rpe1`, `nadig_jurkat`, `nadig_hepg2` all carry **0/300**; `gasperini_k562`, `tian2019_ipsc`, `wessels_thp1` carry 1). CMap's per-(target, line) signatures enter only as the declared secondary of §5.

### 1.2 Targets × independent lines

| number of independent lines | single-cell / DE sources | including CMap |
|---|---:|---:|
| ≥ 2 | **300** | 300 |
| ≥ 3 | 299 | 300 |
| ≥ 5 | 293 | 297 |
| ≥ 8 | 236 | 265 |
| ≥ 12 | 12 | 149 |
| median per target | **8** (min 2, max 17) | **11** (min 4, max 38) |

Median **10 sources** and **33 cross-line pairs** per target (min 1, max 182). The distribution is tight — 88 targets sit at exactly 8 lines — because the four wide sources (Orion ×2, K562-GWPS, KOLF2.1J, CD4 ×3) cover almost the whole panel and everything else is thin.

### 1.3 Leak surface

* **The challenge lines** are A = JURKAT, B = HELA, C = CAL33 (exp 041's CCLE match, AGENTS.md §2.1). No source in the map is one of them: `nadig_jurkat` is Jurkat but carries **0/300** panel targets, so it never enters; HeLa and CAL-33 have no genetic screen on disk; the CMap curated genetic table contains none of the three for panel targets. **The map cannot leak the challenge contexts.**
* **Our own PLOLO folds are sources here.** K562, HCT116, HEK293T and KOLF2.1J each contribute one or more sources, and their rows are exactly the truths the folds are scored against. Every prediction test in §4 therefore reports the **LOLO** version, in which every source from the held line is removed (K562 → `k562_gwps`, `perturbmatch`, `jiang_K562`; HCT116 → `orion_hct116`; HEK293T → `orion_hek293t` **and** `xucao_hek293`, HEK293 being HEK293T's parental line; KOLF2.1J → `kolf2` **and** `feng_ipsc`, whose HipSci pool contains `kolf_2`; H1 → `gse295214_h1`, `h1_2025`). `C_raw` is reported only to show the size of the leak (+0.488 vs +0.387).
* **Sources that are inside the prediction.** `k562_gwps`, `perturbmatch`, `orion_hct116`, `orion_hek293t` (the whitened base), the three CD4 conditions and the two H1 banks (the residual) are all components of `cd4svd5`. §4.4 measures the score with all nine removed.
* CMap's genetic signatures for panel targets cover 15 shRNA lines and 20 CRISPR-KO lines, **none of them a fold line** — so the CMap score is inherently leave-the-held-line-out (its LOLO and raw values are identical to every digit).

---

## 2. The statistic (item 2), frozen before computation

### 2.1 Definition

Common gene axis **G**: the intersection over the 26 sources of their gene symbols, intersected with the challenge's 18,533 → **2,994 genes**, in the challenge's own order. For source *s*, its rows for the panel targets it carries, restricted to G, minus the mean over exactly those rows (each source centred on its own across-target mean; no whitening, no scaling — the statistic is a cosine).

* **`C_raw(t)`** = mean over all pairs of sources from **different cell lines** that both carry *t* of `cos(R_s[t], R_s'[t])`.
* **`C_within(t)`** = the same over **same-line** pairs (10 exist: three K562 sources, three A549 sources, three CD4 states) — reproducibility, not conservation.
* **`C_rel(t)`** = mean over pairs where both sources have a measured per-target split-half reliability of `cos / sqrt(r_s · r_s')`, `r` Spearman–Brown-corrected to full depth, pairs dropped where either `r ≤ 0.02`.
* Intervals: paired bootstrap over that target's own pairs, B = 2,000. Null: target-shuffled pairings that keep the pair structure and gene sets exactly, 200 draws.

### 2.2 The distribution over the 300 targets

| statistic | mean | median | p10 | p90 | share ≤ 0 | n |
|---|---:|---:|---:|---:|---:|---:|
| **`C_raw`** cross-line (315 pairs, median 33/target) | **+0.0263** | **+0.0199** | −0.0028 | +0.0676 | **16.0 %** | 300 |
| `C_within` same-line (10 pairs) | +0.3274 | +0.3213 | — | — | — | 296 |
| target-shuffled null (200 draws) | −0.00034 | −0.00023 | — | — | — | 300 |
| **`C_rel`** reliability-corrected (4 sources, 6 pairs) | +0.2536 | **+0.2280** | −0.2088 | +0.6923 | 22.5 % | 227 |

The null's spread over targets is 0.00133, so the median target's `C_raw` of +0.020 is **15 null SDs** above chance; the per-target 95 % null band averages [−0.0128, +0.0122], so a single target's score is resolved only when it is well outside that band (which is why the pair-bootstrap band is carried per target in the CSV). 459 of the reliability-corrected pairs were dropped at the `r ≤ 0.02` floor.

Two readings the numbers force:

* **Conservation across cell state is ~16× smaller than reproducibility within a line** (+0.020 vs +0.327). The same-line pairs are not a clean replicate control (`k562_gwps` and `perturbmatch` are two analyses of the same cells, the three CD4 states share donors), so +0.327 is an upper bound on "same line, different experiment" — but the gap is far too large to be an artefact of that.
* **Most of the cross-line gap is measurement, not biology.** Dividing by the measured √(r·r′) lifts the median from +0.023 (on the four-source sub-map) to **+0.228** — i.e. of the direction a perfectly-measured 150-cell reference would show, roughly **a quarter is shared across cell lines**, and roughly a fifth of targets share nothing (22.5 % ≤ 0).

### 2.3 What the 2,994-gene common axis costs (a diagnostic, in the frozen plan)

On the three production donors, the pair cosine computed on the common axis tracks the cosine on the donors' own 7,543-gene shared axis at Spearman **+0.852 / +0.875 / +0.960**, with almost the same medians (+0.033 vs +0.032; +0.033 vs +0.031; +0.047 vs +0.043). But the *prediction* is measurably weaker on the narrow axis: 262's two-donor agreement reads ρ +0.403 with `pds` on the donors' 7,543-gene axis and +0.306 on the 2,994-gene common axis. **Every comparison in §4 is therefore made on the same 2,994-gene axis for predictor and baseline alike**, so the comparisons are fair — but the absolute ρ values are all larger on a richer axis, and §4.6 shows the narrow axis slightly *flatters* the map relative to the two-donor baseline rather than penalising it.

### 2.4 Top and bottom of the map (by `C_raw`)

**Top 30** (most conserved): TAF4, SIN3B, SYVN1, HMGN3, HMGXB3, MARCKSL1, DCP2, PTMS, EPHB2, ANP32B, ACLY, CSNK2A1, TBCK, GIGYF2, GMEB2, ZNF507, FOXN2, PQBP1, MED15, SND1, CARM1, ZNF714, UPF3B, PNP, PCGF1, SFXN4, FXYD5, SCAF8, HMGN4, S100A11.

**Bottom 30** (most state-dependent, or least measurable): RAB11FIP5, RBM4B, ZNF786, TXNDC16, TDRD7, SAP30L, PSMB9, IFI27L2, SNRK, SPATA20, NICN1, RFK, PLCG1, CHKB, ZKSCAN3, SNN, ZNF324, DVL3, NT5DC1, DTNBP1, ZNF251, ZNF768, SH2D3A, ZNF174, MAPK7, PGLS, ZNF627, CYTH2, NFE2L1, ABCD1.

Decile cuts: `C_raw` ≤ −0.0028 (bottom) and ≥ +0.0676 (top).

---

## 3. Replicate validation of the score itself (item 3)

### 3.1 Strict split-half — the primary test

Four sources have independent halves on disk: K562-GWPS, Orion HCT116 and Orion HEK293T (271's gemgroup/batch-parity half rows, built from 267's integer pseudobulk sums) and DLD-1 (its two low-MOI replicates, `dld1_low1` / `dld1_low2`). The map was rebuilt **from half A of those four and nothing else**, and again from half B; the two share no cells. Six cross-line pairs, median 3 per target, on a 7,543-gene axis.

| | value |
|---|---|
| Spearman(`C_A`, `C_B`) over the 295 targets where both are defined | **+0.273 [+0.158, +0.388]** |
| median `C_A` / `C_B` | +0.0140 / +0.0120 |
| declared secondary: the full 26-source map with the halved sources' halves swapped and the other 22 sources shared | +0.810 [+0.742, +0.866] — an **upper bound**, the shared sources guarantee most of it |
| DLD-1 replicate-pair reliability (its own known answer) | median +0.0863 on 67 targets |

**Reading.** A conservation score built from three noisy pairs at 70–95 cells per half reproduces at ρ +0.273 — positive, its interval well clear of zero, but weak, exactly as F1's reliability law predicts for a statistic assembled from 150-cell responses. It is a *lower* bound on the production score, which averages 33 pairs at full depth.

### 3.2 Disjoint cell lines — the second, independent test (POST-HOC in construction, same statistic)

A stronger question for a conservation statistic: does conservation measured on one set of cell lines predict conservation measured on a **disjoint** set? The 19 lines were split so that no line appears on both sides.

| split | ρ(score on group 1, score on group 2) | n |
|---|---|---:|
| declared split (G1 = K562, HEK293T, HEK293, KOLF2.1J, cardiomyocyte, A549, HAP1, MCF7, melanoma, Calu-3; G2 = the rest) | **+0.435 [+0.340, +0.522]** | 289 |
| 50 random balanced line splits | median **+0.394**, p10–p90 [+0.292, +0.466], **50/50 positive** | ≥ 30 each |

**This is the finding that makes the map a real object**: a target's conservation is a property of the target, not of which lines happened to measure it.

---

## 4. Does conservation predict transfer? (item 4)

Outcomes are pre-existing and untouched: `cd4svd5`'s per-target `pds` and `cov|base_cos_truth` on the four PLOLO folds (262), and `INCUMBENT|pds_eff` on h1-test (273). Pooling is by within-fold rank (262's own convention); intervals are paired target bootstraps, B = 2,000.

### 4.1 The prediction table

| predictor | K562 fold | HCT116 fold | HEK293T fold | kolf2 fold | pooled (3 cancer folds) | AUC rest > anti-aligned |
|---|---:|---:|---:|---:|---:|---:|
| `C_raw`, all 26 sources (**leaky**, the held line is a source) | +0.508 | +0.497 | +0.469 | +0.366 | **+0.488** [+0.432, +0.540] | 0.756 |
| **`C_lolo`, every source from the held line removed** | +0.411 | +0.390 | +0.365 | +0.278 | **+0.387** [+0.326, +0.444] | **0.694** |
| `C_rel` reliability-corrected (**leaky**) | +0.419 | +0.476 | +0.497 | +0.165 | **+0.461** [+0.398, +0.522] | 0.833 |
| `C_rel_lolo` reliability-corrected, held line out | +0.210 | +0.265 | +0.373 | +0.165 | **+0.290** [+0.208, +0.379] | 0.659 |
| baseline: 262's donor–donor agreement (the two-line version) | +0.207 | +0.322 | +0.373 | +0.179 | **+0.307** [+0.239, +0.372] | 0.652 |
| baseline: CD4 condition count | +0.089 | +0.044 | +0.085 | −0.037 | +0.056 [−0.010, +0.121] | 0.530 |
| baseline: target expression in the challenge control | +0.052 | +0.046 | +0.053 | +0.060 | +0.051 [−0.015, +0.118] | 0.517 |
| control: measured reliability of the non-held sources | +0.200 | +0.238 | +0.390 | +0.210 | +0.283 [+0.215, +0.346] | 0.574 |
| control: response magnitude of the non-held sources | +0.158 | +0.236 | +0.219 | +0.203 | +0.207 [+0.144, +0.269] | 0.549 |
| control: number of pairs (coverage) | +0.235 | +0.184 | +0.171 | +0.322 | +0.197 [+0.130, +0.262] | 0.575 |

AUC > 0.5 means the covariate separates the anti-aligned group in the expected direction (262's own AUC convention is inverted; 0.694 here corresponds to 262's 0.306).

### 4.2 Is it more than "the response is measurable"?

| comparison (pooled, 3 cancer folds, 850 (target, fold) rows) | paired Δρ | partial ρ(`C_lolo`, `pds` \| the other) |
|---|---|---|
| `C_lolo` − donor–donor agreement | **+0.066 [+0.004, +0.133]** | +0.257 [+0.191, +0.326] |
| `C_lolo` − measured reliability | +0.104 [+0.026, +0.183] | +0.347 [+0.284, +0.405] |
| `C_lolo` − response magnitude | +0.181 [+0.101, +0.258] | +0.353 [+0.291, +0.413] |
| `C_lolo` − number of pairs | +0.191 [+0.117, +0.273] | +0.356 [+0.294, +0.416] |
| `C_lolo` − `C_rel_lolo` | +0.103 [+0.014, +0.195] | +0.320 [+0.236, +0.402] |

`C_lolo` beats every control, and every partial correlation stays above +0.25. Against the two-donor baseline the margin is real but small and its interval nearly touches zero.

### 4.3 The same, read as groups (the interpretable form)

| `C_lolo` band | n (target, fold) | mean `pds` | share anti-aligned | share of the total `pds` deficit | median `C_lolo` |
|---|---:|---:|---:|---:|---:|
| bottom decile | 85 | **0.647** | **36.5 %** | 17.5 % | −0.0059 |
| 10–30 % | 170 | 0.706 | 27.1 % | 29.2 % | +0.0026 |
| 30–50 % | 170 | 0.794 | 19.4 % | 20.5 % | +0.0110 |
| 50–70 % | 170 | 0.822 | 12.4 % | 17.7 % | +0.0234 |
| 70–90 % | 170 | 0.901 | 7.6 % | 9.9 % | +0.0432 |
| top decile | 85 | **0.895** | **4.7 %** | 5.2 % | +0.0786 |
| all | 850 | 0.799 | 17.4 % | 100 % | +0.0171 |

The bottom 30 % of rows by conservation carry **46.8 %** of the deficit; the anti-aligned share runs 36.5 % → 4.7 % monotonically across the six bands. For comparison, 263 screened ≈ 130 build-time arbiter cells and its best bin reached mean `pds` **0.685**; the conservation bottom decile reaches **0.647** with one pre-frozen statistic and no screen. The oracle group that knows the held line still sits at 0.295 (R263 §3), so the gap to the ceiling is unchanged in kind.

### 4.4 h1-test (85 Flex-400 targets, the reliability-resolving reference)

Only 5 of h1-test's 85 targets are on the 2026 panel, so the map was rebuilt on **its** targets with the identical construction (same 26 sources, same 2,994-gene axis, same centring rule applied to the rows each source contributes to that map). Coverage there is good: K562-GWPS 82/85, Orion HCT116 84, Orion HEK293T 84, DLD-1 28, IGVF cardiomyocyte 26, TeloHAEC 15.

| predictor | ρ with `INCUMBENT` per-target `pds` | AUC rest > anti-aligned (11/85) |
|---|---|---|
| `C_raw` (H1 sources included — leaky) | +0.330 [+0.113, +0.514] | — |
| **`C_lolo`** (both H1 banks removed) | **+0.239 [+0.020, +0.435]** | **0.634** |
| baseline: donor–donor agreement | +0.267 [+0.067, +0.447] | 0.555 |
| baseline: target expression | +0.103 [−0.108, +0.320] | — |
| control: number of pairs | +0.198 [−0.035, +0.398] | — |
| paired Δρ (`C_lolo` − donor–donor) | **−0.028 [−0.166, +0.106]** | — |
| `C_lolo` terciles → mean `pds` | bottom 0.804 · middle 0.900 · top 0.924 | — |

**On the challenge's own chemistry and depth the map predicts transfer, monotonically — and adds nothing over the two donor lines already in the pipeline.**

### 4.5 The mechanism, and what happens when the incumbent's own sources are removed (POST-HOC)

| conservation score | sources | median pairs/target | ρ with `pds` | ρ with base·truth cosine |
|---|---:|---:|---|---|
| `C_lolo` (held line out) | 24 | 25 | +0.387 [+0.326, +0.444] | +0.388 [+0.328, +0.447] |
| `C_out` (held line out **and** every source the incumbent uses) | 16 | 3 | **+0.198 [+0.127, +0.263]** | +0.225 [+0.155, +0.297] |
| held line out, CD4 pairs only removed | 21 | 10 | +0.370 [+0.309, +0.430] | +0.384 [+0.323, +0.443] |

`C_lolo` predicts the base·truth cosine as well as it predicts `pds` — which is the whole mechanism, since `pds` is a threshold in that cosine at ≈ 0.02 (F4). And with every source the prediction itself uses removed, a 16-source, 3-pair score still reads +0.198: the signal is not an artefact of the prediction's components agreeing with each other, but roughly half of `C_lolo`'s size is. The CD4 residual is not the driver.
`C_out` bands: bottom decile mean `pds` 0.710 / 29.2 % anti-aligned, top decile 0.914 / 4.2 %.

### 4.6 The gene axis, and a split-panel confirmation (POST-HOC, `wide276.py`)

The frozen 2,994-gene axis is the conservative choice (§2.3). Rebuilding every pair on **its own two sources' shared genes** (median 9,517 per pair, min 4,185, max 18,533; each source centred once on its own full official-gene set) gives a score that ranks targets almost identically (ρ +0.811 with the frozen `C_raw`) but on a slightly lower and less negative scale (mean +0.0201, median +0.0143, only 8.3 % ≤ 0 — the narrow axis over-represents targets with a sign flip). Every prediction rises together:

| predictor | wide per-pair axis | frozen 2,994-gene axis |
|---|---|---|
| `C_lolo` (held line out) | **+0.457 [+0.401, +0.510]** | +0.387 [+0.326, +0.444] |
| baseline: donor–donor agreement | **+0.396 [+0.331, +0.456]** | +0.307 [+0.239, +0.372] |
| `C_out` (held line out and every incumbent source) | +0.253 [+0.184, +0.316] | +0.198 [+0.127, +0.263] |
| paired Δρ (`C_lolo` − donor–donor) | **+0.052 [−0.011, +0.117]** | +0.066 [+0.004, +0.133] |

The donor–donor baseline on the wide axis reads +0.396 against 262's own +0.381 — a third confirmation of KA-3. **On the wide axis the map's advantage over the two donor lines is not established**; the conservation signal itself is stronger everywhere.

Split-panel confirmation — the 300 targets cut in two by `sha256(name) mod 2`, so the prediction is read on two disjoint halves of the panel rather than once:

| half | rows | `C_lolo` (frozen axis) | `C_lolo` (wide axis) | donor–donor (wide) |
|---|---:|---|---|---|
| 0 | 407 | +0.392 [+0.306, +0.474] | +0.454 [+0.374, +0.525] | +0.367 [+0.273, +0.455] |
| 1 | 443 | +0.380 [+0.295, +0.453] | +0.458 [+0.374, +0.530] | +0.417 [+0.325, +0.494] |

The prediction is the same in both halves to within 0.012 — it is not carried by a subset of targets.

---

## 5. What conservation is made of (item 5, descriptive)

| covariate | bottom decile (median) | top decile (median) | whole panel | ρ with `C_raw` |
|---|---:|---:|---:|---:|
| CCLE expression breadth (share of 1,600 lines with TPM > 1) | 0.997 | 1.000 | 0.999 | +0.198 |
| CCLE median expression (TPM) | 3.60 | 4.61 | 4.32 | +0.200 |
| CCLE expression SD across lines | 0.823 | 0.773 | 0.778 | −0.030 |
| K562-essential membership | **0/30** | **0/30** | **0/300** | undefined (constant) |
| share of the response energy in the production whitener's **top-decile** directions (the shared programme) | 0.255 | 0.277 | 0.261 | **+0.234** |
| `C_within` (same-line reproducibility) | 0.278 | 0.382 | 0.321 | **+0.404** |
| number of cross-line pairs (coverage) | 32 | 41 | 33 | +0.245 |

* **DepMap co-essentiality could not be computed**: the only DepMap material on disk (`data/external/depmap26q1/extracted/processed/`) is a single-line (HCT15) guide-level LFC/FDR matrix, not the cross-line gene-effect matrix a co-essentiality module needs. Nothing was downloaded. The nearest on-disk proxy is the essentiality flag, and it is **constant at zero across the whole panel** — the panel is non-essential by design (AGENTS.md §2.1), which also removes essentiality as a confound.
* **Expression breadth is uninformative by construction**: every panel target is expressed in ≥ 99.7 % of CCLE lines, because "sufficiently high basal expression of the target gene" was an inclusion criterion (AGENTS.md §2.1). The residual signal is in *level*, not breadth (ρ +0.200), and it is weaker than conservation's correlation with `pds`.
* **The conserved targets are, weakly, the ones the shared programme captures** (ρ +0.234): the top decile puts 27.7 % of its response energy in the whitener's top-decile directions against 25.5 % in the bottom decile. The effect is real but small — conservation is not simply "this target drives the common programme".
* **The strongest correlate of conservation is same-line reproducibility** (ρ +0.404), which is the honest statement of the entanglement in §4.2: a target whose response repeats within a line is also the one whose response repeats across lines.
* **A gene-family signal, POST-HOC.** KRAB/C2H2 zinc-finger targets (35/300 by symbol: `ZNF*`, `ZKSCAN*`, `ZSCAN*`) have a median `C_raw` of **+0.0044** against +0.0223 for the rest (Mann–Whitney p 1.1e-05), and make up 23.3 % of the bottom decile against 6.7 % of the top. Their measured reliability is also lower (+0.041 vs +0.079, p 6.0e-03), and after the reliability correction the conservation gap shrinks to +0.130 vs +0.234 at p 0.096 — so **the zinc-finger deficit is mostly a weak, poorly-measured response, and a residual state-dependence that is not established at n = 28**. The top decile is dominated by chromatin/basal-transcription machinery and RNA/metabolic housekeeping (TAF4, SIN3B, CARM1, MED15, PCGF1, PQBP1, SND1, ACLY, CSNK2A1, DCP2, UPF3B).

### 5.1 Two external cross-checks (declared secondaries)

**LINCS/CMap, the same statistic on a different modality.** Cross-line cosine on the 978 measured landmark genes, per reagent, with the same centring: median **+0.163** across 15 shRNA-consensus lines (113 targets) and **+0.091** across 20 CRISPR-KO lines (82 targets) — far above the single-cell map's +0.020, as L1000 consensus signatures should be. But it carries **no information about transfer**: pooled ρ with `pds` **+0.013 [−0.110, +0.126]** (shRNA) and **+0.125 [−0.004, +0.259]** (CRISPR-KO), and it barely correlates with the single-cell map itself (ρ −0.032 and +0.175 with `C_raw`). Because none of our fold lines has CMap genetic signatures for panel targets, this score is inherently leave-the-held-line-out — its LOLO and raw values are identical to every digit. **This confirms 270's negative in a second construction** (270 used top-100 sign agreement; this is a cosine on the residual space) and extends the closure: *cross-line agreement of L1000 genetic signatures does not measure the conservation that matters for this transfer*, on 82–113 targets and two reagents.

**GPSAdb bulk cross-study sign stability (270's frozen `SS_cross`).** On the 40 targets where 270 could compute it, ρ with `C_raw` = **+0.350 [−0.011, +0.631]** — the same magnitude as 270's own correlation with `pds` (+0.363), and the two resources broadly agree on which targets are conserved (270's high group AKAP1/ANP32B/BCL6/CARM1/PAPOLA/SCAF8/EPHB2 contains four of this map's top 30; its anti-aligned four contain ZKSCAN3, in this map's bottom 15). At n = 40 the interval touches zero, so this is a convergence worth noting, not a validation.

---

## 6. Reading, and the stage-2 recommendation

**What is established.** A per-target conservation statistic over 19 cell lines exists, is reproducible in two independent senses (disjoint cells ρ +0.273; disjoint *lines* ρ +0.435 with 50/50 random splits positive), and predicts per-target transfer under leave-the-held-line-out on both instruments (PLOLO ρ +0.387, AUC 0.694; Flex-400 ρ +0.239, AUC 0.634), beating reliability, response magnitude, coverage, CD4 DE count and target expression, and separating the anti-aligned group 36.5 % → 4.7 % across its deciles. The scientific numbers that outlive the competition question are the distribution itself: **cross-line direction agreement of +0.020 median (16 % ≤ 0) against +0.327 within a line, rising to +0.228 (22.5 % ≤ 0) once the references' measured unreliability is divided out.**

**What is not.** The 17 extra lines beat the two the pipeline already has by ρ +0.066 [+0.004, +0.133] on PLOLO's frozen axis, +0.052 [−0.011, +0.117] on PLOLO's wide axis and −0.028 [−0.166, +0.106] on the Flex-400 reference — **one of three readings clears zero, and it is the narrowest-axis one**. Under programme rule 7 that is a measured answer on the "independent contexts" scaling axis: **on this panel, at this reference depth, the marginal information in a 3rd–19th cell line is close to nothing**, because each new line's response is measured at reliability 0.03–0.09 and the two deep donors already span most of what the map can see. Half of the LOLO signal is the prediction's own components agreeing (`C_out` +0.198 vs `C_lolo` +0.387; +0.253 vs +0.457 on the wide axis).

**Stage 2 — the recommendation.** Not a target-class *dose* on this evidence, and not another arbiter.

1. **Competition track: do not build a dose from this map yet, and do not spend a slot on it.** The mechanism forbids it: `C_lolo` predicts the base·truth cosine (ρ +0.388), i.e. it says *where the whitened donor consensus points away from the held line*, and 262/263 established that the residual cannot repair that group (the anti-aligned deficit is the base, and every observable arbiter that tried to act on it lost `pds`, a wrong pick costing −0.02). A monotone per-target dose `g(C_lolo)` is the only action shape this statistic supports; it is amplitude-changing, so it faces the +0.0089 gate, and nothing here estimates its `pds` gain. **If it is ever priced, the correct prior is that the map buys between −0.03 and +0.07 of ρ over a covariate the pipeline already computes** — which is not a lever.
2. **The one competition-facing thing worth measuring next** is the opposite framing: because the bottom 30 % of rows carry 46.8 % of the deficit and are identifiable *before* the held line is seen, the question is whether a **different rendering** (not a different donor, not a different dose) on those rows can move them across the 0.02 knee — for example the un-deflated CD4 residual, which R273/F4 measured as +30 % `cos_agg` against +24 % deflated on the reliable reference. That is a one-variable arm with a mechanism, and this map supplies the row selector. It needs its own experiment and its own pricing; nothing here prices it.
3. **Scientific track: this map is the deliverable, and it should be extended on the depth axis, not the context axis.** The measured answer above says more lines at 150 cells add little; the reliability-corrected map (+0.228 median) says the ceiling is far higher than the raw map shows. The stage-2 experiment that would move the science is therefore **more cells per condition in a few lines**, not more lines — precisely the trade programme rule 7 asks to be measured, now measured on one of its three axes. The immediate cheap extensions: (a) recompute `C_rel` with reliabilities for the wide sources that have replicate structure and no half-bank yet (KOLF2.1J, the CD4 conditions, Myllia's two channels, Feng's donor lines) so the corrected map covers more than 227 targets; (b) preregister the ZNF-family state-dependence claim and test it on held targets — the raw effect is large (p 1.1e-05) and the reliability-corrected residual is exactly the confirmable part; (c) hand `C_lolo` to 269's router as its feature, where it is the first covariate that beats the labels' own noise floor.
4. **Closure at this construction.** Cross-line agreement of **L1000 genetic signatures** is closed as a predictor of per-target transferability in this pipeline — 270 closed the sign-agreement construction on 82–113 targets, and this experiment closes the residual-cosine construction on the same data (ρ +0.013 / +0.125, and ρ ≤ +0.175 with the single-cell map). The single-cell/DE cross-line map is **not** closed: it reproduces, it predicts, and it is the object stage 2 should deepen.

---

## 7. Files

* `RESULTS.md` (this file); `work/DEFINITION-276.md` (frozen statistic, sha `4e90f8eb…`) and `work/DEFINITION-276.sha256`.
* `inv276.py` → `work/inventory276_sources.csv`, `work/inventory276_targets.csv`, `work/inventory276.json`.
* `lib276.py` (the source registry, the common axis, the residual loader, the pair rules, the leave-the-held-line-out map).
* `build276.py` → `work/rows276.npz` (26 × 300 × 2,994 residual rows), `work/cons276.npz` (per-pair cosines, `C_raw`, `C_within`, pair bootstrap band, shuffled null), `work/build276.json`.
* `rel276.py` → `work/rel276.npz`, `work/rel276.json` (the four half-banks, `C_rel`, the strict split-half and the upper bound).
* `predict276.py` → `work/predict276.npz`, `work/predict276.json` (per-fold outcomes, the prediction table, the paired and partial comparisons, the h1-test map).
* `analyse276.py` → `work/analyse276.json`, **`work/conservation276.csv`** (the per-target deliverable), the disjoint-line splits, the CMap secondary, the decile descriptions.
* `tables276.py` → `work/tables276.md`; `extra276.py` → `work/extra276.md`; `wide276.py` → `work/wide276.md`, `work/wide276.npz` (the POST-HOC robustness reads: the incumbent-source removal, the wide per-pair gene axis, the split-panel confirmation).
