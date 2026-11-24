# 280 — The biological donor as a replicate unit: **a 10x Flex perturbation response in primary human CD4 T cells reproduces across four donors at a median cosine of +0.054, not +0.24 — the Flex-400 regime is a property of the H1 reference, not of Flex chemistry or of cell count**

**The two lines the mandate asks for.**

> **DONOR IS A REAL UNIT: YES, MARGINALLY — and it is not the thing that limits the measurement.**
> Preregistered primary (Stim48hr, the one run-clean condition): `Δ_donor = median cos_A − median cos_C` =
> **+0.0033 [+0.0004, +0.0072]**, n = 241, paired target bootstrap — the interval excludes 0, so the rule passes,
> but the effect is **1/30 of the cell-sampling term** (`D − A` = +0.0988 [+0.0918, +0.1058]), it does **not**
> clear zero in the other two conditions (Rest +0.0019 [−0.0015, +0.0065]; Stim8hr +0.0030 [−0.0002, +0.0057]),
> and at the matched ~350-cell scale the donor is **indistinguishable from the guide construct**
> (`GUIDE − DONOR` = +0.0007 [−0.0528, +0.0235], n = 97).
> **across-donor median +0.0539 [+0.0427, +0.0620]** (Stim48hr, 297-target panel, 3 disjoint pairings,
> ~250 cells/side, 28.7 % of targets ≤ 0); **within-donor at matched cells (the guide split, donors identical on
> both sides, 347 cells/side) +0.0484 [−0.0003, +0.0666]**.
>
> **FOUR-DONOR MEAN BEATS ONE AT MATCHED CELLS: YES — by +0.0067 cosine.**
> With 4 donors the largest honest leave-donors-out training set is **three**, so the measured quantity is the
> three-donor mean, and the preregistration says so. Stim48hr, matched total UMI budget:
> **Δcos(k=3 − k=1) = +0.0067 [+0.0040, +0.0108]**, n = 253 (Rest +0.0071 [+0.0045, +0.0105];
> Stim8hr +0.0055 [+0.0038, +0.0120]) — all three conditions clear zero. At **full** cells the same contrast is
> +0.0129 [+0.0098, +0.0176], so **52 % of the value of three donors is replicate averaging and 48 % is simply
> more cells**.
>
> **And the headline is estimator-independent.** The study's own DESeq2 fits computed *inside* each donor pair
> (`GWCD4i.DE_stats.by_donors.h5mu`, read by HTTP range request, nothing downloaded) give
> **+0.109 / +0.092 / +0.110** — on the 56–68 panel targets that pass DE eligibility inside a two-donor half.
> On *exactly those* targets my statistic reads **+0.154 / +0.141 / +0.139**. **What moves this number is
> which targets you keep, not which estimator you use: DE-eligibility selection is worth +0.10 of cosine, and
> even the selected best fifth of the panel (+0.14) stays below H1's +0.243 over all 85 of its targets.**

**Date** 2026-09-06, 06:48–07:55 UTC (09:48–10:55 Ankara). **Preregistration** `PREREGISTRATION.md`,
sha256 `7273ef3c48ad2ff6676c8f6fd55509b90ad254ea3b13c28bc2d6e1e819353c5b`, hashed at **06:57:35Z**
(`work/prereg280.sha256`) — the file's own header says "07:05Z", which was my estimate written a few minutes
before the hash; the hash timestamp is the binding one, and **no expression value had been read at that point**
(only the S3 key listing, the readme, the sample/QC tables, and the pseudobulk `.obs`/`.var`/`X.indptr`).
**No model, no build, no submission, no `vcc`, no `cell-eval2`; nothing written under `submissions/` or `data/`;
no production file touched.** `uv run python`, `PYTHONDONTWRITEBYTECODE=1` throughout. Peak RSS **5.42 GB**
(`/usr/bin/time -l` on `q1_280.py`); extraction 1.66 GiB. Peak disk in `work/` **46.0 GB**, cap 60 GB.

---

## 0. Known-answer rows

| row | required | measured | verdict |
|---|---|---|---|
| **KA-0** (extraction integrity) | every extracted CSR row's gene sum equals the pseudobulk `obs.total_counts` | **17,495 / 17,495 exact, max\|Δ\| = 0** | **PASS** |
| **KA-1** | pseudobulk `n_cells` summed over a panel target's guides across the 4 donors == `DE_stats.n_cells_target` on the local `GWCD4i.DE_stats.h5ad` (`c355f535…`) | **854 / 854 (target, condition) rows exact, max\|Δ\| = 0** — with the **all-guides** filter | **PASS** (the preregistration guessed `keep_for_DE`; that filter reads 455/854, max\|Δ\| 409. The DE table's `n_cells_target` counts *every* guide's cells) |
| **KA-1b** | the same for `DE_stats.n_guides` | **854 / 854 exact, max\|Δ\| = 0** — with the **`keep_for_DE`** filter (all-guides 795/854) | **PASS**. Together KA-1+KA-1b pin the bookkeeping exactly: cells count all guides, `n_guides` counts the DE-eligible ones |
| **KA-2** | pooled four-donor residual row vs the on-disk `cd4_compact` bank. Preregistered bar: median per-target cosine ≥ 0.50 on ≥ 400-cell targets **and** > 0 for ≥ 95 % | Rest **+0.351**, Stim8hr **+0.382**, Stim48hr **+0.385** (n = 269/270/268 deep targets, 9,828 shared genes); **> 0 for 100.0 %**; q10 +0.27 / +0.29 / +0.29 | **HALF-PASS — see §0.1.** The 95 % condition passes at 100 %; the median bar is **missed on all genes and cleared on the well-expressed half** |
| **KA-3** | join to `conservation276.csv` reproduces `C_raw` | 297 / 297 CD4-present panel targets carry a row, `C_raw` finite on all 297, max\|Δ\| = 0 | **PASS** (regression check) |
| **KA-4 / P5** | label-shuffled null ≈ 0 ± 0.02 | mean **−0.0023**, median −0.0043, sd 0.0945 (Stim48hr, D1+D2\|D3+D4, n = 253) against a real median of +0.0421 | **PASS** |
| **KA-5** | pooled-4-donor cells per (panel target, condition) reconciles with 278 §1.1 | pseudobulk (all guides) **Rest 500 / Stim8hr 476 / Stim48hr 481**; DE table 512 / 490 / 495.5 — the residual is the denominator (297 present targets vs the DE table's 283/287/284), and KA-1 shows the per-row agreement is exact | **PASS** |
| **KA-6** (cross-estimator) | my across-donor cosine should track the **authors' own** `donor_correlation_all_mean` (Pearson of per-gene DESeq2 log2FC across disjoint donor pairs), a completely independent estimator computed by the study | Spearman **+0.700** (Rest, n = 49), **+0.852** (Stim8hr, n = 57), **+0.697** (Stim48hr, n = 61) | **PASS** |

### 0.1 KA-2 diagnosed — the gap is the estimator, not a bug

The bank's `delta` is `ln(2) × DESeq2 shrunken log2FC`; mine is a `log1p(5e4·CPM)` pseudobulk difference. Two
checks separate "different functional" from "wrong pipeline", and both come out on the side of the functional:

* **on-target sign** — my pooled four-donor delta at each target's **own** gene (a column excluded from every
  scored statistic) is **median −0.769 (Rest) / −0.779 (Stim48hr), negative on 98.2 % / 98.9 % of the 283/284
  panel targets**. The pipeline recovers CRISPRi knockdown of the intended gene.
* **expression stratification** — `log1p(CPM)` and a log-fold-change coincide only where CPM is not small. The
  cosine to the bank rises monotonically with baseline expression and **crosses the preregistered 0.50 bar**:

  | genes kept (by median DESeq2 `baseMean`) | 9,828 (all) | 4,914 (top 50 %) | 1,966 (top 20 %) | 983 (top 10 %) |
  |---|---:|---:|---:|---:|
  | Rest, median cosine | +0.353 | **+0.510** | +0.527 | +0.537 |
  | Stim48hr, median cosine | +0.391 | **+0.529** | +0.525 | +0.521 |

  The median norm ratio `|bank| / |mine|` is 3.14–3.48, i.e. the two estimators differ mostly in how hard they
  down-weight low-expression genes. **Verdict: the preregistered bar was set too high for a cross-estimator
  comparison on `log1p(CPM)`'s gene weighting; every diagnostic the bar was meant to protect against
  (sign error, wrong axis, wrong control pool) is separately cleared.** Recorded as a half-pass, not as a pass.

---

## 1. The inventory and the acquisition plan (mandate item 1)

### 1.1 Every cell-level and aggregate object, with sizes, from the bucket listing

`s3://genome-scale-tcell-perturb-seq/marson2025_data/` — anonymous HTTPS, **MIT licence** (stated on the CZI
Virtual Cells Platform dataset page and in `marson2025_data/data_sharing_readme.md`), 57 keys, 1,843.66 GB total.

| object | bytes | what it is |
|---|---:|---|
| `D1_Rest.assigned_guide.h5ad` | **142,828,875,662** | cell-level, donor × condition, `.obs` = `lane_id, n_genes_by_counts, total_counts, pct_counts_mt, top_guide_UMI_counts, guide_id, perturbed_gene_name, perturbed_gene_id, guide_type, PuroR, guide_group, low_quality` |
| `D2_Rest` / `D3_Rest` / `D4_Rest` | 140,574,804,978 / 138,149,946,574 / **118,609,904,492** | **the four-donor Rest set alone = 540.2 GB** |
| `D*_Stim8hr` (4) | 154.8 / 172.8 / 142.7 / 134.6 GB | |
| `D*_Stim48hr` (4) | 156.2 / 137.6 / 137.8 / 159.2 GB | |
| **`GWCD4i.pseudobulk_merged.h5ad`** | **44,566,657,140** | **278,684 pseudobulks × 18,129 genes, one per (guide, donor, condition)**; `X` **CSR, uncompressed**, `data` f8 / `indices` i8 / `indptr` i8, 2,782,896,823 nnz; `.obs` carries `10xrun_id, donor_id, culture_condition, guide_id, guide_type, perturbed_gene_name/id, n_cells, total_counts, keep_{min_cells,effective_guides,total_counts,for_DE,test_genes}` |
| `GWCD4i.DE_stats.by_guide.h5mu` | 29,424,424,894 | per-guide DESeq2 fits, 2 modalities |
| **`GWCD4i.DE_stats.by_donors.h5mu`** | **16,866,278,447** | **per-donor-pair DESeq2 fits**, one modality per pair |
| `GWCD4i.DE_stats.h5ad` | 16,786,240,107 | **already on disk** (`experiments/199-def-runbook/work/`, `c355f535…`) |
| 12 `metadata/marson_*.jsonld` | 5,781–5,793 each | per-file provenance |
| 26 `suppl_tables/*` | 1.2 kB – 45.3 MB | incl. `sample_metadata.suppl_table.csv`, `QC_summaries_per_sample_lane.csv`, `sgrna_library_metadata`, `guide_kd_efficiency` |

**GEO GSE314342 / SRA SRP643211** carries a *different* object: `GSE314342_RAW.tar` =
**171,154,268,160 B** = 284 CellRanger `*_sample_filtered_feature_bc_matrix.h5` files, one per
(10x run × lane × donor × condition), 0.04–0.76 GB each, **plus 52 bulk-RNA-seq CSVs** for the arrayed
validations. Per condition: **Rest 94 files / 54.58 GB**, Stim8hr 94 / 58.32 GB, Stim48hr 96 / 58.25 GB.
These carry **no published guide assignment** — they are the raw CellRanger outputs — so the GEO route would
mean re-deriving the assignment. **Rejected.**

### 1.2 The four-donor replicate unit at the least GB — and the confound that decides the primary condition

Four donors, from `sample_metadata.suppl_table.csv`:
**D1 `CE0008162`** (29 F, Hispanic) · **D2 `CE0010866`** (22 F, African American) ·
**D3 `CE0008678`** (23 F, African American) · **D4 `CE0006864`** (34 M, Asian). All `GEMX_flex_v1`, Ultima.

| condition | D1 | D2 | D3 | D4 | donor ⟂ 10x run? |
|---|---|---|---|---|---|
| Rest | R1 | R1 | R2 | R2 | **NO** — {D1,D2}\|{D3,D4} *is* the run split |
| Stim8hr | R1 | R1 | R2 | R2 | **NO** |
| **Stim48hr** | R2 | R2 | R2 | R2 | **YES — all four donors in one 10x run** |

**This is why the mandate's expected "ONE condition (Rest) × 4 donors" is the wrong choice, and it was
declared in the preregistration before any value was read: `Stim48hr` is the only condition in which the donor
is not aliased with the processing batch, so it is the primary.** (Measured afterwards: the run adds nothing
detectable — see P4 in §2.3 — so the Rest and Stim8hr numbers are usable too.)

### 1.3 What was actually acquired

**Nothing was downloaded at cell level.** The plan taken is a **byte-range extraction from
`GWCD4i.pseudobulk_merged.h5ad`**: `X` is uncompressed CSR, so the 17,495 rows that matter — every targeting
pseudobulk of a panel target (6,477) plus every non-targeting pseudobulk (11,018) — are individually
addressable. The file was pulled with 8 parallel HTTP range requests into `work/parts/p0…p7`
(**never concatenated**, to avoid a second 44.57 GB copy on a disk that was 95 % full; the parts were read as
one virtual file by `work/parts_file.py`), the 17,495 rows extracted to a dense `float32` bank, and the parts
deleted.

| item | value |
|---|---|
| source object | `https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data/GWCD4i.pseudobulk_merged.h5ad` |
| bytes transferred | **44,566,657,140** (exactly the listed size), sha256 of the assembled stream **`fd2b8c21d357f8699ec34e2d5ebc1639612c27a0147a9ca94d4983822d93247e`** |
| peak disk in `work/` | **46.0 GB** (parts 44.57 + bank 1.27 + extras 0.25), against the 60 GB cap; **1.4 GB retained** |
| terms | **MIT licence**; no registration, no auth |
| **retained** `bank280_X.npy` 17,495 × 18,129 f32, 1,268,667,548 B | `44b419ff1e150ef28066a224f53e8b961be96956c32204f5835585f33092754b` |
| `bank280_rows.npy` | `fdf0c060c3d4eade749c3057461b342b075ae8b826e6ab91ff9c1aa1da39d2a3` |
| `obs280.npz` (the full pseudobulk `.obs` + `.var` + `X.indptr`) | `852022b1a9e3f5264f6e9be7f7d05a5ca830ef6b9bd202fc6d325854223493cd` |
| `axis280.npz` (gene-axis map) | `5b3b2611496aca2ec625895d12c848ab2fdcbf57d1fb68860d7ba54843456e7c` |
| `extra280_X.npy` (500 random **non-panel** Stim48hr targets, 3,502 rows — the centring robustness check) | `8195cb26ca8c227bc46690fce7c127409182494a6177aae8be3b3fdbd8f36f14` |
| `de_stats_panel280.csv` (panel rows of the local DE table) | `4e44404ec9c2e6d728ff78d2298f22790988b6a97ce1708e6461226df0826494` |

**The declared secondary was reached, and cost 0 bytes of disk.** `GWCD4i.DE_stats.by_donors.h5mu`
(16,866,278,447 B) was **never downloaded**: its `layers/log_fc` datasets are **contiguous and uncompressed**,
so the 190 panel rows of each of the 6 donor-pair modalities were pulled by HTTP range request —
**1,034 requests, 1.08 GB transferred, 0 B written to disk**, 1,185 s. Results in §2.7. (Downloading it *and*
the pseudobulk would have been 61.4 GB, over the cap; the range route avoided the choice entirely.)

### 1.4 The design, measured from the `.obs`

* **297 / 300 panel targets** present (absent: **EEF1A2, EPHB2, FZD2**).
* Cells per **(panel target, donor, condition)**: **median 101**, q10 16, q90 294 — *one donor is a ~100-cell
  measurement*. Pooled over four donors: **median 482** (KA-5).
* Guides per (target, donor, condition): 2 for 2,969 of 3,508 cells (all guides); 2 for 2,348 of 3,324
  (`keep_for_DE`).
* Non-targeting control pool per (donor, condition): **827–980 NTC guides, 68,441–88,354 cells**.
* Gene axis: 17,771 of the official 18,533 symbols map 1:1 into the CD4 object (0 duplicated symbols);
  all 300 panel genes map and are zeroed; **17,471 scored columns**.

---

## 2. Q1 — how reproducible is a deep, Flex-chemistry response across human donors?

Response object (preregistered §2, `truth267.py` arithmetic): `p = log1p(5e4·S/ΣS)`, control `b` = the same
transform of the **side's own donors'** NTC pool, `delta = p − b`, residual `r = delta −` the column mean over
the 297 panel targets built the same way, panel positions zeroed, cosine on the 17,471 scored columns.
Primary guide filter `keep_for_DE`; all-guides reported as sensitivity where it matters.

### 2.1 The headline: donor halves, ~250 cells/side

| condition | pairing | n | median cosine | 95 % (paired target bootstrap) | mean | ≤ 0 | cells/side |
|---|---|---:|---:|---|---:|---:|---|
| Rest | D1+D2 \| D3+D4 **(run-confounded)** | 251 | +0.0507 | [+0.0375, +0.0628] | +0.0626 | 28.7 % | 230/265 |
| Rest | D1+D3 \| D2+D4 | 251 | +0.0480 | [+0.0312, +0.0726] | +0.0651 | 31.5 % | 234/246 |
| Rest | D1+D4 \| D2+D3 | 251 | +0.0454 | [+0.0331, +0.0625] | +0.0542 | 35.1 % | 254/248 |
| **Rest** | **POOLED** | 753 | **+0.0478** | **[+0.0387, +0.0573]** | +0.0606 | 31.7 % | |
| Stim8hr | D1+D2 \| D3+D4 **(run-confounded)** | 260 | +0.0368 | [+0.0248, +0.0452] | +0.0568 | 33.8 % | 231/248 |
| Stim8hr | D1+D3 \| D2+D4 | 260 | +0.0373 | [+0.0213, +0.0517] | +0.0580 | 33.8 % | 213/244 |
| Stim8hr | D1+D4 \| D2+D3 | 260 | +0.0447 | [+0.0280, +0.0632] | +0.0613 | 33.8 % | 234/240 |
| **Stim8hr** | **POOLED** | 780 | **+0.0397** | **[+0.0305, +0.0458]** | +0.0587 | 33.8 % | |
| **Stim48hr** | D1+D2 \| D3+D4 | 253 | +0.0421 | [+0.0241, +0.0652] | +0.0674 | 32.8 % | 242/264 |
| **Stim48hr** | D1+D3 \| D2+D4 | 253 | +0.0538 | [+0.0390, +0.0651] | +0.0662 | 28.9 % | 232/270 |
| **Stim48hr** | D1+D4 \| D2+D3 | 253 | +0.0599 | [+0.0458, +0.0691] | +0.0720 | 24.5 % | 246/253 |
| **Stim48hr (PRIMARY)** | **POOLED** | **759** | **+0.0539** | **[+0.0427, +0.0620]** | +0.0685 | **28.7 %** | |

**P1 is falsified, in the direction that matters.** The preregistered band was **[+0.08, +0.22]**, chosen
because these are 10x GEM-X Flex measurements — the challenge's own chemistry — at ~250 cells per side.
The measured value is **+0.0539**, *below* the band. Placed next to F1's table:

| reference | chemistry | split unit | cells/half | targets | median split-half residual cosine | ≤ 0 |
|---|---|---|---:|---:|---:|---:|
| H1 `h1-test` (R273 §2, F1) | **10x Flex** | lane parity | 200 | 85 | **+0.243** | **0 %** |
| H1 depth-matched (R273 §2) | 10x Flex | random cells | 80 | 85 | +0.114 | 0 % |
| **GWCD4i CD4 T, Stim48hr (280)** | **10x GEM-X Flex** | **biological donor** | **~250** | **253** | **+0.054** | **28.7 %** |
| **GWCD4i CD4 T, Stim48hr (280)** | 10x GEM-X Flex | **guide construct**, donors held identical | 347 | 97 | **+0.048** | 40.2 % |
| Orion HEK293T (R267 §3) | Orion | batch parity | 95 | 294 | +0.028 | 34 % |
| Orion HCT116 (R267 §3) | Orion | batch parity | 70 | 290 | +0.028 | 32 % |
| K562-GWPS (R267 §3) | 10x 3′ | gemgroup parity | 81 | 269 | +0.047 | 5.6 % |

**The finding.** *A 10x-Flex perturbation response measured in primary human CD4 T cells at ~250 cells per
half sits at the Orion end of F1's ladder, not the Flex end.* Chemistry is held fixed and cell count is
2.6–3.6× H1's; the reproducibility is 4.5× lower. **F1's "the Flex-400 reference resolves every target's
direction" is a property of the H1 reference — its cells and its perturbation responses — not of the Flex
platform and not of the 400-cell budget.** 278 §"one qualifier" argued that "chemistry and the recipient's own
biology carry part of the gap"; 280 measures the split, and it is the biology.

### 2.2 The preregistered 2×2 — where the disagreement actually lives (~50 cells/side, common counting budget)

Every pair below is built from the same (donor, guide) pseudobulks, both sides thinned to a per-target budget
of half the smallest available total, and centred on a **common** subset of targets present in all four donors
(so that a different-donor pair is not penalised by a different centring — which would have biased `Δ_donor`
toward the hypothesis).

| class | pair | contains | Rest (n=232) | Stim8hr (n=239) | **Stim48hr (n=241)** |
|---|---|---|---:|---:|---:|
| **D** | same donor, same guide, disjoint binomial count halves | counting noise only, same cells | +0.1145 [+0.1014, +0.1256] | +0.1320 [+0.1215, +0.1413] | **+0.1202 [+0.1121, +0.1290]** |
| **A** | **same donor, different guide** | + construct + cells | +0.0229 [+0.0190, +0.0278] | +0.0178 [+0.0145, +0.0201] | **+0.0214 [+0.0184, +0.0254]** |
| **B** | **different donor, same guide** | + donor + cells | +0.0251 [+0.0210, +0.0296] | +0.0178 [+0.0138, +0.0221] | **+0.0211 [+0.0166, +0.0259]** |
| **C** | different donor, different guide | + donor + construct + cells | +0.0211 [+0.0149, +0.0256] | +0.0149 [+0.0120, +0.0177] | **+0.0182 [+0.0132, +0.0217]** |

| contrast | Rest | Stim8hr | **Stim48hr (PRIMARY)** |
|---|---:|---:|---:|
| **`Δ_donor` = A − C** (donor, construct held constant) | +0.0019 [−0.0015, +0.0065] | +0.0030 [−0.0002, +0.0057] | **+0.0033 [+0.0004, +0.0072]** ✅ |
| `D − B` (everything but the construct) | +0.0894 [+0.0777, +0.0985] | +0.1142 [+0.1040, +0.1237] | +0.0991 [+0.0923, +0.1063] |
| `D − A` (everything but the donor) | +0.0915 [+0.0803, +0.1010] | +0.1141 [+0.1042, +0.1237] | +0.0988 [+0.0918, +0.1058] |
| `B − C` (construct, donors already differ) | +0.0040 [+0.0009, +0.0086] | +0.0030 [+0.0003, +0.0064] | +0.0030 [−0.0005, +0.0064] |

**Read this table as a variance decomposition.** Splitting the *same cells* by counts reproduces at **+0.12**.
Changing **anything** that changes *which cells* you look at — a different construct in the same person, the
same construct in a different person, or both — collapses it to **+0.018 … +0.021**. So at ~50 cells per side,
**80–85 % of the disagreement between two measurements of one perturbation is cell sampling**, and the donor
and the guide construct each add ≈ +0.003, indistinguishably from one another (A ≈ B to three decimals in every
condition). **P2 passes at the very bottom of its predicted [+0.005, +0.06] band, in 1 of 3 conditions.**

### 2.3 The donor test at the scale that matters (~350 cells/side) — added analysis

`Δ_donor` is measured at ~50 cells/side, which is not the regime the headline lives in. The eight
(donor, guide) pseudobulks of a target can be re-partitioned 4|4 in two ways with **identical cells on both
sides**: by **donor** (2 donors × both guides per side) or by **guide** (all 4 donors × one guide per side,
so the donors are identical on the two sides and only the construct differs).

| Stim48hr, `keep_for_DE`, n = 97 | median | 95 % | cells/side | ≤ 0 |
|---|---:|---|---:|---:|
| **DONOR split** (donors differ) | +0.0477 | [+0.0194, +0.0722] | 347 | 32.0 % |
| **GUIDE split** (donors identical, construct differs) | +0.0484 | [−0.0003, +0.0666] | 347 | 40.2 % |
| **`GUIDE − DONOR`** — donor is a real unit if > 0 | **+0.0007** | **[−0.0528, +0.0235]** | | |

All guides, n = 177: DONOR +0.0469, GUIDE +0.0494, `GUIDE − DONOR` **+0.0026 [−0.0261, +0.0238]**.
In **Rest** the difference is *negative and excludes zero* (`GUIDE − DONOR` = −0.0466 [−0.0677, −0.0108],
n = 91): there, switching the guide construct costs **more** than switching the donor — consistent with F1's
finding that the construct is a real unit the challenge truth (one construct per target) never averages over.

**Conclusion on Q1's pass rule.** The preregistered rule passes, narrowly and in one condition. The honest
statement is: **the biological donor is a real but very small unit of variance in this assay — the same size as
the guide construct, and roughly 1/30 of the cell-sampling term at the depths available.** Person-to-person
variation is *not* what stops a GWCD4i-derived response from being reproducible.

### 2.4 P4 — the run confound is not detectable

| condition | median(run-confounded pairing) − median(mean of the 2 run-balanced pairings) | n |
|---|---:|---:|
| Rest (donor ≡ run) | +0.0084 [−0.0087, +0.0252] | 251 |
| Stim8hr (donor ≡ run) | +0.0023 [−0.0155, +0.0163] | 260 |
| Stim48hr (all four in one run) | −0.0117 [−0.0265, +0.0074] | 253 |

**P4 is not supported**: the pairing that aliases donor with the 10x run does **not** read lower. The
processing batch contributes nothing measurable on top of the donor here, which is why the Rest and Stim8hr
rows in §2.1 can be quoted alongside the primary.

### 2.5 P3 — reproducibility tracks 276's conservation map

| condition | Spearman(across-donor cosine, `C_raw`) | Spearman(·, `C_rel`) | median cosine by `C_raw` tercile (low → high) |
|---|---:|---:|---|
| Rest | **+0.155** (p = 0.014, n = 251) | +0.098 (p = 0.18) | +0.0405 / +0.0460 / +0.0499 |
| Stim8hr | **+0.295** (p = 1.2e-06, n = 260) | +0.144 (p = 0.046) | +0.0127 / +0.0327 / **+0.0767** |
| Stim48hr | **+0.295** (p = 1.8e-06, n = 253) | +0.197 (p = 0.0067) | +0.0291 / +0.0324 / **+0.1005** |

**P3 passes** — the point estimates sit inside the predicted [+0.15, +0.40]. A target whose response is
conserved across 19 cell lines (276) is also the target whose response reproduces across four human donors.
This is a **fourth independent surface** for F8 (276's map reproduces across cells, across disjoint line sets,
under PLOLO, and — now — across biological donors of a primary cell type that is in none of those line sets).

### 2.6 Robustness

* **Centring** (panel-centred vs centred on 500 random **non-panel** targets, Stim48hr donor halves):
  +0.0421/+0.0418, +0.0538/+0.0491, +0.0599/+0.0589. The residual centring is not driving anything.
* **Null** (KA-4): −0.0023 ± 0.0945 against a real +0.0421.
* **Cross-estimator** (KA-6): Spearman +0.70 / +0.85 / +0.70 with the study's own DESeq2 donor-pair correlation.
* **Two ways of pooling the three pairings.** §2.1's POOLED row is the median over all `3 × n` cosines;
  averaging the three pairings **per target first** gives +0.0473 [+0.0357, +0.0593] (Rest),
  +0.0341 [+0.0201, +0.0472] (Stim8hr), **+0.0471 [+0.0332, +0.0643]** (Stim48hr). Every per-target analysis
  below (§2.5, §2.7, §5) uses the per-target average; the two agree to ±0.007.

### 2.7 The headline replicated in the authors' own estimator (the declared secondary)

`GWCD4i.DE_stats.by_donors.h5mu` holds a **separate DESeq2 fit inside each of the six donor pairs**
(`CE0006864_CE0008162` … ), i.e. exactly §2.1's construction with DESeq2 log2 fold changes instead of a
`log1p(CPM)` pseudobulk difference. Same residual centring, same official non-panel gene axis (10,273 authors'
genes → the scored subset).

| condition | D1D2\|D3D4 | D1D3\|D2D4 | D1D4\|D2D3 | **POOLED** | ≤ 0 | n |
|---|---:|---:|---:|---:|---:|---:|
| Rest | +0.1210 | +0.1244 | +0.0822 | **+0.1093 [+0.0673, +0.1349]** | 12.5 % | 56 |
| Stim8hr | +0.0767 | +0.0983 | +0.0988 | **+0.0924 [+0.0635, +0.1120]** | 10.6 % | 66 |
| **Stim48hr** | +0.1274 | +0.1148 | +0.0965 | **+0.1100 [+0.0799, +0.1314]** | 8.8 % | 68 |

**The `n` is the point.** A modality only contains a (target, condition) row if it passed DE eligibility
**inside that donor pair**, so a disjoint pairing retains only **56–68 of the 297 panel targets** — the
strongest fifth. On **exactly those** targets my pseudobulk statistic reads *higher*, not lower:

| condition | n (eligible in both pairs) | authors' DESeq2 | my pseudobulk, same targets | difference | my pseudobulk, all 251–260 targets |
|---|---:|---:|---:|---:|---:|
| Rest | 49 | +0.1104 [+0.0652, +0.1448] | +0.1544 [+0.0799, +0.1931] | −0.0440 [−0.0792, +0.0113] | +0.0473 |
| Stim8hr | 57 | +0.0956 [+0.0699, +0.1378] | +0.1408 [+0.1127, +0.1885] | −0.0452 [−0.0887, −0.0160] | +0.0341 |
| Stim48hr | 61 | +0.1115 [+0.0832, +0.1528] | +0.1387 [+0.1062, +0.1898] | −0.0272 [−0.0727, −0.0009] | +0.0471 |

Per-target Spearman between the two estimators: **+0.702 / +0.864 / +0.753**.

**Three conclusions.**
1. **The headline is not an artefact of the `log1p(CPM)` estimator.** The study's own DESeq2 donor-pair fits
   give the same picture; where they differ, mine is the *more* reproducible of the two.
2. **DE-eligibility selection is worth about +0.10 of cosine** (+0.107 / +0.107 / +0.092). The ~20 % of panel
   targets whose response is strong enough to be fitted inside a two-donor half reproduce across donors at
   **≈ +0.14**; the panel as a whole reproduces at **≈ +0.05**. Any number quoted from this dataset must say
   which population it is over — a reading F1's "reliability is a target property" predicts and 262's group
   structure requires.
3. **Even the selected fifth stays below H1.** +0.14 (best fifth of the panel, CD4) against **+0.243** (H1's
   median over *all* 85 of its targets). The gap in §2.1 survives the most favourable selection available.

---

## 3. Q2 — replicates vs cells (the depth/replicate axis, programme rule 7)

Leave-one-donor-out: truth = the held donor's own residual row (all its cells); prediction = the mean of the
residual rows of a size-`k` subset of the other three, averaged over all `C(3,k)` subsets and all 4 held
donors. **With four donors the largest honest training set is three**, so the mandate's "four-donor mean" is
not estimable under leave-donors-out and the preregistration replaced it with the `k ∈ {1,2,3}` curve.
`matched` = each training donor thinned to `1/k` of the per-target UMI budget, so the **total** counts behind
the prediction are equal at every `k`.

| condition | mode | k=1 | k=2 | k=3 | **Δcos(3−1)** | Δcos(2−1) |
|---|---|---:|---:|---:|---:|---:|
| Rest | full | +0.0264 | +0.0378 | +0.0436 | +0.0172 [+0.0128, +0.0237] | +0.0114 [+0.0067, +0.0146] |
| Rest | **matched** | +0.0230 | +0.0283 | +0.0301 | **+0.0071 [+0.0045, +0.0105]** | +0.0053 [+0.0028, +0.0081] |
| Stim8hr | full | +0.0265 | +0.0324 | +0.0373 | +0.0108 [+0.0071, +0.0185] | +0.0059 [+0.0028, +0.0102] |
| Stim8hr | **matched** | +0.0219 | +0.0255 | +0.0274 | **+0.0055 [+0.0038, +0.0120]** | +0.0036 [+0.0007, +0.0087] |
| **Stim48hr** | full | +0.0291 | +0.0354 | +0.0421 | +0.0129 [+0.0098, +0.0176] | +0.0063 [+0.0038, +0.0102] |
| **Stim48hr** | **matched** | +0.0226 | +0.0270 | +0.0294 | **+0.0067 [+0.0040, +0.0108]** | +0.0043 [+0.0016, +0.0072] |

**P6 passes** (positive, interval excludes 0, all three conditions), below its predicted [+0.01, +0.08] band.
**P7 passes** (the full-cells curve is monotone in `k`). n = 251 / 260 / 253 targets present in all four donors.

**The exchange rate between a replicate and a cell.** Full-cells Δcos(3−1) = +0.0129; matched-budget
Δcos(3−1) = +0.0067. **52 % of the value of going from one donor to three is genuine replicate averaging;
48 % is just three times as many cells.** By 276's conservation tercile (Stim48hr, matched): +0.0113 /
+0.0068 / +0.0113 — the replicate gain is present in every tercile and is not a property of the conserved
targets alone, even though their absolute reliability differs by 2.5× (k=1 cosine +0.0142 / +0.0199 / +0.0356).

---

## 4. The depth curve — would more cells reach the Flex-400 regime? (added analysis)

Stim48hr, both splits thinned to a ladder of counting budgets (`x` = cell-equivalents per side of the guide
split; the donor split's own sides carry ~250 cells at frac = 1):

| fraction of counts | cell-equivalents/side | **GUIDE split** (donors identical) | **DONOR split** | donor ≤ 0 |
|---:|---:|---:|---:|---:|
| 0.0625 | 21.7 | +0.0032 [−0.0005, +0.0135] | +0.0126 [+0.0088, +0.0157] | 24.9 % |
| 0.125 | 43.4 | +0.0135 [+0.0007, +0.0273] | +0.0146 [+0.0113, +0.0188] | 26.1 % |
| 0.25 | 86.8 | +0.0159 [−0.0038, +0.0338] | +0.0263 [+0.0182, +0.0312] | 27.7 % |
| 0.5 | 173.5 | +0.0346 [−0.0082, +0.0517] | +0.0356 [+0.0290, +0.0459] | 26.9 % |
| **1.0** | **347.0** | **+0.0484 [−0.0003, +0.0666]** | **+0.0471 [+0.0332, +0.0643]** | 28.1 % |

Two things fall out.

1. **The two curves are the same curve.** Changing the donor and changing the construct cost the same amount at
   every depth — the second, independent confirmation of §2.2's `A ≈ B`.
2. **Depth does not close the gap to H1.** Spearman–Brown from the measured +0.0471 gives **+0.054 at 400
   cells/side, +0.102 at 800 and +0.222 at 2,000** — 400 cells, the challenge's own budget, buys essentially
   nothing here. Reaching H1's measured +0.243 needs a **6.49× budget: ~2,250 cell-equivalents per side,
   ~4,500 cells per (target, condition) — 9.4× GWCD4i's median 482**, on a screen that already loaded 22 M
   cells over 11.5 k targets. **The mandate's depth axis is answered: for this cell type, cells are not the
   binding constraint at any budget a genome-scale screen can reach.**

And the descriptive statistics say the same thing directly: **Spearman(across-donor cosine, cells per target)
is −0.03 to +0.07 (p ≥ 0.26 in all three conditions), while Spearman(across-donor cosine, number of DE genes)
is +0.38 / +0.49 / +0.59.** What makes a target's response reproducible across donors is **how big the response
is**, not how many cells were spent on it.

---

## 5. Descriptive tail (post-hoc, no pass rule)

`work/desc280_per_target.csv` carries every panel target's across-donor cosine per condition joined to 276's
map and the DE table.

**Most donor-VARIABLE** (lowest across-donor cosine): Stim48hr — `ZNF784 −0.183, TWF2 −0.154, ZNF670 −0.115,
SFXN4 −0.111, RBM4B −0.103, STT3A −0.101, BCKDK −0.090, DTNBP1 −0.088, DBN1 −0.083, ZNF22 −0.082`;
Rest — `HSBP1 −0.162, STT3A −0.152, SFXN4 −0.148, HMGN3 −0.141, ZNF627 −0.140, ZFP62 −0.137`;
Stim8hr — `SFXN4 −0.261, RBM4B −0.160, TWF2 −0.125, ZNF32 −0.121, TMF1 −0.114`.
**The tail is enriched for C2H2 zinc fingers** (ZNF784, ZNF670, ZNF627, ZNF496, ZNF32, ZNF22, ZFP62) — the same
family 279 is testing for state-dependence, reached here from a completely different direction (donors, not
states). `SFXN4`, `STT3A`, `TWF2` and `RBM4B` are donor-variable in two or three conditions.

**Most donor-REPRODUCIBLE**: `ACLY +0.527, PCGF1 +0.502, PLCG1 +0.484, MED13 +0.464, SUPT7L +0.445,
STAT6 +0.442, DNTTIP1 +0.436, CAND1 +0.392, SIN3B +0.381, MED15 +0.366` (Stim48hr) — Mediator, SAGA/ATAC,
Sin3, and the two canonical CD4 signalling nodes (`PLCG1`, `STAT6`).

**Does donor variability track 276's low-conservation tail?** Yes, in all three conditions: median across-donor
cosine in 276's bottom-30 % `C_raw` vs the rest — **+0.0405 / +0.0487** (Rest), **+0.0132 / +0.0467**
(Stim8hr), **+0.0373 / +0.0582** (Stim48hr).

Other Spearman correlates of the across-donor cosine (Stim48hr): `n_total_de_genes` **+0.585**,
authors' `guide_correlation_all` **+0.476**, `C_raw` +0.295, `shared_programme_share` +0.144,
`ontarget_effect_size` −0.042 (n.s.), **`n_cells_target` −0.028 (n.s.)**.

---

## 6. What this establishes, what it does not, and what it costs the programme

**Established (scientific track).**
1. **The replicate ceiling of the challenge's own chemistry in a primary human cell type is +0.054 at ~250
   cells/side, with 29 % of panel targets at or below zero.** Four biological donors, three cell states, the
   scorer's own residual geometry, a validated null and a cross-estimator check.
2. **F1's Flex-400 regime does not generalise.** Same chemistry, more cells, 4.5× lower reliability. Any
   argument of the form "Flex + 400 cells ⇒ a resolving reference" must now name the cell type.
3. **Person-to-person variation is a real but marginal unit** (+0.003 cosine, same size as the guide
   construct, ~1/30 of cell sampling). The programme's reliability law (F1) does not need a donor term.
4. **Replicates beat cells, but only just**: at a fixed sequencing budget, three donors beat one by +0.0067
   cosine; a bit over half of the three-donor advantage survives matching for cells.
5. **F8's conservation map holds on a biological-donor surface** (ρ +0.155/+0.295/+0.295, and +0.340/+0.290 in
   the authors' own estimator on Rest/Stim8hr), the first surface in the programme that is neither a cell line
   nor a split of one experiment.
6. **Depth is not the binding constraint in CD4 T cells** — reaching H1's reliability would need ~9× the
   cells; response magnitude, not cell count, predicts across-donor reproducibility.
7. **The result is estimator-independent** — the study's own donor-pair DESeq2 fits reproduce it (§2.7), and
   **DE-eligibility selection, not estimator choice, is what moves the number** (+0.10 of cosine).

**Not established.** Nothing here is a competition claim, and nothing was priced. The cosines are on a
pseudobulk `log1p(CPM)` estimator that agrees with the study's DESeq2 rows at +0.35–0.39 overall (+0.51–0.53 on
well-expressed genes) — a different functional, so absolute values are not interchangeable with the bank's.
Four donors is a small replicate set: every donor statement rests on n = 4 people (2 pairings of which are
run-balanced), the 2×2 lives at ~50 cells/side, and the matched-scale test has n = 97–177 targets and wide
intervals. `Stim48hr` is the only run-clean condition. The within-donor comparator is a **guide** split, not a
cell split — cell-level data would cost 540 GB for one condition and was correctly out of scope.

**What it costs the programme.** 278 recommended extending on the depth axis. **280 says the depth axis in
GWCD4i is closed**: the only panel-covering deep source we have is a cell type in which a doubling of cells
buys +0.005 of reliability, and the four-donor replicate structure — the thing 267 wanted and could not get —
turns out to be worth +0.0067 cosine at matched budget. A stage 2 that wants a *resolving* reference for these
targets needs a different **cell type**, not more cells or more donors of this one. The remaining acquisition
worth naming is `GWCD4i.DE_stats.by_donors.h5mu` (16.87 GB, MIT) if the authors' own donor-pair DESeq2 rows are
wanted as a second estimator; it is not needed for any claim above.

---

## 7. Files

`PREREGISTRATION.md` (sha256 `7273ef3c…`, hashed 06:57:35Z) · `work/prereg280.sha256` ·
`work/kept_sha256.txt`, `work/src_sha256.txt` — **scripts**: `remote_h5.py` (HTTP-range file object),
`parts_file.py` (multi-part virtual file), `probe280.py`, `obs280.py`, `axis280.py`, `inv280.py`,
`extract280.py`, `extra280.py`, `verify280.py`, `lib280.py`, `ka280.py`, `ka2_280.py`, `ka2b_280.py`,
`ka6_280.py`, `q1_280.py`, `q1b_280.py`, `q2_280.py`, `depth280.py`, `desc280.py`, `probe_mu.py`,
`mu280.py`, `dl280.sh` —
**tables**: `q1_halves_keepDE.csv`, `q1_2x2_keepDE.csv`, `q1b_matched_scale_{keepDE,allguides}.csv`,
`q1_conservation_keepDE.csv`, `q2_curve_keepDE.csv`, `q2_terciles_keepDE.csv`, `depth280_curve.csv`,
`mu280_donor_pairs.csv`, `mu280_pertarget_{Rest,Stim8hr,Stim48hr}.csv`,
`desc280_per_target.csv`, `de_stats_panel280.csv`, `ka280_obs.json`, `ka2_280.json`, `ka6_280.json`,
`mu280.json`, `q1_280_keepDE.json` — **arrays**: `bank280_X.npy`, `bank280_rows.npy`, `obs280.npz`,
`axis280.npz`, `extra280_X.npy`, `q1_pertarget_*.npz`, `q1_2x2_*.npz`, `q1b_*.npz`, `q2_*.npz`, `ka2_*.npz` —
**metadata**: `meta/` (S3 readme, sample metadata, per-lane QC, table-name map), `geo_filelist.txt`,
`geo_suppl_listing.html`, `s3_list_1.xml`, `rss_q1.txt`, `extract280.log`, `q1_280.log`, `q2_280.log`,
`mu280.log`.
