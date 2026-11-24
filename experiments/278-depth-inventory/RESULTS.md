# 278 — The DEPTH axis, inventoried: **214 of the 300 panel targets already have a ≥ 400-cell single-cell measurement in a somatic cell on this disk today** — but 152 of them ONLY in primary CD4 T cells, and only **62** in an immortalised somatic cell line

**The one line the mandate asks for.** *How many of the 300 panel targets can reach a ≥ 400-cell somatic
measurement with data that exists today, and from where?* — **214 / 300, from files already on disk, with no
download at all**: **189 through the genome-scale primary-CD4-T screen (152 of them through nothing else)**
(`GWCD4i`, Zhu/Dann/Marson,
**10x GEM-X Flex v1 on Ultima — the challenge's own chemistry**, median **512 / 490 / 495.5** cells per panel target
in Rest / Stim8hr / Stim48hr), 62 through an immortalised cell line (Orion HEK293T 33, Altair A549 21, Orion HCT116 9,
Xu HEK293 5, Frangieh melanoma 4, K562 3+1, Calu-3 1), and a further 14 only if iPSC-derived cardiomyocytes count
(IGVF, 56 panel targets at ≥ 400 cells) — **228** on the widest somatic reading, **231** counting H1/iPSC too.
**No public dataset that can be downloaded today raises that number.** The one that would — X-Atlas/**Pisces**
(HepG2, iPSC, and **Jurkat resting/active**, i.e. challenge line A) — has **24.3 kB in its repository**,
re-verified 2026-09-06.

**The finding that matters for the programme.** On this panel the DEPTH axis and the CONTEXT axis are
**confounded, and not by a little**: the median panel target's deepest measurement in a somatic cell is **582
cells**, but strip out primary CD4 T cells and it falls to **255**; strip CD4 out of the ≥ 400 count and it falls
from 214 to **62**. Only **13** panel targets have a ≥ 400-cell measurement in **two different somatic contexts**
once CD4's three conditions are treated as one source (183 with CD4's three counted separately). So "more cells per
condition" is available for most of the panel — in exactly one cell type, the one furthest from the challenge's
cancer lines. F1's reliability law says a 400-cell Flex reference resolves every target's direction (median
split-half residual cosine +0.243 at 200/200 cells, **0 / 85 ≤ 0**; Spearman–Brown 0.330 measured at 200/200,
0.46–0.50 fitted at 400, √ρ ≈ 0.68–0.71) against 0.05–0.09 for the 150–200-cell donors the pipeline runs on; this
inventory says that regime is reachable for 71 % of the panel **only by accepting primary CD4 T cells as the
context**.

**One qualifier that F1 itself insists on, and it cuts against a naive depth argument.** F1's headline is *"and the
difference is not depth"*: depth-matched at 80 cells per half the H1 Flex reference still reads **+0.114** (0 % of
targets ≤ 0) against Orion's **+0.028** at 70–95 cells per half (32–34 % ≤ 0). So cells alone do not buy Flex-like
reliability — chemistry and the recipient's own biology carry part of the gap. That makes the CD4 rows more
interesting than a bare cell count suggests (they are Flex, on Ultima, at ~500 cells) and makes "just take more
cells in an Orion-like assay" a weaker plan than the reliability law read on its own would imply.

**And the deepest data we own is already inside the incumbent — as a residual, not as a response estimate.**
The CD4 banks are the source of the `cd4svd5` residual, whose leading direction is deflated out (`VCC_CD4_DEFLATE_K=5`,
READ-FIRST §2). The ≥ 400-cell measurements are being spent on a hub-removal operation, not on a per-target response
whose reliability is exploited. That, not acquisition, is the first thing a DEPTH stage-2 should look at.

**Date 2026-09-06, 06:18–06:50 UTC (09:18–09:50 Ankara).** Metadata and public-catalogue work only.
**No build, no `vcc`, no `cell-eval2`, no model, no submission, nothing under `submissions/` or `data/`, no
production file touched, no count matrix loaded, no download of any kind** (every byte read was already on disk;
the only network traffic was web search / page fetch of public catalogue pages). `uv run python` throughout,
`PYTHONDONTWRITEBYTECODE=1`. `find data scripts src docs submissions -newer work/.marker-start` returns **nothing**;
inside `experiments/` only the sibling 279 agent's own files.

> **Memory, honestly: `scan278_h5ad.py` peaked at 9.40 GiB, above the mandate's ~4 GB cap** (it holds one obs
> column as a 2 M-element Python-object array while `np.unique`-ing it; `GasperiniShendure2019_atscale` and
> `SunshineHein2023` are where it climbs). It stayed under the machine-wide ~10 GB rule and ran for 45 s.
> Every other stage peaked under 1 GiB (`build278.py` 0.78 GiB). The fix, if this is re-run, is to convert each
> obs column to a categorical code array before the unique — I did not re-run it because the scan is cached in
> `work/h5ad_counts/` and re-running would burn the budget for no new number.

---

## 0. Known-answer rows (all read before any new number was believed; `work/ka278.txt`)

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** an independent obs scan of `ReplogleWeissman2022_K562_gwps.h5ad` reproduces 267's per-target cell counts | exact | **272 / 272 targets, max \|diff\| = 0.0** | PASS |
| **KA-2** per-source panel coverage vs 276 §1.1 (`inventory276_sources.csv`) | 18 sources | **15 / 18 exact**; the 3 residuals explained: `altair_a549` 21 raw rows → 20 after the 225/276 builder drops 3 all-NaN rows (COQ10B, HSD17B12, NEK1); `myllia_a549` 15 raw → 14 after the same rule; `feng_ipsc` 166 (gw non-fitness cell metadata) + 31 (gw fitness) = 197 union vs 182 in the bank, which 099 built from the figshare **LFC-by-gene** table, a different (larger) object than the two cell-metadata files | PASS |
| **KA-3** R267 §1's cells-per-fold-target quantiles (min / q25 / median / q75 / max) | (23,127,175,226,1745) · (20,91,153,227,648) · (21,118,203,316,819) | **exact on K562 and HCT116**; HEK293T **(21,118,202.5,316,819)** — 267 printed the midpoint 202.5 as 203 (n = 294 is even), every other value identical | PASS (a display rounding) |
| **KA-4** R270 §4's Feng statements | gw median 22 cells/target; 0 panel targets with ≥ 20 cells in ≥ 2 donor lines; targeted screen 5 panel targets at 12–17 donor lines with ≥ 50 cells | **22**; **0**; **5 targets at 12 / 15 / 15 / 15 / 17 lines** | PASS |
| **KA-5** the h1-test harness is the 2026 reference regime | 85 targets × exactly 400 cells + 18,400 non-targeting; ~20 k UMI | **85 targets at exactly 400, 18,400 non-targeting, median 19,997 UMI/cell** (measured from `X`'s row sums) | PASS |
| **KA-6** CD4 `guide_n` summed over a target's guides == `DE_stats.n_cells_target` | exact | **283 / 283 panel targets in Rest, max \|diff\| = 0.0** | PASS |

---

## 1. Item 1 — on disk, per panel target × line: cells per condition

`work/depth278.csv` is the deliverable: **3,169 rows**, one per (panel target, source), columns
`target, source, line, state, chemistry, unit, cells, control_cells, guides, batches, depth_source, note,
duplicate_of`. `work/depth278_per_target.csv` collapses it to the 300 targets with the deepest somatic /
somatic-or-derived / any-state measurement and its source. `work/depth278_summary.json` carries the thresholds.
Every `cells` value is read from an obs table, a supplementary table shipped with the dataset, or a QC table
shipped with the dataset — **never from a count matrix**. Sources whose on-disk object is a published DE / lfc /
z-score matrix with no cell count are carried with `cells` empty and `depth_source = not-on-disk`, so gaps are
visible rather than silently absent.

`state` classifies the recipient: **somatic** = a differentiated cell (immortalised line or primary),
**derived** = a somatic cell type differentiated from a pluripotent source (iPSC-cardiomyocyte, iPSC-neuron),
**pluripotent** = hESC / iPSC. The mandate's "somatic" count is the first of these.

### 1.1 The depth table by source (measured cells per panel target)

| source | line | state | panel targets | median cells | max | **≥ 400** | ≥ 200 | ≥ 100 | control cells | guides/target |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `cd4_stim48hr` | CD4 T, Stim48hr | somatic (primary) | 284 | **495.5** | 3,997 | **177** | 252 | 272 | 313,989 (Rest) | median 2 |
| `cd4_stim8hr` | CD4 T, Stim8hr | somatic (primary) | 287 | **490** | 4,299 | **173** | 252 | 277 | — | median 2 |
| `cd4_rest` | CD4 T, Rest | somatic (primary) | 283 | **512** | 4,099 | **173** | 251 | 272 | 313,989 | median 2 |
| `igvf_cardiomyocyte` | cardiomyocyte (iPSC-derived) | derived | 87 | **514** | 1,584 | **56** | 72 | 82 | — | — |
| `orion_hek293t` | HEK293T | somatic | 294 (300) | 202.5 | 819 | **33** | 150 | 235 | 219,071 | 284×1, 10×2 |
| `altair_a549` | A549 | somatic | 21 (20 in the bank) | **631** | 946 | **21** | 21 | 21 | 45,320 | 4 |
| `h1_2025` | H1 hESC | pluripotent | 20 | 400 | 400 | 20 | 20 | 20 | 18,399 | — |
| `orion_hct116` | HCT116 | somatic | 290 (300) | 153 | 648 | **9** | 95 | 203 | 165,931 | 280×1, 10×2 |
| `xucao_hek293` | HEK293 | somatic | 5 | 837 | 1,256 | 5 | 5 | 5 | — | — |
| `feng_targeted` | iPSC (HipSci, 19 donor lines) | pluripotent | 5 | 1,544 pooled | 1,923 | 5 | 5 | 5 | 8,241 | — |
| `frangieh_melanoma` | melanoma (patient) | somatic | 5 | 602 | 1,203 | 4 | 5 | 5 | — | — |
| `k562_gwps` | K562 | somatic | 269 (272) | 175 | 1,745 | **3** | 95 | 223 | 75,328 | 261×1, 8×2 |
| `sunshine_calu3` | Calu-3 | somatic | 6 | 225 | 700 | 1 | 3 | 6 | 3,536 | — |
| `tian_ipsc_neuron_crispra` | iPSC-neuron (**CRISPRa**) | derived | 4 | 250 | 592 | 1 | 3 | 4 | — | — |
| `norman_k562` | K562 (**CRISPRa**) | somatic | 3 | 241 | 516 | 1 | 3 | 3 | — | — |
| `tian2019_neuron` / `tian2019_ipsc` | iPSC-neuron d7 / iPSC | derived / pluripotent | 1 / 1 | 3,289 / 2,083 | — | 1 / 1 | 1 / 1 | 1 / 1 | — | — |
| `tian_ipsc_neuron` | iPSC-neuron | derived | 6 | 152 | 239 | 0 | 1 | 6 | — | — |
| `myllia_a549_ch2` / `ch1` | A549 | somatic | 15 / 15 | 94 / 85 | 230 / 190 | 0 | 1 / 0 | 6 / 6 | — | — |
| `wessels_thp1` | THP-1 | somatic | 1 | 78 | 78 | 0 | 0 | 0 | — | — |
| `igvf_telohaec` | TeloHAEC | somatic | 65 | 67 | 147 | 0 | 0 | 4 | — | median 8 |
| `feng_gw_fitness` | iPSC (HipSci, 24 lines) | pluripotent | 31 | 50 | 102 | 0 | 0 | 2 | not identifiable* | — |
| `feng_gw_nonfitness` | iPSC (HipSci, 34 lines) | pluripotent | 166 | 22 | 68 | 0 | 0 | 0 | not identifiable* | — |
| `gasperini_k562_lowmoi` / `_atscale` / `papalexi_thp1` | K562 / THP-1 | somatic | 3 / 1 / 1 | 18 / 1 / 4 | 27 / 1 / 4 | 0 | 0 | 0 | — | — |

Medians for the three fold banks are over the targets the bank keeps (≥ 20 cells); the parenthesised figure is the
number of panel targets with **any** cells in the file. `*` the two Feng genome-wide cell-metadata files carry
**281,334 / 218,616 unassigned cells** (71 % / 68 % of the file) and only **12 / 36** `NonTarget_` calls, so their
control set is not identifiable from this metadata and `control_cells` is left empty rather than guessed; the
targeted screen has **8,241** `NonTarget_` cells. **Excluding the two CRISPRa sources changes nothing**: the somatic
counts stay 214 / 285 / 296 (the two CRISPRa targets are also ≥ 400 in CD4); only the cell-line-only ≥ 400 count
moves 62 → 61.

Sources with panel coverage but **no per-target cell count on disk** (the on-disk object is a published effect
matrix): `kolf2` (273 targets, per-target z-scores), `perturbmatch` (251 — a reanalysis of the *same* K562 cells as
`k562_gwps`, so never counted twice), `dld1_lowmoi` (67, two low-MOI lfc replicates), `gse295214_h1` (64, a
pseudobulk database h5ad), `jiang_{A549,BXPC3,HAP1,HT29,K562,MCF7}` (9 each, per-line log2FC), `igvf_hcasmc`
(sceptre DE, no counts). `nadig_jurkat`, `nadig_hepg2`, `k562_essential`, `rpe1` and Broad-16 carry **0 / 300**
panel targets and contribute nothing to this axis (Broad-16 is 90 non-panel targets × 16 lines at a median of
42–129 cells per KO, R268 §1 — the same 100-cell regime).

### 1.2 The threshold counts (the deliverable)

Deepest **single** source per target (duplicated cells never double-counted):

| population | ≥ 400 | ≥ 300 | ≥ 200 | ≥ 150 | ≥ 100 | ≥ 50 | median deepest |
|---|---:|---:|---:|---:|---:|---:|---:|
| **somatic (immortalised lines + primary CD4 T)** | **214** | 262 | **285** | 294 | **296** | 300 | 582 |
| somatic + iPSC-derived (adds the IGVF cardiomyocytes) | 228 | 268 | 288 | 295 | 297 | 300 | 619.5 |
| any state (adds H1 / iPSC) | 231 | 269 | 288 | 295 | 297 | 300 | 630 |
| **somatic, primary CD4 T EXCLUDED** (cell lines + the patient melanoma) | **62** | — | **210** | — | **288** | — | 255.5 |
| the three production donors only (K562, HCT116, HEK293T) | 34 | — | 198 | — | 286 | — | 232 |

**Every one of the 300 panel targets has at least 50 cells somewhere.** **15 targets have no somatic measurement of
any depth ≥ 200** (`work/gap_lt200_somatic.csv`): ABCD1, ANKRD52, CAPRIN2, E2F4, EPHB2, EZH2, KIF21B, LENG1, MYPOP,
NICN1, PI4KB, PIN1, PRMT7, REV1, STK3. **86 targets are below 400 in every somatic source.**

### 1.3 Replication of the depth, not just its existence

| question | answer |
|---|---:|
| targets with ≥ 400 cells in ≥ 2 somatic **sources** | 183 |
| targets with ≥ 400 cells in ≥ 2 distinct somatic **contexts** | 183 |
| …**excluding CD4's three culture conditions** (i.e. two genuinely different cell types) | **13** |

The 183 is almost entirely "the same target measured deeply in CD4 Rest and CD4 Stim8hr". On the reading that
matters for a conserved-vs-state-dependent question — *two different cell types, both deep* — the panel has
**13 targets**, not 183.

---

## 2. Item 2 — the public record beyond disk

Method: catalogue-level checks of every resource that could plausibly hold a ≥ 400-cell single-cell measurement of
panel targets in a somatic cell, plus a re-verification of the two release watches 270 opened. **No expression
matrix was downloaded.** A gene-by-gene GEO sweep of the 300 symbols was *not* run and I state why below (§2.2).

| resource | somatic context(s) | panel coverage | depth (cells / perturbation) | chemistry | downloadable? terms | does it raise the 214? |
|---|---|---|---|---|---|---|
| **`GWCD4i` primary CD4 T** (Zhu, Dann, … Marson; bioRxiv 2025.12.23.696273, now *Cell* S0092-8674(26)00929-3; GEO **GSE314342** / SRA **SRP643211**; processed on the CZI Virtual Cells Platform) | primary CD4 T, 4 donors × {Rest, Stim8hr, Stim48hr} | **283–287 / 300** | **median 543 / 526 / 548 (Rest / Stim8hr / Stim48hr) over all 11,526 targets; 68.5 % of (target, condition) rows at ≥ 400**; 20,509,715 perturbed cells in the DE table, ~22 M total | **10x GEM-X Flex v1, Ultima** (`sample_metadata.suppl_table.csv`) | yes, **MIT licence**; cell-level `D*_*.assigned_guide.h5ad` per donor × condition | **it IS the 214** (already on disk as DE stats + the `cd4_compact` bank) |
| **X-Atlas/Orion** (Xaira; bioRxiv 2025.06.11.659105) | HCT116, HEK293T | 290 / 294 | **median ≥ 140**; 8 M dual-guide cells released out of ~30 M barcodes | 10x GEM-X 5′, median 16 k UMI | released; **fully on disk** (7,943,468 cells = the released set) | no — we already hold all of it |
| **X-Atlas/Pisces** (Xaira, announced 2026-03-17 with X-Cell) | HCT116, HEK293T, **HepG2**, iPSC, **Jurkat Resting**, **Jurkat Active**, iPSC multi-diff — 25.6 M cells, 7 genome-scale CRISPRi screens | unknown (genome-scale ⇒ ~all 300) | unknown | — | **NO. Repo is 24.3 kB, CC BY-NC-SA 4.0, "coming soon"** (re-verified 2026-09-06; 270 saw the same on 09-05) | **the only candidate that could — and it does not exist yet.** Jurkat = challenge line **A** |
| **Replogle 2022** (K562 GWPS / essential, RPE1) | K562, RPE1 | 269 / 0 / 0 | median 171–175 | 10x 3′ | on disk in full (1,989,578 cells) | no |
| **KOLF2.1J Perturbation Cell Atlas** (Nat Biotech 2026) | KOLF2.1J **iPSC** (pluripotent) | 273 | > 2.5 M cells / 11,692 genes ⇒ ~200–215 mean; 88 % of perturbations > 100 cells; > 5 k UMI/cell | 10x | on disk as per-target z-scores | no (and pluripotent) |
| **Feng 2025 HipSci** | iPSC, 34 donor lines (pluripotent) | 182 | **median 22** cells/target genome-wide; targeted screen 5 panel targets at median 86 cells per (target, line) | 10x | count data figshare 27989294, 7.93 GB, MIT | no — structurally cannot reach 400 |
| **VIPerturb-seq** (Bradu/Satija, bioRxiv 2026.02.12.705613; Zenodo 18460279) | **K562** (+ HEK293, HAP1 pilots) | genome-wide CRISPRi (GuEST-List) | not stated per perturbation; benchmarked against GWPS | **10x Flex v2**, probe-based, combinatorial indexing | Zenodo, Seurat `.rds` | no — K562 is already a fold line |
| **Jiang 2025** | A549, BxPC-3, HAP1, HT-29, K562, MCF7 × 5 cytokine states | 9 | on disk as per-line log2FC only; no cell counts | — | zenodo 14518762 | no |
| **Nadig 2025** | HepG2, Jurkat | **0 / 300** (essential library) | — | 10x | on disk (14 GB) | no |
| **Broad-16** | 16 cancer lines | **0 / 300** | 42–129 cells/KO median (R268 §1) | 10x | on disk | no |
| **IGVF** (TeloHAEC, cardiomyocyte, HCASMC) | TeloHAEC, iPSC-cardiomyocyte, HCASMC | 65 / 87 / — | TeloHAEC median 67; **cardiomyocyte median 514, 56 at ≥ 400** | 10x CRISPRi | on disk | +14 targets, but only if iPSC-derived cardiomyocytes count as somatic |
| **LINCS/CMap 2020** | 38 lines | 148 | **bulk L1000 consensus signatures — not single cells** | L1000 | on disk (17 GB) | no (out of scope: not single-cell) |
| **GPSAdb 2.0** | 1,063 lines | 96 | **bulk RNA-seq** | — | CC BY 4.0 | no (not single-cell) |
| **Tahoe-100M** | 50 cancer lines | **chemical only, 0 genetic** | 100 M cells | 10x | on disk (pseudobulk) | no — noted only, as the mandate asks |

### 2.1 The two release watches, re-verified today
* **X-Atlas/Pisces** `huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces`: **24.3 kB total**, no h5ad/parquet,
  dataset viewer reports no supported data files, CC BY-NC-SA 4.0. README's stated upload plan: all Orion dual-guide
  cells; HepG2 + iPSC cells restricted to a **200-gene validation set** or NTC; **and all dual-guide cells from
  Jurkat Resting and Jurkat Active**. Unchanged in substance since 270 checked on 2026-09-05.
* **X-Cell** remains a model release (CC BY-NC-SA 4.0), not a data release — off the plan, as 205 and 270 found.

### 2.2 Why no gene-by-gene GEO sweep, and what replaces it
The panel was selected (AGENTS.md §2.1) for *non-essential* genes with clean guide assignment and adequate
representation in existing public datasets — i.e. it is a **genome-scale-screen panel, not a famous-gene panel**
(READ-FIRST §7: panel ∩ essential list = 0 on A/B/C). A targeted screen deep enough to give one of these genes
400 cells would have to have chosen it deliberately, and the evidence that essentially none did is already
measured, not assumed: across **all 50 scPerturb files** scanned here at the obs level, the only one with more than
6 panel targets is K562-GWPS (272); every other file tops out at **6** (Tian iPSC-neuron 6, Sunshine Calu-3 6, Xu
HEK293 5, Frangieh 5+5, Tian CRISPRa 4, Norman 3, Gasperini low-MOI 3, Wessels 1, Papalexi 1+1, **and 35 of the 50
carry none at all**), and 270's independent audits found **0 / 300** in Broad-16 and
Nadig, **9 / 300** in Jiang, **6 / 300** in PerturbFate. The panel's depth therefore comes from *how many cells a
genome-scale screen loaded*, not from anyone targeting these genes on purpose — which is exactly why `GWCD4i`
(22 M cells over ~11.5 k targets ⇒ ~500 cells each) dominates the table and why nothing smaller can compete.

---

## 3. Item 3 — the chemistry / depth regime of the challenge's own references

**Primary source, on disk: `docs/reference/vcc2026-metrics-BRIEF.md` §0 (and `-FULL.md` lines 35–36, 653, 686).**

> "Per context: 300 target constructs, one per gene, plus a pooled non-targeting control drawn from 46 constructs.
> **Every perturbation is downsampled to exactly 400 cells and to a median depth of 20,000 UMI per cell**, on a
> common axis of 18,533 genes. Each context is scored independently against its own reference."

So the 400 is a property of **the reference itself**, not only of the submission format (the number of *predicted*
cells is explicitly unconstrained). The regime is: **10x Flex CRISPRi, exactly 400 cells per perturbation, median
20,000 UMI/cell, 18,533 genes, 46 NTC constructs, six cancer lines (A = JURKAT, B = HELA, C = CAL33 by exp 041's
CCLE match).** Arc's public material confirms the platform: CRISPRi + 10x Flex + Ultima UG100, chosen after a
pilot comparison of 3′ / 5′ / Flex on matched H1 perturbations.

**Which public data shares that regime:**

| set | Flex? | cells per perturbation | median UMI/cell | panel targets | verdict |
|---|---|---:|---:|---:|---|
| **A/B/C references (2026)** | yes | **exactly 400** | **20,000** | 300 (held) | the target regime |
| **`h1-test`** `data/harness/h1-test.h5ad` | yes (2025 VCC H1) | **exactly 400** (85/85) | **19,997 measured** | **5 / 85** | **the only object on disk that reproduces the regime exactly** — and it touches 5 panel targets |
| 2025 VCC H1, undownsampled (`vcc2025/{train,validation,test}`) | yes | median **1,045 / 1,090 / 1,141**, max 4,760 | ~54,000 | **25 / 300** | same chemistry, ~2.6× deeper, wrong (pluripotent) context |
| **`GWCD4i` CD4 T** | **yes — GEM-X Flex v1 on Ultima** | median **512 / 490 / 495.5** on the panel | not read (no cell-level file on disk) | **283–287** | **the only large-panel set that shares the challenge's chemistry AND reaches 400 cells** |
| X-Atlas/Orion | no (GEM-X 5′) | median ≥ 140 | 16,000 | 290 / 294 | close on UMI, 2.9× short on cells |
| Replogle GWPS | no (10x 3′) | median 175 | — | 269 | short on both |
| KOLF2.1J | no | ~200 | > 5,000 | 273 | short on both |

The h1-test harness is therefore **regime-exact and panel-empty**, and `GWCD4i` is **regime-matched and
panel-full but context-distant**. There is no object that is both.

---

## 4. Item 4 — the cost of depth

The gap to close: **86 panel targets below 400 cells in every somatic source**, of which **15 are below 200**.

### (a) By download

| candidate | what it would add to the ≥ 400 somatic count | GB | terms | verdict |
|---|---|---:|---|---|
| `GWCD4i` cell-level h5ads (12 files, donor × condition, ~22 M cells) | **+0** — the ≥ 400 rows are already derivable from the on-disk `DE_stats.suppl_table.csv` + the `cd4_compact` bank. What it *would* add is the **per-donor and per-guide split** 267 wanted and could not get (4 biological donors × 3 states is a real replicate unit, unlike the 8–10 two-construct targets in Orion) | 22 M cells ⇒ far over this mandate's 2 GB cap; per-file size not published | **MIT**, CZI Virtual Cells Platform | the only download worth a stage-2 mandate — **for replicate structure, not for depth** |
| **X-Atlas/Pisces** HepG2 / iPSC / **Jurkat Rest+Active** | potentially +the whole panel in 2–4 new somatic contexts at genome-scale depth, incl. challenge line A | unknown | CC BY-NC-SA 4.0 (non-commercial, share-alike) | **NOT AVAILABLE** (24.3 kB) — a watch item, not a plan |
| X-Atlas/Orion (figshare 29190726 / HF) | +0 — 194 measured `pass_guide_filter` at **100.0 % of 7,943,468 cells**, and 3,409,169 + 4,534,299 = 7,943,468 is exactly the released set | 0 | — | nothing to gain |
| Replogle GWPS "raw" (figshare, 61.3 GB) | +0 — same 1,989,578 cells as the 8.8 GB scPerturb h5ad on disk (KA-1) | 61.3 | CC BY 4.0 | nothing to gain |
| Feng count data (figshare 27989294) | +0 — median 22 cells/target genome-wide; 86 cells per (target, donor line) in the targeted screen | 7.93 | MIT | cannot reach 400 at any pooling |
| VIPerturb-seq (Zenodo 18460279) | +0 for new contexts (K562) | — | Zenodo | no new context |
| Tahoe-100M | +0 (chemical) | 1.2 on disk | — | noted only, as instructed |

### (b) By re-extraction from files already on disk

**First, the mandate's specific question: did `build_donor_table.py` cap cells per target? NO — verified by
reading the code.** `scripts/build_donor_table.py:111` takes `min_cells: int = 20` and uses it only as a **floor**
(`keep_names = [n for n, c in ... if n == control_label or c >= min_cells]`, line 125); `group_sums` then sums
**every** cell of every kept group (lines 130, 158: `cells = [int(counts[...]) for n in names]` is a record, not a
selection). `scripts/build_orion_donor.py` is the same with one extra filter, `pass_guide_filter == 1` (line 73).
R267 §1's claim that the full counts were used therefore stands, and KA-1 confirms it end-to-end.

| lever | measured | gain in the ≥ 400 somatic count |
|---|---|---:|
| relax `--min-cells 20` | it is a floor, not a cap; relaxing it only *adds* targets with < 20 cells (6 in HEK293T, 10 in HCT116, 3 in K562) | **0** |
| relax Orion's `pass_guide_filter` | 194 measured the pass rate at **1.0000** on 7,943,468 cells | **0** |
| re-extract the 2025 VCC H1 at full depth instead of the 400-cell `h1_2025` bank | 25 panel targets at median **1,008** cells (max 3,642) vs 20 targets at exactly 400 → **+5 targets, 2.5× depth** | **0 somatic** (H1 is pluripotent); it would make the *reference-regime instrument* deeper, not the somatic map |
| **re-extract IGVF cardiomyocyte from the raw data** | the on-disk object is a **significant-hits-only** DE table (median **6** non-zero genes per panel row) although **56 panel targets have ≥ 400 cells** and 82 have ≥ 100 | **the largest re-extraction prize on disk**: full profiles for 56 ≥ 400-cell targets in a somatic (iPSC-derived) context. Cost: the IGVF raw matrix is not on disk (only the two DE CSVs, 2.1 MB) |
| re-extract IGVF TeloHAEC | the raw h5ad (214,449 cells, 1.47 GB) **is** on disk and already fully used; median 67 cells/target | **0** |
| re-extract `gse295214_h1` per-target cells | the on-disk object is a 15-dataset **pseudobulk database** (20,345 rows; `Replogle2022_K562_gwps` 6,608, `PertOrg` 4,363, `Pan2025_hESC` 2,354, …) with no cell counts in obs | 0 (and pluripotent) |

### (c) By pooling guides / replicates we currently split

| lever | measured | ≥ 400 somatic |
|---|---|---:|
| current: deepest **single** source per target | — | **214** |
| pool independent sources **within a line** (A549 ×3, K562 ×3 — the only lines with > 1 cell-counted source) | A549 15 shared targets, pooled median 195; K562 6 shared targets, pooled median 305.5 | **215** (+1 A549, +2 K562, overlapping the CD4 set) |
| pool **CD4's three culture conditions** into one row | 291 targets, pooled median **1,486** cells (min 28, max 12,395); 270 of them ≥ 400 on CD4 alone | **277** — but this conflates cell **state**: 276 treats `cd4_rest` / `cd4_stim8hr` / `cd4_stim48hr` as three separate contexts and the production builder pins all three as separate banks (`make_submission_pm_cd4.py:173–208`). It is a depth gain bought by destroying the only within-line state contrast the panel has |
| pool **guides** | 261 / 280 / 284 of the panel targets are **one dual-guide construct** in K562 / HCT116 / HEK293T (R267 §1); only 8–10 have two | **0** — there is nothing to pool |
| pool DLD-1's two low-MOI replicates | already pooled in `dld1_lowmoi`; no cell counts on disk in either form | 0 |
| pool Myllia's two A549 channels | already pooled by the compact-donor builder (equal-weight mean); 15 targets, pooled median 179 | **0** — still under 400 |
| pool across **lines** | that is the CONTEXT axis, which 276 measured (17 extra lines beat the two production donors on 1 of 3 readings) — it is not a depth lever | n/a |

**Summary of the cost analysis.** There is **no cheap download and no re-extraction that moves 214**. The three
things that are actually available are: (i) the CD4 cell-level files, worth taking for their *donor* replicate
structure rather than their depth; (ii) an IGVF cardiomyocyte re-extraction, worth 56 deep targets in one more
somatic context if the raw data can be found; (iii) waiting for X-Atlas/Pisces, which is the only prospect that
would put the whole panel at genome-scale depth in HepG2 and Jurkat — and which has been "coming soon" since
2026-03-17.

---

## 5. What this changes for the programme (reading, not a new measurement)

1. **The DEPTH axis is not empty — it is already spent.** 71 % of the panel has a ≥ 400-cell somatic measurement
   on this disk, in the challenge's own chemistry, and the pipeline already consumes it: the CD4 banks are the
   `cd4svd5` residual source (all three conditions pinned as separate banks, `make_submission_pm_cd4.py:173–208`),
   with rank-5 deflated away (READ-FIRST §2, §5). Nobody has to acquire depth to test
   the reliability law's promise on this panel; they have to **use the depth they own as a response estimate**,
   with its reliability, rather than as a whitened residual with its top direction removed. That is a concrete,
   fundable stage-2 question and it needs no data at all.
2. **But the depth is one context, and it is the wrong one.** 13 panel targets — not 183 — have ≥ 400 cells in two
   genuinely different somatic cell types. Any claim of the form "at 400 cells the conserved component is X" is,
   on this panel, a claim about primary CD4 T cells unless it is restricted to those 13. 276's reliability-corrected
   conservation map (median +0.228) rests on 4 sources / 6 pairs at exactly this limitation.
3. **The CD4 axis is 10,113 of 18,533 genes** (276 §1.1). A deep row on a narrow axis is not a deep row on the
   scored axis; that restriction should be priced into any DEPTH argument before it is priced into a build.
4. **"More cells" is not the same lever as "the Flex regime".** F1's depth-matched row (H1 Flex at 80 cells per
   half reads +0.114 with 0 % of targets ≤ 0, against Orion's +0.028 at 70–95 cells per half with 32–34 % ≤ 0)
   says a large part of the reference-quality gap is chemistry and recipient, not cell count. A DEPTH stage-2 that
   simply asks for more cells in the sources we already have is therefore priced against the wrong law; the
   inventory's own answer is that the Flex-regime rows we own are CD4's, and they are already in the build.
5. **The one instrument that is regime-exact touches 5 panel targets.** 273's h1-test-with-union-off reproduces the
   scorer's `pds_cosine` to 5.6e-17 and sits at exactly 400 cells / 19,997 UMI — but 5 of its 85 targets are on the
   2026 panel. Any depth conclusion carried from h1-test to A/B/C crosses both a context gap and a panel gap.
5. **A negative worth writing down**: the acquisition question the mandate was designed to precede — *should we buy
   or generate deeper data?* — is answered **no, not yet**, on measured grounds. Nothing purchasable today raises
   214; the only thing that would is unreleased; and the deepest data already on disk is not being read as a
   response estimate.

---

## 6. Files

* `RESULTS.md` (this file); `MANDATE.md`.
* **`work/depth278.csv`** — the deliverable: 3,169 rows × 13 columns, one per (panel target, source).
* `work/depth278_per_target.csv` — 300 rows: deepest somatic / somatic-or-derived / any measurement, its source and
  line, the number of measured sources, and the number of somatic lines at ≥ 200 cells.
* `work/depth278_summary.json` — threshold counts and the per-source summary.
* `work/gap_lt200_somatic.csv` — the 15 targets with no somatic measurement ≥ 200.
* `scan278_h5ad.py` → `work/h5ad_counts/*.json` (obs-only scans of 46 scPerturb h5ads + Nadig ×2, Altair, Myllia ×2,
  IGVF TeloHAEC), `work/h5ad_index.json`, `work/scan278_h5ad.log`.
* `build278.py` → the three CSVs above plus `work/{igvf_telohaec,igvf_cardiomyocyte,altair_a549,feng_*}_panel.csv`.
* `work/ka278.txt` (the six known-answer rows), `work/cd4_detail.txt`, `work/cost278.txt`, `work/lines278.txt`,
  `work/verify278.txt` (the re-checks of every number quoted above: CD4 medians, Orion/K562 in-bank vs all-rows
  medians, Feng control labels, the CRISPRa-excluded thresholds, and the 189 / 62 / 37 / 152 decomposition of the 214).
* Reused, not rebuilt: 267's `work/inv/<line>-targets.csv` (the three folds' per-target cells, guides, batches),
  276's `work/inventory276_sources.csv` (the source list and panel coverage), 268's R268 §1 (Broad-16 cells/KO),
  270's §4 (the on-disk dataset audit).
