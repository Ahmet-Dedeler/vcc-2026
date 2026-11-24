# 166 — donor acquisition: the panel-covering universe is exhausted, and we already hold all of it

**Date:** 2026-09-02 · deterministic (seed `20260902`) · **no `cell-eval2` was run, no submission
was built, nothing outside `experiments/166-donor-acquisition/` was written.** Total network
download: **~640 KB** (one 28 KB guide list, one TF list, two dataset READMEs, API JSON). Disk
added: **101 MB**, both of it the two donor tables built here. Peak RSS ~3.5 GB.

---

## Verdict

**NO-GO on the 1/k route. No unexploited panel-covering cell line exists in public data.**

I searched Hugging Face, GEO (258 series across two queries), Europe PMC (229 papers),
BioStudies/ArrayExpress, the Arc Virtual Cell Atlas, the IGVF data portal (6,902 measurement
sets, 374 construct libraries), and every h5ad already on this disk (64 files, 118 GB of Orion
included). Every genome-scale human Perturb-seq screen that exists is already in `data/donors/`
or already logged as rejected. The one outstanding candidate, **X-Atlas/Pisces, is an empty
repository** — 3 files, `used_storage: 0`, unchanged since 2026-03-25.

**Two new donor tables were nevertheless built and verified**, because they are the best
remaining *training* assets for experiment 161 even though neither can serve at inference:

| built | tissue | targets | panel | shared targets with held donors | verdict |
|---|---|---:|---:|---:|---|
| `donors/nadig_hepg2` | **liver** (hepatoblastoma) | 2,160 | **0/300** | 2,143 k562_gwps · 1,947/1,955 Orion · 2,086 Jurkat · 2,039 RPE1 | **training-only, high value** |
| `donors/frangieh_melanoma` | **melanocyte** (patient melanoma) | 235 | **5/300** | 202/203 Orion · 156 k562_gwps | training-only, marginal |

**Does this move the effective-donor count? No — not by one unit.** 159's `1/k` law counts lines
that carry the panel *at inference*. Both new lines carry 0 and 5 of the 300. Production stays
at 3.1 effective donors. What changed is the training substrate for the residual model: an 11th
line, a **new tissue**, and ~1,900–2,400 shared targets against four existing donors — an order
of magnitude denser than Jiang 2025's 171.

---

## 0. Known-answer rows — three of them, all reproduced

Rule 10. Every instrument here carries a configuration whose answer was already known.

| instrument | configuration | expected | measured | status |
|---|---|---|---|---|
| `census.py` | `k562_gwps` panel overlap | 272/300 (`build_orion_donor.py`) | **272** | OK |
| `census.py` | `rpe1`, `k562_essential`, `nadig_jurkat` | 0/300 (099 §1.1, brief) | **0, 0, 0** | OK |
| `independence.py` | `k562_gwps × orion_hct116` | 159 §4: +0.0321 / +0.0020 / +0.0321 | +0.0323±0.0006 / +0.0020±0.0001 / +0.0322±0.0006 | OK |
| `independence.py` | `k562_gwps × orion_hek293t` | +0.0289 / +0.0016 / +0.0286 | +0.0302±0.0008 / +0.0017±0.0001 / +0.0298±0.0009 | OK |
| `independence.py` | `orion_hct116 × orion_hek293t` | +0.0306 / +0.0016 / +0.0300 | +0.0297±0.0009 / +0.0019±0.0005 / +0.0291±0.0010 | OK |
| `verify.py` | permuted-label recomputation | must NOT match | cosine −0.077…+0.085 | OK |

`independence.py` **aborts** if the 159 rows do not come back; the first run did abort, which is
how the note in §6 below was found.

---

## 1. The structural fact: panel coverage is a library-size gate, and there is no middle

Computed over every donor table on disk (`census.py`, and the table below is from
`perts.txt` directly, so it is the *post-QC* library, not the paper's advertised one):

| donor | targets | valid symbols | panel | panel % |
|---|---:|---:|---:|---:|
| `orion_hek293t` | 17,977 | 17,078 | **294** | 98.0% |
| `orion_hct116` | 17,767 | 16,891 | **290** | 96.7% |
| `kolf2` | 273 | 273 | **273** | 91.0% |
| `k562_gwps` | 9,787 | 9,437 | **269** | 89.7% |
| `h1_2025` | 255 | 255 | 20 | 6.7% |
| `frangieh_melanoma` **[new]** | 235 | 206 | 5 | 1.7% |
| `k562_essential` | 2,021 | 1,821 | **0** | 0.0% |
| `rpe1` | 2,204 | 1,997 | **0** | 0.0% |
| `nadig_jurkat` | 2,270 | 2,059 | **0** | 0.0% |
| `nadig_hepg2` **[new]** | 2,160 | 1,967 | **0** | 0.0% |

Two regimes and nothing between them:

1. **Genome-wide (≥9,787 targets) → 90–98% of the panel.** Four screens, no exceptions.
2. **Common-essential libraries → exactly 0.0%.** Four independent screens (Replogle K562-ess,
   Replogle RPE1, Nadig Jurkat, Nadig HepG2), all zero, because **the VCC panel excludes
   common-essential genes by construction.** This is not a coverage accident that a bigger
   essential library would fix — it is a disjointness. Any future DepMap-essential screen, in
   any tissue, will read 0/300.
3. `kolf2` is the only member of a third class: a library *deliberately built against the panel*
   (273 targets, 273 on the panel). Only someone building for this competition produces that.

**So the ≥100 bar admits exactly two kinds of source: a genome-wide screen, or a panel-aligned
one.** That is what makes the search finite and it is why the negative result is decisive rather
than merely unlucky.

A corollary that prices a whole family in one number: **75 of the 300 panel targets are
transcription factors** (Lambert human TF list, 1,639 TFs; base rate would predict 26 — the
panel is 2.9× TF-enriched because it excludes common-essentials and ZNFs are largely
non-essential). So *any* TF-atlas, however complete, tops out at **75/300 = marginal band**, and
that ceiling applies to Joung 2023, the fibroblast CRISPRa atlas, and every future TF screen. No
TF atlas needs to be evaluated individually again.

---

## 2. Every source considered

Panel overlap computed by me from the actual file or its guide/target list — never from an
abstract. `[disk]` = already downloaded; `[new]` = examined for the first time in this repo.

### 2a. Genome-scale screens — the only class that can clear the bar

| source | lines | size | panel overlap | verdict |
|---|---|---|---:|---|
| Replogle 2022 K562-gwps `[disk]` | K562 | 8.8 GB | **269–272/300** | **SHIPPED** |
| X-Atlas/Orion HCT116 `[disk]` | HCT116 | 118 GB (pair) | **290/300** | **SHIPPED** |
| X-Atlas/Orion HEK293T `[disk]` | HEK293T | — | **294/300** | **SHIPPED** |
| **X-Atlas/Pisces** `[new]` | HepG2, iPSC, Jurkat rest+active, 16 contexts, 25.6M cells | **0 bytes** | n/a | **UNAVAILABLE.** HF repo `Xaira-Therapeutics/X-Atlas-Pisces` created 2026-03-25, last modified 2026-03-25, `used_storage: 0`, 3 files (`.gitattributes`, `LICENSE.md`, `README.md`), README says "**Coming Soon**". 5 months stale. And when it lands, its **HepG2 and iPSC arms are a 200-gene validation set**, not genome-wide — only the two Jurkat arms would be. CC-BY-NC-SA-4.0. |
| GSE314342 CD4+ T `[disk]` | primary CD4, 4 donors × 3 stim | 465 MB compact | **~284/300** | have; consumed as a post-whitening residual (062). The 4 donors are 4 individuals of **one cell type**, not 4 lines. |
| Feng 2025 iPSC `[disk]` | 24–34 HipSci iPSC lines | 63 MB | **182/300** | have; **neutral, +0.00007 (059)**. Per-line split is dead on depth: 104,094 assigned cells / 2,254 targets / 24 lines ≈ **2 cells per (line, target)** against our donors' 49–203. |
| Replogle RPE1, K562-ess `[disk]` | RPE1, K562 | 2.8 GB | **0/300** | training-only (common-essential) |
| Nadig Jurkat `[disk]` | Jurkat | 9.4 GB | **0/300** | training-only |
| **Nadig HepG2** `[disk, unbuilt]` | **HepG2 — liver** | 5.6 GB | **0/300** | **BUILT HERE** — §3 |
| GSE249595 Jurkat `[disk]` | Jurkat | — | 299/300 targets but **375-gene targeted readout, median 13 guides/cell** | rejected (107) |
| KOLF2 `[disk]` | KOLF2.1J iPSC | 30 MB | 273/300 | have; loses after whitening (057), anti-aligned with its own blob (159 §4) |
| VCC2025 H1 `[disk]` | H1 hESC | 32 GB | 20/300 (donor table) | have; **leaks 85/85 into h1-test** |
| GSE295214 / Pan `[disk]` | multi | 4.1 GB | 64/300 | consumed (117/118) |

### 2b. Below the bar — rejected with the number that killed it

| source | lines | targets | panel | why rejected |
|---|---|---:|---:|---|
| **McFaline-Figueroa sciPlexGxE** GSE225775 `[new]` | glioblastoma × drugs | **494 (kinome)**, 3,165 guides | **23/300** | Below the 40 floor. Computed from `GSM7056149_..._gRNA_sequences.txt.gz` (28 KB — metadata-first). Also genetic × chemical, so the control is drug-dose-confounded. **This closes the last open row of `AGENTS.md` §5.3.** |
| **Broad 16-line Perturb-seq** (bioRxiv 2026.08.24.746802) `[new]` | **16 diverse cancer lines** | **100** | ≤100/300, realistically ~0 (DepMap-essentiality-selected) | Library structurally too small. Data not yet posted (no GEO record under its authors). **Worth re-checking for 161** — 16 lines × 100 shared targets is a dense line×target matrix. |
| **Enhancer CRISPRi CRC** (Wellcome Open Res 2026) `[new]` | colorectal | **12,117 enhancers**, 35,139 guides | **0/300** | Targets enhancers, not genes. Zero gene-symbol targets by construction. |
| **Joung 2023 TF atlas** `[disk]` | hESC | 1,836 TF ORFs | ≤75/300 (§1) | CRISPRa **overexpression**, not knockdown; hESC duplicates `h1_2025`'s cell type. |
| Fibroblast TF CRISPRa (Nat Genet 2025) `[new]` | fibroblast, 2 types | 1,836 TFs | ≤75/300 (§1) | Same 75-TF ceiling; CRISPRa modality. Marginal band at absolute best. |
| ProPer-seq (bioRxiv 2026.02.04) `[new]` | primary T cells | 3,550 sgRNA×effector combos | targeted probe panel | Synthetic-biology readout, not a transcriptome-wide donor. |
| `arcinstitute/Perturb-Sapiens` `[new]` | whole-organism | cytokines + drugs | n/a | **AI-predicted**, not measured. Stack model output. |
| `arcinstitute/State-Parse-Filtered` `[new]` | PBMC | — | n/a | Cell-type/donor splits, not a perturbation screen. |
| Frangieh 2021 melanoma `[disk]` | melanoma | 248 | **5/300** | **BUILT HERE as training-only** — §3 |
| Xu 2023 `[disk]` | HEK293 | 203 | 5/300 | HEK293 ≈ HEK293T, already held → no independence |
| Tian 2021 CRISPRi `[disk]` | iPSC-neuron | 184 | 6/300 | 437 control cells; new tissue but far below bar |
| Sunshine 2023 `[disk]` | Calu-3 lung | 183 | 6/300 | below bar |
| Norman 2019, Dixit 2016, Adamson 2016, Schraivogel, Gasperini, Papalexi, Wessels, Datlinger, sciPlex, and 40 more `[disk]` | — | ≤1,729 | **0–6/300** | all below bar — full table in `work/census.json` |
| Jiang 2025 `[disk]` | 6 lines | — | **9/300** | closed as a donor (072). **But it is a real TRAINING asset**: 6 lines, 171+ shared non-panel targets. |
| TeloHAEC `[disk]` | endothelial | — | 65/300 | NO-GO as a direct donor (105/135) |
| HCASMC `[disk]` | smooth muscle | — | 27/300 | closed |
| A549 Altair / Myllia `[disk]` | A549 | — | 20 / 14 | closed (119–121) |

### 2c. Whole registries closed with one number each

| registry | swept | result |
|---|---|---|
| **IGVF data portal** `[new]` | 6,902 MeasurementSets, **329 Perturb-seq**, 374 ConstructLibrarySets, 106 guide libraries | **The largest gene-targeting guide library in the entire consortium is 24 genes.** Library selection criteria are 216 sequence variants, 41 accessible regions, 35 TSS, 15 cCREs — only 10 "genes" and 12 "TF genes". Max conceivable panel overlap **24/300**. This kills IGVF outright, including its 222 tempting WTC-11-derived **hepatocyte / cardiomyocyte / neuron / NSC** contexts. |
| **GEO** `[new]` | 142 human Perturb-seq/CROP-seq GSEs (2024+) and 116 CRISPRi/a single-cell GSEs (2025+), ranked by sample count | Largest are GSE314342 (CD4, have), GSE264667 (Nadig, have), GSE295214 (have), GSE249595 (rejected). **No genome-scale screen in a line we do not hold.** |
| **Europe PMC** `[new]` | 229 papers, Perturb-seq + atlas/cell-lines, 2025-06 → 2026-12 | Only new multi-line resources are Pisces (empty) and the Broad 16-line 100-gene set. |
| **Hugging Face** `[new]` | `perturb-seq`, `crispri`, `crispra`, `single cell perturbation`, `virtual cell`, + authors `Xaira-Therapeutics`, `arcinstitute`, `tahoebio`, `altoslabs`, `recursionpharma` | Only Orion (have) and Pisces (empty). |
| **Arc Virtual Cell Atlas** `[new]` | repo contents + 12 most recent commits | Contains only `scBaseCount`, `tahoe-100M`, `virtual-cell-challenge`. **No new perturbation atlas** since the 2026-01-12 scBaseCount release. |
| **BioStudies / ArrayExpress** `[new]` | 106 hits | PMC mirrors only; one small alternative-promoter set (E-MTAB-14567). |
| Zenodo | attempted | WAF-blocked (HTTP 403 on the API, twice, with and without UA). Jiang 2025 — this repo's one Zenodo asset — is already held. **Not swept; see §7.** |

### 2d. The strongest external cross-check I found

`altoslabs/primeflow-vcc-datasets` (HF, updated 2026-08-11) is an independent, well-resourced
group's curated list of *every dataset worth training a VCC model on*. Their inventory:

> `cd4t_perturbseq` (D1–D4 × Rest/Stim8hr/Stim48hr) · `feng24` · `jiang24` · **`mcfaline23_gxe`** ·
> `nadig24_hepg2` · `nadig24_jurkat` · `replogle20` · `replogle22_k562` · `replogle22_rpe1` ·
> `vcc_train`

**Their universe is identical to ours**, with one item we had never priced — `mcfaline23_gxe`,
which I priced at 23/300 above — and one we had never built — `nadig24_hepg2`, which I built.
Two independent searches converging on the same closed set is stronger evidence that the set is
closed than either search alone.

---

## 3. What was built, and its verification

### 3.1 `donors/nadig_hepg2` — liver, 79 MB

Built with the repo's own `scripts/build_donor_table.py`, where `nadig_hepg2` was **already
registered and had simply never been run**. `delta (2160, 9623)`, 4,976 control cells, median 49
cells/target (min 20, max 1,213), effect-norm median 14.43.

### 3.2 `donors/frangieh_melanoma` — melanocyte, 22 MB

`build_extra.py`, written here because the scPerturb mirror stores `.X` as **CSC**, which the
shipped row-block streamer cannot read; this walks gene blocks instead. Everything downstream of
the group sums is construction-identical to the shipped builder.
`delta (235, 23712)`, 15,361 control cells, median 175 cells/target, effect-norm median 8.39.

Restricted to `perturbation_2 == "Control"`. Frangieh pools three arms (Control 57,627 / IFNγ
87,590 / Co-culture 73,114); mixing them would put the stimulation response into the control
profile and into every effect vector.

### 3.3 Verification (`verify.py`) — and why the obvious check is not the one that matters

| | `nadig_hepg2` | `frangieh_melanoma` | shipped `k562_gwps` / `orion_hct116` / `orion_hek293t` |
|---|---|---|---|
| `Donor.load` | OK | OK | — |
| `sum(expm1(b_ctrl))` rel. err | 9.12e-09 | 3.94e-09 | 6.22e-09 / 7.64e-09 / 1.71e-09 |
| `sum(expm1(b_ctrl+delta))` rel. err, median | 4.87e-09 | 4.04e-09 | 5.36e-09 / 3.14e-09 / 2.61e-09 |
| — max | 2.78e-08 | 1.59e-08 | 3.16e-08 / 1.87e-08 / 1.70e-08 |
| **independent recomputation, cosine** | **1.000000** | **1.000000** | — |
| — max abs diff | 4.40e-07 | 4.70e-07 | — |
| **permuted-label control, cosine** | **−0.005 … +0.017** | **−0.077 … +0.085** | — |

**The mass invariant on its own is nearly vacuous here** and is reported only comparatively. The
builder writes `delta = prof − b_ctrl`, so `sum(expm1(b_ctrl + delta)) == TS` holds by
construction in exact arithmetic; the check only measures float32 round-trip error. What it *is*
good for is the comparison: both new tables sit in the same 2–9e-09 band as three tables we
already trust, so they are built to the same contract. AGENTS.md §11.1 logs two prior cases of a
control that was invariant by construction — this is not a third, because it is not the check
being relied on.

**The check that carries the weight** is the independent recomputation: five random targets per
table, pseudobulked straight from raw cells through `anndata`'s own backed row-slicing rather
than the builder's group-sum path. A wrong cell→group assignment, an off-by-one in the target
index, a mismatched duplicate-symbol collapse or the wrong control group all break this, and
**none of them break the mass check**. It returns cosine 1.000000 at float32 precision.

**And its negative control**: the same recomputation with target labels permuted collapses to
cosine ≈ 0.00–0.09. Asked the required question — *what would this read if the thing being
tested were absent?* — the answer is "≈ 0, not ≈ 1". The check can fail.

---

## 4. Independence

### 4.1 Basal profile (`independence.py`, 5,959 genes shared across 13 lines)

Raw cosines are ~0.83–0.94 for *every* pair because the housekeeping core dominates, so the
mean-removed number is the one that discriminates.

| line | ctx_A | ctx_B | ctx_C | (raw, ctx_A) |
|---|---:|---:|---:|---:|
| `h1_2025` | **+0.3702** | **+0.4674** | **+0.5397** | +0.9214 |
| `kolf2` | −0.1671 | −0.0812 | −0.1018 | — |
| `orion_hek293t` | −0.1992 | −0.3523 | −0.3317 | +0.8951 |
| `orion_hct116` | −0.4370 | −0.5091 | −0.3835 | +0.8802 |
| `k562_gwps` | −0.5055 | −0.6521 | −0.6106 | +0.8764 |
| **`nadig_hepg2`** | **−0.5516** | **−0.5185** | **−0.4269** | +0.8613 |
| **`frangieh_melanoma`** | **−0.5680** | **−0.4667** | **−0.4049** | +0.8315 |

Both new lines are as basally distant from the three challenge contexts as `k562_gwps` is —
i.e. as distant as our single most productive donor. `nadig_hepg2` is also near-orthogonal to
two shipped donors (`orion_hek293t` −0.024, `kolf2` −0.005), which is the profile of a genuinely
new view. Incidentally, **`h1_2025` is the only donor we hold that is positively aligned with
A/B/C** — a fact I have not seen recorded and which is not exploited anywhere.

### 4.2 Effect transfer, and the confound that had to be removed first

The raw table said `nadig_hepg2 × orion_hct116` reads resid-same **+0.1091** against
`orion_hct116 × orion_hek293t` at **+0.0291** — a 3.7× advantage that would make the new line a
far better transfer partner than our two shipped Orion lines are for each other.

**That comparison is confounded and must not be quoted.** HepG2 can only be scored on HepG2's
targets, which are the Nadig common-essential set; the Orion–Orion pair is scored on 3,000
targets drawn from ~17,700 mostly non-essential ones, and 159 §5.3 already notes that
common-essential genes transfer best. `confound.py` re-runs the Orion pair on **exactly HepG2's
1,917 targets**, and on a size-matched random draw of its own non-HepG2 targets:

| pair / target set | n | same | diff | RATIO | resid same | resid null |
|---|---:|---:|---:|---:|---:|---:|
| `orion_hct116 × orion_hek293t` on **HepG2's 1,917 targets** | 1,917 | +0.0760 | +0.0058 | 13.08 | **+0.0787** | +0.0035 |
| `orion_hct116 × orion_hek293t` on **1,917 non-HepG2 targets** | 1,917 | +0.0207 | +0.0023 | 8.93 | **+0.0198** | −0.0004 |
| `nadig_hepg2 × orion_hct116` | 1,917 | +0.1057 | +0.0226 | 4.67 | **+0.1073** | +0.0055 |
| `nadig_hepg2 × orion_hek293t` | 1,917 | +0.0661 | +0.0097 | 6.79 | **+0.0674** | +0.0004 |

> **Target-set effect ×3.98. Line effect ×1.36.**

Common-essential targets transfer ~4× better than non-essential ones between the *same* two
lines — that is the dominant term, and it independently confirms 159 §5.3 with a controlled
measurement. On identical targets HepG2 is a **modestly** better partner than the Orion pair is
for itself (+0.1073 vs +0.0787, ×1.36), and against HEK293T it is actually **worse** (+0.0674 vs
+0.0787). **HepG2 is a comparable-quality donor, not a superior one.** Without this control the
experiment would have shipped a 3.7× claim that is really 1.36×.

Full pair table (unrestricted target sets, for the record):

| line A | line B | n pert | same | diff | RATIO | resid same | resid null |
|---|---|---:|---:|---:|---:|---:|---:|
| `h1_2025` | `nadig_hepg2` | 61 | +0.1025 | +0.0090 | 11.37 | +0.0958 | +0.0003 |
| `nadig_hepg2` | `orion_hek293t` | 1,955 | +0.0662 | +0.0099 | 6.70 | +0.0678 | +0.0035 |
| `k562_essential` | `nadig_hepg2` | 1,828 | +0.1416 | +0.0239 | 5.94 | +0.1465 | +0.0043 |
| `frangieh_melanoma` | `orion_hek293t` | 203 | +0.0086 | +0.0015 | 5.59 | +0.0066 | +0.0020 |
| `k562_gwps` | `nadig_hepg2` | 2,143 | +0.1216 | +0.0235 | 5.18 | +0.1218 | +0.0047 |
| `nadig_hepg2` | `nadig_jurkat` | 2,086 | +0.1105 | +0.0217 | 5.09 | +0.1085 | +0.0017 |
| `nadig_hepg2` | `orion_hct116` | 1,947 | +0.1069 | +0.0233 | 4.59 | +0.1091 | +0.0011 |
| `nadig_hepg2` | `rpe1` | 2,039 | +0.2090 | +0.0874 | 2.39 | +0.1960 | +0.0121 |
| `frangieh_melanoma` | `k562_gwps` | 156 | +0.0117 | −0.0072 | n/a | +0.0181 | +0.0019 |
| `frangieh_melanoma` | `orion_hct116` | 202 | +0.0076 | −0.0008 | n/a | +0.0085 | +0.0005 |
| `frangieh_melanoma` | `rpe1` | 40 | +0.0427 | −0.0025 | n/a | +0.0488 | −0.0177 |
| `frangieh_melanoma` | `nadig_jurkat` | 43 | +0.0345 | −0.0041 | n/a | +0.0414 | +0.0002 |
| `frangieh_melanoma` | `nadig_hepg2` | 39 | +0.0315 | −0.0021 | n/a | +0.0371 | +0.0000 |

`nadig_hepg2` beats its permutation null by 10–100× on every pair — the transfer is
unambiguously real. **`frangieh_melanoma` is weak**: its decoy cosine goes *negative* in 6 of 7
pairs (making RATIO meaningless), its resid-same is +0.007 to +0.049 on 36–203 targets, and its
strongest pair (`rpe1`, +0.0488) has a null of −0.0177 on 40 targets. It is a Cas9 knockout
screen at MOI>1 in a co-culture-derived line; treat it as a low-weight training pair only.

---

## 5. A note on 159's line-pair matrix, found by the known-answer gate

The gate aborted on first run. Investigating rather than loosening it: 159's `donorpairs.py`
caps a pair at 3,000 targets by a random draw, and **the same draw also sets the gene axis**
(the sampled targets' own coordinates are dropped). So each published row is one realisation.
Measured over 12 independent draws, `orion_hct116 × orion_hek293t` gives same **0.0285 ± 0.0011**
(range 0.0264–0.0305) — and 159's published **0.0306** is at the very top of that range.
Evaluated on **all 17,698** shared targets instead, the same pair reads **+0.0231**, because
dropping 17,698 target coordinates leaves 20,886 genes instead of 35,584.

Nothing in 159's conclusions changes — the ordering and the ratios are robust — but its
line-pair values are single draws with a ±0.001–0.002 sd and a subsample-size dependence, and
should be quoted with that. This is consistent with the caveat `docs/CURRENT_STATE.md` already
carries about 159's magnitudes. My gate is now statistical (5 draws, tolerance
`max(3·sd, 0.0015)`), which still catches a real bug — those miss by ~10× — without pretending a
2-sd draw is a reproduction.

---

## 6. GO / NO-GO

**NO-GO on donor expansion as an inference lever. The 1/k route is closed on data availability.**

The public genome-scale human Perturb-seq universe is: Replogle (K562, RPE1), Nadig (Jurkat,
HepG2), Orion (HCT116, HEK293T), CD4-T GSE314342, Feng iPSC, KOLF2, VCC2025 H1, GSE249595,
GSE295214. **We hold every one.** Reaching 159's "half the ceiling at ~10 independent lines"
would need six more genome-wide screens in new tissues; six exist only inside X-Atlas/Pisces,
which is an empty repository, and even Pisces would contribute only its two Jurkat arms
genome-wide (its HepG2 and iPSC arms are a 200-gene validation set) — and we already hold a
Jurkat line.

**GO, narrowly, on `donors/nadig_hepg2` as a training input to experiment 161.** New tissue
(liver), 11th line, 1,917–2,143 shared targets with four held donors, verified against the
contract, transfer real at 10–100× its null. It moves the effective-donor count by **zero**,
because it carries 0/300 of the panel — exactly the asymmetry 160 §"the unlock" identified.

**Weak GO on `donors/frangieh_melanoma`** — melanocyte is a genuinely new lineage and the table
is verified, but the transfer signal is marginal and the decoy level is negative. Use it as a
low-weight pair or not at all.

The honest strategic read is the one 159 and 160 already reached, now closed on the data side
as well: **the remaining path is the cross-line residual model (161), not another donor.**

---

## 7. What I did NOT measure, and what is surface-limited

- **No metric was computed.** No `pds`, `mse`, `fid`, `jac`, `nmae`, `reach`; no `cell-eval2`;
  no submission; nothing added to production. Per the brief's rule 4, deciding whether to ship
  these tables is a separate experiment's job. Nothing here should be read as a score claim.
- **`fid`, `jac`, `reach` are untouched by construction** — rule 9. Every number here is a
  profile-space pseudobulk statistic, none has `n_real` in a denominator or union.
- **Zenodo was not swept.** Its API returned HTTP 403 (WAF) on both attempts. Jiang 2025, this
  repo's only Zenodo asset, is already held, and the GEO/EPMC/HF sweeps would surface any
  genome-scale screen that also has a Zenodo mirror — but Zenodo-only releases are a residual
  gap. It is the one stone I left unturned.
- **The Broad 16-line dataset has no located accession.** 100 genes puts it below the panel bar
  regardless, but its 16-line structure is valuable to 161 and its data location should be
  re-checked in a few weeks.
- **`nadig_hepg2` is the shallower of the two Nadig screens**: 4,976 control cells and median 49
  cells/target, against Jurkat's 12,013 and 86. Per-target reliability is correspondingly lower;
  159 §1.3 says a 150-cell reference is only 5–12% reliable on the target-specific component, so
  49 cells is thin. This is a real caveat for any model that weights HepG2 rows equally.
- **The `frangieh_melanoma` pair statistics rest on 36–203 shared targets**, near 159's MIN_PERT
  floor of 30. They are noisy and should not be quoted as levels.
- **Panel-overlap counts differ by one or two** between `census.py` (raw h5ad labels, e.g.
  `k562_gwps` 272) and the donor tables (post `min_cells ≥ 20` filter, 269). Both are correct
  for what they count.
- **X-Atlas/Pisces was checked on 2026-09-02.** It could land at any time and would be the
  single most valuable acquisition available. Worth a weekly `hf datasets info` poll — that is a
  10-second check with a very large payoff.

---

## 8. Files

| file | what it does |
|---|---|
| `census.py` | §1, §2 — metadata-only census of all 64 on-disk h5ads; reads `.obs` codes and `.var` index through h5py, never `.X` |
| `shared.py` | §2 — shared-target counts between candidate lines and the seven held donor tables |
| `build_extra.py` | §3.2 — CSC-aware donor builder with a stratum filter |
| `verify.py` | §3.3 — contract check, independent recomputation, permuted negative control |
| `independence.py` | §4 — basal cosines incl. challenge contexts, 159-style line-pair transfer, the known-answer gate |
| `confound.py` | §4.2 — the target-set vs line-quality control |
| `donors/nadig_hepg2/`, `donors/frangieh_melanoma/` | the two built tables (79 MB, 22 MB) |
| `work/*.json` | every number above |

```bash
cd experiments/166-donor-acquisition
uv run python census.py data/external                  # ~3 min, the survey
uv run python shared.py
uv run python ../../scripts/build_donor_table.py --donor nadig_hepg2 --outdir donors
uv run python build_extra.py --source frangieh_melanoma
uv run python verify.py                                # known-answer + negative control
uv run python independence.py                          # aborts if 159 does not reproduce
uv run python confound.py
```
