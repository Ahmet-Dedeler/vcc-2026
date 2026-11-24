# RESERVE R3 — the single-cell perturbation datasets on disk whose outcomes no agent has inspected

**Declared 2026-09-06 by experiment 292 (metadata only).** Nothing in this file was produced by reading an
expression value. Every number is a cell count, a label, a set intersection or a sha256. `X`, `layers`, `obsm`,
`varm` and `obsp` were never opened; the h5ad files were read with `h5py` through `obs`, `var` and `uns` only.

**RESERVE R3 DECLARED: 8 datasets, 8 human genetic-perturbation sets, 174 perturbations in the held map
(pooled; 152 in the 7-file core), 1 new line; NOTHING COMPUTED BEYOND METADATA.**

---

## 0. What a reserve is here, and why this one is different from R1 / R2

`docs/RESEARCH-PROGRAMME-2026-09.md` rule 2: a reserve is data whose **outcomes** no agent has inspected, scored
only for a preregistered claim, once per claim, every access logged. Both earlier reserves are spent — R1
(`orion_hek293t` half-B, claim C1, and relabelled 2026-09-06 13:40Z as split-half robustness because the half-B
cells sat inside the HEK293T truth 276 had already evaluated against) and R2 (four Broad-16 lines, claim C2).
291's two surfaces (Broad-16's 12 working lines, GWCD4i's published K562 comparison) were never reserves and have
now been looked at once with 290's features.

R3 is stronger than R1 in exactly the way R1 failed: **none of its cells is inside any bank, any map row, any
donor table or any fold used by the programme**, and no derived object of theirs exists on disk. What HAS been
read is their `obs` label vocabulary and cell counts, by five census scripts that never open `X`
(`E166/shared.py`, `E218/census166_copy.py`, `E225/census225.py`, E270's obs scan, `E278/scan278_h5ad.py`) and by
this experiment. Metadata reads leave the outcomes untouched; that is the property the access rule protects.

---

## 1. Known-answer rows (all before any reserve number)

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** every scPerturb file's `n_cells` and `n_labels` reproduce 278's independent obs-only scan | equality on all files with a 278 record | **50/50 files, exact** (`work/scan292_ka.json`) | PASS |
| **KA-2** 279's held gene list and 290's feature-store gene list are one object | same length, same order | **17,106 genes, element-wise identical** — so the mandate's overlaps (a) and (c) are ONE number, reported once | PASS |
| **KA-3** the 300 panel targets read from `conservation276.csv` equal `data/challenge/pert_counts.csv` | equality | **300, equal** | PASS |
| **KA-4** the sealed HCASMC asset is byte-identical to what 126 acquired | `b87452c2…` (`docs/DEF-RUNBOOK.md`:296) | **`b87452c2db1f508f5013bc15dff8088828ab931b35492d77ed6aea569fb5d9b9`**, 1,531,381,676 B — sha256 only, no decompression, no header read | PASS |

---

## 2. The reserve

Eight h5ad files, 1,575,574,586 B, 427,073 cells, 85,859 control cells. `sha256` in full in
`work/sha256.txt`; `path` + digest below.

### 2.1 Core (7 files) — human, CRISPR, transcriptome-wide, ≥ 50 cells per gene, own controls

| # | file (`data/external/scperturb/…`) | sha256 | bytes | line / cell type | modality | genes ≥ 50 cells (in held map / in panel) | cells per gene (median, q10) | control cells | guides per gene | replicate column |
|---|---|---|---:|---|---|---|---|---:|---:|---|
| 1 | `AdamsonWeissman2016_GSM2406681_10X010.h5ad` | `e70fcd49808cab8d724de8d5a332940911206e1c8ef44cc7b568d048ed795c85` | 471,286,951 | K562 (repeat line) | CRISPRi | **82** (80 / 2: `STT3A`, `SYVN1`) | 583, 253 | 10,007 | 1 | — |
| 2 | `DatlingerBock2017.h5ad` | `2ee8b95ef1b2567ce5cbf83fc3a20bf4318f32dc007e1228f064fe428125f7da` | 39,085,166 | Jurkat (repeat line) | CRISPR-KO (CROP-seq) | **28** (28 / 0) | 150, 74 | 1,320 | **3** | `replicate` (6) |
| 3 | `DatlingerBock2021.h5ad` | `9f8c35058d71442a904391dabf6574bbb24ae31ed00cef5d4e4ed498828f3592` | 33,600,780 | Jurkat (repeat line) | CRISPR-KO (scifi-RNA-seq) | **20** (20 / 0) | **1,702**, 1,062 | 4,497 | **2** | `sample` (384) |
| 4 | `DixitRegev2016.h5ad` | `a448af6aa250a6cdca472d26bcba5aec857381d8aeb03024b9358aebc2d6d23e` | 309,282,071 | K562 (repeat line) | CRISPRi (Perturb-seq) | **10** (8 / 2: `E2F4`, `ELK1`) | **2,435**, 1,928 | 10,247 | **2.5** | `guide_id` (1,729) |
| 5 | `PapalexiSatija2021_eccite_RNA.h5ad` | `03e9602d124c261281936717d5795445805e6e704a52b0310bcde94225aa0929` | 147,215,278 | THP-1 (repeat line) | CRISPR-KO (ECCITE) | **23** (22 / 1: `IFNGR2`) | 662, 294 | 2,386 | **4** | `guide_id` (107) |
| 6 | `PapalexiSatija2021_eccite_arrayed_RNA.h5ad` | `7f8b06d92955940c453b7bec5dec879d64c1b2e2667aa5b6f7bb4bbc6665d976` | 52,339,395 | THP-1 (repeat line) | CRISPR-KO, **arrayed** (one gene per well) | **7** (7 / 0) | 726, 288 | 2,009 | 1 | `guide_id` (14) |
| 7 | `ShifrutMarson2018.h5ad` | `122c31ae44177f79b2d0441d8d10bf47028af0d6e536543e8c734048f8eecf29` | 253,843,410 | **primary human CD8⁺ T cells — NEW context** | Cas9 RNP KO | **18** (17 / 1: `STAT6`) | 566, 404 | **39,813** | **2** | `replicate` (2 donors) |

### 2.2 Restricted (1 file) — usable, but not study-independent

| # | file | sha256 | bytes | line | modality | genes ≥ 50 (held / panel) | cells per gene | controls | restriction |
|---|---|---|---:|---|---|---|---|---:|---|
| 8 | `TianKampmann2019_day7neuron.h5ad` | `e5883ffef02c0da955b8d6f1dd3c92ee1dc4d73f43d6a1ed9562444097fca5df` | 268,921,535 | iPSC-neuron day 7 (repeat line) | CRISPRi | **24** (23 / 1: `RFK`) | 1,314, 299 | 15,580 | **same study** as the map bank `tian2019_ipsc` (Tian 2019). It adds a cell STATE, not an independent study; a claim using it must say so. |

### 2.3 What each contributes

- **Adamson 10X010** is the reserve's volume: 82 genes at 583 cells, an independent 2016 CRISPRi study in a line
  the map already holds four times over. It is the single best row for "does the map's `C_held` predict how well a
  gene's response transfers into an unseen K562 study?"
- **Datlinger 2017 / 2021** are the reserve's construct axis: 3 and 2 guides per gene on the same TCR-signalling
  library, in the same line, from the same laboratory five years apart — and **all 20 of 2021's genes are inside
  2017's 28** (measured), so the pair is an exact same-line, same-library, different-chemistry replication on 20
  genes. Programme point 2 (a second construct per target beats more cells of the first) has no on-disk test bed
  outside CD4; this is one.
- **Dixit 2016** is the reserve's depth: 2,435 cells per gene and 2.5 guides, 10 genes.
- **Papalexi ECCITE + arrayed** are the reserve's guide-assignment control: 4 guides per gene in the pooled screen
  and **all 7** arrayed genes inside the pooled 23 (measured), run one gene per well with no guide-calling step,
  in the same paper. Whether a response
  agrees between the pooled and the arrayed measurement bounds how much of F11's donor × guide effect is
  assignment noise. THP-1 is a map line, but through `wessels_thp1`, which is **Cas13 knockdown** — so Papalexi is
  a new study AND a new modality in that line.
- **Shifrut 2018** is the only NEW context: primary human CD8⁺ T cells from two healthy donors, Cas9 RNP
  knockout, 18 genes at 566 cells with 39,813 control cells. The programme's other primary-T data is GWCD4i
  (CD4, CRISPRi, pooled donors) and is used. This is a different cell type, a different perturbation chemistry
  and a different laboratory.
- **Tian day-7 neuron** adds the differentiated state of a map bank's own iPSC screen.

### 2.4 Comparability with the map, measured

The map's response axis is 279's 2,739 genes. Measured `var` overlap (a `var`-only read, `work/varaxis292.json`):
Adamson 2,732 · Datlinger2017 2,732 · Datlinger2021 2,636 · Dixit 2,731 · Papalexi 2,725 · Papalexi arrayed
2,723 · Shifrut **2,739** · Tian day7 **2,739**. Every reserve file can carry a response vector on ≥ 96 % of the
axis the map is written on, so no re-definition of the response space is needed to score a claim.

---

## 3. The frozen access rule

> **R3 access rule (frozen 2026-09-06).** The eight files listed in §2 are a reserve. **Outcomes — responses,
> pseudobulks, differential expression, cosines, conservation values, agreement statistics, reliabilities, any
> feature derived from an expression value — may be computed on them exactly ONCE, for a preregistered claim
> written into `docs/CONFIRMATION-LEDGER.md` before the first byte of `X` is read, by an experiment the
> coordinator names in that preregistration.** Metadata reads (`obs`, `var`, `uns`, file digests, cell counts,
> label vocabularies) are free and unlimited. **No exploratory look**: no "quick check", no pilot, no
> hyper-parameter chosen on the reserve, no second pass after seeing the first. The preregistration must name the
> statistic, the pass rule, the baselines, the controls and the known-answer rows; the result is appended to the
> ledger whatever it is, and the file is then SPENT and never used again for a claim.
> A file may be spent alone: spending Adamson does not spend Shifrut. Each row of §2.1/§2.2 is spent
> independently, and the ledger records which.

Two riders, from what went wrong with R1 and C2a:

1. **Independence is a property of the cells, not of the moment the test was frozen** (R1's relabelling). A claim
   on R3 must state that no derived object of the file exists on disk — true today, and checkable by the
   `evidence` column of `work/inventory292.csv`.
2. **Depth-condition the pass rule** (programme rule 9). The reserve's cells-per-gene span 150 (Datlinger 2017) to
   2,435 (Dixit); its q10 spans 74 to 1,928. A threshold frozen without stratifying or covarying on cells (and on
   UMIs per perturbation, which is NOT in this metadata and must be measured inside the single access) will fail
   the shallow files for the reason C2a failed, and the source papers' own reliability statistics should be read
   first.

---

## 4. What the reserve CANNOT confirm

- **It cannot measure cross-line conservation of its own.** Only **3 genes** (`EGR1`, `IRF1`, `SOCS1`) appear at
  ≥ 50 cells in two distinct reserve studies, and **0** in three. The reserve's usable outcome is per-gene
  *agreement with the map's response for that gene measured in other lines* — the shape of 291's C3 / C4 — not a
  conservation score computed inside the reserve.
- **Modality differs from the map's dominant CRISPRi.** Four of the eight files are Cas9 knockout (Datlinger ×2,
  Papalexi ×2, Shifrut). A knockout and a knockdown of the same gene are not the same perturbation; a claim must
  either restrict to the CRISPRi files (Adamson, Dixit, Tian day7 — **115** genes ≥ 50 cells, 110 in the held map,
  5 on the panel) or the Cas9-KO files (**69** genes, 67 in the held map, 2 on the panel), or carry modality as a
  declared covariate over the pooled 174.
- **Tens of targets per set → pooled claims only.** Per-file n is 7–82. A per-file interval on 7–28 genes is
  uninformative (291's C3 failed condition 1 at n = 82 with a ±0.22 interval). The reserve supports ONE pooled
  claim over 152 (core) or 174 (with Tian) genes, with the file as a block in the bootstrap — not eight claims.
- **Panel coverage is 6 targets** (`STT3A`, `SYVN1`, `E2F4`, `ELK1`, `IFNGR2`, `STAT6`; 7 with `RFK`). The reserve
  says nothing directly about the challenge panel and must never be used to select or tune anything competition-facing.
- **It cannot test F10 / F1's reliability laws as designed.** Only Datlinger 2017/2021, Dixit, Papalexi and
  Shifrut carry a second construct per gene, and only Shifrut carries a biological-donor unit (2 donors) — a
  split-half reliability law over cells is computable, but a donor law is not.
- **It cannot confirm anything about the incumbent's board score.** Its lines are not A/B/C and its targets are
  not the panel; the competition track gains nothing from spending it.
- **It is not a substitute for the sealed HCASMC asset**, which remains sealed under its own older rule
  (`docs/DEF-RUNBOOK.md` §5.-3): 27 panel targets, competition-facing, one shot after the D/E/F freeze. R3 does
  not authorise touching it and this experiment did not.

---

## 5. Inventory over every dataset on disk

Full table: `work/inventory292.csv` (104 data rows: 50 scPerturb h5ad files, 43 rows of 276's source inventory,
11 aggregate rows for the remaining collections). **USED 62 · UNUSED 42.** "USED" means an experiment computed an
outcome from it (a bank, a map row, a feature, a score); the `evidence` column carries the bank name and the file
and line of the code that reads it. "UNUSED" means no code on disk reads it for anything but a census.

### 5.1 The 28 conservation-map banks and the rest of the used set

The 28 banks are `E279/work/build279.json` `sources`; 276's inventory adds `tian_ipsc_neuron_crispra`,
`norman_k562`, `dld1_low1/2`, the three `cd4_*` conditions, the two `jiang_*` aggregates, the three `cmap_*`
tables, `kolf2` and `feng_ipsc`. Twelve scPerturb files feed these banks:

| scPerturb file | bank | reader |
|---|---|---|
| `ReplogleWeissman2022_K562_gwps` / `_K562_essential` / `_rpe1` | `k562_gwps`, `k562_essential`, `rpe1` | `scripts/build_donor_table.py:40,46,52` |
| `FrangiehIzar2021_RNA` | `frangieh_melanoma` | `E166/build_extra.py` |
| `SunshineHein2023`, `TianKampmann2021_CRISPRi`, `TianKampmann2021_CRISPRa`, `TianKampmann2019_iPSC`, `XuCao2023`, `NormanWeissman2019_filtered`, `GasperiniShendure2019_lowMOI`, `WesselsSatija2023` | `sunshine_calu3`, `tian_ipsc_neuron`, `tian_ipsc_neuron_crispra`, `tian2019_ipsc`, `xucao_hek293`, `norman_k562`, `gasperini_k562`, `wessels_thp1` | `E225/build225.py:60-67` |

### 5.2 The 38 UNUSED scPerturb files, and why 30 of them are not in the reserve

| reason | files | note |
|---|---|---|
| **in the reserve** | 8 | §2 |
| mouse | `CuiHacohen2023`, `LaraAstiasoHuntly2023_{exvivo,invivo,leukemia}`, `LiangWang2023`, `SantinhaPlatt2023`, `SchiebingerLander2019_{GSE106340,GSE115943}`, `WeinrebKlein2020`, `GehringPachter2019` | 10. `CuiHacohen2023` is a symbol-matching false positive: its 31 "gene" labels are cytokine TREATMENTS (`IL2`, `IFNG`, …), not genetic perturbations |
| chemical | `AissaBenevolenskaya2021`, `ChangYe2021`, `LotfollahiTheis2023`, `SrivatsanTrapnell2020_sciplex{2,3,4}`, `ZhaoSims2021` | 7 |
| not a transcriptome | `FrangiehIzar2021_protein` (24 proteins, and the SAME cells as a used bank), `PapalexiSatija2021_eccite_protein` (4), `PapalexiSatija2021_eccite_arrayed_protein` (4) | 3 |
| enhancer / non-gene targets | `GasperiniShendure2019_atscale`, `GasperiniShendure2019_highMOI` (multi-guide combination labels, high MOI), `XieHon2017` (chromosome coordinates) | 3 |
| targeted readout | `SchraivogelSteinmetz2020_TAP_SCREEN__chromosome_{8,11}_screen` | 2 — K562 CRISPRi with real gene targets (9 and 10 at ≥ 50 cells, 4 guides each), but 4,191 / 3,185 measured genes cover only **28 / 30** of the map's 2,739 response-axis genes. No comparable response vector; excluded on a measured number, not a guess |
| too small / shallow | `AdamsonWeissman2016_GSM2406675_10X001` (7 genes, pilot lane of the reserve's own study), `AdamsonWeissman2016_GSM2406677_10X005` (median 3 cells per gene), `McFarlandTsherniak2020` (1 gene, `GPX4`; the rest is drugs across 209 lines) | 3 |
| identity unresolvable on disk | `JoungZhang2023_atlas`, `JoungZhang2023_combinatorial` | 2 — see §5.3 |

### 5.3 The one recoverable exclusion

`JoungZhang2023_atlas.h5ad` (5,804,454,475 B, sha256
`8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22`) is **1,145,823 hESC
cells over 3,369 TF ORF perturbations** — by far the largest untouched human genetic-perturbation set on this
disk. Its `obs.perturbation` is a bare numeric ORF id and no column of `obs`, `var` or `uns` maps it to a symbol;
the mapping lives in the paper's supplementary table, which is not on disk and which this mandate may not
download. Two further caveats before anyone spends effort on it: the modality is **ORF overexpression**, whose
sign is not the map's knockdown sign, and hESC repeats the map's `H1` line. **Recommendation: a coordinator may
authorise a ≤ 100 KB metadata download to resolve the ids; that is a metadata read and does not spend the
reserve. Nothing else about the file may be read until it is declared.**

### 5.4 The non-scPerturb collections

`hcasmc-igvf` — **UNUSED and already sealed** by 126 under `docs/DEF-RUNBOOK.md` §5.-3 (27 panel targets, one shot
after the D/E/F freeze). This experiment computed its sha256 only (KA-4) and did not decompress it, read its
header, or derive any value-dependent statistic. It is **not** part of R3 and R3's rule does not govern it.
`hspc_crispri500` — UNUSED, mouse HSPC DE table, ineligible. `gse249595` — a guide FASTA, one supplementary xlsx
and a GEO file listing; **no expression matrix on disk**, nothing to score. Everything else under
`data/external/` is either a used perturbation source (`orion`, `cmap`, `broad16line`, `tahoe`, `nadig2025`,
`perturbmatch`, `dld1-lowmoi`, `altair-a549`, `myllia-a549`, `jiang2025`, `kolf2`, `feng-ipsc`/`feng2025`,
`gse295214`, `igvf-*`, `cd4t`, `cd4_cell2026`, `vcc2025`, `replogle2022` placeholders) or annotation, model
checkpoints and code (`depmap26q1`, `ccle_broad_2025`, `regulon`, `omnipath-dorothea`, `slim`, `sigrecon`,
`state-pretrained`, `state-vcc`, `xcell`, `perturbation-decomposition`).

---

## 6. Suggested first claim (NOT preregistered here; the coordinator writes it)

The reserve's shape fits one pooled, out-of-map generalisation test of the two live claims at once, over the 174
genes (152 core-only) it holds at ≥ 50 cells:

> For each reserve gene, compute its response in the reserve file and its agreement with the map's consensus
> response for that gene measured in the 28 banks' other lines. Test (i) that 279's `C_held` predicts that
> agreement (F8 off-map, 291's C4 shape) and (ii) that 290's `gen_pc20` predicts it with the SIGN F13 claims
> (broad, spectrally central responses more conserved), with the annotation ridge `M_ann`, magnitude, the gene's
> own DE count and a shuffled-gene null as the baselines, the file as the bootstrap block, cells-per-gene as a
> declared covariate, and modality (Cas9 KO vs CRISPRi) as a declared stratum.

That is one preregistration, one access, all eight files, and it is the only use of R3 that answers a question
the write-up actually asks. Anything smaller wastes it.

---

## 7. Process notes and one breach

- No build, no `vcc`, no `cell-eval2`, no download, nothing written under `data/` or `submissions/`. All writes
  are inside `experiments/292-reserve-declaration/` except the R3 entry appended to `docs/CONFIRMATION-LEDGER.md`.
  Sibling 285 was neither read nor written.
- **Peak RSS breach, recorded.** The mandate capped this run at ~6 GB; `work/scan292.py` peaked at **11.5 GiB**
  (`work/scan292_ka.json`), because it materialised whole `obs` string columns as NumPy object arrays for the
  two largest files (`JoungZhang2023_atlas`, 1,145,823 cells; `GasperiniShendure2019_atscale`, 203,496 distinct
  labels). It completed, `memory_pressure` stayed healthy (752,283 free pages after) and no sibling job failed,
  but the cap was exceeded and the script must be run per-file in a subprocess if it is ever re-run. Everything
  after that step (the label→gene resolution, the inventory, the reserve totals) reads 278's and 292's JSON
  artefacts only and peaks under 0.1 GiB.
- `work/parse292.py` carries the FROZEN per-dataset label grammars (`Tcrlibrary_<GENE>_<n>`, `ES.sg<n>.<GENE>`,
  `p-sg<GENE>-<n>`, `<GENE>g<n>`, `<GENE>_p<plasmid>`, …). They were written after reading label vocabularies
  only. A resolved symbol counts as a gene perturbation iff it is in 279's 17,106 held genes ∪ the 300 panel.
- Files: `work/scan292.py`, `work/scan292.json`, `work/scan292_ka.json`, `work/parse292.py`, `work/perts292.json`,
  `work/inventory292.py`, `work/inventory292.csv`, `work/reserve292.json`, `work/sha256.txt`,
  `work/sha256_scperturb.txt`, `work/sha256_other.txt`, `work/varaxis292.json`, `work/varoverlap292.json`,
  `work/reserve_totals292.txt`.
