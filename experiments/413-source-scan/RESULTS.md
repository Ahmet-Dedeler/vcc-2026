# 413 — Source scan, 2026-09-22: no new panel-covering cancer LINE is public; the best new item is a second genome-wide K562 CRISPRi screen (VIPerturb-seq, 284/300 at ≥20 cells), and the only pending line that matters is still X-Atlas/Pisces Jurkat (empty repo)

**Date:** 2026-09-22, ~19:30–20:35 UTC. **Mandate:** `MANDATE.md` (PI, 19:35Z). **Hygiene:** no builder, no `vcc`, no
`cell-eval2`, nothing under `submissions/`, nothing retained outside `experiments/413-source-scan/` (one transient exception below). `uv run python`
throughout. Peak RSS < 1 GB (obs-only h5py reads; no `.X` was ever loaded). Transient downloads 2.6 GB
(Perturb-ME MOI1 h5ad 1.43 GB and DeSCOPE PC9 h5ad 1.13 GB, each read for `.obs` only and **deleted**; the A549
head-to-head DE table, 732 MB, was streamed through `awk` and never stored). Retained: **75 MB** in `work/`
(largest: GSE345058 cell annotations 40 MB, HGNC table 17 MB). Nothing > 2 GB was downloaded. One deviation from the
write-only-here rule: a 34 KB 10x web page was saved once to `/tmp/tenx.html` and deleted immediately. The RDS reader
(`rdata`) ran as an ephemeral `uv run --with` overlay; the shared venv is unchanged.

---

## 0. Verdict in ten lines

1. **Known answer reproduces exactly:** `perts.txt ∩ panel` = **269 / 290 / 294 / 0** for
   `k562_gwps / orion_hct116 / orion_hek293t / nadig_jurkat`, both literally and after HGNC alias resolution (§1).
2. **No public dataset released since 2026-08-15 adds a panel-covering new cell LINE.** The 1/k route is still
   closed on data availability (166's verdict stands for a fourth time, after 225, 270 and 306).
3. **Best new item (missed by 225, which only searched GEO): VIPerturb-seq** (Satija lab; Zenodo 18460279,
   CC BY 4.0, released 2026-02-14). Genome-wide K562 CRISPRi screen read out with fixed-cell probes (Flex-style,
   the challenge's chemistry): **300/300 panel in the manifest, 284/300 at ≥ 20 cells** (median 47 cells per panel
   target), 906,837 cells incl. 7,949 `NO-TARGET` controls. It is **K562 again**, so it adds a *line* only for the **30** panel targets absent from our
   K562-GWPS donor. Five of those are among the 10 panel targets that currently have ≤ 2 held lines:
   **ABCD1, EPHB2, KIF21B, NICN1, SPATA20**. 17.3 GB of Seurat RDS; every file is > 2 GB (none downloaded).
4. **Best new cancer-line item: GSE345058** (CZ Biohub multimodal atlas, GEO public 2026-08-26). A549 lung
   adenocarcinoma, dox-Cas9 knockout, 1,000 genes, 606,075 cells, **median 565 cells/target**. It is deep, but only
   **21/300** panel targets, because 550 of its 1,000 genes are in the common-essential (Nadig) library. It is a
   *training* asset (877–929 targets shared with each held donor), not a donor.
5. **Perturb-ME** (A375 melanoma, genome-scale knockout, figshare 2026-08-23): library 300/300, but the perturbed
   cells come only from HLA-sorted gates. **24/300 panel targets reach ≥ 20 cells** (median 11). Not a donor.
6. **Criterion (b) is empty in public data.** HeLa, CAL33/HNSCC and Jurkat have **no** transcriptomic screen of
   ≥ 30 panel targets anywhere: GPSAdb bulk knockdowns give HeLa **4**, Jurkat **1** and HNSCC **0** panel targets. A second
   independent compendium (the LLM-built GEO atlas, 5,350 KD/KO signatures) gives **4 / 1 / 1**. The only
   Jurkat perturbation screens are Nadig (0/300, held), Datlinger (0/300, known) and Song 2025 (GSE247599,
   10 genes, **0/300**, new here). CMap's BICR6 (HNSCC) × 82 panel targets is on disk and was closed by 176.
7. **X-Atlas/Pisces is still an empty repo** (`usedStorage 0`, 3 files, last modified 2026-03-25, checked today).
   Xaira's own `docs/data.md` confirms the planned public subset includes the **full Jurkat Resting and Jurkat
   Active screens (10,872 / 10,878 perturbations, 2.8 M cells each)**. It remains the single most valuable pending
   item (see §5).
8. **Illumina Billion Cell Atlas** (CRISPR on/off × 20,000 genes × 200–250 lines) is **commercial alliance only**
   (AstraZeneca, Lilly, Merck, Formation Bio and two more; July 2026 expansion). No public release. The CZI Billion Cells Project's only
   public CRISPR release is the CD4 atlas (held). Tahoe's genetic arm (300 M cells via Parse) is unreleased; the
   Tahoe/Arc/Biohub 120 M-cell project is chemical.
9. **Registries closed with one number each** (§3): Arc Virtual Cell Atlas — no new data (last commit 2026-08-28,
   README). scPerturb — v1.4 of 2024-08-26 is still current. LINCS — no post-2020 build. Hugging Face — 1,648 matching datasets, zero
   new genetic screens in cancer lines. EBI Perturbation Catalogue inventory (merged 2026-09-18, 30 datasets) —
   everything held except Song 2025 (0/300) and a primary-Treg set.
10. **Recommendation (§6):** build nothing new into the donor bank for the 1/k law. Worth doing, in order:
    (1) extract VIPerturb-seq K562 as a **second K562 replicate plus 30-target gap-filler** (cloud or a one-off
    ~11 GB download); (2) keep the Pisces poll running and pre-register a Jurkat build for the day it lands;
    (3) optionally add GSE345058 A549 (40 MB metadata now; 11–25 GB counts) as a deep epithelial *training*
    line for the cross-line residual model. Everything else is below the bar.

---

## 1. Known-answer row (required first)

`work/known_answer.py` and `work/known_answer2.py` (literal and HGNC-resolved, through the same `overlap.py`
used for every candidate below):

| donor | perturbations | panel (literal) | panel (HGNC-resolved) | expected | status |
|---|---:|---:|---:|---:|---|
| `k562_gwps` | 9,787 | **269** | 269 | 269 | PASS |
| `orion_hct116` | 17,767 | **290** | 290 | 290 | PASS |
| `orion_hek293t` | 17,977 | **294** | 294 | 294 | PASS |
| `nadig_jurkat` | 2,270 | **0** | 0 | 0 | PASS |

The first version of the resolver mapped only the source side and read 289/293 for Orion, one short. The cause is a
panel symbol that HGNC has since renamed (`TMEM104` → `SLC38A12`). The resolver now maps both sides, and all four
rows reproduce. Every "panel" number below is the **literal** intersection. The resolved number is identical in
every row of §2.

---

## 2. The table — every candidate examined, with exact panel overlap

"Library" = targets in the downloadable perturbation list. "≥ 20 cells" = targets with at least 20 assigned cells,
the donor builder's threshold. `[new]` = not in any prior sweep (045/166/225/270/306). Panel counts were computed
here from the files; none is quoted from a paper.

| # | dataset | URL | released | line(s) | modality | #perts | **panel: library / ≥ 20 cells** | size | license / access |
|---|---|---|---|---|---|---:|---|---|---|
| 1 | **VIPerturb-seq genome-wide** (Bradu … Satija, bioRxiv 2026.02.12.705613) `[new — 225 found no GEO record; data is on Zenodo]` | zenodo.org/records/18460279 | 2026-02-14 | **K562** | CRISPRi (GuEST-List, 3 sg/gene), fixed-cell split-probe readout | 18,885 targeted genes in manifest (+ `NO-TARGET`, 7,949 cells); 17,490 at ≥ 20 cells | **300 / 284** (≥ 50 cells: 134; in the 3.6 GB "filtered" object: 124). Median 47 cells per panel target | 17.3 GB (binA/B/C 3.6/3.8/2.9 GB, filtered 3.6 GB, vimentin-sorted 1.9 GB, cell-line-mixing pilot 1.6 GB) | CC BY 4.0, open |
| 1b | VIPerturb "custodian export" (third party, Z. Peng) `[new]` | figshare 33719593 / 33732421 | 2026-09-14 | K562 | 20-cell means from the *filtered* object | 6,168 | **118** (v3 full-depth subset: 835 targets, **29**) | 470 MB / 421 MB | CC BY 4.0. **Its control is a pseudo-control drawn from perturbed cells**: its control pattern matched only 66 cells, although the filtered object carries all 7,949 `NO-TARGET` cells (manifest flag). It is therefore not a clean donor. It is also evidence that someone else is building a pipeline on this source |
| 2 | **GSE345058** CZ Biohub multimodal perturbation atlas (Liu … Leonetti, bioRxiv 2026.06.01.728087) `[new]` | ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345058 | GEO public 2026-08-26 | **A549** (lung adenocarcinoma) | dox-inducible Cas9 **knockout**, CROP-seq, 4 batches | 1,000 (4,000 sg; 211 NTC sg) | **21 / 21** (≥ 50 cells: 21). Median 631 cells per panel target; 606,075 cells, median 26,103 UMI | h5ad 25.0 GB; raw mtx 11.0 GB; cell annotations 40 MB (downloaded) | GEO, public |
| 3 | **Perturb-ME** (Wang, Gu, Frangieh …, bioRxiv 2026.08.17.745330) `[new]` | figshare 33264225 | 2026-08-23 | **A375** melanoma | genome-scale Cas9 KO; FACS-sorted HLA-low / HLA-high + NTC channel; Perturb-CITE-seq | ~19,340 detected (library 19,563) | **300 / 24** (≥ 50: 1). Median 11 cells per panel target; all perturbed cells sit in sorted gates | 3.36 GB (MOI1 h5ad 1.43 GB, read and deleted) | CC BY 4.0 |
| 4 | DLD-1 GSE337988 **Med/High** MOI arms (we hold Low1/Low2 as `data/external/dld1-lowmoi`) | GEO GSE337988 | 2026-07-10 | DLD-1 | CRISPRi, Zim3-dCas9, 5,450-target sublibrary 2 | 5,450 | **67 / 67**: the *same* 67 as Low1/Low2 (header-only range reads) | 155 MB per lfc matrix (+ se / z matrices) | GEO. Two more replicate arms at higher MOI, **zero new targets**. The genome-wide library's other sublibraries are not deposited |
| 5 | Head-to-head CRISPRko vs CRISPRi Perturb-seq (bioRxiv 2026.07.04.736492) `[new]` | zenodo 21136747 / 21282294 | 2026-07-02/09 | K562, **A549** | Cas9 KO and CRISPRi benchmark libraries | 25 genes (A549 d7) | **1** (AKT2) | 0.3–0.7 GB per table | CC BY 4.0 |
| 6 | DeSCOPE PC9 (bioRxiv 2026.04.13.718147) `[new]` | zenodo 21622013 | 2026-07-27 | **PC9** (lung adeno) | CRISPR (single + "+ctrl" pairs) | 50 | **3 / 3** (DEK, E2F4, SND1); 57,441 cells | 1.13 GB (read, deleted) | CC BY 4.0 |
| 7 | Song 2025 Jurkat HIV latency Perturb-seq (Nat Cell Biol 2025) `[new; surfaced by the EBI inventory]` | GEO GSE247599 | 2024-08-22 | **Jurkat** | Cas9 KO; DMSO / PMA-I GFP± sorted | 10 genes | **0** | 634 MB tar | GEO |
| 8 | 10x Genomics GEM-X Flex 1M K562 CRISPRi aggregate `[new to this repo]` | 10xgenomics.com/datasets/16-plex_GEM-X_Flex_1M_human_K562_CRISPR_aggregate | 2024-12 | K562 | CRISPRi, **Flex** probes, 16-plex | 847 genes (+ 280 lncRNAs; 2,692 "Ignore" guides) | **15** | feature reference 3.4 MB (downloaded) | 10x terms |
| 9b | Gastric-cancer dcPerturb-seq (bioRxiv 10.1101/2025.04.16.649236 v2) `[new]` | — | preprint v2 2026-09-03 | GES1 + **5 gastric cancer lines** | CRISPR dcPerturb-seq × 4 epigenetic drugs | 226 | unmeasurable: no accession found | — | — |
| 9 | ASD cortical-neuron CRISPRi (bioRxiv 2026.09.02.748925) | — | preprint 2026-09-04 | hESC-derived neurons (not cancer) | CRISPRi Perturb-seq | 1,408 | unknown: no accession found (306 also found none) | — | — |
| 10 | **X-Atlas/Pisces** | huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces | **not released** | Jurkat Rest / Active (full screens), HepG2 / iPSC (200-gene validation subset), Orion re-release | CRISPRi FiCS | Jurkat 10,872 / 10,878 | **unmeasurable** (no file; no gene list published) | 0 B today; ~400 GB class when full | CC BY-NC-SA 4.0 |
| 11 | **Illumina Billion Cell Atlas** | illumina.com (BioInsight) | commercial | 200–250 lines | CRISPRi/a × ~20k genes | ~10 M perts | — | ~20 PB | **paid alliance only** |
| 12 | Tahoe 300 M (Parse GigaLab; "genetically or chemically perturbed") | — | announced 2025-12, unreleased | tumour pools | genetic/chemical | — | — | — | — |
| 13 | Tahoe **EmeraldBay** `[new, out of scope]` | huggingface.co/datasets/tahoebio/EmeraldBay | 2026-06-05 | 52 cancer lines | **drugs** (91 incl. combinations) | — | 0 (chemical) | 57.7 GB | CC BY 4.0 |

**Already held or already closed and re-confirmed, not re-reported:** K562-GWPS, Orion HCT116/HEK293T, CD4
GSE314342, Nadig HepG2/Jurkat, Replogle RPE1/K562-essential, Feng iPSC, KOLF2, VCC2025 H1, GSE295214, Jiang 2025,
Broad-16 (figshare 33273600, 0/300), DLD-1 Low1/Low2, Myllia/Altair A549, PerturbFate A375 (6/300), CMap/LINCS 2020,
GPSAdb, scPerturb, IGVF, Frangieh, the 2026 CD4 reanalysis preprints, and the enhancer-CRISPRi CRC screen
(0/300 by construction).

---

## 3. Registries swept — each closed with one number

| registry | what was run (script / file) | result |
|---|---|---|
| **Hugging Face** | 32 keyword searches, 500 each, by last-modified (`hf_scan.py` → `hf_scan.json`); 30+ org listings (Xaira, Arc, Tahoe, Altos, Recursion, CZI/Biohub, Genentech, Broad, theislab, …) | **1,648** matching datasets. New since 06-01 and on-topic: `sunweiwei/ai4sci-virtual-cell-zeroshot-data` (2026-09-22, 192 GB: **a rebuild of Replogle + Nadig + Orion + VCC2025 for a VCC-2026-style zero-shot task**, nothing new), `nttduc/perturb_transfer_data` (Frangieh + Feng + sciplex3 re-hosts), gene-embedding tables, UCE embeddings, Tahoe EmeraldBay (chemical). **Zero new genetic screens in cancer lines.** Pisces: `usedStorage 0`. X-Cell model: 3 files, last modified 2026-03-16 |
| **GEO** (E-utilities; `geo_scan.py`) | 6 query families (Perturb-seq/CROP-seq; CRISPRi/a + single-cell; CRISPR screen + RNA-seq; HeLa + CRISPR; HNSCC lines + CRISPR; Jurkat + CRISPR), human GSE, PDAT 2025-06-01 → 2026-09-30 | 336 unique series; 123 single-cell CRISPR series listed with seen/new flags in `geo_sc_candidates.txt`. **Only GSE345058 (A549, 21/300) is a new cancer-line screen with > 100 targets.** HNSCC + CRISPR, 2026-08→: **0 series** |
| **Europe PMC** (papers + preprints; `epmc_scan.py`) | 4 title/abstract query families, 2025-10-01 → 2026-09-22 | 264 unique hits read by title. New data-bearing items are all in §2 (VIPerturb, GSE345058, Perturb-ME, head-to-head, ASD neurons). Everything else is methods, primary cells, mouse or non-coding screens |
| **bioRxiv API** full crawl (`biorxiv_crawl_par.py`) | every preprint 2026-08-01 → 2026-09-22, keyword-filtered | **10,601 scanned, 100 unique hits**; 2 new items (gastric dcPerturb-seq, LLM GEO atlas), see §3b |
| **figshare** (`registry_scan.py`) | 7 queries, published ≥ 2026-05-01 | 113 items. On-topic: VIPerturb export (1b), Perturb-ME (3), Broad-16 (held), HSPC Perturb-seq (**mouse**), Feng iPSC (held), CRC enhancer-screen table (non-coding) |
| **Zenodo** (API, page size ≤ 25) | 5 phrase queries, created ≥ 2026-04-01 (`zenodo_scan.txt`) | New on-topic: VIPerturb (1), head-to-head (5), DeSCOPE PC9 (6), **CIPHER** curated compendium (21729034; re-hosts CD4 / Orion / K562 plus lab sets: CRISPRa K562/RPE1 and fibroblast controls, none new in a cancer line), perturbmatch and sigrecon (held), GBM43 PDX in-vivo Perturb-seq (small), an A375 NF1 base-editing set (one gene) |
| **BioStudies / ArrayExpress** | same 7 queries | 156 fuzzy hits ≥ 2026-05-01, **0 transcriptomic perturbation screens** (closest: S-BIAD3248 CROPseq-multi, imaging; E-MTAB-16438, a fitness screen) |
| **scPerturb** | Zenodo versions of 7041849 | **v1.4, 2024-08-26 is still latest** |
| **Arc Virtual Cell Atlas** | `gh api` commits since 2026-08-01; repo contents | one README commit (2026-08-28). Contents still `scBaseCount`, `tahoe-100M`, `virtual-cell-challenge` |
| **CZI Virtual Cells Platform / Billion Cells** | site search, press | the only CRISPR dataset is the CD4 atlas (held; *Cell* 2026-08-28) |
| **EBI Perturbation Catalogue** | PR #380 inventory `datasets.tsv` (merged 2026-09-18; `dl/ebi_perturbseq_datasets.tsv`) | 30 datasets: all held or known except `song_2025_jurkat_hiv` (0/300, row 7) and `arce_2025` (primary Treg/Teff). **An independent curated inventory lands on our universe** |
| **LINCS / CMap** | web | no build after LINCS 2020; the 2026 Zenodo item is a repackaging |
| **10x / Parse / Scale** | web + 10x feature reference | no new public genome-scale screen; the 10x K562 Flex set is row 8 (15/300) |

### 3b. bioRxiv full crawl, 2026-08-01 → 2026-09-22

`biorxiv_crawl_par.py` (the serial `biorxiv_crawl.py` was replaced mid-run by a 6-worker version with the identical
keyword filter, `biorxiv_kw.py`, when the API slowed to ~10–12 s per call): **10,601 preprints scanned, 0 failed
pages, 107 keyword hits = 100 unique DOIs** (`biorxiv_hits_2026-08-01_2026-09-22.json`). Triage by title and
abstract: 29 flagged by a non-human keyword filter (bacteria, mouse, yeast; a few are human studies that
mention mice), about 40 methods, benchmarks or models, and the rest are screens with
fitness, imaging, enhancer or MPRA readouts. The data-bearing human items are all in §2 already: VIPerturb v2
(posted again 2026-08-20), GSE345058's preprint, Perturb-ME, Broad-16 (held), and the ASD neurons. Two were new to this scan:

- **Gastric-cancer dcPerturb-seq** (bioRxiv 10.1101/2025.04.16.649236, v2 2026-09-03): **226 GC-related genes ×
  6 gastric lines** (GES1 normal; SNU719, NUGC3, LMSU, SNU1750, HGC27) × 4 epigenetic drugs + DMSO, 625,866
  cells. These are genuinely new epithelial-carcinoma lines, but **no data accession was found** (GEO search on
  gastric + Perturb-seq/CRISPR screen returns no matching series; the full text served did not include a data
  statement). Panel overlap is **unmeasurable**. A GC-driver library of 226 genes is unlikely to exceed the 30-target
  bar, but that is an expectation, not a count. Worth re-checking when the paper is published.
- **LLM-built GEO bulk perturbation atlas** (Soul & Young, bioRxiv 2026.08.08.743502; perturbMatchData, Zenodo
  21792179): 6,802 signatures (3,571 KD, 1,779 KO, 1,452 OE) from 4,453 GEO series. Metadata read from
  `soulj/perturbMatchPaper/data/perturbation_metadata.rds` (`perturbmatch_atlas.py`). Knockdown or knockout
  signatures reach **76 panel targets** in any line. By context line: **HeLa 4** (EZH2, IFNAR2, ZBTB48, ZKSCAN3),
  **Jurkat 1**, **HNSCC 1** (EZH2). This is a second, independent bulk compendium that confirms §4.

---

## 4. Criterion (b): HeLa, CAL33 and Jurkat, and the epithelial carcinomas

| context | best public transcriptomic perturbation data | panel targets |
|---|---|---:|
| **B = HeLa** (cervical adenocarcinoma) | GPSAdb bulk knockdown/knockout groups in HeLa / HeLa S3 (242 groups, 182 genes) | **4** (DHX35, EAPP, EZH2, ZKSCAN3) |
| | LLM-built GEO atlas (5,350 KD/KO signatures) HeLa | **4** (EZH2, IFNAR2, ZBTB48, ZKSCAN3) |
| | any HeLa Perturb-seq / CROP-seq | **none found** (GEO, EPMC incl. full text, bioRxiv, HF, Zenodo, figshare, EBI inventory) |
| **C = CAL33** (HNSCC) | GPSAdb HNSCC lines (CAL27, FaDu, Detroit562, SCC-4/9/13/25/68, UM-SCC-1, …: 48 groups, 35 genes) | **0** |
| | CMap L1000 CRISPR-KO in **BICR6** (HNSCC) | 82 (on disk; landmarks only; closed by 176 inside the pooled `cmap_xpr` arm) |
| **A = Jurkat** (T-ALL) | GPSAdb Jurkat | **1** (EZH2) |
| | Nadig (held) / Datlinger 2017 (known) / Song 2025 (new) | **0 / 0 / 0** |
| | **X-Atlas/Pisces Jurkat Rest + Active** | unmeasurable, **unreleased** |
| epithelial carcinomas | DLD-1 (held, 67), **A549 GSE345058 (21, new)**, A549 Altair/Myllia (20/14, held), PC9 (3, new), MCF7/HT29/BxPC3 via Jiang (9, held), Calu-3 (6, held) | — |

GPSAdb context counts: `work/gpsadb_context_lines.py` on 270's saved metadata. CMap: `work/cmap_context_lines.py`.
**No HeLa, CAL33 or Jurkat perturbation set with ≥ 30 panel targets exists in public data.** For A/B/C that is a
measured absence, not a search gap. For D/E/F, the nearest epithelial-carcinoma sets are DLD-1 (held) and A549
(GSE345058, new but 21/300).

---

## 5. What each item would actually add, per panel target

`work/per_target_k.py`: the number of held lines (K562-GWPS, Orion HCT116, Orion HEK293T, CD4 at ≥ 20 cells in any
condition) covering each panel target is **4 for 255 targets, 3 for 35, 2 for 9, 1 for 1** (C5orf22). The ten
targets with ≤ 2 lines are ABCD1, C5orf22, CDK8, EPHB2, GLYR1, IL10RB, KIF21B, NICN1, SPATA20 and UNK.

| candidate | panel hits | adds a *line* for | of which currently ≤ 2 lines |
|---|---:|---:|---|
| VIPerturb K562 (≥ 20 cells) | 284 | **30** (targets absent from our K562-GWPS donor) | **5**: ABCD1, EPHB2, KIF21B, NICN1, SPATA20 |
| GSE345058 A549 | 21 | 21 | 1: C5orf22 |
| DLD-1 (held; reference) | 67 | 67 | 1: C5orf22 |
| Pisces Jurkat (if released) | unknown; library ~10.9k expressed genes | likely most of the panel | would take k from 4 → 5 on most targets **and** is the context-A line itself |

For the other 254 panel targets, VIPerturb is a **second independent K562 measurement** (a different lab,
library, guides and readout chemistry) of a line we already hold. Under 159's `1/k` law it does not increase *k*.
What it can do is denoise the K562 row and give a per-target reliability estimate for K562: two K562 screens
disagreeing on a target flags it. Its median 47 cells per target is shallower than K562-GWPS, so the gain is
bounded; price it before building (414's harness, not here).

Training value (shared targets with held donors, `work/shared_targets.json`): the A549 1,000-gene library shares
**877 / 927 / 929** targets with K562-GWPS / HCT116 / HEK293T and 550 with Nadig-Jurkat's essential set. That is
denser than Jiang (171 shared) and deeper per target than Nadig HepG2 (565 vs 49 median cells). It is a Cas9
*knockout* screen, however, and the KO class has priced weakly before (Frangieh 166/225, `cmap_xpr` 176).

---

## 6. Ranked recommendation for the PI

| rank | action | what it buys | cost | risk / caveat |
|---|---|---|---|---|
| **1** | **VIPerturb-seq K562 → a second K562 donor table** (pseudobulk from `genome_wide_binA/B/C.RDS`; the 7,949 `NO-TARGET` cells are in every bin and in the filtered object; control label is `NO-TARGET`, not `non-targeting`) | the only public K562 rows for **ABCD1, EPHB2, KIF21B, NICN1, SPATA20** (5 of the 10 thinnest panel targets) plus 25 more K562 gap-fillers; a replicate K562 measurement for 254 targets (reliability weight / denoise); Flex-style probe readout matching the challenge chemistry | 10.2 GB (three bins; R or `rds2py` to read Seurat v5) — do it on the cloud box or as one scheduled download; **nothing > 2 GB was fetched here** | K562 is not a new line (k unchanged for 270 targets); 47 cells/target is shallow; guides are CRISPRi with a different effector. Price on PLOLO before shipping. The third-party 470 MB export is **not** a substitute: its control regex missed `NO-TARGET`, so it uses a pseudo-control drawn from perturbed cells, and its rows are 20-cell means |
| **2** | **Keep the Pisces poll** (`data/pisces-watch.log`, last check 17:00Z today: `data_files=0`) and **pre-register the Jurkat build now** | Jurkat Rest/Active are the context-A line and, per Xaira's `docs/data.md`, will be released as **full screens** (~10.9k perturbations, 2.8 M cells each) | 0 now; ~100–150 GB when it lands | CC BY-NC-SA 4.0. Jurkat perturbation data is legal training data, but it *is* context A, so any A/B/C gain from it must not be read as D/E/F transfer. It moves D/E/F only through k 4 → 5 |
| **3** | **GSE345058 A549 as a training line** for the cross-line residual / key models (not a donor) | a deep (565 cells/target) epithelial carcinoma with 877–929 shared targets | 40 MB annotations (have); counts 11.0 GB mtx or 25.0 GB h5ad | Cas9 KO with DSB artifacts; essential-heavy library; 21 panel targets only |
| 4 | DLD-1 Med/High arms | replicate reliability for the 67 DLD-1 targets | 2 × 155 MB (+ se matrices) | higher-MOI stress (the paper reports it); zero new targets |
| — | Perturb-ME, PC9, head-to-head A549, Song Jurkat, 10x K562 Flex, EmeraldBay | — | — | below the bar: 24 (sorted gates) / 3 / 1 / 0 / 15 / chemical |
| — | Illumina Billion Cell Atlas | would be decisive (200+ lines) | commercial membership | not public; nothing to do but note it |

**Bottom line for the 1/k law:** the number of public panel-covering independent *lines* is unchanged today (K562,
HCT116, HEK293T, CD4; plus the panel-fitted KOLF2/Feng, which project to zero in October per 225). The next
increment arrives only if Xaira uploads Pisces or Illumina opens its atlas.

---

## 7. What I did not do, and the limits of this scan

- **No panel number for Pisces or the Illumina atlas.** Neither has a public perturbation list.
- **VIPerturb's 284 is a manifest count** (`cell_count ≥ 20` per gene across bins A/B/C). It was not recomputed
  from cells, and guide-level knockdown QC is not applied (the authors' filtered object keeps 6,725 genes, 124 on the panel).
- **GSE345058's 21 is from `cell_annotations.csv.gz`** (every labelled cell). QC is the authors'; there is no
  separate knockdown/KO-efficiency filter.
- **Perturb-ME's 24 counts cells in the HLA-low and HLA-high gates.** Those gates are selected on a phenotype, so even
  those 24 are biased samples.
- Keyword sweeps can miss a deposit that never uses the words Perturb-seq / CROP-seq / CRISPRi / CRISPR screen
  / perturbation atlas in its title or abstract. Four independent routes (GEO, EPMC, bioRxiv full crawl, repository
  keyword search) plus two curated inventories (EBI's, Altos') all converge on the same universe.
- One web search was blocked by a safety filter (a query about the Perturb-ME immune phenotype). The dataset itself was
  then fully characterised from its figshare record.
- The ASD-neuron preprint's data location was not found (as in 306). It is not a cancer line.

---

## 8. Files (`work/`)

| file | content |
|---|---|
| `known_answer.py`, `known_answer2.py`, `known_answer*.json` | §1 |
| `overlap.py` | literal + HGNC-resolved panel intersection (both sides mapped) |
| `hf_scan.py`, `hf_inspect.py`, `hf_scan.json` | Hugging Face sweep |
| `geo_scan.py`, `geo_scan_*.json`, `geo_sc_candidates.txt`, `geo_older.txt` | GEO sweep |
| `epmc_scan.py`, `epmc_*.json` | Europe PMC sweep |
| `biorxiv_crawl.py` (serial, superseded), `biorxiv_crawl_par.py`, `biorxiv_kw.py`, `biorxiv_hits_*.json` | bioRxiv full crawl |
| `perturbmatch_atlas.py/.json`, `dl/perturbmatch_atlas_metadata.*` | LLM GEO atlas (§3b) |
| `epmc_fulltext_context*.json` | Europe PMC full-text search for HeLa / Jurkat / HNSCC Perturb-seq (§4) |
| `candidates.csv` | §2 as a machine-readable table |
| `registry_scan.py`, `registry_scan.json`, `zenodo_scan.txt` | figshare / Zenodo / BioStudies |
| `viperturb.py`, `viperturb_overlap.json` | row 1 |
| `figshare_derivatives.py`, `figshare_derivatives_overlap.json` | rows 1b, 3 (library, model) |
| `perturbme_obs.py`, `perturbme_obs_overlap.json` | row 3 (cells per target) |
| `gse345058.py`, `gse345058_overlap.json` | row 2 |
| `dld1_arms.py`, `dld1_arms_overlap.json` | row 4 |
| `small_sets.py`, `small_sets_overlap.json`, `pc9_overlap.json` | rows 5, 6 |
| `song_jurkat_overlap.json`, `tenx_k562_flex_overlap.json` | rows 7, 8 |
| `gpsadb_context_lines.py/.json`, `cmap_context_lines.py/.json` | §4 |
| `per_target_k.py/.json`, `shared_targets.json` | §5 |
| `dl/` | the small source files (metadata, libraries, manifests, READMEs); 70 MB |
