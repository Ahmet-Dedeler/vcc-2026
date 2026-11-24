# 301 — downloaded sources (URL, licence / terms, sha256)

All downloads are in this directory. **Total downloaded: 116,142,365 B = 116.1 MB** (including the 434,429 B PMC
landing page read for the licence statement), plus the two derived promoter FASTAs built from Ensembl REST
(11,731,720 B) — **grand total 127.9 MB**, under the mandate's 200 MB cap.
Nothing was written outside `experiments/301-metadata-acquisitions/`.

## 1. Replogle et al. 2022 — the VCC 2025 library's protospacer source

Arc's own account of the 2025 challenge library (Arc Institute, "Behind the Data of the Virtual Cell
Challenge", <https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge>) states that the dual-guide
construct design uses "protospacer sequences and cloning strategy from (Replogle, Bonnar, et al. 2022)".
That paper's Table S1 is the library table.

| file | URL | licence / terms |
|---|---|---|
| `replogle2022_mmc1.xlsx` (Table S1: "Full sgRNA content of the K562 day 8 genome-scale library, K562 day 6 essential-wide library, and RPE1 day 7 essential-wide library") | <https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc1.xlsx> | Replogle et al., *Cell* 185(14):2559–2575.e28 (2022), doi 10.1016/j.cell.2022.05.013, PMC9380471. **Open access, CC BY 4.0** (stated on the PMC record: "This is an open access article under the CC BY license"). |
| `replogle2022_mmc2.xlsx` (Table S2: per-perturbation summary statistics incl. percent knockdown) | <https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005979-mmc2.xlsx> | same |

| file | bytes | sha256 |
|---|---:|---|
| `replogle2022_mmc1.xlsx` | 1,340,967 | `1e8d67490d562a48c03a641fc6170119548b19bc6dfe64a7cca1f6c3b6d506ae` |
| `replogle2022_mmc2.xlsx` | 3,000,043 | `364748fe04e9a9cdbcf239d8012c94cd17f72880b43a84817e7505bdb78e7177` |

`pmc9380471.html` (the PMC landing page, fetched only to read the licence statement and enumerate the
supplementary files) is also in this directory.

## 2. Horlbeck et al. 2016 — hCRISPRi-v2, the parent library, and its TSS annotation

Supplementary file 3 fills guide ids absent from Replogle's K562 day-8 library (the 14 non-targeting pairs) and
supplies the rank profile that shows the missing 20 A/B pairs are not reconstructible. **Supplementary file 2 is
the TSS annotation the library was designed against and is the basis of the accepted promoter window set** (see
§3 and `RESULTS.md` §2).

| file | URL | licence / terms |
|---|---|---|
| `horlbeck2016_supp3_hCRISPRi-v2.xlsx` (Supplementary file 3: "Library composition of hCRISPRi-v2 and hCRISPRi-v2.1"; 209,078 sgRNAs) | <https://cdn.elifesciences.org/articles/19760/elife-19760-supp3-v2.xlsx> | Horlbeck et al., *eLife* 5:e19760 (2016), doi 10.7554/eLife.19760. **CC BY 4.0** (eLife default). |
| `horlbeck2016_supp2_TSS.xlsx` (Supplementary file 2: "TSS annotations for hg19 and mm10 genomes"; 20,883 human rows, FANTOM-CAGE-derived P1 / P2 / P1P2 TSSs) | <https://cdn.elifesciences.org/articles/19760/elife-19760-supp2-v2.xlsx> | same |

| file | bytes | sha256 |
|---|---:|---|
| `horlbeck2016_supp3_hCRISPRi-v2.xlsx` | 31,182,240 | `afcc9b4c1f1ce2302a0f360c99239663caf6443644e12512b888b1de2bbde6bf` |
| `horlbeck2016_supp2_TSS.xlsx` | 2,790,673 | `2808fde8c583e23d01a4b2cd07fa15a90d49472cbe22c0d8bbb99f2f9eea68a7` |

## 3. Ensembl — promoter windows (item 2)

Two window sets were built. **The second is the accepted one**; the first is kept as the record of an attempt its
known-answer row rejected (22/68 vs 68/68 on-target recovery — `RESULTS.md` §2).

| file | URL | licence / terms |
|---|---|---|
| `biomart_genes.tsv` (86,411 human genes: Ensembl id, symbol, chromosome, start, end, strand, biotype) | `https://www.ensembl.org/biomart/martservice?query=…` (query in `biomart_genes.xml`), GRCh38, Ensembl current release | Ensembl data are released with **no restrictions** (Ensembl terms of use: "no restrictions on use of Ensembl data"). |
| `promoters_tss1000_grch38.fa.gz` — **REJECTED set**: TSS ± 500 bp around the Ensembl gene-level start/end, one window per measured symbol | built by `fetch_promoters.py` from `https://rest.ensembl.org/sequence/region/human?coord_system_version=GRCh38` (POST, 50 regions per request) | same |
| `promoters_cage_tss_hg19.fa.gz` — **ACCEPTED set**: ±500 bp around the Horlbeck 2016 CAGE TSS peaks (§2), hg19/GRCh37 | built by `fetch_promoters_hg19.py` from `https://grch37.rest.ensembl.org/sequence/region/human` (POST, 50 regions per request) | same |

| file | bytes | sha256 |
|---|---:|---|
| `biomart_genes.tsv` | 4,972,614 | `da62f6905a37dc5fa65be0e4fc6e653c345bbb62c7cf700a9df7024ee4785787` |
| `promoters_tss1000_grch38.fa.gz` (17,569 windows) | 5,451,999 | `6896a0f06c0c7b42c18d830d0b3c47a6fdbbde355998b2fb2412d09a3e5c0310` |
| `promoters_cage_tss_hg19.fa.gz` (18,387 windows, 16,936 genes, 20,602,846 bp) | 6,279,721 | `149a5a7038e62269adb30693d3be10e1aa6502b16c20cceac554fcc73639940e` |

FASTA header format: `>SYMBOL|<ENSG or P1/P2/P1P2/ENST>|chrom|strand|TSS|window_lo|window_hi|len`; the sequence is
always the **forward genomic strand** of `[window_lo, window_hi]`, and `strand` (+1 / −1) gives the sign convention
for the recorded offset-from-TSS. The two sets are on different assemblies (GRCh38 / GRCh37) — do not mix them.

## 4. Joung et al. 2023 — the TF-ORF id table and the per-cell TF assignment (item 3)

| file | URL | licence / terms |
|---|---|---|
| `joung2023_mmc2.xlsx` (Table S1A: the TF ORF library — `Name` = `TFORF####`, RefSeq gene name, isoform id, RefSeq/Gencode id, 24-nt ORF barcode sequence; 3,550 ORFs) | <https://ars.els-cdn.com/content/image/1-s2.0-S0092867422014702-mmc2.xlsx> | Joung et al., *Cell* 186(1):209–229.e26 (2023), doi 10.1016/j.cell.2022.11.026. **Open access, CC BY 4.0.** |
| `joung2023_mmc3.xlsx` (Table S2: differentiation statistics per `TFORF####-SYMBOL`) | <https://ars.els-cdn.com/content/image/1-s2.0-S0092867422014702-mmc3.xlsx> | same |
| `joung2023_mmc4.xlsx` (Table S3: NMF / correlation cluster order per `TFORF####-SYMBOL`) | <https://ars.els-cdn.com/content/image/1-s2.0-S0092867422014702-mmc4.xlsx> | same |
| `GSM6719970…GSM6719979_210322_TFmap_S01…S10.csv.gz` (per-cell TF-barcode assignment for atlas sublibraries S01–S10; columns: SHARE-seq cell barcode `R1.xx,R2.xx,R3.xx`, 24-nt ORF barcode, read count) | `https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6719nnn/GSM<id>/suppl/GSM<id>_210322_TFmap_S<nn>.csv.gz` | NCBI GEO series **GSE217460**. NCBI/GEO data are **public domain in the US** (NLM/NCBI data usage policy); the series carries no additional restriction. |
| `geo_GSE217460_filelist.txt` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217460/suppl/filelist.txt> | same |

| file | bytes | sha256 |
|---|---:|---|
| `joung2023_mmc2.xlsx` | 3,124,062 | `ffe5f0f4af7b6fb70a0dbc3cea3f94cedf1f0537cd2455953ad3d236bda9a7f6` |
| `joung2023_mmc3.xlsx` | 265,648 | `4d199de0c44cc958c88f857f622ca7f442d4fe6aa280d51e3001f333984a8884` |
| `joung2023_mmc4.xlsx` | 133,877 | `6ecdee87e75424d5e3b16ea31b83eba88acb753c9c9197d4c057ff944a0b052b` |
| `geo_GSE217460_filelist.txt` | 3,343 | `9d13e06a18486d0aa956962dd4604e40b9ac8b97f59fced97970d718d4fbe723` |
| `GSM6719970_210322_TFmap_S01.csv.gz` | 6,721,959 | `1abc53f41f85cf662f39365e3b3b132e7fbe8e063887ca0517c23aeb8a625e7a` |
| `GSM6719971_210322_TFmap_S02.csv.gz` | 6,577,831 | `6d2fec764b9da55f34b83290b988cf6988d1831a155ffadf8cff9ce7a75680b7` |
| `GSM6719972_210322_TFmap_S03.csv.gz` | 5,925,334 | `9e0f9845c97460bf1e71cf41f246210913cd7923d16ae52121710655cd0bbc09` |
| `GSM6719973_210322_TFmap_S04.csv.gz` | 6,840,179 | `f5c1d5a3770d3544e16a363b36d0f01d4558bf9b6a6e6c0e1e6ae3d3bd54ac1b` |
| `GSM6719974_210322_TFmap_S05.csv.gz` | 7,306,913 | `56943d2a9cb86774a4fa75fc531ab0d5321095fa8b65ababc107d56cccb62d2e` |
| `GSM6719975_210322_TFmap_S06.csv.gz` | 7,408,569 | `81b869ba4e1ee9b6c84662df8dd2a3e6ec5423235d85ca8ff53b6471a47276b3` |
| `GSM6719976_210322_TFmap_S07.csv.gz` | 7,305,091 | `749207c0677a9574a2cfc61bd2428e8c2e7187c60cc615477e784cf42ff0ea4b` |
| `GSM6719977_210322_TFmap_S08.csv.gz` | 6,572,647 | `516004e06ad795486801b25f9aa30c975b1dac81753cce629024a69efb8ef76b` |
| `GSM6719978_210322_TFmap_S09.csv.gz` | 7,173,457 | `b377862b5ba8823efaeb7adca8d93dce0ef1865f3d29b6acd2ba9b68c36a9f08` |
| `GSM6719979_210322_TFmap_S10.csv.gz` | 7,062,489 | `ce573ba6822ef7b15c52f8249876d35b7ff6f8ccc30be2d6f2d7154c0d4ea7d0` |

## Local files read (not downloaded)

Read-only, `obs` / `var` / `uns` and coordinates only:

- `data/external/vcc2025/{train,validation,test}/adata_*.h5ad` — `obs/guide_id` categories and `var/_index`.
- `data/external/vcc2025/gene_names.csv` — the 18,080 measured symbols.
- `data/external/scperturb/JoungZhang2023_atlas.h5ad` — `obs` (`cell_barcode`, `perturbation` codes and
  categories, `batch`, `louvain`, `cell_line`, `nperts`), `var` keys, `uns` keys. **`X` was never opened.**
  sha256 `8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22`, 5,804,454,475 B.
- `data/external/scperturb/JoungZhang2023_combinatorial.h5ad` — not opened;
  sha256 `ce5188bd3f924f9c79fea6d0b2fb49e48a71d0cf6a2a5010ee9a9bd6e3289a01`, 834,073,300 B.
- `experiments/298-reviewer-checks/work/seed298.csv`, `seed298_plan.md` — the 34 arms and the recipe.
- `experiments/293-challenge-construct-ceiling/work/run293_*.npz` — key listing only (see `RESULTS.md` §4).
- The scPerturb curation notebook `dataset_processing/notebooks/JoungZhang2023.ipynb` was read from
  <https://raw.githubusercontent.com/sanderlab/scPerturb/master/dataset_processing/notebooks/JoungZhang2023.ipynb>
  (MIT-licensed repository) into the session scratchpad, not into this directory.
