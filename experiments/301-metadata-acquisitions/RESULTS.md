# 301 — three metadata acquisitions, and the seed scan they unlock

**SPACERS: FOUND (Replogle et al. 2022 *Cell* Table S1, `mmc1.xlsx`, tab `TabA_K562_day8_library`; 14 non-targeting
pairs filled from Horlbeck et al. 2016 *eLife* Supplementary file 3) — 328 / 348 guide_ids, 34 / 34 seed-scan arms;
PROMOTERS: FETCHED (Ensembl GRCh37 REST on Horlbeck 2016 Supplementary file 2's CAGE TSS annotation, 18,387 windows
over 16,936 genes; a first GRCh38 gene-level-TSS attempt was built and REJECTED by its known-answer row);
JOUNG IDS: 3,136 / 3,369 resolved (99.48 % of the atlas's 1,145,823 cells; 1,647 distinct TF symbols);
SEED SCAN: 4 / 13 failing arms with a candidate off-target (working arms 4 / 21) at the primary reading —
12-nt seed, CRISPRi active window, distal — Fisher p = 0.679; no failing-vs-working difference at any seed length.**

2026-09-07, ~50 min inside the 2.5 h box. No `vcc`, no `cell-eval2`, no build, no submission. Nothing written
outside `experiments/301-metadata-acquisitions/` except the R4 candidate note appended to
`docs/CONFIRMATION-LEDGER.md` (§3). `JoungZhang2023_atlas.h5ad` was opened for `obs` / `var` / `uns` only — **`X`
was never read**. R3's eight files were not touched. Downloads **116.1 MB** + two derived promoter FASTAs
(11.7 MB) = **127.9 MB**, each with URL, licence/terms and sha256 in `work/sources301.md`. Peak RSS well under the
6 GB cap (largest object: 512 k-row CSV dictionaries and a 1.15 M-element code array).

---

## 1. The 2025 VCC sgRNA library table — FOUND

**Where it is.** Arc's own account of the 2025 library ("Behind the Data of the Virtual Cell Challenge") says the
dual-guide constructs use "protospacer sequences and cloning strategy from (Replogle, Bonnar, et al. 2022)". That
paper's **Table S1** — `1-s2.0-S0092867422005979-mmc1.xlsx`, tab `TabA_K562_day8_library`, 11,294 rows — is the
table: columns `unique sgRNA pair ID`, `gene`, `transcript`, `ensembl gene id`, `sgID_A`, `targeting sequence A`,
`sgID_B`, `targeting sequence B`. The article is open access under **CC BY 4.0** (PMC9380471). Naming matches the
2025 `obs/guide_id` scheme exactly: a guide id `SOX4_P1_A|SOX4_P1_B` is Replogle's `gene`=`SOX4`,
`transcript`=`P1` row, positions A and B; `P1`/`P2`/`P1P2` and the `ENST…` forms are hCRISPRi-v2 CAGE-TSS
annotations (Horlbeck et al. 2016 Supplementary file 2). **Non-targeting ids match verbatim**: Replogle's TabA
carries 585 rows with `sgID_A|sgID_B` = `non-targeting_00005|non-targeting_01689`, which is literally the first
2025 NTC guide id.

**The name-match check** (`work/build_spacers.py` → `work/vcc2025_spacers.csv`). The three splits' `obs/guide_id`
categories were read (obs only): train 189, validation 85, test 136 → **348 distinct guide ids**, 317 gene-targeting
+ 31 non-targeting.

| route | resolved |
|---|---|
| Replogle TabA by `gene`+`transcript` | 297 / 317 targeting |
| Replogle TabA by exact `sgID_A\|sgID_B` | 17 / 31 non-targeting |
| Horlbeck 2016 hCRISPRi-v2 by sgID | 14 / 31 non-targeting |
| **total** | **328 / 348** |

All 656 resolved spacers are 20 nt and pure ACGT. **All 34 arms of the 17 two-construct targets resolve**, so item 4
is unblocked at full coverage.

**The 20 that do not resolve** are gene-targeting constructs absent from Replogle's K562 genome-scale library —
`AASS, ANTXR1, ATP1B1, CAMSAP2, CNN3, COL6A1, CXCL12, DZIP3, KDR, LAD1, NREP, PMEL, RAB3B, SH3BP4, SNCA, SRC,
SV2A, TGFBR2, THY1, WFS1`, i.e. genes not expressed in K562 but expressed in H1. Their **individual spacers exist**
in hCRISPRi-v2 (20/20 gene-transcript keys present) but the **A/B pair cannot be reconstructed**: across the 10,271
Replogle pairs whose two sgIDs are both in hCRISPRi-v2, only 3,019 (29 %) are the rank-1 + rank-2 pair; the rest are
spread over (1,3) 1,263, (1,4) 724, (2,3) 624, … Replogle applied a selection rule that is not the hCRISPRi-v2
selection rank, so guessing a pair here would be fabrication. (1,206 / 10,271 pairs also differ in at least one
protospacer string between the two tables — hCRISPRi-v2 records a 5′ cloning G "whether or not it is present in the
genome"; Replogle's table is the authoritative one for the actual library and is what `vcc2025_spacers.csv` uses.)

**Where else I looked, and what was there:** `github.com/ArcInstitute/arc-virtual-cell-atlas`
`virtual-cell-challenge/README.md` — count matrices and `pert_counts_*.csv` under
`gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/`, **no guide library table**; the Arc "Behind the
Data" page — the Replogle citation, no URL or accession for a library table; the 2025 *Cell* paper
S0092-8674(25)00675-0 — reachable, no library table of its own; the three 2025 h5ads and `data/challenge/` (298
already showed 0/410 ACGT strings). The pointer, not the table, is what Arc published; the table is Replogle's.

## 2. Promoter sequences — FETCHED (second attempt; the first was rejected by its known-answer row)

**Attempt 1, REJECTED.** `work/fetch_promoters.py`: Ensembl BioMart GRCh38 gene coordinates + strand (86,411
genes), TSS = `start_position` on +, `end_position` on −, window TSS ± 500, sequence from
`rest.ensembl.org/sequence/region/human?coord_system_version=GRCh38` → `work/promoters_tss1000_grch38.fa.gz`,
17,569 windows, 5,451,999 B, sha256 `6896a0f06c0c7b42c18d830d0b3c47a6fdbbde355998b2fb2412d09a3e5c0310`.
**Known-answer row KA-1 failed: 22 / 68.** Only 22 of the 68 arm spacers matched their own target gene's window at
19 nt + NGG. The gene-level Ensembl boundary is frequently far from the CAGE peak the library targets (both ARPC2
arms, both CAST arms, both KDM2B arms, … all missed), so this window set is wrong for the question and is kept only
as the record of the rejected attempt.

**Attempt 2, ACCEPTED.** `work/fetch_promoters_hg19.py`: the TSS annotation the library was *designed against* —
Horlbeck et al. 2016 **Supplementary file 2**, "TSS annotations for hg19", 20,883 human rows
(`gene, transcript, chromosome, strand, TSS source, Primary TSS 5'/3', Secondary TSS 5'/3'`) — restricted to the
18,080 symbols measured in the 2025 h5ads, window = `[min(peak coords) − 500, max(peak coords) + 500]`, sequence
from `grch37.rest.ensembl.org/sequence/region/human`. → **`work/promoters_cage_tss_hg19.fa.gz`, 18,387 windows over
16,936 genes, 20,602,846 bp, 6,279,721 B, sha256
`149a5a7038e62269adb30693d3be10e1aa6502b16c20cceac554fcc73639940e`.** Header:
`>SYMBOL|P1|P2|P1P2 or ENST|chrom|strand|TSS|lo|hi|len`; sequence is the forward genomic strand, `strand` gives the
sign convention for offsets. **KA-1: 68 / 68.** 1,144 of the 18,080 measured symbols have no hCRISPRi-v2 TSS
annotation and are absent from the scanned universe (mostly non-coding / symbol drift); this is a stated coverage
limit of the scan, not of the download.

## 3. JoungZhang 2023 ORF id → symbol — RESOLVED, 3,136 / 3,369

**Why it was blocked.** scPerturb's own curation notebook
(`dataset_processing/notebooks/JoungZhang2023.ipynb`, cell 11) renames the GEO h5ad's `TF` column to
`perturbation` while reading `obs` column-wise through raw h5py, with the comment
`# TODO: Perturbation here is just a number. Find out what it means!!!`. The numbers are the **category codes** of a
categorical whose category array was left behind in the source file's `uns`; the atlas on disk therefore carries
codes `'0'…'3368'` and no vocabulary.

**How it was resolved** (`work/joung_map.py`, output `work/joung_orf_map.csv`, meta `work/joung_map_meta.json`).
Joung's per-cell TF-barcode assignment is public as GEO **GSE217460** `…_210322_TFmap_S<nn>.csv.gz` (one row per
SHARE-seq cell barcode: `R1.xx,R2.xx,R3.xx`, 24-nt ORF barcode, read count). Ten sublibraries S01–S10 were
downloaded (66.9 MB) and joined to the atlas's `obs.cell_barcode` (`R1.xx,R2.xx,R3.xx,P1.yy-…`) by the R-triple,
after identifying which `P1.yy` group is which sublibrary. That identification is itself measured, not assumed:

| sublibrary | best `P1` group | purity | runner-up |
|---|---|---:|---|
| S01 | P1.22 | 0.632 | P1.32 0.154 |
| S02 | P1.30 | 0.630 | P1.32 0.154 |
| S03 | P1.38 | 0.667 | P1.32 0.162 |
| S04 | P1.46 | **1.000** | P1.32 0.153 |
| S05 | P1.54 | **1.000** | P1.32 0.148 |
| S06 | P1.62 | 0.646 | P1.54 0.164 |
| S07 | P1.07 | **1.000** | P1.32 0.145 |
| S08 | P1.15 | **1.000** | P1.07 0.212 |
| S09 | P1.23 | **1.000** | P1.31 0.162 |
| S10 | P1.31 | **1.000** | P1.32 0.152 |

(Purity = fraction of joined cells whose ORF barcode equals the modal barcode of their perturbation code. Six of the
ten separate at 1.000; the four ~0.63 groups are the three P1 groups holding ~97 k cells, i.e. more than one
sublibrary pooled. Every runner-up is ≤ 0.21.)

ORF barcodes were then mapped to `TFORF####` + gene symbol through **Joung Table S1A** (`mmc2.xlsx`, 3,550 ORFs;
CC BY 4.0). Codes with ≥ 10 joined cells and ≥ 80 % modal purity give **2,146 anchors**, and the anchors are
**strictly increasing in TFORF number with 0 violations over 2,145 consecutive pairs** — confirming that the codes
are the rank of the label `TFORF####-SYMBOL` in sorted order. Codes between two anchors were filled only where the
number of missing codes exactly equals the number of library ORFs strictly between the two anchor ORFs (forced, no
choice); everything else is left `AMBIGUOUS`.

- **3,136 / 3,369 codes resolved** (2,146 ANCHOR + 990 INTERPOLATED), **233 ambiguous** (median 8 cells each,
  5,917 cells total).
- **99.48 % of cells covered**: 1,139,906 / 1,145,823.
- **Independent confirmation the ordering is right**: Table S1A's last two entries are `TFORF3549` = **GFP** and
  `TFORF3550` = **mCherry**, and the map lands them on codes 3367 (49,539 cells) and 3368 (78,274 cells) — the two
  largest categories in the file, i.e. exactly the two controls, predicted by the ordering rather than assumed.
- 1,647 distinct TF gene symbols (median 1 ORF isoform per symbol, max 19); 1,018,010 perturbed cells.

## 4. T7 — the seed scan

`work/seed_scan.py` (scan + known-answer rows) and `work/seed_contrast.py` (final reading). Recipe = 298's
`seed298_plan.md` steps 1–7, with steps 2 and 4 supplied by items 1 and 2.

**Known-answer rows.** KA-2: the failing/working split read back from `seed298.csv` is **13 / 21**, as 298 recorded.
KA-1: with the own-target exclusion switched off, all **68 / 68** arm spacers recover their intended site (3′ 19 nt
exact + NGG) inside their own gene's window — this validates the spacer table, the promoter windows, the strand
handling and the PAM rule together, and it is the row that rejected the GRCh38 window set (22/68).

**Method.** PAM-anchored: enumerate NGG sites in each window and its reverse complement, test the L nt 5′ of each
against the 130 spacers' seeds for L ∈ {12, 11, 10, 9}; record gene, orientation, PAM and strand-aware
offset-from-TSS. Groups: 13 failing arms, 21 working arms, and — as a background group — the **31 non-targeting
constructs' 62 spacers**, scanned identically. 4,084 raw hits.

**Two filters applied before reading.** (a) *CRISPRi active window*: keep only hits with offset in
[−50, +300], the window in which dCas9-KRAB represses. (b) *Distal*: drop hits whose candidate promoter window lies
within 10 kb of any window of the arm's own gene — these are the arm's own locus re-found under a neighbour's
symbol. This removes 24 hits over 3 construct-gene pairs, including the only two 12-nt CRISPRi-window hits for
`ZNF581_P2`, which land on **ZNF580** at **0 bp** window separation (adjacent paralogues sharing a promoter region,
not an off-target). Without this filter the headline would have been 5/13 and would have been wrong.

**Result.**

| seed L | failing | working | non-targeting | Fisher OR (F vs W) | p |
|---:|---|---|---|---:|---:|
| **12** | **4 / 13** | **4 / 21** | 4 / 31 | 1.89 | **0.679** |
| 11 | 9 / 13 | 15 / 21 | 11 / 31 | 0.90 | 1.000 |
| 10 | 12 / 13 | 19 / 21 | 23 / 31 | 1.26 | 1.000 |
| 9 | 13 / 13 | 21 / 21 | 30 / 31 | — | 1.000 |

(CRISPRi window, distal. Mean candidates per construct at L = 12: failing 0.38, working 0.29, non-targeting 0.13;
Mann-Whitney p = 0.473. Over the whole ±500 window the L = 12 row is 8/13 vs 11/21, Fisher p = 0.728. At L ≤ 10 the
scan saturates on background — the non-targeting group reaches 30/31 — so only L = 12 and L = 11 carry information,
and neither separates the groups.)

**The candidate list, L = 12, CRISPRi window, distal:**

| arm | class | candidate off-target(s) | offset from TSS |
|---|---|---|---|
| `CAST_P1` | FAIL | UTP18 | +259 |
| `GRB10_P2` | FAIL | NEUROG3, SNX8 (both spacers) | +274; +291 / +299 |
| `KDM2B_ENST00000377071.4` | FAIL | WRN | −24 |
| `KLF10_P1` | FAIL | ST3GAL3 | −27 |
| `DPF2_P1` | work | CEP120 | +247 |
| `IGF1R_P1` | work | IFITM5 | +114 |
| `MED15_P1` | work | AFF4, TMCC1, ZSWIM5 | +106 / +92 / +208 |
| `TMSB4X_P2` | work | GNG10 | +141 |

**Stopped at the candidate list, as the mandate directs.** 293's saved arrays are per-target summary scalars, not
per-construct residual rows: `experiments/293-challenge-construct-ceiling/work/run293_{train,validation,test}[_thin*].npz`
contain, per cell-count regime `{25,50,100,200}`, the keys `<n>|{a,b,d,bfix,a2,x}_{cos,pds}`, `<n>|{a,b,d}_own`,
`<n>|targets` (4–5 two-construct targets per split) and `<n>|panel_{cos,pds,targets}` (~100 targets) — every array
is one number per target, 18–28 KB per file. **No gene-column residual row for any construct is saved**, so the
expression test (plan step 6) cannot read a stored row. It is rebuildable without any download from
`agg293.py` + `lib293.py` over the 2025 h5ads; that is the next run's work, not this one's.

**Reading.** This is an existence-proof design at n = 13 vs 21, and it produces no existence proof of the
hypothesis: failing constructs are **not** enriched for promoter seed matches at any seed length, on either window
definition, and the largest apparent signal (ZNF581 → ZNF580) is a shared-locus artefact. Four failing arms carry a
distal 12-nt seed + NGG inside a candidate gene's CRISPRi window — but so do four working arms, at a rate the
Fisher test cannot separate (p = 0.679), and both sit close to the non-targeting background (4/31). This is
consistent with 296's direct measurement (failing constructs' responses are reproducible but carry no shared vector
signature) and with the published 0.8–4.7 % off-target rate for a validated CRISPRi library: a 13-arm scan is
underpowered to estimate that rate, and it did not luck into a case. **293's reading stands — the two-construct term
is knockdown failure, not seed off-targeting** — and T7 can now be written as a measured negative rather than an
untestable claim.

## Artefacts

| path | what |
|---|---|
| `work/sources301.md` | every download: URL, licence/terms, bytes, sha256 (112.9 MB total) |
| `work/vcc2025_spacers.csv` | 348 guide ids → protospacer A/B, sgIDs, gene, transcript, Ensembl id, route |
| `work/promoters_cage_tss_hg19.fa.gz` | 18,387 CAGE-TSS ± 500 windows, GRCh37 (the accepted set) |
| `work/promoters_tss1000_grch38.fa.gz` | the rejected GRCh38 gene-TSS set (KA-1 22/68), kept as the record |
| `work/joung_orf_map.csv` | 3,369 codes → `TFORF####`, symbol, isoform, RefSeq/Gencode id, n_cells, evidence |
| `work/joung_map_meta.json`, `work/joung_map_work.json`, `work/joung_anchors_S01.json` | sublibrary↔P1 purities, anchors, resolved map |
| `work/seed_hits301_cage.csv`, `work/seed_scan301_cage.json` | 4,084 raw seed+PAM hits and the per-arm summary |
| `work/seed_contrast301.json` | the filtered contrast, both window definitions, all four seed lengths |
| `work/seed_L12_crispri_window.json` | the 18 L=12 CRISPRi-window hits with locus separation |
| `work/seed_{hits,scan}301_grch38_REJECTED.*` | the rejected first run's output, kept as the record |
| `work/{build_spacers,fetch_promoters,fetch_promoters_hg19,joung_map,seed_scan,seed_contrast}.py` | the code |
