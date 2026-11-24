# 813 — Lane B4, step 1: inventory of public perturbation data the project lacks

**Answer (2026-10-16 06:40Z).** One new item adds a fourth held-out line covering every lab panel:
**IGVF TeloHAEC genome-wide CRISPRi Perturb-seq**, released on the IGVF portal on **2026-10-12**, four days ago.
- It targets 19,276 genes with 114,885 guides. The A/B/C panel is **300/300 in its library** and L1–L6 are **298–300/300**.
- In the resting rep-1 file alone, **298/300 panel targets have ≥ 100 cells** (median 353 cell-incidences).
- Cell-level h5ad files are public. No account is needed (licence: portal says "no certificate"; it is IGVF consortium data).
- Caveats: it is not a cancer line (hTERT-immortalised aortic endothelium), the design is MOI ~5 (median 4 guides per cell), and
  the data are shallow (median 5,452 UMI per cell).

The second find is cell-level **DLD-1 CRISPRi (GSE337988)**, a fourth cancer line. Production holds only its two low-MOI
LFC tables (2,827 gene rows). The deposit also has the processed cell matrices: 1,196,592 low-MOI QC cells × 36,601 genes.
But it overlaps each lab panel at only **65–82** targets, and DLD-1 is MSI colorectal, the same as HCT116.

**Release watches:** X-Atlas/Pisces is still empty, re-read at 06:19Z.

This is step 1 only. No data matrix was downloaded and nothing was spent on the cloud. Steps 2–3 wait for the lab lead's
GO.

## 0. Scope and hygiene
- 2026-10-16, 06:19–06:40Z. Web search, catalogue APIs (GEO E-utilities, the GEO FTP listing, the IGVF portal API, the Hugging
  Face API, the Zenodo API, the GitHub API) and small metadata files only. No `vcc`, no `cell-eval2`, no git, no build, nothing
  written outside `experiments/813-public-data/`. `window800.sh` printed **CLEAR** at 06:34Z.
- **Metadata range reads (declared, so the "no data download" line is exact).** Three helper scripts read HDF5 metadata from
  remote files through HTTP range requests, each with a hard byte cap. None reads an expression matrix:
  - `work/h5peek813.py`: file structure plus obs/var/X shapes;
  - `work/obsq813.py`: quantiles of three obs QC columns;
  - `work/cov813.py`: the guide-name list plus the guide-assignment index array.

  | target file | bytes fetched |
  |---|---:|
  | IGVF IGVFFI2522NCFR (66.3 GB) | 19.2 MB + 36.7 MB + 16.0 MB + 40.1 MB |
  | GSE337988 `*_assays.h5`, three files | ~2 MB each |

  Small catalogue files went to `work/meta/` (8 MB): the GEO sample/series text, guide manifests and 44 Cell Ranger
  `metrics_summary` files (387 B each). Every file's sha256 is in `work/SHA256SUMS-813.txt`.
- **Local compute:** sub-second aggregations only. The range-read runs were network-bound: under 1 min each and well under
  1 GB RSS.

## 1. Known-answer rows (read before any new number was believed)

| row | required | measured | verdict |
|---|---|---|---|
| KA-1 Pisces state vs production's watcher (`data/pisces-watch.log` last line 2026-10-16T03:20:38Z `data_files=0`) | no data files | HF API 06:19Z: files `.gitattributes, LICENSE.md, README.md`, usedStorage 0, lastModified 2026-03-25T22:58:18Z, cc-by-nc-sa-4.0 | PASS (unchanged since 270/306/413/482c) |
| KA-2 DLD-1 processed cell counts vs the published QC total (116 README: 1,869,846 cells after QC) | exact | Low 1,196,592 + Med 399,359 + High 273,895 = **1,869,846** (`meta/tree_GSE337988_*`) | PASS, exact |
| KA-3 DLD-1 sublibrary-2 ∩ A/B/C panel vs 413/482c | 67 | **67** (GEO guide manifest, 4,955 genes; and 116's promoter symbols, 4,715) | PASS |
| KA-4 TeloHAEC guide-assignment internal check | nnz = mean guides/cell × cells | 6,280,216 = 4.7010 × 1,335,923 (6,280,216.1) | PASS |
| KA-5 `cov813.py` (file) reproduces the first inline run | identical JSON values | 0 differing keys (`telohaec_rep1_coverage.json` vs `_v2.json`) | PASS |

## 2. Inventory

"Fold" means a held-out line for the 801 engine: a line kept out of every donor mix, whose own rows are the truth. "Donor"
means a new line for 811's model or for production. "Panel overlap" is against the A/B/C panel (`data/challenge/pert_counts.csv`)
and against the lab's L1–L6 draws (`experiments/801-lab-folds/work/L*/bundle/pert_counts.csv`). VERIFIED means read from the
record or the file today. EST means an estimate; its basis is stated.

### 2.1 New or newly usable (ranked)

| # | dataset | source / licence | line, state | modality, chemistry | targets | cells per target | depth | size | panel overlap | what it adds |
|---|---|---|---|---|---|---|---|---|---|---|
| **1** | **TeloHAEC genome-wide CC-Perturb-seq** (Engreitz lab, IGVF), released **2026-10-12** | IGVF portal. Merged: analysis set IGVFDS8166LJMC, file IGVFFI3974SAVD. Per condition × rep: IGVFDS2068OODO (resting r1, IGVFFI2522NCFR), IGVFDS6017INUS (resting r2), IGVFDS9554DAPM / IGVFDS4026PWBE (oxLDL r1/r2), IGVFDS5371SLCS / IGVFDS2589XTPD (IL-1β r1/r2). Downloads are presigned S3 (us-west-2) via `api.data.igvf.org/.../@@download/`. Data-use summary "no certificate". DOIs `10.65695/IGVFDS…`; no paper yet | TeloHAEC: hTERT-immortalised human aortic endothelium, female, **non-cancer**; 3 conditions (resting, oxLDL 7 d, IL-1β 7 d) × 2 biological reps | CRISPRi dCas9-KOX1-KRAB, lentiviral, **MOI ~5** (median 4 guides/cell, mean 4.70; 10 % carry 0); 3′ scRNA ("CC-Perturb-seq"), guide UMI counts and pipeline guide calls in `obsm` | **19,276 genes**, 114,885 guides (median 6/gene), 1,149 NTC guides; 1,960 K562-essential genes in the library (VERIFIED, `cov813.py`) | resting rep 1 alone: per gene q10/q50/q90 = **183 / 311 / 476** cell-incidences; A/B/C panel median 353 (VERIFIED). Rep 2 + 2 stimulated conditions not counted yet. EST: resting total ~550/gene (rep-2 file is 0.67× rep 1's size) | median **5,452 UMI**, 2,590 genes per cell (rep 1, VERIFIED) — ~¼ of Arc's ~20 k | resting r1 1,335,923 cells (66.3 GB); r2 44.5 GB; oxLDL 41.6 + 38.7; IL-1β 35.2 + 38.2; merged 199 GB; X float64 CSR, 62,757 genes | **A/B/C 300/300** in library, 298 ≥ 100 cells; **L1 298, L2 299, L3 299, L4 300, L5 298, L6 299** in library; ≥ 100 cells 295–299 per panel (VERIFIED, rep 1) | **A fourth fold line for every lab panel** (L1–L6 and A/B/C; L7/L8 at the lab lead's discretion), at near-full coverage. Not in any production donor mix (`docs/october-decisions.json` has no TeloHAEC; the on-disk `data/external/igvf-telohaec` is the older 2,167-target pilot IGVFFI8645YMGX, 214,449 cells, used only by 122/135). Also a new lineage for 811 (three conditions = a state axis). |
| **2** | **DLD-1 CRISPRi, cell level** (GSE337988; Genentech: Yeung, …, Mayba, Xie; "Joint analysis of multiply perturbed cells…", public 2026-07-10) | GEO, no licence stated. Processed: `GSE337988_sublib2_processed_objects_{Low,Med,High}_{assays.h5,se.rds}`; per channel: `filtered_feature_bc_matrix.h5` + `crispr.umi_correct.h5` + hashing | DLD-1 colorectal adenocarcinoma (MSI, like HCT116) | CRISPRi, dox/Shield-1-inducible Zim3-dCas9, 10x 3′ GEM-X v4 + sgRNA capture + MULTI-seq; arms MOI 0.3 / 3 / 10 | sublibrary 2: 4,955 genes, 6,451 promoter perturbations (VERIFIED, GEO guide manifest) | EST, low arm: 1,196,592 QC cells (VERIFIED) ÷ 6,451 promoters × 0.5–0.8 singly assigned ≈ **90–150 per promoter**, ~120–190 per gene. Guide calls sit in `se.rds` colData (not range-readable). 116's Low1/Low2 author-LFC split-half cosine median 0.070 (on the 2,827 DE rows) | median 6.5–10.7 k UMI per low-arm channel (44 `metrics_summary` files, VERIFIED) | low processed 4.3 GB + 0.64 GB rds (dense int32, chunks 1000×1000); low per-channel filtered h5 5.8 GB; whole RAW.tar 26.5 GB | A/B/C **67**; **L1 65, L2 69, L3 82, L4 81, L5 80, L6 73** (VERIFIED). Arc-like pool for DLD-1-specific panels: **2,421** (sublib2 ∩ axis ∩ donor union − K562-essential − all 6,824 prior panel targets); of those, K562-GWPS carries 666 | **A fourth cancer fold**, either thin (~75 targets per L panel, wide CIs) or as 2 new DLD-1-specific Arc-like panels. Not in the October donor mix; it is only a same-line registry entry (`scripts/def_sameline.py`). Independence caveat: DLD-1 and HCT116 are both MSI colorectal lines. |
| 3 | **KOLF2.1J atlas, cell level** (Mali lab, Nat Biotech 2026) | figshare+ 10.25452/figshare.plus.27261219, CC BY 4.0 (218); HF mirror `Boom5426/iPSC-CRISPRi-PerturbSeq` (updated 2026-10-01, 267.6 GB, `KOLF_Pan_Genome_QC_Filtered.h5ad` + Feng counts; mirror of a board team) | KOLF2.1J iPSC (pluripotent) | CRISPRi, 10x 3′ | 11,739 | ~200–215 mean (278) | > 5 k UMI | 250.6 GB (pan-genome 189.4) | 273 (278) | Turns the kolf2 fold from per-target z-scores into cells, so a **second DE surface** (six members from emitted cells) on a non-H1 line becomes possible. 218 said NO-GO on the Mac budget; on EC2 it is a ~$10–20 job. |
| 4 | GWCD4i primary CD4 T (Marson; GSE314342; CZI VCP) | MIT | primary CD4 T × 3 conditions | CRISPRi, **10x Flex** | ~11.5 k | ~500 | — | 1.8 TB cell level | 283–287 | Already 805's DE-CD4 surface. Listed for completeness; not this lane's. |

### 2.2 Checked, nothing to add (with the number that closes each)

| resource (checked today unless dated) | result |
|---|---|
| **X-Atlas/Pisces** (HF `Xaira-Therapeutics/X-Atlas-Pisces`) | 3 files, usedStorage 0, last modified 2026-03-25 (KA-1). X-Cell model repo: last push 2026-03-16. Production's watcher keeps polling. |
| **GEO sweep** (E-utilities, `gds`, human/any, published 2026-08-15 → 10-16, three Perturb-seq/CRISPRi/CROP-seq queries; `meta/geo_esummary_0815_1016.json`) | 35 series. The only cancer-line screen of size is GSE345058 (A549 Cas9 KO, 1,000 genes, known from 413/482c: median residual 0.58, 21/300). Everything else is enhancer/variant maps (CAD GSE347104, GSE328800, GSE347395), Cas12a DDR combinatorial (GSE343308), mouse/E. coli, ORF screens (GSE344535). |
| **IGVF portal, all Perturb-seq measurement sets** (`meta/igvf_perturbseq.json`, 7,724 sets) | Besides #1: Gersbach K562 CRISPRi of accessible regions (enhancers, 2026-07-27); Gersbach / Hon / Huangfu WTC11- or HUES8-derived TF Perturb-seq (~2,000 TFs; the TFP3 project `github.com/IGVF/tf_perturb_seq`; pluripotent-derived, TF-only); Norman Hs27 CRISPRa (TFs, 2026-08-24); Norman RPE1 multiome CRISPRi (ATAC readout, 2026-09-01); Quertermous HCASMC Parse Perturb-seq (on disk as `hcasmc-igvf`). None is a genome-scale cancer line. |
| **Arc**: HF org `arcinstitute`; `arc-virtual-cell-atlas` repo | PIE_* datasets (2026-10-05) are PIE's sources and splits, already read by 700. The atlas repo's last commit is 2026-08-28 (README): no new perturbation data. |
| **Hugging Face** dataset search (10 queries, 200 newest each; `meta/hf_search.json`) | New perturbation sets since 08-15 are derivatives of data on disk: `sunweiwei/ai4sci-virtual-cell-zeroshot-data` (Replogle/Nadig/Orion re-splits, 190 GB); `Emulated-Inc/sci-perturbation-response-new-gene` (K562-GWPS summaries); `windmirror/perturb_datasets` (Adamson etc.); `syntensor/cellsim-crispri-screens` (a benchmark harness). |
| **Zenodo** (6 queries, newest 50; `meta/zenodo_search.json`) | Nothing new and relevant. VIPerturb-seq record 18460279 unchanged since 2026-02-14 (K562 Flex; 413: median 47 cells per panel target). scPerturb v1.4 (2024-08-26) is still current. |
| **CZI** Billion Cells Project / Virtual Cells Platform | The only public CRISPR release is still the CD4 atlas (#4). |
| **Illumina Billion Cell Atlas** | Commercial alliance only (413). |
| **Tahoe genetic arm / Parse GigaLab** | No public CRISPR release. Parse's public GigaLab set is a cytokine screen; Tahoe's 300 M-cell project is unreleased. |
| **Broad 16-line Perturb-seq** (bioRxiv 2026.08.24.746802, figshare 33273600) | = `data/external/broad16line` (on disk; 100 genes, 0/300 panel, Cas9 KO). |
| **Krogan-lab MDA-MB-468** (`github.com/y-doctor/krogan-reports`, public reports) | Breast cancer line, CRISPRi, 113 genes × DMSO/vorinostat, ~10 k cells per sample. Too small; no deposited counts. |
| **Weissman / Replogle / Satija / Marson 2026** | Nothing beyond what is on disk (VIPerturb K562, GWCD4i). Nadig HepG2/Jurkat are essential-only (0/300). |

## 3. What wiring them would buy (input to step 2; no PROPOSAL.md until the GO)

- **TeloHAEC (#1) is the strongest checker addition.** Every BULK-W panel gets a fourth held line with ~300 targets. That makes
  24 working folds instead of 18, and BULK-ABC gets a fifth fold.
  - The line is not in any donor mix, so PLOLO needs no refit for it beyond excluding it from operator, key and composition
    fits. It was never in them.
  - It is the first fold line from a non-epithelial, non-haematopoietic, non-pluripotent lineage. That is the sort of line
    shift D/E/F may bring.
- **What has to be checked before it counts (step 3):**
  1. Knockdown: the target's own CPM in carriers vs non-carriers, per guide.
  2. Power: split-half cosine by guide halves, read against 278's F1 law. At ~5.5 k UMI and MOI 5, expect weaker per-target
     reliability than Orion.
  3. **615's blob test with the control definition stated.** In an MOI-5 design the "control" is either cells with 0 guides
     or NTC-only cells (rep 1 has 64,464 NTC incidences), or every non-carrier. Every carrier also carries ~3.7 other
     perturbations, so the target pseudobulk is shrunk toward the screen mean.
  4. Retrodict before trusting it, as any new fold must: it must reproduce the BULK-ABC known answers (cap 6 +, 702 −) in sign
     before it is read.
- **Rough cost of #1:** the resting reps only (111 GB), streamed from IGVF S3 (us-west-2) to a us-east-1 worker.
  - Per-guide pseudobulk sums by CSR row blocks with the `obsm/guide_assignment` CSR. No matrix is ever densified.
  - Instance: m7i.2xlarge-class with ~250 GB gp3, 3–4 h, about $3–6.
  - Adding the two stimulated conditions: about $10–15 in all.
  - Nothing comes to the Mac except per-target tables (< 1 GB).
- **DLD-1 (#2):** a fourth cancer line, but thin on the existing panels (65–82 targets each). Its best use is two DLD-1-specific
  arc800 panels drawn from its 2,421-target pool.
  - Cost: ~5 GB of processed objects; the colData sits in `se.rds` (read with `rdata` via `uv run --with`, as 270 did); a
    streamed pass over the 1000×1000-chunked dense int32 matrix. About 1–2 h on EC2, ≤ $5.
  - Its knockdown and power can be read from the same pass. 116's split-half cosine (0.070 on author LFCs) warns that power
    may be marginal.
- **KOLF2.1J cells (#3)** is a DE-side addition (a second, unworn DE surface besides DE-H1 and 805's DE-CD4). Lower priority:
  it is pluripotent, and it is a ~190 GB pull.
- **Donor side:** none of these is a production donor candidate today, so there is no DOSSIER.md. TeloHAEC is a new lineage
  for 811's leave-one-line-out model. As a production donor it would need its own held-out price (225/238: sources saturate).

## 4. What is left
- Step 2 (PROPOSAL.md for #1 and #2, possibly #3) and step 3 (EC2 fetch, checks, fold wiring) wait for the lab lead's GO
  after the judge. Cloud spend so far: **$0** (lane budget $50 before a GO).
- Unverified, to settle first on EC2:
  - TeloHAEC rep-2 and stimulated-condition cell counts;
  - the per-gene cells after a single-target or NTC-control rule;
  - DLD-1's per-promoter singly-assigned cell counts.

## 5. Files (sha256 in `work/SHA256SUMS-813.txt`)
- `work/h5peek813.py` `b7df91a1…`, `work/obsq813.py` `bdba9aca…`, `work/h5tree813.py` `ee147d2f…`, `work/cov813.py` `bb56c8ad…`
- `work/meta/telohaec_rep1_coverage_v2.json` `a89101bd…` (from `cov813.py`; same values as the first inline run, `a307fcfa…`)
- `work/meta/peek_IGVFFI2522NCFR.json` `6408d956…`, `work/meta/obsq_IGVFFI2522NCFR.json` `96dae243…`,
  `work/meta/telohaec_guide_names.txt` `ba0e0918…`
- `work/meta/tree_GSE337988_Low_assays.json` `687e20bb…`, `work/meta/tree_GSE337988_MedHigh_assays.json` `9b3c7444…`,
  `work/meta/GSE337988_{series_brief,gsm_brief,filelist}.txt`, guide manifest `…crispr_demux_index_map_prod.csv.gz` `debaf4f1…`
- `work/meta/geo_esummary_0815_1016.json`, `work/meta/igvf_*.json`, `work/meta/hf_search.json`, `work/meta/zenodo_search.json`

## 6. Data terms (R7 of VERDICT-1; read 2026-10-16 09:18Z, before any data byte)
- **Source:** the IGVF External Resource Sharing Policy, https://igvf.org/external-resource-sharing-policy/ ("Last updated
  22_01_06"). Page as fetched: `work/meta/igvf_external_resource_sharing_policy.html` (sha256 `f04a7c23…`; it embeds dynamic
  markup, so its bytes can change on a refetch). Extracted text: `work/meta/igvf_external_resource_sharing_policy.txt`
  (sha256 `fc17237bc9abf0a7a0906431cdff4ac73baeee21567810a090985120787ac793`).
- **Operative sentence (verbatim):** "External users may freely download, analyze and publish results using any Impact of
  Genomic Variation on Function (IGVF) resources without restrictions. This applies to all data, analyses, software and
  protocols, and includes no grace period for IGVF Consortium members…"
- **Reading:** our use is allowed: download, analysis as a held-out evaluation line, and, if it ever happens, a model input
  or a published result. Contacting the producing lab is optional. Citation is requested, in three forms: the
  Consortium's latest publications; an acknowledgement of the IGVF Consortium and the Engreitz lab; and the accession
  numbers (IGVFDS…/IGVFFI…). Those go into any write-up that uses the line. The portal's file record agrees: its
  data-use summary is "no certificate" (no controlled access).
- **Gate status:** PASS on this reading. Stage 0 re-checks that the operative sentence is still present before the first
  fetch. If it has gone or changed, STOP with no byte fetched.

## 7. Stage 0a: KA-SIM, KA-SG bounds, A/B choice (2026-10-16 09:37–09:58Z; GO of VERDICT-2, G1)
**Answer:** `freeze813.json` written. sha256 **`2f6e5169bcfacb859c9ece156e5b168558ffb99cc6c6647c4eb956b55e71ceba`**.
- S3: `s3://<your-s3-bucket>/output/lab/813/freeze813.json` (3,485 B).
- Local: `work/s0a/freeze813.json`.
- Estimator: **A** (B −0.020 lower in median cosine).
- KA-SG power verdict: "informative" (naive outside its bounds in 200/200 replicates).
- **KA-SIM: FAIL by the frozen rule.** The offset condition (ii) reads **0.200 > 0.10** on design 3 (the injected guide-load offset);
  design 1 reads 0.081 (pass).
- So stage-0 gate 4 fails as run. Per C3 nothing is re-tuned and nothing goes further. **Stage 0b is NOT to be launched on this
  freeze.**
- Below are three defects that I found in my own implementation when I read the numbers, written up before any real row exists.
  The lab lead and the judge decide what follows.

**Run.** `cloud/jobs-s0a.txt`, run `20261016T093745Z-s0a813`.
- r7i.xlarge on-demand, 0.328 h, **actual $0.094** (cap $0.32).
- DONE 0; job wall 17:54 (smoke + full); peak RSS 7.4 GB.
- `labcloud.sh status` after it shows no 813 instance. The one running Project=lab instance is de803-b9, another lane's.
- Read: metadata only, 180.6 MB of range reads of IGVFFI2522NCFR (66.3 GB): guide names, the guide assignment, n_guides,
  total_counts, pool and barcodes. X was never opened.
- Code: `sim813.py` `b2c8a20e…`, `kasg813.py` `b94bf49e…` (the frozen sha), `work/h5peek813.py` `b7df91a1…`.
- Inputs: pilot b_ctrl `cbc0f2cc…`, K562 delta `0c942228…`, guide names `cd8e5780…`, assignment indices `bec179d9…`.
- Stage-0 targets (100, seed 813): sha `93cc51a3…`. Axis (3,000 genes): sha `4a81b3e0…`. Halves by `sha1('813:'+barcode)[0]&1`.
- Files: `work/s0a/{freeze813.json, detail813.json (de1698ae…), job1.log (62d9a062…), *-smoke.json}`.

**Numbers (R = 200 for multi and inject, 50 for twin and B).**

| item | value | rule | verdict |
|---|---|---|---|
| model | 1,335,923 cells, 19,276 targets (9,579 with a K562 truth row), 269/300 panel targets with truth, remainder mass 0.737, median ‖δ*‖ 4.18 | — | — |
| (iii) noise-free code check (twin, k=1 vs zero-guide, expected sums) | min cos **1.0000000** | > 0.9999 | PASS |
| A median cos(δ̂, δ*), multi design | **0.765** | — | — |
| B (first-order deconvolution) median cos | 0.745 (B − A = −0.020) | B only if ≥ +0.02 | **A** |
| twin median cos | 0.460 | — | (see defect 2) |
| (i) A / twin | 1.66 | ≥ 0.9 | PASS (not like-for-like; defect 2) |
| (ii) ‖mean(A − δ*)‖ / median‖δ*‖ | design 1 **0.081**, design 3 **0.200** | ≤ 0.10 both | **FAIL** (design 3) |
| KA-SG β bounds (central 95 %, unbiased design) | [0.0968, 0.1046] | — | (see defect 3) |
| KA-SG dsh bounds | [7.96, 8.25] | — | (see defect 3) |
| power: naive (k ≥ 2 vs zero-guide) outside bounds under injection | 200/200 = 1.00 | ≥ 0.8 → informative | "informative" |
| A inside the bounds under injection | 0.045 | (diagnostic) | — |

**Three defects in what ran (my implementation; found on reading, not tuned).**
1. **Lumped strata.** `sim813.py` caps the stratum at n_guides ≥ 12 (KCAP = 12). PROPOSAL-2 §2 says
   strata = (n_guides_total, pool, rep), with no cap.
   - The fraction of carriers grows with k, so inside the lumped 12+ bin, carriers have more guides than non-carriers. Any
     effect that scales with guide load then leaks into A.
   - That is design 3's residual offset of 0.200. The injection is deliberately strong: 0.5 × median‖δ*‖ PER GUIDE, and the
     n_guides tail runs to 189.
   - The same mechanism would leak a real guide-load technical effect. So this is a real property of the capped estimator, and
     of any estimator that coarsens k.
   - The preregistered (uncapped) A was not what ran. I have not re-run it.
2. **The twin does not match incidences.** It keeps one random guide per cell, so each target has ~4.7× fewer carrier
   incidences than in the multi-guide design. PROPOSAL-2 specified the same incidences per target.
   - So (i)'s ratio of 1.66 compares unequal depths and does not answer its question.
   - A matched-incidence twin is analytic: one virtual single-guide cell per incidence, each with its cell's library size.
3. **KA-SG's statistics are dominated by low-count log bias, so their bounds do not transfer to real data.**
   - Single-guide rows have very few counts: ~10 % of cells, so ~10 cells or ~0.1 M UMIs per target, split in halves. On a
     5-CPM gene that is well under one expected count per half.
   - `log1p(TS·x/N)` of such counts carries a large deterministic (Jensen) bias. It is shared across targets and across both
     halves.
   - So on an UNBIASED design, β centres at **0.10** (not 1) and dsh at **8.1** (not 0). Their narrow bounds encode the
     simulation's own count regime and effect sizes.
   - On real data, with different effect sizes and expression, an unbiased A would fall outside them. The "informative" power
     verdict is technically true (naive moves far more) but not usable.
   - **G5 consequence, stated:** KA-SG as frozen cannot serve as a real-data gate. If stage 0 proceeds without a valid KA-SG,
     only G0-BLOB guards against k-dependent offsets on real data, and the folds must carry that caveat in their registry entry.
   - A usable form would compute both rows on genes with enough expected single-guide counts per half (or in count space, not
     log1p), and calibrate on the matched-incidence twin. That is a new statistic, so it needs a new freeze.

**What I recommend (for the lab lead; nothing done).** Either:
- a corrected 0a (≤ $0.15): uncapped strata as preregistered, the matched-incidence twin, and KA-SG either redefined in
  count space on well-covered genes or dropped (then G5's caveat applies). It would be frozen as a new `freeze813-v2.json`,
  with the defects above recorded as the reason, judged as a revision; or
- stop here and record the negative: "a multi-guide stratum estimator that coarsens guide load leaks load-dependent effects,
  0.20 of a median effect per the injected test".

I lean to the corrected 0a. Defects 1 and 2 are implementation departures from the preregistered design, not results about
the data, and the data are untouched. Stopped; no stage-0b launch.

**Citation (G6, IGVF terms):** IGVF Consortium; Engreitz lab (Stanford); TeloHAEC genome-wide CC-Perturb-seq, analysis
set IGVFDS2068OODO, file IGVFFI2522NCFR (resting rep 1, metadata only so far).

## 8. Stage 0a v2 (RULING-0a; ONE corrected run, no third). Written before the run, 2026-10-16T10:02:34Z
- **Negative finding recorded (RULING-0a item 2): coarsening guide count leaks load-dependent offsets.**
  - A stratum-matched estimator that lumps n_guides ≥ 12 into one stratum (v1, freeze `2f6e5169…ceba`) left a shared offset of
    **0.200 × median‖δ*‖** under an injected per-guide load effect. The uncoarsened design reads 0.081 with no injection.
  - Mechanism: a cell's chance of carrying a given target grows with its guide count, so inside a lumped bin carriers have
    more guides than non-carriers.
  - It binds every later estimator in this lane: no coarsened k.
- **Code diff v1 → v2** (`work/s0a/sim813-v1-to-v2.diff`, sha256 `7d825fca29dfa884a7e502bbfbfaedf5c112479679e3cc29044c890330d5704f`;
  v1 kept as `work/s0a/sim813-v1.py` `b2c8a20e…`). **v2 `sim813.py` sha256 `270b87abc3d28060d0030409c1ae86706ca6ccf3f67e4fc737d534b0c9465c00`.**
  The diff is limited to the following:
  - (a) **Uncapped strata:** n_s = (max n_guides + 1) × n_pool, kbin = n_guides_total (PROPOSAL-2 §2). B's stratum
    composition matrix becomes sparse because uncapped strata make it 19 k × 1,330. This is the same arithmetic, needed to
    run (a) inside the budget.
  - (b) **The matched-incidence twin:** one virtual single-guide cell per (cell, target) incidence, carrying that cell's
    library size, pool and half. Cells with no target guide stay zero-guide cells. Each target therefore has exactly its
    multi-guide incidence count (PROPOSAL-2).
  - The freeze output per RULING-0a item 3: `"kasg": {"status": "uninformative (G5)"}` with no bounds. The raw KA-SG numbers
    stay only in the detail file, and nothing reads them. The file names are now `freeze813-v2.json` / `detail813-v2.json`.
  - Unchanged: seeds, axis, stage-0 targets, R (200 / 200 / 50 / 50), the injection, estimators A/B/S/naive, `kasg813.py`
    `b94bf49e…`.
- **Bars unchanged:** (i) A / twin ≥ 0.9; (ii) offset ≤ 0.10 on designs 1 and 3; (iii) noise-free min cos > 0.9999; B only
  if ≥ +0.02.
- **KA-SG (G5): uninformative at this depth.**
  - Consequence: on real data only G0-BLOB guards against k-dependent offsets.
  - The TeloHAEC folds carry that caveat in their registry entry.
- **If v2 fails any bar:** recorded and STOP.
- **Run:** r7i.xlarge on-demand, metadata range reads only, smoke run first. Expected ~25 min, ≈ $0.11. The launcher's
  minutes cap is 40 ($0.22 cap); a run past ~34 min would exceed the ruling's $0.15, and I will report it if so.

### 8.1 Result of v2 (run 2026-10-16 10:02–10:30Z): KA-SIM FAILS bar (ii) again. Stage 0 ends here (no third 0a)
- **`freeze813-v2.json` sha256 `ed11ef339102d95bd5dbcaae44f56fcfcc3279539be3c5cb7cdcf6b33e9cc5a3`.**
  - S3: `s3://<your-s3-bucket>/output/lab/813/freeze813-v2.json` (3,539 B).
  - Local: `work/s0a/freeze813-v2.json`. Detail `detail813-v2.json` `93c17df7…`; log `job1-v2.log` `aeef0b19…`.
  - Code `sim813.py` `270b87ab…` (as written above), `kasg813.py` `b94bf49e…`.
- **Run:** `20261016T100242Z-s0a813v2`, r7i.xlarge on-demand, 0.449 h, **actual $0.129** (≤ the ruling's $0.15).
  - DONE 0; job wall 24:45; peak RSS 8.0 GB.
  - Metadata range reads only (172 MB of 66.3 GB); X never opened.
  - `labcloud.sh status` after it shows no 813 instance; the only running Project=lab instance is de803-b9, not mine.

| bar | v1 (capped, unequal twin) | **v2 (preregistered)** | rule | v2 verdict |
|---|---|---|---|---|
| (i) A / matched-incidence twin, median cos | 1.66 (not like-for-like) | **0.971** (A 0.763, twin 0.786) | ≥ 0.9 | PASS |
| (ii) shared offset, design 1 (no injection) | 0.081 | **0.084** | ≤ 0.10 | PASS |
| (ii) shared offset, design 3 (injected guide load) | 0.200 | **0.143** | ≤ 0.10 | **FAIL** |
| (iii) noise-free min cos | 1.0000000 | **1.0000000** | > 0.9999 | PASS |
| A vs B | B −0.020 | B 0.743 vs A 0.763 (−0.019) | B only if ≥ +0.02 | A |
| KA-SG | bounds (invalid; defect 3) | **status "uninformative (G5)"**, no bounds | — | G5 caveat |
| strata | 13 × 7 (capped) | 190 × 7 = 1,330 (uncapped) | — | — |

**Reading.**
- Uncapping the strata cut the load leak from 0.200 to 0.143, but not under 0.10.
- With no injection, A is fine: cos 0.763, within 3 % of an ideal single-guide experiment with the same incidences.
- Its offset of 0.084 is the expected −μ̄ (one random guide's mean effect), which is shared by every target.
- So the estimator recovers per-target direction but does not hold a strong load-dependent effect out of the target rows.

**Probable mechanism (a hypothesis; not tested, since a third 0a is not allowed).**
- In the uncapped design, high-guide-count strata (the tail reaches 189) hold few cells, and carriers are enriched there.
- Where a (stratum, half) cell group has carriers but no non-carriers, its carriers stay in P_t with no matching background
  (w = 0).
- Where non-carriers are few, the background mixture is noisy.
- The injection scales as exp(k · τ), with ‖τ‖ = 0.5 × median‖δ*‖ per guide. So exactly those high-k carriers carry the
  largest load effect, and they leak it into P_t.
- The bar tests that extreme regime by design. A real load effect is unknown and could be much smaller, but the frozen bar
  is what binds.

**Consequences (binding, recorded):**
1. **Stage 0 stops at 0a. Stage 0b is not launched, and no TeloHAEC expression byte is fetched.** The multi-guide truth line
   is NOT validated by its preregistered simulation gate. Under the bars as frozen, TeloHAEC genome-wide does not become a
   lab fold line.
2. **Negative findings, for the estimator class (they bind any later attempt here):**
   - coarsened k (v1) leaks 0.200;
   - uncoarsened stratum matching (v2) still leaks 0.143 × median effect under a per-guide load of 0.5 × median effect;
   - first-order co-occurrence deconvolution (B) neither helps direction (−0.02 cos) nor removes the offset (0.084).
   - A design-3 pass would need an estimator that also handles the sparse high-k tail. That would be a new construction,
     so a new proposal.
3. **What did pass, and stays useful:**
   - per-target direction recovery at matched incidences (0.971 of a single-guide experiment) without a load effect;
   - the noise-free code identity;
   - the terms gate (§6).
   - If someone later re-proposes this line with a k-tail rule fixed in advance (e.g. drop carriers whose stratum-half has
     fewer than N non-carriers, preregistered), the v2 simulation is the ready known answer. That is for the judge, not
     this run.
4. **Lane spend:** stage 0a v1 $0.094 + v2 $0.129 = **$0.223** of the stage-0 cap of $5. Nothing else was launched.

**Citation (G6):** IGVF Consortium; Engreitz lab (Stanford); TeloHAEC genome-wide CC-Perturb-seq, analysis set
IGVFDS2068OODO, file IGVFFI2522NCFR (metadata only).
