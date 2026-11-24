# 356 — Data scout 2026: nothing new to acquire today; two items to watch, two zero-cost held assets to analyse

Status: **COMPLETE** (web research + local inventory reads only; no downloads, no builders, no submissions).

## Known-answer reproduction (rule 3, own code)

Re-read `experiments/323-current-family-validation/work/metrics-cap15.json` (`raw` block):

| member | re-read value | record (343) | match |
|---|---|---|---|
| nmae | 0.95083646 | 0.950836 | exact |
| pds | 0.88963585 | — (info) | — |
| reach | 0.25434726 | — (info) | — |
| mse | 0.97494624 | — (info) | — |
| fid | 0.49755452 | — (info) | — |
| jac | 0.07475983 | — (info) | — |

NMAE matches 343's headline to 6 decimals. Proceeding.

## What was measured

Web search/fetch (all URLs accessed 2026-09-10) for 2026 single-cell perturbation releases,
Jurkat/HeLa/CAL33 perturb transcriptomes, DepMap 2026 + line embeddings, and open-weights
transfer models; cross-checked every candidate against the local `data/external/` inventory
(35 entries) and small local reads (cmap `siginfo_beta.txt`, depmap26q1 metadata, dir listings).
Memory: JSON/CSV reads only, peak RSS < 1 GiB. No file over 2 GB touched or downloaded.

## §1 — Single-cell perturbation datasets released/updated in 2026

| candidate | facts (URL, date, size, license, format, cells, perts, lines, modality) | held? | verdict |
|---|---|---|---|
| X-Atlas/Pisces (Xaira) | Preprint bioRxiv 10.64898/2026.03.18.712807, posted 2026-03-20; 25.6M perturbed cells, 7 genome-scale CRISPRi screens (HCT116, HEK293T, HepG2, iPSC, Jurkat Resting, Jurkat Active, iPSC Multi-Diff), median 25,478 UMIs, 78.7% KD; license CC-BY-NC-4.0 (preprint) / CC-BY-NC-SA-4.0 (HF). HF `Xaira-Therapeutics/X-Atlas-Pisces` verified 2026-09-10: **24.3 kB total, "Coming Soon"** — full data NOT uploaded | orion (= Pisces HCT116/HEK293T) held; rest not | WATCH — highest future value, not acquirable today |
| Zhu et al., Cell 2026 (Aug 28; bioRxiv Dec 2025) | 22M primary CD4+ T cells, 4 donors × Rest/Stim8h/Stim48h, genome-wide CRISPRi, probe-based; **MIT license**; CZI Virtual Cells Platform v1.0 (processed, 2025-12-22) + cell-level h5ad + pseudobulk + DESeq2 DE stats | **HELD** (`cd4t/` + `cd4_cell2026/` suppl tables verified on disk) | not new — confirms held bank is the current SOTA CD4 atlas |
| Moonen et al., bioRxiv 2026.03.09.710372 | 4.1M primary CD4+ T cells, genome-wide + 1,032 CREs; CC-BY-4.0 | not held | secondary: same lineage as held Zhu; only if a second CD4 donor set is needed |
| You et al. 100M-cell genome-wide Perturb-seq atlas (P588 abstract, 2026) | CRISPRi + CRISPRa + shRNA, "dozens of cell types", Illumina emulsion kit; conference abstract only (sciencedirect P588) | no | no public data found; note §5 flag re leader |
| Tahoe–Arc–Biohub 120M-cell / 225k-interaction set (announced 2026-01-12) | **drug** perturbations (Mosaic), open-source promised, first shared internally | tahoe (100M) held | not released; chemical, not genetic |
| Tahoe 300M GigaLab project (announced 2025-12-11) | upcoming | no | not released |
| Parse–Graph immune atlas (announced 2026-01-20) | planned, hundreds of millions of patient immune cells | no | not released |
| PerturbAI 8M mouse-brain in vivo atlas (2026-03-17; HF + NVIDIA blueprint) | mouse neuronal, thousands of perts | no | wrong species/tissue; no value for A/B/C |
| Perturb-DBiT spatial (Nat Biotech, 2026-06-11) | mouse tumors, ≤288 sgRNAs | no | wrong modality/scale |
| VIPerturb-seq (Bradu/Satija; K562 + HAP1 demo) | platform paper, ~12k-cell demo | K562 held; HAP1 in held jiang2025 | no new line |
| Feng iPSC atlas, Cell Genomics (held `feng2025`/`feng-ipsc`); Nat Biotech 2026 iPSC atlas | same programme | **HELD** | not new |
| scPerturb / PerturBase (NAR 2025, 122 sets/5M cells) | aggregators of public GEO data | scperturb held | mining source only; no new line identified |

## §2 — Jurkat / HeLa / CAL33 specifically

| line | finding | source |
|---|---|---|
| Jurkat (A) | GSE249595 genome-wide CRISPRi (rest+stim) — **held** (`gse249595/`) | local inventory |
| Jurkat (A) | Pisces Jurkat Resting + Jurkat Active full screens — **not yet public** (see §1) | HF verified 2026-09-10 |
| Jurkat (A) | **256 JURKAT shRNA knockdown signatures already on disk** in held cmap (`siginfo_beta.txt`: 15,197 JURKAT sigs = 14,079 compound + 846 vehicle + 256 `trt_sh`, 256 distinct hairpins; measured 2026-09-10 with own code) | local, NEWLY SURFACED — 1 hairpin/gene, 978 landmarks, pilot-scale; zero acquisition cost |
| HeLa (B) | 26,170 HELA sigs in held cmap — **all compound** (`trt_cp` 24,771 + vehicle); zero knockdown | local, measured 2026-09-10 |
| HeLa (B) | no HeLa Perturb-seq with transcriptome readout found anywhere (2026 or earlier) | web, negative |
| CAL33 (C) | DepMap member (ACH-000518, Cellosaurus CVCL_1108); only viability CRISPR (Allevato Brunello screen, no transcriptome) + observational HNSCC patient scRNA (GSE234933, GSE103322) | web |
| CAL33 (C) | **no CAL33 perturbation transcriptome of any modality found** (sc or bulk KD) | web, negative |

## §3 — Line covariates: DepMap, CCLE, embeddings

- **DepMap latest = 26Q1** (announced 2026-04-01; +25 models, Olink proteomics 161 lines). No 26Q2 data release exists (only a schedule Q&A thread). 26Q1 is **held** (`depmap26q1/`), but the local extract is a 16-line subset (apollo31/32 batch), so the three context lines must be pulled from the portal or the HF mirror (`ChanghaoKan/crispr-depmap`: CRISPRGeneEffect_26Q1 441 MB + Model.csv 697 kB; license undeclared on the mirror — prefer the portal, CC-licensed).
- **All three lines are DepMap/CCLE members**: CAL-33 = **ACH-000518** (Cellosaurus, verified); HELA_CERVIX = **ACH-001086** (DepMap portal URL, verified); JURKAT present in CCLE expression + DepMap CRISPR (Harmonizome, verified; exact ACH to be read from Model.csv, not guessed here). CCLE expression held (`ccle_broad_2025/`, 312 MB TPM matrix).
- **2026 "cell line embedding" resources**: none new found. The live prior stack (ESM-2, STRING, GenePT, DepMap, JUMP-Cell Painting, scGPT) is exactly X-Cell's cross-attention set — all public/held (GenePT vectors quarantined per 333).

## §4 — Open-weights transfer models (2026)

| model | weights? | license | inference hw | verdict |
|---|---|---|---|---|
| X-Cell (Xaira, Mar 2026; 55M Mini → 4.9B Ultra; zero-shot line transfer claimed) | **NO** — HF repo verified 2026-09-10 via API: README + PNG only (273 kB), "coming soon"; Ultra never promised openly | would be CC-BY-NC-SA-4.0 | Mini: 8 GB VRAM | cannot acquire; Mini was trained on Replogle–Nadig (public) only, so even release adds no Pisces information at the small end |
| STATE v1.0 (Arc, Cell) | held (`state-pretrained/`, `state-vcc/`); local `jurkat_0.99/` + `hepg2_0.99/` prediction dirs exist | non-commercial | — | held; rules note it is "not able to win" |
| Tahoe-x1, TranscriptFormer | open weights not verified | — | — | not verified |

## Ranked "what to acquire" (≤ 5)

1. **X-Atlas/Pisces Jurkat Resting/Active screens (watch HF repo).** Same-line genome-wide CRISPRi for context A; HepG2/iPSC add new lines for October. Expected value: highest of any known item. Action: re-check monthly; no download possible today.
2. **JURKAT shRNA L1000, 256 signatures (held cmap — analyse, don't download).** Zero-cost same-line KD labels for context A. Caveats: 978 landmarks, 1 hairpin/gene. Suggested follow-up experiment (not this one): intersect the 256 genes with the 300-target panel and test sign concordance vs our design effects. No license issue (public LINCS, already on disk).
3. **Full DepMap 26Q1 Model.csv + CRISPR/expression rows for ACH-000518 / ACH-001086 / JURKAT (portal download, < 1 GB, license-clean).** Gives verified IDs + dependency/expression covariates for all three context lines. Small, clean, unblocks any covariate arm.
4. **Moonen 4.1M CD4 atlas (bioRxiv Mar 2026, CC-BY-4.0).** Second CD4 donor set adjacent to Jurkat lineage; complements held Zhu. Caution: likely multi-GB — size must be confirmed in PLAN before any pull (none pulled here).
5. **LINCS CRISPR-KO expansion audit (held cmap + L2S2 catalogue).** Marino et al. 2025 report 7,508 KO genes × ≥10 lines / 1.678M signatures; verify whether the held cmap build contains the KO subset and whether HELA/JURKAT have KO (not just compound/shRNA) coverage. Analysis-only.

## Licensing red flags

- Xaira CC-BY-NC-SA-4.0: NC is fine for the challenge (rules §7 permit proprietary data with rights; training use is allowed), but **ShareAlike** could encumber derived artifacts — check before ingesting if it ever releases.
- Zhu CD4 MIT: clean. L1000/GEO/DepMap-portal: clean. HF mirror datasets with undeclared licenses (e.g. `ChanghaoKan/crispr-depmap`): prefer the primary portal.

## Verdict against the PLAN kill criterion

The kill criterion fired in full: no candidate offers a *currently downloadable*, license-clean,
new-line genome-wide CRISPR perturbation transcriptome, same-line A/B/C labels beyond what is
held, or open transfer-model weights. Verdict: **NEGATIVE — nothing to acquire today.**
The two watch items (Pisces Jurkat, Tahoe–Arc–Biohub release) and two zero-cost held analyses
(JURKAT shRNA-256, DepMap three-line pull) are the complete output. No download was made
(none under 2 GB and license-clean was needed).

## Limits

Web-keyed: HF "Coming Soon" pages and press releases change; Pisces/X-Cell status is as of
2026-09-10. H1-panel numbers are not quoted here (no computation run). JURKAT shRNA-256 utility
is untested (978 landmarks, single hairpins). DepMap ACH for JURKAT not verified to the digit.
CAL33 "no perturb transcriptome" is a search negative, not a proof.

## Exact commands to reproduce

```bash
# known answer
python3 -c "import json;m=json.load(open('experiments/323-current-family-validation/work/metrics-cap15.json'));print(m['raw'])"
# JURKAT/HELA cmap audit
python3 -c "import pandas as pd;s=pd.read_csv('data/external/cmap/siginfo_beta.txt',sep='\t',usecols=['cell_iname','pert_type','pert_id'],low_memory=True);print(s[s.cell_iname.isin(['HELA','JURKAT'])].value_counts(['cell_iname','pert_type']))"
# depmap subset inventory
python3 -c "import csv;print([r['cell_line'] for r in csv.DictReader(open('data/external/depmap26q1/extracted/metadata/cell_line_metadata.csv'))])"
# X-Cell weights check: GET https://huggingface.co/api/models/Xaira-Therapeutics/X-Cell (siblings == [.gitattributes, README.md, png] => no weights)
```

Peak memory < 1 GiB; wall time ~1 h (web research dominated).
