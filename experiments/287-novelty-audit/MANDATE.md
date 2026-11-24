# 287 — Novelty audit: every finding F1–F12 (and the C1 / C2 readings) against the ORIGINAL dataset papers and the perturbation-transfer literature — what the sources already report, what is new here, what they contradict. Reading and citation only, ≤ 3 hours. No compute beyond small checks, no build, no submission.

You are a research subagent. Work in `experiments/287-novelty-audit/`. **Edit nothing outside it** (the coordinator
integrates your audit into `docs/FINDINGS-2026-09.md`). You MAY use web search / web fetch for papers, preprints, GEO /
CZI dataset pages and supplementary tables; write what you fetch (URL, DOI, date accessed) into `work/refs.md`. No
`vcc`, no `cell-eval2`, nothing under `data/` or `submissions/`; `uv run python` only for small checks (e.g. reading a
supplementary table). **Time box: 3 hours.** Three compute experiments (284, 285, 286) run in parallel — do not read or
write their directories, and do not start heavy jobs.

## Why (third external review, 2026-09-06, item C.5)
"Compare each claimed finding against the original dataset papers before calling it novel." The programme's scientific
deliverable is a write-up of which perturbation responses are conserved across cell states, which depend on the
recipient, and what a reference can resolve. Before any of it is called a contribution, each finding needs the prior
report it extends or repeats, cited, and any contradiction stated.

## Inputs on disk
`docs/FINDINGS-2026-09.md` (F1–F12 with the 13:40Z amendments — read the amended wording, not the original headlines);
`docs/READ-FIRST.md` §6; `docs/CONFIRMATION-LEDGER.md`; the RESULTS.md of 267, 268, 271, 272, 273, 276, 279, 280, 281,
282, 283 (their §"what this establishes" / "limitations" sections state each claim's scope); on-disk texts:
`data/external/xcell/xcell_biorxiv_full.txt`, `xcell_preprint.txt`, `pisces_README.md` (X-Atlas / Orion, Pisces);
`data/external/cd4_cell2026/data_sharing_readme.md` and `experiments/280-gwcd4i-donor-replicates/work/meta/` (GWCD4i
readme, sample metadata, guide-efficiency and QC tables); `data/external/gse295214/` (Pan 2025 hESC; `media-*.xlsx`
supplements); `experiments/268-broad16-context-residuals/work/inventory.json` and `data/external/broad16line/metadata`
(the 16-line set's provenance); `experiments/270-data-reconnaissance/RESULTS.md` (the dataset map).

## Papers to check (fetch what is not on disk)
GWCD4i (Zhu, Dann, … Marson; GEO GSE314342; the genome-scale CD4 T-cell CRISPRi paper) — its OWN numbers on donor
reproducibility (`donor_correlation_all_mean`), guide reproducibility (`guide_correlation_all`), activation-state
dependence, and what it says about construct effects; X-Atlas / Orion (Xaira; HCT116 / HEK293T) — replicate and batch
reproducibility, cross-line conservation; Replogle 2022 (K562 / RPE1 genome-wide Perturb-seq) — split-guide
reproducibility, K562-vs-RPE1 conservation of responses, the essential-vs-non-essential contrast; Nadig 2025
(Jurkat / HepG2, "transcriptome-wide characterization of genetic perturbations") — shared vs cell-type-specific
responses; the Broad 16-line Perturb-seq paper (identify it from 268's metadata); scPerturb (Peidli 2024); Pan 2025 hESC
(GSE295214); the VCC 2025 / 2026 challenge papers and the published 2025 leaderboard analyses; the transfer / baseline
literature (Ahlmann-Eltze 2024–2025 "deep learning-based predictions of gene perturbation effects do not yet outperform
simple linear baselines", Wu 2024 PerturBench, Csendes 2025, Bendidi 2024, GEARS, scGPT, STATE / Arc Institute); and for
286's estimator, the precedents for cross-validated / cross-product covariance and PCA denoising (cvPCA — Stringer 2019;
Gavish–Donoho; split-half PCA; Marchenko–Pastur noise floors in single-cell).

## Protocol
For EACH finding F1–F12, C1, C2 (and F3's "the whitener learns its own noise", F4's `pds` knee, 280's donor × guide
2×2, 282's construct effect, 279/281's family axis):
1. **Our claim** as stated after the amendments, with its number, unit and data.
2. **Prior reports** on the same question: paper, section / figure / table, and a quotation of at most 15 words in
   quotation marks with attribution (copyright: never more; paraphrase the rest).
3. **Verdict:** NEW (no prior report) / KNOWN (already reported, cite) / EXTENDED (known qualitatively; our number, unit
   or dataset is new) / CONTRADICTED (which prior result, on what data, and why the two disagree — depth, statistic,
   estimator, population).
4. **What would settle a contradiction** (a measurement on disk, or not available).
Be adversarial: prefer "KNOWN" when a paper reports it, and say when our statistic is not comparable to theirs.

## Deliverable
`NOVELTY-AUDIT.md`: the table (finding · our claim · prior reports · verdict · notes), the reference list with URLs /
DOIs and access dates (`work/refs.md`), a section "what the write-up can claim as new, and with what caveats", a section
"contradictions to resolve before publication", and one line: **NEW: n, EXTENDED: n, KNOWN: n, CONTRADICTED: n (list)**.
No edits outside this directory. Then stop.
