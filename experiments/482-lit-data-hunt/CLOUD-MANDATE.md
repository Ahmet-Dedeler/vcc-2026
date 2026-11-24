# 482c — Literature and data hunt with verification (Cursor Cloud VM, Grok 4.7 extra-high)

You are a research teammate running on your own cloud VM, working on a branch of the private Virtual Cell Challenge 2026
repo. Read `AGENTS.md` and `docs/FINDINGS.md` §0–§2 first. Write ONLY inside `experiments/482c-cloud-lit-data/`
(create it): `RESULTS.md` plus any small scripts you used. Commit only that directory to your branch; do not open a PR,
do not touch other files, do not merge anything. Budget: up to ~8 hours. Commit progress to RESULTS.md as you go.

## The task

Arc Institute's Virtual Cell Challenge 2026: predict CRISPRi (dCas9-KRAB) Perturb-seq responses (10x Flex chemistry,
400 cells per target, 18,533 genes) for 300 target genes in human cancer cell lines never seen perturbed, given only each
line's non-targeting control cells. The final test (Oct 22) uses three new cancer lines and a new 300-target panel.
Validation lines are Jurkat, HeLa and CAL-33. Public perturbation data may be used to train models, and a public screen
of the SAME cell line as a test line is our single largest lever (it would switch on our "same-line" mode).

Find, and VERIFY by opening the data (you may download to the VM, up to ~30 GB in total):
1. **Public single-cell CRISPR perturbation datasets in human cancer lines** (Perturb-seq, CROP-seq, TAP-seq, Flex-based
   screens), especially 2025–2026 releases, genome-scale or ≥ 1,000 targets, CRISPRi, and lines NOT listed below.
   Search GEO, ArrayExpress, Zenodo, figshare, Hugging Face, CZ CELLxGENE, Arc Virtual Cell Atlas, scPerturb, bioRxiv /
   arXiv supplements, company releases (Xaira, Illumina, Parse, Scale, Tahoe, 10x, CZI Billion Cells, Genentech).
   For each top candidate, open the files and report: cell line(s), modality, number of targets, cells per target,
   controls, chemistry, and a quality check: the fraction of targets whose own gene is knocked down in its cells
   (on-target knockdown rate).
2. **How Arc chose the 2026 lines and targets**: the challenge paper (Cell 2026, "Virtual Cell Challenge 2026:
   Benchmarking zero-shot generalization across cellular contexts"), Arc blog posts, talks, FAQs. Anything that narrows
   which lines the final three are.
3. **Methods of the 2025 Virtual Cell Challenge winners** and 2025–2026 cross-cell-line perturbation-prediction methods
   with public code that beat simple transfer or linear baselines on held-out cell lines.

Already known and on disk (skip unless there is a NEW release): Replogle 2022 K562 genome-wide and essential, RPE1;
Nadig 2025 HepG2 and Jurkat (essential); X-Atlas/Orion (Xaira) HCT116 and HEK293T; X-Atlas/Pisces (empty so far);
genome-scale CRISPRi in primary CD4+ T cells (2026); KOLF2.1J iPSC atlas; Feng 2025 iPSC atlas; Jiang 2025 (A549, HT29,
MCF7, BxPC3, HAP1); Broad 16-line Cas9 knockout screen (Boyle 2026); VIPerturb K562; Pan 2025 H1 hESC (GSE295214);
VCC 2025 H1 data; Tahoe-100M (drugs); DLD-1 low-MOI; Myllia A549; Sunshine Calu-3; Wessels THP-1; CMap/LINCS; IGVF
cardiomyocyte, TeloHAEC and HCASMC screens; scPerturb collection; Arc State and Stack.

## Output (`experiments/482c-cloud-lit-data/RESULTS.md`)

- A ranked shortlist (top 10): exact URL(s), line(s), modality, targets, cells, chemistry, release date, licence, total
  size, the verification you did (what you opened, the knockdown rate), and one line on why it could help.
- The full inventory below it. Mark each fact VERIFIED (you opened the data or read the primary source) or REPORTED.

## Rules

- Read-only on the web: do NOT create accounts, post, comment, email, submit forms, open issues or PRs, or message
  anyone. Nothing about this repo or its method may appear anywhere public.
- Do not commit data files; commit only RESULTS.md and small scripts (< 1 MB) in your directory.
- Do not run `vcc` or `cell-eval2`; do not edit files outside your directory.
