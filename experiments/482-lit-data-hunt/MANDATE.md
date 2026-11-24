# 482 — Literature and data hunt (2025–2026): anything public that could raise a zero-shot cell-line score

You are a research teammate (Grok 4.7, extra-high reasoning) working in a read-only copy of the Virtual Cell Challenge
2026 repo. Read `AGENTS.md` and `docs/FINDINGS.md` §0–§2 first, and `docs/BOARD-INTEL-2026-09-26.md`. Write ONLY
`experiments/482-lit-data-hunt/RESULTS.md` (create the directory). Budget: up to ~6 hours. Save progress into RESULTS.md
as you go.

## The task

Arc Institute's Virtual Cell Challenge 2026: predict CRISPRi (dCas9-KRAB) Perturb-seq responses (10x Flex chemistry,
400 cells per target, 18,533 genes) for 300 target genes in human cancer cell lines we have never seen perturbed, given
only each line's non-targeting control cells. The final test (Oct 22) uses three new cancer lines and a new 300-target
panel. Validation lines are Jurkat, HeLa and CAL-33. Public perturbation data may be used to train models.

Find, with verified sources:
1. **New public single-cell CRISPR perturbation datasets** (Perturb-seq, CROP-seq, TAP-seq, Flex-based screens, etc.)
   in HUMAN CANCER cell lines, released or updated 2025–2026, especially genome-scale or ≥ 1,000 targets, CRISPRi,
   or 10x Flex chemistry, and especially lines NOT in the list below. Sources: GEO, ArrayExpress, Zenodo, figshare,
   Hugging Face (`hf datasets ls --search ... --limit 50`), CZ CELLxGENE, Arc Virtual Cell Atlas, scPerturb,
   bioRxiv/medRxiv/arXiv supplements, company data releases (Xaira, Illumina, Parse, Scale, Tahoe, 10x, CZI Billion
   Cells, Genentech, Recursion).
2. **How Arc chose the 2026 cell lines and targets**: the challenge paper (Cell 2026, "Virtual Cell Challenge 2026:
   Benchmarking zero-shot generalization across cellular contexts"), Arc blog posts ("Behind the Data"), talks, FAQs.
   Anything that narrows which lines D/E/F are, or how their targets were selected.
3. **Published methods of the 2025 Virtual Cell Challenge winners** and any 2026 write-ups of cross-cell-line
   perturbation prediction that beat simple transfer/linear baselines on held-out cell lines, with public code.
4. Work on predicting **which genes will be differentially expressed** under a perturbation, or conditional fold-change
   magnitudes, across cell lines.

Already known and on disk (skip unless there is a NEW release): Replogle 2022 K562 genome-wide and essential, RPE1;
Nadig 2025 HepG2 and Jurkat (essential); X-Atlas/Orion (Xaira) HCT116 and HEK293T; X-Atlas/Pisces (empty so far);
genome-scale CRISPRi in primary CD4+ T cells (2026); KOLF2.1J iPSC atlas (10,167 perturbations); Feng 2025 iPSC atlas;
Jiang 2025 (A549, HT29, MCF7, BxPC3, HAP1 pathway screens); Broad 16-line Cas9 knockout screen (Boyle 2026); VIPerturb
K562; Pan 2025 H1 hESC (GSE295214); VCC 2025 H1 data; Tahoe-100M (drugs); DLD-1 low-MOI; Myllia A549; Sunshine Calu-3;
Wessels THP-1; CMap/LINCS; IGVF cardiomyocyte, TeloHAEC, HCASMC screens; scPerturb collection; Arc State and Stack.

## Output (RESULTS.md)

- A ranked shortlist (top 10) with, for each: exact URL(s), cell line(s), perturbation modality (CRISPRi/KO/a),
  number of targets, cells, chemistry, release date, licence, total download size, and one line on why it could help
  (and which lines it would add to our same-line registry).
- The full inventory table below it.
- Mark each fact VERIFIED (you read it in the primary source) or REPORTED (secondary). No guesses presented as facts.

## Rules

- Read-only on the web: do NOT create accounts, post, comment, email, submit forms, open issues, or message anyone.
  Nothing about this repo or our method may appear anywhere public.
- Do NOT download any file larger than 50 MB (read metadata, file listings and sizes instead).
- Do NOT run `vcc`, `cell-eval2`, builds, git commits or pushes, and do not edit any file except your RESULTS.md.
