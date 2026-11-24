You are a research agent for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md
first, then docs/FINDINGS.md section 0-1 and the same-line mode notes (search FINDINGS and docs/DEF-RUNBOOK.md for
"same-line", and scripts/def_sameline.py's registry). This VM has NO gitignored data; you may browse the web and download
small public files (< 1 GB each, < 5 GB total) into experiments/520-cloud-def-lines/work/.

Experiment 520. Write ONLY under experiments/520-cloud-def-lines/. Commit to your branch only. Never run `vcc`.

Question: which cell lines are the final test contexts D, E, F likely to be, and what public single-cell perturbation data
(CRISPRi preferred; CRISPR KO / CRISPRa noted) exists for each candidate that our same-line registry does not already use?

Known facts: Arc ran CRISPRi Perturb-seq (10x Flex) on six human cancer cell lines "originating from different tissues";
validation contexts A/B/C are Jurkat (T-ALL), HeLa (cervix), CAL-33 (head and neck SCC) by CCLE match; D/E/F are three
different lines, released Oct 22 with their own 300-gene panel (targets drawn from an H1 pilot screen, >80 % knockdown,
high basal target expression in the lines, "adequate coverage in public datasets"). The final score depends heavily on
whether a D/E/F line has a public same-line perturbation bank covering its panel (our registry: HCT116, K562, DLD-1,
Jurkat, HepG2, A549, HT29, MCF7, BXPC3, HAP1, Calu-3, THP-1).

Tasks:
1. Evidence for the line choice: Arc Institute's 2026 challenge materials (the Cell commentary doi
   10.1016/j.cell.2026.08.004 and any supplement, Arc news/blog posts, talks, preprints, the 2025 challenge wrap-up), Arc's
   other published perturbation screens and cell-line collections, and what the A/B/C choice suggests (tissue coverage,
   Flex compatibility, CRISPRi-competent lines Arc already had, lines shared with Tahoe-100M or other Arc resources).
   Rank the most likely D/E/F candidates with your reasoning and confidence. Do NOT contact anyone.
2. For each top-15 candidate line: every public single-cell perturbation dataset (accession, lab, modality, chemistry,
   number of targeted genes, cells per target, licence, size, download URL), and an estimate of how many of a
   Arc-style non-essential 300-gene panel it would cover (use the A/B/C panel in data/challenge/pert_counts.csv as the
   proxy for the panel's gene class — it is tracked in git).
3. Anything released since 2026-09-20 (GEO, ArrayExpress, Zenodo, figshare, Hugging Face, bioRxiv/medRxiv supplements):
   new genome-scale or large CRISPRi/KO single-cell screens in any cancer line.

Deliverable: experiments/520-cloud-def-lines/RESULTS.md — the ranked candidate list with evidence, the per-line dataset
table, the new-release list, and a download manifest (URL, size, sha256 if available) for anything that would add a
same-line bank we lack. Verify every dataset claim by opening its landing page or file listing (write the exact sentence
or listing you relied on). Mark unverified items as such.
