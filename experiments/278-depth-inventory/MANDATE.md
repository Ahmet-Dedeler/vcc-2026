# 278 — The DEPTH axis, inventoried: for every panel target, the deepest measurement of its response available per cell line on disk or in the public record, and what a ≥ 400-cell somatic measurement would cost. Metadata and public-catalogue work only. ≤ 3 hours. No download > 2 GB, no build, no submission.

You are a research subagent. Work in `experiments/278-depth-inventory/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/` or `data/` (anything fetched goes to `work/`, ≤ 2 GB total, metadata
tables only). `uv run python`. Peak RSS under ~4 GB. You may use web search / fetch for public catalogues. **Time box:
3 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you.

## Why (R276 §6.3, F1, F8 in `docs/FINDINGS-2026-09.md`)
The programme has now measured two of its three scaling axes on this panel: more CONTEXTS at 150 cells add little
(276: 17 extra lines beat the two production donors on one of three readings), and model CAPACITY / knobs around the
incumbent are exhausted (READ-FIRST §5). The third axis — DEPTH, cells per condition — is where the reliability law
(F1: split-half reliability 0.05–0.09 at 150–200 cells, a third of Orion targets without a measurable direction; Flex
at 400 cells resolves every target) says the information is. The reliability-corrected conservation map (276: median
+0.228 vs the raw map) says the ceiling is far above what 150-cell responses show. **Before anyone proposes acquiring
or generating deeper data, inventory what depth exists: per panel target and line, on disk and in the public record.**

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; `docs/FINDINGS-2026-09.md` F1, F8; 267's `RESULTS.md` §1 (the on-disk per-target
cell counts for K562-GWPS, HCT116, HEK293T); 270's `RESULTS.md` §4 (the on-disk datasets, contexts, units); 276's
`RESULTS.md` §1 (the inventory of (dataset, line) per target — REUSE it; add depth); `scripts/scperturb_files.tsv`
(the scPerturb catalogue on disk and its sizes); `data/external/*/` metadata files (do not load count matrices; read
obs tables / metadata only); AGENTS.md §2.1 (the panel's design).

## What to establish
1. **On disk, per panel target × line: cells per condition** (perturbed cells; and controls), guides per target,
   batches — from the obs tables of every genetic-perturbation set on disk (the three production donors, Replogle
   K562/RPE1, Nadig HepG2/Jurkat, Jiang's 6 lines, Feng iPSC, HSPC, KOLF2.1J, DLD-1, IGVF sets, Myllia A549, CD4
   conditions, Broad-16's 90 non-panel targets for context). Table: target × line → cells; summary: for each target,
   the deepest somatic measurement on disk; how many panel targets have ≥ 400 cells in ANY somatic line, ≥ 200, ≥ 100.
2. **The public record beyond disk** (catalogues only: scPerturb's full list, GEO/SRA searches for the panel's gene
   symbols with "Perturb-seq / CRISPRi / CRISPR screen single-cell", the 2025/2026 VCC data releases, X-Atlas/Orion's
   published depth, Replogle 2022's full GWPS depth vs what we extracted, Tahoe (chemical, note only)): which panel
   targets have a ≥ 400-cell single-cell perturbation measurement in a somatic line anywhere, in which study, at what
   chemistry, and whether the data are downloadable and under what terms. Do not download expression matrices.
3. **The chemistry / depth regime of the challenge's own references**: what is public about the 2026 A/B/C references
   (Flex, 400 cells, ~20 k UMI — from `AGENTS.md` §2 and `docs/reference/`), and which public sets share it (h1-test does).
4. **Cost of depth**: for the deepest gap (panel targets with no somatic measurement ≥ 200 cells), what it would take to
   close it — (a) by download (which datasets, GB, terms), (b) by re-extraction from files already on disk that we
   sub-sampled (did `build_donor_table.py` cap cells per target? R267 §1 says the full counts were used — verify), (c) by
   pooling guides / replicates we currently split. Report each as a table with numbers.

## Deliverable
`RESULTS.md`: the on-disk depth table (per target × line, summarised) with the ≥ 400 / ≥ 200 / ≥ 100 counts; the public-
record table; the reference-regime statement; the cost table; and one line: **how many of the 300 panel targets can
reach a ≥ 400-cell somatic measurement with data that exists today, and from where**. Save the per-target table as
`work/depth278.csv`. No model, no submission. Then stop.
