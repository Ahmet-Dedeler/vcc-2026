# 301 — Three metadata acquisitions the write-up is waiting on, and the one analysis they unlock: (1) the VCC 2025 sgRNA library table (protospacer sequences), (2) promoter sequences for a seed scan, (3) the JoungZhang 2023 TF-ORF id → symbol table; then T7's seed scan of the 2025 library's 13 failing constructs if (1) and (2) land. ≤ 2.5 hours. No expression value of any reserve is read. No build, no submission.

You are a research subagent. Work in `experiments/301-metadata-acquisitions/`. **Edit nothing outside it** except
appending, if item 3 succeeds, a "Reserve R4 candidate (metadata only)" note to `docs/CONFIRMATION-LEDGER.md` in R3's
format. Downloads go to `work/`, **≤ 200 MB total**, each with URL, licence / terms and sha256 recorded in
`work/sources301.md`. `uv run python`; peak RSS under ~6 GB. **Time box: 2.5 hours.** No `vcc`, no `cell-eval2`,
nothing under `data/` or `submissions/`. **Reserve rule:** `JoungZhang2023_atlas.h5ad` may be read for `obs` / `var` /
`uns` only — never `X`; R3's eight files are SPENT and are not touched.

## Why
`docs/FINDINGS-2026-09.md` v2.5 and the manuscript (`experiments/294-writeup/PAPER.md`) carry three items blocked on
metadata that is not on disk (R298 §4; R292 §5.3): (T7) the 2025 challenge files hold protospacer NAMES
(`GENE_P1_A|GENE_P1_B`), not sequences, and no genome or promoter sequence exists on this machine — so whether a
failing construct's reproducible response is a seed-match off-target cannot be tested; (R4) the JoungZhang 2023 hESC
TF-ORF atlas (1.1 M cells, 3,369 ORFs) is the largest untouched human genetic-perturbation set on disk but its
`obs.perturbation` is a bare numeric ORF id no on-disk column resolves. Both need ≤ 100 KB–50 MB of public metadata.

## Items
1. **The 2025 VCC sgRNA library table.** Look where the 2025 release documents its guides: the Arc Virtual Cell Atlas
   (`github.com/ArcInstitute/arc-virtual-cell-atlas`, `virtual-cell-challenge/`), the Arc "Behind the Data" page, the
   2025 Cell paper's supplementary tables (S0092-8674(25)00675-0; the publisher blocks plain fetches — try the
   supplementary-file URLs, PMC, or the Arc data bucket), GEO / the challenge's S3 bucket. The target is a table
   mapping the 2025 `guide_id` names (e.g. `SOX4_P1_A`) to protospacer sequences (≈ 410 constructs × 2 spacers). If
   found: sha256, terms, and a check that its names match the three splits' `obs/guide_id` categories (`data/
   external/vcc2025/{train,validation,test}/*.h5ad`, obs only). If not found after a genuine search, say exactly
   where you looked.
2. **Promoter sequences for the seed scan** (only if item 1 succeeded): TSS ± 500 bp for protein-coding genes from
   Ensembl BioMart (the `upstream_flank` / `downstream_flank` sequence query, GRCh38, current release) or a
   comparable public source, ≤ 50 MB gzipped; strand-aware; keep gene symbol + Ensembl id.
3. **The JoungZhang 2023 ORF id → symbol table.** The paper's supplementary table (Joung et al., Cell 2023, "A
   transcription factor atlas of directed differentiation") or the GEO series metadata (GSE217460 / GSE217812 — verify)
   or the scPerturb curation notes; ≤ 100 KB. If found: map every `obs.perturbation` id in the atlas's `obs` (read
   `obs` only) to a symbol, report coverage, and append the R4 CANDIDATE note to the ledger — files, sha256, cell
   counts, modality (ORF overexpression; sign not comparable to knockdown), the same once-per-claim rule as R3,
   and what it cannot confirm. Do not declare it a reserve for any claim; the coordinator does.
4. **T7's seed scan** (only if items 1 and 2 landed): for the 13 failing arms of the 2025 library's 17 two-construct
   targets (`experiments/298-reviewer-checks/work/seed298.csv`; plan `work/seed298_plan.md`), scan each failing
   construct's two spacers (PAM-proximal 9–12 nt seed, NGG PAM, either strand) against every promoter window;
   for each candidate off-target gene, read the failing construct's residual response at that gene's column from
   293's saved rows (`experiments/293-challenge-construct-ceiling/work/` — state which array; if the per-construct
   residual rows are not saved, stop at the candidate list) and report whether the candidate is knocked down
   (lfc below the library's own-gene |lfc| median). Null: the 21 non-failing arms scanned the same way. Report
   candidates per arm and the failing-vs-working contrast; n = 13 vs 21, so an existence proof at best.

## Deliverable
`RESULTS.md`: what was found and where (items 1–3, with sha256 and terms), the guide-name match check, the ORF
coverage, the R4 note (or why not), the seed scan (or the exact blocker); `work/sources301.md`; one line: **SPACERS:
FOUND / NOT FOUND (source); PROMOTERS: FETCHED / SKIPPED; JOUNG IDS: n / 3,369 resolved; SEED SCAN: n / 13 failing
arms with a candidate off-target (working arms n / 21) — or NOT RUN (blocker)**. Then stop.
