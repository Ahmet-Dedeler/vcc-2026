# 413 — source scan: new public perturbation data for cancer lines (mandate from the PI, 2026-09-22 19:35Z)

**Current state (read these first; they supersede docs/READ-FIRST.md):** `experiments/410-pertscale-revisit/BOARD-RESULT-410b.md`
(new best 0.2366), `experiments/410-pertscale-revisit/RESULTS.md`, `docs/NEXT-SLOT-DECISION.md` (last two sections).
Our predictions average donor perturbation responses from three production lines (K562-GWPS, X-Atlas/Orion HCT116 and
HEK293T) plus CD4 T cells and H1; direction quality follows a 1/k law in the number of independent lines (AGENTS.md §8.1,
exp 159). More panel-covering lines are the most valuable thing that could exist.

**Your territory (others: 414 re-prices bulk arms on PLOLO, 415 builds a key — do not do their work):** find public data
released or updated since 2026-08-15 — and anything the earlier sweeps missed — with CRISPR(i/a/ko) perturbation
transcriptomics (Perturb-seq, CROP-seq, bulk RNA-seq screens, L1000 XPR) in HUMAN CANCER or transformed cell lines that
(a) cover ≥ 30 of our 300 panel targets (`data/challenge/pert_counts.csv`, column target_gene), or (b) are in or near
HeLa (cervical adenocarcinoma), CAL33 (head-and-neck squamous), Jurkat (T-ALL) — our contexts A/B/C (exp 041) — or other
epithelial carcinoma lines likely for the Oct-22 panel. Check X-Atlas/Pisces status, Arc Virtual Cell Atlas releases,
scPerturb updates, GEO / ArrayExpress / Zenodo / figshare / Hugging Face (`hf datasets ls --search ... --limit 100`),
CZI, Illumina's Billion Cells Project (any public release?), and 2026 bioRxiv/arXiv preprints with data links.

**Read first so you do not re-report:** experiments/045-data-scout, 166-donor-acquisition, 225-source-sweep,
270-data-reconnaissance RESULTS.md; `ls data/external/` (what we already hold).

**Known answer (do first):** reproduce panel overlap for sources we hold from their perturbation lists:
`data/donors/{k562_gwps,orion_hct116,orion_hek293t,nadig_jurkat}/perts.txt` ∩ panel = 269 / 290 / 294 / 0.

**Deliverable:** `experiments/413-source-scan/RESULTS.md` — a table (dataset, URL, release date, cell lines, modality,
#perturbations, EXACT panel overlap when a gene/perturbation list ≤ 200 MB is downloadable, total size, license/access)
and a ranked recommendation of what the PI should download and build into a donor bank. Do not download anything
> 2 GB (list it for the PI instead). No builder, no vcc, no cell-eval2; write only in experiments/413-source-scan/.
