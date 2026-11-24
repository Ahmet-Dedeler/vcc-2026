# 293 — The construct ceiling in the challenge's OWN regime: on the 2025 H1 dataset (Arc's dual-guide CRISPRi library, H1 hESC, 10x Flex, 400 cells per construct), how well do two constructs of the same gene agree, against two halves of one construct's cells, at matched cells — on every target that carries two constructs (train, validation, test splits)? Preregistered, ≤ 3 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/293-challenge-construct-ceiling/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/`; no download. `uv run
python`, streaming / row-blocked reads only (the 2025 files are 11–15 GB), peak RSS under ~10 GB, one file at a time.
**Time box: 3 hours.** One sibling (292) runs in parallel; do not read or write its directory. Both reserves are
SPENT (`experiments/267-*/work/reserve/`; Broad-16 JVE127 / MG63 / SLR23 / UMRC3) and are not needed here.

## Why (F14, `docs/FINDINGS-2026-09.md` v1.9; R285 §6, R273 §2)
285 measured, on one matched design over five public CRISPRi datasets, that the guide construct is the dominant
irreducible term of a perturbation response: two constructs of one gene agree at ρ_guide 0.09 (K562-GWPS, 748
perturbations), 0.17 (Jurkat), 0.41 (HepG2), 0.02–0.06 (CD4) of the same-construct cell split, while batches are free
and a different study on the same line keeps ρ 0.79. The 2026 challenge's truths are 400 cells of ONE construct per
target in six lines, zero-shot, and the 2026 panel shares only 25 of 300 targets with the 2025 H1 library — so no
source on disk carries the 2026 constructs and the construct ceiling binds every entry. **What is that ceiling in the
challenge's own regime — Arc's dual-guide design (two guides per vector, validated for knockdown), H1, 10x Flex, 400
cells?** 273 read it on the 5 two-construct targets of the test split (P1 vs P2: median +0.014 against +0.243 for lane
halves, R273 §2) — n = 5. The training and validation splits carry more (8 in training: TMSB4X, KLF10, KDM2B, PRCP, …).
This experiment measures all of them on one matched design and states the ceiling with its interval.

## Inputs (metadata first, then the cells of the two-construct targets only)
`data/external/vcc2025/train/adata_Training.h5ad` (14.8 GB; obs `batch`, `guide_id`, `target_gene`; 151 targets, 189
guide categories, **8 targets with two constructs**), `validation/adata_Validation.h5ad` (count its two-construct
targets), `test/adata_Test.h5ad` (11.4 GB; 101 targets, 136 categories, **5**: SOX4, GRB10, MED15, PCCB, DPF2) and its
2026-faithful subsample `data/harness/h1-test.h5ad` (52,400 cells; 273's object — the 5 test targets' construct halves
are in R273 §2 and `experiments/273-flex-reference-instrument/truth273.py`; reuse its residual-space and lane-split
code). Construct naming: `GENE_P1_A|GENE_P1_B` vs `GENE_P2_A|GENE_P2_B` (two promoters) or two transcript ids;
the common form `GENE_P1P2_A|GENE_P1P2_B` is ONE construct. Non-targeting controls carry their own guide ids.
Also `experiments/285-reliability-benchmark/lib285.py` (the matched-n nested ladder, thinning, the two readouts:
residual cosine and 218's scorer statistic; the ρ construction) — reuse verbatim where it fits.

## Protocol (write `PREREGISTRATION.md`, hashed into `work/PREREGISTRATION.sha256`, before any cosine)
1. **Inventory (metadata only):** every target with ≥ 2 constructs in each split, cells per construct, batches per
   construct; the controls per split. State n (expected 13 + validation's).
2. **Known answers:** 273's lane-parity 200/200 median on h1-test (+0.242920) at 0.0e+00 with your code; 273's
   5-target construct row (median +0.0135; per target DPF2 −0.03, MED15 −0.10, GRB10 +0.01, PCCB +0.03, SOX4 +0.35)
   reproduced from `h1-test.h5ad`; the union of two halves reproduces the full row (integer identity).
3. **The matched design (frozen; 285's construction):** per two-construct target, at n ∈ {25, 50, 100, 200} cells per
   side where both constructs allow (8 draws averaged): (a) same construct, random cell halves — the cell-sampling
   reference; (b) construct P1 vs construct P2; (d) lane halves within one construct; residual space = each side's
   profile minus the mean over that split's targets on the same side (159 / 273 convention; build the centring from a
   matched-n sample of every target of the split, never from full-depth rows — R285 §3's artefact), target-gene columns
   removed; the split's own non-targeting cells as the control for the scorer's statistic. **ρ_guide** = (b) / (a) at
   the same n, per target and pooled (median), paired target bootstrap (20,000 draws); Spearman–Brown extrapolation of
   (a) to 400/400 and the implied between-construct agreement at the challenge's budget; the scorer's per-target `pds`
   (218's `_per_target`, the split's other targets as distractors) for (a) vs (b). Also the on-target knockdown of each
   construct (its own gene's column, excluded from every scored statistic) and whether the weaker construct explains
   the disagreement (285 §6b's stratification, here as a per-target covariate).
4. **The challenge-regime ceiling, stated:** the pooled ρ_guide with its interval; the SB-extrapolated between-
   construct cosine at 400 cells; the fraction of a 400-cell truth's reproducible variance that is construct-specific
   (1 − ρ_guide); how it compares with 285's public-library values (0.09 / 0.17 / 0.41) and with 273's 5-target read.
5. **Descriptive (post-hoc, labelled):** which targets' constructs agree (SOX4-like) and which do not (MED15-like), and
   whether knockdown strength or promoter choice (P1 vs P2 vs transcript) separates them.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the inventory; the KA rows; the per-target table (cells per construct,
(a), (b), (d) at each n, ρ_guide, on-target knockdown per construct, `pds` (a) vs (b)); the pooled ceiling with
intervals and the SB extrapolation; the comparison with 285 / 273; `work/construct293.csv`; one line: **CHALLENGE-
LIBRARY CONSTRUCT CEILING (H1, Flex): ρ_guide … [lo, hi] on n targets; between-construct cosine at 400/400 (SB) … vs
same-construct …; scorer `pds` (a) … → (b) …; knockdown explains … of it**. No model, no build, no submission. Then stop.
