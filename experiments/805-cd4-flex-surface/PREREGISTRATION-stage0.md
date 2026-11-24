# 805 stage 0 — the power gate on GWCD4i (preregistered before any cell-level expression value is read)

Written 2026-10-16 ~02:00Z by the 805 worker. Mandate: `MANDATE.md` (lab lead [c61d4c]). Only metadata has been read so far:
the published supplementary tables on disk (`data/external/cd4_cell2026/suppl_tables_2026-08-24/`) and the S3 listing.
A layout probe (`probe805.py`) reads HDF5 structure, obs, var and indptr only, never expression values.

## What is already known from the published tables (DESeq2 pseudobulk, all ~550 cells, 10 % FDR; NOT our test)
Rest, `ontarget_significant` and `n_cells_target` ≥ 400 and not K562-essential: 4,543 targets. Their published
`n_downstream`: q10/q25/q50/q75/q90 = 0 / 0 / 2 / 6 / 52; 14 % have ≥ 20. A/B/C's own panel targets in Rest (283; 0
essential) have median 4 (163 eligible by the same filters: median 3, 21.5 % ≥ 20). **So the expected outcome of the
mandate's random draw is a FAIL for target-selection reasons alone**: a random genome-wide CRISPRi target in resting
T cells mostly does nothing downstream. H1's 89 is on Arc's strength-stratified panel. That is why a secondary readout
on Arc-selected targets is fixed here too (below); it does not change the gate verdict.

## Cells (both target sets, one rule)
- Source: `s3://genome-scale-tcell-perturb-seq/marson2025_data/D{1,2,3,4}_Rest.assigned_guide.h5ad`, read on EC2.
- Usable cell: `low_quality` False and exactly one assigned guide (`guide_id` not "multi-guide").
- Controls: usable cells whose `guide_type` is the non-targeting class.
- A target's cells: usable cells whose guide maps (sgRNA library table `target_gene_id`, the curated target; fallback
  `perturbed_gene_id` of obs when a guide is absent from the table) to the target's Ensembl ID (`target_contrast` of
  `DE_stats.suppl_table.csv`), AND whose guide has `signif_knockdown` True in Rest (`guide_kd_efficiency.suppl_table.csv`;
  Arc's panel keeps clean, strongly knocked-down constructs).
- Donor mixture: π_d = donor d's share of all usable control cells (the four files). Controls: 18,400 = Σ_d n_d with
  n_d = largest-remainder(π_d × 18,400), drawn without replacement per donor. Each target: 400 cells with the same
  quotas (largest-remainder π_d × 400); a donor short of its quota gives its shortfall to the other donors in
  proportion to their spare cells (recorded per target). Seed 805 (numpy default_rng([805, k]) per group k).
  Rationale: four donors (possibly of different sex) must not enter as a target-vs-control composition difference.
- Gene axis: the official 18,533 genes (`data/challenge/gene_names.csv`) matched by symbol to the file's
  `var.gene_name`; genes absent from the file are zero columns; a symbol present twice in the file takes the first
  column (count recorded). Overlap recorded.
- Depth: per-cell totals T_c on the official axis over all selected cells (targets + controls of BOTH sets);
  f = min(1, 20,000 / median T_c); every count is binomially thinned with probability f (seed 805), Arc's "median depth
  of 20,000 UMI per cell" by a single global factor. f and the post-thinning median are recorded. If f = 1 the truth is
  shallower than Arc's and that is stated next to the gate.

## Target sets
- **P (primary, decides the gate; the mandate's rule):** pool = Rest rows of `DE_stats.suppl_table.csv` with
  `ontarget_significant` True, `n_cells_target` ≥ 400, symbol on the official axis, symbol not in
  `data/donors/k562_essential/perts.txt`. Order = `np.random.default_rng(805).permutation(sorted(pool))`; walk it and
  take the first 100 targets with ≥ 400 usable cells under the rule above.
- **S (secondary, reported, does NOT decide):** every target of `data/challenge/pert_counts.csv` (A/B/C's panel,
  Arc-selected) that passes the same pool filters and has ≥ 400 usable cells. Purpose: separates "the line has no DE
  power at Arc depth" from "random targets do nothing". S's targets are A/B/C's, so S can never be part of 805's surface
  (stage 1 excludes them); reading S's truth chooses nothing that is later priced on the surface.

## Statistics
- DE: `experiments/614-k562-real-emitter/ref614.py` on the truth (evaluate323's reproduction of cell-eval2's DE:
  Wilcoxon vs the 18,400 controls, BH, the CPM gate on controls). P and S are each built as their own truth with the same
  18,400 controls (two `ref614.py` runs, so BH families are per target as in the scorer).
- n_real per target = q < 0.05 on the DE universe, the target's own gene excluded (ref614's count).
- Knockdown: the target's own gene, mean CPM in its cells / mean CPM in controls ("remaining"), and its lfc and q from
  the reference (when the gene is in the DE universe).
- Blob test: 615's split-half share (`blob615.old_surface`, TS 5e4, target genes excluded, per-label random halves),
  on the P truth with its single control pool.

## Known-answer rows (run on the same instance, same code, BEFORE reading P or S; abort the screen if either fails)
- **KA-1 (DE engine):** `ref614.py` on `data/harness/h1-test.h5ad` with its 85 targets must reproduce 614's
  `work/ref_h1ka/reference.json` `de_sha256` = 1bb8c76d1177047569e156af5a7c40e48f8fedce0ce3cf85b7d82998c5485e99 exactly
  (and so its n_real median 89).
- **KA-2 (blob statistic):** `blob615.old_surface` on h1-test: split-half share 0.0487 and plain 0.0376 (615 RESULTS),
  |Δ| ≤ 0.002 (the same seed and code: expected exact).

## The gate (mandate, verbatim numbers), on P
**PASS iff median n_real ≥ 20 AND ≥ 60 % of the 100 targets have n_real ≥ 10.** (H1: 89 and 68 %; K562: 5; Orion clean: 1.)
Also reported on P: the knockdown distribution and the blob split-half share (must be ≤ 0.10; a share > 0.10 is a
FAIL of the surface whatever the power).

## If P fails on Rest: Stim48hr, once, by the same rule
Identical rule with `culture_condition == "Stim48hr"` and the `D{1..4}_Stim48hr` files (guide efficiency from Stim48hr,
seed 805, P and S re-drawn by the same procedure, KA rows re-checked on the same instance). Then stop and report
whatever it reads. No third condition, no re-draw, no threshold change.

## Reported with no rule attached
n_real quantiles (0/10/25/50/75/90/100) for P and S; S's gate statistics; DE DOWN share; per-target n_real vs the
published `n_downstream` (Spearman); the axis overlap; f; library medians of targets and controls; donor quotas and
shortfalls; the 400-cell availability (how many pool targets had to be skipped).

## Amendment 1 (2026-10-16 02:25Z, after KA-1 failed on the first launch; no CD4 expression value had been read)
Run `20261016T021043Z-s0-805` stopped at KA-1 as this file required (logs: `work/s0_ka_fail1/`). On the x86 EC2 instance
ref614.py on h1-test gave **n_real per target identical to 614's ref_h1ka for all 85 targets** (median 89, 58 ≥ 10) and
KA-2 passed (split-half 0.048693, plain 0.037613), but the DE parquet's sha256 was 796a7e29… instead of 1bb8c76d…, and the
npz sha differed too. A byte-level sha across two platforms (Apple arm64, where 1bb8c76d… was made, vs x86_64) is the
wrong known answer: libm and BLAS last-bit differences change the float bytes without changing a single call.
**KA-1 is replaced by:** (a) n_real_by_target identical for all 85 targets to 614's `ref_h1ka/reference.json`, AND
(b) the DE frame (`reference_de.parquet`) has exactly the rows (target, feature) of 323's
`experiments/323-current-family-validation/work/reference_de.parquet` (the frame 1bb8c76d… hashes), and max |Δ| over
p_value, p_adj and log2_fold_change ≤ 1e-9, AND (c) the BH call set (p_adj < 0.05) is identical. Nothing else changes.

## Amendment 2 (2026-10-16 02:35Z, after the second launch stopped at KA-1; still no CD4 expression value read)
Run `20261016T021703Z-s0b-805` (logs `work/s0_ka_fail2*`): rows identical, **n_real identical for all 85 targets, the
BH call set identical (916,215 tests)**, log2FC max |Δ| 2.2e-16; but p-values differ in 521 of 916,215 tests by 1e-9 …
1.18e-4 RELATIVE (max abs 7.4e-6; 287 of them in COASY, then SETDB1 48, MED19 38), the rest ≤ 7.5e-12 relative. The
pattern (a few targets, every p range) is the Wilcoxon tie structure moving with last-bit differences of normalised
values between arm64 and x86_64; it never moves a call at the scorer's α. The 1e-9 bound of Amendment 1 was set without
measuring this. **KA-1 is now:** same rows; n_real_by_target identical; the call set (p_adj < 0.05) identical;
max |Δ log2FC| ≤ 1e-9; p_value and p_adj relative |Δ| ≤ 1e-3 for every test and ≤ 1e-9 for ≥ 99.9 % of tests. The
gate's statistic (n_real, a count of calls) is exactly reproduced, which is what this screen depends on.

## Amendment 3 (2026-10-16 02:40Z; launch 3 stopped in `select` after reading obs/var only)
The four donor files do not share one var axis (D3's differs from D1's). Each file is now mapped to the official axis by
its own symbols, and only genes present in ALL four files keep counts (a gene some donor lacks is zero for every cell),
so a donor's missing gene can never become a target-vs-control difference. Per-donor overlap, the common and union counts
are recorded in select.json. Nothing else changes.

## Amendment 4 (2026-10-16 02:52Z; written after launch 4's `select.json` (obs counts only) and BEFORE any gate number)
Launch 4 (`20261016T023623Z-s0d-805`) selected P and S under the shortfall rule: ~170 of its ~260 targets could not fill
their donor quotas (shortfalls up to 244 of 400 cells: many targets are carried by only some of the four donors). With
one pooled control set, such a target's donor mix differs from the controls', so donor-specific genes (sex, HLA, eQTL
genes) would be called DE and inflate n_real. **The rule is now: a target qualifies only if every donor holds at least
its quota (n_d ≥ quota_d for all four donors); no shortfall transfers.** P walks the same seed-805 permutation further;
S keeps the A/B/C targets that qualify. Everything else unchanged. Launch 4 is allowed to finish as a DIAGNOSTIC only
(its gate numbers never decide anything): comparing its shortfall vs no-shortfall targets measures the donor artefact.
