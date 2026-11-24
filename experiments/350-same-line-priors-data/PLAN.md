# 350 — Same-line priors data inventory (PLAN)

## Question
Contexts A/B/C are Jurkat / HeLa / CAL33 (041). We hold no perturbation training data
from those lines, but public same-line perturbation transcriptomics may exist — most
concretely Nadig 2025 Jurkat Perturb-seq. 341 showed a line's OWN non-panel targets can
fit a calibration map (+0.033 weighted sign, +0.07 cosine; rejected under a retrieval
guard). On the board, nmae/reach/fid are half the score and 343 shows the deficit is
per-gene direction information at the design stage. This experiment is DESCRIPTIVE ONLY:
inventory what same-line data exists per context, and measure what it could supply
(common response vector, per-gene DE propensity, calibration-map target counts).
No arm, no gate, no score claim.

## Inputs (paths)
- `data/external/nadig2025/jurkat_raw_singlecell.h5ad` (9.37 GB) and
  `data/external/nadig2025/hepg2_raw_singlecell.h5ad` (5.61 GB) — dense float32 X,
  obs label col `gene` (categorical, 2,394 levels incl `non-targeting`), var symbols in
  `var/__categories/gene_name` (bytes, decode). Code path reused:
  `scripts/build_donor_table.py` DONORS dict (`label_col='gene'`,
  `gene_col='gene_name'`, `control='non-targeting'`); prebuilt donor tables at
  `data/donors/nadig_jurkat/` and `experiments/166-donor-acquisition/donors/nadig_hepg2/`.
- `data/challenge/pert_counts.csv` (300 panel targets), `data/challenge/gene_names.csv`
  (18,533 axis), `data/challenge/context_A.h5ad` (18,400 controls, CSR).
- H1 design: `experiments/323-current-family-validation/work/effects.npz['cap15']`
  (85×18,533, 319 composition design effects), `per-target-cap15.csv`,
  `metrics-cap15.json` (known-answer source).
- Donor comparison: `data/donors/k562_gwps/delta.npy` (+`genes.txt`) if present.
- cmap: `data/external/cmap/{cellinfo_beta.txt,siginfo_beta.txt,geneinfo_beta.txt,genetic_sigs.parquet,level5_beta_trt_sh.gctx}`.
- CCLE/DepMap: `data/external/ccle_broad_2025/ccle_cache.npz` (+`data_mrna_seq_tpm.txt`),
  `data/external/depmap26q1/`.

## Estimator / procedure
1. **Inventory (h5py only, no X load):** n_cells/n_genes from `X.shape`; obs/var column
   lists; perturbation label cardinality (`obs/__categories/gene`, incl control);
   control cell counts from `obs/gene` codes; gene id type (`var/gene_id` Ensembl +
   symbol categories; note HSPA14-class duplicates if seen); overlaps: panel∩bank
   targets, bank-var∩challenge-axis, bank-targets∩challenge-axis (bytes-decoded).
2. **Common vectors (chunked row-block streaming, float64 accumulators):**
   For each bank: codes for control + seeded-random 300 non-panel targets
   (`rng = np.random.default_rng(350)`); stream X in 8,192-row blocks; accumulate group
   sums (301×G). CPM = 1e6·sums/row-total; lfc = log2((CPM_pert+1e-9)/(CPM_ctrl+1e-9))
   (exact `exact_de323.de_from_counts` lfc formula). Restrict to challenge-axis genes
   via symbol→index map. Common vector = mean over 300 targets. Report: entry-level
   fraction lfc<0; per-target fraction with mean lfc<0; fraction of common-vector
   entries <0; top 25 up / top 25 down genes by signed value. DE propensity per gene =
   fraction of 300 targets with |lfc|>0.5; Spearman(propensity, context-A mean CPM)
   on shared genes (context A: backed CSR read, chunked mean of CPM).
   Peak target ≤ ~1 GB (one dense block 8192×9624×4 B ≈ 315 MB).
3. **Comparisons (cosine + sign agreement):** cos(Jurkat common, H1 cap15 row-mean
   restricted to shared genes); cos(Jurkat common, HepG2 common); cos(H1 design,
   HepG2 common). Donor-bank mean only if a `delta.npy` mean is directly loadable
   (state the space mismatch: donor deltas are log1p(TS·p) natural-log pseudobulk
   deltas vs our log2 CPM ratios — direction comparison with caveat, not a score).
4. **Web + local checks:** HELA in cmap `genetic_sigs.parquet` (count, pert_type
   split trt_sh vs oe, gene axis from `geneinfo_beta.txt`/gctx); HeLa Perturb-seq /
   CROP-seq / ECCITE / shRNA compendia literature; CAL33 / CAL-33 / HNSCC
   (SCC-25, FaDu, Detroit 562, Cal27) knockdown/knockout expression literature;
   CCLE/DepMap basal rows for CAL33 (verify sample ids in `ccle_cache.npz`).
   Cite URLs + access dates (2026-09-10).
5. **Outputs:** `work/jurkat_common.npz` + `work/hepg2_common.npz` (gene names on the
   challenge-axis intersection + vector in the log2-CPM-ratio space + target list +
   meta), `work/inventory.json`, RESULTS.md tables, every load command.

## Known-answer reproduction (before new measurement)
Own code recomputes from 323's saved files: mean of `per-target-cap15.csv` rows with
`metric == 'de_wilcoxon_lfc_nmae'` must equal `metrics-cap15.json` raw nmae
0.9508364630325662 (343's value), and `effects.npz['cap15']` shape must be
(85, 18533). If either fails, stop.

## Controls / falsifiers
- HepG2 is the control line (NOT A/B/C): if cos(Jurkat common, HepG2 common) ≈
  cos(Jurkat common, H1 design), the "common vector" is a generic CRISPRi axis, not
  line-specific information — kills any "Jurkat common vector fixes context A"
  proposal at the design stage.
- Sign agreement is reported alongside cosine; a high cosine driven by a few large
  genes with chance signs elsewhere is stated as such (top-50 inspection).
- H1-panel numbers are development evidence, not A/B/C forecasts (rule 5).

## Kill criterion
This is an inventory: it succeeds as a clean descriptive table even if every cosine
is ~0. The idea it would kill: "same-line common vector / propensity is
line-specific enough to supply design-stage direction" — dead if Jurkat–HepG2 common
cosine ≥ Jurkat–H1-design cosine, or if both cross-line cosines are ≈0.

## Budget
~3 h: mostly single-pass streaming (≈15 GB reads, minutes) + web. Backed/h5py only;
no densification; scratch deleted after.
