# 270 — Declared SECONDARY feature: cross-LINE sign stability of a target's knockdown signature in LINCS/CMap 2020

Written 2026-09-05 ~16:45Z, BEFORE any CMap signature values were read for the panel (only the curated signature
metadata `data/external/cmap/genetic_sigs.parquet` — sig_id, line, pert_type, gene — and the gctx shapes were read).
This is an ADDITION to the primary GPSAdb feature (`FEATURE-DEFINITION.md`, sha f8a6bd39…): same statistic, a
different, denser resource (many lines, one platform). Labelled secondary; it does not change the primary claim.

## Inputs (all on disk, nothing downloaded)
- `data/external/cmap/level5_beta_trt_sh.gctx` (MODZ level-5 z-scores; rows = signatures, columns = 12,328 genes) —
  only the `trt_sh.cgs` signatures (consensus across shRNAs, one per gene per line).
- `data/external/cmap/level5_beta_trt_xpr.gctx` — `trt_xpr` (CRISPR knock-out) signatures; several per (gene, line);
  the per-(gene, line) signature is the plain mean of their z-vectors (stated here; MODZ not recomputed).
- `data/external/cmap/geneinfo_beta.txt` for the landmark flag; the 300-panel symbols as in `pert_counts.csv`
  (`cmap_name` exact match).

## Definitions (fixed before computing)
For target t and a set of per-line signature vectors z_l (l = cell line):
1. `S_l` = top **N = 100** genes by |z_l| among the **978 landmark genes** (directly measured), excluding t itself.
2. For every pair of distinct lines (a, b): `U_ab = S_a ∪ S_b`; `c_ab` = fraction of U_ab with sign(z_a) = sign(z_b)
   (zero counts as disagreement); pairs used only if |U_ab| ≥ 20 (always true here except degenerate signatures).
3. `SS_cmap_sh(t)` = mean of c_ab over line pairs of the shRNA consensus signatures (primary variant of this file).
   `SS_cmap_xpr(t)` = the same over line pairs of the CRISPR-KO per-line means.
   `SS_cmap_method(t)` = the same statistic between the shRNA consensus and the CRISPR-KO mean **within the same
   line**, averaged over lines that have both (cross-method, same context).
   Sensitivity rows (declared): all 12,328 genes instead of landmarks only; N = 50.
Chance is 0.5; recorded per target: number of lines (sh, xpr), number of pairs.

## Tests (identical to the GPSAdb file; run by `test270.py` on `work/feature270_cmap.csv`)
Spearman with 262's mean cancer-fold `cd4svd5` `pds` (paired target bootstrap B = 2,000, permutation p), AUC for
rest vs anti-aligned (majority rule over cancer folds), ρ with kolf2 apart, ρ with the pipeline-independent
held-vs-donor truth agreement. Decision rule as in the primary file (both intervals must exclude chance).
Controls: the number of lines per target (a popularity proxy) is reported against `pds` as well; if the feature
tracks `pds` no better than the line count, it is not a stability signal.
