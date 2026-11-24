# 404 — A better neighbourhood for the block-sign tilt: similarity metric, source of neighbours, K

## Hypothesis
The shipped kNN tilt (399: 25 nearest K562-GWPS perturbations by full-row Pearson, median
top-1 r ≈ 0.11) leaves signal on the table: full-row Pearson is diluted by ~8k mostly-quiet
genes plus a shared common-response direction. A denoised similarity (query-focused gene
subset, common-response removal, rank-based, or PC projection) and/or neighbours found in a
recipient-like bank should raise held-out block-sign accuracy by ≥ +0.005 (unweighted) on
Jurkat AND HepG2 with bootstrap CI above zero, without losing RPE1.

## Exact inputs (paths)
- Banks (`delta.npy` = pert − control, `perts.txt`, `genes.txt`, `b_ctrl.npy`):
  `data/donors/{k562_gwps,k562_essential,nadig_jurkat,rpe1,orion_hct116,orion_hek293t}`,
  `experiments/166-donor-acquisition/donors/nadig_hepg2`.
- GWPS PCs: `data/donors/k562_gwps/svd_U.npy` (9787×800), `svd_S.npy` (8248,), `svd_Vt.npy` (800×8248).
- Frame: copy of `experiments/398-knn-tilt-lolo/lolo398d.py` (LOLO, 370 block rule = top-2200
  by held bank's own control among GWPS-shared genes minus own targets; truth = held delta
  sign; votes = six-bank family minus recipient; shrink-to-marginal; per-target paired bootstrap).
- Known answers: `experiments/398-knn-tilt-lolo/work/orion_m/results398.json`
  (kNN25_shrink400: Jurkat 0.57525 / HepG2 0.58579 / RPE1 0.59881;
  marginals 0.54171 / 0.56410 / 0.59702) and 399 PREREG §Evidence table.

## Estimator
- Recipients L ∈ {nadig_jurkat, nadig_hepg2, rpe1}; vote banks = six-bank family minus L;
  shrink 400 (sweep 200/800 for the winner); targets = L perts present in GWPS.
- Similarity variants for neighbour search in GWPS (K=25 first):
  (a) Pearson, all 8248 genes (baseline); (b) Pearson on the query row's top-500 |delta|
  genes (per-query subset, re-centred/re-normalised); (c) cosine on rows z-scored per gene
  across GWPS (removes the common response); (d) Spearman (average-tie ranks per row);
  (e) Pearson on rows projected onto GWPS top-50 PCs (`U[:,:50]*S[:50]`, re-normalised).
- Source variants (best similarity + baseline similarity): neighbours by GWPS rows (S1) vs
  by a recipient-like bank's rows over GWPS∩S genes, pert overlap only
  (S2: Jurkat→nadig_hepg2 [same study/protocol], HepG2→orion_hct116 [per brief example],
  RPE1→orion_hek293t [epithelial]; swaps Jurkat→k562_essential, HepG2→hek293t, RPE1→hct116
  if time) vs union of the two top-25 lists (≈50 votes, compared against K=25 and K=50
  single-source arms to hold vote count fair).
- K ∈ {10,25,50,100} for the best similarity; shrink ∈ {200,400,800}.
- Controls: marginal; random-K (same K, seed 20260916+1).
- Endpoints: unweighted block accuracy (fid reads this) + |delta|-weighted (nmae proxy);
  paired per-target bootstrap, 10,000 resamples, Δ vs baseline (a) with 95% CI;
  median top-1 neighbour similarity per metric (diagnostic for the "0.11" claim).

## Controls / falsifiers
- Random neighbours score below marginal (re-checks 398: gain = neighbourhood, not votes).
- Union arms are only credited vs K-matched single-source arms.
- A similarity that wins by re-ranking but whose median top-1 sim is lower is flagged, not trusted.

## Kill criterion
No variant clears +0.005 unweighted on BOTH Jurkat and HepG2 with CI above zero (and no
significant RPE1 loss) → neighbourhood stays as shipped; report the negative.

## Logistics
`NUMBA_CACHE_DIR` not needed (numpy only, no `src/`/`scripts/` imports). mmap all
`delta.npy`; peak extra dense ≈ 1× 9787×8248 float32 (323 MB, Spearman/zscore, freed after
use). One job at a time; baseline ≈ 12 s, full sweep budgeted < 60 min.
