# 357 — Jurkat shRNA truth proxy (PLAN)

## Hypothesis
Context A is Jurkat (041). LINCS holds 256 JURKAT shRNA knockdown signatures
(worker 356). If ≥8 of those genes intersect panel A targets, each signature is a
noisy same-line, same-target proxy for the true response direction in Jurkat on the
978 landmark genes — a truth proxy for exactly what 343 says we lose on
(per-gene direction on the held-out line at the design stage).

## Inputs (exact paths)
- `data/external/cmap/siginfo_beta.txt` (465 MB; read chunked, never whole into RAM twice)
- `data/external/cmap/level5_beta_trt_sh.gctx` (11.7 GB; h5py, read only needed columns)
- `data/external/cmap/level5_beta_trt_xpr.gctx` (6.5 GB; inventory only)
- `data/external/cmap/geneinfo_beta.txt` (landmark map, `pr_is_lm == 1`)
- `data/challenge/pert_counts.csv` (panel A/B/C targets — NOTE: single panel file;
  `data/challenge/context_{A,B,C}.h5ad` + manifest checked for per-context panels)
- `data/challenge/gene_names.csv` (18,533 challenge axis)
- `experiments/323-current-family-validation/work/reference.npz` (`targets`: 85 H1 eval targets)
- `experiments/323-current-family-validation/work/effects.npz::cap15` (319-composition
  design rows on H1 panel, fallback if no A-panel design array exists)
- A-panel design rows: search `experiments/319-cap-variant/work/build*/`,
  `experiments/312-normrel-arm/work/build*/` for 300×18,533 arrays; read-only read of
  `scripts/make_submission_pm_cd4.py` for the artifact name (never run it).
- H1 truth for step 5: `experiments/323-current-family-validation/work/reference.npz`
  (`lfc`), `reference_de.parquet`.
- Known answer: `experiments/323-current-family-validation/work/metrics-cap15.json`
  (raw pds 0.8896358543417369, mse 0.9749462418923823, nmae 0.9508364630325662,
  reach 0.25434726355678233, fid 0.4975545182292535, jac 0.07475983298131939)
  and/or 343's recomputed mean NMAE 0.950836 over 58 eligible targets.

## Estimator
1. Inventory with pandas chunks: JURKAT `trt_sh` rows (sig_id, pert_iname, cell_iname,
   pert_time, pert_dose, distil_cc_q75/tas if present); intersect pert_iname with
   panel targets and H1-85 targets. Count JURKAT `trt_xpr`. Same for HELA `trt_sh`.
2. h5py: open level-5 gctx, read row (gene) metadata once, column (sig) metadata once,
   then read ONLY the columns for intersecting sig_ids, rows with `pr_is_lm == 1`.
   Map landmark symbols → challenge axis.
3. Locate A design rows (see above). If absent: report exact recipe needed, do NOT
   restore the 12.6 GB h5ad; fall back to H1-panel cap15 design rows for the
   same targets where they exist (context-independent consensus rows — verify this
   claim from builder code before using; otherwise H1-calibration only).
4. Per intersecting target: |z|-weighted sign agreement
   agree = Σ|z|·1[sign(design)==sign(z)] / Σ|z| over landmark genes (z≠0),
   vs (i) 100 wrong-target design rows, (ii) cross-target common vector,
   (iii) constant all-down vector; fraction of landmarks with |z|>2.
   Pool across targets + 1,000-draw target bootstrap.
5. H1 calibration: same agreement of Jurkat shRNA signature vs 323 H1 reference lfc
   on the same targets' landmark genes; compare with our H1 design-row agreement
   to the same truth (is the same-line proxy comparable to Perturb-seq truth?).

## Controls / falsifiers
- Wrong-target design rows (100 draws): proxy must beat the null distribution.
- Common-vector and all-down controls: a "better than chance" claim must survive them
  (094/155: reference pools skew down; all-down is a strong null).
- |z|>2 fraction: if the shRNA signal itself is weak, agreement is noise.
- H1 calibration (step 5): if Jurkat-shRNA↔H1-truth agreement ≈ chance, the proxy
  cannot adjudicate our A rows even if step-4 numbers look good.
- Known-answer reproduction before any new measurement (rule 3).

## Kill criteria (from mandate)
- <8 panel-A targets intersect → proxy too thin; report count and stop.
- If no A design array on disk → report exact recipe, no h5ad restore; continue only
  with what exists on disk (H1-panel rows + H1 calibration).
- A clean negative (agreement ≈ wrong-target null) is a good result.

## What would kill the idea
Agreement of our A rows with Jurkat shRNA signs within the wrong-target null band,
or shRNA↔H1-truth agreement at chance (proxy invalid), or <8 intersects.

## Budget / memory
~3 h. `memory_pressure` ≥50% free before >4 GiB steps; cap 12 GiB RSS;
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`; chunked reads; delete scratch.
No downloads (all local). No reserve paths. No builders run.
