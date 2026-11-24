# 378-audit-a-arm-chain — PLAN

## Hypothesis
The 375 ctx-A arm artifact and its supporting claims (372 calibration, 374 retrieval, 373 patch
inertness) are correct as stated in the assignment. The audit re-derives each from scratch.

## Exact inputs (paths)
- Audited artifact: `experiments/375-ctx-a-arm/work/ctx_a_blend_rows_b035.npz` + `.json`
- 372 artifact: `experiments/372-jurkat-calibration-artifact/work/ctx_a_calibrated.npz` + `.json` + `fit372.json`
- Official axis: `data/challenge/gene_names.csv` (col `gene_name`), panel: `data/challenge/pert_counts.csv` (col `target_gene`)
- Donor banks: `data/donors/nadig_jurkat/{delta.npy,genes.txt,perts.txt}`, same for `k562_gwps`, `orion_hct116`, `orion_hek293t`
- Frame defs (read AFTER blind numbers for items 3/4): `experiments/341-known-line-calibration/core341.py`,
  `experiments/334-*/core334.py`, `experiments/341-known-line-calibration/selection.json`,
  `experiments/334-*/selection.json`; audited scripts `build372.py`, `run374.py`, 375 docs — read only after own numbers.
- Builder: `scripts/make_submission_pm_cd4.py` (HEAD) vs `experiments/373-blend-rows-patch/work/patched/make_submission_pm_cd4.py`
- Known answer: `experiments/323-current-family-validation/work/metrics-cap15.json` raw members
  (fid 0.4975545182292535, jac 0.07475983298131939, mse 0.9749462418923823, nmae 0.9508364630325662,
  pds 0.8896358543417369, reach 0.25434726355678233); 343 NMAE decomposition (0.951, mag-ratio 0.205,
  sign-acc 0.588/0.600) — reproduced by re-reading, not re-running the 142 s scorer (budget).

## Estimators (own code, from scratch)
1. Integrity: `numpy.load`, list keys/shapes/dtypes/finiteness, scalar strings, sha256 of file bytes;
   gene-axis membership vs official 18,533 list; intersection of `genes` with 300 panel symbols must be empty;
   `targets` ⊆ panel, |unique|==259; duplicates check.
2. Equality: load both npz, `max|375.matrix − 372.prediction_raw|` and vs `372.matrix`; exact + digest.
3. Ridge (blind): common axis = ∩(4 bank genes) ∩ canonical?? — approximation: ∩(4 banks) minus
   training targets minus 300 panel symbols; unit-RMS row normalisation; stack 3 donor deltas
   (3·n_train rows) as X, Jurkat delta rows as Y on same gene axis; ridge λ = n_rows (6018-class);
   solve (X'X+λI)W = X'Y per split; even/odd target split on Jurkat perts; weighted sign gain =
   mean over eval targets of (|y|-weighted sign agreement of same-line-corrected vs diagonal baseline)?
   — Will pin definition after blind first pass: report own gain + paired-target CI (paired t / percentile
   bootstrap 10k). Pass if within ±0.005 of +0.042.
   Note: exact endpoint definition lives in 341/core; blind pass uses |y|-weighted sign accuracy
   gain (same − diag), then read core341 and re-run exactly.
4. Retrieval (blind): cosine-rank across targets; for β∈{0,0.35} blend prediction_raw with consensus
   under fixed-norm rule, score retrieval on 341 frame; report Δ(0.35−0). Read core334/core341
   definition first (assignment explicitly permits), implement independently.
5. Patch: `diff -u` + `git diff --no-index`; classify each hunk: definition / conditional-on-flag
   (`os.environ.get("VCC_CTX_A_PRIOR")` guard) / config-key-conditional; flag any unconditionally
   executed behaviour change; verify inertness when unset by code inspection (no builder run per rules).

## Controls / falsifiers
- Item 1/2: any key/shape/dtype/sha mismatch, any panel symbol in `genes`, any target outside panel → MISMATCH.
- Item 3: shuffled-label control (expect gain ≤ 0); diagonal baseline sanity; CI must exclude 0 for a real gain;
  kill if |own − 0.042| > 0.005 after exact-definition rerun.
- Item 4: β=0 must reproduce the unblended baseline; kill if β=0.35 reduces retrieval (Δ<0 beyond noise).
- Item 5: any unconditional hunk → FAIL (not inert).

## What kills the idea
Any MISMATCH on 1/2/5 blocks the slot; a miss on 3/4 downgrades the supporting claim (not the artifact).

## Budget / memory
~1.5 h. `memory_pressure` free 78% at start; cap 12 GiB RSS; NUMBA/OMP threads 4; npz/donor matrices
are small (259×4951, banks ≤ ~2000×~7000); no .h5ad densification; delete scratch.
No downloads, no builders, no commits, no submits.
