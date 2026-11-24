# 358 — sign-classifier-lolo (PLAN)

## Hypothesis

The SIGN of a held-out line's response of gene g to knockdown of target t is
predictable from label-free features — donor-line consensus sign and agreement,
the recipient's basal expression of g and of t, donor-vs-recipient basal shift,
plus a regulon prior — by a shallow classifier whose parameters are shared
across lines, trained only on the other five lines' observed responses (never on
the recipient's). This is worker 353's proposal 2. It differs from 341 (which
used the recipient's own perturbation labels) and from 332/340 (which regressed
magnitudes and were judged by retrieval).

## Exact inputs (paths)

- Donor banks (read-only, mmap): `data/donors/k562_gwps`,
  `data/donors/orion_hct116`, `data/donors/orion_hek293t`, `data/donors/rpe1`,
  `data/donors/nadig_jurkat`,
  `experiments/166-donor-acquisition/donors/nadig_hepg2` — files
  `delta.npy` (response), `b_ctrl.npy` = log1p(50000·p) control profile
  (verified in `scripts/build_donor_table.py:150` and
  `scripts/build_orion_donor.py:112`), `perts.txt`, `genes.txt`, `meta.json`.
- Frame: `experiments/334-text-conditioned-transfer/selection.json` (six lines,
  1,024 targets, 512 coordinate genes; targets ∩ genes = ∅ by construction).
- Target blocks: `order(targets, "334-split|")` even/odd halves (same salt and
  halves as 334/341, so target leakage matches the record).
- Optional prior: `data/external/omnipath-dorothea/dorothea-abc.tsv` (7.4 MB;
  TF→target edges with stimulation/inhibition flags). Included only as
  static per-(t,g) indicators; if parsing exceeds 15 min it is dropped and the
  omission recorded.
- Baselines to reproduce (341 RESULTS.md "Measured result"): source-mean
  weighted sign .589398, cross-line diagonal weighted sign .589398, own cosine
  .195532/.198735, retrieval .790095/.794614 (six lines, 1,024 targets,
  512 coordinates, pooled over 6×1024 target rows).

## Estimator (own pipeline, no 323 run, no builder)

Per outer held line L (6) × target block b (2; train = other half, test = this
half), mirroring 341's leakage regime (donor test-half responses never touched):

- Response matrix: unit-RMS directions `direction = raw / sqrt(mean(raw²))`
  per (line, target), float64, own reimplementation of `rms_direction`.
- Train rows: all (inner recipient D ∈ 5 donors) × (512 train targets) ×
  (512 genes) = 1,310,720 rows. Consensus/agreement for a D-row are computed
  over the other 4 donors (excluding L and D) — leakage-matched; including D
  would leak the label into the consensus. Agreement is stored as a fraction
  so train (/4) and test (/5) share one semantic.
- Features per (recipient R, target t, gene g), all label-free w.r.t. R:
  1. `cons` — mean unit-RMS response over donors; 2. `sign(cons)`;
  3. `|cons|`; 4. `agree_frac` — donor sign-agreement fraction;
  5. `basal_g` — R's `b_ctrl[g]`; 6. `basal_t` — R's `b_ctrl[t]`
  (109–138/1024 targets missing from `genes.txt` in 4 banks; line-median
  imputed + `miss_t` indicator = feature 7);
  8. `donor_mean_g` — mean donor `b_ctrl[g]`; 9. `diff_g = basal_g − donor_mean_g`;
  10–14. (prior, DoRothEA ABC: 7.4 MB, loads 0.1 s) `is_tf(t)`, direct edge
  t→g presence + net sign (50/524,288 pairs; 17 signed), `log1p(shared
  regulator count)` (4.1% nonzero), edge g→t presence + sign (82 pairs;
  20 signed). 15 features total; control B uses 1–4, diagnostic B2 adds 10–14.
- Label: `sign(R's observed unit-RMS response)`; sample weight `|response|`
  (334 defines no noise floor, so all 512 coords, weighted — 343's functional).
- Fit 1 (baseline): cross-line diagonal `cons_5 × gain` (gain = 334/341
  per-gene diagonal on train targets) — the number to beat.
- Fit 2: L2 logistic regression (standardised, lbfgs) on all features.
- Fit 3: `HistGradientBoostingClassifier` on all features
  (max_iter=150, lr=0.05, max_leaf_nodes=15, max_depth=4, min_samples_leaf=64,
  l2=10, no early stopping — the classifier mirror of core334.model).
- Delivered object for endpoints: `sign_pred × |cons_5|` (magnitude-preserving,
  §8.1e-style: primary sign endpoints depend only on signs; retrieval guard
  differs from diagonal only on flipped signs, so the guard is fair).
- Endpoints per held line (assembled over both blocks = all 1024 targets):
  |response|-weighted sign accuracy (primary), top-64 weighted sign
  (core341.sign_endpoints reimplemented), retrieval guard (core334.endpoints
  reimplemented). Pooled mean over 6×1024 + per-line means; paired target
  bootstrap (10,000) + crossed line/target interval (core341.interval
  reimplemented, seed 3580909).

## Controls / falsifiers

- Control A (recipient carries nothing): fit-2 features with gene-indexed
  recipient columns (`basal_g`, `diff_g`) permuted by one fixed derangement
  over the 512 genes per (outer, block), applied identically in train and
  test; `donor_mean_g` untouched; `basal_t` untouched (target-indexed —
  disclosed as a limit).
- Control B (donor-only): logistic on features 1–4 only.
- Control C (structure memorisation): logistic on full features with labels
  shuffled across the 512 train targets within each training line (one fixed
  derangement per (outer, block), as in 341's `label_permutation`).

## Kill criterion (mandate-verbatim)

DEAD if (fit 2 or 3 pooled weighted sign − diagonal) ≤ +0.005, OR control A is
within its interval of the fit (recipient features carry nothing), OR control C
explains ≥ half the gain. ALIVE if the gain is ≥ +0.015 with control A resolved
below it — then also report the retrieval change so the PI can price the trade.
No A/B/C forecast will be made; H1-frame and six-line-frame numbers are
development evidence only; no retrieval/cosine/proxy gain will be called a
score gain.

## Budget

≤ 5 h wall; ≤ 12 GiB RSS (expected peak < 3 GiB: direction array 25 MB,
per-fit float32 feature block ~70 MB, freed per block); threads
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4` + `threadpool_limits`. No downloads.
Scratch cached in `work/` (basal-target index, prior edge sets); deleted if large.
