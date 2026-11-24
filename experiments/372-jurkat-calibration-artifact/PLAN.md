# 372 — PLAN: Jurkat same-line calibration artifact for context A

2026-09-10 UTC. Assignment 372-jurkat-calibration-artifact. Numerical + engineering; ≤12 GiB; threads capped. Write-only in this directory.

## Hypothesis
346 priced a same-line ridge map (donor consensus → recipient response, fitted on recipient's own public non-evaluation targets, blended at fixed row norm) on H1: design-stage weighted sign +0.0348, H1-local overall +0.0071 at β=0.5 with pds cost −0.020 raw. The same construction for context A (= Jurkat, 041) using the public Nadig Jurkat bank yields a production-usable artifact: a calibrated Jurkat prediction for the A/B/C panel targets on a covered gene axis, with documented counts, λ, poison test, and held-out sign gain. This is an engineering artifact, not a score claim; no emission, no build, no submission.

## Exact inputs (paths)
- Recipient bank: `data/donors/nadig_jurkat/` — `delta.npy` (2,270 × 8,881), `genes.txt`, `perts.txt`, `b_ctrl.npy`. Raw cells `data/external/nadig2025/jurkat_raw_singlecell.h5ad` only if 346 needs a raw-derived quantity (it does not; recorded here as unused).
- Donor banks for X (same three production donors as 346 §"Deviation"): `data/donors/k562_gwps`, `data/donors/orion_hct116`, `data/donors/orion_hek293t` (delta/b_ctrl/genes/perts).
- Panel exclusions: `data/challenge/pert_counts.csv` (300 symbols), `data/challenge/gene_names.csv` (18,533 canonical order).
- Estimator primitives reused verbatim (import, not copy): `experiments/341-known-line-calibration/core341.py` (`fit`, `predict`, `sign_endpoints`), `experiments/334-text-conditioned-transfer/core334.py` (`rms_direction`, `source_blocks`, `diagonal`, `derangement`, `endpoints`). Loader pattern from 346 `fit346.py::bank` (exact-name row/column selection + scalar-index fixture).
- Known-answer record: 346 `work/fit.json` (1,965 training / 5,275 covered / 80 test), 346 RESULTS.md design sign 0.600185→0.634995, 323 `work/metrics-cap15.json` raw members, 341 Jurkat weighted-sign gain +0.03380218 vs diagonal.

## Recipe (346 exact, deviations stated)
1. **Gene axis rule (346 `build_axes` verbatim, minus Pan mask):** 346 covered = (donor genes ∩ canonical ∩ Pan-measured-mask) − training − panel/eval union. Jurkat bank has no `measured_mask.npy`; the axis is (3 donors ∩ Jurkat ∩ canonical) − training − 300 panel symbols. Document: common 4-bank ∩ canonical size, how many training symbols and how many panel symbols intersected the axis and were dropped, final covered size.
2. **Training targets:** Jurkat perts ∩ 3-donor perts, minus all 300 panel symbols. Assert Jurkat ∩ panel = 0. Report count (expected ≈2,006; NEXT-METHOD-SCREEN quoted 1,986 — record both, do not tune to match).
3. **Consensus definition:** 346 `rms_direction` unit-RMS raw directions; `baseline = direction[donors].mean(0)`; query = baseline[test]. Same.
4. **Ridge:** `core341.fit(x, y)` verbatim, λ = number of rows (3 × n_train), on RMS-scaled features; x = three stacked leave-one-out donor consensuses on training targets, y = Jurkat response tiled 3×. Same. Report λ, n_rows, n_features, normal-equation residual, subset-linearity identity (<1e-12).
5. **Test/query targets:** panel targets with a full 3-donor consensus, in `pert_counts.csv` order. DEVIATION — assignment says "100 A targets"; the panel is 300 shared across A/B/C (manifest.json), 259 have full 3-donor consensus, 0 are in the Jurkat bank. There is no buildable "100". Ship all predictable panel targets in panel order and document the deviation; do not invent a 100-subset.
6. **Matrix definition (ambiguity resolved here):** assignment's `matrix` ("the calibrated Jurkat prediction") + `beta:0.5` + `mode:'blend_rows'` + `norm_rule:'fixed_norm_on_covered'` is ambiguous between raw map output and blended rows. Ship `matrix` as the β=0.5 fixed-norm blend of donor-consensus direction and calibrated prediction (346 `effects346.blend` verbatim, with design row replaced by the donor-consensus row since no A design exists), scaled to the consensus norm; also store `prediction_raw`, `consensus`, `consensus_norm` as extra keys so the builder can reconstruct either. Document.
7. **Genes:** covered axis symbols, panel-target genes dropped (state how many of the 300 intersected the axis). No NaN allowed.

## Controls / falsifiers
- **Known-answer first:** re-run 346's three-donor axis counts with own code (expect 1,965 / 5,275 / 80 test / 5 uncovered AASS,CNN3,COL6A1,HMGCS1,MED23). If mismatch, stop.
- **Poison test:** refit with Jurkat rows for the query (panel) targets replaced by 1e6 noise (they do not exist in the bank; the assembly holds zeros there) → every gain/scale/beta and every prediction byte-identical. Proves panel exclusion.
- **Shuffled-label map:** `derangement(n_train, 346 SEED)` on Jurkat training labels, same architecture; evaluated on Jurkat's own held-out half (below). Must be ≈0 or negative gain.
- **Within-bank held-out sign check:** split the 2,006 training targets into halves (seed 3720910, sorted order split evens/odds for determinism + one random-half replication). Fit ridge on half 1, predict half 2's donor consensus, score `core341.sign_endpoints` weighted sign vs Jurkat truth; compare `same_line_full` vs `crossline_diag` (346's `P.diagonal` on inner blocks). Report gain to compare with 341's +0.0338 on Jurkat. This is a diagnostic, not a score.
- **Sanity gates (abort on fail):** no NaN/Inf; per-row cosine(calibrated, consensus) moderately positive (median in (0.1, 0.99), not ~1, not ~0); covered fraction reported; bank shas recorded; panel overlap asserted 0.

## Outputs
- `work/ctx_a_calibrated.npz`: `targets` (panel targets with consensus, panel order), `genes` (covered axis), `matrix` (N×G float64 blended), `prediction_raw`, `consensus`, `consensus_norm`, `beta=0.5`, `mode='blend_rows'`, `context='A'`, `norm_rule='fixed_norm_on_covered'`.
- `work/ctx_a_calibrated.json`: bank shas, λ, counts, poison result, held-out sign gain, cosine summary, coverage fraction, deviations.
- `RESULTS.md` (mandatory sections), incremental.

## Budget / safety
≤4 h. One ridge fit (≈6k rows × ≈5k features, ~30 s, ≤10 GiB like 346) + one held-out replication + artifact write. `memory_pressure | tail -1` ≥50% before fit; `NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`; h5py/backed reads, never densify single-cell matrices; delete scratch; no git, no submit, no builder, no reserve, no download (all inputs on disk).
