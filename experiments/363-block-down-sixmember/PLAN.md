# 363 — Block-down six-member pricing (PLAN)

## Hypothesis
361 showed our 2,200-gene CPM block (84.6% of our calls) has delivered sign
accuracy 0.508 vs the reference — a coin flip — and that forcing it all-down at
the *delivered-lfc* stage buys exactly the line's block down-tilt (H1: 51.8% →
+0.0081 raw fid), with nmae −0.0083 scaled and reach +0.018 scaled (noisy), and
pds/mse UNPRICED. The PI measured the same tilt on donor banks' own
high-expression genes: Jurkat (= context A) 0.540, K562 0.511, H1 0.516,
HCT116 0.507, HEK293T 0.502. If the tilt mechanism transfers, the fid gain on
A could be ~4x H1's. But an all-down block is a COMMON component across all 85
predictions, and 155 found a uniform down-shift cost pds 4x what it bought on
reach. This experiment prices, with actual emitted counts through the 323
evaluator, what the design-stage sign intervention costs pds and mse, and
checks whether fid/nmae/reach move as 361's exact arithmetic predicts after
emission.

## Exact inputs (read-only; never written)
- `experiments/323-current-family-validation/work/effects.npz` (`cap15`
  85x18533 on `genes` 18533-axis; `targets` 85 H1 panel) — design source.
- `experiments/323-current-family-validation/work/reference.npz` (`genes`
  18080-axis, `targets`, `keep` 10779-universe mask, `ref_lin` 10779 basal CPM,
  `q`/`lfc` reference DE, `b_ctrl`, `libraries`, `sums`, `jk`, `counts`).
- `experiments/323-current-family-validation/work/stats-cap15.npz`
  (delivered `q` 85x10779 for call-set membership check), `metrics-cap15.json`
  (known-answer row), `per-target-cap15.csv`, `h1-anchors.json` (H1-local
  scaling), `reference_de.parquet` (used inside scorer via `E.score`).
- Emission/scoring code imported read-only: `run323` (builder, map_axis),
  `evaluate323` (load_reference, block_stats, score, MEMBERS, TS),
  `report323` (mse_bootstrap, ci, scaled), `vcc2026.rc_emit`
  (rc_hybrid_block, CallSet(400,2200,1.9)), `vcc2026.donor_transfer`
  (emit_counts). No builder/script edits.

## Estimator
Per target, on the 10779-universe exactly as 361 (`fid361.py:157-165`):
head = top-400 universe genes by |D| where D = cap15 mapped universe-by-name
(2 universe genes absent from the effects axis zeroed); block = top-2200 by
`ref_lin` among non-head universe genes. Applied to the 18533-axis cap15 row
by gene name (absent genes are no-ops):
- `cap15` — unmodified (known answer).
- `block_down` — design[block] = -|design[block]|.
- `block_random` — same block, iid 50/50 random sign x |design| (seed 20260910),
  the "re-randomising the block" control.
- `block_up` — design[block] = +|design[block]| (fid mirror control).
Magnitudes preserved elementwise (|modified| == |cap15| on the block;
off-block bit-identical). Head/block membership is |.| -invariant so the
partition cannot shift under the intervention. Each arm then mapped
18533->18080 via `run323.map_axis` and emitted+scored with 323's exact
parameters (builder `R.builder()`, make_rngs(0) fresh per arm, min_cpm 10,
400 cells, spread block tail 1.0 keep_top 8000 keep_spread 0.25 lo 0 hi 0.60,
dispersion 0.1 overdispersion 4.0, key = where(gate, effect, 0),
jitter 1000+i, call_set_first CallSet(400,2200,1.9), emit_counts) and
`E.block_stats`/`E.score`. Outputs (`stats-<arm>.npz`, `per-target-<arm>.csv`,
`metrics-<arm>.json`) written ONLY to `experiments/363-block-down-sixmember/work/`.
Analysis: six raw members per arm; H1-local scaled overall with 323 anchors
(mse clipped [0,1], nmae floored -6, mean/6); paired 10,000-target bootstrap
(seed 20260909) of each arm minus cap15 per member with mse as whole-panel
ratio (claim cap recomputed per resample, report323.mse_bootstrap identity
gated at 1e-10); per-target fid delta vs per-target block down-tilt
(extrapolation table for the PI); design-stage flip fractions.

## Controls / falsifiers
1. Known answer: `cap15` arm must reproduce `metrics-cap15.json` EXACTLY —
   emitted-count digest, panel-total digest, all six raw metrics, mse legs —
   plus stats arrays and per-target rows (342's hash-checked pattern). If it
   fails, stop; changed arms are not interpreted.
2. Block reconstruction check: per target, verify the reconstructed 2,200-gene
   block against `stats-cap15.npz` q<0.05 membership (361: whole block called,
   2200/2200). Report match rate (median/min block-calls-in-block).
3. `block_up` must be the fid mirror of `block_down` on H1 (down-tilt 51.8% ->
   up arm loses ~= what down gains, up to the 49.5% flip-rate asymmetry).
4. `block_random` discriminates mechanism: if it costs pds as much as
   `block_down`, the price is re-randomisation of tiny effects, not the common
   component — report which.
5. Bootstrap identity gates (mse 1e-10, per-target means 1e-12) must pass or
   the contrast table is void.

## Kill criterion (board-worthiness for context A only)
Board-worthy ONLY if H1 pds AND mse raw losses, expressed in scaled-overall
units (Delta_raw/span/6 with H1 spans pds 0.44367, mse 0.96203), are EACH
smaller than half the projected A-side fid gain (+0.03 raw fid on one of three
contexts ~= +0.033 scaled on the panel mean ~= +0.0055 overall). I.e. need
pds loss < ~0.0055 overall AND mse loss < ~0.0055 overall (equivalently raw pds
loss < ~0.0146, raw mse worsening < ~0.0318), else NOT board-worthy even if fid
moves as predicted. State the measured pds/mse price in those units either way.
A clean negative with controls is a good result.

## Compute / memory plan
`memory_pressure | tail -1` >= 50% free before each arm; sequential arms only
(~142 s, ~5 GiB each; cap 12 GiB RSS; NUMBA_NUM_THREADS=4
OMP_NUM_THREADS=4); no full-matrix densification (per-target blocks only);
delete no 323 files; scratch arrays freed per arm with gc. No downloads, no
reserve data, no submissions/builders, no git writes.
