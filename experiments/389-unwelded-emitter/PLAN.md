# 389-unwelded-emitter — PLAN

## Hypothesis

The scorer's DE table (LFC magnitude/sign → nmae/fid/reach/jac) is computed from
per-cell ARITHMETIC CPM (`mean_calc: arithmetic`, target_sum 1e6; direction = sign of
predicted log2FC), while pds/mse come from BULK per-gene sums. The production emitter
(`rc_hybrid_block`) welds every cell to one library total T = S // n, which makes the
arithmetic mean a deterministic function of the per-gene total — the DE channel is
welded to the bulk channel, so every DE-side steering idea pays a pds tax through
the bulk. With K heterogeneous per-cell library classes (levels near the real
controls' 25/50/75% quantiles, total S unchanged), per-gene totals stay pinned
(pds/mse fixed) while the per-gene arithmetic mean becomes steerable by placement:
counts in LOW-library cells raise a gene's arithmetic CPM, HIGH-library cells lower
it. Price: (S1) the lattice/|z|-drift cost of heterogeneity alone; (S2) the gain from
steering the arithmetic LFC on the 2,600-gene call set at 343's nmae-optimal 1.48x.

## Exact inputs (all read-only, never edited)

- Effects: `experiments/323-current-family-validation/work/effects.npz` (`cap15`,
  85x18,533; sha-checked against `work/compose.json` as in run385/run342).
- Reference: `E.load_reference()` from `experiments/323-current-family-validation/work/`
  (`reference.npz`, `reference_de.parquet`, `reference.json`).
- Builder: `R.builder()` = `scripts/make_submission_pm_cd4.py` via `run323`
  (`profile_from_effect`, `emission_spread_vector`, `emission_rank_axis`, `make_rngs(0)`,
  `CallSet(400, 2200, 1.9)`).
- Primitives imported from `src/vcc2026/rc_emit.py` (read-only import):
  `wilcoxon_z`, `design_gene`, `histogram_of`, `assign_cells`, `repair_rows`,
  `phi_tables`, `EME`, `de_from_block`, `CallSet`, `rc_hybrid_block` (K=1 reference).
- Scorer: `evaluate323.block_stats` (= `exact_de323.de_from_counts`, the scorer's
  Wilcoxon on per-cell CPM at each cell's OWN library total) + `evaluate323.score`.
- Anchors for H1-local scaling: `experiments/323-current-family-validation/work/h1-anchors.json`
  + `report323.scaled/mse_bootstrap/ci` (same as 342).
- Class levels: measured by me from `data/harness/h1-test.h5ad` control libraries
  (25/50/75% quantiles); assignment text quotes context-A values for A/B/C only.
- 388 status (2026-09-13): `experiments/388-de-channel-prize/` contains only `work/`,
  no RESULTS.md → S3 has no named estimator; S3 will be skipped, shuffled control runs.

## Estimator

New `rc_hybrid_block_k` in this directory (`k_emit389.py`), same pass 1 and sub-gate
block as `rc_hybrid_block` (same `emit_counts` calls, same rng streams → per-gene
totals bit-identical, asserted), then:

- Cell totals: classes = list of (fraction, raw level); cell counts by largest
  remainder; levels scaled to sum exactly S (integer fix-up like the T+1 remainder
  convention); asserted `block.sum(1) == cell_totals`, `cell_totals.sum() == S`.
- Design: per-class Phi tables `v=log1p(1e6*c/L_k)` (the equal-T count<->CPM
  monotonicity generalises per class; documented in code). Per-gene total split
  across classes (S1: even per-cell counts; S2: greedy split to hit the arithmetic
  target). Within each class the existing `design_gene`/`histogram_of` family scan
  runs unchanged against that class's tables with per-class Fbar target
  `0.5 +- |z*|*sd_global`; sign-pattern choice (2^K, K=3 -> 8) by pooled realised |z|
  via a new grouped `hist_z_k` njit kernel (pairs grouped by exact CPM value via
  integer cross-multiplication `c_a*L_b == c_b*L_a`, because two (count,class) pairs
  can share one CPM value). Realised z re-measured exactly with `wilcoxon_z` on the
  materialised per-cell CPM matrix.
- Placement: `assign_cells` per class (same call sequence as `EME.place`, so K=1
  issues identical RNG draws → bit-identical block), global `repair_rows` unchanged.
- Steering (S2): per-gene arithmetic target `A_t = ref_lin * 2^(s*LFC_design)` with
  `LFC_design = log2(bulkCPM/ref_lin)`, `s` in {1.48, 1.25, 2.0}, applied to
  prescribed-call-set genes only (`CallSet(400,2200,1.9)` membership per target);
  complement genes use the S1 split. Unreachable targets clip; clip rate reported.
- Harness `run389.py` (adapted from `run385.py`/`run342.py::emit_arm`): emits + scores
  arms, writes `work/stats-<arm>.npz`, `work/per-target-<arm>.csv`,
  `work/metrics-<arm>.json`. `report389.py` (adapted from `report342.py`): paired
  10,000-target bootstrap, seed 20260913, H1-local scaling with 323 anchors.

Arms: `cap15` (production `rc_hybrid_block`, S0 gate) → `k1` (new fn, K=1,
None; must match cap15 bit-for-bit, else within 7e-5 seed spread with the RNG cause
named) → `k3plain` (S1) → `k3s148` (+`k3s125`,`k3s200` if time, S2) → `k3shuf`
(shuffled-class control: same library multiset, steering targets deranged).

## Controls / falsifiers

1. S0 known answer: `cap15` must match `work/metrics-cap15.json` raw members
   (fid 0.49755452, jac 0.07475983, mse 0.97494624, nmae 0.95083646, pds 0.88963585,
   reach 0.25434726) + both digests, else stop.
2. Invariant asserts per block: per-gene totals == `rc_hybrid_block`'s; per-cell
   totals == class levels; delivered steered-gene arithmetic mean within tolerance
   of target (logged, clip rate reported); achievable (A_min,A_max) range per gene
   computed analytically (greedy class fill) and reported as distributions.
3. Shuffled-class control must not gain (else the effect is not the LFC channel).
4. pds must be bit-identical across arms sharing totals (bulk-pinned proof); mse
   allowed only jackknife-placement wobble (reported, not gated).

## What kills the idea (from assignment)

- S1 alone loses >= 0.005 H1-local overall with a CI excluding zero → stop after S1,
  mechanism named (lattice error and/or |z| drift and/or BH-count move).
- Reachable |dLFC| at K=3 < 0.1 log2 for typical call-set genes (analytic range
  statistics, computed before the S2 panel) → report ceiling, stop.

## Compute / memory / rules

- ONE harness run at a time (~140 s / ~5 GiB each); `memory_pressure` free >= 40%
  before every heavy step; `NUMBA_CACHE_DIR=experiments/389-unwelded-emitter/work/numba-cache`,
  `NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`; <= 12 GiB RSS.
- Write only in `experiments/389-unwelded-emitter/`; no git, no builder, no submit.
- H1-panel numbers are development evidence, not A/B/C forecasts; no board conversion.
