# 393 — mse shared component (PLAN)

## Hypothesis

`mse`'s signal term X splits into a shared part (mean effect across the
context's 300 targets) and a centred part. In 187's per-context identity

```
MSE_c − 1 = [s²A_c + Kc_c − min(claim_c, s²B_c) − 2sX_c] / ΣD
```

B is the across-perturbation CENTRED sum of squares, so A − B = P·‖d̄‖² is
the shared component's energy and the budget cliff reads only the centred
part. With separate scales c (shared) and s (centred) the numerator becomes
c²·P‖d̄‖² − 2c·X_shared + s²B − min(claim, s²B) − 2s·X_centred + Kc: the
shared scale has its own optimum c* = X_shared / (P‖d̄‖²), independent of the
cliff. If our shared component is poorly aligned with the truth's, c* < 1
and shrinking it improves mse; pds (a retrieval rank) may or may not pay.

## Exact inputs (all read-only, never written)

- KA-i: `experiments/187-joint-amplitude/mse187.py` (C_p-free LS solve on
  a8+a10+a11 → X = 0.019240, ΣD = 0.229141, residuals ≤ 6.1e-5) and its
  inputs `experiments/096-mse-floor/work/{a8,a10,a11,a12}-stats.npz`
  (+ `experiments/096-mse-floor/mech_a12_audit.py`), `experiments/180-def-readiness/work/knee180.json`.
- KA-ii (211): `experiments/211-amplitude-reopened/work/chat211.json`
  (q_c per emitter) and `bx211.json` (A_c/B_c on the shipped artifacts);
  q_c(a15) = 1.000000 / 0.999904 / 1.000000.
- Task 1: shipped 319p2 artifact
  `experiments/319-cap-variant/work/build_p2/the-photocopier-demands-a-raise.h5ad`
  (11.8 GB; obs `context`, `target_gene`; 900 blocks × 400 rows),
  real controls `data/challenge/context_{A,B,C}.h5ad`, gene axis from the
  artifact's own `var/_index`.
- Task 2: board snapshots `data/leaderboard/2026-09-11T1910Z.json`,
  `2026-09-13T0544Z.json` (paths TBD — exact filenames resolved at run
  time; the two dated snapshots named in the mandate), full records
  `data/board-entries/*.json`, build metas `submissions/*.h5ad.meta.json`,
  `scripts/restore_prediction_h5ad.py`, `docs/STORAGE-AUDIT-2026-09-09.md`.
- Task 4: H1 harness `experiments/323-current-family-validation/`
  (`evaluate323.py`, `exact_de323.py`, `work/effects.npz` key `cap15`,
  `work/reference.npz`, KA = `work/metrics-cap15.json`);
  `experiments/385-per-target-callset-size/run385.py` as the exact cap15
  chain template.
- Scorer (installed, never vendored): `cell_eval2.prep.bulk_lognorm_means`,
  `cell_eval2.metrics.delta._row_weights`, `_across_pert_budget`;
  score_mse = clip((0.987862879 − MSE_c)/0.950162566, 0, 1) averaged
  (confirm against 187's `score()`).

## Estimator

1. KA: re-run 187's `solve_ls` on a8+a10+a11 via my own code (import 096's
   `measure` read-only, or re-implement from `mse187.py` logic) and check
   X = 0.019240 ± 2e-5, ΣD = 0.229141 ± 2e-4; reproduce 211's q_c from
   `chat211.json` / recompute claim saturation on a small subsample.
   STOP if either fails (mandate rule).
2. Task 1: stream the 319p2 h5ad per (context, target) block with h5py
   (bincount per block; never densify). Per block: pseudobulk total,
   mu = log1p(5e4 · tot/tot.sum()), d = mu − b_ctrl(context) with b_ctrl
   from the pooled real controls of that context (same construction as
   210's `control_frames`, verified equal). Per context: w =
   `_row_weights`, A_c = Σ w_p(‖d_p‖² − d[p,target_p]²),
   B_c = `_across_pert_budget`, P‖d̄‖² via the same weights, shared share
   (A−B)/A, cos(d̄_A,d̄_B), cos(d̄_A,d̄_C), cos(d̄_B,d̄_C).
3. Task 2: tabulate our scored board entries with full records + build
   metas; for consecutive arms that changed predicted effects, compute
   ΔP‖d̄‖² and ΔB from restorable artifacts (≤ 2 restores, deleted after
   use). Fit X_shared per context if identified; else bracket
   X_shared ∈ [0, P‖d̄‖²] (c* ∈ [0, 1]).
4. Task 3: closed-form Δ raw mse for c ∈ {0, 0.25, 0.5, 0.75, 1, 1.25}
   across the bracket (panel mean over contexts, board units).
5. Task 4 (≤ 2 harness runs, one at a time, gated on free ≥ 45%):
   run 0 = cap15 KA reproduction (must match metrics-cap15.json before any
   arm is scored); run 1 = cap15 with the design's mean-across-85-targets
   component scaled to 0.5 before rendering (and 0.0 only if it fits in the
   same run script without a third run — else dropped and stated).
   Read pds + raw mse legs on H1 (mse score clamps near 1 there).

## Controls / falsifiers

- KA-i/ii must PASS or the whole experiment stops (stated numbers above).
- Task 1 sanity: A_c/B_c on 319p2 must land within a few % of 211's
  rc_callset values (A 0.23643/0.22643/0.23662, B 0.23434/0.22433/0.23448);
  B/A must sit in 0.990–0.999; weight_sum must equal 300/18532 ≈ 0.01619
  if all targets resolve.
- Task 3 bracket endpoints are analytic bounds, labelled as such — not
  measurements of X_shared.
- Task 4 KA run must reproduce the emitted-count digest + panel digest +
  six members within 323 tolerances, else no arm result is reported.

## Kill criterion (from the mandate)

KILL the idea if the bracket's best case gains < 0.01 raw mse, or the H1
pds cost at c = 0.5 exceeds −0.01 raw (i.e. pds drops by more than 0.01).
Either condition alone kills it. A clean negative with controls is the
expected deliverable unless both gates pass.
