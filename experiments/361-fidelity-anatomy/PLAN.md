# 361-fidelity-anatomy — PLAN

## Hypothesis
`fid = k/max(n_pred,n_real)` is dominated by the ~2,200-gene CPM block whose genes are
chosen by basal expression, not effect — so their delivered signs may be ~coin-flip, and any
sign prior with >50% accuracy on those genes buys fid. The question is the price on nmae and
reach (same signs, reference gate) — measured here by exact sign-flip counterfactuals on
saved arrays. pds is NOT priced (needs emission + 323 run).

## Inputs (all read-only, exact paths)
- `experiments/323-current-family-validation/work/stats-cap15.npz` — our q, lfc per
  target×universe gene (85×10,779), sums/jk unused.
- `experiments/323-current-family-validation/work/reference.npz` — reference lfc, q over the
  same universe; `keep`, `genes`, `targets`, `ref_lin` (basal CPM).
- `experiments/323-current-family-validation/work/effects.npz::cap15` (+ `genes`, `targets`) —
  design effects on the 18,533 axis, mapped by gene name to the universe (2 universe genes
  absent from the effects axis → 0, per 347).
- `experiments/323-current-family-validation/work/per-target-cap15.csv` — recorded per-target
  fid/nmae/reach (known answers: fid mean 0.49755452 over 85; nmae 0.95083646 over 58;
  reach mean 0.25434726 over 82 non-NaN).
- `experiments/323-current-family-validation/work/metrics-cap15.json` — recorded means.
- `experiments/323-current-family-validation/work/h1-anchors.json` — local baseline/replicate
  anchors for scaled conversion.
- `experiments/348-reach-ordering/reach348.py` — exact reach walk (copied logic into our
  script, read-only import not required; we re-implement and verify to 1e-9).
- Scorer source for adjudicability/conventions:
  `.venv/.../cell_eval2/metrics/direction.py` (`_components`: scored set =
  {q_pred<0.05}∩{adjudicable}, own target removed; `F=k/max(n_pred,n_real)`,
  n_real=0 & n_pred>0 → k/n_pred) and `evaluate323.py score()` (prepare_de, α=0.05).

## Estimator (exact arithmetic, no emission, no 323 run)
1. Known-answer: implement fid per target from (q_pred, lfc_pred, lfc_ref, q_ref) with
   adjudicable = finite & nonzero ref lfc, own-target excluded; require max abs diff 1e-9
   vs recorded per-target fid. Recompute nmae per target (gate = ref-significant, finite ref
   lfc, own target removed, ≥10 genes else NaN) and reach per target (scorer ordering:
   sig-first, q_pred, p_pred, |lfc| desc, gene name; walk purity≥0.9 deepest k*/n_real) and
   require match to recorded means (nmae 1e-9, reach 1e-9 per target).
2. Decomposition: per-target n_pred, n_real, k; share of our calls that are
   reference-significant vs not + k-accuracy on each; head (top-400 |design|) vs CPM block
   (top-2200 ref_lin among non-head) vs other-called split of our calls (sign accuracy,
   contribution to F); reference sign tilt (frac ref-lfc<0) among block calls per-target and
   pooled, among all adjudicable genes, among ref-significant genes.
3. Counterfactuals (ONLY block-call signs change; head signs, call set, magnitudes, ordering
   keys unchanged; lfc'_block = ±|lfc| per rule):
   a. all-down; b. LOO reference-common-vector sign (truth ceiling); c. LOO design-common-vector
   sign (design-time); d. down on top-half block by basal CPM, design sign elsewhere;
   e. current. Recompute exact fid/nmae/reach per target; mean; paired 10,000-target
   bootstrap (seed 20260909, resample targets with replacement, per-member valid-cohort
   means like report342) of Δ vs current; H1-local overall Δ = (Σ scaled Δ)/6 with
   mse/pds/jac carried at 0 (stated, not measured).
   - LOO common vector: for target i, common_j = mean over the other 84 targets of the
     source row (ref lfc for b; design effect for c); sign(common_j); exact-zero → keep
     current sign (ties carry no information; count them).
   - nmae recompute: numerator Σ|lfc'−l| over the unchanged gate; denominator unchanged.
   - reach recompute: ordering keys (q,p,|lfc|) unchanged → same order; only match flags change.

## Controls / falsifiers
- Known-answer gate (§1) must pass before any counterfactual is read; if adjudicability or
  own-target handling differs, read `evaluate323.py`/`direction.py` until it matches, else
   STOP and report non-reproduction.
- Reach walk verified per-target to 1e-9 against recorded values (same bar as fid).
- Bootstrap identity: resample weight-1 row reproduces point means exactly.
- Sanity: counterfactual (a) block sign-flip rate and pooled tilt consistency (all-down
  accuracy on block ≈ pooled frac-down).

## Kill criterion (from mandate)
Dead if no design-time counterfactual (a, c, d) gains ≥ +0.008 raw fid with nmae and reach
not resolved negative (paired 95% intervals entirely harmful = resolved negative).
If one passes, specify the emitter change (sign of designed z on CPM-block genes) for a PI
323 arm; pds unpriced.

## Resources
Peak < 2 GiB (all arrays ≤ 85×10,779 float64, ~7 MB each); wall time minutes. Threads capped
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`. Writes only into `experiments/361-fidelity-anatomy/`.
No downloads, no sealed data, no builders.
