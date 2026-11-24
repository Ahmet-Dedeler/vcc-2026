# 347 — PLAN: sign-confidence strata on the H1 development panel

## Hypothesis
343 measured design weighted sign accuracy 0.600 on the reference gate, flat across
|lfc| and basal-CPM quintiles. A stratum with weighted accuracy ≥ 0.65 AND ≥ 10% of
|l|-mass would let 188's identity buy NMAE (scale up where e > ½, shrink elsewhere —
a surgical version of what 317 tested globally and lost 28–41:1). If no such stratum
exists among all honest pre-emission orderings, confidence-weighted amplitude is closed.

## Inputs (read-only, saved arrays only; no 323 run, no builder)
- `experiments/323-current-family-validation/work/reference_de.parquet` (gate truth: target,
  feature, p_adj, log2_fold_change)
- `.../work/stats-cap15.npz` (delivered per-target p/q/lfc on the 10,779-gene universe)
- `.../work/reference.npz` (universe axis `genes[keep]`, `ref_lin` basal CPM)
- `.../work/effects.npz::cap15` (DESIGN effect matrix, 85×18,533, mapped to universe by gene name)
- `.../work/metrics-cap15.json` + `experiments/343-*/work/{summary,signs}-cap15.json` (known answers)

## Estimator
Pooled gate-pair table over NMAE-eligible targets (gate = q<.05, finite nonzero lfc,
gene≠target, ≥10 genes/target; expect 58 targets). Per pair: ref l, delivered lh/p/q,
design ef. Stratum stats: n, mass share Σ|l|/Σ|l|_all, weighted + unweighted sign accuracy
(DESIGN primary, delivered secondary), delivered magnitude ratio Σ|lh|/Σ|l|. Exact per-target
NMAE recompute with delivered-lfc amplitude proxy: in-stratum ×{2,3} × out {×1, ×0.5};
report mean-over-targets NMAE and Δ vs base (~0.9508).

Strata: (a) design-|effect| deciles, per-target-within-gate and pooled-global; (b) delivered
called vs not + delivered-p pooled deciles; (c) design sign vs leave-one-out cross-target
common vector (mean of other 84 design rows): agree/disagree/≈0, where ≈0 = |common| below
the pooled 10th percentile of |common| (threshold reported); (d) donor-sign agreement —
EXPECTED SKIP: stage.npz stores only consensus arms + CD4 aggregates (verified: no per-donor
rows; donors `del`eted in run323.py:240; cd4_votes ∈ {0,3} only, a carried-by-all-or-none
flag, not a graded agreement count); will verify and state precisely; (e) gene class
RPL*/RPS*/MT-*/HIST* vs rest + reconstructed emitter 2600-set (head = top-400 universe genes
by |design effect| per target — the |H| key per rc_emit.py:839; rest = top-2200 by ref_lin
among non-head) vs outside; (f) (a)-top-3-deciles ∩ (c)-agree.

## Controls / falsifiers
1. Known-answer first: mean NMAE 0.950836; delivered w-acc 0.588 (0.638 called / 0.552
   uncalled); designed w-acc 0.600 (0.649 / 0.558); gate mass in own set 0.346; magnitude
   ratio 0.205. Must match to ~1e-6/1e-3 or STOP and report.
2. Decile sanity: pooled deciles must carry ~10% mass each only if |l| uniform — they won't;
   report mass share explicitly so a "high-accuracy, no-mass" decile can't pass the gate.
3. (c) circularity guard: common vector is leave-one-out (other 84 targets), so the target's
   own row never votes. Common-sign ≈ design-sign agreement is expected to be high (shared
   programmes); the question is whether AGREE strata carry accuracy, not whether they exist.
4. No retrieval/cosine/proxy numbers presented as score gains; H1 only, not A/B/C forecasts.

## Kill criterion (from mandate)
No stratum with ≥ 10% mass reaches weighted DESIGN sign accuracy ≥ 0.65 → confidence-weighted
amplitude is CLOSED; report. If one exists, report exact NMAE gain of the rescaling
counterfactual and flag for a 323 run (do NOT run it).
