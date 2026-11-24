# 384 — PLAN: do PromoterAI-derived per-gene priors beat θ_g / expression?

Hypothesis: per-gene regulatory/constraint features (promoter sensitivity,
constraint, promoter class, expression breadth) order per-gene sign accuracy on
held-out-line DE and add reach beyond 189's θ_g and control expression.
Falsifier: no feature beats control expression's held-out decile spread by
+0.02, or best feature prices Δreach < +0.01 raw (348 walk, paired bootstrap).

Inputs (read-only): 323/work {stats-cap15, reference(.npz+.parquet),
per-target-cap15.csv, metrics-cap15.json, effects.npz}; 189/work/prior189.npz
(θ_g); 348/reach348.py (imported walk conventions, not copied).
Downloads → data/ (all public, <100 MB total, licences noted in RESULTS):
gnomAD v2.1.1 LoF-by-gene (LOEUF/pLI), GTEx v8 median-TPM (tissue breadth),
UCSC CpG islands + GENCODE TSS (CpG-promoter flag). PromoterAI precomputed is
licence-gated by email (README) → fails the ≤30-min rule; documented, not used.

Estimator:
1. Known answers: 348 O1=0.254347/O4=0.280745 via imported reach348 walk;
   189 θ_g shape/tau check.
2. Per-gene H1 sign accuracy a_g = frac delivered-sign==ref-sign over 348's
   pool (ref q<0.05, minus own target); support floor reported.
3. Per feature: Spearman vs a_g, vs θ_g, vs log basal CPM; split-target
   held-out decile spread (order fixed, accuracy on held-out half) vs ≥200
   feature-permutation null; OLS a_g ~ θ_g + expr + feature (ΔR², t).
4. Best feature as reach ordering key: 348 walk_reach verbatim, Δ vs O1 and
   O4, 10k paired-target bootstrap CI.

Kill: max spread − expression spread < +0.02, or Δreach < +0.01 → NEGATIVE.
Controls: permutation nulls, shuffled-feature ordering, expression-only key,
θ_g-only key. Threads ≤4, RSS ≤6 GiB, one job at a time.
