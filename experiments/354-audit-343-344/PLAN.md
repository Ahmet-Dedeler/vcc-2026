# 354 — Independent audit of 343 (NMAE decomposition) and 344 (null-key killcheck)

## Hypothesis
The numeric claims in `experiments/343-nmae-decomposition/RESULTS.md` and
`experiments/344-null-key-killcheck/RESULTS.md` are re-derivable from the stated inputs
with independently written code. Each claim is scored MATCH / MISMATCH with both values.

## Inputs (read-only, never modified)
- 343: `experiments/323-current-family-validation/work/{stats-cap15.npz, reference.npz,
  reference_de.parquet, effects.npz, metrics-cap15.json, per-target-cap15.csv}`.
  `reference.npz` is 3.0 GiB (contains 10779×18400 float64 control matrices I do not need);
  load with `mmap_mode='r'` and touch only `genes/targets/lfc/q/keep/ref_lin`.
- 344: `data/challenge/context_A.h5ad` (215 MB) via h5py, manual CSR→CSC, ≤ 6 GiB.

## Estimator (343, all from scratch with numpy/pandas only)
- Gate per target: universe genes with `q<0.05`, finite nonzero `lref`, gene symbol ≠ target.
  Eligible iff ≥ 10 genes. NMAE_t = Σ|l̂−l|/Σ|l|; mean over eligible targets.
  Reference = `reference.npz` lfc/q (verify bit-equal to parquet on one target first).
  Ours = `stats-cap15.npz` lfc/q. Design = `effects.npz::cap15` mapped by gene name
  (genes absent from the effects axis → design 0).
- Magnitude ratio, |lfc|-weighted sign accuracy: computed POOLED and as mean-of-per-target;
  the variant matching the claim decides the reported convention (state both).
- Counterfactuals per target then averaged: (a) sign(lref)·|l̂|; (b) least-squares scalar
  a=⟨l̂,l⟩/⟨l̂,l̂⟩; (c) grid search gain over {0.25,…,6.0} step 0.25 (try step 0.1/0.5 if mismatch).
- Strata: quintiles of `ref_lin`; called = our q<0.05; UP = lref>0; mass = Σ|lref| shares.
- Claim 6 bins: quintiles of `ref_lin` over universe genes (pre-check: 20th pct = 15.73,
  min = 5.00 → bottom quintile should read 5.0–15.7, an independent edge check).

## Estimator (344, context A only)
- Control CPM per gene = mean over 18,400 cells of counts/(lib/1e6); universe = genes > 5 CPM.
- Hand-rolled Wilcoxon rank-sum z (average ranks, tie correction) on log1p(CPM),
  per construct (400 cells) vs the pooled 18,400 control cells (self included, literal reading
  of "against the pooled 18,400 controls"; sensitivity run vs the other 18,000 if time permits).
- BH step-up α=0.05 per construct (hand-rolled); count constructs with any hit; median hits.
- Key = mean |z| over 46 constructs per universe gene; hand-rolled Spearman vs control CPM.

## Controls / falsifiers
- Known-answer gate: reproduce 323's recorded raw NMAE 0.9508364630325662 and the 58-row
  per-target mean before trusting any downstream number; reproduce parquet==npz equality.
- 344 has no recorded local artifact to match; the cross-check is internal consistency
  (universe size → hit counts → key correlation each depend on the previous step).
- This audit FAILS (reports MISMATCH, does not "fix") if any claim differs beyond stated
  tolerance after exhausting documented convention variants (pooled vs mean, grid step).

## Kill criterion
A clean MATCH on claims 1–3 (the load-bearing NMAE/sign numbers) plus claim-4 pool size
validates the audit's own plumbing; any MISMATCH there stops dependent strata claims and is
reported with the likely cause (convention vs arithmetic). Small Wilcoxon numerics on 344
are acceptable per the mandate and reported as deltas, not matches.

## Budget
~1.5 h. 343 script seconds; 344 a few minutes single-threaded. Peak target < 6 GiB RSS.
