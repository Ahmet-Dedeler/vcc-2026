# 403 — PLAN: does neighbourhood denoising of donor rows help the BULK channel?

## Hypothesis
A K562-GWPS donor row `r_t` is mostly sampling noise per target (F1: split-half
residual cosine ≈ 0.05 at ~80 cells/half). 397/398 showed the SIGN of the 25
nearest GWPS perturbations generalises cross-line. If the shared (transferable)
component dominates the neighbourhood mean, then shrinking `r_t` toward its
neighbourhood mean should raise both own-cosine to a held line's truth and
retrieval rank — i.e. the same idea would pay pds/mse through the bulk. If the
neighbourhood mean is mostly shared noise / generic programme, shrinkage blurs
identity: own cosine may rise while retrieval falls (the 341 regime) — which is
NOT a bulk candidate.

## Inputs (exact paths)
- Source: `data/donors/k562_gwps/{delta.npy,perts.txt,genes.txt}` (9787×8248).
- Recipients (truth = their own delta rows):
  `data/donors/nadig_jurkat`, `data/donors/k562_essential`, `data/donors/rpe1`,
  `experiments/166-donor-acquisition/donors/nadig_hepg2` (same 5-file layout).
- Reference numbers: `experiments/371-consensus-block-prior-lolo/RESULTS.md`
  (own-target negativity table), `experiments/398-knn-tilt-lolo/RESULTS.md`
  (Jurkat marginal 0.5441 — context only, different endpoint, not reproduced).
- Neighbourhood convention copied from `experiments/398-knn-tilt-lolo/lolo398.py`:
  Pearson correlation of FULL GWPS rows (8248 genes), self excluded.
- pds form from `docs/reference/vcc2026-metrics-BRIEF.md` §1:
  PDS_p = 1 − k_p/(n−1), k_p = #{d<d_self} + ½(#{d==d_self}−1); d=1 on zero norm.

## Estimator
Per recipient R: shared gene axis (intersection, order of R); queries = targets
present in BOTH banks. Vectors = delta rows on shared genes, with the QUERY
target's own gene column dropped (zeroed) from estimate, truth and every gallery
row (per-query correction, exact, no approximation).
- (a) raw: GWPS row r_t.
- (b) meanK: mean of K nearest GWPS rows (Pearson, self excluded), K∈{10,25,50}.
- (c) shrink: (1−w)·r_t + w·meanK, w∈{0.25,0.5,0.75} (9 arms).
- (d) control: mean of 25 RANDOM GWPS rows (self excluded), 3 seeds.
- (e) control: shrink toward the GWPS grand mean (all 9787 rows), same w grid.
Endpoints per recipient, per arm:
  (i) retrieval = mean over queries of 1 − k/(N−1), gallery = ALL R rows
      (N = full bank size), cosine distance, pds tie rule;
  (ii) own cosine = mean over queries of cos(estimate, truth);
  (iii) cos_agg = pooled cosine ΣE·T/sqrt(Σ|E|²Σ|T|²) (207/159 convention).
Uncertainty: paired 10,000-target bootstrap (seed 403) of each arm minus (a),
for (i) and (ii). Also recorded: shared-axis size, n queries, median top-1/top-25
neighbour Pearson (neighbourhood-strength QC).

## Controls / falsifiers
- (d) must sit at/below (a) on retrieval (398: random-25 −0.020 on the sign task).
  If random means beat raw, the instrument is broken (norm/centering artefact).
- (e) separates "shrinkage toward anything helps" from "the neighbourhood helps".
- k562_essential (same line, other study) is the upper-bound recipient; the three
  cross-line recipients carry the verdict.
- 341-regime tripwire: Δcos>0 with Δretrieval<0 (either CI-resolved) = NOT bulk.

## Kill criterion (from the assignment)
If NO shrinkage arm (c) improves retrieval on ≥2 of the 3 cross-line recipients
with bootstrap CI strictly above zero WHILE own cosine does not fall (Δcos CI
not resolved-negative, i.e. lower bound ≥ −0.002), the route is DEAD for the
bulk — report NEGATIVE. A cosine gain with a retrieval loss is the 341 regime,
not a candidate, even if retrieval is flat.

## Cost / safety
mmap source, float32, chunked matmuls; est. peak < 3 GiB (cap 8 GiB RSS);
single process, NUMBA untouched (no src/ imports), OMP/MKL ≤ 4 threads.
Known-answer job first (seconds); main run est. < 10 min. One heavy job at a
time; production build (399) has priority — nice +n10, no BLAS oversubscription.
