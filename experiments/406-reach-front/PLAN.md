# 406 — Where does our reach walk break, and can a neighbourhood-aware ORDER key move the front?

## Hypothesis
reach walks the reference's significant genes in our confidence order and pays the deepest
≥90%-pure prefix. On H1 the walk breaks at k* ≈ 14 on average (388). Two order levers are
closed: design–marginal-prior agreement boost (390 §8, monotone loss) and promoting confident
block genes to the front (397b, null). Untested: the ORDER WITHIN THE HEAD (the first 400
genes by |design|), re-ranked by a per-target neighbourhood confidence. 397b showed that on
head genes with |P−0.5| ≥ 0.2, P(right|agree) ≈ 0.65 vs P(right|disagree) ≈ 0.46 — so down-weighting
or demoting disagreeing head genes at the front could deepen k* without touching any sign
(fid/nmae/jac unchanged by construction).

## Exact inputs (paths)
- Frame + reach walk VERBATIM from `experiments/397-knn-tilt/knn397.py` (everything above the
  `# ---- kNN tilts` marker, exec'd, not copied): `cur` (delivered order), `walk_reach`,
  `head` (top-400 by |E1| = cap15 design), `block`, `cb`, `sig/adjud/pool`, `n_conf`,
  `SG0` (delivered signs), `M0`, `truth=sign(lrm)`, `E0/E1`, `GIDX/LREF/elig`.
- Shipped-sign construction VERBATIM from `experiments/397-knn-tilt/knn397d.py` (six-bank
  family, K=25, shrink 400, family-marginal default; signs on called block only): banks
  `data/donors/{k562_gwps,k562_essential,nadig_jurkat,orion_hct116,orion_hek293t}` +
  `experiments/166-donor-acquisition/donors/nadig_hepg2`.
- Key definition: |z| := |E1| (cap15 design magnitude — the quantity that defines the head;
  the 323 analytic frame carries no z-statistic, so this substitution is documented, not hidden).
  conf := |Ps−0.5| with the SHIPPED shrunk Ps (six-bank, shrink 400); NaN (no vote anywhere) → 0.
  agree := (shipped kNN sign == delivered sign), defined on every gene via the fallback chain.

## Estimator
- Signs fixed at SG399 everywhere (descriptive + all variants) → fid/nmae/jac identical by
  construction; only `walk_reach` is re-run with re-permuted `cur`.
- Within-head variants re-permute ONLY head genes among the slots they already hold in the
  delivered order (stable, delivered-rank tiebreak): (a) delivered order; (b) key
  |E1|·(1+β·agree), β∈{0.5,1,2}; (c) key |E1|·(1+β·conf·agree−β·conf·disagree), same β;
  (d) demote head genes with (disagree & conf≥0.2) to end of head, stable; (e) deranged
  control: per-target permutation of conf within head genes, re-run (cβ) and (d); plus (b,β=1)
  with deranged agree.
- Endpoint: reach per variant; Δ vs (a) with paired 10,000-target bootstrap CI (seed 406,
  resample the 82 scored targets); break-rank distribution per variant.
- Descriptive (under SG399 + delivered order): per-target break rank b = k*+1 (None if the walk
  never dips below 0.9), breaking gene's head/block/uncalled class, |E1|, conf, agree,
  correctness; P(10)/P(20)/P(50); table by n_conf tercile.

## Controls / falsifiers
1. Known answers before any new number: delivered reach 0.25435 (85-target mean, 82 scored)
   and E0 block precision 0.5083 (390). Mismatch → stop.
2. Construction cross-checks (not gates, vintages differ): SG399 block precision ≈ 0.5338 and
   reach(SG399, delivered) ≈ 0.2519 (399 prereg addendum: 0.2544−0.0025).
3. Chunked Orion voter must reproduce verbatim `votes()` EXACTLY (int64 counts) on all small
   banks fully and on a 200-row Orion subset before use.
4. Deranged control (e) must read near zero; a variant that wins only with deranged conf too is
   an order-statistic artefact, not signal.
5. Memory: Orion banks stay mmap'd, streamed in ≤500-row chunks; cap 8 GiB RSS; no writes
   outside `experiments/406-reach-front/`.

## Kill criterion (from assignment)
No variant clears +0.005 reach with CI above zero and the deranged control near zero → the
order stays; the descriptive table is the deliverable regardless.
