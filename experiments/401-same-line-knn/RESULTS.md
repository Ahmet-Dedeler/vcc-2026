# 401 — Same-line votes: on Jurkat, the neighbourhood tilt with Jurkat-only votes reaches 0.607 vs 0.581 for the six-bank family (split-half)

Status: **COMPLETE** (PI-run, 2026-09-16 05:52Z, `split401.py`, `work/results401.json`; three seeded target splits, sd ≤ 0.001).
Nadig Jurkat targets split in halves; eval-half targets take their 25 nearest GWPS perturbations; votes from (a) the six-bank
family (399's construction, shrink 400), (b) the FIT-half Jurkat rows only (same line, disjoint targets), (c) both; truth = the
eval target's Jurkat delta sign on Jurkat's top-2,200 block (own targets excluded). Units: unweighted block accuracy.

| prior | accuracy |
|---|---|
| six-bank marginal (399's default) | 0.5417 |
| six-bank kNN25, shrink 400 (= 399) | 0.5806 |
| Jurkat-only marginal (370's same-line tilt) | 0.5807 |
| **Jurkat-only kNN25, shrink 50** | **0.6065** |
| Jurkat-only kNN25, shrink 400 | 0.5855 |
| six banks + Jurkat, kNN25, shrink 400 | 0.5882 |

## Reading
Same-line rows carry the same neighbourhood signal on top of a much stronger marginal: +0.026 over 399's construction on Jurkat
truth. Context A is Jurkat; B/C have no same-line bank. Candidate arm (402): per-context per-target prior — A from Jurkat-only
kNN25 shrink 50, B/C unchanged from 399. Needs the builder to read `psign_by_target_A` when present (small patch; apply after the
399 build releases the builder). Expected: +0.026 × 0.85 (block share of the call set) / 3 contexts ≈ +0.007 raw fid ≈ +0.004
overall, if challenge A behaves like Nadig Jurkat (319pA's bulk-channel failure was a pds cost, which this channel does not pay).
