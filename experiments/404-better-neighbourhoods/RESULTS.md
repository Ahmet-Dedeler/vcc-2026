# 404 — A better neighbourhood for the block-sign tilt: no similarity or source beats Pearson-on-GWPS by more than noise

Status: **NEGATIVE** (COMPLETE; 2026-09-16; `run404base.py`, `run404sim.py`, `run404src.py`).

## Known-answer reproduction (numbers, match)
Copied `experiments/398-knn-tilt-lolo/lolo398d.py` verbatim (only the output path changed)
and re-ran: kNN25_shrink400, six-bank votes — Jurkat **0.57525**, HepG2 **0.58579**,
RPE1 **0.59881**; marginals 0.54171 / 0.56410 / 0.59702. Match to the recorded
`work/orion_m/results398.json` and the 399 PREREG table: **exact to all 4 decimals shown
(5th decimal identical in JSON), including every shrink arm and CI**. Proceeded.

## What was measured (inputs, estimator, controls)
- LOLO frame of 398d unchanged (370 block rule, truth = held-bank delta sign, votes =
  six-bank family minus recipient, shrink 400 unless varied, endpoints unweighted =
  fid-proxy and |delta|-weighted = nmae-proxy). Recipients: nadig_jurkat (2,252 targets),
  nadig_hepg2 (2,143), rpe1 (2,187).
- Similarity variants for GWPS neighbour search: (a) Pearson all 8,248 genes [baseline,
  bit-identical code path]; (b) Pearson on query's top-500 |delta| genes; (c) cosine on
  per-gene z-scored rows; (d) Spearman (average-tie ranks); (e) Pearson on top-50-PC
  projection (`svd_U/S`). K ∈ {10,25,50,100}, shrink ∈ {200,400,800}.
- Source variants: neighbours by GWPS rows vs by a recipient-like bank's rows (shared
  genes, GWPS-overlap perts: Jurkat←hepg2 / ←ess, HepG2←hct116 / ←hek293t,
  RPE1←hek293t / ←hct116) vs union of both top-25 lists (≈43–49 votes, judged against
  K-matched GWPS K=50). All source arms re-scored on the identical target subset.
- Statistics: paired per-target bootstrap, **10,000 resamples**, Δ vs baseline (a)K25
  with 95% CI. Controls: family marginal, all-down, random-25/50 (fresh seeds).

## Results (units: block-sign accuracy, unweighted; Δ vs shipped (a)K25/sh400)

Similarity @ K=25, shrink 400:

| arm | Jurkat Δ [95% CI] | HepG2 Δ [95% CI] | RPE1 Δ [95% CI] |
|---|---|---|---|
| (b) top-500 Pearson | −0.0008 [−0.0012,−0.0004] | +0.0005 [+0.0002,+0.0009] | +0.0019 |
| (c) z-scored cosine | **+0.0010** [+0.0007,+0.0014] | **+0.0013** [+0.0009,+0.0017] | +0.0014 |
| (d) Spearman | +0.0003 [+0.0001,+0.0006] | −0.0001 [−0.0004,+0.0001] | −0.0013 |
| (e) PC50 Pearson | −0.0019 | +0.0014 [+0.0008,+0.0019] | +0.0038 |
| random-25 (control) | −0.0549 | −0.0555 | −0.0579 |

Best K/shrink per line (Δ vs (a)K25/sh400, unw): Jurkat `c_kNN100` +0.0039 / `c_sh200`
+0.0049 (0.5801); HepG2 `e_kNN25` +0.0014 (nothing else > +0.0013); RPE1 `e_kNN10`
+0.0075 / `a_kNN10` +0.0058. Weighted endpoint mirrors unweighted throughout
(e.g. Jurkat c_sh200 +0.0069 w; HepG2 c_kNN25 +0.0021 w).

Neighbour source (subset Δ vs subset GWPS-Pearson K25, unw):

| pair (n tgts) | S-Pearson25 | S-z25 | S-top500 | union-P | union-Z | GWPS-K50 ref |
|---|---|---|---|---|---|---|
| Jurkat←hepg2 (2,072) | −0.0063 | −0.0085 | −0.0057 | +0.0016 | +0.0017 | **+0.0023** |
| HepG2←hct116 (1,935) | −0.0073 | −0.0064 | −0.0051 | +0.0006 | +0.0009 | −0.0010 |
| RPE1←hek293t (1,981) | −0.0138 | −0.0139 | −0.0121 | −0.0034 | −0.0034 | −0.0042 |
| Jurkat←ess (1,903) | +0.0007 [−0.0,+0.0015] | −0.0015 | +0.0001 | +0.0031 | +0.0027 | +0.0021 |
| HepG2←hek293t (1,943) | −0.0130 | −0.0124 | −0.0117 | −0.0007 | −0.0004 | −0.0010 |
| RPE1←hct116 (1,974) | −0.0071 | −0.0062 | −0.0052 | −0.0017 | −0.0013 | −0.0042 |

Median top-1 neighbour similarity (diagnostic): (a) 0.33–0.34, (b) 0.59–0.60,
(c) 0.29–0.30, (d) 0.28–0.29, (e) 0.82–0.83 — higher nominal similarity (b, e) does
NOT predict accuracy; per PLAN §falsifier these are flagged, not trusted.

## Verdict against the PLAN.md kill criterion
**Kill criterion triggered → NEGATIVE.** No variant clears +0.005 on both Jurkat and
HepG2 with CI above zero. Best joint similarity (c) is +0.0010/+0.0013; best union is
Jurkat-only (+0.0031, and only +0.001 over K-matched GWPS-K50). The neighbourhood
construction stays as shipped (Pearson-on-GWPS, K=25, shrink 400).

## Secondary readings for the PI (not score claims)
1. (c) z-scored cosine is the only similarity ≥ baseline on all three lines
   (+0.0010/+0.0013/+0.0014 @K25) — real but 4–5× below the gate; candidate only as a
   no-cost swap inside a future arm, never as its own A/B.
2. K/shrink optima are line-dependent: Jurkat prefers K≥50 + shrink 200 (c_sh200
   0.5801, +0.0049 — the largest single-line move seen, still short of the joint gate);
   RPE1 prefers K=10 + shrink 800; HepG2 sits at shipped K=25/sh400. Relevant to 402
   (Jurkat-only votes for context A): a Jurkat-tuned K50/sh200/c-sim variant is the one
   priced follow-up, but it must be validated on a Jurkat split-half (401 frame), not LOLO.
3. Surrogate-bank neighbour search is anti-informative (−0.005…−0.014) except same-type
   Jurkat←k562_essential (≈ neutral, +0.0007, CI touches zero). Cross-study rows are
   worse finders even when the study matches (Jurkat←hepg2 −0.0063).
4. Context check: the brief's "median top-1 r = 0.11" reproduces exactly on the
   challenge panel (269/300 targets in GWPS: median **0.1066**, p10 0.068, p90 0.238),
   while LOLO targets sit at 0.34. LOLO gains are measured on much cozier
   neighbourhoods than the panel will offer — discount transfer accordingly.

## Limits
- LOLO targets (≈2,000 bank-overlap perts) differ from panel targets (top-1 sim 0.34
  vs 0.11); nothing here prices the panel directly.
- (d) Spearman ties ~40–50% exact-zero delta entries at average rank — variant is
  valid but zero-heavy by construction.
- Bootstrap resamples targets, not genes; CIs ignore gene-level correlation.

## Exact commands to reproduce
```
python3 experiments/404-better-neighbourhoods/run404base.py   # known answer (~30 s)
python3 experiments/404-better-neighbourhoods/run404sim.py    # similarity × K × shrink
python3 experiments/404-better-neighbourhoods/run404src.py    # neighbour sources
```
Outputs: `work/{base,sim,src}/results{398,404sim,404src}.json`. No `src/`/`scripts/`
imports (numpy/scipy only); all banks mmap'd.

## Peak memory and wall time
Base 30 s wall; source sweep 188 s wall (measured); similarity sweep ~10 min wall
single-process (uninstrumented estimate). Dense peak bounded by 2× 9787×8248 float32
+ 10k×2252 int64 ≈ 0.8 GiB plus mmap; system free stayed ≥ 54% throughout, no scratch
files left behind.
