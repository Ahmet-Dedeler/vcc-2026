# 403 — Neighbourhood denoising buys cosine and sells retrieval on every line: DEAD for the bulk

Status: **NEGATIVE** (kill criterion fails on the retrieval leg 0/3; textbook 341
regime — cosine gain with retrieval loss — on all four recipients).

## Known-answer reproduction (own code, before any new number)

371's bank-QC table (`experiments/371-consensus-block-prior-lolo/RESULTS.md` §Bank QC),
own-target delta negativity among perturbations whose gene is on the bank axis:

| bank | reproduced cov | reproduced frac<0 | recorded (371) | match |
|---|---|---|---|---|
| k562_gwps | 7617/9787 | 0.99028 | 0.9903 | exact |
| k562_essential | 1835/2021 | 0.99673 | 0.9967 | exact |
| nadig_jurkat | 2045/2270 | 0.99707 | 0.9971 / mandate 0.99707 | exact |
| nadig_hepg2 | 1958/2160 | 0.99387 | 0.9939 | exact |
| rpe1 (not in 371) | 1928/2204 | 0.99741 | — | new QC row |

Shapes asserted in code: `D.shape == (len(perts), len(genes))` on all five banks.
Neighbour code path verified against `knn397.py:78-79` (centred full-GWPS-row
Pearson, self excluded) — identical construction.

## What was measured (inputs, estimator, controls)

- Source rows: K562-GWPS `delta.npy` (9787×8248, mmap). Recipients (truth = own
  delta rows): nadig_jurkat (2270×8881), nadig_hepg2 (2160×9623), rpe1
  (2204×8749), k562_essential (2021×8563, same-line upper bound).
- Per recipient: shared gene axis (recipient order); queries = targets present in
  both banks (≈2000 each); the query target's own gene column dropped from
  estimate, truth and every gallery row (exact per-query correction; 1800–1979
  queries had the gene on axis, remainder no-drop).
- Arms per query from GWPS rows only: (a) raw; (b) mean of K nearest GWPS rows
  (Pearson, self excluded), K∈{10,25,50}; (c) (1−w)·raw + w·meanK,
  w∈{0.25,0.5,0.75}; (d) mean of 25 RANDOM rows, 3 seeds; (e) shrink toward the
  GWPS grand mean, same w grid.
- Endpoints: (i) retrieval = mean 1−k/(N−1), pds tie rule (brief §1), gallery =
  ALL recipient rows; (ii) mean own cosine; (iii) pooled cos_agg = ΣE·T /
  sqrt(Σ|E|²Σ|T|²). Paired 10,000-target bootstrap (seeds 403/1403) for each arm
  minus raw on (i) and (ii).

## Results (units: fractions; CIs 95% paired-target bootstrap, arm minus raw)

Raw levels — the regime is healthy (retrieval far above 0.5, random control at
chance, §below):

| recipient | shared genes | queries | raw retrieval | raw owncos | raw cosagg | med top-1 / top-25 Pearson |
|---|---|---|---|---|---|---|
| nadig_jurkat | 7352 | 2252 | 0.8421 | 0.1607 | 0.2052 | 0.341 / 0.220 |
| nadig_hepg2 | 7424 | 2143 | 0.8023 | 0.1292 | 0.1659 | 0.341 / 0.220 |
| rpe1 | 7094 | 2187 | 0.7860 | 0.1328 | 0.1553 | 0.330 / 0.212 |
| k562_essential | 8202 | 2003 | 0.9025 | 0.2758 | 0.3722 | 0.321 / 0.212 |

Shrinkage grid (K=25; K=10/50 in JSON, same pattern — larger K loses more
retrieval): every cell retrieval-negative, every cosine-positive, all resolved:

| arm | Jurkat Δretr / Δcos | HepG2 Δretr / Δcos | RPE1 Δretr / Δcos | Ess Δretr / Δcos |
|---|---|---|---|---|
| mean25 | −0.137 [−0.146,−0.129] / +0.044 | −0.147 / +0.040 | −0.150 / +0.058 | −0.133 / +0.052 |
| shrink25 w.25 | −0.039 [−0.043,−0.036] / +0.026 | −0.044 / +0.022 | −0.043 / +0.029 | −0.042 / +0.038 |
| shrink25 w.50 | −0.072 [−0.078,−0.067] / +0.046 | −0.080 / +0.039 | −0.081 / +0.053 | −0.074 / +0.065 |
| shrink25 w.75 | −0.104 [−0.111,−0.097] / +0.053 | −0.114 / +0.047 | −0.115 / +0.065 | −0.104 / +0.071 |
| grand w.25 | −0.012 [−0.013,−0.011] / +0.002 | −0.016 / +0.003 | −0.014 / +0.006 | −0.016 / +0.002 |
| random25 (3 seeds) | retr 0.496–0.503 / cos −0.115 | retr ≈0.50 / −0.070 | retr ≈0.498 / −0.041 | retr ≈0.500 / −0.215 |

(Full 19-arm × 4-recipient table with CIs: `work/results403.json`; per-target
arrays: `work/per_target403.npz`.) Cos_agg rises 15–30% relative on every
neighbourhood arm (e.g. Jurkat shrink25 w.5: 0.205→0.265) — the aggregate-cosine
view alone would green-light this; retrieval vetoes it.

## Verdict against the PLAN.md kill criterion

**FAIL — record the negative.** Criterion: a shrinkage arm with retrieval CI
above zero on ≥2/3 cross-line recipients and own cosine not falling. Measured:
**0/3** — all 9 shrinkage arms lose retrieval on all 3 cross-line recipients
with CIs strictly below zero (mildest: −0.039 [−0.043,−0.036] at w=0.25);
cosine rises everywhere (+0.02…+0.07, all resolved). Controls behave: random
means sit at exactly chance retrieval (0.496–0.503, instrument healthy) and
below raw on cosine; grand-mean shrink also loses retrieval, so shrinkage per
se costs identity and the neighbourhood costs ~3× more per unit cosine. Same
pattern on same-line k562_essential (retrieval 0.903→0.860 at w=0.25). This is
the 341 regime the mandate names: **a cosine gain with a retrieval loss is NOT
a bulk candidate**. No K/w dose rescues it. Do not build.

## Limits

- Pure analysis on donor delta rows (pseudobulk-level responses), no whitening,
  no amplitude calibration, no emitted counts, no scorer — per mandate.
- Gallery is the full recipient bank (~2000–2270 rows), not a 300-way challenge
  panel: absolute retrieval levels are instrument readings, not board forecasts.
- Queries = targets present in both banks; the gene-drop applies where the
  target is on the shared axis (88% of queries), skipped otherwise.
- Neighbourhood from GWPS rows only (397/398 convention); a six-bank vote was
  not tested — it answers the sign task, and porting it to row-denoising needs
  design choices outside this mandate.
- Neighbour-strength QC differs by query population: median top-1 Pearson
  0.32–0.34 here vs 397's 0.11 on its 82 H1 targets (same code path — target-set
  property, not a bug; GWPS pert names are all unique, no sibling rows).
- Memory incident, disclosed: the first two runs peaked 10.3–12.4 GiB RSS
  (above the 8 GiB cap) — a `(nT, K, Gs)` fancy-index gather materialised up to
  3.3 GiB transient per K. Rewrote as per-query accumulation; definitive run
  peaks 3.82 GiB. Machine stayed ≥54% free throughout. Result JSON md5-identical
  across the pre- and post-fix runs (`a796cd3a…`), so no number below depends on
  the fix.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 40%
NUMBA_CACHE_DIR=experiments/403-bulk-knn-denoise/work/numba-cache \
  nice -n 10 uv run python experiments/403-bulk-knn-denoise/run403.py
```

Reads only: the 5 donor banks above. Writes only `work/`: `results403.json`,
`per_target403.npz`, `run.log`. Deterministic (seeds 403/1403, stable argsort).

## Peak memory and wall time

`/usr/bin/time -l`, definitive run: **peak RSS 3.82 GiB** (cap 8 GiB), **wall
89 s** (jurkat 24 s / hepg2 21 s / rpe1 22 s / essential 22 s), single process,
OMP/MKL 4 threads, `nice -n 10`. Well inside the 2 h budget.
