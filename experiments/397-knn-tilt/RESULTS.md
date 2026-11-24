# 397 — Neighbour-conditioned block tilt on H1: block precision 0.530 vs CLEAN4 0.522 (analytic, scorer truth)

Status: **COMPLETE** (PI-run, 2026-09-16 05:27Z, 9 s; `knn397.py`, `work/results397.json`). Frame copied from 390's
`consensus390.py`; known answers reproduced exactly (delivered block precision 0.5083, CLEAN4 0.5223 = 390 §5).

## Construction
For each of the 82/85 H1 targets present in K562-GWPS: neighbours = top-K GWPS perturbations by Pearson correlation of the
GWPS response rows (self excluded; median top-1 correlation 0.11 — the rows are noisy); P_down per gene = fraction of down
votes over the neighbours' rows (GWPS alone, or GWPS + K562-essential + Nadig Jurkat rows of the same perturbations);
optional shrinkage toward the marginal with prior weight 50; signs on the called CPM block only, |E0| magnitudes (392's
construction); head untouched. Controls: random neighbour sets of the same K (3 draws), the marginal tilts, CLEAN4.

## Results (raw Δ vs delivered E0; paired 4,000-target bootstrap; H1-local overall estimate at the board slopes)
| arm | block precision | Δfid | Δnmae (lower better) | Δreach | overall est. |
|---|---|---|---|---|---|
| CLEAN4 marginal (= 392) | 0.5223 | +0.0116 [+0.0060, +0.0177] | +0.0017 | +0.018 | +0.0095 |
| GWPS+ess+Jurkat marginal | 0.5225 | +0.0118 | +0.0010 | −0.014 | +0.0038 |
| kNN10 gw+ess+jk shrink50 | **0.5313** | **+0.0184 [+0.0120, +0.0251]** | −0.0021 | −0.004 | +0.0101 |
| kNN25 gw+ess+jk shrink50 | 0.5300 | +0.0174 [+0.0113, +0.0245] | −0.0008 | −0.010 | +0.0081 |
| kNN50 gw+ess+jk | 0.5293 | +0.0169 [+0.0105, +0.0236] | −0.0018 | +0.002 | +0.0103 |
| kNN200 gw+ess+jk shrink50 | 0.5286 | +0.0165 | −0.0014 | +0.008 | +0.0110 |
| random-25 GWPS (3 draws) | 0.508–0.510 | +0.0003…+0.0016 | +0.001…+0.005 | +0.008…+0.035 | — |
| random-500 GWPS (3 draws) | 0.5175 | +0.0078 | +0.0013…+0.0020 | +0.004…+0.010 | — |
Random neighbour sets converge to the GWPS marginal (0.5175) as K grows and never reach the kNN values: the gain is the
neighbourhood. nmae moves from a small harm (CLEAN4 +0.0017) to a small gain (−0.002). Reach is inside its noise everywhere.

## Reading
A target-specific, denoised tilt beats the target-agnostic tilt on H1 by +0.007…+0.009 block precision (≈ +0.006 raw fid over
392's construction). H1 is the line where every prior has been weakest; 398 prices the same construction leave-one-bank-out.

## 397b — what else the neighbourhood confidence buys (analytic, `knn397b.py`, 2026-09-16 05:45Z)
- **Block precision is graded by confidence |P − 0.5|:** [0, .05) 0.505 (delivered 0.509); [.05, .1) 0.524; [.1, .2) 0.551;
  [.2, .3) 0.587; [.3, .5) 0.605 (delivered 0.51 throughout). Only 10–18 % of block calls sit in the reference gate.
- **Head genes:** where the neighbourhood is confident (|P − 0.5| ≥ 0.2, 9 % of head calls) the delivered sign is right 0.65–0.68
  when it agrees with the kNN sign and only 0.46 when it disagrees — a real signal, but 784 disagreeing calls on 85 targets
  (≈ +0.2 % head precision if flipped). Not an arm today.
- **Reach order:** promoting confident block genes to the front of the order does nothing (τ 0.25: −0.0016 [−0.0088, +0.0044];
  higher τ touches ≤ 1 gene per target). With 399's signs and the delivered order H1 reach reads 0.2441 vs 0.2544 with E0 signs
  (inside H1's ±0.04 reach noise; the board's 392 also read reach −0.004 where H1 said +0.018). Closed.
- **Confidence-targeted magnitudes:** ×2 on block genes with |P − 0.5| ≥ 0.2 is flat (+0.0003), ×4 harms (+0.0075 nmae),
  ×8 harms more; the gate's block genes carry small true |lfc|. Magnitudes stay dead. Closed.

## 397c — neighbourhood signs on UNCALLED genes (`knn397c.py`; nmae and reach read them, fid does not)
Uncalled non-head genes by kNN confidence: [.1, .2) n=108k kNN precision 0.545 vs delivered 0.505 (on the reference gate 0.642
vs 0.530); [.2, .3) 0.571 vs 0.507 (gate 0.726 vs 0.549); [.3, .5) 0.575 vs 0.509 (gate 0.810 vs 0.560). Replacing the
delivered sign there (vs 399's signs, |E0| magnitudes): conf ≥ 0.2 → **nmae −0.0015 [−0.0024, −0.0007] resolved**, reach 0;
conf ≥ 0.1 → nmae −0.0025 [−0.0057, +0.0012], reach +0.0015 [0, +0.0046]. 390 §7 closed the MARGINAL prior on uncalled genes
(reach −0.022); the neighbourhood version is positive. Needs an emitter extension (the prior mode steers the called block only,
and uncalled genes have no class capacity reserved) — a later variant, not part of 399.

## 397e — head flips where the shipped neighbourhood (six-bank, shrink 400) disagrees (`knn397e.py`)
With the heavy shrink the neighbourhood rarely disagrees with a delivered head sign at confidence ≥ 0.2 (245 of 30k head calls,
0.8 %), and the precision swing there is noise-sized (0.465 → 0.535 on 245 calls; other thresholds go the other way). Closed:
no head arm from this construction.
