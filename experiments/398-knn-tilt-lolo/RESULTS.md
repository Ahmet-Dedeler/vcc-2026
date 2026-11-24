# 398 — Neighbour-conditioned tilt, leave-one-bank-out: beats the marginal tilt on every held-out bank, +0.035 on Jurkat

Status: **COMPLETE** (PI-run, 2026-09-16 05:30–05:40Z; `lolo398.py`; `work/run.log` four banks, `work/run2.log` five banks
with HepG2 held out and voting). 371's LOLO framework, 370's block rule (top-2,200 by the held bank's own control among genes
shared with GWPS, the bank's own targets excluded), truth = the held bank's delta sign, unweighted (fid) and |delta|-weighted
(nmae proxy) accuracy, per-target paired bootstrap CI for kNN minus marginal.

| held-out bank (targets) | all-down | marginal tilt (vote banks) | kNN25 shrink50 | Δ vs marginal, unweighted [95 % CI] | Δ weighted |
|---|---|---|---|---|---|
| K562-essential (2,003) | 0.5115 | 0.5526 | 0.6136 | +0.0610 [+0.0578, +0.0642] | +0.088 |
| Nadig Jurkat (2,252) | 0.5073 | 0.5441 | 0.5795 | **+0.0354 [+0.0332, +0.0376]** | +0.051 |
| Nadig HepG2 (2,143) | 0.5296 | 0.5625 | 0.5799 | +0.0174 [+0.0157, +0.0192] | +0.022 |
| H1-2025 (237) | 0.5170 | 0.5259 | 0.5350 | +0.0091 [+0.0054, +0.0129] | +0.016 |
| RPE1 (2,187) | 0.5603 | 0.5666 | 0.5710 | +0.0044 [+0.0024, +0.0065] | +0.000 |
K = 10 / 25 / 50 / 100 within ±0.004 of each other on every bank; shrinkage toward the marginal never hurts and rescues RPE1
(pure kNN10 there −0.0003, shrunk +0.0072). Random neighbour sets of the same K score BELOW the marginal on every bank
(e.g. Jurkat random-25 −0.020), so the gain is the neighbourhood, not the vote count.

## Reading
The neighbourhood tilt is a strictly better block-sign predictor than the marginal tilt on every public line tested, with
the gain largest on lines close to K562 (Jurkat = context A) and smallest on the epithelial RPE1. Ships as 399
(`VCC_DE_SIGN_PRIOR=knn25sq80`, artifact `experiments/399-knn-tilt-arm/work/knn25s50_prior_challenge.npz`).

## 398c — the Orion (raw HCT116 / HEK293T, 17.7k + 18.0k perturbations) banks in the vote family (`lolo398c.py`, `work/orion/`)
| held-out | marginal, CLEAN4-family votes | marginal, + Orion votes | kNN25 shrink50, + Orion votes (Δ vs that marginal) |
|---|---|---|---|
| RPE1 (epithelial) | 0.5666 | **0.5970** | 0.5880 (−0.0090 [−0.0111, −0.0071]) |
| Nadig HepG2 | 0.5625 | 0.5641 | 0.5839 (+0.0198 [+0.0181, +0.0216]) |
| Nadig Jurkat | 0.5441 | 0.5417 | 0.5814 (+0.0397 [+0.0372, +0.0420]) |
Reading: the pooled Orion rows carry a strong **epithelial marginal tilt** (+0.030 on RPE1) that CLEAN4 lacks; on RPE1 the
25-neighbour vote (m = 50) is then WORSE than that marginal, i.e. the shrinkage weight should grow with the strength of the
marginal, while Jurkat keeps its +0.04 neighbourhood gain. Contexts B (HeLa) and C (CAL33) are epithelial carcinomas.
Tomorrow's one-variable candidate on 399: the Orion-inclusive vote family (with m raised to ~200) for B/C, CLEAN4 family for A.
Not part of 399 (registered before this ran).
