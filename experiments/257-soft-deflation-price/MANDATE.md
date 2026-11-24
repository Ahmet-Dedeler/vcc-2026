# 257 — Price a SOFTER CD4 deflation on top of `cd4svd5`, on the transfer factors the board just measured. Measurement only. DO NOT BUILD, DO NOT SUBMIT.

You are a research subagent. Work in `experiments/257-soft-deflation-price/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. PLOLO fold loads are the heavy step: one
fold at a time, peak RSS under ~10 GB, never alongside another ≥ 15 GiB job (check `ps` for builders / `vcc prep`;
none should be running). Record peak RSS per script.

## Why (2026-09-05 00:55Z board point, `experiments/249-cd4-deflate/BOARD-RESULT.md`)
`cd4svd5` (`VCC_CD4_DEFLATE_K=5`: the CD4 residual centred over its carried rows and its top-5 left-singular
directions REMOVED before sizing, 246's `lib246.deflate(mode="svd", k=5)` resolved variant) is the incumbent at
**0.2160944**, +0.0042634 over `h1union`. The board paid **`pds` +0.00545 overall (raw +0.01492 = 1.22× PLOLO's
+0.012263)** and charged **`mse` −0.00080 (raw +0.00454: board cosine −1.735 % against PLOLO's aggregate −0.4617 %,
i.e. a cosine transfer factor of 3.76 for this "re-shape occupied rows" class)**; `nmae`/`fid`/`reach`/`jac` were flat
as priced. So the arm netted +0.0043 instead of its +0.0086 upside because the projector spends cosine at 3.8× the
local read. **Question: does a SOFTER removal — shrink the top-k component by s < 1 instead of removing it, or a
different k — buy more `pds` per unit of board cosine, on top of what `cd4svd5` already ships?** 255 already showed
k = 10 (full projector) is dead on these factors (`pds` +0.00054 on 1/3 folds, aggregate −1.08 % more cosine).

## READ FIRST
`docs/READ-FIRST.md` (v2.3), `experiments/249-cd4-deflate/BOARD-RESULT.md` (§2 for the factors), `experiments/237-gate-validation/BOARD-RESULT.md`,
`experiments/243-transfer-factor/RESULTS.md` §4.4 and §5.2 (the class rule; read cosine MEASURED, never the floor),
`experiments/246-direction-hypotheses/lib246.py` + `sweep246b.py` + `price246.py` (the deflation code and 246's
price recipe), `experiments/249-cd4-deflate/RESULTS.md` §2 and §9, `experiments/255-cd4-deflate-k10/RESULTS.md` §2.1,
`experiments/162-retrodicting-validator/` (PLOLO, the LOO σ 0.001494 raw), `AGENTS.md` §8.1b (gates), §11.1 rule 10.

## The instrument and the price recipe (fix these before measuring; put them in a `price257.py`)
On 162's PLOLO, three cancer folds, the `resolved` construction (the deflated residual re-sized by `solve_cd4_fraction`,
exactly as the builder does it — 249 §2 / 254 KA 3 explain why the resolved variant is the package's object):
- `pds` leg, overall = 2.2138 × (k_pds × Δ`pds`_3cancer vs `cd4svd5`) / 6, with **k_pds = 1.22** (point) and 0.85
  (conservative, 0.7× of it).
- `mse` leg, overall = −1.051 × Δraw / 6, with Δraw = −1.996 × Δcos_board (197's identity slope at these amplitude
  terms, 249 §9) and Δcos_board = **k_cos × rel_cos_agg(PLOLO, vs cd4svd5) × 0.10547** (cd4svd5's measured board
  cosine: h1union's 0.107332 × (1 − 0.01735)); **k_cos = 3.76** for a re-shape of occupied rows (point; also report at
  2.15, the macro-based read, and at 1.36, the old re-dose value, so the sensitivity is visible). If a candidate's
  delivered amplitude moves > 0.2 % (report it via 246's delivered ratio), add 243's amplitude leg and say so.
- DE-three (`fid`, `reach`, `jac`): 0.24 × Δcos_board × (their 243 §4.3 propagation), and `nmae` by the 216/188 exact
  identity (229's `six229.py` computes it; it was VALIDATED at scale tonight, BR237 — treat its magnitude as a
  lower bound on the loss, ×1.0 point, ×1.5 conservative).
- **Retrodiction rows, must pass before any candidate is priced:** (R1) `cd4svd5` vs `h1union` through this recipe
  reads **+0.0043 ± 0.0006 overall, `pds` +0.0055 ± 0.0006, `mse` −0.0008 ± 0.0003** (the board point);
  (R2) `svd k=5 s=1` reproduces 246's `B2_cd4_svd5_resolved` Δ`pds` **+0.012263** and aggregate **−0.4617 %** to the
  digit (its `price246.json` row); (R3) `svd k=10 s=1` vs `cd4svd5` prices negative on the point factors (255 §2.1's
  increment: `pds` +0.00054, aggregate −1.08 %).

## The screen
Add a `shrink` variant to 246's deflation (in YOUR copy under this directory, `lib257.py`; do not edit `lib246.py`):
`R_c − s·U_k U_kᵀ R_c` for the top-k centred directions, support re-zeroed, then resolved sizing. Grid:
**s ∈ {0.25, 0.5, 0.75, 1.0} × k ∈ {3, 5, 8}**, plus **per-direction shrinkage** where each of the top-5 directions is
shrunk alone at s = 1 (five cells: which direction carries the `pds`, which carries the cosine cost). For every cell,
on the three cancer folds: Δ`pds` vs `cd4svd5` (per fold, 3/3 sign count, LOO σ), rel `cos_agg` and macro vs
`cd4svd5`, energy removed, delivered-amplitude ratio, then the six-member price at point and conservative factors.
Controls: the energy-matched random subspace at the best cell's energy (246's `rand`, one seed), and 246's
leave-top-3-out retention for the best cell.

## Deliverable
`RESULTS.md`: the three retrodiction rows first (abort and report if R1 or R2 fail); the full grid as one table
(increment over `cd4svd5`, point and conservative price); the best cell; **whether ANY cell clears +0.0019 overall
over `cd4svd5` at the conservative factors** (the shape-only gate; a softer removal is shape-only if amplitude stays
within 0.5 %); the builder flag such a cell would need (the builder has only `VCC_CD4_DEFLATE_K`, a full projector;
name the flag and the one-line change, do not make it); and one line **CLEARS / DOES NOT CLEAR**. If it clears, add
`PREREGISTRATION-DRAFT.md` in 249's form (six-member prediction, five falsifiers). A well-supported negative is a
valid result. No model name, no description, no submission. Then stop.
