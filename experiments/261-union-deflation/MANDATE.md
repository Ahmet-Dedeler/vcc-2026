# 261 — Deflate the H1-UNION residual's target-shared directions before sizing, on top of `cd4svd5`: the CD4 mechanism's direct analogue, priced on the factors that mechanism just measured. Measurement only, ≤ 3 hours. DO NOT BUILD, DO NOT SUBMIT.

You are a research subagent. Work in `experiments/261-union-deflation/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold at a time, peak RSS under ~8 GB
(a sibling, 260, may still be finishing at < 1 GB). **Time box: write up after ~3 hours of wall clock even if the grid
is incomplete.**

## Why (the one hypothesis left with a mechanism AND a board-validated analogue)
`cd4svd5` — the CD4 residual centred over its carried rows with its top-5 target-shared left-singular directions removed
before sizing (246/249) — is the incumbent at 0.2160944: **+0.00545 overall on `pds` (raw +0.01492, 1.22× PLOLO)** for
**−0.00080 on `mse`** (cosine cost 3.76× PLOLO's aggregate read), `experiments/249-cd4-deflate/BOARD-RESULT.md`. The
mechanism (P249 §2): a residual that enters UN-whitened keeps its target-shared hub — the component that makes
perturbations look alike — which the whitener strips from the donors and nothing strips from a residual. **The H1
union (`gse295214_h1` + `h1_2025`, `VCC_H1_UNION_G=0.28`, ~84 carrying rows) enters the same way and has never been
deflated.** 257 showed the CD4 projector's `pds` is direction-specific (direction 1 carries 75 %) and 259 showed the
union's directions buy `pds` on top of `cd4svd5` when their magnitude cost is removed. Question: **does removing the
union residual's top-k target-shared directions before its carrying-row sizing buy `pds` at a cosine cost that prices
positive on the re-shape factors, on top of `cd4svd5`?** One variable. A negative closes the analogue.

## READ FIRST
`docs/READ-FIRST.md` (v2.8; §3 pricing rule, §4 instruments), BR249 §2–§3, `experiments/246-direction-hypotheses/lib246.py`
(`deflate` — apply it to the UNION residual object, centred over ITS carrying rows, support re-zeroed; `mode="svd"`,
`"rand"`, the shuffled-projector and split-half controls as 246 §4 ran them for CD4), `experiments/257-soft-deflation-price/`
(`lib257.py`, `price257.py`, `sweep257.py`: the resolved construction, the price recipe, the fold caches — byte-copy
`work/plolo/`), `experiments/249-cd4-deflate/RESULTS.md` §2–§3 (how the deflated residual is re-sized: CD4 by
`solve_cd4_fraction`; the union by the carrying-row rule `solve_carrying_row_fraction` at g — deflate BEFORE sizing,
exactly as 249 did for CD4), `experiments/259-union-without-igvf/price259.py` (the union path through
`apply_frozen_cd4_residual` on the `cd4svd5` base), `AGENTS.md` §8.1b, §11.1 rule 10.

## Construction (three cancer folds, `kolf2` apart, per-fold rows, LOO σ)
On the `cd4svd5` base (CD4 deflated k = 5, f 0.34, the shipped stack), the H1-pair union residual `U` on its carrying rows:
`U_c = U − mean_rows(U)`, remove the top-k left-singular directions of `U_c` (k ∈ {1, 2, 3, 5, 8}), support re-zeroed,
then the carrying-row sizing at g = 0.28 (the RESOLVED variant — the sizing rule sees the deflated object, as the
builder would). Also k = 3 and 5 at g ∈ {0.20, 0.35} (the dose may re-optimise once the hub is gone). Controls at the
best k: energy-matched random subspace (246 `rand`, one seed), shuffled projector, split-half projector (fit on one
half of the union's rows, applied to the other), leave-top-3-out `pds` retention, and the per-target share of the top 3.
**Known-answer rows first:** (R1) k = 0 reproduces `cd4svd5`'s rows from 257/258/259 at 0.0e+00 on every fold and on the
production frame; (R2) the price recipe retrodicts BR249 (+0.0043 ± 0.0008; `pds` +0.0055, `mse` −0.0008) and 257's
best cell; (R3) the same `deflate` code applied to CD4 at k = 5 reproduces 246's `B2_cd4_svd5_resolved` (+0.012263,
−0.4617 %) — the code path is the proven one.

## Price — COHERENT PESSIMISTIC / POINT / OPTIMISTIC (READ-FIRST §3)
Legs as 257 §1.3 with the bands: `pds` k ∈ {0.85, 1.22} (re-shape class; pick the end that lowers the total for the
cell's sign), cosine k_cos ∈ {1.36, 3.76} (the union's rows are occupied rows; 3.76 is the measured re-shape value —
pessimistic for a cosine loss), DE-three rate ∈ {0.0425, 0.24}, `nmae` identity ×1.5 on a loss / ×0 credit on a gain
(pessimistic) and ×1.0 (point), amplitude leg charged when the delivered ratio moves > 0.2 %. **CLEARED = pessimistic
≥ +0.0019 over `cd4svd5` AND leave-top-3-out ≥ 50 % AND `pds` > 0 on 3/3 folds at ≥ 2 LOO-σ AND the random and
shuffled controls below half the arm.**

## Deliverable
`RESULTS.md`: R1–R3 first; the k × g table with all three scenarios; the controls; the best cell in full; the builder
flag it would need (analogous to `VCC_CD4_DEFLATE_K`: name it `VCC_UNION_DEFLATE_K`, describe the one-place change
next to `deflate_cd4_residual`, do not make it); one line **CLEARS / DOES NOT CLEAR**; if it clears,
`PREREGISTRATION-DRAFT.md` in 249's form with five falsifiers priced on the three scenarios. A clean negative closes
the analogue and is a valid result. No model name, no description, no submission. Then stop.
