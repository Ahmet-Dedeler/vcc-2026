You are starting experiment 183 in ~/Code/virtual-cell-challenge-2026. Work only in `experiments/183-unique-boost/`. **No builds, no `vcc`, no `cell-eval2`, no submission.** Peak RSS ≤ 12 GB. `uv run python`.

READ FIRST: `AGENTS.md` §4.2 (whitening / hubness), §11.1 rules 1–13 and the 181 amendments; `docs/CURRENT_STATE.md` 2026-09-03 00:00 UTC `abc-knee` result AND the 181 audit REOPEN item; `experiments/156-delivered-geometry/RESULTS.md` (centroid / rank-k / row-gram already screened in delivered space — do not redo them); `experiments/170-gamma-mse/RESULTS.md` (pds peaks at gamma 0.50–0.55; 0.65 is the mse/ρ knee); `experiments/169-*/RESULTS.md` if present (sharpening); `experiments/173-sign-repair/RESULTS.md` (un-whitening rider, +0.002, closed).

## Why this experiment exists

Audit 181: **"production sits at the pds-optimal point; only new information can help" is WEAK** — every input to that claim was an un-whitening move. 170's criterion for a remaining `pds` lever is: **norm-preserving, raise X not A** (unique / discriminable energy, not amplitude), and it is unexplored outside sign repair.

This experiment tests ONE such lever, not a family.

## The lever: unique-component boost (then renormalize)

For each predicted effect row `h` (already whitened, as production ships):

```
c = mean_p h_p                  # predicted-effect centroid in the same space
u = h - α * proj_c(h)           # strip a fraction of the shared direction
h' = u * (‖h‖ / ‖u‖)            # preserve the row's L2 exactly
```

α = 0 is identity (known-answer). α = 1 is full anti-share then renormalize (AgingLab's public move, different operator). Sweep α ∈ {0, 0.25, 0.5, 0.75, 1.0}.

This is NOT un-whitening (the operator stays put). It is NOT amplitude (row norms are identical). It raises the share of each row that is unlike the others.

## What to do

1. Load the production-shaped predicted effects for the **challenge frames** (PLOLO / 162's instrument) and for **h1-test** (k562_gwps + orion_hek293t only — no h1_2025). Reuse existing frozen effect matrices; do not rebuild.
2. Apply the transform in the space the scorer reads (delivered / challenge frame), not in design space. 156 exists because that distinction is real.
3. Known-answer row: α = 0 must reproduce the incumbent pds to numerical noise. If it does not, abort.
4. Price `pds` on PLOLO (challenge frame) and on h1-test. Do not quote absolute Orion pds. Do not run `cell-eval2`.
5. Also report: mean pairwise cosine of the rows (should fall as α rises), mean cosine-to-truth if you have truth, and whether `‖h'‖ == ‖h‖`.
6. If the best α beats incumbent by ≥ 0.003 raw pds on BOTH PLOLO and h1-test, write a six-member price sketch (mse will move because the direction changed; nmae / DE are profile-class and unpriced — say so). Otherwise NO-GO.
7. Write `RESULTS.md`. Negative is a result. Do not modify production.

## Territory

182 owns pretrained State (new information). You own this one transform. Do not wander into emission, alpha, or donors.
