# 274 — The union's DIRECTIONS without its MAGNITUDE cost: design and price an `nmae`-neutral rendering of the six-source union on top of `cd4svd5`. The one competition candidate with a measured positive direction. ≤ 5 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/274-nmae-neutral-union-renderer/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold per process, peak RSS
under ~10 GB. **Time box: write up after ~5 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you (pessimistic
pricing, strong baselines, paired bootstraps, name the unit).

## Why (BR237 §3, R258 §6, R259 §0, R273 §0)
On the board, ARM6 (the six-source union at g = 0.35 replacing the H1 pair) **bought `pds` +0.00248 and `mse` +0.00088
overall and paid `nmae` −0.00716** — the union's directions are worth +0.0034, its magnitude on the Wilcoxon-significant
genes costs twice that (raw `nmae` 0.934 → 0.960). 258 priced an `nmae`-at-zero rendering of the same directions on
`cd4svd5` at **+0.0014 pessimistic / +0.0034 point** (k_pds 0.7 / 2.48), and 259 confirmed the cost is
`igvf_cardiomyocyte`'s while the pair + DLD-1 is `nmae`-cheap. Nobody has designed the renderer. The `nmae` identity
(216/188) is now board-validated (BR237: predicted −0.004783, realised −0.00716, same sign), so it can price a rendering
that changes magnitudes. **Question: is there a rendering that keeps the union's DIRECTION contribution on `pds`/`mse`
while holding the per-gene MAGNITUDES on the scored DE set at the incumbent's, and does it clear the pessimistic
gate over `cd4svd5`?** 273 adds a fact to respect: the challenge truth resolves directions (reliability 0.24–0.5), so
direction gains are real currency there.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; BR237 (`experiments/237-gate-validation/BOARD-RESULT.md`) §1–§3; 258's
`RESULTS.md` §0–§6 and `price258.py` (the union path on the `cd4svd5` base, the identity's pools, the price recipe —
REUSE; byte-copy its `work/plolo/` caches); 259's `price259.py` (`union2` path, if a second-residual variant is
tried); `experiments/229-stack4-build/six229.py` (the identity: `nmae_wald89`, the top-89 pools, `NMAE_SPAN`, 4-fold
pooling — KA-1 must reproduce ARM6's leg −0.004783); `experiments/188-magnitude-calibration/RESULTS.md` (what
L1-optimal magnitude calibration did and why it failed — do not repeat it; it recalibrated EVERYTHING, this renders
only the union's contribution); `experiments/228-theta-gain-plolo/RESULTS.md` (per-gene gain, dead on `nmae`);
`scripts/make_submission_pm_cd4.py` around the union merge and the preserve step (`apply_frozen_cd4_residual`,
`prepare_h1_union`, the merged add, the global preserve) — name the one place a renderer would enter; READ-FIRST §3.

## Constructions (freeze in `render274.py` before measuring; all on `cd4svd5` = CD4 deflated k = 5 re-solved + the H1 pair at 0.28)
R0 **baselines**: the incumbent (`cd4svd5`); the six-source union at g = 0.35 as ARM6 built it (KA: its identity leg
   −0.004783 reproduces; its PLOLO `pds`/cosine rows match 258); the `nmae`-at-zero bound 258 priced (the ceiling
   for any renderer).
R1 **magnitude-matched union (per gene, pooled set)**: compose base + CD4 + union as ARM6, then for every target rescale
   the effect's entries on the identity's top-89 pool (the genes the identity prices — and separately on the
   incumbent's own top-89) so that their per-gene |magnitude| profile equals the INCUMBENT's on that pool (rank-matched
   or ratio-matched: try both), leaving the direction and the off-pool entries as composed.
R2 **direction-only union**: take the ARM6 composition's per-target DIRECTION (unit vector) and the INCUMBENT's per-target
   norm and per-gene magnitude spectrum (sorted |values|), i.e. re-render the union arm with the incumbent's amplitude
   structure — the rotation without the magnitude change.
R3 **the cheap subset**: R1/R2 applied to the pair + `dld1_lowmoi` (259's `nmae`-cheap subset) at g ∈ {0.28, 0.35}, and to
   the six-source union at g ∈ {0.20, 0.28, 0.35}.
R4 **second-residual variant** (259's construction B, pair kept at 0.28 + `union2` of the added sources) rendered as R2.
Controls: R1/R2 with the union's direction REPLACED by the incumbent's own (must reproduce the incumbent — the
renderer's identity check); the union at its own magnitudes (R0); energy-matched random directions rendered as R2.

## Measurement (three cancer folds, per fold and pooled; `kolf2` apart)
Through the fold instrument: `cos_agg`, macro, `pds` (per target; the scorer's statistic as 262/273), the `nmae` identity
on its three pools, top-89 precision, delivered amplitude (must stay within 0.5 % of the incumbent's — a renderer that
moves amplitude is priced in the amplitude-changing class), paired target bootstraps, leave-top-3-out retention, by
262 group. **KA-2:** the recipe retrodicts BR237 by member (−0.0036 ± 0.001 with 258's factor set) and BR249.

## Price (READ-FIRST §3): coherent pessimistic / point / optimistic with the new-source bands for the union legs (k_cos
0.94 / 0.96, k_pds 0.7 / 2.48 by the sign of the leg, DE-three 0.0425 / 0.24, identity ×1.5 on a loss and ×0 on a gain).
CLEARED = pessimistic ≥ +0.0019 over `cd4svd5` AND `pds` > 0 on 3/3 folds with the interval above 0 AND retention ≥
50 % AND amplitude within 0.5 % AND every control below half.

## Deliverable
`RESULTS.md`: KA-1/KA-2; the renderer definitions; the table of every cell (six members priced three ways); the
best cell in full; the production change it would need (a renderer stage after the union merge and before the
preserve — one place, named, NOT made); one line **CLEARS / DOES NOT CLEAR**; if it clears, `PREREGISTRATION-DRAFT.md`
in 249's form with five falsifiers and the known-answer rows a build would need (flag-off byte identity with
`aba27758…`; ARM6's identity leg; the renderer's identity check). A negative — "the union's direction cannot be
rendered without its magnitude cost" — closes the last measured positive direction on this base and is a valid
result. No model name, no submission. Then stop.
