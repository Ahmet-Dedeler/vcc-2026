# 277 — Render the identifiable LOW-CONSERVATION rows differently: one-variable arms selected by 276's leave-the-held-line-out conservation map, priced on PLOLO and on the Flex-400 reference. ≤ 5 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/277-conservation-class-rendering/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold per process, peak RSS
under ~10 GB. **Time box: write up after ~5 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you (pessimistic
pricing, controls, paired bootstraps, name the unit, refit inside the split).

## Why (R276 §5–§6, R273 §4, R262 §5, F8 in `docs/FINDINGS-2026-09.md`)
276 built a per-target conservation statistic over 19 lines that reproduces (ρ +0.27 / +0.44), predicts transfer with
the held line left out (ρ +0.39 on PLOLO, AUC 0.69 for the anti-aligned group), and was confirmed on the reserved
HEK293T half-B truths (claim C1, ρ +0.361 [+0.257, +0.458]). **The bottom 30 % of rows by conservation carry 46.8 % of
the `pds` deficit and are identifiable BEFORE the held line is seen.** Every per-target *selector* and *dose* failed
(262–265, 271): the residual cannot repair a row whose whitened base points away from the held line. 276's stated
question is different: **can a different RENDERING of those rows — not a different donor, not a different dose — move
them across the 0.02 knee?** 273 measured two cosine effects on the reliable reference that are invisible on PLOLO:
the un-deflated CD4 residual reads +30 % `cos_agg` against +24 % deflated, and the soft noise floor +4.4 %. Neither is
priced. One variable per arm: the rendering applied to the low-conservation class, everything else `cd4svd5`.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 276's `RESULTS.md` §4–§6 and `work/conservation276.csv` (columns `C_lolo_<fold>`
— the held-line-left-out score per fold; `C_lolo_kolf2`; there is no H1-specific LOLO column: for h1-test compute the
score from sources that exclude H1 lines exactly as `build276.py` does, and say so); 262's `diag262.py` (groups, per-
target arrays); 273's `lib273.py` / `instrument273.py` (the exact scorer statistic on h1-test, union OFF, CD4 re-extracted
for its 85 targets — REUSE); 257's `lib257.py` (deflation on/off), 272's `op272.py` (the soft floor c = 1), 246's
`lib246.py` (`compose`, `score_effects`); `scripts/make_submission_pm_cd4.py` (where a per-target class mask would enter:
the CD4 deflation is applied to the whole residual in `deflate_cd4_residual` — a class-conditional variant would need a
row mask there; name it, do not make it); BR249 (the deflation's board point); READ-FIRST §3.

## Constructions (freeze in `render277.py` before measuring; the CLASS is the bottom 30 % of C_lolo per fold, declared once)
V0 `cd4svd5` (KA: reproduces 262's rows at 0.0e+00).
V1 **un-deflated CD4 on the low-conservation class**: for rows in the class, the CD4 residual enters WITHOUT the k = 5
   deflation (as `h1union` had it); all other rows deflated as `cd4svd5` (one re-solve of the CD4 fraction for the
   mixed residual — state how, and report the delivered amplitude).
V2 **the opposite**: deflated on the class, un-deflated elsewhere (the mirror control that must lose if V1's mechanism is
   real).
V3 **no CD4 residual on the class** (the residual removed where the base is unreliable; the base alone there).
V4 **the soft noise-floor operator (272 A′, c = 1) applied to the class only** (rows whitened by the damped operator, the
   rest by the incumbent operator) — if the composition allows a per-row operator choice cheaply; else state why not.
V5 **class boundary sensitivity** for V1 only: bottom 20 % and 40 % (reported, not selected on).
Controls for V1: a RANDOM row set of the same size (3 seeds); the bottom 30 % by 262's donor-donor agreement (the
pipeline's own proxy — V1 must beat it to justify the map); the top 30 % by conservation (must not gain).

## Measurement
PLOLO three cancer folds (`kolf2` apart): `cos_agg`, macro, `pds`, per-target `pds`, by class and by 262 group, the
`nmae` identity, delivered amplitude, paired target bootstraps, leave-top-3-out retention. **Flex-400 (h1-test, union
OFF, 273's instrument)**: the same for V1–V3 vs the union-OFF incumbent, with the class from the H1-excluded score.
**KA:** 262's `cd4svd5` rows; 273's incumbent rows; BR249 retrodiction with the recipe.

## Price (READ-FIRST §3)
Coherent pessimistic / point / optimistic; the deflation's own class factors (re-shape: k_cos 3.76 for a cosine loss /
1.36 point, k_pds 0.85 / 1.22) for V1–V3, the operator bands for V4; identity ×1.5 on a loss; amplitude leg if > 0.2 %.
CLEARED = pessimistic ≥ +0.0019 over `cd4svd5` AND `pds` > 0 on 3/3 folds with the interval above 0 AND retention ≥ 50 %
AND every control (random set, agreement-selected set, mirror) below half AND the Flex-400 reading does not contradict
the sign.

## Deliverable
`RESULTS.md`: KA rows; the class definition and its size per fold; the arm table (six members, three scenarios, per
class and total); the controls; the Flex-400 re-read; one line **CLEARS / DOES NOT CLEAR**; if it clears,
`PREREGISTRATION-DRAFT.md` in 249's form (five falsifiers; the production variable: a per-target class mask consumed by
`deflate_cd4_residual` and the C_lolo table it needs frozen as an artifact with a sha — named, not made). A clean
negative — "the low-conservation rows cannot be moved by rendering either" — closes 276's competition-facing framing.
No model name, no submission. Then stop.
