# 262 — DIAGNOSTIC: where does `cd4svd5` lose `pds`, target by target, and which structural property explains it? No pricing, no build, ≤ 2 hours.

You are a research subagent. Work in `experiments/262-where-pds-is-lost/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold at a time, peak RSS under ~8 GB.
**Time box: 2 hours wall clock, then write up.** This experiment PRICES NOTHING and proposes no candidate as cleared;
its output is a ranked diagnosis and at most three one-variable hypotheses, each with the cheap screen that would test it.

## Why
`pds` is the member where first place is furthest ahead (`docs/READ-FIRST.md` §1: our scaled 0.77809 vs the leader's
~0.826) and the only member tonight's arms bought. Every source arm measured since the board points puts its `pds`
on **three or four targets per fold** (R257 §0.5–0.6, R258 §6, R259 §5, R261 §6), and the CD4 projector's `pds`
lives in one direction (R257). We do not know WHICH targets `cd4svd5` discriminates badly, or why. The scorer's `pds`
is per-perturbation (`docs/reference/vcc2026-metrics-BRIEF.md`; AGENTS.md rule 8: never pool across the axis the
metric refuses to average across), so the loss has an address.

## READ FIRST
`docs/READ-FIRST.md` (v3.0), `docs/reference/vcc2026-metrics-BRIEF.md` (the exact `pds` definition the scorer uses),
`experiments/162-retrodicting-validator/plolo.py` (the instrument; per-target `pds` arrays exist — 257/258/259/261
saved `pertarget*.npz` per fold), `experiments/257-soft-deflation-price/` (`lib257.py`, `sweep257.py`; byte-copy its
`work/plolo/` caches), `experiments/246-direction-hypotheses/lib246.py` (`compose`, `score_effects`, coverage helpers),
`experiments/243-transfer-factor/RESULTS.md` §4.3 (k_pds by class), `AGENTS.md` §2.1 (how the panel was built:
`n_real`, essentiality), `scripts/make_submission_pm_cd4.py` around the donor registry (which donors carry which target).

## What to measure (three cancer folds; `kolf2` apart; the `cd4svd5` construction as 257/261 build it, k = 0 of 261 == R1)
1. **Per-target `pds` for `cd4svd5`** on every fold (the scorer's definition, from the instrument's per-target arrays):
   distribution, the bottom and top deciles by name, fold-to-fold agreement of the bottom decile (is it the same
   targets, or fold-specific?).
2. **For every target, the structural covariates** available without new data: number of production donors carrying
   it and which; CD4 carried (complete / partial / absent) and its CD4 row norm; H1-union carried; the prediction's
   cosine to the fold truth; the truth's energy (‖Δ‖) and the prediction's ‖Δ‖ ratio; the target's essentiality flag
   as AGENTS.md §2.1 defines it; the prediction's nearest-neighbour cosine to OTHER targets' predictions (how
   confusable it is) and the truth's nearest-neighbour cosine to other truths (how confusable it SHOULD be).
3. **Which covariate explains the loss**: rank the covariates by how well they separate the bottom decile from the rest
   (AUC or a rank correlation with per-target `pds`, per fold and pooled by rank); report effect sizes, not a fitted
   model. In particular: is low `pds` a property of (a) the prediction being confusable with other predictions
   (the residual-hub story, which is what the deflation fixed — measure how much of that remains), (b) low coverage /
   weak donor support (nothing to predict from), (c) small true effects (nothing to discriminate), or (d) magnitude
   error at correct direction? Quantify the share of the fold's `pds` deficit (vs a perfect-rank ceiling, state it)
   attributable to each group.
4. **The `cd4svd5` step itself, per target**: which targets gained `pds` from `h1union` → `cd4svd5` and which lost;
   whether the losers share a covariate (this is the direction the next arm must not repeat).

## Deliverable
`RESULTS.md`: the known-answer row (the `cd4svd5` construction reproduces 261's R1 / 257's `k5_s1` at 0.0e+00); the
per-target tables (bottom/top deciles with covariates, per fold); the covariate ranking with effect sizes; the deficit
decomposition; and **at most three one-variable hypotheses**, each stated as a mechanism, the covariate group it
targets, the builder variable it would use (existing or named), and the PLOLO screen that would test it in ≤ 2 hours
with its known-answer row and controls. **No price, no "clears", no preregistration.** A result of "the loss is
diffuse and no covariate explains it" is valid and useful. No model name, no description, no submission. Then stop.
