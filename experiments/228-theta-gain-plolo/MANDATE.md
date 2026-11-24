# 228 — Price 224's per-gene reliability gain on ALL SIX MEMBERS through PLOLO. This is 224 §10's specified successor.

You are a research subagent. Work in `experiments/228-theta-gain-plolo/`. **No `vcc`, no submission
build, no production change. No `cell-eval2` while another agent is running one — check first.**
Peak RSS ≤ 12 GB. `uv run python`. Five other agents are live (221, 224 finishing, 225, 226, 227);
write only inside your own directory.

## What you inherit

**Read `experiments/224-learned-transfer/RESULTS.md` in full before anything else, especially §10.**
224 asked whether any trained model beats the donor mean on held-LINE transfer and found that
**exactly one rung does — the smallest**:

> A **per-gene diagonal** (one scalar per gene, ~5,070 parameters), fitted by inner
> leave-one-training-LINE-out on `(donor consensus → held-out line)` pairs, beats the equal-weight
> donor mean by **+5.16% relative transfer cosine** on a 7-line surface (95% CI `[+4.63, +5.72]`,
> **7/7 held lines positive**) and by **+8.9% at production's `k = 3` on the 300-gene challenge
> panel**, at **Δ`pds` +0.0118**, **on top of production's own whitening `gamma = 0.65`** (+9.63%
> there, more than it adds to the unwhitened consensus). Pre-registered gate was +2.9%.

**Every higher-capacity rung lost**, two of them in this project's classic ways (132's
same-target-signal replacement; 195's clear-on-perturbations-fail-on-lines). The 25M-parameter
operator bought cosine by **inflating the decoy** (−0.032 `pds`). **The slope in `k` is negative**,
so capacity is not the binding constraint and **no GPU was rented**.

**The mechanism is identified and it is already ours.** The fitted gain **is 189's per-gene
cross-line reliability prior** — `r(gain, theta) = +0.81 / +0.90`, `r` with 224's own reliability
statistic **+0.94**, and gains fitted for different held lines correlate at **r = 0.95**.
**189 measured this in the ORDERING channel, where it bought nothing. Nobody had put it in the
DIRECTION channel.** All four nulls pass: gene-shuffle **−7.2%**, moment-matching **−3.5%**,
label-derangement **−25.7%**, wrong-recipient better than the arm. Leak audit clean.

**224 states loudly that none of this justifies a slot**, because it is an **effect-level**
measurement on donor tables — not run through `profile_from_effect` + the scorer read-back, so per
156 it measures an object the board never sees, and it is not a PLOLO number. **Your job is to turn
it into one.**

## The task, as 224 §10 specifies it

1. **Wire the gain in as a per-gene multiplier on the production effect matrix.**
2. **Fit it leave-one-donor-LINE-out on NON-PANEL targets only.** 224 §7: the panel surfaces read ~0
   from **power, not substance** — 7,312 gains cannot be fitted from 73–120 panel targets; fitting on
   **1,758 non-panel targets** and evaluating on **232 panel targets** of a held line recovers
   **+8.3% ± 0.3 across five seeds**. Keep that discipline; it is also what makes the arm
   **October-robust** (see §5 below).
3. **Price it on ALL SIX MEMBERS through PLOLO in the challenge frame**, three cancer folds, with
   `kolf2` reported **separately** (it is pluripotent, gains 3–5×, and inflated a coordinator
   prediction on 2026-09-04 badly enough to force a band reissue).

## The two things 224's surface could not see — these are the crux, do not skip them

1. **The `VCC_ALPHA` interaction.** Rescaling per-gene magnitudes changes the optimal *global*
   amplitude, and **211 found `claim` already at 100.0000% of its cap** (permitting −0.12%). So the
   gain may be partly spending amplitude we do not have. **Sweep the gain jointly with `VCC_ALPHA`**;
   report the gain at fixed α **and** at re-optimised α, and say how much of it is amplitude.
2. **The DE call set moves.** The gain changes per-gene effect sizes, therefore which genes clear
   significance, therefore `fid` and `jac`. **Price them; do not carry them at zero** — that is
   §8.4's entire lesson, and it is why a direction arm priced `mse`-only is understated ≈**2.51×**.

## Rule 13 and the toll
224 expects a **profile-class** arm, so the `reach` toll should be **0** per 169. **Do not assume
it — run the test.** `AGENTS.md` §11.1: the **per-gene panel-total digest** is the test, and it is
single-valued. Emission-law profile-CHANGING costs **−0.0146 scaled**; profile-class costs 0. Report
the digest.

## Falsifiers, all mandatory
1. **Re-run 224's four nulls through YOUR instrument** — gene-shuffle, moment-matching,
   label-derangement, wrong-recipient. If they do not reproduce their signs, your wiring is wrong.
2. **Leak audit.** Fit on non-panel targets only; state every donor↔surface overlap. `h1_2025`
   covers **85/85** of h1-test's targets and is inside the shipped union. `debench.Surface.pds()`
   must **never** size an arm (188 Phase 3).
3. **Leave-top-3-out retention on cosine and `pds` separately.** On `h1union` the board realised
   ~**100%** of the cosine prediction and only **70%** of the `pds` one.
4. **Marginal over the CURRENT shipped stack** — three donors + CD4 (`f = 0.34`) + the H1 union
   (`g = 0.28`), board **0.2118310**. A gain that duplicates what the stack already supplies is
   worth zero (208's rider rule).

## Reporting
**Do not convert PLOLO to board points with a scalar** — the factor is **not identified**: 1.384 on
`cd4f034`, 1.000 on `h1union`, matched `cos_agg` functional. Quote PLOLO units and **name the fold
set, the functional, the gate and the tag**. 221 owns the reconciliation. 224's own §10 estimate
(≈+0.0029 overall) uses §8.2's exchange rate and is explicitly *"if the effect-level gain survives"*
— **treat it as a hypothesis to test, not a number to reproduce.**

## §5 — October, and why this arm may be the best one we have
The panel is **replaced 2026-10-22 with 300 targets disjoint from today's**, and **only the last
entry after that date is prize-eligible** (`DEF-RUNBOOK` §0, §5). A gain fitted on **non-panel
targets** is **panel-independent by construction** — unlike every candidate in 225's source ladder,
which is selected against today's panel. **State this explicitly and quantify it if you can**
(e.g. refit on a disjoint synthetic panel via `def-bundle` and show the gain survives). That
property may be worth more than the headline number.

## Deliverable
`RESULTS.md`: the six-member PLOLO price at fixed and re-optimised α, the digest/toll result, the
four reproduced nulls, retention, the marginal over the shipped stack, the October statement, and a
one-line verdict. **"It does not survive the scorer read-back" is a complete and valuable result** —
156 exists precisely because effect-level gains can vanish there.

**Do not build a submission package and do not choose a model name.** If it clears, the coordinator
decides the slot.

## Territory
221 calibration · 224 (finishing) learned transfer at effect level · 225 sources-as-donors ·
226 the operator's fitting sample · 227 `nmae` and emission. **You own the PLOLO pricing of 224's
gain.** Note 226 also touches the operator: **it varies the operator's fitting ROWS, you apply a
per-gene multiplier to the effect.** Different objects — do not re-fit the operator.
