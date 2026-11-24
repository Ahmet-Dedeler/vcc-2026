# 265 — The proper statistical check of 263's `kolf2` majority result (and 263's two-donor cells): paired target-level uncertainty and a permutation test over drop assignments. Analysis of saved arrays only, ≤ 90 minutes. No fold rebuilds, no build, no submission.

You are a research subagent. Work in `experiments/265-majority-permutation-test/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. This is arithmetic on arrays 263/264 already
saved; peak RSS under ~4 GB. **Time box: 90 minutes.**

## Why
263 reported the three-donor majority rule on `kolf2` at +0.0121 `pds` "8 σ" and 264 called its random control
+0.0068 "56 % of the arm" from three seeds. Both statements outrun the measurement: the "σ" divides a per-fold Δ`pds`
by 162's 0.001494, which is the leave-one-arm-out RMS of PLOLO→board prediction error, not a per-fold sampling σ; and
three random seeds cannot establish what fraction of a gain is generic. An outside re-analysis of the saved per-target
scores (GPT, 2026-09-05) gives majority − incumbent **+0.01209 with a target-bootstrap 95 % interval ≈ [−0.0034,
+0.0280]** and majority − mean(random drops) **+0.0053 ≈ [−0.0115, +0.0228]**. Settle it properly, so the decision
"is a correctly re-fitted three-donor validation worth building?" rests on a real test.

## READ FIRST
`docs/READ-FIRST.md` (v3.4, §3: the uncertainty rule), 263's `RESULTS.md` §7 and `lib263.py` / `sweep263.py` /
`analyse263.py` (where the per-target `pds` arrays and the per-arm caches live under `work/`), 264's `RESULTS.md` §0.8 and
its `work/` (the drop-one random seeds on `kolf2`, the matched cancer-fold construction), 262 KA-3 (the per-target rank
statistic recomputed from the saved prediction × truth cosine matrices equals the instrument's — **so a target's `pds`
depends only on ITS prediction row against all truths**; verify this on 263's arrays before relying on it).

## The test
1. **Per-target additivity (KA).** For `kolf2`, precompute per-target `pds` for each of the three possible single-target
   choices (drop donor j for target t, norm-matched as 263's `:nm`) — 3 × 240 rows from the saved cosine matrices or
   by re-scoring single-target changes; verify that the majority arm's total equals the sum of its per-target choices
   (0.0e+00 against 263's `majority:raw:all:nm` row) and likewise for the incumbent and the three random-drop seeds. If
   additivity fails, say so and stop at the bootstrap (step 2) with the full-arm evaluations only.
2. **Paired target-level bootstrap** (≥ 10,000 resamples of targets): majority − incumbent; majority − each random-drop
   seed and − their mean; the majority's excess over the always-drop-`k562`, -`hct116`, -`hek293t` rows. Report the
   95 % intervals and the fraction of resamples ≤ 0.
3. **Permutation over drop assignments** (needs additivity): ≥ 10,000 uniformly random drop assignments, and ≥ 10,000
   FREQUENCY-PRESERVING shuffles (permute which target gets which dropped donor, keeping each donor's drop count exactly
   as the majority rule's); the majority's `pds` against each null: p-value and the null's mean and 95 % band; the same
   for the τ = 0 gate (129 targets). Also report the oracle drop's position in the same nulls (it must be extreme).
4. **The cancer folds, 264's matched construction**: the same bootstrap for (b) − (a) and (b) − (c) per fold, and the
   frequency-preserving permutation where additivity holds; and 263's two positive reverse cells (`cd4d:raw:tau0:rev:nm`,
   `union:raw:tau0:rev`) and the best two-donor arbiter cell, bootstrapped.

## Deliverable
`RESULTS.md`: the additivity KA; a table of every comparison with point estimate, bootstrap 95 % interval, permutation
p-value (uniform and frequency-preserving), null mean and band; and one line answering: **does the majority rule's
excess over frequency-preserving random drops survive at p < 0.05 on `kolf2`?** If yes: state what a correctly
re-fitted three-donor validation would need (operator re-fit including the third donor, 244's recipe, ~30 min per fold)
and recommend it or not. If no: H2's majority variant closes on its own numbers. No price, no preregistration, no
model name, no submission. Then stop.
