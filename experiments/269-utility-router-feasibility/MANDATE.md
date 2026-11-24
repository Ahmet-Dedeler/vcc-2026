# 269 — Selective prediction: can a small model learn WHEN to act on a target — keep the incumbent row, switch donor, or apply a correction — under leakage-safe held contexts AND held target families? Feasibility stage, ≤ 4 hours. Measurement only.

You are a research subagent. Work in `experiments/269-utility-router-feasibility/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. Peak RSS under ~6 GB (siblings
267/268/270 run alongside). **Time box: write up after ~4 hours.** Read `docs/RESEARCH-PROGRAMME-2026-09.md` first; its
rules bind you (name the unit, refit inside the split, strong baselines, shuffled controls, paired bootstraps,
report the two generalisation settings SEPARATELY).

## Why
263–265 showed that per-target donor selection has a large oracle ceiling (+0.046 `pds`) and that no pairwise-agreement
arbiter captures it; that closed those arbiters, not a learned model that could use target biology, response reliability
and training across contexts. 066 (`experiments/066-leakage-safe-target-gate/RESULTS.md` §60–70) is the discouraging
prior: a nested gate captured 6–9 % of an oracle's headroom and even a same-target lookup did not transfer HEK293T →
HCT116. So this is a FEASIBILITY test, not an assumed breakthrough. **The formulation is utility, not pick accuracy:**
predict how much a specified action would improve or harm the target's `pds`, act only where the predicted gain is
confidently positive, and ABSTAIN = keep the incumbent row (never "predict zero").

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 066's `RESULTS.md`; 262's `RESULTS.md` §4–§5 and its per-target covariate arrays
(`work/cov262-*.npz`, `diag262.py`: coverage, CD4 carry and norms, union carry, donor-donor agreement, norms, nearest-
neighbour cosines, CD4 DE count); 263's per-target arrays (`work/cov263-*.npz`, `identity263-*.npz`: per-target `pds`
of the incumbent, of each always-drop-j / choose-j arm, of the oracle) and 264's (`work/cov264*.npz`, the matched
construction); 257's per-target arrays (`work/pertarget257*.npz`: `cd4svd5` vs `h1union` per target — the deflation
on/off action); `data/external/depmap26q1/extracted/` (gene-effect matrix: a target's co-essentiality profile is its
biology feature; check what is there) and `data/external/omnipath-dorothea/`, `data/external/regulon/` (pathway /
regulon membership, if usable); `AGENTS.md` §2.1.

## Stage 1 (this mandate)
1. **Labels (frozen action set), per (target, context)** on the four PLOLO contexts (K562, HCT116, HEK293T, KOLF2 apart):
   u₀ = 0 (keep the incumbent row); u_j = Δ`pds` of choosing donor j (norm-matched) for that target alone (263's
   choose-j arrays; note 265: per-target `pds` is not exactly additive through the pipeline — use the per-target arrays
   as the labels and say so); u_defl = Δ`pds` of the deflation on vs off for that target (257's arrays). State how many
   (target, context) rows exist and the label distributions.
2. **Features (build-time-observable only; NOTHING from the held line's truth):** target biology — the DepMap gene-effect
   profile of the target reduced to ≤ 20 PCs fitted on the TRAINING targets only, plus pathway / regulon membership
   counts if on disk; response reliability — donor-donor agreement, donor row norms, CD4 DE count, CD4 / union carry,
   nearest-neighbour cosine among predictions; recipient basal programmes — the held line's control profile projected
   on ≤ 10 programmes fitted on the training lines' controls. List every feature.
3. **Splits, reported separately:** (A) leave-one-CONTEXT-out over the three cancer contexts (train on two, test on the
   third; KOLF2 as a fourth training context only, never a test); (B) held TARGET FAMILIES: cluster targets by DepMap
   co-essentiality into ~10 families fitted on all targets' profiles (allowed: profiles are not labels), hold out whole
   families within each context. (C) both at once (held context AND held family) — the honest number.
4. **Models, small:** ridge and a shallow gradient-boosted regressor (≤ 50 features, depth ≤ 3) predicting each action's
   utility; decision = the action with the highest predicted utility if its predicted gain exceeds a margin τ (sweep
   τ, report the curve), else keep. **Baselines:** keep everything; every fixed action; the oracle (ceiling); SHUFFLED
   FEATURES (rows permuted within context); shuffled TARGET labels; 066's nested gate if reproducible. Metric: the
   decision-level Δ`pds` on the held rows (sum of realised utilities of the chosen actions, per context, with a paired
   target bootstrap), and the fraction of the oracle's headroom captured; also pick accuracy for reference only.
5. **Failure analysis:** for the held rows where the model acted and lost, which features drove it; for the targets the
   oracle gains most on, whether ANY feature separates them (AUC) — this is the scientific output even on a negative.

## Deliverable
`RESULTS.md`: labels and features tables; the three split settings with baselines and bootstraps; the τ curve; the
failure analysis; one line per setting: **BEATS FIXED + SHUFFLED on held data: YES / NO (captured … % of the oracle)**,
and the go / closure for stage 2 (more contexts from 268's data, or closure of this construction). No model name, no
submission. Then stop.
