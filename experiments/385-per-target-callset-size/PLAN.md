# 385 — per-target call-set size: PLAN

Hypothesis: n_real,p (reference significant genes per target, scorer's gate) is
predictable from build-time observables, so K_p = clip(c·n̂_p, 400, 4474) beats
constant K = 2,600 (= CallSet(400, 2200, 1.9)) on jac+fid.
Inputs: 323 work/{reference.npz, effects.npz::cap15, stage.npz (pre-whitening
donor-consensus norm), h1-anchors.json, metrics-cap15.json}; 360
work/CRISPRGeneEffect_subset.csv (DepMap essentiality); rc_emit.CallSet.
Estimator: (1) known-answer: own cap15 emit+score via 380's call sequence,
exact digest match. (2) per-target n_real vs observables: ‖effect_p‖ row norm,
#{g: CPM-gate ∧ |effect|>τ} for τ ∈ {0.05, 0.1, 0.2, 0.5}, pre-whitening row
norm (stage.npz uniform arm), DepMap |GeneEffect| of target; Spearman +
LOO OLS of log(n_real+1). (3, gated) price two c values vs constant K and a
shuffled-K control; paired-target bootstrap (report380 pattern); H1-scaled via
323 anchors. Per-target count realised as CallSet(n_head_p, n_cpm_p, 1.9) with
n_head_p = min(400, K_p), rest by CPM — protected 400-|H| head kept inside.
Controls/falsifiers: shuffled-K control (same K multiset, random targets);
192 §3.3 prior (Pearson +0.24 → expect kill); mse_legs + pds must not move
(call set preserves panel totals by construction).
Kill: best predictor Spearman < 0.4 → report NEGATIVE, stop, no pricing runs.
Kill: rule fails shuffled-K control on jac+fid by ≥ +0.005 H1-scaled → NEGATIVE.
Budget: ≤10 GiB, threads 4, ≤4 harness runs, one at a time; memory gate first.
