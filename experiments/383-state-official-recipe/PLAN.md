# 383 — STATE official recipe: PLAN (30 lines max)

Hypothesis: the official STATE recipe (model=state, ESM2 pert features, 40k steps,
infer from step-20000) yields a checkpoint WITH perturbation identity (held-line
cosine-rank pds >= 0.7) whose H1 deltas are decorrelated enough from cap15 that a
50/50 blend raises truth cosine. Prior state_sm (106/113: pds 0.50) and one-hot
pretrained (182: 0/300 coverage) do not price this recipe.
Inputs: notebook data/external/state-vcc/official_vcc_notebook.ipynb (sha in 106
PROVENANCE.json); repo data/external/state-vcc/state @9bbfe78 (read-only); support
set https://storage.googleapis.com/vcc_data_prod/datasets/state/competition_support_set.zip
(8.7 GB download > 2 GB — noted here; goes to experiments/383-state-official-recipe/data/;
verify sha256 f1d5fd56...507f7c vs 106 bundle-manifest.json); harness data/harness/h1-test.h5ad;
challenge data/challenge/context_{A,B,C}.h5ad + pert_counts.csv + gene_names.csv.
Estimator: (a) recipe training per notebook cell 14, model=state, MPS via 106 mps patch
logic copied into my dir; (b) identity = cosine-rank pds on recipe's held hepg2 test
perturbations; (c) H1 six members via 323 score() frame reused read-only, outputs in my dir.
Controls/falsifiers: KNOWN ANSWER FIRST — recompute score() from 323 stats-cap15.npz and
match experiments/323-current-family-validation/work/metrics-cap15.json raw to 1e-6; if no
match, stop. Download hash must match manifest or stop. Throughput probe (200 steps) before
committing to a multi-hour run.
Kill criteria: (i) held-line pds <= 0.55 after the recipe fraction I can train -> STOP at
deliverable 3, no A/B/C inference (106's trap). (ii) measured full-state MPS rate projects
step-20000 beyond ~6 h -> report rate + largest feasible run, no cloud (per stop rules).
(iii) memory_pressure free < 40% or pgrep def_pipeline alive -> wait, one heavy job at a time.
Deliverables in order: RESULTS.md recipe+provenance; checkpoint (or measured infeasibility);
identity number; H1 six members vs cap15 raw; state_delta_{A,B,C,h1test}.npz; blend cosines.
