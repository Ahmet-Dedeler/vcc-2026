# 351 — PLAN: board forensics + 2025 winners + proposals

## Hypothesis
- A: The public board's full snapshot history contains a readable signal about HOW the top
  teams (esp. leader "ViRtUaL cElL builders" / PromoterAI) gained: which members carried
  ≥+0.005 steps, whether nmae gains pair with mse losses (our 317 amplitude trade) or not,
  and where our profile (best 0.21866) is anomalous vs the top-30 median.
- B: The 2025 VCC winners' public write-ups/code name a mechanism (model class, training
  data, cell/count generation, DE-metric handling) that suggests a cheap, not-yet-closed
  experiment on our H1 harness (323).
- C: At most three concrete, falsifiable H1-harness experiments follow; "nothing new" is
  acceptable.

## Exact inputs (paths)
- Board snapshots: `data/leaderboard/*.json` (69 files, 2026-08-21 → 2026-09-09T2137Z).
  Schema per snapshot: `{"entries": [...], "num_entries": N}`; per entry: `id`, `team_name`,
  `team_id`, `model_name`, `submission_date`, `status`, `score_avg`, `score_pds/mse/nmae/fid/reach/jac`,
  raw `pds_cosine`, `expr_mse_unbiased_capped_norm`, `de_wilcoxon_*`.
- Own-entry ground truth: `data/board-entries/*.json` (9 files incl. `xuDelEik7zGWMXdNu0R4.json`
  = 319 at 0.21865973944346428).
- Ledger: `docs/BOARD-LEDGER.md`; handoff gap table: `docs/SPRINT-HANDOFF-2026-09-09.md`
  (leader 0.2630424866 − us 0.2186597394 = 0.0443827472; nmae +0.02282, reach +0.01868…).
- 323 harness known answers: `experiments/323-current-family-validation/work/metrics-cap15.json`
  (`raw.nmae 0.9508364630325662`, `raw.reach 0.25434726355678233`, etc.).
- 343/344: `experiments/343-nmae-decomposition/RESULTS.md`, `experiments/344-null-key-killcheck/RESULTS.md`.
- Prior-experiment titles for Part C: `experiments/*/RESULTS.md` (grep).
- Web: search/fetch only, no downloads >2 GB (none planned).

## Estimator / procedure (Part A, local, no harness runs)
1. Load all 69 snapshots with `json` stdlib only (<2 MB each, serial, RSS < 500 MB).
   Dedupe entries by `id` (keep first-seen date, last-seen snapshot, full latest record).
   Sort each team's entries by `submission_date`; team's BEST = max `score_avg`.
2. Known-answer reproduction BEFORE new measurement (own code):
   (a) 319's six scaled members + overall from `data/board-entries/xuDelEik7zGWMXdNu0R4.json`
       vs the same id's row in the latest snapshot containing it — must match to ≤1e-9;
   (b) 323 `raw.nmae 0.950836` re-read from `metrics-cap15.json` (reproduces 343 §"recomputed mean NMAE");
   (c) leader-minus-us gap decomposition recomputed from the two best rows — must match the
       handoff table to ≤1e-6 overall-points.
   If any fails: stop and report.
3. Top-30: rank teams by best `score_avg` in the LATEST snapshot (2026-09-09T2137Z, 840 entries);
   emit best-entry six members + overall, median of each member over the 30, and HIGH/LOW profile
   flags vs median. Locate "Ahmet Dedeler" best (expect 0.21866, rank ~13).
4. Consecutive pairs: within each team, order entries by `submission_date` (drop equal-timestamp
   ties); Δ = later − earlier on overall + six members. Flag improvements ≥ +0.005, record carrying
   members (largest Δ contributions = Δmember/6). Quantify Pearson corr(Δnmae,Δmse),
   corr(Δreach,Δpds), corr(Δnmae,Δreach) over ALL consecutive pairs (scaled members), plus the same
   restricted to improvers (Δoverall>0) and to top-30 teams. Report n with pairwise-complete observations.
5. Leader trajectory: all ids with team_name "ViRtUaL cElL builders", ordered by date: model_name,
   date, overall + six members. Note test1…PromoterAI_precomputed sequence.
6. Jaccard: list every best-entry with `score_jac` > ours (0.0286846490); show their other members.
7. Outputs: `work/entries.csv` (one row per distinct id), `work/top30.csv`, correlation table in
   RESULTS.md (also `work/correlations.json`, `work/improvements.csv`, `work/leader_trajectory.csv`).

## Controls / falsifiers
- Dedupe check: distinct-id count ≥ latest snapshot's 840; every latest-snapshot id present in deduped set.
- Sum check: `mean(six scaled)` == `score_avg` to ≤1e-9 on emitted rows (spot-check 30/30).
- Date sanity: first-seen ≤ last-seen; submission_date parses as UTC ISO.
- Correlation guard: report Spearman alongside Pearson (rank-robustness); report n; do NOT call a
  board correlation a causal lever or a score forecast.
- Scope guard: board numbers are A/B/C validation evidence; H1-panel numbers are development
  evidence (343). Never convert a retrieval/cosine/proxy gain into a score gain.

## What kills the idea (per part)
- A: kill = dedupe or known-answer reproduction fails (then PARTIAL, report only what reproduces).
- B: kill = no public 2025 winner write-up/code found after documented search (then "nothing new" path).
- C: kill = every candidate maps onto a CLOSED experiment (312/317/227/088/232/234/344/342/318/229/331–341
  list in assignment) with no material difference; then verdict "Nothing new" is the honest result.

## Budget / resources
- ~3 h total. Part A ~60–90 min (script + checks). Part B ~60 min (search/fetch, ≤15-word quotes).
  Part C ~30 min (grep + falsifier design). No builder/scorer/submit/cloud. Peak RSS target <1 GiB.
