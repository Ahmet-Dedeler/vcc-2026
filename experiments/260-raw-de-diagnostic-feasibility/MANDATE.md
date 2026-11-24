# 260 — FEASIBILITY ONLY: a private raw-count diagnostic with the scorer's REAL Wilcoxon DE sets, on K562-GWPS. 90 minutes. No harness build, no scoring, no model change.

You are a research subagent. Work in `experiments/260-raw-de-diagnostic-feasibility/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2` scoring run, nothing under `submissions/` or `data/harness/`. `uv run python`. Peak
RSS under ~6 GB (a sibling, 259, is running PLOLO folds). **Time box: 90 minutes wall clock, then write up.** You MAY
download `ReplogleWeissman2022_K562_gwps.h5ad` (8.8 GB, `scripts/scperturb_files.tsv` row 50, via
`scripts/download_scperturb.sh` or the URL there) into `work/` if — and only if — the per-target cell counts cannot be
read from anything already on disk; read it with backed/`h5py` access, never load `X`.

## The question
Every local `nmae` number we have is a PROXY: 216/188's exact identity uses fixed top-89 gene pools per target
(`POOL_N` = 89) where the board scores the reference's actual Wilcoxon-significant sets (`docs/reference/vcc2026-metrics-BRIEF.md`).
Tonight the identity was validated at scale on the board (ARM6: predicted −0.004783, realised −0.00716,
`experiments/237-gate-validation/BOARD-RESULT.md`), so it is a conservative bound, not a miss. The only local surface
with real held-out truth and real DE sets is `data/harness/h1-test.h5ad` (`scripts/make_harness.py`: 400 cells per
perturbation, ~20,000 UMI/cell median, scored by the real `cell-eval2` on the cloud, `scripts/cloud_score.sh`), and it is
LEAKED for every current arm (`data/donors/h1_2025` covers 85/85 of its targets, `AGENTS.md` §"TWO TRAPS"). **Is a
K562-GWPS-based diagnostic feasible as a NON-leaked surface for one narrow question — does the proxy miss a
material DE penalty on one frozen comparison — and what would it cost?** Decide GO / NO-GO with a time estimate;
do not start building it.

## READ FIRST
`docs/READ-FIRST.md` (v2.7) §4 (instruments), `AGENTS.md` §"TWO TRAPS IN h1-test" and the harness paragraphs around
lines 540–600, `docs/HARNESS-DEPTH-CONFOUND.md` and `scripts/check_harness_depth.py` (why a depth-mismatched harness
produces 90 % negative DE calls), `scripts/make_harness.py` (the 2026-faithful downsampling recipe), `scripts/cloud_score.sh`
and `experiments/221-*/RESULTS.md` §0 (the real scorer on the cloud, $1.48 per score, bit-exact), `docs/reference/vcc2026-metrics-BRIEF.md`
(what `nmae`/`fid`/`reach`/`jac` need from the reference), `scripts/make_submission_pm_cd4.py` around the donor registry
and PerturbMatch (lines ~840–1010: which production inputs derive from K562), `data/donors/k562_gwps/meta.json`.

## What to establish (numbers, with the file each came from)
1. **Availability**: what of K562-GWPS is on disk (the compact donor is pseudobulk effects, not cells); whether raw
   counts require the 8.8 GB download; its median UMI/cell and gene coverage against the challenge's ~20,000 and 18,533.
2. **Panel overlap and cell counts**: for the 300 A/B/C panel targets (`data/challenge/pert_counts.csv`), how many are
   in K562-GWPS, and their per-target cell counts (how many reach 150; how many 400); whether `make_harness.py`'s
   400-cell / 20k-UMI recipe can be met, and what a 32–40-target diagnostic would look like.
3. **The leak-audit surface, enumerated**: every production input that is K562-derived — the `k562_gwps` donor, the
   `k562_essential` union-registry source, the union whitening operator (fit on which donors? `data/donors/whitener_union_pruned.json`),
   PerturbMatch's K562-slot replacement, anything else — and what "a model with every K562-derived input removed" would
   be. State plainly whether that is still the production model or a different one (then the diagnostic calibrates
   the instrument, not the shipped prediction).
4. **Cost**: wall-clock and machine steps for: download, harness build (`make_harness.py`-style with
   `check_harness_depth.py` acceptance), a K562-free rebuild of ONE comparison (two predictions), two cloud scores
   (`cloud_score.sh`, ~$1.50 each), and the reading. Where the two predictions' `nmae` on the real DE sets would be
   compared with the identity's proxy.
5. **Verdict**: GO (with the frozen comparison you would preregister and the total time) or NO-GO (which step makes it a
   major project: extraction, the exclusion audit, depth, or coverage). GPT's prior from outside the repo: many targets
   have ~150 cells, few 400; K562 pseudobulk has already shaped development, so it is not an untouched holdout.

## Deliverable
`RESULTS.md` with the five sections above, every number sourced, and one line **GO / NO-GO**. No harness, no scorer
run, no builder change, no submission, no model name. Then stop.
