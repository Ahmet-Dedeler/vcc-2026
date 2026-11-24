# 285 — The apples-to-apples reliability benchmark: which source of disagreement between two measurements of one perturbation dominates at the challenge's budget — cell sampling, UMI depth, the guide construct, the donor, the batch, or the study — measured on MATCHED targets and response strata, in response agreement AND the scorer's own statistic, on the cell-level data on disk plus one bounded CD4 extraction. Measurement only, ≤ 7 hours. No build, no submission.

You are a research subagent. Work in `experiments/285-reliability-benchmark/`. **Edit nothing outside it** and write
nothing under `data/`; any download goes to `work/` (hard cap **20 GB on disk at any moment**, **40 GB transferred in
total**; `df -h` before each part; delete what you finish with). `uv run python`, backed / streaming reads only, peak RSS
under ~10 GB; any pass ≥ 8 GiB one at a time, and only after `memory_pressure` reports ≥ 40 % free (hold otherwise — two
sibling experiments, 284 and 286, run in parallel; do not read or write their directories). **Time box: 7 hours.**
`docs/RESEARCH-PROGRAMME-2026-09.md` binds you (name the unit; preregister before looking; strong baselines). Both
reserves are SPENT: never read `experiments/267-*/work/reserve/` or the Broad-16 lines JVE127 / MG63 / SLR23 / UMRC3.

## Why (third external review, 2026-09-06)
The programme's reliability numbers come from different designs that do not line up: 267 (batch-parity and random
cell halves at 70–95 cells, Orion / K562-GWPS), 273 (lane halves at 200/200, H1 Flex), 280 (donor halves at ~250 and
guide halves at 347 cells, CD4; its "depth curve" thinned COUNTS from pseudobulks, not cells — retracted 13:40Z), 282
(donor × guide × state at ~66–107 cells). No design on disk puts the units side by side on the SAME targets at the SAME
cells per side, and none separates cells per condition from UMIs per cell. **The deliverable is not another universal
reliability number. It is the decomposition, per source, of which disagreement a better estimator or more measurement
can actually remove**: cell sampling and UMI depth are removable by measuring more; the guide construct, the donor and
the study are not removable by more cells of the same construct. Across-donor agreement is also not the challenge's
within-reference split-cell repeatability; this benchmark measures both and says which is which.

## READ FIRST
`docs/READ-FIRST.md` §4 (instruments) and §6 (the 13:40Z corrections); `docs/FINDINGS-2026-09.md` F1, F10, F11;
267's `RESULTS.md` §1–§3 and its scripts (`inv267.py`, `half267.py`, `truth267.py`: the streaming pseudobulk passes,
the residual space, the batch / random / guide splits; its cached half-banks under `experiments/267-…/work/half/<line>/`
and the reliability rows `work/half/truth-<fold>.npz` `rel|*`); 273's `lib273.py` / `truth273.py` (the scorer-exact
per-target `pds` statistic — reproduces cell-eval2 to 5.6e-17 — and the h1-test halves by lane, random cells and
guide); 280's `RESULTS.md` §1 (the GWCD4i objects, sizes, obs columns, the HTTP-range reader `work/remote_h5.py`,
`work/parts_file.py`), §2.2 (the 2×2 by donor × guide at ~50 cells), §4 (the count-thinning curve to correct), and its
bank (`work/bank280_X.npy` 17,495 × 18,129 f32 = every (guide, donor, condition) pseudobulk of the 297 panel targets +
every non-targeting pseudobulk; `work/obs280.npz` keys `donor, cond, gtype, run, pert, guide, n_cells, total, keep_*`);
268's `work/halves/<line>.npz` and `build_halves*.py` (Broad-16 guide / cell halves on 12 working lines, 90 targets);
`scripts/build_donor_table.py` (the response object).

## The sources (inventory from metadata first; then freeze the plan; then measure)
| source | on disk | units available |
|---|---|---|
| H1 hESC, 10x Flex (`data/harness/h1-test.h5ad`, 85 targets × 400 cells, 48 lanes, 5 two-construct targets) | yes | same-guide cell halves at n ≤ 200; lanes; guide (5) |
| K562-GWPS (`data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad`, 8.8 GB, streams; gemgroups, guide ids; 269 panel targets, 23–1,745 cells, median 175; 8 two-construct targets) | yes | cell halves at n ≤ ~90 (deep targets more); gemgroups; guide (8) |
| K562 essential screen (`data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad`, 1.5 GB) — the SAME line, an independent experiment (different library, different run) | yes | **different study, same line**, on the targets shared with GWPS (≈ 2,000 essential genes, none on the panel — state the overlap) |
| Orion HCT116 / HEK293T (`data/external/orion/data/*.parquet`, 332 files, 118 GB; 267's inventory) | yes; 267's half-banks cached | batch halves, random halves (cached); a raw pass only if the cell ladder needs it (267: 11–13 GiB, 17–22 min per line) |
| Nadig 2025 Jurkat, HepG2 (`data/external/nadig2025/*_raw_singlecell.h5ad`, 14 GB; obs `gem_group`, `sgID_AB`, `gene`) | yes | cell halves; gemgroups; guides per target (inventory) |
| Broad-16, 12 working lines (268's halves) | yes | guide halves, cell halves, 90 targets; lines |
| GWCD4i CD4 T (280's bank: (guide, donor, condition) pseudobulks) | bank only | guide × donor × state; UMI thinning; **no cells** → the bounded extraction below |

## The bounded CD4 cell-level extraction (frozen plan before any byte; the review's "scope a bounded extraction")
Probe ONE object with 280's `remote_h5.py`: `D2_Stim48hr.assigned_guide.h5ad` (137.6 GB; Stim48hr is the run-clean
condition, 280 §1.2). Read `obs` (`perturbed_gene_name`, `guide_id`, `lane_id`, `total_counts`, `low_quality`) and the
layout of `X` (encoding attrs; CSR `indptr` if present; chunking / compression). **Proceed only if `X` is uncompressed
CSR with range-addressable rows** (as the pseudobulk object was): plan the rows of every cell whose target is a panel
target plus a non-targeting sample of ≤ 2× those cells, estimate the bytes from `indptr`, and extract them if the
estimate is ≤ 25 GB for the donor (parts deleted as read; the kept bank ≤ 3 GB with sha256). A second donor
(`D3_Stim48hr`, 137.8 GB) only if the first came in under 15 GB and ≥ 2 hours remain. Otherwise stop at the plan (sizes,
layout, the exact reason) and use the bank-only CD4 units. Record bytes transferred, requests, sha256 of what is kept.

## Frozen statistics (write `PREREGISTRATION.md`, hash it into `work/PREREGISTRATION.sha256` before measuring)
- **Response object:** 267's / 273's residual space — per view, `log1p(50,000 · S/ΣS)` minus the view's OWN controls,
  centred on the view's mean over the set's targets (159's panel-centred convention), on the shared gene axis with the
  target columns zeroed. The **units** (each a separate row, never pooled across sources):
  (a) **same guide, random cell halves** at MATCHED cells per side n ∈ {25, 50, 100, 200, 400} where available
  (without replacement, 8 draws averaged) — the term removable by more cells; (a′) **same cells, thinned counts** at
  fractions {1, ½, ¼, ⅛} at fixed n (280's binomial thinning) — the term removable by depth, reported apart from (a);
  (b) **different guides, same target, same line/donor**, at matched n; (c) **different donors, same guide** (CD4);
  (d) **different batches / gemgroups / lanes** (267's parity, 273's lanes) at matched n; (e) **different study, same
  line** (K562 GWPS vs essential, shared targets); (f) different lines (276's statistic, one reference row).
- **Matching:** targets binned by full-depth response magnitude (row norm decile; DE-gene count where a table exists);
  unit contrasts reported on the SAME targets where the design allows, otherwise reweighted to a common stratum
  histogram; medians and shares ≤ 0 per stratum.
- **Two readouts per cell:** the residual cosine between the views; and the scorer's own statistic — per-target `pds`
  with 273's exact recipe (view A as the prediction, view B as the truth, every target of that source as a distractor,
  midrank) plus the aggregate cosine. Say where they disagree.
- **Uncertainty:** paired target bootstrap (20,000 draws) within a source; leave-one-donor-out / leave-one-lane-out /
  leave-one-gemgroup-out sensitivity where the unit exists; report the number of cells screened (multiplicity).
- **The depth law, corrected:** the measured cell-ladder curve (a) per source against Spearman–Brown from the smallest
  n and from the largest — the replacement for 280 §4. State, per source, the n at which (a) reaches (b), (c), (d):
  beyond that n more cells of one construct cannot help.
- **Known answers (abort if any fails):** 267's batch-split medians +0.0471 / +0.0276 / +0.0282 from its cached
  `truth-<fold>.npz` at 0.0e+00; 273's h1-test +0.2429 (lane parity, 200/200) and 159's 0.330039 at 0.0e+00; 280's
  Stim48hr donor-half median +0.0539 from the bank; the union of two halves reproduces the source's full row (267 KA-1).

## Order of work (the cheap cells first; the extraction in the background)
h1-test ladders and lanes → K562-GWPS ladders, gemgroups, guides, and the GWPS-vs-essential study row → CD4 bank units
(guide / donor / state, thinning) and 280's 2×2 re-read at matched n → Broad-16 guide vs cell halves → Nadig → Orion
raw ladder only if ≥ 2 hours remain and the gate holds. Start the CD4 layout probe early; run the extraction (if
admissible) while the cheap cells compute.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the inventory; the benchmark table (rows = source × unit × n; columns =
n targets, cells/side, UMIs/cell, median residual cosine, share ≤ 0, mean `pds`, `cos_agg`) with the stratum-matched
contrasts and intervals; the depth law per source; the CD4 extraction record; `work/bench285.csv` (per source, unit, n,
target); one line: **AT MATCHED CELLS (n = …): same-guide cell split … / thinned counts … / different guide … /
different donor … / different batch … / different study …; THE DOMINANT REMOVABLE TERM AT 400 CELLS: cells / depth /
neither; DEPTH LAW: measured … vs Spearman–Brown …**. No model, no build, no submission. Then stop.
