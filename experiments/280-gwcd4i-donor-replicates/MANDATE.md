# 280 — The replicate unit the programme lacks: GWCD4i's cell-level files (4 biological donors × 3 states, 10x Flex, ~540 cells per target). Feasibility and a preregistered analysis plan, with a bounded acquisition if the sizes allow. ≤ 4 hours. No build, no submission.

You are a research subagent. Work in `experiments/280-gwcd4i-donor-replicates/`. **Edit nothing outside it** and write
nothing under `data/` — any download goes to `work/` (hard cap **60 GB total on disk at any moment**; check
`df -h` before each file; delete what you finish with unless the plan below keeps it). `uv run python`, backed / streaming
reads only, peak RSS under ~10 GB. You may use web fetch for GEO / SRA / the CZI Virtual Cells Platform metadata.
**Time box: 4 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you (name the unit; preregister before looking).

## Why (R278 §4, F1, F8, R267 §3)
278 found that the panel's deep measurements are almost all in **GWCD4i** (Zhu, Dann, … Marson; GEO GSE314342, MIT
licence): 283–287 of the 300 panel targets at a median of ~540 cells per (target, state), **10x GEM-X Flex** — the
challenge's own chemistry — across **4 biological donors × {Rest, Stim8hr, Stim48hr}**. On disk we have only its DE
table and the pooled `cd4_compact` bank. The programme's reliability law (F1) rests on random / batch / guide halves;
**no fold on disk has a BIOLOGICAL replicate** (person-to-person). Two questions only this data can answer:
(Q1) how reproducible is a deep, Flex-chemistry perturbation response across human donors of the same primary cell
type — the true replicate ceiling of the challenge's regime; (Q2) does averaging four donors beat one (the value of
replicates vs cells — programme rule 7's depth/replicate axis), per target and by conservation class (276's map).
Both are scientific results; both bound what a Flex-400 reference can resolve. Neither touches the competition build.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 278's `RESULTS.md` §1 (the GWCD4i row: files, chemistry, licence, cells) and §4;
`docs/FINDINGS-2026-09.md` F1 and F8; 267's `truth267.py` (how halves and their controls become effect rows in the
instrument's residual space — reuse the object); `scripts/build_donor_table.py` (the response object; `min_cells` is a
floor, not a cap — R278 §4b); `experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad` and
`data/external/cd4t/cloud-compact/cd4_compact` (what is on disk; their provenance in `docs/DEF-RUNBOOK.md` §"CD4");
276's `work/conservation276.csv` (`C_raw`, `C_rel`, the classes).

## Stage 1 (this mandate)
1. **Inventory from metadata only:** GEO GSE314342 / SRA SRP643211 / CZI platform: the per-file list of the cell-level
   `D*_*.assigned_guide.h5ad` objects (donor × condition), sizes, formats, what obs columns they carry (donor, guide,
   condition, batch/lane), and terms. Which subset gives the four-donor replicate unit for the panel targets at the
   least GB — expected: ONE condition (Rest) × 4 donors. Write the acquisition plan with exact URLs, sizes and the
   scratch/disk budget BEFORE downloading anything.
2. **Preregister (hash a `PREREGISTRATION.md` before any cell is read):** the two questions as frozen statistics:
   (Q1) per panel target, the across-donor reliability of the response — cosine of donor-half effect rows (donors
   {1,2} vs {3,4}, all three pairings) in 267's residual space, with a cell-count-matched WITHIN-donor split as the
   comparator (so donor-to-donor disagreement is separated from cell sampling); predictions: median across-donor
   cosine, its relation to the within-donor split, and to 276's conservation class (state the direction you expect
   and the pass rule for "donor is a real unit": across-donor < within-donor at matched cells, interval excluding 0);
   (Q2) the response from the four-donor mean vs one donor vs two donors, scored by cosine to the held-out donors'
   consensus (leave-donors-out), per target: the replicate-vs-cells curve at matched total cells.
3. **Acquire the planned subset if it fits the cap** (else stop after the plan and say so). Verify checksums where
   given; record sha256 of what you keep.
4. **Compute Q1 and Q2 on the 283–287 panel targets present**, with the known-answer row: the pooled four-donor Rest
   response reproduces the on-disk `cd4_compact` Rest bank's rows for those targets (state the max |Δ| and why any
   difference exists — e.g. guide filtering, the platform's processing).
5. **Descriptive (post-hoc, labelled):** which targets are donor-variable (candidates for genotype- or state-
   dependent responses), and whether donor variability tracks 276's low-conservation tail.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the inventory and plan; what was acquired (files, GB, sha256, terms); the KA
row; the Q1 table (across-donor vs within-donor reliability, per target summary, by class) and the Q2 curve; the
descriptive tail; one line: **DONOR IS A REAL UNIT: YES / NO (across-donor median …, within-donor at matched cells …);
FOUR-DONOR MEAN BEATS ONE AT MATCHED CELLS: YES / NO (Δcos …)**. If the download does not fit the cap, the deliverable
is the plan with sizes and the coordinator decides. No model, no build, no submission. Then stop.
