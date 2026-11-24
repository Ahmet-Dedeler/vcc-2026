# 298 — The reviewer-facing checks the audits named, on SAVED outputs only: (1) the single-file jackknife of claim C5; (2) 290's hallmark-only ridge against the random-matched-set baseline; (3) F14's study-vs-construct contrast reported in both statistics (reliability ratio and raw correlation) on the same K562 pair; (4) a seed scan of the 2025 library's failing spacers against promoters (T7). ≤ 2.5 hours. No new data access, no build, no submission.

You are a research subagent. Work in `experiments/298-reviewer-checks/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`; no download except, for item 4, a ≤ 50 MB public
promoter-sequence table if none is on disk (cite; sha256). `uv run python`. Peak RSS under ~6 GB. **Time box: 2.5 h.**
**Reserve rule:** R3 is SPENT — item 1 re-analyses 295's SAVED per-gene outcomes (`work/reserve295.csv`) and must not
open any R3 file's `X`; the same for every other item: saved tables and arrays only, no raw expression pass except the
K562-GWPS pass item 3 needs (a public bank, not a reserve).

## Why
The audits (287, 297) and the write-up (294 §4.5) name four checks a reviewer will ask for. None needs new data.
1. **C5's file-block interval [+0.153, +0.566] rests on 8 blocks and one file (Adamson) carries 78 of 174 genes.**
   Leave-one-file-out: recompute C5a (ρ(`C_held`, `A_g`)) and C5b (ρ(`gen_pc20`, `A_g`)) dropping each file in turn,
   with the gene bootstrap, from `experiments/295-reserve-r3-claim/work/reserve295.csv`; also Adamson alone and
   everything-but-Adamson. State whether both verdicts hold without Adamson.
2. **F13's ablation (297 row 1c, NEW-as-ablation) reports the DROP from removing hallmark identity, not the FIT of
   identity alone.** From 290's saved feature table (`experiments/290-remainder-response-features/work/feat290.csv`)
   and its exact fold structure (`pred290.csv` `component_fold` / 284's components), refit the 290 ridge with (i)
   the nine hallmark projections + `pe_prolif_abs` ONLY, (ii) the same with random matched gene sets (20 draws;
   reuse 290's `ctrl290.py` construction), (iii) spectral columns only (`gen_pc5`, `gen_pc20`, `gen_cos`,
   `gen_abscos`, `breadth`, `logmag`), against the same target (`remainder_half2`) with 290's component block
   bootstrap. Report (i) − (ii) with its interval — the number that settles 297's contradiction with Wang 2026's
   Hallmark contrast.
3. **T8: F14's "a different study in the same line agrees better than a second construct" quotes a reliability
   ratio (0.79); Replogle 2022's same-line cross-experiment number is a raw correlation (0.50).** On the SAME
   K562-GWPS vs K562-essential shared perturbations 285 used (its `k562_*_study_*.npy` arrays and the pass code in
   `experiments/285-reliability-benchmark/work/k562_analyse.py`), report both statistics side by side: the
   reliability ratio at n = 25 and 50 (285's), the raw per-perturbation cosine and Pearson at full depth, and the
   same two raw statistics for the construct pair on 285's 748 two-construct perturbations — so the inversion is
   stated in one common unit as well as in the ratio.
4. **T7: what is a failing construct's reproducible response?** For the 2025 library's 34 two-construct arms
   (`experiments/293-challenge-construct-ceiling/work/construct293.csv`; spacer sequences from the 2025 files' obs
   `guide_id` if they encode protospacers, else from the library metadata on disk — state the source; if no spacer
   sequence is on disk, stop item 4 at the plan and say so), run Hartman 2026's recipe as far as the data allow:
   a 9–12 bp seed match of each failing spacer (|lfc| < 1.0) against annotated promoter windows (TSS ± 500 bp from
   `data/external/broad16line/metadata/gene_locations_biomart.csv` or the HGNC/Ensembl tables on disk; genome
   sequence is NOT on disk — if a promoter-sequence table cannot be fetched under the cap, report the match logic and
   stop), and for each failing construct whether its reproducible residual response shows the predicted off-target
   gene's knockdown in its own column. Report the count of failing constructs with a candidate seed off-target.
5. Known-answer rows before each item: C5a +0.5016 / C5b +0.3709 re-read from `reserve295.csv`; 290's primary cell
   +0.317385 from its store; 285's K562 study row (+0.0605 within / +0.0479 across at n = 25; ρ 0.791); 293's
   16-target ρ_guide +0.075.

## Deliverable
`RESULTS.md`: the four items with their KA rows, tables and intervals; `work/jackknife298.csv`, `work/ablation298.json`,
`work/t8_298.json`, `work/seed298.csv` (or the plan); one line: **C5 WITHOUT ADAMSON: C5a … / C5b … (holds / fails);
HALLMARK-ONLY − RANDOM-MATCHED: … [lo, hi]; STUDY vs CONSTRUCT in raw cosine: … vs …; FAILING SPACERS WITH A SEED
OFF-TARGET: n / 34 (or: not testable on disk)**. No claim changes are made here; the coordinator amends the record.
Then stop.
