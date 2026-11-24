# 382 — Night audit: all five headline numbers reproduce (one mapping trap found and cleared)

Status: **COMPLETE**. Every item is a MATCH at its stated tolerance. Two apparent
mismatches encountered during the work were my own errors, both diagnosed and fixed
with evidence below — neither touches an audited record.

## Known-answer reproduction (own code, before any new number)

| quantity | mine | record | match |
|---|---|---:|---|
| 323 cap15 raw fid | 0.4975545182292535 | `metrics-cap15.json` same | exact (0.0) |
| 323 cap15 raw nmae | 0.9508364630325662 | 0.9508364630325663 | 1.1e-16 (1 ulp summation order) |
| context_A shape / median lib / gene head | 18400×18533 / 20109.0 / TSPAN6,TNMD,DPM1 (h5ad var == CSV) | AGENTS.md §2, `gene_names.csv` | exact |
| Jurkat own-target negativity | 2045 rows, 0.99707 | 370 mandate (~99%) | exact to 5 dp |
| K562 own-target negativity | 7617 rows, 0.99028 | 370's check 7617/0.99028 | exact |
| A-block cutoff / mapped genes | 130.43 CPM / 8,284 | 370: 130.4 / 8,284 | match |
| 345 gate rebuild | 58 eligible targets | 345 (58) | exact |

## What was measured (inputs, estimator, controls)

All code written from the spec (`docs/reference/vcc2026-metrics-BRIEF.md` §§3/4/7)
and the assignment text. Audited scripts never opened. Scripts `audit1.py`–`audit5.py`,
outputs `audit1.json`–`audit5.json`, all in this directory.

1. **363**: fid/nmae from `stats-{cap15,block_down,block_random}.npz` + 323
   `reference.npz` (universe column order assumed shared — validated by the match
   itself). fid: own q<0.05 ∩ adjudicable (ref lfc non-NaN, nonzero), target excluded,
   k = sign-match via product>0 (directionless own call = miss), 0/0 dropped.
   nmae: ref gate (q<0.05, finite ref lfc, target out), ≥10 genes, denom>0 finite,
   nonfinite own lfc → 0. Gate cross-check vs `reference_de.parquet`.
2. **346**: scaled=(raw−b)/(r−b) from 346 `metrics-*.json` + 323 `h1-anchors.json`;
   mse clamped to [0,1] after scaling (verified: cap15 raw 0.9749 → −0.0044 → 0.0);
   overall = mean of six. Compared to `results.json` full-precision contrasts and to
   the assignment's 4-decimal headlines separately.
3. **370**: A mean-CPM from `context_A.h5ad` via h5py CSR scan, 2,000-row chunks;
   block = top-2,200 A-CPM among Jurkat-mapped genes; P(down) fit on even target rows,
   held-out eval on odd rows (and swapped); eval = block minus donor-axis top-400
   |delta| head; exact zeros dropped.
4. **Board**: deltas from `score_avg` in `data/board-entries/`; each `score_avg`
   checked against the mean of its six `score_*`. Entry identity verified by
   file note, not filename assumption (see trap below).
5. **345**: K562 all-row mean (`delta.npy` mmap, 9,787 rows) mapped by gene name to
   the H1 universe; |lfc|-weighted sign accuracy per target on the parquet gate
   (58 eligible), macro + pooled. Sign convention verified first (own-target 0.99028).

## Results

| # | item | recomputed | recorded | tol | verdict |
|---|---|---|---|---|---|
| 1 | 363 cap15 fid / nmae | 0.49755452 / 0.95083646 | same / same+1.1e-16 | 1e-9 | MATCH |
| 1 | 363 block_down fid / nmae | 0.49391986 / 0.95476594 | exact / +1.1e-16 | 1e-9 | MATCH |
| 1 | 363 block_random fid / nmae | 0.49139852 / 0.95885208 | exact / exact | 1e-9 | MATCH |
| 2 | 346 same_ridge_b05 − cap15 | +0.00707185 | contrast +0.00707185 | 1e-6 | MATCH |
| 2 | 346 same_ridge_b1 − cap15 | −0.01494242 | contrast −0.01494242 | 1e-6 | MATCH |
| 2 | 346 diag_cross_b1 − cap15 | −0.01226534 | contrast −0.01226534 | 1e-6 | MATCH |
| 3 | 370 prior / all-down (even/odd mean) | 0.57655 / 0.51457 | 0.5766 / 0.5147 | ±0.003 | MATCH |
| 4 | 319p−319 / 319pA−319p / 254−cd4svd5 / 319d3−319 | +0.00115924 / −0.00127936 / +0.00133253 / −0.00026860 | +0.0011593 / −0.0012794 / +0.0013325 / −0.000269 | headline rounding | MATCH |
| 4 | score_avg == mean of six (6 entries) | max diff 2.8e-17 | — | 1e-12 | MATCH |
| 5 | K562 all-mean weighted acc, macro / pooled | 0.5267 / 0.5474 | PI 0.527 / 0.547; 345 said 0.275 | — | MATCH PI; 345 refuted |

Units: raw metric units (items 1, 5: fractions/raw NMAE), H1-local scaled overall
points (item 2), board overall points (item 4), fractions (item 3).

Notes: pds/mse cannot be recomputed from 363's saved files (no design effects, no
emitted counts) — stated, not attempted. Item 1 nmae diffs of 1.1e-16 are the known
1-ulp `nanmean` summation-order effect (same as 346's known-answer record). All 85
targets fid-defined; 58 nmae-eligible in every arm.

## Discrepancies encountered (both mine, both resolved)

1. **Item 4 mapping trap.** I first mapped "319" to `Q519I2pLWRCW7PAEx5ul.json`
   (score 0.2177189) and got 319p−319 = +0.00210 and 319d3−319 = +0.00067 (wrong sign).
   That file's note says "312 normrel exp −1.5" — it is the 312 arm, not 319. The true
   319 entry is `xuDelEik7zGWMXdNu0R4.json` ("319 cap variant, suppression-only",
   0.21865973944346428). With the correct file all four deltas match headlines to
   ≤4.1e-07 (rounding). Lesson for future audits: identify board entries by note
   content, never by filename guess.
2. **Item 3 head-definition bug.** My first version took the top-400 |delta| *within
   the block* (always removing exactly 400, n_eval 2.04M) instead of the donor-axis
   top-400 intersected with the block (n_eval 2.39M, matching 370's 2,393,189). That
   gave prior 0.5646 / all-down 0.5106 — outside tolerance. After the fix: 0.57655 /
   0.51457, matching to ~1e-4. The n_eval mismatch (2.04M vs 2.39M) was the diagnostic.

## Verdict against the PLAN.md kill criterion

No audited record fails. 363's fid/nmae, 346's scaled Δs, 370's prior accuracies, the
four board deltas (with 319 = `xuDelEik7zGWMXdNu0R4`), the score_avg==mean-of-six
identity, and the PI's 0.527/0.547 withdrawal values all reproduce from raw files
with independent code. 345's 0.275 does not reproduce (0.5267 macro on covered gate
genes) — the withdrawal stands.

## Limits

- Item 1 validates 363's arithmetic from saved stats, not the stats themselves
  (emission/Wilcoxon upstream); the 363 known-answer record already ties those to 323.
- Item 3 is one deterministic even/odd split vs 370's 5-random-split mean; agreement
  to 1e-4 is better than the ±0.003 gate and better than split noise needs to be.
- Item 5 evaluates covered gate genes only (62.7% of gate mass by gene count is in
  the K562 bank); counting absent genes as wrong would lower the number — the PI's
  0.527 uses the same covered-gene convention (exact agreement confirms it).
- H1-panel numbers are development evidence, not A/B/C forecasts.

## Exact commands to reproduce

```bash
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/382-night-audit/audit1.py
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/382-night-audit/audit2.py
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/382-night-audit/audit3.py
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/382-night-audit/audit4.py
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/382-night-audit/audit5.py
```

Audited scripts were not read before numbers were final and are not needed to
reproduce this audit.

## Peak memory and wall time

`memory_pressure` read 80–81% free before the heavy runs; cap was 12 GiB RSS.
Largest live objects by construction: K562 delta mmap (323 MB, paged), context_A
CSR chunk scan (~100 MB live), stats/reference arrays (<100 MB) — peak RSS stayed
an order of magnitude under cap (not instrumented sample-by-sample). Wall: audits
1/2/4/5 in seconds-to-~2-min; audit3 ~3 min per run (two runs: bug + fix). Total
well inside the 2 h budget. No scratch files over 1 MB retained (`audit*.json`,
KBs each); no downloads.
