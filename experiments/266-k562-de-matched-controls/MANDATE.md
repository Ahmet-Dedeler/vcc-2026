# 266 — 260's NO-GO revisited narrowly: do batch-matched controls change K562-GWPS's real Wilcoxon DE-set sizes on the panel? A multi-seed pilot on the rows 260 already extracted. ≤ 60 minutes. No harness, no scoring, no build, no submission.

You are a research subagent. Work in `experiments/266-k562-de-matched-controls/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2` run, nothing under `submissions/` or `data/harness/`. `uv run python`. Peak RSS under
~6 GB — 260's `extract_and_de.py` peaked at 8.25 GiB on float64 CPM temporaries; use float32 in place. **Time box: 60
minutes.** A sibling (265) runs light array analysis alongside.

## Why
260 (`experiments/260-raw-de-diagnostic-feasibility/RESULTS.md` §2c, §5) found that the scorer's own Wilcoxon recipe on
the actual K562-GWPS cells gives the panel targets a median of **5** significant genes and that only 31 % clear
`nmae`'s 10-gene gate — the NO-GO's whole basis. But its control pool was **globally pooled across 267 gemgroups**
(`work/extract_and_de.py` around line 38) whose per-gemgroup median UMI runs 2,854–22,130 (260 §1). Depth
heterogeneity inside the control pool inflates within-group variance and can suppress a rank test's calls, so the low
DE counts may partly be a control-pool artefact rather than the panel's biology. **Question: with controls drawn from
the SAME gemgroups as each target's perturbed cells (and depth-matched), how many significant genes do the same
targets get?** If the median moves from 5 toward the board's ~89-gene regime, 260's NO-GO is wrong; if it stays in
single digits, the NO-GO stands on a cleaner measurement.

## READ FIRST
260's `RESULTS.md` §1–§2c and `work/extract_and_de.py`, `work/measure_gwps.py`, `work/candidates_de_power.csv`,
`work/de_power_summary.json`, `work/extract_rows.npz` (the seed-0 extraction: which rows, what is cached);
`docs/reference/vcc2026-metrics-BRIEF.md` (the Wilcoxon recipe, the CPM > 5 gate, the 10-gene `nmae` gate);
`scripts/check_harness_depth.py` (the depth-ratio acceptance band); `docs/HARNESS-DEPTH-CONFOUND.md`.

## The pilot (same targets, same recipe, different controls)
1. Take 260's 50 best-covered candidate targets (`candidates_de_power.csv`). For each, the perturbed cells' gemgroup
   composition (from `obs.gemgroup`/`batch`); draw controls **from those gemgroups in proportion**, matched in count to
   260's design (and a second draw depth-matched by UMI quantile within gemgroup). Three seeds each.
2. Re-run 260's Wilcoxon recipe exactly (its gate, its BH threshold, its cell cap) on (a) the pooled controls (reproduce
   260's per-target significant-gene counts to the digit — the known-answer row), (b) gemgroup-matched controls, (c)
   gemgroup + depth-matched controls. Report per target and as medians / the fraction clearing the 10-gene gate, with
   the seed spread.
3. Diagnostics: the control pool's UMI CV under each scheme; the per-target depth ratio (must stay inside 0.95–1.05); and
   for the five best-covered targets the overlap of the DE sets across schemes (are the "new" genes the same genes at
   a looser bar, or genuinely different?).

## Deliverable
`RESULTS.md`: the KA row (pooled reproduces 260), the three-scheme table, the medians and gate fractions, and one line:
**does batch-matched control selection change 260's coverage verdict (median significant genes ≥ 20 and ≥ 60 % of
candidates clearing the 10-gene gate would reopen it)? YES / NO.** No harness, no scorer run, no model name, no
submission. Then stop.
