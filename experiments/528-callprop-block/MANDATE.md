# 528 — The DE call set's COMPOSITION: a call-propensity block instead of the top-CPM block (worker mandate)

Written 2026-09-29 11:00Z by the coordinator. You are a research worker (Claude Opus 5.5, max effort). Read `AGENTS.md`,
then `docs/FINDINGS.md` §0–§2, then `experiments/517-flex-programme/RESULTS.md` §7 ("Per-gene detection (dropout)" and
`callprop517.py`), `experiments/521-strength-nreal/RESULTS.md` (call-set size, §5 global K, the regimes f and cmp),
`experiments/513-reach-walk-anatomy/RESULTS.md` (walk513.py reproduces the production emitter's H1 DE members exactly),
`experiments/192-predicted-callset/RESULTS.md`, and 441's PS3 (search `experiments/441*` for "PS3": re-choosing the block
by prior confidence lost jac from power). Your directory: `experiments/528-callprop-block/` (write ONLY there).

## The question
Our DE call set is `CallSet(400, 2200, 1.9)`: the 400 genes with the largest |design row| plus the 2,200 highest-CPM genes
in the context's controls (the "CPM block", 85 % of our calls). Its SIZE is being read on the board (523: K 1,600 / 2,000,
10-06). Its COMPOSITION has not been revisited since the block was chosen. 517 found that per-gene excess dropout (the
Poisson-expected detection rate given each control cell's library, minus the observed one) is a reliable per-gene
property shared across Flex contexts (r 0.51–0.69 among A/B/C), carries no sign information, but predicts which genes the
reference calls: at matched expression the most dropout-prone third is in h1-test's DE sets 2–3× as often, and a 2,200-gene
block ranked by a fitted call propensity (log CPM + excess dropout, fitted on half the targets) captures 36.7 % of the
held-out targets' gate genes against 33.9 % for the top-CPM block (per-target mean 0.446 vs 0.390), with no loss in the
block's median |reference lfc|. At our call-set size jac ≈ overlap / n_pred, so that is worth up to ~+0.006 raw jac on H1,
IF the emitter can make the swapped-in (lower-count) genes significant; fid's k also counts sign matches on real DE genes,
which are more predictable than non-DE ones. Price it for real.

## Objective
1. **Known answers first.** Reproduce (a) 517's callprop517 capture numbers (36.7 % / 33.9 %) exactly, and (b) the production
   emitter's h1-test DE members (walk513.py or 435's production run, `stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse`)
   exactly, before reading anything new.
2. **Deliverability.** For block genes swapped in by propensity, does the K-class emitter (src/vcc2026/rc_emit_k.py, the
   CallSet path) actually make them significant under BH at α 0.05 (they are lower-count on average)? Measure the realised
   call rate of designed calls by count stratum on the real emitter, and the jac/fid consequence of calls that fail to
   materialise (441 PS3's failure mode).
3. **Price on H1 with the real production emitter, all six members**, ≥ 3 emission seeds and a target bootstrap: the full
   swap (propensity-ranked 2,200), partial swaps (replace the bottom 10 / 25 / 50 % of the CPM block by the top propensity
   genes not already in it), and the propensity key built from the recipient's OWN controls only (no H1 truth in the fit
   for the held-out targets; fit coefficients on public data or on half the H1 targets and apply to the other half). Report
   each at the board regimes 521 uses (f 1.0 / 0.5 / 0.35 and cmp) with 521's conversion to board units.
4. **Transfer checks.** The block is context-specific by construction (each context's controls). Check the propensity
   feature's stability across A/B/C (from their controls) and whether the swapped genes differ materially between contexts;
   if a public cancer-line surface with cells exists where it can be checked (521's cancer folds: K562 / HCT116 / HEK293T
   thinned truths; note their technical blob, AGENTS §6.5), price it there too as relative direction only.
5. If positive: a builder hook, default off and bit-identical when unset (e.g. `VCC_RC_CALLSET_BLOCK=cpm|prop:<artifact>:<sha>`),
   as a patch file in your dir (NOT applied) with its unit tests, a known-answer design (unset rebuilds 446bj2z's output
   `ec3a5ae6` byte for byte), the per-context A/B/C block artifacts computed from the challenge controls, the October form
   (computed on the day from D/E/F's own controls: what def_day would need), and draft board bins for a DE-only arm vs
   446bj2z (pds/mse bit-identical; the primary read on fid + jac overall, as 523 does, since the swap changes the emitter's
   pass-2 placement draws: reach carries the re-roll noise).

## Territory (disjoint from the other work)
- YOURS: which genes are in the DE call set (the block's composition), and whether the emitter can deliver them.
- NOT yours: the call set's size (523, on the board 10-06/07), sign sources (517: programmes and the uncalled-gene hook),
  the bulk channel, production code, builds, slots.

## Rules
- Memory: a production build chain (`experiments/527-arc-hook-arm/build527.sh`) runs until ~14:00Z 2026-09-29. Until its
  process has exited, keep EVERY process ≤ 4 GB RSS and check `vm_stat` free+inactive ≥ 33 GiB before any step expected above
  2 GB. After it exits: ≤ 10 GB per process, one heavy step at a time, free+inactive ≥ 23 GiB. The 517 worker also runs
  real-emitter H1 jobs (~7.5 GB each, `experiments/517-flex-programme/drive517.sh`): never run a heavy step while one of its
  `run517.py` processes is alive. Report peak RSS (`/usr/bin/time -l`).
- No builds, no `vcc`, no `cell-eval2` CLI, no git, no edits outside your directory. Never start long-lived shell commands
  whose arguments name production scripts (build chains watch process lists).
- Known-answer row in every screen; target-bootstrap CIs; negative results are results. Finish with `RESULTS.md` (numbers,
  the recommendation, and for a positive: the patch, tests, KA design, artifacts and draft bins) and a final message
  summarising the numbers.
