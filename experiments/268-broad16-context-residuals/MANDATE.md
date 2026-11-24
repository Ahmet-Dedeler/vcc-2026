# 268 — Multi-context biology on Broad-16: the reproducible context-residual map. Which part of a perturbation response is shared across lines, which part is line-specific and REPRODUCIBLE, and is the line-specific part programme-level and related to basal state? Feasibility stage, ≤ 4 hours. Measurement only.

You are a research subagent. Work in `experiments/268-broad16-context-residuals/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`, nothing written under `data/` (unzip into `work/`).
`uv run python`. Peak RSS under ~6 GB (siblings 267/269/270 run alongside). **Time box: write up after ~4 hours.**
Read `docs/RESEARCH-PROGRAMME-2026-09.md` first; its rules bind you (experimental unit, refit inside the split, strong
baselines, correct- vs wrong-recipient controls, shuffled labels, paired bootstraps).

## Why
The programme's organising question is which responses are conserved across cell states and which depend on the
recipient. The only multi-line genetic-perturbation set on disk with many lines is Broad-16 (`data/external/broad16line/`,
`data.zip` + `metadata.zip`, 78 MB): 90 non-control targets, 87 complete across 16 contexts. Experiment 111 tested ONE use
of it (a DepMap context-program operator, NO-GO, `experiments/111-depmap-context-program-gain/RESULTS.md`) and 224 found
that context-conditioned heads on our donors contributed nothing (`experiments/224-learned-transfer/RESULTS.md` §193);
neither established that a line-specific response component is absent — they showed two constructions did not
predict it. **First measure whether line-specific residuals are REPRODUCIBLE; only then ask whether they are
predictable.**

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 111's `RESULTS.md` and its scripts (how the data was parsed, what the units are);
224 §"context conditioning" (what was tried); `experiments/262-where-pds-is-lost/RESULTS.md` §5 (the anti-aligned group
— our donors' consensus points away from the held line for ~17 % of targets: the phenomenon 268 should be able to see
in a 16-line set); `docs/reference/vcc2026-metrics-BRIEF.md` (why per-perturbation direction matters).

## Stage 1 (this mandate)
1. **Inventory first, and DECLARE THE RESERVE:** parse `metadata.zip` / `data.zip` (pseudobulk or cells? guides?
   replicates? batches? per target × line); write the unit table. Then name **four lines as the programme's reserved
   confirmation lines** (pick them by a fixed rule — e.g. every fourth line alphabetically — before looking at any
   residual), exclude them from EVERY fit and selection below, and list them in `RESULTS.md` §0 so the coordinator can
   record them in `docs/CONFIRMATION-LEDGER.md`.
2. **Variance decomposition, per target, on the 12 working lines:** response = shared (cross-line mean) + line-specific
   residual + noise. Estimate the noise floor from replicate / guide / batch halves where the data has them (say which
   unit), else from a cell-level split if cells exist, else state that the floor is unavailable and use the between-
   guide difference if guides exist. Report, per target and pooled: the fraction of response variance that is shared,
   line-specific-and-reproducible (split-half cosine of the residual > the noise floor), and noise. **This table is the
   deliverable's core.**
3. **Is the line-specific part programme-level?** Unfold the (target × line) residual tensor into a (target·line) × gene
   matrix; SVD; how many components capture the REPRODUCIBLE part (compare against the split-half reproducibility of
   each component: fit components on half A, project half B). Report the effective rank and the top components' gene
   loadings in words (not a full table).
4. **Does it relate to basal state?** For each line, basal programme activities from its control profile (a small set:
   the top principal programmes of the 12 control profiles, fitted INSIDE the training lines). Leave-one-line-out:
   predict the held line's residual for each target from (a) target-only (the shared response: the strong baseline),
   (b) recipient-only, (c) target + recipient basal programmes (linear, few parameters), with the controls:
   **correct recipient vs WRONG recipient (another line's basal state) vs shuffled line labels.** Report per-target cosine
   of the predicted residual to the held line's reproducible residual, paired bootstrap. If (c) does not beat (a) and
   the wrong-recipient control, context learning is not demonstrated on this data — say so plainly.
5. **The bridge to our problem (small):** which of the 300 VCC panel targets are among the 90, and for those, is the
   line-specific-reproducible fraction larger for 262's anti-aligned group than for the rest (on our held lines)?
   Report the overlap count and the comparison if ≥ 10 targets overlap; otherwise state the overlap only.

## Deliverable
`RESULTS.md`: the unit table and the four reserved lines; the variance decomposition table; the programme-level result
(effective rank of the reproducible residual); the LOLO prediction table with the three controls; the panel bridge;
one line: **REPRODUCIBLE LINE-SPECIFIC COMPONENT: YES / NO (fraction …); PROGRAMME-LEVEL: YES / NO (rank …);
BASAL-STATE-PREDICTABLE under LOLO: YES / NO**, and what stage 2 would be (or the closure at this construction's scope).
No model name, no submission. Then stop.
