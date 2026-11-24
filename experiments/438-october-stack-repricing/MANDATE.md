# 438: re-price the September bulk stack on October-like panels, and make that a one-command step for Oct 22. Mandate from the PI, 2026-09-23 17:00Z

**Why.** The bulk stack we ship was priced on the A/B/C panel by 418/420's PLOLO machinery:
- the per-target key (410, clip 0.3–3);
- sharpening q 1.25 (417a; q 1.5 in 433k);
- the composition stack H1 ×2 + CD4 ×1.5 (417e);
- the per-context re-level.

On Oct 22 the panel changes:
- 300 new, essential-free targets (238, 436);
- K562-GWPS coverage between 146 and 300 of 300 depending on the draw (238: `efp` 146, `efr` 300);
- the H1 union's coverage collapses (h1_2025 20 → 3; gse295214_h1 64 → 34, DEF-RUNBOOK §5.-2).

DEF-RUNBOOK says to re-price the H1 dose "as a first-class step on Oct 22, not as an afterthought", but no tool does
that. 432 proved the stack RUNS on synthetic D/E/F. Nobody has checked whether its WEIGHTS are still right there.

**Read first:**
- `AGENTS.md` (all);
- `docs/FINDINGS.md` §0–§2;
- `experiments/420-bulk-stack/RESULTS.md` (your machinery, and its known answers) and `lib420.py`;
- `experiments/418-bulk-composition/RESULTS.md`;
- `experiments/238-essential-free-panel/RESULTS.md` (the October-like synthetic panels, `efp` / `efr`, and how they
  were drawn);
- `docs/DEF-RUNBOOK.md` §5 (coverage-drop responses);
- `experiments/432-def-rehearsal-new-knobs/RESULTS.md` §6–§8 (the day-of sequence this will plug into).

**Territory.** Only `experiments/438-october-stack-repricing/`.
- No git, no `vcc`, no `cell-eval2`, no builder runs, no edits to `scripts/`, `src/` or other experiments. Import
  their code read-only.
- ≤ 10 GB RSS per process, strictly serial, behind the memory gate (`memory_pressure` ≥ 40 % free, free+inactive
  ≥ 20 GiB).

**Tasks.**
- **(a) A panel-agnostic pricer.** Wrap 420's pricing so that the panel (the target list) is an argument, read from a
  bundle's `pert_counts.csv`, instead of `data/challenge`. Everything that depends on the panel must follow it:
  - the PLOLO folds' target set;
  - the axis (union of the production donors' genes minus the panel);
  - the donor coverage;
  - the key.
  Keep 420's arithmetic.
- **(b) Known answer.** On the real A/B/C panel your wrapper must reproduce 420's stored numbers exactly, for the stack
  rows (`PP|K3+sh`, `b:h1x2|K3+sh`, `a:cd4x1.5|K3+sh`, the stack) on all four folds. Abort if it does not.
- **(c) Re-price on `efp` and `efr`** (both essential-free; `efp` is the pessimistic-coverage draw).
  - Price on the same cancer folds, vs the October analogue of 417a:
    - the stack's two halves and the full stack;
    - the stack neighbours (H1 ×1 / ×1.5 / ×2 / ×2.5; CD4 ×1 / ×1.25 / ×1.5 / ×2);
    - sharpening q 1.25 / 1.5 on the stack;
    - the key at clip (0.3, 3) vs 420's LODO key at clip (0.5, 3).
  - Report per fold, the 3-cancer-fold mean with a bootstrap CI, and the energies (ρ vs the cliff), as 420 did.
- **(d) The answer.** Do the A/B/C-optimal weights (H1 ×2, CD4 ×1.5, q 1.25–1.5, key 0.3–3) still hold on
  October-like panels?
  - If not, which way do they move, by how much, and what does the October stack lose by keeping the A/B/C weights?
  - Say plainly whether the H1 ×2 half survives the H1 union's coverage collapse.
- **(e) The day-of command.**
  - One command that, given the day's D/E/F bundle, prints the re-priced stack table and the recommended weights.
  - Its expected runtime and RSS.
  - Where it plugs into 432 §8: which dump recipe (`VCC_H1_UNION_G`, `VCC_CD4_TARGET_G`, q) it would change.
  - Write it as text for DEF-RUNBOOK; do not edit the runbook.

**Deliverable:** `experiments/438-october-stack-repricing/RESULTS.md`:
- the known answers;
- the tables;
- the direct answer to (d);
- the day-of command.
