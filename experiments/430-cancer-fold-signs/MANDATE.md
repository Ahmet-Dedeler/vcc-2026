# 430: which composition gives the best DE-gene SIGNS on held-out CANCER lines? Mandate from the PI, 2026-09-23 06:20Z

**Read first:**
- `docs/FINDINGS.md` §0–§2.
- `experiments/427-confident-magnitude/RESULTS.md` §427e. On H1, the unwhitened equal-weight mean of six non-H1
  donors has better sign accuracy on the reference's DE gate than our delivered LFC: 0.578 vs 0.561. Blended with ours
  and given 426's MDE add, it would buy nmae −0.0245 analytically, against −0.0098 for the plain MDE.
- `experiments/417-bulk-sharpen/BOARD-RESULT-417c.md`. The board read unwhitened (γ 0) carriers as WORSE on nmae
  (+0.0089) while H1 said better. H1 is hESC; A/B/C are cancer lines.
- `experiments/422-carrier-sign-only/RESULTS.md`, `experiments/416-de-carrier-menu/RESULTS.md`.
- 414 / 418 / 420 (PLOLO folds, the composition machinery: `experiments/418-bulk-composition/`,
  `experiments/420-bulk-stack/lib420.py`, 415's `sens415.compose_pp`).

**Question.** Is "unwhitened donor-mean signs beat our whitened composition" an H1 artefact, or does it hold on held-out
cancer lines? The answer decides whether a sign change on the carriers is worth a board slot.

**Territory.** Only `experiments/430-cancer-fold-signs/`.
- No git, no `vcc`, no `cell-eval2`, no builder or `src/` edits, no H1 emission runs (429 owns the H1 harness today).
- Analysis only, ≤ 8 GB RSS, behind the memory gate (`memory_pressure` ≥ 40 % free and free+inactive ≥ 20 GiB). A
  ~20 GB production build may be running.

**Tasks.**
- **(a) Folds.** Use the PLOLO cancer folds the repo uses: held-out line L ∈ {`orion_hct116`, `orion_hek293t`,
  `k562_gwps`}, plus kolf2 as a non-cancer reference. Take the panel targets (A/B/C panel, `data/challenge/pert_counts.csv`)
  that L covers.
  - The truth is L's own donor δ row. It is pseudobulk, and there are no p-values, so define a **gate proxy** per target:
    - (i) the top-K genes by |δ_L|, K ∈ {50, 200, 1000};
    - (ii) |δ_L| ≥ τ, τ ∈ {0.1, 0.2}.
  - Also build an L1 proxy of nmae on each gate: Σ|l̂ − δ_L| / Σ|δ_L|, with l̂ the candidate scaled to its own
    L1-optimal global k (report k).
- **(b) Candidates for the sign/lfc**, every one built WITHOUT L:
  1. the production-style composition for that fold: the 418/420 "PP" composition, whitened, i.e. the same
     construction our carriers use, with L held out;
  2. its γ 0 (unwhitened) variant, if the machinery allows;
  3. the equal-weight mean of all other donor tables (`data/donors/*`: k562_gwps, k562_essential, rpe1, nadig_jurkat,
     orion_hct116, orion_hek293t, minus L; never `h1_2025`);
  4. the CLEAN4 per-gene prior sign (`experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz`, recomputed
     without L's bank if L is one of its banks);
  5. blends 0.5·(1) + 0.5·w·(3), w ∈ {1, 1.5, 2}, each row scaled to (1)'s row norm first.
- **(c) Report**, per fold and averaged over the three cancer folds (bootstrap CI over targets):
  - sign accuracy on each gate proxy, pooled and per-target mean;
  - the L1 nmae proxy.

  Include a known-answer row: candidate (1) must reproduce 420's stored fold quantities for the same composition where
  they overlap (e.g. its per-target cosine with L's row), or say exactly why it cannot.

**Deliverable:** `experiments/430-cancer-fold-signs/RESULTS.md`:
- the known answer;
- the table;
- a direct answer: does the unwhitened donor mean (or a blend) beat the whitened composition on cancer-line DE-gene
  signs, and by how much relative to H1's +0.017;
- what this predicts for a carrier-sign arm on the board, in words and with a sign.
