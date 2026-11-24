# 184 — Tahoe-100M as the cross-line transfer substrate

You are a research subagent. Work ONLY in `experiments/184-tahoe-transfer/`.
No `vcc`, no `cell-eval2`, no submission, no writes outside your directory (downloads go to
`data/external/tahoe/` — that is the ONE exception, and nothing else under `data/` may be touched).
Peak RSS ≤ 12 GB. Disk budget for downloads: **≤ 120 GB**, and stop and report if you need more.
`uv run python` always.

## READ FIRST
`AGENTS.md` §5.3 (data inventory), §11.1 (method rules 1–13). Then:
- `experiments/159-transfer-ceiling/RESULTS.md` — the numbers this experiment exists to attack.
- `experiments/166-donor-acquisition/RESULTS.md` — why no more CRISPRi lines exist.
- `experiments/161-crossline-residual/RESULTS.md`, `160-control-subspace`, `174-coexpression`,
  `175-donor-pathway-match`, `104-context-heldout-credit-adapter` — the five failed attempts to
  predict the line-specific residual. Know exactly what each tried so you do not repeat it.

## The problem, stated precisely

Production predicts a **context-independent** effect matrix: for target `t` it averages the
whitened effect of `t` over three donor lines (`k562_gwps`, `orion_hct116`, `orion_hek293t`) and
applies the same vector to challenge contexts A, B and C.

159 measured, on real data:
- only **~30% (κ² ≈ 0.34)** of a target's target-specific response energy is **shared** across cell
  lines; ~70% is line-specific;
- production sits at **31.1% (HCT116) / 28.5% (HEK293T)** of the *infinite-donor* context-independent
  ceiling, and the ceiling is approached as **1/k** in the number of independent lines — half of it
  needs ~10 lines against our 3.

166 then proved **no further panel-covering CRISPRi line exists in public data**. So the
context-independent route is data-capped, and the line-specific 70% is where the score is.

Five experiments tried to predict that 70% from the recipient's basal profile and all failed. **Every
one of them had at most 3–4 cell lines to fit the transfer on.** That is a sample-size regime in which
a null is uninformative about the mechanism.

## The new information: Tahoe-100M

`https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M` — ~100M cells,
**1,100+ drugs × ~50 cancer cell lines**, single-cell RNA-seq, released by Arc. **This project has
never touched it** (one passing mention in `AGENTS.md`; zero experiments).

It is chemical, not genetic, so it cannot serve as a donor for the 300 panel targets. **That is not
what it is for.** It is the only dataset in existence that measures **the same perturbation in ~50
different cell lines**, which is exactly the substrate needed to estimate a **cross-line transfer
operator**: how a perturbation response in line `L1` maps into line `L2` as a function of the two
lines' basal profiles. Our CRISPRi data cannot estimate that operator; Tahoe can.

Note also that a competitor is already mining it: `PertMind` (team `luka`, arXiv 2608.16419, RL over
Tahoe-100M) holds **the board's best raw `reach`, 0.2908**.

## What to do, in order. Stop at the first gate that fails and write it up.

### Phase 0 — known-answer row (mandatory, rule 10)
Before touching Tahoe, reproduce **one** number from 159 on our own CRISPRi donors using their
frozen artifacts: the same-target cross-line per-target cosine **+0.1881** (HCT116↔HEK293T, 156's
frozen exact-space truth) or the corrected all-cells **+0.0317 same / +0.0006 different**. If you
cannot reproduce a 159 number, abort and say so — your instrument is not calibrated.

### Phase 1 — acquisition (cheap first)
Find the smallest sufficient representation. You need a **(perturbation × cell line × gene)
pseudobulk tensor plus per-line basal (DMSO/vehicle) profiles**, not 100M cells. Check, in order:
Arc's repo for precomputed pseudobulk or metadata; the HF dataset (`hf datasets info ...`) for a
pseudobulk/summary config; the plate-level parquet shards. Report exact file list and sizes BEFORE
downloading more than 20 GB. If only cell-level data exists, stream and aggregate — never
materialise the matrix.

Deliverable: `tensor.npz` with drugs, lines, genes, the pseudobulk log-CPM tensor, cell counts per
cell, and the per-line control profile. Record shapes, a sha256, and how many (drug, line) cells are
populated.

### Phase 2 — THE FEASIBILITY MEASUREMENT (this is the real deliverable)
**Inside Tahoe alone**, with ~50 lines, answer: *is a perturbation's response in a held-out line
predictable beyond the mean of the other lines?*

Protocol — leave-one-LINE-out over all lines with adequate coverage:
1. Baseline `M`: predict the held line's response to drug `d` as the **mean** response to `d` over
   the other lines (this is production's estimator, transplanted).
2. Arm `T`: predict it as a **basal-conditioned** transform of the other lines' responses. Start with
   the simplest defensible family and say why: e.g. a per-gene multiplicative gain that is a linear
   function of the recipient's basal log-expression of that gene; and a low-rank bilinear operator
   `Δ_L2 ≈ Δ_L1 W(b_L1, b_L2)`. Do not build a deep model in this phase.
3. Endpoints, per (drug, held line), macro-averaged over lines: **per-gene sign accuracy** on the
   held line's own significant genes, and **cosine** to the held line's true response.
4. Controls that must run: a **line-label permutation null** (shuffle which basal profile is paired
   with which held response) and a **drug-label permutation null**. An arm that does not beat both
   nulls is not a result.
5. Report the **learning curve in the number of training lines** (`k` = 3, 5, 10, 20, 40). This is
   the number that decides the whole direction: if the arm's advantage over `M` is flat in `k`, the
   five CRISPRi failures were not sample-size limited and this route is closed for good — a valuable
   terminal negative. If it grows in `k`, we have found the missing information.

**Gate:** arm `T` beats `M` by ≥ **+0.01 absolute sign accuracy** on held lines, macro, beating both
nulls. Below that, write the negative and stop — do NOT proceed to Phase 3.

### Phase 3 — only if Phase 2 passes: transfer to CRISPRi
Apply the operator learned on Tahoe to our CRISPRi donor effects and score cross-line transfer on
surfaces with truth: hold out `orion_hct116`, `orion_hek293t`, and `h1-test` (donors
`k562_gwps` + `orion_hek293t` only for h1-test; `data/donors/h1_2025` covers 85/85 h1-test targets
and is **pure leakage** — it must never enter an h1-test screen). Endpoints: per-gene sign accuracy
`q` on the reference's DE pool (incumbent **0.537**) and raw PDS. Report both, with a permutation
null. Do not build a submission.

## Reporting
`RESULTS.md` with: the known-answer row, the acquisition inventory, the Phase-2 table including the
learning curve in `k`, the nulls, and a plain verdict. **A well-supported negative is a full
result** — and given how much of the project's remaining plan rests on the line-specific residual
being unreachable, a clean 50-line null is worth as much as a win.

## Territory
185 owns the emission law and the DE call set. 186 owns competitive intel and literature. 187 owns
amplitude/`mse`/`nmae`. You own Tahoe and cross-line transfer. Do not touch production code.
