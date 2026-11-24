# 313 — Implementation parity: does production's per-target weighting compute the SAME object 308 measured? Shared-fixture test over PerturbMatch, gene masks, missing donors and near-zero rows.

You are a research subagent. Work in `experiments/313-normrel-parity/`. **Edit nothing outside it.** No `vcc`, no
`vcc submit`, no `cell-eval2`, **no build**, nothing under `submissions/` or `data/`, no download. `uv run python`.
**Peak RSS under 4 GB — a heavy 20 GB build is running concurrently, so stay light and do NOT load the production
donor banks.** **Time box: 2.5 hours.**

## Why this is the priority
Entry `Q519I2pLWRCW7PAEx5ul` scored **0.21771888483638538** (Δ **+0.0016245** over `cd4svd5`) carrying
`VCC_NORMREL_EXP=-1.5`. An external review (GPT 6 Astra) then found — and the coordinator verified — that
**production is not implementing the object 308 tested**:
- **308 / `experiments/263-donor-selection-arbiters/lib263.py` (~line 76):** PerturbMatch is blended into the K562
  donor row **before** the per-target weighting and the norm match.
- **Production (`scripts/make_submission_pm_cd4.py`, `apply_perturbmatch_k562_replacement`, ~line 1550):** the rows
  are weighted and norm-matched **first** (`build_v9_pre_whiten_effects`, which the coordinator patched), and the
  PerturbMatch correction is then applied using `k562_share` computed from the **uniform** weights, because the
  patch changed only the pre-whiten function.

**These operations do not commute.** Until they are shown to agree, the board point validates the shipped
implementation and NOT 308's local calibration — which may be part of why the board's `pds` came in 2.13× the local
raw prediction while `nmae` came in with the wrong sign.

## What to do
1. **Read both implementations** and state the two orderings precisely, in one paragraph each, with line
   references: 308's (`experiments/308-kd-weighted-donors/lib308.py`, `sweep308.py`; and 263's `lib263.py` for the
   PerturbMatch blend and the `:nm` norm match) and production's (`build_v9_pre_whiten_effects` and
   `_normrel_pre_whiten_effects` as patched, then `apply_perturbmatch_k562_replacement`).
2. **Build ONE shared synthetic fixture** — small (tens of genes, tens of targets, 3 donors), fully deterministic,
   saved to `work/fixture313.npz` with its seed — and run BOTH implementations on it. The fixture must exercise:
   (a) PerturbMatch present and absent; (b) a gene mask where a donor covers only part of the axis;
   (c) a target missing from one donor, and from two; (d) a near-zero donor row (norm ≈ 0, the division guard);
   (e) a donor row far above its own panel median (the case the arm is supposed to act on).
3. **Report the divergence**, per case: max |difference|, the relative row-norm difference, and whether the
   per-target weights themselves agree. Say plainly whether they are the same object.
4. **If they differ, characterise it**: is production's ordering a *different but defensible* estimator, or a bug?
   Which one does 308's +0.00259 [+0.00000, +0.00527] actually describe? Does the difference plausibly account for
   an `nmae` sign flip? Answer from the fixture, not from argument.
5. **Do NOT change production.** Today's entry is the standing result and is preserved. If a parity-correct variant
   is warranted, specify it as a SEPARATE arm for a future slot with its own env value — never as a silent fix.

## Known-answer rows (before any conclusion)
- With the exponent set to 0 both implementations must reproduce the **uniform** result exactly (max |diff| 0).
- With PerturbMatch absent and all donors present, state whether the two orderings coincide analytically; prove it
  on the fixture to ≤ 1e-12 or show they do not.
- Your reading of production's `k562_share` must be demonstrated, not asserted: show the value it takes under the
  patched path for a target where the normrel weights are far from uniform.

## Deliverable
`RESULTS.md`: the two orderings; the fixture; the per-case divergence table; the verdict sentence —
**"production's per-target weighting IS / IS NOT the object 308 measured (max |diff| …; the difference is
…)"** — and, if they differ, the exact specification of the parity-correct arm (env name, value, the one function
to change) WITHOUT applying it. Plus `work/fixture313.npz` and `work/parity313.json`. Then stop.
