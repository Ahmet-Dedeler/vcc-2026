# 226 — Does the whitening operator get better with MORE FITTING ROWS? The one operator question nobody asked.

You are a research subagent. Work in `experiments/226-operator-rows/`. **No `vcc`, no submission
build, no production change, no `cell-eval2` while another agent is running one — check first.**
Peak RSS ≤ 12 GB (193 ran two 8.8 GB jobs and had to serialise them; do the same). `uv run python`.

## The question

The shipped operator (`data/donors/whitener_union.npz`, sha `e2f65ecf…`) was fitted **once**, on
**2026-08-22**, on a matrix of **(33,787 × 17,942)** — roughly **n ≈ 2p**. A covariance/whitening
operator at n ≈ 2p is barely determined, and every regularisation result we have is a way of coping
with that.

**Two experiments closed the operator's *shape* and neither touched its *sample size*:**
- **193** closed the **regulariser** — "the ridge is optimal, the geometry is CLOSED" (γ = 0.65,
  λ = 0.0).
- **212** closed the **exponent**, by a ratio argument and therefore instrument-independently.

**Nobody has asked whether fitting it on more rows helps.** And this is not a static question any
more: **branch A fired on 2026-09-04** — `h1union` scored **+0.0069, rank 3, all six members up** —
so the source stack is now *expected to grow*, and **every source we add is also more fitting data
for the operator.** 222's falsifier 1 confirms the operator's sha was **unchanged** when H1 was
added: H1 entered as a residual source and contributed **zero rows** to the operator.

**If more rows help, the effect compounds with 225's source ladder and we have been leaving it on
the floor every time we add a source.**

## READ FIRST
- `AGENTS.md` **§8.1, §8.2, §8.4** (the corrected **2.51×** multiplier; **+11.6–13.6%** amplitude-floor
  and **+2.5–2.9%** shape-only bars) and **§4.1c** (raw `mse` is the board's cosine readout).
- **`AGENTS.md`'s new warning: the PLOLO→board cosine factor is NOT identified** — 1.384 on
  `cd4f034`, 1.000 on `h1union`, matched functional. **Report in PLOLO units and name the fold set
  and the functional (`cos_agg` vs macro, which gate, which tag). Do not convert with a scalar.**
- `experiments/193-whitener-shrinkage/RESULTS.md` and `experiments/212-gamma-cosine/RESULTS.md` —
  **what is already closed. Do not re-open γ, λ, or the exponent.** Your variable is **rows**.
- `scripts/build_union_whitener.py` — the rebuild path and its provenance note.
- `docs/DEF-RUNBOOK.md` **§0.1**: the operator is a **function of the panel** and fails *silently* on
  a new one — 98.8% of coordinates paired with the wrong gene while every check passes. 212 priced
  that trap at **−0.0429 overall on pds alone.** Everything you build must survive it.

## Work

### Phase 0 — known-answer rows, before anything (mandatory)
1. **Re-fit the operator on today's exact inputs and reproduce sha `e2f65ecf…`.** `DEF-RUNBOOK` §4
   lists this as a known-answer row; if you cannot reproduce it, stop — the rebuild path is wrong
   and nothing after it counts.
2. Reproduce **193's** γ-vertex and **212's** exponent result on the shipped operator.

### Phase 1 — the row curve
Re-fit the operator at increasing row counts (subsample the existing rows: 25%, 50%, 75%, 100%) and
measure PLOLO transfer cosine at each. **This gives the slope in `n` on data we already hold, with
no new sources**, and it is the cheapest possible read on whether the operator is sample-limited.
**If the curve is flat from 50% to 100%, the operator is already saturated and this experiment ends
here — report that.** That is a complete and valuable result.

### Phase 2 — only if Phase 1's curve is still rising at 100%
Add rows from sources **not currently in the operator**: H1 (`gse295214_h1`, `h1_2025`), CD4, and
anything 223/225 have built. Re-fit and measure. Report cosine **vs number of fitting rows**.

### Phase 3 — falsifiers, all mandatory
1. **Leak audit, and it is the dangerous one here.** `data/donors/h1_2025` covers **85/85** of
   h1-test's targets. An operator fitted on H1 rows and evaluated on h1-test is **leaking**. Use
   **PLOLO donor folds** in the challenge frame; state every source↔surface overlap. `debench.
   Surface.pds()` must **never** size an arm (188 Phase 3).
2. **Row-shuffle null** — refit on the same number of rows drawn with gene labels permuted. A gain
   that survives is measuring sample size, not structure.
3. **Wrong-recipient null** — 184 died here (arm +0.0665, control +0.0667, null +0.0663).
4. **Leave-top-3-out retention** on both cosine and `pds` separately (222 §9: 14.6% vs 72.7% on the
   same arm — the two channels have opposite fragility).
5. **Panel-identity check** — assert the refitted operator's axis matches the panel it will be used
   with, and demonstrate the failure mode deliberately once so the assertion is proven to fire.

## Deliverable
`RESULTS.md`: the reproduced sha and 193/212 rows, the **cosine-vs-rows curve**, the falsifiers, and
a one-line verdict. **"The operator is saturated in `n`" is a complete result** and would close the
last untested dimension of the geometry.

If it is NOT saturated, say what fitting-row count the curve implies and **what that is worth at
§8.4's bar** — and note explicitly whether the gain compounds with 225's ladder or duplicates it.

## Territory
**221** owns cloud scoring and the local→board calibration. **224** owns learned transfer.
**225** owns sources-as-donors. **You own the operator's fitting sample, and nothing else** — do not
price sources as donors, do not fit models, do not touch the builder or `submissions/`.
