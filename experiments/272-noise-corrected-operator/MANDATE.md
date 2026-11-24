# 272 — The whitening operator has learnt its own training noise (R271 §1.2). Test a NOISE-CORRECTED operator: damp the eigen-directions the donors' sampling noise dominates, and cross-fit as the check. Measured on PLOLO, priced. ≤ 6 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/272-noise-corrected-operator/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/`, `data/donors/` or `data/`. `uv run python`. **Heavy
steps are serial and gated**: an operator refit is an 18k × 18k eigendecomposition (`scripts/build_union_whitener.py`
inlines `Whitener.fit`; ~8–12 GiB, 10–30 min); run at most ONE at a time, only when `memory_pressure` ≥ 40 % free AND
`vm_stat` free+inactive ≥ 20 GiB AND no other builder / `vcc prep` (check `ps`); everything else under ~8 GB. **Time box:
write up after ~6 hours**; do the cheap construction on all three folds BEFORE any refit.

## Why (R271 §1.2, R267 §3, R262 §4.2)
Each PLOLO fold's operator (`experiments/062-cd4-post-whiten/cloud-output/whiteners/held_<line>.npz`) is fitted on the
two donors' FULL-depth rows. 271 whitened the same donors' batch-parity HALF rows through it and found them
**anti-correlated: median cos −0.13 … −0.24, 97–100 % of targets ≤ 0**, against a deranged null of 0.00 ± 0.01, while
the HELD line's half rows through the same operator read +0.03 … +0.07. The operator's amplified directions are the
ones in which each fitted row's two halves cancel — it has learnt the sampling-noise realisation of its training rows
(the donor rows are ~5–9 % signal energy at their depth, R271 §1.1). The production operator `union_pruned` is fitted
the same way on the three production donors. **Hypothesis: damping the eigen-directions that the donors' measured
sampling noise dominates raises the transfer cosine on every target (the systematic component), not per target;
on the marginal band +0.01 cosine buys +0.095 `pds` (R262 §4.2).** The operator-hyperparameter family (γ, λ, exponent,
rows: 193 / 212 / 226, READ-FIRST §5) moved the board by ≤ +0.0004 and is closed — this is a different construction
(a measured noise covariance, not a knob), and it must be priced as such.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; R271 §1.1–1.2 and `diagwht271.py` (the anti-correlation measurement — reproduce it
as your KA-0), `halfrows271.py` (the donors' half-banks, from 267's `work/half/<line>/bank-*`); R267 §1 (how halves were
built); `scripts/build_union_whitener.py` (the fit: covariance over which rows, the γ = 0.65 suppression, the pruning,
`axis_sha256`, the RSS guard) and `experiments/162-retrodicting-validator/plolo.py` (`Operator.apply`, `HELD_OPS`, the
fold construction); `experiments/170-*/RESULTS.md` (how γ was chosen: an end-to-end V with its vertex at 0.65 — the
noise correction must be compared against that tuned baseline, not against an untuned operator); READ-FIRST §3
(pricing), §5 (operator family closures); 257's `price257.py` / 261's `price261.py`.

## Constructions (frozen in `op272.py` before measuring)
A. **Diagonal Wiener damping in the FITTED eigenbasis (cheap; all three folds first).** Let the fold operator have
   eigenvectors V and eigenvalues λ (of the covariance it was fitted on; recover them from the .npz or recompute the
   projections). Estimate the donors' sampling-noise covariance from the half differences: for each donor row t,
   d_t = (row_A − row_B)/2 (a noise realisation at full depth up to the half-vs-full factor — state the factor you use
   and why; 271 measured the half-bank norms at 1.37–1.39× full), N = mean_t d_t d_tᵀ over the donors' carried rows.
   Per direction i: n_i = v_iᵀ N v_i, s_i = clip((λ_i − n_i) / λ_i, 0, 1). The damped operator multiplies direction i's
   weight by s_i^κ, κ ∈ {0.5, 1, 2} (κ = 0 is the incumbent operator — the KA). Report the spectrum of s_i (how many
   directions are noise-dominated, how much whitened energy they carry, and whether they are the top or the bottom
   of the spectrum — 271 §1.2 says the amplified ones).
B. **Cross-fitted operator (one fold, the check).** Refit the fold operator on half-A rows only (same recipe, same γ,
   same pruning), apply to half-B rows and to the held line; the donors' half-B rows through it must NOT be
   anti-correlated (this is the test that the anti-correlation is in-sample noise learning, not a property of
   whitening). Then the symmetric fit and the average of the two operators as the "cross-fitted" operator for
   predictions. Do this on `orion_hct116` (the fold whose donors have the lowest reliability) unless the memory gate
   refuses; state the wall clock and peak RSS.
C. **The production operator, diagnostic only:** the same s_i spectrum for `data/donors/whitener_union_pruned.npz`
   from the THREE donors' half-banks (267 built all three). No production change; report what fraction of its
   whitened energy sits in noise-dominated directions.

## Known-answer rows
KA-0: 271's anti-correlation numbers reproduce (median cos(WA, WB) per donor per fold to 3 decimals). KA-1: κ = 0
reproduces `cd4svd5`'s fold rows at 0.0e+00 (from 271/263's caches). KA-2 (construction B): the standard full-row refit
on the chosen fold reproduces the cached held operator (eigenvalues to 1e-6 relative, `axis_sha256` equal) BEFORE any
half-row fit is trusted. KA-3: the price recipe retrodicts BR249 (+0.0043 ± 0.0008).

## Measurement (three cancer folds for A; the one fold for B; `kolf2` apart)
Through `cd4svd5`'s full construction (CD4 deflated k = 5 re-solved, union at 0.28) with the modified operator in place
of the fold's: `cos_agg`, macro, `pds`, per-target `pds`, by 262 group and reference-reliability tercile, delivered
amplitude, the `nmae` identity; paired target bootstraps; leave-top-3-out retention. **Controls:** damping by a SHUFFLED
noise spectrum (s_i permuted), damping by the same amount on the OPPOSITE end of the spectrum, and the tuned-γ
baseline (the incumbent) — the arm must beat all three. Price on the coherent-pessimistic / point / optimistic rule
(the operator class is unmeasured on the board beyond the ≤ +0.0004 hyperparameter points — say so; use the re-shape
band as pessimistic for a cosine LOSS and the new-source band for a GAIN). CLEARED = pessimistic ≥ +0.0019 over
`cd4svd5` AND `pds` > 0 on 3/3 folds with the interval above 0 AND retention ≥ 50 % AND every control below half.

## Deliverable
`RESULTS.md`: KA rows; the s_i spectra (A on three folds, C on production); B's anti-correlation test and its fold
result; the measurement and price tables with controls; one line **CLEARS / DOES NOT CLEAR**; if it clears,
`PREREGISTRATION-DRAFT.md` in 249's form with the production construction (a damped `union_pruned` — the file it would
replace, `data/donors/whitener_union_pruned.npz`, and the registry pin it would need; NOT made) and the five falsifiers.
A negative — "the operator did learn its noise, and un-learning it does not move the transfer" — is a valid result
that closes this construction and leaves the finding. No model name, no submission. Then stop.
