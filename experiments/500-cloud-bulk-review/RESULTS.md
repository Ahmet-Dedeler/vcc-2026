# 500 — bulk channel (pds, mse)

No silent pds bug. The one mse formula error is real and points the expensive way, and on the exponents 448 already measured it sits inside the re-level loop's own 1% band, so it would not have changed the factors that rehearsal shipped. Expected cost of sitting on the wrong point instead of the EV point is about 1e-5 overall on a context as expansive as 448's F. Not a validation slot.

pds rank is not being taxed by a gene-order mix-up, a log-space mix-up, a double context scale, a jackknife that disagrees with cell-eval2, or ballast that fails to land the bulk profile. Those were checked by running the production functions on synthetic inputs (no `data/challenge` on this VM). Details and the negatives are below.

## 1. CONFIRMED — re-level targets ρ\* with exponent 2 after it has measured a different exponent

`experiments/448-def-rehearsal-q175/relevel_iter.py:29` sets

```text
rho_star = exp(2 * sd * Φ^{-1}(q)),  q = 0.31/2.09,  sd = 0.01  →  0.9793
```

and line 66 aims the secant at that constant: `f = f1 * (rho_star / rho1) ** (1/e)`. The same constant is hardcoded in the not-converged fallback at `scripts/def_day.sh:156`.

The hedge ntclevel states (`experiments/417-bulk-sharpen/ntclevel417.py:11-13`) is a place in ln-scale: `s = s_cliff * exp(sd * Φ^{-1}(q))`, with slopes −0.31 below the cliff and +1.78 above, per unit scale. Below the cliff `ρ ∝ f^e`, so that place is `ρ*(e) = exp(e * sd * Φ^{-1}(q))`. Exponent 2 is correct only when the render is quadratic. 448 measured it is not (`experiments/448-def-rehearsal-q175/RESULTS.md`: e = 1.89 / 1.90 / 2.30 from f=1 to lvl1, and 2.71 on F near the top). When e>2 the code sits closer to the cliff than its own hedge; that is the side with slope +1.78.

Reproduction, both on 448's own ntclevel json (factors from that RESULTS table):

```bash
uv run python experiments/448-def-rehearsal-q175/relevel_iter.py \
  experiments/448-def-rehearsal-q175/work/efr/448-lvl0.ntclevel.json:D=1,E=1,F=1 \
  experiments/448-def-rehearsal-q175/work/efr/448-lvl1.ntclevel.json:D=1.2084,E=1.2562,F=1.2871
# stock: D 1.2216, E 1.2718, F 1.2452   (e = 1.892 / 1.897 / 2.302, target 0.9793)

uv run python experiments/500-cloud-bulk-review/relevel_iter_proposed.py \
  <same two arguments>
# proposed: D 1.2223, E 1.2725, F 1.2435  (rho*(e) = 0.9804 / 0.9804 / 0.9763)
```

F's factor moves by −0.14%. On all four 448 builds the proposed script still converges to the factors that were shipped (`D=1.2218,E=1.2721,F=1.2514`), because F's ρ 0.9784 is 0.68% above `ρ*(e=2.741)=0.9718` and the loop's tolerance is 1%.

Arithmetic of the mis-centering, from `repro_bulk.py` (`work/repro.json`, key `relevel`), loss = 0.31·|u| below the cliff and 1.78·u above, u the ln-scale error, sd 0.01:

| context | e | ρ in the code | ρ*(e) | scale, code / EV | excess expected raw mse |
|---|---|---|---|---|---|
| D, f=1→lvl1 | 1.89 | 0.9793 | 0.9805 | 0.9994 | 8.7e-6 |
| E | 1.90 | 0.9793 | 0.9804 | 0.9995 | 7.2e-6 |
| F, f=1→lvl1 | 2.30 | 0.9793 | 0.9763 | 1.0014 | 4.7e-5 |
| F, local top | 2.71 | 0.9793 | 0.9721 | 1.0027 | 2.0e-4 |

Member: mse only (pds is scale-invariant). Direction: worse on an expansive context, because the bias is toward the cliff. One context's mse is 1/18 of overall and the scaled member moves about 1/(r−b) per unit raw, with r−b ≈ −0.95 (brief §2, b≈0.989, r≈0.036). 2.0e-4 raw on F is about −1e-5 overall. A synthetic render through `profile_from_effect` at q=1.75 plus cell-eval2's `_across_pert_budget` gave e = 1.79–1.81 (`work/repro.json`), same side of 2 as 448's D/E.

ntclevel's first step still uses `cliff = sqrt(claim/B)` (`ntclevel417.py:125`), the same e=2 assumption. That is what put 448's F lvl1 past the cliff (ρ 1.0568); `def_day.sh` already refuses to ship a past-cliff build. The secant is what remains in the shipped factor, and that is what the patch changes.

Proposed fix, not applied: `relevel_iter_proposed.py` (diff in `patch_relevel_rho.diff`) and `patch_def_day_fallback.py` for the `python3` heredoc at `def_day.sh:148`. The fallback keeps the inlined `Φ^{-1}` so it still runs on system python3, which is what `def_day.sh` calls.

## 2. CONFIRMED mechanism, no emitted-profile change — sharpen's restored L2 counts coordinates the renderer clips

`bulk_sharpen` (`scripts/make_submission_pm_cd4.py:4153`) restores the gated non-panel L2 of `sign(h)|h|^q`, and `profile_from_effect` (`:4246`) then emits `clip(expm1(b+h), 0) / sum`. Coordinates with `h < −b` contribute to the restored norm and contribute nothing further to the counts (the gene is already at 0).

On the production key `lodo415_lo0.5_hi3` (sha `8c82bef37d0a59b3`, energy constant c=0.6091, 32.3% of rows at the 0.5 floor, 16.3% at the cap 3), the final row scales are 0.305 / 0.609 / 1.827. Fraction of non-panel L2 sitting past the floor, median over 12 synthetic rows (`repro_clip.py` run, printed inline):

| final scale | median excess L2 | worst row |
|---|---|---|
| 0.305 (floor) | 0 | 0.34 |
| 0.609 (median) | 0.025 | 0.73 |
| 1.827 (cap) | 0.196 | 0.42 |

Read-back cosine against the designed sharpened effect at the cap: median 0.977, min 0.832 (`work/repro_followup.json`). That cosine is low because the design vector contains the unrealizable tail, not because the counts miss the realizable part.

A redistribution that moves the clipped L2 onto coordinates that still have room (`repro_clip.py` follow-up, 16 rows) changes the emitted effect by a cosine of 1.0000 median and 0.9999 minimum against the stock render, at scale 1.827. pds does not move. The read-back L2 rises about 1.4% relative (0.942 → 0.956), which the re-level removes as a global scale. Do not patch this expecting a pds gain. 417's board numbers already include this renderer.

`mass_consistent_effect` does not close the worst cap rows either (min cosine stays 0.833 at scale 1.827): the clipped coordinates have zero derivative, which that function's own comment says (`:4136`).

## Negatives (checked, not score bugs)

Reproduction is `uv run python experiments/500-cloud-bulk-review/repro_bulk.py` and `repro_followup.py`. G=4096, 400 cells, libraries lognormal around 20,109, sharpen q=1.75, effects N(0, 0.12) plus a t₃ tail. Absolute jackknives are not Cp46 (different gene count, multinomial rather than biological noise); ratios are the claim.

- **Jackknife.** `block_jackknife` (`make_submission_pm_cd4.py:888`) against `cell_eval2.moments.jackknife_correction` at TS=5e4: relative error 0 on a ballast-shaped block and on a block with a zeroed gene, 1.2e-16 on an equal-depth block. The boost's target is the scorer's C.
- **Saturation.** An un-boosted block of 384 ordinary cells plus 16 multinomial cells at 995,000 has C = 0.365× a 400-cell equal-depth draw of the same profile. After `ballast_rows`, on the floor/median/cap grid, Chat/C₄₀₀ is 1.003–1.486 (`work/repro_followup.json`). The boost is what produces saturation, and on this multinomial reference the runbook's `Chat ≥ C` assumption holds. Biological C_real is larger than multinomial C; that comparison needs the challenge bundle, which is not on this VM. 481 F2 (Cp46's 400-cell assumption) was not re-opened.
- **Ballast lands the profile it was given.** At final scale 0.305, effect-space cosine of the emitted pseudobulk against `profile_bulk` is 0.799 and a fresh multinomial draw from `profile_bulk` itself is 0.800. Median scale 0.936 vs 0.935. Cap 0.983–0.992 vs a noise floor of 0.992. Abundance cosine is 0.99998. The floor-scale effect cosine looks bad because the effect is small relative to pseudobulk sampling, and ballast matches that floor. Ratio mode (`ballast_rows:963`) is the behaviour 410 priced; the alloc comment at `:865` is stale and was not treated as a bug.
- **K sizer.** `ballast_count:958` sizes K from the first 384 libraries and from `min(profile_bulk/profile_de)` over genes above 1e-9. The argmin genes in this draw were CPM 8–40, not dust, and a CPM>5 threshold did not change K. At the cap, ratio 0 (a gene the sharpened bulk drives to 0) pins K at 16. A sweep at fixed carriers, K=1→16, moved effect cosine from 0.961 to 0.989: more ballast matches the bulk better. The DE-channel cost of K=16 is outside this review.
- **Pert-scale energy.** `apply_pert_scale` on `lodo415_lo0.5_hi3` preserves non-panel Frobenius to 1e-15 (67.090237 → 67.090237). The global constant is absorbed by the re-level; the per-row weights are the priced key.
- **Context scale, once.** `profile_from_effect` is called with `bulk_row * BULK_CTX_SCALE[str(context)]` (`:5715`). A missing context name raises KeyError. Sharpen restores L2 of the already-scaled row, so a positive scale survives as amplitude and is not applied twice. RC emitters refuse `effect_scale != 1` (`:5340`), so that knob is not a second scale on this path.
- **Bulk load discards the load-side pert scale, on purpose.** `:6587` records `apply_pert_scale`, `:6620` replaces `effects` with the dump. `def_day.sh:99-103` sets `VCC_H1_UNION_G=0.56` (and the LODO key, if any) only in the dump subshell; sharpen and context scale are applied at emission, so they are not baked twice. The counts follow the dump. The meta still describes the discarded load-side `pert_scale`. Provenance only; `max|loaded-own|` is printed. No score change while the dump is the intended matrix.
- **CD4 log base.** The builder calls `Bank.load(..., require_cd4_scale=True)` (`:2757`). The check requires the natural-log string. Not re-verified against a bank file (donors are not on this VM).
- **Zero repair.** `zero_repair` moved 32 counts on an 8.6e6-count block (`work/repro.json`). Invisible to pds/mse.
- **Existing tests.** `experiments/410-pertscale-revisit/test_ballast410.py` passes. It prints `max|log comp/bulk|` up to 0.13 and does not assert it; the assert allows a jackknife undershoot when K<2. `experiments/417-bulk-sharpen/test_sharpen417.py` was not run: it reads `data/challenge/pert_counts.csv` and `gene_names.csv`. The full `scripts/test_make_submission_pm_cd4.py` suite was not run.

## Re-run

```bash
uv run python experiments/500-cloud-bulk-review/repro_bulk.py
uv run python experiments/500-cloud-bulk-review/repro_followup.py
uv run python experiments/500-cloud-bulk-review/repro_clip.py
uv run python experiments/410-pertscale-revisit/test_ballast410.py
```

Pure Python. Nothing here depends on BSD grep, `shasum`, or the macOS `REPO` path in `def_day.sh`. `data/challenge` was not used.

## Coordinator's disposition (2026-09-27 06:50Z; brought to master from `cloud/500-bulk-review`)
Accepted as recorded; no change before October. The re-level's ρ* uses exponent 2 where 448 measured e = 1.9–2.7. The
review prices the cost at ~1e-5 overall, inside the loop's own 1 % band, and every board-read arm was levelled the same
way, so changing it now would move a calibrated constant for no readable gain. No silent pds or mse bug was found.
