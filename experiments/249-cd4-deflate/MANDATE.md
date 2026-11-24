# 249 — BUILD `cd4svd5`: deflate the CD4 residual's top-5 target-shared directions before sizing it. DO NOT SUBMIT.

You are a build subagent. `experiments/249-cd4-deflate/PREREGISTRATION.md` is binding — read it first, do
not edit it. **NEVER `vcc submit`.** Work in `experiments/249-cd4-deflate/`; you may edit exactly one file
outside it, `scripts/make_submission_pm_cd4.py` (247 and 222 changed it today — read their RESULTS and do
not undo them), and you may build into `submissions/` under the `249-` prefix only. `uv run python`.
**Builds ~20 GiB, `vcc prep` ~25 GiB: strictly serial, never alongside another build/prep** (`ps aux |
grep -E "make_submission|vcc prep|def_pipeline"` first; 244 may be running its rehearsal — wait for it),
gate on `memory_pressure` ≥ 40% free AND `vm_stat` free+inactive ≥ 30 GiB; record peak RSS. Commit the
builder change with its known-answer hashes (`git add <file>`, never `-A`).

## READ FIRST
`docs/READ-FIRST.md`; the preregistration; `experiments/246-direction-hypotheses/RESULTS.md` §3.3–§5 and
`lib246.py` `deflate()` (the construction, verbatim); `experiments/222-h1-union-build/` (`build222.sh`,
`stage222.py`, `repro222.py`, `profile222.py`, `metadiff222.py`, `probe222.py` — reuse them);
`experiments/247-oct22-guards/RESULTS.md` §3 (the builder's current CD4 sizing paths and the byte-identity
references); `experiments/243-transfer-factor/RESULTS.md` §5.2 (measure the board cosine off the
package); `AGENTS.md` §11.1 rule 10, §14.

## The change
`VCC_CD4_DEFLATE_K` (int ≥ 1; unset = inert). When set, after `prepare_cd4_mean`/`freeze_complete_mean3`
and BEFORE `solve_cd4_fraction` (or 247's `solve_cd4_dose`), replace the CD4 residual by
`deflate(res, "svd", k)` exactly as 246 defines it: centre over carried rows (`target_condition_count >
0`), `G = R_c R_cᵀ` on the carried rows, remove the top-k eigen-directions, re-zero outside the support.
Deterministic (`eigh`), no seed. Record `effects.cd4.deflation` = {k, n_carried, energy_removed_frac,
top_k_energy_frac_of_centred, sha256 of `U_k`} and `config.cd4_deflate_k`, written only when on.
`solve_cd4_fraction` itself is not modified. **The H1 union path must not see the flag.**

## Known-answer rows (abort on failure)
1. **Flag unset is byte-inert:** staging `effects_array_sha256` `abeb0749…` (222 `stage-on028.json`, re-verified
   by 247) reproduces; then the full A/B/C package: `.h5ad` **`8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6`**
   on 12,648,807,800 B (240 row 1 / 247 §3.4). If 247's byte-identity build is still in `submissions/`
   or its work dir, hash it instead of rebuilding, and say so.
2. **The staged arm reproduces 246:** through 162's PLOLO on the builder's own staged effects (222's
   `repro222.py` pattern), `B2_cd4_svd5_resolved` on the three cancer folds — macro **−0.81%**, aggregate
   **−0.46%**, Δ`pds` **+0.01226** vs the shipped stack — to the digit; `energy_removed_frac` **0.202**;
   the CD4 frame's solved fraction moves the way 246 §3.3 reports (0.374 → 0.396 in PLOLO's frame; report
   production's).
3. **Energy-matched random control through the same builder path** (one seed): `pds` ≤ +0.001 — proves the
   builder applies a projector, not a norm cut.

## Falsifiers — report each with its number (the preregistration's 4 and 5, plus 222's build gates)
1. Rule-10 audit diff against `submissions/222-h1union-g028.h5ad.meta.json`: exactly `cd4_deflate_k` added,
   the pre-declared CD4-frame consequences, `operator.sha256` unchanged, `alpha` 0.982, `cd4_target_f` 0.34,
   `h1_union_g` 0.28, every `effects.h1_union.*` field unchanged. Print every difference.
2. Per-gene panel-total digest split by CD4-carried row (profile222's scale-corrected form): carried rows
   move in direction; non-carried rows differ by one shared scalar only.
3. Non-vacuity: `‖ΔE‖/‖E‖` against `h1union` per block, zero-block count = 0.
4. Challenge `n_pred` within ±2% of `h1union`'s (probe222 pattern).
5. Delivered amplitude: stream the package's `A_c`/`B_c` (243 `totals243.py`) — within 0.5% of
   `h1union`'s; report the predicted board raw `mse` from 197's identity at the arm's own terms.
6. `vcc prep --dry-run` clean; nnz vs 4.75e9; max counts/cell vs 1e6.

## Deliverable
`submissions/249-cd4svd5-f034.{h5ad,vcc}` with sha256 and bytes; `RESULTS.md` with the three known-answer
rows, every falsifier's number, the config diff verbatim, the builder commit hash, and one line
**BUILDABLE / NOT BUILDABLE**. Do not choose a model name, do not write a description, do not submit.
If the flag cannot be made inert or the staged arm does not reproduce 246, stop and report.
