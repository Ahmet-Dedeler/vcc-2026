# 222 — BUILD `h1union`: a coverage-aware sizing rule + the two H1 sources. DO NOT SUBMIT.

You are a build subagent. The arm is preregistered at
`experiments/222-h1-union-build/PREREGISTRATION.md` — **read it first; it is binding and must not be
edited.** **NEVER run `vcc submit`** — the coordinator submits.

Work in `experiments/222-h1-union-build/`, plus `submissions/` for the package and
`scripts/`+`src/` for the sizing rule. `uv run python`.
**Peak RSS ≤ 22 GB, builds strictly SERIAL** (210 measured 19.8 GiB). Three research agents are
running; check free memory before the emission stage and wait if under 25 GB.

## READ FIRST
1. `experiments/222-h1-union-build/PREREGISTRATION.md` — binding, five falsifiers.
2. `experiments/218-independent-donors/RESULTS.md` — **the science. Its §7 has the mechanism, its
   §7.2 the 12/12 fold×dose table, its §7.3 the four falsifiers it already ran.** You are the
   engineering; do not re-derive its measurements, reproduce the ones it names as known answers.
3. `experiments/210-cd4-build/` — `build210.sh` is your starting script and `RESULTS.md` its gates.
4. `scripts/make_submission_pm_cd4.py` — `solve_cd4_fraction` and the residual staging.
5. `AGENTS.md` §6.1, §8.4, §11.1 **rule 10**, §14.

## The change — one variable, but it is a CODE change, so the bar is higher

**The bug:** `solve_cd4_fraction` sizes a residual by the **median row norm over every panel row**.
A source covering ≤ 50% of the panel has a median row norm of **exactly 0**, so its solved fraction
is **exactly 0**. Pan-2026 H1 (64/300) and `h1_2025` (20/300) are therefore contributing **nothing**
to every build ever made, silently.

**The fix:** make the sizing rule **coverage-aware** — size a residual on the rows it actually
carries, not on all 300. Then add the **union of the two H1 sources** (77/300, 7 shared) as an
un-whitened post-whitening residual at 218's chosen dose.

**Design constraints, non-negotiable:**
- **CD4's own dose, form and solved fraction must not move.** `cd4_target_f` stays **0.34**,
  `effects.cd4.*` must be **untouched** in the built meta. If your fix changes CD4's fraction by any
  amount, it is two variables and you must stop and report.
- The change must be **behind a flag, defaulting to the current behaviour**, so the default path is
  provably unchanged.
- Follow 218's construction exactly: un-whitened, post-whitening, sized on carrying rows. It
  measured `pds` **up** on 12/12 cells; a build that loses `pds` is not its arm.

## Known-answer rows (mandatory, abort on failure)
1. **With the flag off**, rebuild `cd4f034` and require sha256
   **`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`** on
   **12,648,807,800 bytes**, by `shasum` on the file. If the default path moved at all, stop.
2. **Prove it against git, not against the current tree** — resolve the pre-change builder out of
   git and diff behaviour, the way 191 and 202 did. A test that only checks your own new code is
   not a control.
3. Reproduce **218's** own headline on your staged effects before emitting: macro cosine **+7.68%**
   and `pds` **+0.0049** on the three cancer-line folds, and the non-carrying-row bound
   **|Δ| ≤ 0.018 pp**. If your union differs from 218's, these will not reproduce — that is the check.

## Falsifiers — report each with its measured number
1. The preregistration's falsifier 5: **config diff is exactly the pre-declared sizing change plus
   its deterministic consequences**, `operator.sha256` unchanged, `alpha` **0.982**,
   `cd4_target_f` **0.34**, `effects.cd4.*` untouched. **Print every difference, not only the
   expected ones.**
2. **The profile moved on the carrying rows and nowhere else.** Per-gene panel-total digest against
   `210-cd4f034.h5ad`: expect it to **differ** on rows the union carries and to be **identical**
   elsewhere. Report the counts. This is the sharpest available check that the arm is what it claims.
3. **Non-vacuity**: `‖ΔE‖/‖E‖` against `cd4f034`, with a per-block minimum and a count of zero
   blocks. **A zero here is the exact bug you are fixing** — 200, 207 and now 218 have all found a
   silently-zeroed residual, so assert it loudly.
4. **Challenge `n_pred` within ±2%** of `cd4f034`'s (probe against `data/challenge/context_*.h5ad`).
5. **`vcc prep --dry-run` clean**; nnz against the 4.75e9 cap; max counts/cell against 1e6.

## Deliverable
`RESULTS.md`: the two known-answer rows, 218's reproduction, every falsifier with its number, the
config diff verbatim, the package path and sha256, and a one-line **BUILDABLE / NOT BUILDABLE**
verdict. Then stop.

**Do not submit. Do not choose a model name or write a description.**

**If the coverage-aware rule turns out to change CD4's fraction and you cannot separate them, say so
and stop** — a two-variable arm is worth less than an unspent slot, and we have an unspent slot.
