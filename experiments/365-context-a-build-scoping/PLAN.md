# 365 — PLAN: read-only scoping of a one-variable context-A arm on the production builder

## Hypothesis (design-level)
If 346/362/363 show context A (= Jurkat) can be improved from the public Nadig
Jurkat bank, the production change is a single env-gated, unset-inert, per-context
correction applied to context A's design effect rows only, on top of 319p
(`VCC_NORMREL_EXP=-1.5`, `VCC_NORMREL_MODE=cap`, `VCC_CD4_PARTIAL_G=0.34`), reusing
the `context_effects` per-context entry point that `VCC_ALPHA=knee` already
established in `write_submission_from_effects`. One variable, one sitting, ~80 min
decision-to-build.

## Exact inputs (all read-only; nothing outside this directory written)
- `scripts/make_submission_pm_cd4.py` (builder; composition flow `build_v9_pre_whiten_effects`
  → PM replace → whiten → sign-repair → `apply_frozen_cd4_residual` → `build_production_effects`
  → `run()` → `write_submission_from_effects`; per-context entry `context_effects`, knee
  precedent at run() lines ~5325–5382/5433–5438)
- `scripts/def_pipeline.py` (`day_env`, `SHIPPED_ENV`, `assert_build_carries_arms`,
  `def_postbuild.py` config-diff gate, novelty `check_submission_novelty.py`)
- `experiments/319-cap-variant/run319p.sh` (gated build recipe; 319p = cap + partial rows)
- `experiments/254-partial-rows/patch254.py` + `RESULTS.md` §1 (byte-inert flag pattern, KA-1)
- `docs/DEF-RUNBOOK.md` (§2 stages, §3 decision table, §4 known-answer rows)
- Candidate constructions: `experiments/346-same-line-calibration-h1/PLAN.md` (ridge map →
  per-target calibrated matrix = `replace`), `experiments/362-common-vector-oracle-sixmember/PLAN.md`
  (common vector add, scale s = pooled-median |D|/pooled-median |C| = `add_vector`),
  `experiments/363-block-down-sixmember/PLAN.md` (design[block] = −|design[block]| = `sign_block`)
- Bank facts: `data/donors/nadig_jurkat/` (`delta.npy` 2270×8881 float32 log1p(TS=5e4)
  pseudobulk deltas, `genes.txt`, `perts.txt`, `meta.json`)

## Estimator / method
No estimator. Read the six sources above, verify line numbers by direct read, run
only small read-only probes (numpy/pandas, < 1 GiB): bank∩panel overlap, bank gene
coverage, one KA cross-check each against 323 (`metrics-cap15.json` raw members)
and 343 (NMAE 0.950836, sign 0.588/0.600). Write `DESIGN.md` (§1–§5 of the mandate)
so a competent engineer can implement in one sitting; snippets ≤ 15 lines,
illustrative only.

## Controls / falsifiers
- KA-1: nadig_jurkat median row norm recomputed == `meta.json` 12.371221542358398;
  panel overlap == 0/300 (225/235's record); 323 raw nmae == 343's 0.950836.
- Falsifier that kills the `replace` mode on A/B/C: bank∩panel = 0/300 means a
  per-target replacement matrix carries no A/B/C row — stated as a constraint, not
  worked around (no rescaling trick creates rows that do not exist).
- Falsifier for the whole doc: if the insertion point cannot be stated as one
  file:line with the normalisations before/after enumerated, the doc fails its
  mandate ("write it so a competent engineer could implement it in one sitting").

## What result would kill the idea
- `context_effects=None`-when-unset cannot be shown bit-inert (it can — the knee
  path already proves it), or the A-correction cannot be confined to one config key
  plus its consequences under the postbuild diff gate.
- Any finding that the correction must touch `scripts/` shared semantics (whitener,
  preserve, emitter) rather than ride them: that is not a one-variable arm.

## Budget
3 h. No builds, no `vcc submit`, no `def_pipeline.py`, no `make_submission*.py`
invocation, no edits outside `experiments/365-context-a-build-scoping/`.
