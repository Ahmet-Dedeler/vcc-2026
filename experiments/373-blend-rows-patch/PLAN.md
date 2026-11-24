# PLAN.md — 373-blend-rows-patch

## Hypothesis
366's staged `VCC_CTX_A_PRIOR` patch can be extended with a third artifact mode
`blend_rows` — a per-row fixed-norm blend of the named context's post-preserve design
rows toward an artifact matrix (346's rule) — as a chained exact-anchor patch
(366's 4 replacements, then 373's E5–E8), byte-inert when unset, with old-mode
behaviour bit-identical and every malformed blend artifact refused loudly.
Engineering deliverable: `patch373.py` + `test373.py` + README. No score claim.

## Exact inputs (paths)
- Builder HEAD: `scripts/make_submission_pm_cd4.py` (6348 lines; HEAD sha for this
  path `1318cd3c…`; 366 `--check` re-verified 2026-09-10 UTC: all 4 anchors ×1).
- 366's staged patch (read-only): `experiments/366-ctx-a-prior-patch/patch366.py`
  (E1 env+loader, E2 applicators, E3 run() block, E4 config keys),
  `test366.py` (45 tests), `README.md`, `RESULTS.md`.
- Fixed-norm blend rule reference: experiment 346's rule (rescale to ‖r‖ on covered
  genes; cited in mandate).
- Known-answer record: `experiments/323-current-family-validation/work/metrics-cap15.json`
  (reused via 366's test suite, re-run on the combined copy).

## Estimator / construction
`patch373.py` chains: `patch = E5–E8 ∘ patch366.patch` on a byte copy into
`work/patched/` (never writes `scripts/`). `--check` applies the full chain to
HEAD bytes and prints pre/mid/post sha256.
- E5 loader: `CTX_A_PRIOR_MODES += ("blend_rows",)`; mode-aware allowed/missing key
  sets — blend_rows allows exactly `{targets,genes,matrix,beta,mode,context,norm_rule}`
  (so `p_down`/`vector`/`tau`/`alpha`/`scale_rule` are refused here, and
  `targets`/`matrix`/`beta`/`norm_rule` are refused in the two old modes by 366's
  existing unknown-key gate); `targets` must be 1-D unicode (`U`/`S` kind, `O`
  refused), non-empty, duplicate-free; `genes` 1-D, non-empty, duplicate-free;
  `matrix` 2-D float64 of shape `(len(targets), len(genes))`, all finite (NaN/inf
  refused); `beta` scalar finite in (0,1] (`0` refused — see β=0 note); `norm_rule`
  defaults to `'fixed_norm_on_covered'`, any other value refused; `context`
  defaults to `'A'` (366 convention), validated against A–F. Panel-membership of
  targets is NOT knowable at import (official/panel axes unknown) → checked at
  apply time (stated split; README documents it).
- E6 applicator `apply_ctx_a_blend_rows(design, prior, official_genes,
  panel_positions, panel_targets)`: genes map BY NAME (misses raise); genes hitting
  panel-coordinate positions raise; artifact targets outside `panel_targets` raise;
  design must be `(len(panel_targets), len(genes))`; per touched row on covered
  genes C: `r = D[i,C]`, `m = M[k]`; `mn==0 → raise`; `rn==0 → skip (already zero,
  documented)`; `m_scaled = m·(rn/mn)`; `blended = (1−β)·r + β·m_scaled`;
  `‖blended‖==0 → raise`; write `blended·(rn/‖blended‖)` back on C only. Dispatcher
  in `build_ctx_a_context_effects` gains `panel_targets=None` kwarg (default keeps
  all old calls working) and routes `mode == "blend_rows"` to the new applicator.
  Audit: `mode/context/rows/rows_touched (rows written)/covered_genes/mean_cosine
  (over written rows, cos(r,m))/beta/artifact_sha256`.
- E7 run() block: pass `panel_targets=tuple(map(str, spec.perts))` at 366's call
  site; make the `[ctx-a-prior]` print mode-aware (`touched_entries` vs
  `rows_touched` via `.get`).
- E8 config: mode-conditional dict — always `mode/sha256/context`, plus
  `alpha/tau` for the old modes, plus `beta/targets(count)` for blend_rows.

β=0 note: the loader refuses β=0 artifacts (beta ∈ (0,1]). The mandated
"β=0 ⇒ bitwise identical" test therefore exercises the *arithmetic path* with a
hand-built prior dict (`beta=0.0`, bypassing the loader) and proves the rescale is
the identity there (IEEE: `1.0·r + 0.0·m == r` bitwise, `rn/rn == 1.0`).

## Controls / falsifiers
- F1: any anchor count ≠ 1 → abort (366's `replace_once` + 373's; no partial
  application). 373 anchors live in 366-inserted text, so the chain order is load-bearing.
- F2 kill criterion: unset import not inert; any old-mode synthetic expectation
  differs (366's 45 re-run on the combined copy must pass 45/45 unchanged);
  `‖blended‖/‖r‖ − 1` exceeds 1e-12 relative on any touched row; β=0 hand-built
  prior not bitwise identical; any untouched row or non-target context differs by
  even one bit; any mandated refusal imports/applies; `scripts/` contains
  `CTX_A_PRIOR` or `blend_rows` after the run.
- F3: `test373.py` imports `test366.py` BY PATH and re-runs it with `PATCHED`
  overridden to the combined copy (no write into 366's dir); new blend tests use
  an independent recomputation (separate code path, stable-argsort-free).

## What result would kill the idea
Any F2 firing. In particular: if the dispatcher change alters old-mode numerics
(45/45 must be bit-exact — the old path is branched, not edited, except the
`kw.pop`), or if the fixed-norm rescale cannot hold 1e-12 (it is one multiply by
`rn/‖blended‖`; expected error ~1e-16).

## Budget
3 h. Memory: no process over 4 GiB (largest: builder import + 5×50 synthetic
float64 ≈ KB). `memory_pressure` free = 76% at start. No builds, no submits, no
downloads.
