# PLAN.md — 366-ctx-a-prior-patch

## Hypothesis
A pre-staged exact-anchor patch script (254's pattern) can add a one-variable
context-A prior arm (`VCC_CTX_A_PRIOR`) to a COPY of `scripts/make_submission_pm_cd4.py`
with a ≤~150-line diff, byte-inert when unset, applying `sign_block` or `add_vector`
to one context's post-preserve design rows only, with all malformed inputs refused
loudly. Engineering deliverable: patch script + pure tests + README. No score claim.

## Exact inputs (paths)
- Builder HEAD: `scripts/make_submission_pm_cd4.py` (6348 lines; anchors verified by read
  2026-09-10 UTC: env block ends L893–897; `def build_production_effects(` L2943;
  knee tail `del frames`/`gc.collect()` + `del donors, ...` L5426–5429;
  config `cd4_partial_g` line L5518).
- Pattern: `experiments/254-partial-rows/patch254.py` (replace_once + --check).
- Inert-flag discipline: `experiments/254-partial-rows/RESULTS.md` §1.
- Block definition: `experiments/361-fidelity-anatomy/RESULTS.md` (head = top-400 |design|
  per row; block = top-2,200 control-CPM among non-head).
- Loader constants: `RC_CALLSET_HEAD=400`, `RC_CALLSET_CPM=2200` (builder L1173–1174).
- Control CPM source: `load_control_frame(...).b_ctrl` (builder L4159); ranking by CPM ==
  ranking by `b_ctrl` (both monotone in per-gene control totals; ties identical at zero).
- Prior artifacts (for README conversion recipes only, not loaded here):
  `experiments/370-jurkat-block-sign-prior/work/jurkat_block_prior.npz` (sign_block source),
  `experiments/350-same-line-priors-data/work/jurkat_common.npz` (add_vector source).
- Known-answer record: `experiments/323-current-family-validation/work/metrics-cap15.json`
  raw six members + `effects.npz::cap15` / `reference.npz` for the block-size recomputation.

## Estimator / construction
1. `patch366.py`: four exact-anchor replacements on a byte copy of the builder into
   `work/patched/` (+ `--check` verifying each anchor occurs exactly once in HEAD bytes,
   printing pre/post sha256; never writes `scripts/`):
   - E1 env+loader after the `CD4_PARTIAL_G` block: `VCC_CTX_A_PRIOR` unset/blank → None
     (inert); else load `.npz` AT IMPORT, validate structure (required `genes,mode,vector`;
     optional `p_down,tau,alpha,scale_rule,context` with defaults `tau=0.0, alpha=1.0,
     scale_rule='gate_median', context='A'`; mode ∈ {sign_block,add_vector};
     context ∈ {A..F}; shapes agree; all finite; `p_down` ∈ [0,1]; `tau` ∈ [0,0.5];
     `alpha` finite ≠ 0; unknown extra keys → ValueError). Refuse `+ VCC_ALPHA=knee`
     (two per-context mechanisms, mirrors L838–843). Axis-by-NAME mapping and
     zero-panel-support checks happen at apply time (official axis unknown at import).
   - E2 two pure functions before `def build_production_effects(`: `apply_ctx_a_prior`
     (per-row head/block partition via stable argsort; `sign_block`:
     `design[block] = sign*|design[block]|` with sign from the `p_down` tau-rule
     (else −1) or `sign(vector)` (0 → untouched); genes absent from the artifact →
     untouched; `add_vector`: `row += alpha*s*vector`, per-row
     `s = median|design|/median|vector|` over `b_ctrl>0` universe genes; refuses
     empty universe / zero vector median / panel support / gene misses) and
     `build_ctx_a_context_effects` (None when off → writer path bit-identical;
     else `{ctx: corrected, others: same object}` + `effect_audit["ctx_a"]`).
   - E3 run() block between the knee tail and `del donors, ...`: refuse non-A/B/C
     bundle (identity-bound, 041) and unknown context; load that context's
     `ControlFrame`; apply; set `context_effects`; audit; print one `[ctx-a-prior]` line.
   - E4 config: one dict with exactly `ctx_a_prior_mode/sha256/alpha/tau/context`
     when on, absent when off (h1_union_g convention). Postbuild/drive allowance via
     the EXISTING `--allow-config-diff` flag (254 precedent) — no def_postbuild.py /
     def_pipeline.py patch required; documented in README.
2. `test366.py`: imports the patched copy BY PATH (spec_from_file_location on
   `work/patched/make_submission_pm_cd4.py`, `__file__` pointed at the real builder so
   its `ROOT` resolves; env set before exec). Synthetic frame 3 ctx × 5 targets × 50 genes
   (head_k=5, block_k=10): unset→None; sign_block exact-vs-independent-recomputation +
   head/B/C bitwise identical; add_vector exact + B/C identical; 8 malformed artifacts
   refused (missing key, NaN, bad mode, bad context, extra key, p_down range, tau range,
   alpha=0); knee-combo refused at import; panel-support + gene-miss refused at apply;
   `CTX_A_PRIOR` absent from the real builder source (scripts/ untouched).
3. Known-answer reproduction (own code, numpy only): reload 323 `effects.npz::cap15` +
   `reference.npz`, recompute the 361 partition (head top-400 |design|/row; block top-2200
   `ref_lin` among non-head) and require block size exactly 2200 on all 85 targets, plus
   the six raw members equal DESIGN §0's table to 6 dp. If this fails, stop.

## Controls / falsifiers
- F1: any anchor count ≠ 1 → patch aborts (no partial application).
- F2: kill criterion — the idea is dead if: unset import carries any `CTX_A_PRIOR`
  attribute ≠ None; any synthetic B/C byte differs; any malformed artifact imports;
  head entries move under sign_block; or the added-line count exceeds ~150 by >20%.
- F3: real `scripts/` tree must not contain `CTX_A_PRIOR` after all runs (checked in test).

## What result would kill the idea
Any F2 firing, or DESIGN's L5325 anchor (`context_effects = None`) proving to be the only
correct point (it is NOT: inserting there runs before the knee frames load and would be
overwritten by the knee rebuild; the settled point is the knee tail — documented
discrepancy; knee combo refused so both orders coincide when off).

## Budget
~4 h. Memory: no process over 4 GiB (largest: 85×18533 float64 ≈ 13 MB + builder import).
`memory_pressure` free = 73% at start. No builds, no submits, no downloads.
