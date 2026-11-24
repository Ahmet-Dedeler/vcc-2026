# 373 — `blend_rows` mode for `VCC_CTX_A_PRIOR`: chain applies, 100/100 tests pass

Status line: **ENGINEERING_ONLY** (no score claim; staged patch under
`work/patched/`, `scripts/` untouched; PI applies per README §1 when a
blend_rows artifact clears its gate).

## Known-answer reproduction (numbers, match)

| quantity | reproduced | recorded | match |
|---|---|---|---|
| 366-only sha on HEAD bytes | `95913dc0…` | 366 RESULTS `95913dc0…` | exact (chain reproduces 366) |
| 366's 45 tests on the combined copy | 45 pass | 366 RESULTS 45/45 | exact, unmodified rerun |
| 323 raw `pds/mse/jac/nmae/fid/reach` (inside that rerun) | 0.88964/0.97495/0.07476/0.95084/0.49755/0.25435 | 323 via 366 | within 5e-6 |

If any of these had failed, work would have stopped per PLAN.md F-callout. None did.

## What was measured (inputs, estimator, controls)

- Input: `scripts/make_submission_pm_cd4.py` at HEAD `1318cd3c…` (6348 lines).
- Estimator: `patch373.py` — full chain `E5–E8 ∘ patch366.patch` (E5 loader
  branch incl. mode-aware key sets + inline blend parse; E6
  `apply_ctx_a_blend_rows` + `panel_targets` dispatcher; E7 call-site +
  mode-aware print; E8 mode-conditional config), applied to a byte copy at
  `work/patched/`. `--check` runs the chain on HEAD bytes, prints
  pre/mid/post sha256.
- Controls: 366-pattern inert-flag discipline (combined copy imports clean with
  flag unset — proven by the rerun); synthetic frame 5 rows × 50 genes with
  independently recomputed expectations (separate code path sharing only
  `np.linalg.norm` — the op-independent property check is the 1e-12
  fixed-norm test); 20 load refusals + 6 apply refusals + old-mode coexistence;
  E8 line evaluated from the shipped source for both modes; `blend_rows`
  absent from `scripts/`.

## Results tables (units: checks, all must pass)

| group | checks | result |
|---|---|---|
| 366 rerun on combined copy (known-answer + inert + sign_block + p_down + add_vector + 9 malformed + knee/panel/miss/vacuous/context + untouched) | 45 | 45 pass |
| blend exact (bitwise vs recomputation, audit fields, cosine, norm identity 0.00e+00, untouched rows/genes, B same-object, C bitwise, A differs+equals) | 15 | 15 pass |
| beta edges (0 ⇒ bitwise identity via hand-built prior; 1 ⇒ pure direction + norm, 2 rows) | 5 | 5 pass |
| load refusals (p_down, vector, tau, alpha, scale_rule, NaN, inf, beta∈{0,1.5,NaN}, dup targets, dup genes, shape, missing beta, extra key, bad norm_rule, empty, object dtype, 2 reverse) | 20 | 20 pass |
| apply refusals (off-panel target, gene miss, panel support, zero row, panel_targets=None, wrong context) | 6 | 6 pass |
| old-mode coexistence (sign_block identical with new kwarg; no blend keys) | 2 | 2 pass |
| E8 config expression (old key set, new key set, targets count==2) | 3 | 3 pass |
| scripts/ untouched + anchors + loader defaults | 4 | 4 pass |
| **total** | **100** | **100 pass, 0 fail** |

Patch arithmetic: HEAD `1318cd3c…` → 366-only `95913dc0…` (+155) →
366+373 `c6f366ca…` (+116, total +271). No F-falsifier fired.

## Verdict against the PLAN.md kill criterion

**Not killed — ships as a staged patch.** F2 required: unset import clean ✓
(rerun 45/45 bit-exact — the dispatcher branches, the old path is unedited
except `kw.pop`); norm identity worst rel err **0.00e+00** vs 1e-12 gate ✓;
β=0 bitwise ✓; untouched rows/contexts bitwise ✓; every mandated refusal
fires ✓; `scripts/` clean ✓.

## Limits

- The combined copy has never been through a real build: byte-inertness at
  staging and the rule-10 five-key diff are proven only at the unit level; the
  PI's step (d) build must re-prove 254 KA-1 before any flagged run.
- Panel-membership of artifact targets is apply-time, not import-time (the
  loader cannot see `spec.perts`) — a blend artifact with an off-panel target
  fails 29 min into nothing: it fails at the `[ctx-a-prior]` step before any
  rendering, loudly. Validate artifacts with `test373`-style apply checks first.
- A design row exactly zero on the covered genes is skipped (identity there);
  `rows_touched` counts written rows, `artifact_targets` counts artifact rows.
- `assert_build_carries_arms` has no `ctx_a_prior` row (366's limit, inherited).
- No emission, no scorer, no A/B/C forecast is claimed.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50% (was 76%)
uv run python experiments/373-blend-rows-patch/patch373.py --check
uv run python experiments/373-blend-rows-patch/patch373.py
uv run python experiments/373-blend-rows-patch/test373.py   # 100 passed, 0 failed
```

## Peak memory and wall time

`test373.py`: peak self RSS **240 MiB** (`resource.ru_maxrss`), wall **~5 s**
(two builder imports + 100 checks). `patch373.py --check`: < 1 s, < 100 MB.
No process exceeded 4 GiB; synthetic npz live in a `TemporaryDirectory`;
no downloads, no builds, no submits.
