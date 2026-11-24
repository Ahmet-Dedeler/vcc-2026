# 365 — DESIGN: a one-variable context-A arm (`VCC_CTX_A_PRIOR`) on top of 319p

Read-only scoping. No builder, pipeline, or source file was edited, no build was
run. All file:line citations are against git HEAD as read 2026-09-10 UTC
(builder `scripts/make_submission_pm_cd4.py`, 6348 lines).

Base: **319p** = `VCC_NORMREL_EXP=-1.5` + `VCC_NORMREL_MODE=cap` + `VCC_CD4_PARTIAL_G=0.34`
(`experiments/319-cap-variant/run319p.sh:30-33`), fixed `VCC_ALPHA=0.982`.
Recipe: `experiments/319-cap-variant/run319p.sh`. Inert-flag pattern: 254
(`experiments/254-partial-rows/patch254.py`, `RESULTS.md` §1). Config-diff gate:
AGENTS.md rule 10; `scripts/def_postbuild.py`.

---

## 0. Bank facts constraining the design (measured here, read-only, numpy only)

| fact | value | source |
|---|---|---|
| `data/donors/nadig_jurkat/delta.npy` | 2270 × 8881, float32, all finite, log1p(TS=5e4) pseudobulk deltas | header read |
| median row norm recomputed | **12.371221542358398 == `meta.json` to all digits** | KA-1 PASS |
| bank ∩ A/B/C panel | **0 / 300** | KA-1 PASS (agrees with 225 §8.3 / 235) |
| bank genes ∩ official 18,533 axis | **8,284** (10,249 official genes unmeasured) | KA-1 PASS |
| 323 raw members re-read | pds 0.88964, mse 0.97495, jac 0.07476, nmae **0.95084**, fid 0.49755, reach 0.25435 | `experiments/323-current-family-validation/work/metrics-cap15.json` |
| 343 headline re-read | NMAE 0.950836 (== 323 to 6 dp), delivered sign 0.588, designed sign 0.600 | `experiments/343-nmae-decomposition/RESULTS.md` |

**Constraint C1 (hard).** `replace` mode (per-target externally calibrated rows,
346's construction) is **vacuous on the A/B/C panel**: the bank carries none of the
300 targets, so there is no row to blend. This is not tunable — no rescaling creates
rows that do not exist. On A/B/C only target-*independent* modes (`sign_block`,
`add_vector`) can act; `replace` becomes meaningful only on a panel the bank
overlaps (October). The builder must **raise** on an empty carry, never ship a
silent no-op (235's lesson; AGENTS.md §8.1 live-defect note).

**Constraint C2.** Any bank-derived object covers only 8,284/18,533 genes. The
artifact must carry its measured mask; unmapped genes contribute exactly 0 and the
loader must assert the mapping by gene NAME (never by position).

---

## 1. Insertion point

**File:line: `scripts/make_submission_pm_cd4.py:5325` — `context_effects = None` in
`run()`, with the correction applied immediately after the shared `effects` matrix
is built and before `write_submission_from_effects(...)` at line 5433.**

Data flow before the point (all context-INDEPENDENT, all untouched):
`build_v9_pre_whiten_effects` (:1410, donor mix + normrel-cap −1.5) →
`apply_perturbmatch_k562_replacement` (:1598, 50:50 K562 slot) →
`apply_v9_whitening` (:1864, γ=0.65 union_pruned) →
`apply_frozen_cd4_residual` (:2700, CD4 k=5 deflation → `solve_cd4_dose` → merged
add of CD4 + H1-union + partial → **ONE global preserve rescale**, :2854–2860) →
`build_production_effects` (:2943) returns shared `final`, 300 × 18,533.

After the point (per-context, rides existing code): `write_submission_from_effects`
(:4370) already branches per context at :4502–4509 — `effects` for all, or
`context_effects[str(context)]` when given. That branch is `VCC_ALPHA=knee`'s only
entry point (:4396–4399, precedent block :5315–5427). The A-arm sets:

```python
context_effects = {"A": corrected_A, "B": effects, "C": effects}  # ≤15 lines
```

when the flag is on; `None` (today's path, byte-identical) when off.

Normalisations that act on the correction, in order:
1. **The global preserve rescale does NOT act on it** — it fired inside
   `apply_frozen_cd4_residual`, upstream. The correction changes A's delivered
   amplitude. `sign_block` and norm-matched `replace` are amplitude-neutral by
   construction (|row| preserved); `add_vector` is not — note it in the meta and
   price it through 096's claim route (211: A's claim reads at cap, so extra
   amplitude there is pure cost).
2. **A's own CPM gate DOES act on it** — `profile_from_effect` (:3246) zeroes every
   gene ≤ 10 CPM on A's own `b_ctrl` (:3262–3265). Correction entries on gated-out
   genes are silently discarded: restrict `sign_block` to the delivered call
   universe (head + CPM block, §2a) and record the gated-out fraction in the audit.
3. **Renderer + emitter** — `expm1(b_ctrl+gated)` render (:3317), then `rc_callset`
   emission keyed on the gated effect (top-400 |H| ∪ top-2200 control CPM,
   `RC_CALLSET_*` :1161–1185) with A's own control cells. Panel coordinates
   (all 300 target genes) are unscored by `pds` and excluded from every member —
   the loader must refuse panel-coordinate support (mirror :2739–2742 asserts).

Refusals at import (rule 10, one variable): `VCC_CTX_A_PRIOR` + `VCC_ALPHA=knee`
together is two per-context mechanisms — refuse (mirrors :838–843
K-vs-basis and :734–741 sizing/dose refusals). Non-A/B/C bundle (`--expect-contexts`
≠ A,B,C): refuse — this arm is identity-bound (041: A = Jurkat); October re-derives.

---

## 2. Flag + artifact (`VCC_CTX_A_PRIOR=<path>`)

Env parse mirrors :885–896 (`CD4_PARTIAL_G`): unset/blank → `None` (inert);
set → path must exist, else raise. Byte-inertness proof follows 254's KA-1 exactly:
staged `effects.final_sha256` == 319p's AND full `.h5ad` sha == 319p's
(`submissions/319p-cap15-partial.h5ad.meta.json`) with the flag unset through the
patched builder; plus the pre-patch byte copy importing with the flag SET and
carrying no attribute (254 `RESULTS.md` §1 table, `envcheck254.py` 11/11 pattern).

Artifact: one `.npz` + sidecar `.json`, built offline by whoever owns 346/362/363:

```
genes    (G,) unicode — gene symbols on the OFFICIAL axis order for arrays below
targets  (T,) unicode — panel symbols carried (replace mode; empty otherwise)
matrix   (T,G) float64 — calibrated rows, official axis, 0 on unmeasured genes (replace)
vector   (G,) float64 — common vector, official axis, 0 on unmapped genes (add_vector)
sign     (G,) int8 — {-1,0,+1} gene prior, official axis (sign_block)
mode     unicode — one of sign_block | add_vector | replace
measured (G,) bool — bank measured mask on the official axis
meta: bank path + file shas, panel id, builder of artifact, dose defaults
```

Loader validation (all loud raises): `mode` known; shapes consistent with the
bundle's 300 × 18,533; `genes` maps onto the official axis by NAME with zero
misses; every array finite; **zero support on the 300 panel-coordinate columns**;
`targets ⊆ panel`; `replace` with `targets ∩ panel = ∅` raises (C1 — the
A/B/C case); unknown extra keys raise.

Application (all on the post-preserve design matrix, A rows only, B/C bit-identical):
- **(a) `sign_block`** (363): block = top-2200 by **A's own NTC control CPM**
  among non-head universe genes (head = top-400 by |design| — membership is
  |·|-invariant so the partition cannot shift under the intervention, 363 PLAN).
  `design_A[block] = sign[block] * |design_A[block]|`, `0`-prior genes untouched.
  Magnitude-preserving elementwise → row norms, preserve scalar semantics, and
  `mse` claim untouched.
- **(b) `add_vector`** (362): `row_A(t) = row(t) + α·s·c`, `s = pooled-median|D| /
  pooled-median|c|` over measured genes (362's scale rule; recompute at build,
  record in audit), `α` from the artifact sidecar (default 1.0; PI doses).
  Amplitude-changing → record `√(A_A/A_ref)−1` per context (243's stream) in audit.
- **(c) `replace`** (346): for carried targets,
  `row_A(t) = (1−λ)·row(t) + λ·M(t)` then rescale to `‖row(t)‖` (346's fixed-norm
  blend; `λ` from sidecar, default 0.5). Uncarried rows untouched. Raises on A/B/C.

Config keys (the `h1_union_g` convention, :5528–5536 — absent when off):
`ctx_a_prior` (path sha, not the path alone — rule 11), `ctx_a_prior_mode`,
`ctx_a_prior_dose` (α or λ or "sign"), `ctx_a_prior_targets` (carried count).
Effects audit: `effects.ctx_a` block mirroring `effects.cd4_partial` (carried rows,
support entries, gated-out fraction, per-context amplitude ratio, pre/post
`final_sha256`). The one-key diff gate (rule 10) then reads exactly these keys
plus their consequences — same shape as 254's F1 row (config +1, audit changed =
`final_sha256` + preserve echoes, added keys all under one prefix, everything else
pinned).

Pinning: record the artifact's sha256 in config at build time AND verify the file
at load (the `UNION_SOURCE_REGISTRY` pattern, :917–1020 — path AND sha, verified
in `verify_union_source_hashes` :2569). Under the 80-minute window a registry
entry in source is preferred (rule 11, reviewable diff); fallback is config-only
recording with load-time hash check — state which was used in the round RESULTS.md.

---

## 3. Invariants touched / untouched

| invariant | touched? |
|---|---|
| Norm preservation (`output == whitened·scale` on panel/unsupported entries, :2861–2875) | untouched — asserts fire upstream of the insertion point |
| `--expect-contexts A,B,C` + manifest agreement (postbuild :240–244) | untouched |
| Panel-support repin (`perturbmatch_replacement_targets`, `cd4_complete_targets`) | untouched — the prior is not a donor; counts unchanged |
| Frozen input/implementation hashes | builder source edit moves its own bytes (record in commit msg, 254 §1 pattern); `FROZEN_IMPLEMENTATION_SHA256` does not pin the builder itself (:62–103), but `run319p.sh:36-41`-style hash gates must be updated and the new hash cited |
| Postbuild config diff (`ALLOWED_CONFIG_DIFFS` = contexts/targets/expect_contexts, postbuild :81) | expected FAIL on the new keys → run postbuild with `--allow-config-diff ctx_a_prior[,ctx_a_prior_mode,...]` (:452); `def_pipeline.py` hardcodes `--allow-config-diff whitener_id` (:614) so the A-arm build passes the extra value explicitly (or via `--allow-arm-drift` only if deliberately dropping a shipped arm — NOT this case) |
| Novelty ledger (`check_submission_novelty.py check`) | new bytes ⇒ novel by construction; `record` after submit with entry id |
| `assert_build_carries_arms` (driver :171–248) | extend the asked/got table with the `ctx_a_prior*` rows (else the driver cannot see a dropped arm — 239 §2.1 failure class) |
| RNG streams (`make_rngs`, seeds) | untouched — correction is deterministic, pre-RNG |
| Incumbent reference for postbuild | 319p's meta (`submissions/319p-cap15-partial.h5ad.meta.json`), NOT 249/222 — diff must read exactly the §2 keys |

---

## 4. Test plan

Pure unit tests on a tiny synthetic frame (4 targets × 12 genes, 2 fake contexts;
no h5ad, no builder import beyond the new function — 254 `unit254.py` pattern,
44/44 there). Each mode: hand-computed expected array, exact equality.
1. OFF: `context_effects is None`; shared matrix and a canary audit bit-identical
   to the pre-patch function on the same frame (KA-1 in miniature).
2. ON, each mode: output == hand computation; B/C rows == OFF to the last bit;
   panel-coordinate support raises; non-finite raises; gene-name miss raises.
3. `replace` with empty carry raises (C1); `sign_block` with all-zero prior returns
   rows unchanged but still writes the audit (non-vacuity checkable).
4. Config block: keys absent when off, present with exact values when on.
Byte-inert test (254 KA-1 full): patched builder, flag unset → staged
`effects.final_sha256` == 319p's; full gated build → `.h5ad` sha == 319p's
(`run319p.sh` gate pattern: memory ≥ thresholds, no other builders, hash-gate the
builder bytes first). 20-minute dry-run path: unit suite (~2 min) + `def_pipeline.py
--dry-run --stages build` to verify flag threading and the `--allow-config-diff`
value end to end (prints commands, runs nothing) + synthetic-frame integration
through the new function. Honest limit: a full 20 GiB/≈22 min build plus
≈24 GiB prep cannot fit in 20 minutes — the dry path covers everything except
wall-clock build verification, which rides the normal gated build.

---

## 5. Time estimate (decision → upload, single engineer, gate open)

| step | min |
|---|---:|
| Patch builder (env block + loader + 3 applicators + config/audit keys, ~120 lines, exact-anchor style per `patch254.py`) | 25–35 |
| Unit suite + byte-inert staging check | 10 |
| Gated full build (≈1300 s class at 20 GiB, alone; 254 §2.2: 1331 s / 20.0 GiB) | 25 |
| Postbuild + rule-10 one-key diff + novelty `check` | 5 |
| Real `vcc prep` (runbook: 12–26 s / 22–24 GiB; run319p writes the `.vcc`) + upload (13–16 min at ~3.8 MiB/s) | 20 |
| **Total** | **≈ 80–95** |

Fits the ~80-minute window only if the memory gate opens immediately and no
rework: pre-stage the patch as an exact-anchor script (254 pattern) BEFORE the
decision so decision-time work is review + run, not writing. If 346/362/363 land
mixed (e.g. only `sign_block` clears its gate), ship that mode alone — modes are
independent code paths sharing only the loader.
