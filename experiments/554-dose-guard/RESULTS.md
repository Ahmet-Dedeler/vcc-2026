# 554 — the emitter-dose guard (UNAPPLIED patch + tests), and what calibrating it found

Owner: session "Explore contribution areas", at the coordinator's request (2026-10-01 ~09:00Z; due 10-08 12:00Z). Class
(AGENTS §11.1c): **scorer mechanics / delivery** — it changes no prediction; it stops a build that silently lost its DE
design. Written 2026-10-01 ~09:45Z.

## Why
553 found that the K-class emitter (`src/vcc2026/rc_emit_k.py`) can deliver its DE design at dose 0 **without any error**:
its library classes are the context's control quartiles, and on H1 the steering dose (CLEAN4 block prior, coarse MDE, the
uncalled-gene hook that rides on the MDE) collapses to 0 between control-library q75/q25 3.5 and 5, and the placement
crashes by 7.24. Mechanism (rc_emit_k.py:1114–1165): the prior's log2 move is bisected down (`lam`, 10 steps) until every
class keeps a fixed 6,000-count margin; `lam` 0 means no steering at all, and the MDE is skipped whenever `lam` is 0.

## The patch (`patch_guard554.py`; applies AFTER 538's applier; refuses a second or partial apply; atomic writes)
| file | change | default |
|---|---|---|
| `scripts/def_preflight.py` | per context, control-library q25 / q75 / **q75/q25** from the library sizes it already streams, always recorded in the preflight JSON (`library_spread_q75_q25`) | a FAIL row only when q75/q25 > **3.0**; `DEF_DAY_ALLOW_WIDE_LIBS=1` makes it a passing "ALLOWED" row (JSON `library_spread_allowed`) |
| `scripts/def_postbuild.py` | `--rc-audit <jsonl>`: per-context medians of the builder's per-block `prior_lambda`, `mde_lambda`, `call_set_missed`; the dump must belong to the meta (same block count, and its pooled medians equal the meta's `output.rc_callset.*.median`); without `--rc-audit`, the meta's pooled medians (a global collapse only); JSON gains `dose_guard` (values, floors, margins) | a FAIL row only when a median is past its floor (**prior ≥ 0.9, MDE ≥ 0.9, missed ≤ 100**) or the dump is not the meta's; `DEF_DAY_ALLOW_LOW_DOSE=1` makes it a passing "ALLOWED" row |
| `scripts/def_day.sh` (538-patched) | every level build runs with `VCC_RC_AUDIT_DUMP=$W/$N-<tag>.rc_audit.jsonl`; the final build's dump is linked as `$W/${N}_final.rc_audit.jsonl` and passed to def_postbuild as `--rc-audit`; non-final dumps are deleted with their h5ads; both overrides logged at start and written to `def_day.final.env` | always on (diagnostics only) |

**Why def_day.sh too:** the final meta records the emitter's diagnostics only POOLED (min / median / max over all 900
blocks, `make_submission_pm_cd4.py:5099`). A one-context collapse (300 of 900 blocks) moves neither the pooled median nor a
minimum that is already 0 in healthy builds. The per-block record exists only as the builder's `VCC_RC_AUDIT_DUMP` jsonl
(`make_submission_pm_cd4.py:942-944`: "diagnostics only, the emission is untouched. Unset is a no-op"; not recorded in the
config, so no config diff).

**Disk budget of the dump:** one JSON line per block, 900 lines, **1.3–1.6 MB per build** (the seven dumps on disk). A day
writes lvl0 + up to `MAX_LEVELS` (4) level builds, so ≤ ~8 MB transient in `$W`; the non-final ones are deleted with their
h5ads, leaving **~1.6 MB** (`$W/$N-<tag>.rc_audit.jsonl` + its hard link `$W/${N}_final.rc_audit.jsonl`).

## Calibration (`calib554.py` → `work/calib554.json`)
| quantity | data | floor | margin |
|---|---|---|---|
| control-library q75/q25 | every bundle on disk (A/B/C, 180 def-bundle, 199 def-bundle-r, efp, efr, R6 S1, R7a T1 — all byte-identical clones of A/B/C's controls): **A 1.947, B 1.572, C 1.855**; H1 1.575; R7b as built **7.238 (trips)**; R7b floored 2.124 | 3.0 | ≥ 1.05 on every real bundle |
| per-context MDE dose (dumps) | 433a (coarse MDE, production since): A 0.956, B/C 1.0; October rehearsals 432 b0/b1/b2/ef_b1 and 448 efr lvl1: D 0.951–0.958, E/F 1.0; **429a (plain MDE, the known failure): A 0.0, B 0.936, C 0.856** | 0.9 | A/D 0.051–0.058; 429a trips on A and C |
| per-context prior dose (dumps) | 1.0 in every context of every dump (429a included) | 0.9 | 0.1 |
| per-context missed calls (dumps) | medians 4–21.5 | 100 | ≥ 78 |
| pooled medians (no-dump fallback), the coordinator's list | 446bj2z / zcp / zcph, 453dr, 457dr, 460dr, 497dr, 474efp, 497efr, 510efr, def452: prior 1.0, MDE 1.0, missed 6–13 (none trips); **465efr: prior 0.278 (trips)** | 0.9 / 0.9 / 100 | — |

**On the MDE floor** (the coordinator asked whether 0.9 sits too close for A): A/D's healthy median is 0.951–0.958 across
six builds (spread ~0.007), so 0.9 leaves 0.05 (~15× that spread) and still catches 429a's C (0.856), a real delivery
defect of the pre-coarse MDE. Kept at 0.9. If a D/E/F context reads 0.85–0.9, it is a ~10–15 % MDE dose loss (≈ −0.0003
overall): stop, look, and override with `DEF_DAY_ALLOW_LOW_DOSE=1` if nothing else is wrong.

## What calibrating it found (reported 2026-10-01 ~09:30Z; confirmed and recorded by the coordinator, c15e7935)
The board's **programme-prior builds delivered under half their prior** (`rc_callset.prior_lambda`, pooled median over 900
blocks): 458 **0.243**, 458d (the 09-30 HARM read) **0.433**, 461dr **0.441**, 465efr **0.278**; gated forms 458g / 469dr
0.953, 458e 0.991; every current-design CLEAN4 build on the coordinator's list 1.0 (an older emitter, 392's, had a
pooled median 0.901). **On H1 the same programme priors were delivered at λ 1.0 on 83 of 85 targets** (MED23 0.85–0.93 and
PHF12 0.97–1.0, in programme and CLEAN4 runs alike, over 50 diag files; corrected 2026-10-01 after the other peer's audit 593
and the coordinator's re-check, 57633e25)
(462 arcbank/arcfull, 467 gated*, every 517 arm; their diag jsons). Likely mechanism: the programme's balanced up/down
signs (~56 % down vs CLEAN4's 84 %) need class capacity A/B/C's thinner sub-gate mass cannot supply. Coordinator's reading:
458d's HARM **size** is confounded with delivery, its **direction** stands (458g at λ 0.95 still read −0.0049 vs 446b);
October is unaffected (programme off, CLEAN4 at λ 1.0). 446bja (A only) is pooled 1.0 / min 0.243, so its A-context dose
is unknown; a dose-only rebuild with the dump would settle it (listed as optional, not run).

## Tests
- **`test_guard554.py` (fast, temp trees, live files hashed before/after): 5 tests OK.** A: `--check` clean and writes
  nothing, apply, second apply refused, a partial tree refused, an anchor miss writes nothing; `bash -n` and `py_compile`.
  B: `dose_guard` on synthetic metas/dumps (healthy; one context collapsed → one FAIL row naming it; prior and missed
  floors; ALLOW → passing row; a dump that is not the meta's; pooled fallback healthy and collapsed; a non-K-class meta) and
  on the real dumps (433a: no row, A's MDE 0.95605; 429a: FAIL on A 0 and C 0.8555, not B; 433a's dump against 429a's meta:
  "does not belong"). C: `def_day.sh` on 538's stub harness — with every knob unset (and with the October `LODO=1, DE_MAG`
  env) the 554 day makes exactly the 538 day's calls except `VCC_RC_AUDIT_DUMP=$W/<name>.rc_audit.jsonl` on the level
  builds (not on the stack-dump build) and `--rc-audit $W/day_final.rc_audit.jsonl` on def_postbuild; `def_day.final.env`
  identical; the final dump linked, the non-final ones pruned; with both overrides set, both are logged and recorded.
- **`identity554.py` (real data, ~60 s, peak 0.45 GB) → `work/identity554.json`:**
  - def_preflight, base (538 tree) vs patched, knobs unset: **identical text and exit code on all 6 bundles** (masking only
    the elapsed-seconds and "wrote <path>" lines); the JSON gains the spreads.
  - def_postbuild, no `--rc-audit`: **identical text and exit code on 11 of 12 metas**; on **465efr** the patched script adds
    exactly one FAIL row ("pooled: median prior_lambda 0.2778 past its floor 0.9"; 37/41 → 37/42 checks) — the guard
    catching the under-delivered programme prior, not a false alarm. (Every run exits 1 in both versions for unrelated
    rows of these old metas against today's arms table.)
  - Through the CLIs: a bundle with R7b's as-built controls in every context fails preflight (q75/q25 7.238 ×3) and passes
    as "ALLOWED" with `DEF_DAY_ALLOW_WIDE_LIBS=1`; def_postbuild `--rc-audit` on 433a adds no row, on 429a adds the FAIL
    row (A 0, C 0.8555) and, with `DEF_DAY_ALLOW_LOW_DOSE=1`, the passing ALLOWED row.
- **Limitation, stated:** the dump/meta consistency check (block count + pooled medians) catches a stale, truncated or
  mixed dump, not a different build whose pooled medians happen to be identical.
- **Not run here (needs the real builder):** that `VCC_RC_AUDIT_DUMP` leaves the output h5ad byte-identical. The coordinator
  runs it when applying: the inert KA (efp relaxed through the patched def_day, where the dump is now always on) must still
  produce `ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3`.

## Apply (the coordinator, after 538's applier on 10-05, outside 10-05 00:00–18:00Z for anything heavy)
```bash
uv run python experiments/554-dose-guard/patch_guard554.py --root . --check
uv run python experiments/554-dose-guard/patch_guard554.py --root .
uv run python experiments/554-dose-guard/test_guard554.py
uv run python experiments/554-dose-guard/identity554.py
bash -n scripts/def_day.sh
```
Then the inert KA (550's step 4 command) — it now also writes `$W/474efp_final.rc_audit.jsonl`; check the h5ad sha, that
`postbuild.json` has a `dose_guard` record with per-context values (D/E/F ≈ 0.95–1.0), and commit via
`commit_after_chain.sh`. The patch refuses to run on a tree without 538's edits (anchor miss, nothing written).
`calib554.py` / `identity554.py` take their inputs from module-level lists: add 524dr's meta once it exists.
