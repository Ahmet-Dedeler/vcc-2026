# 366 — Pre-staged `VCC_CTX_A_PRIOR` builder patch: applies cleanly, 45/45 tests pass

Status line: **ENGINEERING_ONLY** (no score claim; patch staged under `work/patched/`,
`scripts/` untouched; PI applies per README §1 when a prior clears its gate).

## Known-answer reproduction (own code, `test366.py::test_known_answer`, numpy only)

| quantity | reproduced | recorded | match |
|---|---|---|---|
| 323 raw `pds/mse/jac/nmae/fid/reach` | 0.88963585 / 0.97494624 / 0.07475983 / 0.95083646 / 0.49755452 / 0.25434726 | DESIGN 365 §0 table (0.88964 / 0.97495 / 0.07476 / 0.95084 / 0.49755 / 0.25435) | all within 5e-6 |
| 361 partition on `effects.npz::cap15` | block size exactly 2200 on 85/85 targets | 361 RESULTS "block calls/target exactly 2200" | exact |
| universe→effects mapping | misses exactly {HSPA14-1, TBCE-1} | 361 "2 universe genes absent from the effects axis" | exact set |

If any of these had failed, work would have stopped per PLAN.md F-callout. None did.

## What was measured (inputs, estimator, controls)

- Input: `scripts/make_submission_pm_cd4.py` at HEAD `1318cd3c…` (6348 lines).
- Estimator: `patch366.py` — four exact-anchor replacements (E1 env+loader, E2
  `apply_ctx_a_prior` + `build_ctx_a_context_effects`, E3 run() block at the knee tail,
  E4 five config keys), applied to a byte copy at `work/patched/`. `--check` verifies
  each anchor occurs exactly once in HEAD bytes and prints pre/post sha256.
- Controls: 254-pattern inert-flag discipline (unset ⇒ `None`, helper returns
  `(None, None)`); synthetic frame 3 ctx × 5 targets × 50 genes with independently
  recomputed expectations (separate argsort code path in the test); 9 malformed
  artifacts + knee combo + 4 apply-time refusals; `CTX_A_PRIOR` absent from `scripts/`.

## Results table (units: checks, all must pass)

| group | checks | result |
|---|---|---|
| known-answer (323 members + 361 block) | 9 | 9 pass |
| unset ⇒ inert | 3 | 3 pass |
| sign_block exact (5 head-bitwise + full-matrix + audit + B same-object + C bitwise + A differs) | 10 | 10 pass |
| sign_block p_down/tau exact | 1 | 1 pass |
| add_vector exact + scale audit + B/C bitwise | 3 | 3 pass |
| malformed refused at load (missing key, NaN, mode, context, extra key, p_down, tau, alpha=0, shape, path) | 10 | 10 pass |
| import-time refusals (knee combo, bad mode) | 2 | 2 pass |
| apply-time refusals (panel support, gene miss, vacuous add, wrong context) | 4 | 4 pass |
| scripts/ untouched + patched copy carries 4 anchors | 3 | 3 pass |
| **total** | **45** | **45 pass, 0 fail** |

Patch arithmetic: HEAD `1318cd3c…` → patched `95913dc0…`, **+155 lines** (mandate said
≤ ~150; 3% over — the excess is the audit dict and the four refusal branches, kept
deliberately; no further trim without losing a tested refusal).

## Verdict against the PLAN.md kill criterion

**Not killed — ships as a staged patch.** F2 required: unset import clean ✓; B/C bytes
identical ✓ (B is the same object); every malformed artifact refused ✓; head entries
bitwise under sign_block ✓; +155 lines, inside the 20% tolerance band ✓. No F-falsifier
fired. This is engineering, not evidence for any member: no emission, no scorer, no
A/B/C forecast is claimed.

## Limits

- The patched copy has never been through a real build: byte-inertness at staging
  (`effects.final_sha256` with flag unset) and the rule-10 five-key diff are proven only
  at the unit level; the PI's step (d) build must re-prove 254 KA-1 before any flagged run.
- `add_vector` scale `s` is per-row over `b_ctrl > 0` genes (mandate's wording), not 365's
  pooled median — a flagged build's `effects.ctx_a.scale_*` audit records which ran.
- `replace` mode is not implemented (365 C1: vacuous on A/B/C; the loader refuses it).
- `assert_build_carries_arms` has no `ctx_a_prior` row: a flagged build passes the driver
  silently — the PI must add the row when this graduates (README §3).
- H1-panel numbers above are development evidence; nothing here forecasts A/B/C.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50% (was 74%)
uv run python experiments/366-ctx-a-prior-patch/patch366.py --check
uv run python experiments/366-ctx-a-prior-patch/patch366.py
uv run python experiments/366-ctx-a-prior-patch/test366.py   # 45 passed, 0 failed
```

## Peak memory and wall time

`test366.py`: peak RSS **240 MB** (`/usr/bin/time -l` max resident 251,478,016 B), wall
**~1.2 s**. `patch366.py`: < 1 s, < 100 MB. No process exceeded 4 GiB; no scratch files kept (synthetic npz live in a
`TemporaryDirectory`); no downloads, no builds, no submits.
