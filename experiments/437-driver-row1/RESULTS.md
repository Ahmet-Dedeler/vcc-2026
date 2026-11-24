# 437 — 432's October driver installed and proved: its defaults rebuild 433h byte for byte (row 1 PASS); two more driver traps fixed; T6 guarded in the builder

Status: 2026-09-23 16:40Z, PI.

## What was installed (from `experiments/432-def-rehearsal-new-knobs/proposals/`, after the 433g/433z/433j builds)
| file | before | after | what |
|---|---|---|---|
| `scripts/def_pipeline.py` | `27edc486` | **`5463fe28`** | 432's proposal (details below) plus two 437 fixes |
| `scripts/def_pin_panel_artifacts.py` | `2aab299c` | **`10af47e7`** | 432's `--pert-scale <npz> [--pert-scale-id]` (T2) |
| `scripts/def_pert_scale_key.py` | — | **`0e375b77`** | 432's `key432.py` with `ROOT` fixed for `scripts/` (T1) |
| `scripts/make_submission_pm_cd4.py` | `2c90d9e1` | **`63dfd7d0`** | T6 guard (below) |

432's driver proposal adds:
- `SHIPPED_ENV` with every September knob, and `DEF_ENV`/`ABC_ENV` keys;
- `september_rows` arms guard;
- postbuild allow-list;
- `key` stage.

`experiments/428-zero-gene-audit/run_arm_generic.sh` now pins driver `5463fe28`. Note that `unset` now means
"driver default": a knob an arm wants OFF must be exported empty.

## Known answers
- **Proposal tests before install**, reproducing the worker's results:
  - `test_arms432.py`: 10/10 PASS.
  - `test_pin432.py`: 7/7 PASS.
  - Their record files in 432's `work/` were restored from git afterwards, because the tests write there.
- **Installed-driver arms test** (`test_arms437.py` = 432's test pointed at `scripts/def_pipeline.py`): **9/10**. The
  one "failure" is T1, "the CURRENT guard accepts cd4svd5". It is expected: the current guard is now the September
  one, which refuses cd4svd5 exactly as P3 requires.
- **Key script KA:** `scripts/def_pert_scale_key.py --bundle data/challenge --op data/donors/whitener_union_pruned.npz
  --ka` reproduces the pinned A/B/C key **`d7f76a22…` byte for byte** (6 s).
- **ROW 1 (the point of this experiment).** `run437.sh` exports ONLY the two day-of values: 433h's stack dump
  (`VCC_BULK_LOAD` `bb19f0e5…`) and its re-level (`VCC_BULK_CTX_SCALE=A=1.0088,B=1.05,C=1.0736`). Everything else
  comes from the driver's defaults. It runs the A/B/C stage list through `postbuild`, and the artifact is
  **`5d9e2e380bcff3dc…`, 433h's own artifact sha: PASS**.
  - So the new `SHIPPED_ENV` + `ABC_ENV` is exactly 433h's configuration, and October's defaults cannot silently drop a
    September knob: 432's T4 is closed.
  - Row 1 ran on driver `48b20c4c`. The two later edits (below) change only the pin stages' environment and where the
    restore runs, not the build.
- **T6 guard tests** (`test_t6.py`, import-level): **6/6 PASS**. A/B/C unset, A/B/C with F 0.51 (417's A/B/C dump
  recipe), D/E/F carrying + G, + the driver's F 0.34, CD4 ×1.5 via G: all accepted. CD4 ×1.5 via F under carrying:
  refused.
  - The builder change is import-time only. It raises on (carrying ∧ F ≠ 0.34) and never fires on A/B/C, where sizing
    is unset, so A/B/C builds are unchanged by construction.
  - `test_coarse426`, `test_mde426`, `test_zr428` and `test_zs429` all PASS on `63dfd7d0`.

## Two more driver traps, found by row 1's first attempt and fixed here
- **T3b: a day-of export whose partner is a driver DEFAULT makes the pin stages refuse.**
  - Try 1 exported LOAD/CTX from the start. `repin-frozen` runs the builder with the exports but not the defaults, so
    it saw `VCC_BULK_CTX_SCALE` without `VCC_PERT_SCALE_BALLAST` and refused at import
    (`work/run437.try1-repinfrozen-refused.log`).
  - Fix: `Runner.run(..., bare_env=True)` for `repin-frozen` and `pin-artifacts`. Every `VCC_*` is removed there,
    which is exactly 432's validated sequence of exporting nothing until after the pins, but now whatever the operator
    exported. The log line says `[bare env: every VCC_* removed]`.
- **T10: `--restore-builder-after` only restored inside the `prep` stage.**
  - Any run whose stage list stops before prep (every arm build) or that fails mid-way left the production builder
    PINNED. Row 1 left it at `b750b3c0`, and `run437.sh`'s own explicit restore caught it.
  - Fix: the stage loop is wrapped in `try/finally`, and the restore runs there, but only if the builder on disk still
    equals THIS pin's `builder_sha256_after`, so a stale record never reverts a later committed builder.
  - Tested:
    - (b) a census-only run next to the old record says "not restoring", builder `2c90d9e1` kept;
    - (c) a real `pin-artifacts` run is restored in `finally`, `b750b3c0 → 2c90d9e1` identical.

## Still open for October (432 §6)
- T7: `cliffcheck414.py` is A/B/C-only. Use 432 §4's `ntclevel417` pair on D/E/F.
- T8: LOAD does not check the dump's operator or key.
- T9: synthetic D/E/F are A/B/C's controls relabelled.
- The `SHIPPED_ENV` now carries `VCC_DE_MDE=0.1,0.1,0.3,coarse` before its board read (433a on 09-24). If 433a reads
  "no transfer", remove it from `SHIPPED_ENV` and re-run row 1 against the new shipped arm.
