# 534 — resume mandate: finish Part B if needed, then Part C

You are the same worker that ran Part B, resumed. Nobody can message you mid-run. If something blocks you, write
`experiments/534-generalization-audit/work/partC/BLOCKED.md` with the reason and stop.

## 0. First, Part B

If `experiments/534-generalization-audit/RESULTS-B.md` does not exist yet, finish Part B first: follow
`work/partB/RESUME.md`, under MANDATE-B's rules and windows. Then continue here. If RESULTS-B.md exists, go straight to
Part C.

## 1. Part C: execute `experiments/534-generalization-audit/PREREGISTRATION-C.md` exactly

Question: do three small steps inside every October build's base hold on the October-like panels? The steps are
normrel −1.5, its cap form, and the CD4 partial rows at g 0.5. The preregistration is binding: designs B0–B3 at the
October bulk configuration (`h2c1`, `V1_5_3`, q 1.75); panels `abc`, `real`, `efp`, `efr` with their existing
held-out operators, read only; the known answers in its §3 first; the reading rules in its §4.

- **Reserve R6 is never used in Part C.** Its access rule allows C10a and C10b only.
- No operator fits are needed. Reuse your redirected lib438 copy. Part C's caches and outputs go under
  `work/partC/`.
- The normrel and partial toggles: find where the production values live (`lib418.pre_normrel` via
  `sens415.Donors.pre`; `PARTIAL_G` in lib438; `compose_arm` on the `abc` rules path) and set them in your copy only.
  Prove the toggles are inert at production values (KA) before reading any B1–B3 number.

## 2. What to write

`experiments/534-generalization-audit/RESULTS-C.md`. **§1 must be the decision table the coordinator will cite**, one
row per step:

| step | Δ pooled real/efp/efr (k = 1, board units) | 95 % CI (paired target bootstrap, 4,000 draws, seed 534) | folds positive (of 9) | verdict (HOLDS / FRAGILE / FAILS) | Δ on `abc` | the board's A/B/C read |

Then the per-panel values, the KA table, peak RSS and timing, and one honest sentence per step on what it means for
October. Also write `work/partC/results534C.json`.

The coordinator will add a default-off October knob for any step marked FAILS, citing §1. So the verdicts must follow
§4's rules exactly: no rounding into a verdict, no post-hoc changes to a rule.

## 3. Rules (as in MANDATE-B, unchanged)

- **Write scope:** only under `experiments/534-generalization-audit/` (`work/partC/`, `RESULTS-C.md`, and Part B's
  files if you are finishing B). Never edit `RESULTS.md` (Part A, the owner's), a `PREREGISTRATION-*.md` or a
  `MANDATE-*.md`. No git, `vcc`, `cell-eval2`, builder edits or `submissions/`.
- **Memory gate:** `memory_pressure` ≥ 40 % free, `vm_stat` free+inactive ≥ 27 GiB, and no builder or chain running
  (the same checks as MANDATE-B). 27 GiB, not PREREGISTRATION-B §9's 33, is the owner's amendment of 2026-09-29T15:43Z:
  33 GiB was unreachable with the desktop apps open. Read the floor from `work/partB/GATE-OVERRIDE.json` at every
  poll, exactly as `gate534.py` does. The owner's `swapguard534.sh` (running until 2026-10-03T00:00Z) raises it to
  33 while swap free < 400 MiB (AGENTS §9), which is how the swap condition reaches your gate. Record the amendment
  as the owner's deviation in RESULTS-C.md. Heavy processes strictly serial, each under `/usr/bin/time -l`. Log starts and ends
  to `work/partC/heavy.log`.
- **Windows:** no heavy process starts during 2026-09-29 23:45Z → 2026-09-30 09:45Z, or 2026-10-05 02:00–06:00Z.
  The first window covers the coordinator's build535.sh at 09-30 06:00–~09:30Z: ~20 GB builds whose steps are
  launched by a bash script your chain regex does not name. Also treat any bash process running `build5NN*.sh` as a
  chain. This applies whether you are finishing Part B or running Part C. Stop and write `work/partC/RESUME.md` if you
  would otherwise run into a window. Target: done before 2026-10-03 00:00Z.
- **Names:** `*534*` only. No `chain5NN`, `run_ka_` or `def_day` in any script name or long-lived argv.
- Every number in RESULTS-C.md comes from a log or json in `work/partC/`.
