# 217 — BUILD `cd4f080` and verify every falsifier. DO NOT SUBMIT.

You are a build subagent. The arm is preregistered at `experiments/217-cd4-f080/PREREGISTRATION.md`
— **read it first; it is binding and must not be edited.** **NEVER run `vcc submit`.**

Work in `experiments/217-cd4-f080/` plus `submissions/`. **No code change is needed — one env var.**
If you think a code change is required, stop and say so. `uv run python`.
**Peak RSS ≤ 22 GB, builds strictly SERIAL** (210 measured 19.8 GiB). Two research agents are
running; check free memory before the emission stage and wait if under 25 GB.

## READ FIRST
1. `experiments/217-cd4-f080/PREREGISTRATION.md` — binding, including its five falsifiers.
2. `experiments/210-cd4-build/` — `MANDATE.md`, `RESULTS.md`, **`BOARD-RESULT.md`**. You are
   repeating that build with one number changed; **its `build210.sh` is your script.**
3. `AGENTS.md` **§8.4** (why the dose is being pushed), §6.1, §11.1 rule 10, §14.
4. `experiments/207-cd4-dose/RESULTS.md` §8 — the sweep, and its meta-assertion table.

## The change
Copy `build210.sh`. Change **one line**: `VCC_CD4_TARGET_F=0.34` → `0.80`. Keep
`--emitter rc_callset`, `VCC_ALPHA=0.982` and every other flag byte-identical.
**Set `f`, never `fraction`.** **Do NOT re-solve α.**

## Known-answer row (mandatory, abort on failure)
Rebuild the incumbent first at `VCC_CD4_TARGET_F=0.34` and require sha256
**`7036e89f8a03c82d69e79a5c6690b275adae52dbb3ff72bb6ab23e7ca622bf7a`** on **12,648,807,800 bytes**,
verified with `shasum` on the file. If it does not reproduce, stop.

## Falsifiers — report each with its measured number
1. Preregistration falsifier 5 (config diff exactly 2 fields + the pinned meta values) — fully
   checkable now, and it decides whether this is one variable.
2. **207 §8.2's meta table** at the new dose: `target_f_requested` **0.80**, `fraction_source`
   `VCC_CD4_TARGET_F`, `solved_fraction` ≈ **0.758** (`= 0.80 / 1.0556645369046553`),
   `fully_supported_targets` **276**, `alpha` **0.982**.
3. **The profile moved, and by how much.** Per-gene panel-total digest against
   `210-cd4f034.h5ad` must **differ on 900/900**. Report delivered median row norm and `A` against
   `cd4f034`'s (210 measured 0.99367× and 0.97970× against `rc_callset`) — **the amplitude drift is
   the thing most likely to surprise at this dose, so measure it, do not assume it.**
4. **Non-vacuity**: `‖ΔE‖/‖E‖` against the `f = 0.34` artifact, with a per-block minimum and a count
   of zero blocks. 210 measured 0.235141 against `f = 0.10`; at this dose expect larger. **200 and
   207 both found the CD4 residual can be identically zero while every other check passes.**
5. **Challenge `n_pred` within ±2%** of `cd4f034`'s (probe against `data/challenge/context_*.h5ad`;
   `rc_callset` read 2,599). A large move would mean the call set shifted and the DE members are in
   play beyond what the preregistration predicts.
6. **`vcc prep --dry-run` clean**; nnz against the 4.75e9 cap and max counts/cell against 1e6.

## Deliverable
`RESULTS.md`: the known-answer row, every falsifier with its number, the config diff verbatim, the
package path and sha256, and a one-line **BUILDABLE / NOT BUILDABLE** verdict. Then stop.
**Do not submit. Do not choose a model name or write a description.**
