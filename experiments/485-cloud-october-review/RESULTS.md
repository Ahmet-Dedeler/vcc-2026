# 485 — October day-of review

The day script can ship the wrong build, and it can throw away a safe one. Both are reproduced. The board is dark on Oct 22, so neither is caught by a score.

Apply `patches/def_pipeline.py.patch` and `patches/def_day.sh.patch` together (see `patches/README.md`). The helper and re-level patches go with the shell. Eight findings, all CONFIRMED. Checked on GNU bash 5.2.21; bash 3.2 is not installed on this VM, and the note under each shell finding says what was executed here versus only scanned.

| Rank | Where | What October does |
|---|---|---|
| F1 critical | `scripts/def_day.sh:99-101`, `:113-116`; `scripts/def_pipeline.py:696-704` | A retry ships the previous h5ad (old q, old factors) and records the new factor string. A dump retry can keep the previous npz. |
| F2 high | `scripts/def_day.sh:142-143`; `experiments/448-def-rehearsal-q175/relevel_iter.py:65-71` | A factor outside [0.7, 1.4], or a flat secant, aborts before the safe-build fallback. lvl0 is left on disk and is not linked as `_final`. |
| F3 medium | `scripts/def_day.sh:134-140` | A combination build is labelled CONVERGED whenever the file exists and no rho is above 1, including rho 0.999 (outside the 1% band). |
| F4 medium | `scripts/def_day_helpers.py:25-27`; `scripts/def_day.sh:122-123` | `past_cliff` exits 1 both when every rho is ≤ 1 and when the json is unreadable. The shell keeps the h5ad in both cases. |
| F5 medium | `scripts/def_day.sh:173`, `:181-183` | A non-zero zero-gene count and a failed `vcc prep --dry-run` are logged and the script still prints DONE. |
| F7 medium | `scripts/def_day.sh:44-46`, `:183` | Success says "commit the pins" and the exit trap then restores the builder, so `git diff` is empty. A failed restore still exits 0. |
| F8 medium | `scripts/def_day.sh:112`, `:117`; `scripts/def_pipeline.py:525` | CD4's 275/290 gates have no pass-through, and the factor string and `ntclevel --contexts` are hardcoded `D,E,F`. Both stop the day loudly; they do not ship a wrong file. |
| F6 low | `scripts/def_day.sh:72-78` | If `nctx` cannot read the same-line npz, the script turns same-line mode off and continues. |

Reproductions, from the repo root:

```bash
uv run python experiments/485-cloud-october-review/repro/repro_october_bugs.py
uv run python experiments/485-cloud-october-review/repro/repro_fixes.py
```

The first exits 0 when the bugs behave as claimed. The second exits 0 when the patched copies close them. This run's logs are `repro/last-bugs.log` and `repro/last-fixes.log` (both end `all ... matched the claims`).

## F1 — a non-empty file is treated as a successful build (CONFIRMED)

`build()` runs the driver in a subshell and then only tests `[ -s "$W/$N-$tag.h5ad" ]`. The driver's exit code is discarded. The driver itself refuses to overwrite (`def_pipeline.py:696-697`, `SystemExit`). So a second invocation in the same work dir exits 1 immediately, the shell measures the h5ad that is already there, and `POINTS` records the factor string of the invocation that did not build.

The die text at `def_day.sh:142-143` ("lower q and rerun") happens before the cleanup `rm` at line 166, and it tells the operator to rerun. `chain452.sh` deletes h5ads only after a success, and it uses a fresh work dir, which is why 452's passing rehearsal did not hit this. A rerun after a stop, in the same work dir, does.

The arm guard makes the same hole on a first build. `assert_build_carries_arms` (`def_pipeline.py:703-704`) runs after the h5ad has been published. Its own message says a dropped arm still builds and `vcc prep` still passes. The shell throws that non-zero exit away because the file exists.

The dump is the same shape, and it fires on the first run. `make_submission_pm_cd4.py:6597-6603` writes the npz and `raise SystemExit(0)` before `write_submission`, which is the only meta writer. `r.run` therefore returns, and `assert_build_carries_arms` raises `FileNotFoundError` on a meta path that was never created. Every successful dump makes the driver exit 1. The shell only checks `[ -s "$DUMP" ]`. A later dump that fails before `savez` (memory gate, early crash) leaves the old npz in place, and the shell accepts it. Filenames do not contain q or the LODO flag. `452`'s `driver-dump.log` is not in this checkout; the `FileNotFoundError` was reproduced by calling `assert_build_carries_arms` on a missing meta, and the builder source is the exit-before-meta path. 452's log is consistent with it: setup ends `18:24:05Z`, the dump sha is logged at `18:24:17Z`.

Reproduction (`repro_october_bugs.py`, F1 / F1b / F1c):

```
SHELL_ACCEPTS rc=1 file=BUILT q=1.75 factors=D=1.2000,E=1.2000,F=1.2000  requested=D=1.0500,E=1.0500,F=1.0500
SHIP driver: REFUSING THE BUILD: artifact does not carry the configuration
assert_build_carries_arms on a dump (no meta) raises FileNotFoundError
DUMP_KEPT bytes=17 content=OLD_DUMP_LODO_OFF
```

October: `def_day.final.env` carries the new q from the shell variable. The h5ad meta carries the old q and the old scale. `final.env` and the artifact disagree, and the dark board does not say so.

Fix: `patches/def_pipeline.py.patch` skips the arm assert only when `VCC_BULK_DUMP` is set, no h5ad was published, and the npz is non-empty (it logs that skip). `patches/def_day.sh.patch` deletes that tag's h5ad before the driver, requires rc 0, and dies otherwise. The dump subshell stays, so `VCC_H1_UNION_G` and `VCC_PERT_SCALE` still do not leak; the shell now requires rc 0, a non-empty npz, and a `[bulk-dump] wrote` line in `$W/logs/build.log` (the full stage log, not the driver's 25-line tail). A meta check before the hard link compares `bulk_sharpen_q` and `bulk_ctx_scale` to the factor string and q (absent scale only when every factor is 1).

## F2 — a bad secant step aborts before the safe fallback (CONFIRMED)

The fallback at `def_day.sh:146-163` runs only when the loop finishes with `FINAL` empty. Line 142-143 is `die` when `inbound` rejects a factor, and line 130 is `die` when relevel prints no `VCC_BULK_CTX_SCALE` line. Either exit skips the fallback, so a finished lvl0 is never linked as `_final`.

Three triggers, all reproduced with `experiments/448-def-rehearsal-q175/relevel_iter.py`:

- First step `f = current_f * ntclevel_scale`. With scale 1.55 the script prints `next f 1.5500` and a WARNING, then still prints the `VCC_BULK_CTX_SCALE` line and exits 0. The shell dies: `STOP: next factors D=1.5500,E=1.2000,F=1.2000 leave the builder bound [0.7, 1.4]: lower q and rerun` (rc 1). `FALLBACK_SHIPPED` is absent. Un-levelled rho below about 0.50 pushes the first factor past 1.4 (`rho* / 1.4^2`). 452's lowest was F at 0.591, first factor 1.287, inside the bound, with about 8% of headroom. A quieter line on the day crosses it.
- Flat local exponent. rho 0.80 then 0.801 between f=1 and f=1.15 gives exponent 0.009 and `next f 6730587730.8757`, with the WARNING, and still emits the VCC line. `inbound` then dies.
- Same f, different rho. `log(f1/f0)` is 0. Unpatched relevel raises `ZeroDivisionError` and prints no VCC line. The shell dies with "relevel_iter printed no factors". This needs rho to move between two builds at the same factors; a bit-identical rebuild dedupes via `set()`.

452's own step was replayed from the logged 4-decimal rhos. D's next factor is 1.2069 and F's bracket is 1.2452, matching `experiments/452-def-day-rehearsal/work/efr/def_day.log` lines 29-32. E prints 1.2719 against the log's 1.2718: one ulp from feeding 4-decimal rho into the replay, not a different formula. The patched relevel prints the same stdout as the unpatched one on that input.

The relevel docstring says the loop approaches from below. The code does not clamp. 452 stayed below the cliff because the shell deleted lvl1 (F rho 1.0568), which is the happy path and still works.

October: the day stops with a usable lvl0 on disk and no `_final` link. The recovery string is "lower q and rerun", which in the same work dir is F1.

Fix: `patches/relevel_iter.py.patch` exits 2 with `NO_STEP` and no VCC line when the exponent is non-finite or `|e| < 0.05` (452's exponents are 1.89–2.73; the happy-path stdout is identical). Out-of-bound finite factors are still printed. `patches/def_day.sh.patch` `break`s to `best_safe` on relevel rc ≠ 0, a missing factor line, or `inbound` failure. `best_safe` is the old fallback, with rho* computed from the same q and sd as relevel (`0.31 / (0.31 + 1.78)`, sd 0.01) instead of a second literal 0.9793. Non-finite and rho>1 builds are skipped. The fix check ships lvl0 for the 1.55 proposal and the nearer safe build for the flat pair, and it does not call `build`.

## F3 — a combination build is accepted without its own rho (CONFIRMED)

When relevel says `CONVERGED` but no single build carries every context's factor, the shell builds `lvlC` and, if the h5ad exists, sets `FINAL` and breaks. `past_cliff` deletes only when some rho is above 1. rho 0.999 is under 1 and outside the 1% band of rho* ≈ 0.9793 (the band's upper edge is about 0.989). The script logs `CONVERGED` and the fallback does not run.

Reproduction: a condensed copy of lines 134-140 labels that file `CONVERGED_LABEL lvlC measured rho 0.999 0.999 0.999`.

This is uncommon. The factor string is `.4f`, and every value in [0.7000, 1.4000] round-trips, so `same()` (tolerance 5e-5) finds the build that used that string. lvlC exists for the case where each context converged in a different build. The missing check is still real.

Fix: the shell keeps lvlC only when `inband` says every rho is finite, positive, ≤ 1, and within 1% of rho*. Otherwise it deletes that h5ad and falls through to `best_safe`. `inband` on a 452-lvl3-like rho (0.9795 / 0.9794 / 0.9790) returns 0; on 0.999 it returns 1.

## F4 — an unreadable cliff check keeps the artifact (CONFIRMED)

`past_cliff` returns 0 when any rho is above 1, and 1 otherwise. An exception is an uncaught Python error, which is also exit 1. The shell's `if past_cliff; then rm; fi` keeps the h5ad on exit 1.

Reproduction with the production helper: `safe=1 past=0 crashed=1`, and the shell fragment prints `KEEP safe / DELETE past / KEEP crashed` for a json with no `per_context`.

On its own this does not link `_final`: relevel usually dies later (F2) or another build is chosen. Combined with F1, a retry can ship the undeleted past-cliff file. The trigger is a malformed json; `ntclevel417.py` writes valid json on the success path, so this is narrow.

Fix: the helper returns 2 on a missing, non-finite, or unreadable rho. The shell deletes on 0, keeps on 1, and dies on anything else.

## F5 — the zero-gene audit and the prep dry-run do not stop the script (CONFIRMED)

Line 173 pipes `audit428.py` into `head -1` and never reads the status. Line 181 logs a prep failure as `WARN` and line 183 prints DONE anyway. 452's log line 57 is the clean summary (`0 in 0 blocks; lfc < -5: 7`) and the script continued, which is what we want for the `lfc < -5` count. A non-zero zero count would follow the same path.

Reproduction with a stand-in that prints `24 in 10 blocks` and exits 0, then a prep command that exits 1:

```
pipeline_rc=0
JSON_WRITTEN
WARN prep dry run failed
DONE
```

`pipeline_rc=0` and `JSON_WRITTEN` are the correction to a SIGPIPE theory. Python block-buffers stdout to a pipe, so the sidecar write finishes and the pipeline status is 0. The bug is the missing gate, not a lost json. The summary line is the first line `audit428.py` prints (`experiments/428-zero-gene-audit/audit428.py:64`).

Fix: write the audit to `$W/audit428.txt`, die on rc ≠ 0, and die unless the summary contains `ZERO block total (target excluded): 0 in `. The `lfc < -5` count is logged and not gated (452 shipped with 7; the repair floor is -5). Prep dies on rc ≠ 0.

## F7 — "commit the pins" runs before the trap restores the builder (CONFIRMED)

Line 183 tells the operator to commit the pins. The EXIT trap then runs `def_pin_panel_artifacts.py --restore`. 452's log is the sequence: line 60 DONE, line 61 `restored scripts/make_submission_pm_cd4.py ... (identical)`, line 62 `builder after restore: 1dcb29a3 (want 1dcb29a3)`. After exit, `git diff` on the builder is empty. The pins are regenerable by re-running the pin stages (the 410 key pin, then the LODO pin if it was on). The meta already records the builder hash from build time.

The trap logs both hashes and does not compare them. A restore command that fails, followed by a successful `log`, leaves the script at rc 0. Reproduced: incoming exit 0, restore rc 1, script rc 0, log `builder after restore: bbbb (want aaaa)`.

`exit 7` is not swallowed on this bash when the trap's last command succeeds (reproduced, rc 7). The bug is the success path exiting 0 when the hashes differ, and the message that says to commit a diff the trap has already reverted.

Fix: the trap captures `$?` first, skips restore when `DEF_DAY_KEEP_PINS=1` (so the diff is there to commit), and otherwise exits 1 when the restore command fails or the full sha differs from the pre-pin sha. Default remains restore, which is what 452 checks. An explicit `exit` inside the trap does not re-enter it (three restore invocations for three non-keep runs). The DONE line matches the mode. On bash 5.2, `exit 7` through a successful trap stays 7.

## F8 — CD4 gates and context names stop the day by hand (CONFIRMED)

`def_day.sh` never passes `--cd4-thresholds`. The driver default is `default`. `experiments/047-cd4-direct/extract_cd4.py` raises `ValueError` below 275 targets per condition or 290 in the union. Setup then has no `rc=0` line and the script dies. That is a loud stop, and the recovery on the day is editing the script. Experiment 244's two-invocation workaround is in the runbook body, not in the one-pager command.

`FACT="D=1,E=1,F=1"` and `ntclevel417.py --contexts D,E,F` are literals. Preflight checks `len(contexts)==3` (`scripts/def_preflight.py:189`), not the names. A three-context bundle with other labels gets through setup; lvl0 exports no scale, so the builder's context-key check does not fire; ntclevel then looks for `context_D.h5ad` and dies. It does not auto-ship.

Fix: `DEF_DAY_CD4_THRESHOLDS=relaxed` appends `--cd4-thresholds relaxed` and logs it. Any other non-empty value dies. Unset stays strict. After setup, the script reads `preflight.json`'s `contexts` (three identifier strings) and builds the unit factor string and the ntclevel context list from that. The Python is single-quoted so bash does not brace-expand it.

## F6 — a failed nctx count turns same-line off (CONFIRMED)

`def_sameline.py` is `|| die`. The next line, `nctx`, is not. On rc ≠ 0 the arithmetic test is false and the script logs "same-line mode off" and continues. Reproduced with an object-array npz: numpy raises `allow_pickle`, stdout is empty, the fragment prints `MODE_OFF` and `continued`.

The writer saves `contexts` as a unicode array, and `nctx` reads that (a unicode array of length 2 prints 2). An object array or a corrupt npz is the trigger. Same-line is the H1-priced lever; dropping it silently is a different submission. With `DEF_DAY_PROGRAMME_SOURCE` set, `nctx == 0` still continues unless `DEF_DAY_PRIOR_GATE` is also set (that path dies later at line 106 because `SL` is empty).

Fix: die if `nctx` fails or prints a non-integer. Die if a programme source is set and the count is 0.

## Checked, and not a bug on this path

- The dump subshell does not leak `VCC_H1_UNION_G`, `VCC_PERT_SCALE`, or `VCC_BULK_DUMP` (reproduced: set inside, unset outside).
- `def_day.sh` has no `\|` grep, so the BSD-grep `$`-before-`\|` bug is not in this file. Setup's `rc=0` grep is anchored `^--- \[repin-panel-support\]`. The factor extraction `grep -o 'VCC_BULK_CTX_SCALE=[A-Z0-9=.,]*'` matched 452's positive decimals on the Mac (the log shows the factors that were built). A leading `-` is outside the class; the patched shell treats an empty capture as "no factors" and falls back.
- `"${POINTS[@]}"` is expanded only after lvl0 has appended. The bash 3.2 `set -u` error on an empty array is not reached. Reproduced on bash 5.2: an empty expansion under `set -u` exits 0. The patch does not add `"${EXTRA[@]}"` on an empty array; the CD4 flag is `DRV+=(--cd4-thresholds relaxed)` on an array that already has elements.
- bash 3.2 was not executed. The scan of the patched shell found no `declare -A` and no `\|`. `bash -n` on 5.2 passes. `local`, arrays, `${var:0:16}`, and `pipefail` are 3.2 features; the patch uses those and not `mapfile`, `|&`, or associative arrays.
- The whitener skip is guarded by `verify_existing_operator` (axis sha). The setup log is truncated with `>`.
- The h5ad publish is `os.replace`. `other_builders` matches `ps` output, not `pgrep -f`.
- `.4f` round-trips on [0.7000, 1.4000], so `same()` finds the build that used the string.
- Omitting ntclevel `--f` on D/E/F is the calibration factor f_c = 1, not the re-level factor. rho is measured B over that claim. 452's lvl1 relative scales are not reused as absolute factors; the secant is. The replay matches the log on D and F.
- Sharpening runs after the dump point, so q is not baked into the dump matrix. `VCC_BULK_LOAD` replaces the effect matrix on the level builds.
- rho* printed as 0.9793 and the shell fallback's old literal 0.9793 differ from the full value 0.979343 by 4.3e-5, about 200× smaller than the 1% band. The patch computes it once, in the helper.
- `--min-free-gib 23` is the long-standing rehearsal setting (the driver default is 26). Not a new bug. `free_gib()` also sums speculative and purgeable pages (`def_pipeline.py:328-335`) and hardcodes a 16384-byte page, which matches the M5. On this Linux VM `vm_stat` is missing, so the function returns 0 and a heavy stage refuses. That is fail-closed, and it is not the October machine.
- `REPO=~/Code/virtual-cell-challenge-2026` is the production path. It is wrong on this VM and correct on the day machine.
- A killed build does not leave a partial final h5ad (`os.replace`). F1 needs a completed previous file.

## Fix checks

`repro_fixes.py` runs the patched copies, not the production scripts.

- 452 4-decimal replay: patched and unpatched `relevel_iter.py` stdout are identical, including `next f 1.2069` and `next f 1.2452`.
- Flat exponent exits 2 with `NO_STEP` and no `VCC_BULK_CTX_SCALE` line. Same-f `ZeroDivisionError` does too.
- `past_cliff` exit codes are safe=1, past=0, unreadable=2. `inband` accepts the 452 lvl3 rhos and rejects 0.999. `best_safe` skips rho>1. `meta_matches` accepts 452's scale at q 1.75, rejects q 1.50, accepts a missing scale when every factor is 1, and rejects a missing scale when a factor is 1.2.
- The extracted `build()` from the patched shell replaces a pre-existing h5ad with what the driver just wrote, and a driver that exits 1 removes the pre-existing file and stops the script. A rho of 1.0568 deletes the artifact and returns.
- The extracted loop, on the 1.55 proposal, logs the bound and ships lvl0. It does not build. On the flat pair, relevel's rc 2 ships the nearer safe build.
- The extracted restore trap: sha mismatch exits 1, a failing restore command exits 1, `exit 7` stays 7, `DEF_DAY_KEEP_PINS=1` does not restore. The trap does not re-enter.
- `patch --dry-run -p1` applies all four patches. Production files were not modified.
