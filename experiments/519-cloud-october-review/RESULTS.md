# 519 — October-path review (2026-09-28)

Reviewed the 09-28 axis-row change, 502's pass-through (not yet applied), finish497dr/gate497, the re-armed 0928–1005 slot chains, and `october_recipe.py` with `de_mag` decided. No challenge h5ad, donor table, or harness matrix was read. Nothing here was applied to `scripts/` or to another experiment. `vcc` was not run.

Repros (each exits 0 while the defect is still in the tree):

```text
uv run python experiments/519-cloud-october-review/repro/repro_axis.py
uv run python experiments/519-cloud-october-review/repro/repro_502.py
uv run python experiments/519-cloud-october-review/repro/repro_recipe.py
uv run python experiments/519-cloud-october-review/repro/repro_gate497.py
uv run python experiments/519-cloud-october-review/repro/repro_slots.py
```

All five exited 0 on this checkout. `git apply --check` passed for every patch below; none was applied.

## Findings

### 1. `abc_targets` reads the file the Oct 22 download replaces — blocks October

`scripts/def_postbuild.py:318`

```python
abc_targets = panel_set == {str(t) for t in pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"]}
```

`panel_set` is the bundle's `pert_counts.csv` (`def_postbuild.py:226`). DEF-RUNBOOK copies `data/challenge` aside and then `vcc datasets download controls` overwrites `data/challenge`. The day command's bundle is that tree. After the download, line 318 compares the new target list to itself, so `abc_targets` is true for every D/E/F build whose `--bundle` is `data/challenge`.

Lines 323–326 then require the recomputed column digest and the symbol digest to equal A/B/C's constants. A correct refit whose axis moved is failed. `repro_axis.py` runs lines 318–329 (dedented, `ROOT` pointed at a temp tree, the module's `ABC_*` constants) and gets that false FAIL: `abc_targets=True` and the check named "new contexts on A/B/C's own 300 targets".

The other three shapes, while `data/challenge/pert_counts.csv` is still A/B/C's file:

| Case | What line 318 does |
|---|---|
| New contexts, new targets, axis moved | `abc_targets` false. The "moved off" rows do not fire. No failure from this block. |
| New contexts, A/B/C's own 300 targets, axis is A/B/C's (497dr) | `abc_targets` true. Lines 323–326 pass. This is the row that unblocked 497dr. |
| A/B/C itself | `new_panel` is false. Lines 327–329 require the symbol digest to be A/B/C's. The new branch is not consulted. |

It does not fail open on a wrong operator. Line 285 still requires `operator.axis_sha256` to equal the axis recomputed from this bundle. The repro records A/B/C's symbol digest against a moved recomputation and line 285 fails. When the recorded digest equals the recomputation, line 285 passes. The new rows are absent in that case; they are not a passing substitute for line 285.

Two more fail-closed cases, same block:

- Preserved A/B/C-target bundle after the download (`--bundle data/challenge-abc`, while `data/challenge/pert_counts.csv` is already the new list). `abc_targets` is false and the axis is still A/B/C's, so lines 319–322 fire "IDENTICAL to A/B/C — impossible".
- A panel that shares most A/B/C targets but whose symmetric difference misses the donor union. The axis is `union(donor coverage) − panel`, so that swap leaves the symbol digest equal to A/B/C's. Lines 319–322 call that impossible. The repro shows the set identity on `{IN_UNION}` versus a one-gene swap of `ABCD1` for `NOT_IN_UNION`, then the same false FAIL from the real lines. A partial overlap whose axis did move adds no row.

Lines 323–326 also require `columns_sha256 == ABC_OPERATOR_COLUMNS_SHA`. With the challenge's fixed gene order, the same symbols are the same indices, so this is redundant with the symbol digest. It is the wrong digest to use as a panel identity (the file's own header: index digests collide across panels). The symbol row is the one that carries the answer.

Pinned A/B/C target digest, newline-joined sorted symbols of the current `data/challenge/pert_counts.csv` (300 names): `c9c4c9a69f76afd4507e9a7619fe9925c19ff34e5878c7854493683b23bea5ca`.

Patch: `patches/def_postbuild_abc_targets.patch`. Compare that pin, drop the two "impossible" rows, and on the same-target branch require only the symbol digest. Line 285 stays the binding check. Not applied.

### 2. A failed 502 rider build exits 0 and can leave `def_day.sh` half-patched — blocks October

`experiments/502-october-de-riders/chain502.sh:59-63` and `patch_defday502.py:89`.

`(b)` failure logs and falls through. `set -uo pipefail` does not include `-e`, and the next command is `log "END"`, so the script exits 0. `repro_502.py` runs that tail under bash and gets rc 0. `scripts/commit_after_chain.sh:41-43` then looks for `(b) PASS`, does not find it, and does not commit. The working tree is left at the patched shas. A re-arm's release gate (`chain502.sh:23`) requires `def_day.sh` prefix `e29e793a7964df03`, which no longer matches, so it waits until `202610031200` and exits 1 without restoring.

`patch_defday502.py:89` writes `def_day.sh` and then the helpers. A kill between those two calls leaves one file new. The same release gate does not match, and `back()` never runs because the process is already dead. The repro writes the post-502 script and leaves the helpers old: the prefix is `610db78277239f50`, not `e29e793a7964df03`.

Line 49 (memory never reached 30 GiB before `(b)`) already exits 1, and it also does not restore. Same stuck tree, with a visible rc.

The commit regex is not the bug. `commit502.spec:7-8` is `REQ=\(a\) PASS` and `REQ=\(b\) PASS`. `read -r` keeps the backslashes, and `grep -E` of that pattern matches chain502's log lines `(a) PASS:` and `(b) PASS:` (count 1). The unescaped form `(a) PASS` matches nothing (count 0), because the parentheses are a capturing group around `a`. A successful `(a)` and `(b)` would commit, and only those two files: `EXPECT` is `610db78277239f50` and `20142142e4574b7f`, which are the sha prefixes of `repro/def_day.sh.new` and `repro/def_day_helpers.py.new`.

Patch: `patches/chain502_fail_closed.patch` (restore a half-applied pair from the pre-502 copies before the wait; `(b)` failure logs `END` and exits 1) and `patches/patch_defday502_atomic_and_riders.patch` (write both temps, then `os.replace`). Not applied. A consistent post-502 tree whose `(b)` failed is still not committed and a re-arm still will not retry it; that is the stop, and the rc is now 1.

### 3. `command` prints a line that drops `de_mag`, and exits nonzero — costs score silently if the line is run

`scripts/october_recipe.py:91` and `:120-125`. `docs/october-decisions.json` has `de_mag` DECIDED `dn_k250_zs` (2026-09-28T13:07:13Z). `scripts/def_day.sh` contains no `DEF_DAY_DE_MAG`.

Today, with no flags:

- stdout: `DEF_DAY_LODO=1 DEF_DAY_DE_MAG=dn_k250_zs bash scripts/def_day.sh data/challenge W stem 1.75 0.5 noprog 4`
- stderr: `PROBLEM: PENDING: lodo_hi_pct, sharpen_pt, programme_source, prior_gate, mde_sign` and `PROBLEM: scripts/def_day.sh has no DEF_DAY_DE_MAG`
- exit 2. Pending is tested before the missing-knob exit 3, so the rc hides the knob. The message does not.

`--allow-pending` prints the same command and exits 3, stderr only the missing knob.

Honoring the rc refuses. Executing the stdout line on today's `def_day.sh` drops the rider: `DEF_DAY_DE_MAG` is not a `VCC_*` name, the wipe at `def_day.sh:66` does not export `VCC_DE_MAG`, and nothing else reads it. The repro's bash fragment prints `VCC_DE_MAG=<unset>`.

`installed_problems` is `var not in dd`. A copy of `def_day.sh` with the token only in a comment, plus `--allow-pending`, exits 0 and still prints the rider. The live builder does contain `VCC_DE_MAG`, so the builder half of the check is satisfied; the shell half is the hole.

`scripts/test_october_recipe.py:49-51` still expects `--allow-pending` to exit 0 and to print `DEF_DAY_LODO=1 bash scripts/def_day.sh ...` with no rider. That test fails: `AssertionError: 3 != 0` and the missing-knob problem. Cosmetic. The recipe itself is what the day runs.

502's pass-through is what makes the printed variable do anything. It is not applied. Landing 502 (finding 2's success path) is the real fix for the day; the patch below stops the checker from being fooled and stops a refused command from being pasteable.

Patch: `patches/october_recipe_refuse_stdout.patch` (require `export VCC_DE_MAG="$DEF_DAY_DE_MAG"`, print `# REFUSED` in front of the line when `problems` is non-empty) and `patches/test_october_recipe_live.patch`. Not applied.

### 4. gate497 does not check what chain497 checked before relabel — costs score silently on a re-gate

`experiments/497-october-switches/chain497.sh:58` requires `DONE. Next` in `def_day.log` before relabel. That line is written only after postbuild, the audit, and the driver's dry-run. `finish497dr.sh` goes straight to relabel (`finish497dr.sh:9`). `work/run/def_day.log` has no `DONE. Next`. The 12:18:05Z lines are a manual continuation, including `postbuild PASS (manual re-run` and an explicit note that the driver's D/E/F prep was not re-run. `finish497dr.sh:17` logs `gate497 rc=$?` and does not exit on it. chain497's own gate call also ignores the rc; chain497 at least had the `DONE. Next` gate first.

`gate497.py` never mentions `postbuild` or `DONE. Next`. The zero-gene check (`gate497.py:46-47`) is a regex on the last log line containing `ZERO block total`. The repro's forged line `someone typed: ZERO block total (target excluded): 0 in 0 blocks` matches. The docstring (`gate497.py:8`) says `gate460.py`. Cosmetic.

This package, as it sits, is not a failed postbuild: `work/run/postbuild.json` has no `"pass": false`, `work/run/audit428.txt` line 1 is the zero-gene line, and `work/gates497.log` is a real PASS (cliff, LODO `c998b18c…`, S_P2, config diff exactly the four bulk keys, novelty rc 0, package `6ed67b2a…`). Requiring the string `DONE. Next` would reject that package. The lock that does not is the json and `audit428.txt`.

`chain1003.sh:40`, `chain1004.sh:58`, and `chain1005.sh:60` run `gate497.py --write` and then, as the next statement (`set -e` is off), grep the existing `gates497.log`. `--write` exits 1 without rewriting the log when a gate fails. A stale PASS file still matches. The repro shows a `false` followed by those greps printing `STALE_PASS`, and `false || return 1` printing `REJECTED`.

Patches: `patches/gate497_require_postbuild.patch`, `patches/finish497dr_postbuild_json.patch`, `patches/stale_gate_chain1003.sh.patch`, `patches/stale_gate_chain1004.sh.patch`, `patches/stale_gate_chain1005.sh.patch`. Not applied.

### 5. The 12 h-late 09-28 submits stayed on 09-28. One later chain can release on a stale `END` — latent, costs the 10-01 arm

`arm0928.sh:20-23`: slot 1's 12:00Z cap sits inside `until date >= 202609280001`. A start after 00:01Z never enters the body, so the cap does not run. That is why the 12:11Z re-arm submitted. Slot 2's cap edit is `arm0928.sh:28` (`202609282330`), also inside its until, and there is no post-loop date check (0929 and later have one; `arm0929.sh:33`). Slot 2 submits whenever `CAPTURED own entry` is already in slot 1's log, including after 23:30.

The logs:

- slot 1 re-armed `2026-09-28T12:11:05Z`, SUBMITTING `12:11:29Z`, SUBMIT rc 0 `12:37:12Z`, CAPTURED `q3VuTImh4FRGKXH2waCp` `12:37:36Z`, END.
- slot 2 re-armed `12:11:07Z`, SUBMITTING `12:38:06Z`, SUBMIT rc 0 `13:04:58Z`, CAPTURED `be90smZUDcIV7j5MgQ03` `13:05:17Z`, END.

Both are 09-28 UTC. `scripts/slot_entry_id.py` on those logs returns those two ids. A synthetic log whose last `ARMED slot 1: 453dr` has no following capture returns `none`, so the killed 09-27 arm (pid 72218 / 72220) is not reused.

`chain0929.sh:12-15` releases once the clock is past `202609280005`, `chain0928.log` contains `armed 09-28` (it does, from 09-27, when the chain armed the waiters and exited), and no bash process matches `arm0928.sh|chain0929`'s awk. After the waiters exited at 13:05Z it can arm `arm0929.sh` the same afternoon. `arm0929.sh` sleeps until `202609290001` and the post-loop check refuses a start at or after `202609291200`. The packages it arms are 457dr and 446bd2, not 453dr or 446bj2z, and it re-runs novelty before submit. The 18:00Z break can start `arm0929` while slot 2 is still inside the 23:30 wait. Slot 2 had already captured at 13:05Z, so that overlap did not happen. A capture landing in the last seconds before 23:30 could still start `vcc submit --wait` with no post-loop check and overlap 09-29 00:01Z. One in flight, not a second submit of the same package.

`chain1002.sh:14`, `chain1003.sh:14`, `chain1004.sh:13`, `chain1005.sh:13` wait on `202610012000`, `202610022000`, `202610032000`, `202610042000`. None of them can fire on 09-28. `arm0930.sh:33`, `arm1001.sh`, `arm1002.sh:41`, `arm1003.sh:41`, `arm1004.sh:41`, `arm1005.sh:41` all have the post-loop noon stop.

`chain1001b.sh:13` is the latent one. The until matches `CAPTURED own entry|NOT SUBMITTING|exiting| END$` anywhere in `submit0930s2.log`, with no date floor and no process check. There is no `submit0930*.log` in this checkout, so the 09-28 fire does not release it. A one-line synthetic log ending in ` END`, and a line containing `exiting`, both match `grep -qE`. Releasing early makes `rule1001b.py` see no 09-30 ids and arm a fallback for 10-01. `arm1001.sh` then sleeps until `202610010001`, so the UTC day of the submit is still 10-01; the arm can be the wrong one.

Patch: `patches/chain1001b_wait.patch` (date floor `202609300005`, bash-plus-path process check, drop `exiting| END$`). Not applied. No patch for `arm0928.sh`: both slots have already submitted on 09-28, and the missing post-loop check is what let the recovery submit.

### 6. `riders_match` treats a non-dict `de_mag` as absent — cosmetic

`experiments/502-october-de-riders/repro/def_day_helpers.py.new:96` (and the same line in `patch_defday502.py:58`).

`got_mag` is `None` when `de_mag` is missing, not a dict, or a dict whose `variant` is missing (`{} or {}`). The ON path (`dn_k250_zs`) returns 1 for all three: fail closed. The OFF path (`-`) returns 0 for `{}` and for a string. The current builder writes `de_mag` as a dict with `variant` only when `VCC_DE_MAG` is set, so the October ON path is not this hole. A non-numeric `tau` raises `ValueError` outside the json `try` (loud, rc not 0). `clean4:*` in the shell accepts `clean4:foo`; the builder, once 496 is applied, dies on the float. Today the builder has no `VCC_DE_MDE_SIGN`, and the shell dies at `grep -q VCC_DE_MDE_SIGN` before the export.

Patch: `patches/def_day_helpers_riders_match.patch` together with the riders half of `patches/patch_defday502_atomic_and_riders.patch`. Apply both, or `chain502.sh`'s `cmp` against `repro/*.new` fails. Not applied.

## Checked and found correct

- Export order in `repro/def_day.sh.new`: the `VCC_*` wipe is line 66, the dump subshell closes at line 200, `export VCC_DE_MAG` is line 210. A bash copy of that order leaves the dump subshell with `VCC_DE_MAG` empty and the following level env at `dn_k250_zs`, including when `VCC_DE_MAG` was set before the wipe.
- `scripts/def_pipeline.py:511-513`: `VCC_DE_MAG` is not in `SHIPPED_ENV` or `DEF_ENV`, so a level-build child (`{**os.environ, **resolved}`) keeps it. `bare_env` drops every `VCC_*` before the overlay and does not put it back. Pin stages use `bare_env`. The dump subshell is forked before the export, so the dump does not see the rider even though its build is not `bare_env`.
- `ALLOW` gains `de_mag` / `de_mde_sign` only when the corresponding `VCC_*` is set (`def_day.sh.new:319-320`). `def_day.final.env` gains the same lines only when set (`:308-309`). `riders_match` runs either way: unset requires both keys absent.
- Shell `case` for `DEF_DAY_DE_MAG` is `*[!a-z0-9_]*`. `dn_k250_zs` passes it and the builder's `VCC_DE_MAG` regex. A looser token passes the shell and the builder's `SystemExit` fires. Fail closed.
- `chain502.sh:22` matches `bash` plus the path (`awk '$1 ~ /^bash$/'`). On a synthetic process list it matches `def_day.sh` and `chain511.sh` and does not match `chain502.sh` or a python builder. The deadline is a timestamp (`202610031200`), not an iteration cap. `(a)` failure calls `back()`, which logs `RESTORED`; the committer refuses that word (`commit_after_chain.sh:39`). `grep -E` is used, not `$` before `\|`. No `declare -A`.
- The applied-bytes check is already in the chain: `cmp` against `repro/*.new` before the tests, and a mismatch calls `back()`.
- 497dr's recorded gate log matches a package that was prepped (dry-run and real), novelty rc 0, and whose postbuild json has no failing check. The defect is that a later run is not locked to that json.
- Slot chains 0929–1005 do not re-submit 453dr or 446bj2z. `slot_entry_id.py` keys off the last `ARMED` of the tag. 1002–1005 cannot release on 09-28. 0929/0930/1001+ arm scripts refuse a late wake after that day's 12:00Z.
- `chain1001b.sh`'s `grep -qE` is the prescribed form. A local `grep -E ' END$'` anchors (one match on a two-line sample). This VM's grep is not BSD grep; the chains that were read use `grep -qE`, which is the form AGENTS 6.2 requires, so the `$`-before-`\|` bug is not in these waits.

## Patches (not applied)

| File | What it changes |
|---|---|
| `patches/def_postbuild_abc_targets.patch` | Pin `ABC_TARGET_SET_SHA`, drop the impossible rows, symbol-only same-target row |
| `patches/chain502_fail_closed.patch` | Restore a half-applied pair; `(b)` failure exits 1 |
| `patches/patch_defday502_atomic_and_riders.patch` | Temp files then `os.replace`; riders_match rejects a non-dict |
| `patches/def_day_helpers_riders_match.patch` | The same riders_match change in `repro/def_day_helpers.py.new` |
| `patches/october_recipe_refuse_stdout.patch` | Require the export line; comment out a refused command |
| `patches/test_october_recipe_live.patch` | Live-state test expects rc 3 and the missing knob |
| `patches/gate497_require_postbuild.patch` | Read `postbuild.json` and `audit428.txt`; docstring name |
| `patches/finish497dr_postbuild_json.patch` | Refuse to relabel when `postbuild.json` has a failing check |
| `patches/stale_gate_chain1003.sh.patch` | `return 1` when gate497 `--write` fails |
| `patches/stale_gate_chain1004.sh.patch` | same |
| `patches/stale_gate_chain1005.sh.patch` | same |
| `patches/chain1001b_wait.patch` | Date floor, process exit, no whole-file `END$` / `exiting` |

---

## Coordinator's disposition (2026-09-28 13:55Z, on master)
All five repros were run locally first and demonstrate their defects (`repro_axis`, `repro_502`, `repro_recipe`,
`repro_gate497`, `repro_slots`: "ALL DEMONSTRATED").
1. **Applied (critical): `def_postbuild_abc_targets.patch`.** The coordinator's own 09-28 fix read the live
   `data/challenge/pert_counts.csv`, which the Oct 22 download replaces, so every correct D/E/F build would have failed
   postbuild. Now the A/B/C target set is a pinned digest (`c9c4c9a6…`, recomputed independently here) and the
   same-target branch checks the symbol axis only; line 285 remains the binding operator check. Verified: 497dr 42/42;
   verdict lists identical to HEAD on 446bj2/data/challenge (34) and 474efp/ef-bundle (41).
2. **Applied, rewritten:** `chain502.sh` now reverts a partial 502 pair at start (both the pre- and post-502 pairs are
   accepted), and a `(b)` failure or a pre-`(b)` memory stop RESTORES the pre-502 files (and deletes the rehearsal
   h5ads) instead of leaving an applied-but-unproven pass-through. `patch_defday502.py` writes both temps, then
   `os.replace`: output bytes unchanged (== repro/*.new, 0644, test_riders502 9/9). The `riders_match` content change
   (finding 6) was NOT applied: it is cosmetic (the ON path already fails closed) and would break the tested byte check.
   chain502 and its committer were killed and re-armed (13:50Z).
3. **Applied:** `october_recipe.py` checks for the export line (not the token) and prints a refused command as
   `# REFUSED …`. The test was rewritten to accept both states (before 502: rc 3 + the problem; after: rc 0 with the
   rider) instead of the patch's today-only expectation. 6/6 pass.
4. **Applied:** `gate497.py` refuses on a missing/failing `postbuild.json` and reads the zero-gene audit from
   `audit428.txt` (negative control: the 20:20Z failing json gives rc 1; today's passing json rc 0). Chains 1003/1004/1005
   return 1 when `gate497.py --write` fails; all three killed and re-armed (13:51Z). `finish497dr.sh` is a finished
   one-off and was not patched.
5. **Applied:** `chain1001b.sh` waits for 09-30 00:05Z and for the 09-30 waiters' processes to exit; the whole-file
   `END`/`exiting` matches are gone. Killed and re-armed (13:52Z).
