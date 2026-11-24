# 591 — run_arm_generic.sh's prep gate gets a pressure-based branch (UNAPPLIED patch + stubbed tests, ALL PASS); the build gate stays as it is, on purpose

2026-10-01, session "Explore contribution areas [57577c]", at the coordinator's request (its message, 13:3xZ). Nothing
live was edited: `experiments/428-zero-gene-audit/run_arm_generic.sh` is unchanged (sha `c58b95533782c597`). Class:
operations.

## The problem
`gate_ok` opens only at memory_pressure free ≥ 40 % AND vm_stat free+inactive ≥ 23 GiB (plus: no other builder, no
`vcc prep`, disk ≥ 20 GiB). With ~18 GiB of other apps' pages in the compressor, free+inactive sits under 23 GiB while
memory_pressure reads ~60 % free and swap is unused. `wait_gate` then gives up after 240 × 30 s (exit 3); 549's prep was
run by hand on 09-30.

## The change (`patch_gate591.py`)
- `gate_ok <mode>` with two memory branches; every other condition is unchanged:
  - `fi23`: pressure ≥ 40 % AND free+inactive ≥ 23 GiB. Today's gate, both modes.
  - `pressure55`: pressure ≥ 55 % AND swap free ≥ 400 MiB (`sysctl -n vm.swapusage`, M/G/K parsed). **Prep only.**
- **Why not the build too.** Invocation 1 runs `def_pipeline.py … --min-free-gib 23`. Before its `repin-panel-support`
  and `build` stages, `def_pipeline` re-checks free+inactive+speculative+purgeable ≥ 23 GiB and raises
  `SystemExit("… refusing")`, i.e. it fails instead of waiting. A shell gate opened on pressure would start a run that
  def_pipeline then kills: today's 2-hour wait becomes a failed build. Giving the build the same branch means the same
  alternative inside `def_pipeline.py`. That file is a pinned driver (run_arm_generic's `case` list of driver shas, and
  538's applier edits it on 10-05), so it is left to the coordinator. The prep (`vcc prep`, invocation 2) has no second
  check, which is exactly where 549 stalled.
- The gate line keeps its prefix and field order (`gate: free …% free+inactive … other builders … vcc preps … disk …`)
  and appends `swap free … MiB  mode …  memory branch fi23|pressure55|none`.
- Call sites: `wait_gate build` before invocation 1, `wait_gate prep` before `vcc prep`; `wait_gate` passes its mode and
  the exhaustion message names it.
- `--apply` refuses while a bash process runs `run_arm_generic.sh` (`--force` overrides). It writes by atomic rename, so
  a running bash would keep reading the old inode anyway. `--check` writes nothing and runs `bash -n` on the result.
- Diff: 24 changed lines; patched sha `e023aebc5cd17950`.

## Tests (`test_gate591.sh`, macOS /bin/bash 3.2, BSD awk/grep; 34 s; ALL PASS)
The gate functions are cut from the ORIGINAL and the PATCHED script text and run against stub readers on PATH (`memory_pressure`, `vm_stat`, `sysctl`,
`uv` for the other-builders probe, `ps`, `df`).
- **R1–R3, 600 gate calls, 0 failures.** The grid is pressure 39/40/54/55/81 % × free+inactive 22.9/23.0/30.0 GiB × swap
  0 / 399 / 400 MiB / 1.5 GiB × {clean, another builder, a vcc prep, disk 19.9, disk 20.0} × {build, prep}.
  - R1: build mode equals the original's exit code at every point.
  - R2: prep mode equals the original everywhere except the intended branch (pressure ≥ 55, swap ≥ 400 MiB, all other
    conditions clean, fi23 closed). So it is identical wherever free+inactive ≥ 23 and pressure ≥ 40.
  - R3: every gate line keeps the prefix and names its mode and the expected branch.
- **R4.** `bash -n` passes; exactly one `wait_gate build` sits directly before invocation 1 and one `wait_gate prep`
  directly before `vcc prep`; no bare `wait_gate` call remains.
- **R5 (live readers, no stubs).** On this Mac at 13:29Z the line reads: `gate: free 78%  free+inactive 26.2 GiB  other
  builders 0  vcc preps 0  disk 98.1 GiB  swap free 1622 MiB  mode prep  memory branch fi23` (rc 0). So the real sysctl
  output parses.

## To apply (coordinator)
`uv run python experiments/591-prep-gate/patch_gate591.py --check`, then `--apply` when no `run_arm_generic.sh` runs,
then `bash experiments/591-prep-gate/test_gate591.sh`. Before `--apply` the test patches the live text in memory; after
it, the test runs the live (patched) file against `original_run_arm_generic_c58b9553.sh`, a byte copy of today's
script whose sha it checks. The patch refuses a second application. Commit the script with its
new sha. Note for anything that pins `run_arm_generic.sh`'s sha: 476's chain only logged it (`c58b9553…`); nothing found
asserts it.
