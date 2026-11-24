# 535 — 528's call-propensity block in production, and the arm 446bj2zcp (results)

Preregistration: `PREREGISTRATION.md` (2026-09-29 16:45Z; slot amended to 10-01 s2 at 09-30 06:38Z, bins unchanged).

## 1. Production change (committed `06120021`, pushed)
`build535.sh`, run 2 (run 1 failed closed at 06:00Z: test_patch528 re-patches the production emitter and found it already
swapped in; the script now runs 528's tests before either file is swapped).

| step | result |
|---|---|
| emitter | `src/vcc2026/rc_emit_k.py` 1658ae2c → **0f56209d** (== 528's tested `k_emit528.py`) |
| builder | `scripts/make_submission_pm_cd4.py` 39c1dc66 → **4daae237** (pins the new emitter) |
| 528's tests | test_emitter528v2 26/26, test_emitter528 23/23, test_patch528 22/22 |
| builder suite | 56 passed, 0 failed, 1 skipped |
| KA535a | 446bj2z's env, knob unset → `ec3a5ae6` byte for byte |
| KA535b | 446bj2z's env + the INERT artifact (`e62ee765`, CPM-order key) → `ec3a5ae6` byte for byte |
| KA535c | 446bj2zmz's env (Arc hook), knob unset → `ea5aaa92` byte for byte |

The committed message first said "builder 2973c48b": that is test_patch528's ROOT-pinned test copy. The production file is
4daae237; the unpushed commit was amended to say so (09-30 10:48Z).

## 2. The arm 446bj2zcp
446bj2z's env + `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ec…`.

First build (09-30 08:04–08:37Z, `428/work/build_446bj2zcp.nopkg-0930/`): config diff vs 446bj2z exactly
`['de_sign_emitter_sha256', 'rc_callset_block']` (PASS), zero-gene audit 0 genes in 0 blocks, CLIFF OK; artifact
`callprop535.h5ad` 11,969,322,168 bytes, sha256 `0b3733d5e6e0b2f0…`.

## 3. Incident (09-30): the package was never made, then the prediction was deleted
- `run_arm_generic.sh`'s prep gate needs ≥ 20 GiB free disk after the build. The data volume had **17.7 GiB free** (99 %
  full: 1.8 TiB used), so the gate polled for two hours and logged `GATE never opened` at 10:37:18Z. No `vcc prep`, no
  package, no novelty line.
- `build535.sh`'s last line then ran `rm -f …/callprop535.h5ad` unconditionally. It was meant to drop the source once the
  `.vcc` existed (AGENTS §10). With no package, it deleted the only copy.
- `arm1001cp.sh` then found no package and armed no slot 2 (10:38Z), as designed.

Fix: (1) freed 52 GiB under AGENTS §10 (`git clean -xdf` on 417's 09-23 build/dump trees, 428a/g/h and 399's failed
build; production's work-level dump `.npz` files unchanged, `bb19f0e5` / `00ccfcba`); (2) `rebuild535.sh` rebuilds the arm
with build535's env on the committed builder/emitter, requires the new artifact to equal `0b3733d5` byte for byte (else it
writes `464/work/NOGO-1001-s2`), and leaves the h5ad in place until the arm is scored; (3) `arm1001cp.sh` now takes the pid
file and log of the build it waits for, re-armed on `work/rebuild535.pid`. Rules added to AGENTS §10 (disk headroom; a
cleanup step must be conditional on the package).

## 4. Rebuild (09-30 10:51–11:32Z)
- `rebuild535.sh` on builder 4daae237 / emitter 0f56209d (== HEAD), build535's env, disk 75.9 GiB at release.
- **The rebuilt artifact equals the first build's `0b3733d5e6e0b2f0…` byte for byte** (11,969,322,168 bytes): the build is
  deterministic across runs, and the package is the arm that passed the gates at 08:36Z.
- Gates again: `[446bj2zcp] PASS: exactly ['de_sign_emitter_sha256', 'rc_callset_block']`, zero-gene audit 0, CLIFF OK
  (identical to 446bj2z's bulk: a DE-only knob), novelty rc 0.
- Package `submissions/446bj2zcp-…-zs250-cp.vcc`, sha256 `fa556a05477048e85f5418bf2c2545b407160634d5fe37b606113c89cfc449ea`.
- `arm1001cp.sh` armed **10-01 slot 2** at 11:32:53Z (`464/arm1001.sh 2 446bj2zcp …`, pid 74867, name #23 "I see problems and
  I solve them.", description empty); it waits for slot 1's (460dr) capture. The source h5ad stays until the arm is scored.
- Reader: `experiments/464-slot-1001/read1001.py cp <entry>` (written 11:11Z, before any 10-01 entry).

## 5. Board read
10-01 slot 2, after 460dr's capture. Bins in `PREREGISTRATION.md`.
