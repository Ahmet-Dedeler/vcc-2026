# 663 — idall3 build, the 10-05 slot-2 fill

Coordinator, 2026-10-04. idall3 = 498p2g (the best scored entry, JWSghc0o7MN5D6o8Iadc, 0.2749190) + the UNGATED v3 identity
hook (`VCC_DE_MDE_SIGN=file:<604 v3 artifact, λ 0.5>:6cd84dc2…:0`, no KEYGATE), per
`experiments/612-key-gated-hook/PREREGISTRATION-draft-IDALL3.md` (FINAL, cf3f8004). Class: transfer (a DE sign prior). It reads against
Δ_496a, which is 10-05's slot 1. 618 predicts a TIE.

## Scripts
- `build_idall3.sh`: runs `run_arm_generic.sh idall3` with 498p2g's round-2 environment. It expects a config diff of exactly
  `de_mde_sign` + `de_sign_emitter_sha256` and ends `READY idall3 <sha>` or `NOT READY <reason>` in `work/build_idall3.log`.
- `cutoff_idall3.sh`: at 22:30Z, kills the build only if the driver never started. The driver did start at 18:39Z, so the
  cutoff is moot.
- `fill1005s2.sh`: when chain1005 leaves slot 2 empty and the build is READY, arms `496/arm1005.sh 2 idall3`. Cap 10-05 04:00Z.
- `prep_fallback663.sh` (armed 19:56Z, pid 34190): see below.

## Build record
- 19:23Z: PASS (config diff exactly as expected); zero-gene audit 0 blocks; CLIFF OK `{"A": 0.9494, "B": 1.025, "C": 1.0908}`,
  equal to 498p2g's (the hook does not touch the bulk).
- From 19:23:58Z, the prep gate stayed closed. The iOS Simulator held ~17.5 GiB in the compressor, leaving free+inactive at
  18.6–19.7 GiB (fi23 needs 23) and swap free at 386 MiB (pressure55 needs 400, at memory_pressure free 58–60 %).
  `wait_gate` gives up after 240 × 30 s, at about 21:24Z.
- 20:17:51Z: the gate opened on fi23 (free+inactive 23.3 GiB, pressure free 66 %). vcc prep rc 0 (peak footprint 24.1 GB);
  package `dad0e10eb85edd71a2ce8e7d539974651c096f1464abcf52e2254aa6ab4ee4af` (2,773,841,920 B), novelty rc 0.
  `READY idall3` at 20:18:46Z. fill1005s2 armed slot 2 at 20:19:09Z (pid 41181). The fallback exited without acting.

## Prep fallback (prep_fallback663.sh)
It acts only if `runidall3.log` says "GATE never opened (prep)":
1. It kills `fill1005s2` (command line verified) and takes over its role.
2. It waits for the same gate, with the same function and thresholds. At 10-05 01:30Z it adds one logged deviation: it preps
   when memory_pressure free is ≥ 50 % and no other build or prep is running. `vcc prep` takes about 25 s with a ~24 GB peak
   footprint (446bj2zcph: 25.06 s, 24.0 GB). 549 was prepped by hand on 09-30 the same way. If the gate is still closed at
   01:45Z, it gives up, so the prep never collides with recipe_dr524's 02:00Z build.
3. It runs run_arm_generic's post-gate steps unchanged, and also requires prep rc 0 and "prep complete".
4. It applies build_idall3.sh's READY criteria. On READY it arms slot 2 under fill1005s2's own conditions (chain1005 armed s1
   with s2 `<none>`; no live arm1005 slot-2 waiter).
