# 456 — the zs rider pre-built on every 09-27 base (for the 09-28 slot 2), and the 09-28 automation

The 09-28 rule (NEXT-SLOT-DECISION, fixed 01:40Z 09-26, before any 09-27 entry) sends the zs pass
(`VCC_DE_MAG=dn_k250_zs`, a pure-nmae post-pass) on the best 09-27 base. These builds make every branch ready with no
build on the day. zs bins: `433-coarse-mde-arm/PREREGISTRATION-0926.md`.

| arm | = | config diff vs its base | CLIFF | package | novelty |
|---|---|---|---|---|---|
| 446bz (built 09-24) | 446b + zs | exactly `de_mag` | OK 0.9603 / 1.0282 / 1.0785 | `8020a8d5…` | rc 0 |
| **446bj2z** (03:42Z) | 446bj2 + zs | exactly `de_mag` | OK 0.9602 / 1.0281 / 1.0784 | `f252b815…` | rc 0 |
| **446br2z** (04:21Z) | 446br2 + zs | exactly `de_mag` | OK 0.9601 / 1.0280 / 1.0783 | `d5bc76f0…` | rc 0 |

All three zs variants are built and gated, so 09-28's slot 2 is ready in every branch. 446br2z's zero-gene audit is 0
(7 cells below −5).

- 446bj2z's zero-gene audit is 0 (18 cells below −5, 1 above +3).
- Incident, 02:32Z: the first launch waited for `" END$"` in chain455.log. That log quotes 453's `…END` line, so chain456
  started beside 455's build. run_arm's other-builders gate parked it; it was killed before pinning anything and
  re-armed to wait for chain455's own last line plus its exit. The rule is now in AGENTS 6.2.
- `chain0928.sh` (pid 53630) applies the rule after the 09-27 second slot and arms `arm0928.sh` s1/s2. `read0928.py`
  reads the entries.
