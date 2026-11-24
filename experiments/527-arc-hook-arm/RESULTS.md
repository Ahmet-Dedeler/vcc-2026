# 527 — 446bj2zmz built (the Arc-sourced uncalled-gene hook on 446bj2z): READY, not yet scored

Built 2026-09-29 10:45–12:40Z by `build527.sh` (preregistration `PREREGISTRATION.md`, 10:44Z, before any 09-30 entry).

| step | result |
|---|---|
| patch_builder517.py on f5d5b92f | builder 39c1dc6682b46644 (committed 67df0421 by commit_after_chain.sh) |
| builder suite | 56 passed, 0 failed, 1 skipped |
| 517's test_patch517.py (copies, ROOT pinned) | ALL 26 PASS |
| **KA-1** (446bj2z's env, knob unset) | output `ec3a5ae6…` == 446bj2z byte for byte |
| **KA-496** (496a's env, `clean4:0.05`) | output `e21ff752…` == 496a byte for byte |
| arm config diff vs 446bj2z | PASS: exactly `de_mde_sign`, `de_sign_emitter_sha256` |
| hook audit per context | A: 6,386 of 9,929 gated genes carry a sign (2,646 down), 86 weak targets on CLEAN4 (2,550 genes); B: 6,271 signed (2,533 down); C in the build log |
| zero-gene audit | 0 in 0 blocks; lfc < −5: 18, lfc > +3: 1 (identical to 446bj2z and 496a) |
| cliff | CLIFF OK (bulk untouched) |
| novelty | rc 0 |
| package | `submissions/446bj2zmz-…-zs250-hkmz05le1.vcc` sha256 `ed90e3b34eb75e3315c7a45ce6f73f27c96fb3b8e44021064c9749a3116017cf` |

Slot: 10-01 slot 2 if 458d reads null (not harm) or is missing (464/rule1001b.py, amended 10:56Z); otherwise the 10-07 rule.
Read with the bins in `PREREGISTRATION.md` (primary Δnmae raw vs 446bj2z; pds/mse must be bit-identical).
