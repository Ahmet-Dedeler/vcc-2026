# 502 — October pass-throughs for the DE riders (prepared 2026-09-27 07:35Z; not applied)

**Why.** `def_day.sh` wipes every `VCC_*` at start and exports only its own knobs, so a DE rider that transfers on A/B/C
does not reach October unless the day script passes it. Two are in flight: zs (`VCC_DE_MAG=dn_k250_zs`, 09-28 read as
446bj2z) and the uncalled-gene sign hook (`VCC_DE_MDE_SIGN=clean4:0.05`, 10-05 read as 496a). On H1 they stack exactly
(496z).

**Patch** (`patch_defday502.py`, on 497's def_day e29e793a + helpers 60b2a777):
- `DEF_DAY_DE_MAG=<variant>` and `DEF_DAY_MDE_SIGN=clean4:<tau>`: validated before setup. Each is exported for the level
  builds only, after the dump, and written to `def_day.final.env`.
- A new helper `riders_match` checks the final meta before the `_final` link: `de_mag.variant` and `de_mde_sign.{source,tau}`
  must equal the request, and must be absent when not requested.
- 5b's allowed diff gains `de_mag` / `de_mde_sign` only when set. 499 F2 already allows `de_sign_emitter_sha256` once
  496p's emitter pin differs from the incumbent's.
- Both unset: the same exports, argv and checks as 497's script. `riders_match` then requires both keys absent, which
  every current build satisfies.

**Checks at prep time:** `bash -n`; `test_riders502.py` 9/9 (absent/present/mismatch cases); 497's helper tests 9/9 on
the 502 helpers; 485's `repro_fixes.py` PASS on the 502 copies (`repro/`).

**To apply** (after a rider reads "transfers"): patch the committed 497 files, run an inert known answer (efp relaxed
== ac49d319 with both unset), then a rehearsal with the rider on (a 497dr-like build scorable against 497dr).
