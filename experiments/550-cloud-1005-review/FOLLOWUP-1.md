Follow-up 1 (2026-09-30 21:40Z, from the session that launched you; worker 537's notes). Add these to your review:

G. `experiments/550-cloud-1005-review/snapshots/537-PREREGISTRATION-draft-CPH-20260930T1635Z.md` is 537's draft preregistration
   for the 549 arm (446bj2zcph). Check 549's build path against it too. Its falsifiers: config diff vs 446bj2z exactly
   `rc_callset_block`, `de_mde_sign`, `de_sign_emitter_sha256`; vs 446bj2zcp (`snapshots/446bj2zcp.h5ad.meta.json`) exactly
   `de_mde_sign`; pds and mse equal to 446bj2z's to <= 1e-6 raw. Say which of these a local gate enforces before a slot and
   which only the board read would catch (see item E: write the bulk-identity check if none exists).
H. 537's RESULTS §5.4 (snapshot) recommends that in October `DEF_DAY_MDE_SIGN` be conditional on `DEF_DAY_CALLSET_BLOCK=prop`,
   with CLEAN4 (`clean4:0.05`), not Arc's programme, as the default hook source (its P9/P12: CLEAN4's zero-information floor on
   the block is ~0, Arc's is -0.002 ... -0.003). On the patched tree (530 + 538's re-anchored 528 + 538 extras): can
   `docs/october-decisions.json` + `scripts/october_recipe.py` express "hook on only with the block, source clean4", and
   does anything refuse or warn on the risky combinations (hook without the block; the Arc source without an explicit
   decision)? If not, write the smallest patch script (not applied) and a test.
