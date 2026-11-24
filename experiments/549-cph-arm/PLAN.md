# 549 — arm 446bj2zcph (call-propensity block + CLEAN4 uncalled-gene hook), plan (2026-09-30T15:12Z)

- **What:** 446bj2z + `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ec…` (446bj2zcp's
  block) + `VCC_DE_MDE_SIGN=clean4:0.05` (496a's hook). Expected config diff vs 446bj2z: `rc_callset_block`, `de_mde_sign`
  (+ `de_sign_emitter_sha256`).
- **Why:** 537 v3 §0 B / §3.2 (real production emitter, H1, 3 seeds): block + hook +0.0067 / +0.0059 / +0.0059 / +0.0059 vs
  446bj2z (f 1 / 0.5 / 0.35 / cmp); block alone +0.0051 / +0.0048 / +0.0045 / +0.0053; the hook signs the block's
  swapped-in genes (0.508 sign accuracy unsteered) and cuts its leaked calls 20.4 → 8.3 per target.
- **Class (AGENTS §11.1c):** scorer mechanics (the block is computed from each context's own controls) plus a transfer-class
  prior (CLEAN4 is cross-line; H1 is held out from it).
- **Launch:** by hand, ONLY if 446bj2zcp reads ≥ PARTIAL on 10-01 (`464/read1001.py cp …`), or NULL/NEG with Δjac raw
  ≥ +0.002 and Δfid < 0 (537 §0 C). `nohup bash experiments/549-cph-arm/build549.sh > experiments/549-cph-arm/work/build549.nohup 2>&1 &`
  It pauses 537's driver between jobs (SIGSTOP) and always resumes it (SIGCONT, EXIT trap).
- **Slot:** 10-07 at the earliest, by the 10-07 rule (written before any 10-06 entry). Bins: 537 §5 (to be copied into a
  PREREGISTRATION.md here before any 10-07 entry exists).
