# 537 — DRAFT board bins for 446bj2zmzcp: the call-propensity block with the ARC-sourced uncalled-gene hook, read against 446bj2zcph (not armed; for the coordinator)

Written 2026-09-30 by the 537 worker from real production-emitter H1 runs (RESULTS §3.2, P4), before any 10-01 entry exists
and before the arm is built. Rule draft v2 (19:25Z, after P12 below) reads it on 10-08 (after CPHK20 when that applies) when
the block reads ≥ PARTIAL (or is rescued) on 10-01 AND 496a reads ≥ PARTIAL on 10-05, against 446bj2zcph (10-07 s1): one
question — the hook's SOURCE at a fixed block and a fixed delivery. (v1 had it on 10-07 s2.)

## The arm
- **446bj2zmzcp** = 446bj2z + `VCC_RC_CALLSET_BLOCK=prop:experiments/528-callprop-block/work/callblock528_abc.npz:04fb24ecf81582f0af49b4e44b4dfb7b35def66425b20afd6e834c0d907bf0dd`
  + `VCC_DE_MDE_SIGN=file:experiments/517-flex-programme/work/arch1meanz_ne_hook_abc.npz:165c6b78a16de2990ed87ca0e3436ab8c4023f7dddcdb30cf00b8f5eebda9f6b:0.05:0.05:le1`
  (527's hook: Arc's H1 programme as the mean z over its non-essential rows, on the design-strong targets; CLEAN4 τ 0.05 on
  the weak third by 468's strength measure).
- Both knobs are in the committed builder 4daae237 / emitter 0f56209d (535's KA535c rebuilt 446bj2zmz on it).
- **Expected config diff vs 446bj2zcph: exactly `de_mde_sign`**; vs 446bj2z: `rc_callset_block`, `de_mde_sign`,
  `de_sign_emitter_sha256`. pds and mse equal 446bj2z's to ≤ 1e-6 raw.

## Why (H1, the production emitter, 3 seeds; the Arc source attenuated to r 0.37, predict517's reading of 458d)
| comparison | f 1 | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| block + Arc hook (r 0.37) − block + CLEAN4 hook (the arm vs 446bj2zcph) | −0.0019 | −0.0008 | +0.0007 | +0.0006 |
| its members (fid / jac / reach / nmae raw) | +0.0006 / −0.0002 / −0.0133 / −0.0013 | +0.0006 / −0.0001 / −0.0068 / −0.0008 | +0.0006 / −0.0001 / +0.0000 / −0.0014 | +0.0006 / −0.0001 / +0.0003 / −0.0010 |
| block + Arc hook (r 0.37), total vs 446bj2z | +0.0048 | +0.0051 | +0.0065 | +0.0065 |
| block + CLEAN4 hook (446bj2zcph), total vs 446bj2z | +0.0067 | +0.0059 | +0.0059 | +0.0059 |
| block on the Arc-hook base: r 0.37 (r 0, 528) | +0.0062 (+0.0085) | +0.0054 (+0.0074) | +0.0053 (+0.0075) | +0.0061 (+0.0082) |

Seeds (f 0.35, arm vs 446bj2zcph): +0.00061 / +0.00075 / +0.00068. The Arc source is better on the two robust members (fid
+0.0006, nmae −0.0008 … −0.0014 raw) and worse on reach at f ≥ 0.5 (−0.007 … −0.013 raw: a small-gate term, the lottery 458d
paid on the board).

## D/E/F: the Arc source's value is steep in its attenuation (P12, 3 seeds, added 19:25Z)
Arc − CLEAN4 hook on the block (f 1 / 0.5 / 0.35 / cmp), by how much of the Arc programme's sign information survives:

| Arc source | vs the CLEAN4 hook on the block | members at f 0.35 (fid / jac / reach / nmae raw) |
|---|---|---|
| r 0 (H1's own information; 528) | +0.0062 / +0.0064 / +0.0061 / +0.0062 | +0.0041 / −0.0004 / +0.0001 / −0.0146 |
| r 0.37 (A/B/C's, 458d) | −0.0019 / −0.0008 / +0.0007 / +0.0006 | +0.0006 / −0.0001 / +0.0000 / −0.0014 |
| no sign information (signs permuted among its 6,773 covered genes) vs CLEAN4's real hook | −0.0042 / −0.0045 / −0.0030 / −0.0028 | |
| (the same, both sources permuted) | −0.0021 / −0.0035 / −0.0018 / −0.0023 | −0.0016 / −0.0016 / −0.0023 / −0.0010 |

With no sign information the Arc hook on the block LEAKS (78 unsuppressible calls per target vs 20 for the block alone; its
coverage is 2.7× CLEAN4's) and costs −0.0016 … −0.0034 against the block alone, where CLEAN4's permuted hook costs ≈ 0.
The value crosses zero near r ≈ 0.40, just past A/B/C's measured 0.37. For October the Arc source therefore has an
expectation ≈ 0 against CLEAN4 at A/B/C-like attenuation and a −0.003 downside if D/E/F attenuate more, so the A/B/C read
must show MORE than the r 0.37 prediction before the source changes. The bins below are raised accordingly (v2, 19:25Z).

## Primary statistic, noise and predictions
- **Primary = Δ(fid + nmae) in overall units vs 446bj2zcph** (σ ≈ 0.00005 overall; RESULTS §2), reach reported.
- Board prediction: fid +0.0006 raw × 0.07 … 0.67 (sign-prior class) and nmae −0.0008 … −0.0014 raw × T (T = 496a's
  board factor) → **Δ(fid + nmae) overall +0.0001 … +0.0007**; reach −0.013 … 0 raw × (1 … −0.76), i.e. anywhere in
  −0.0025 … +0.0019 overall.

## Bins (primary Δ(fid + nmae) overall vs 446bj2zcph; v2 after P12)
- **≥ +0.0010 and Δreach raw ≥ −0.005: ARC.** The source beats its r 0.37 prediction by ~2× or more (A/B/C attenuate less than
  458d implied), which leaves D/E/F a margin: October's hook source is Arc's (`DEF_DAY_MDE_SIGN=hookbank:h1_2025_full:0.05:0.05:le1`
  + `DEF_DAY_HOOK_ESTIMATOR=meanz_ne`, 530's bundle), with a rehearsal first.
- **+0.0003 … +0.0010, or ≥ +0.0010 with Δreach raw < −0.005: A/B/C-ONLY.** The read matches r ≈ 0.37: the source earns its
  signs on A/B/C, but at that attenuation its D/E/F expectation is ≈ 0 against CLEAN4 with a −0.003 downside (P12). October
  keeps CLEAN4; the Arc source may ride on A/B/C-only entries.
- **< +0.0003: CLEAN4** (the simpler, cross-line source) stays October's hook.
- **Falsifiers:** pds / mse ≠ 446bj2zcph's by > 1e-6 raw (stop); |Δjac| > 0.001 raw (a call leak: the Arc hook leaked
  ~590 calls per target on H1's gated-programme base, not on this base — compare `call_set_unsuppressible`).
