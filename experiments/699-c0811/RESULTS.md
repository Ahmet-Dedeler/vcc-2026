# 699 — C0 (R-ABC2 + the lab's dumphook811c0) on def_day, A/B/C: READY c0811 52410f55 (readiness [1efb7d], 2026-10-20)

Preregistration: `PREREGISTRATION.md` (788a66e2; the lab's prediction quoted in 290e3d88: pred 0.2763696, HARM ≤ 0.2758696 or
< 0.2734849). Board read (10-20 s2): **ABOVE** (below).

## Built (10-20 02:50 → 07:04Z)
- **GO-c0811** 02:49:41Z (6d324b66): hook `experiments/811-lolo-model/dumphook811c0.sh` 297c3832 (hook811.py c5a5614a; 811's c0_l25,
  λ 0.25, NO swap707), committed at e603fe73 and clean; λ-0 switch HOOK811_LAMBDA=0; KA-on target 96204332.
- **Readiness's KAs** on 697's stack dump 2dfe327d (content == R-ABC2's own e8401242):
  - KA-off: λ 0 returns the input byte for byte (2dfe327d; the npz_same content check also passes);
  - KA-on: the lab's C0 dump 96204332 byte for byte;
  - 13 s.
- **Recipe:** the frozen decisions (04be2af3) resolve to R-ABC2's recipe exactly; + DEF_DAY_DUMP_HOOK=dumphook811c0.sh is the only
  addition. Tree == frozen before and after; builder restored e99bbd4a.
- **def_day:** 690g ran the hook at 02:51:51Z and validated its dump (file 43f35a77). Content checks: C0's stack dump == 697's, and the
  dump def_day used == the KA-on output (96204332).
  - lvl0 rho 0.8762 / 0.8141 / 0.7805 (R-ABC2 0.8673 / 0.7994 / 0.7766; C1 0.8870 / 0.8289 / 0.7885). lvl1 at 1.0572 / 1.0968 /
    1.1201 → 0.9730 / 0.9689 / 0.9890: E just outside ±1 %.
  - **CONVERGED at lvl2:** D 1.0572 / E 1.1030 / F 1.1201, final rho 0.9730 / 0.9794 / 0.9890.
  - Postbuild 42/42, zero-gene 0, builder peak 15.7 GiB.
- **Package** `submissions/699c0811-october-recipe-abc.vcc`, sha256 `52410f557cda181feab0bfe2f0fbdae919953c4a00b933fd2740db370eb93e5a`
  (output ab0c2d79); novelty rc 0.
- **gate699, ALL GATES PASS on the first run** (no false FAIL; def_metadiff's atol floor in force):
  - the diff vs R-ABC2 is exactly `bulk_ctx_scale` + `bulk_load_sha256` (REAL: array bulk differs); sameline is provenance-only;
  - adopted == [DEF_DAY_DUMP_HOOK]; 690g hook lines OK.
  - **READY c0811 52410f55** at 07:03:56Z, in c0811_699.log (the coordinator's 1020s2 ready-log) and READY-science.log; evidence
    build_c0811/runc0811.log.
- The live run's KA-on and dump-check lines say "C1 dump": this is the stale wording of the frozen copy (fixed in the source, f6450e98).
  It means C0's dump 96204332.
- The queue continues: W1 (673) started def_day at 07:04:20Z, right after C0.

## Board read (10-20 s2, entry 4rq3qHIuv8bHAxpTkUeu, captured 07:29:25Z by the coordinator's waiter; recorded by readiness)
- **overall 0.2780624 vs pred 0.2763696: D = +0.0016928 → ABOVE** (> 0.2768696). **Not HARM on either line** (> 0.2758696, > 0.2734849):
  the swap-OFF X2 route (SIXA Amendment 1) is not closed by its A/B/C condition. It is the team's highest board overall so far.
- Vs R-ABC2 (nnHQfcZHsT3as9EafA3i, 0.2739849): Δoverall **+0.0040775** (pred +0.0023847 = 0.00269 ÷ ρ 1.128, so 1.7× the deflated
  increment and 1.5× the undeflated).

  | member | scaled Δ vs R-ABC2 | C1 − R-ABC2H (for comparison) |
  |---|---|---|
  | pds | −0.000270 (raw −0.000074) | −0.000588 |
  | mse | **+0.023141** (raw −0.021983) | +0.017726 |
  | fid | +0.000010 | −0.000008 |
  | reach | +0.000603 | −0.000023 |
  | jac | −0.000013 | +0.000003 |
  | nmae | +0.000993 | +0.000119 |

- **Bulk (Δpds + Δmse)/6 = +0.0038118**, almost all mse. **DE +0.0002657**, about 1 re-roll sd.
- **C0 vs C1 (same 811 blend, without vs with 707's swap):** C0 0.2780624 vs C1 0.2777649, so **C0 − C1 = +0.0002975**. In bulk,
  where these members are exact: C0 − R-ABC2 = +0.0038118, and C1 − R-ABC2 = R-ABC2H's +0.0006458 + C1's +0.0028564 = +0.0035022.
  So the no-swap arm is +0.00031 higher in bulk. **On A/B/C, once 811 is applied, the swap adds nothing and is slightly negative**: 811's
  many-line consensus already carries the swap's strong-tier direction gain. Facts only; which branch ships is SIXA's rule
  (708's day-of swap decision plus the lab's DEFo/DEF reads).
- **Class:** transfer (811 without the swap). Tries on the board for 811: 2 (C1, C0).
