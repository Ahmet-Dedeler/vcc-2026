# 697 — R-ABC2H: R-ABC2 + 707's strong-tier hook, on def_day, on A/B/C (readiness [4592ee])

Preregistration: `PREREGISTRATION.md` (9056de70, before R-ABC2's read). Board read (10-19 s2): **COMPOSES** (below).

## Board read (10-19 s2, entry ndt8swt7ILYn3qMJYZ9a; read695 --rabc2h by the coordinator [223d39]; recorded here by readiness)
- **overall 0.2748936 vs pred 0.2750179: D = −0.0001244**, inside COMPOSES (0.2745179, 0.2755179). The swap composes with the
  October recipe on A/B/C, so the 6a route (R0 + 707's swap, switched on by 708's day-of rule on D/E/F) has its A/B/C
  composition read. Nothing changes for R0: the hook stays off in the frozen decisions file.
- Vs R-ABC2 (nnHQfcZHsT3as9EafA3i, 0.2739849): Δoverall **+0.0009086** (707 on ALLSW: +0.0010330).

  | member | scaled Δ | |
  |---|---|---|
  | pds | −0.002351 | raw −0.001041 |
  | mse | +0.006226 | raw −0.005922 (lower is better) |
  | fid | +0.000144 | |
  | reach | +0.000076 | |
  | jac | +0.000016 | |
  | nmae | +0.001341 | |

- **Bulk, in overall units:** (Δpds + Δmse)/6 = **+0.000646**, vs 707's primary +0.000842 on ALLSW. These two members reproduce
  exactly, so θ's base takes about 0.0002 off the swap's bulk gain, a small real interaction. It stays above 707's TRANSFERS line
  for the bulk. The swap trades pds for mse, as on ALLSW.
- **DE members:** (Δfid + Δreach + Δjac + Δnmae)/6 = +0.000263, about 1 re-roll sd (≈ 0.00025): the hook moves the bulk, so pass 1
  re-draws.
- **Class:** transfer (707's strong-tier direction). Its held-out reads are 707's and the lab's. This board read says only that
  it composes with the frozen recipe on A/B/C. Whether it ships on D/E/F is 708 AMENDMENT 2's day-of rule (wait708 →
  decision708-d708def.txt); the coordinator decides the 6a build from that file.
- C1's prediction (698) resolves from this read. The lab filled it in at 02:18Z: pred = 0.2748936 + 0.000888 = **0.2757814**
  (`811/PREREG-C1-ABC.md`; 698 RESULTS).

## Built (10-17 11:22 → 15:33Z) and re-gated (16:21Z)
- R-ABC2's resolved recipe (== 695's recipe695.json: R-ABC + θ 0.25) + `DEF_DAY_DUMP_HOOK=experiments/707-strong-tier-direction/dumphook707.sh`,
  added by the driver. docs/october-decisions.json keeps dump_hook off for R0.
- **The hook's first def_day run:** 690g ran it on the θ stack dump (`dump_h1x2_cd4x1_s0.25.npz`) and def_dumphook_check validated it;
  the replacement dump af27efbd became VCC_BULK_LOAD.
- Re-level: lvl0 0.879 / 0.815 / 0.783 (R-ABC2 0.867 / 0.799 / 0.777); lvl1 0.973 / 0.969 / 0.992, all safe. CONVERGED at lvl2,
  D 1.0555 / E 1.1019 / F 1.1121, final rho 0.973 / 0.979 / 0.979. Postbuild 42/42, zero-gene 0, novelty rc 0.
- Package `submissions/697rabc2h-october-recipe-abc.vcc`, sha256 `a43a64d072014f2a34723eb965c6f3c34d4686164a86b68f3e401178189a7dd8`.

## The gate's false FAIL and the fix
gate697 (15:33Z) failed the content-aware diff vs R-ABC2: it listed bulk_sharpen_pt next to the hook's keys. The ONLY difference
inside that record was the two tercile `edges` at the last float digit (0.20105594997123832 vs …835; 0.6530179660732803 vs …802).
Counts, carried flags, the key and its sha, moved_vs_tercile: all identical. So the hook's rewrite of the dump perturbs the round-off
of the key-ratio edges and nothing else. `scripts/def_metadiff.py` now compares float config values to 1e-12 relative (default;
`--float-rtol 0` = exact). Test M9: a 1-ULP change is not a difference, 0.1 % still is, rtol 0 is exact; 9/9. The other gates' tests
still pass (694 11/11, 695 6/6, 697 3/3). `ready697r2.sh` re-gated the finished build (package sha unchanged, novelty re-run):
ALL PASS, diff vs R-ABC2 = exactly bulk_ctx_scale + bulk_load_sha256 (the hook's keys), `READY rabc2h a43a64d0…` 16:21:56Z.

## Prediction (frozen in the prereg)
pred = R-ABC2's board read (10-18 s1) + 707's Δoverall +0.0010330; bins ±0.0005. BELOW ⇒ 6a ships without the swap. R0 is unaffected.
**Resolved 10-18 02:2xZ:** R-ABC2 read 0.2739849 (nnHQfcZHsT3as9EafA3i) ⇒ **pred = 0.2750179**; COMPOSES (0.2745179, 0.2755179);
BELOW ≤ 0.2745179; ABOVE ≥ 0.2755179. Bulk check: 707's primary Δ(pds+mse) on ALLSW was +0.0008417 (overall units).
