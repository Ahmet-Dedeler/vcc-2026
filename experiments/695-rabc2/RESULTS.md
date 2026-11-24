# 695 — R-ABC2: the October recipe after the 705 / 707 reads, on def_day, on A/B/C (readiness [4592ee])

Preregistration: `PREREGISTRATION.md` (d688fa90, Amendment 1 f82f931c, resolved inputs d349c231 frozen before GO). Board read
(10-18 s1): **COMPOSES**, θ stays for R0 (below).

## Board read (10-18 s1, entry nnHQfcZHsT3as9EafA3i, captured 01:40:13Z by the coordinator's waiter; recorded 02:2xZ by readiness)
- **overall 0.2739849 vs pred 0.2740445: D = −0.0000596**, inside COMPOSES (−0.0005, +0.0005). θ 0.25 stays decided in
  docs/october-decisions.json; nothing reverts.
- Vs R-ABC (PooEbFClgTLMQoXiOCQP, 0.2743318): Δoverall −0.0003469 (705 on ALLSW: −0.0002874).

  | member | scaled Δ | |
  |---|---|---|
  | pds | −0.001575 | raw −0.000669 |
  | mse | −0.000356 | raw +0.000331 |
  | fid | −0.000038 | |
  | reach | −0.000270 | |
  | jac | +0.000020 | |
  | nmae | +0.000137 | |

- **Bulk, in overall units:** (Δpds + Δmse)/6 = −0.000322. 705's primary on ALLSW was −0.000240. These two members reproduce
  exactly, so the extra −0.000082 is a real (small) interaction of θ with R-ABC's base. It is still far from 705's HARM line (−0.0005).
- **DE members:** (Δfid + Δreach + Δjac + Δnmae)/6 = −0.000025, inside the re-roll sd (≈ 0.00025). θ moves ballast K, so pass 1
  re-draws.
- **Class:** transfer (θ's held-out read is 705/522). The A/B/C board read says θ composes with the October recipe as built by def_day;
  it costs ~0.0003 on A/B/C by design.


## Built (10-17 01:32 → 05:37Z)
- Recipe = R-ABC's env + `DEF_DAY_LODO_STRENGTH=0.25` (705 CLEARED; 707 WEAK, so no hook; K stays 2,600). def_day.sh on 180's def-bundle
  with 453's operator 66610be7, the tree with def_day.sh 11ff102e (690i + AGENTS floor; 690k R-KA proof) and builder e99bbd4a.
- Package `submissions/695rabc2-october-recipe-abc.vcc`, sha256 `54e6c5cf1ca066e1506c7abaeec2ddb399d6be0f684ab0a35ec2607d4c165ce0`;
  novelty rc 0.
- Re-level: lvl0 rho 0.867 / 0.799 / 0.777 (R-ABC 0.872 / 0.802 / 0.773). lvl1 at 1.0626 / 1.1069 / 1.1230 overshot F (1.0038) and was
  deleted as past the cliff. CONVERGED at lvl2, D 1.0626 / E 1.1132 / F 1.1105, final rho 0.973 / 0.9795 / 0.9774. Postbuild 42/42, zero-gene 0.

## θ's set-value rehearsal on A/B/C (the first def_day run with the rider)
def_day logged "LODO strength rider: theta 0.25"; the sidecar records θ 0.25; the dump is `dump_h1x2_cd4x1_s0.25.npz`; tercilec carried the
LODO key (lodo_lo0.5_pct_def, 294 carried). The content-aware diff vs R-ABC is exactly bulk_ctx_scale, bulk_load, bulk_load_sha256 and
bulk_sharpen_pt. That is the signature predicted on 10-16: θ rewrites only the dump's key, and pert_scale does not move. gate695 ALL PASS.

## Prediction and bins (frozen before the build)
pred = 0.2743318 − 0.0002874 = **0.2740445**. COMPOSES (−0.0005, +0.0005); ABOVE ≥ +0.0005; BELOW ≤ −0.0005 (≤ 0.2735445) reverts θ for R0.
