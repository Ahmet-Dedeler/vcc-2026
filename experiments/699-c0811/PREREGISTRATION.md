# 699 — C0: R-ABC2 + the lab's dumphook811c0 (811 on the UNSWAPPED dump), on def_day, on A/B/C (tag c0811)

Written 2026-10-20 02:3xZ by readiness [1efb7d], before the lab's C0 hook exists and before any C0 build or read. Requested by the
coordinator [3c7a7e] under SIXA-DECISION Amendment 1 (02:35Z): if 708 turns the swap OFF on D/E/F, X2 = R0 + 811 = C0, and that
increment had no A/B/C board read. **The prediction is the lab's**, written before the build (to be quoted here and in RESULTS).
HARM (the coordinator's, fixed now): D ≤ −0.0005 vs the lab's pred, OR C0 < R-ABC2 − 0.0005 = **0.2734849**. Either closes the
swap-OFF X2 route.

## What is built (only if the lab reports READY; no C0 def_day build starts after 10-21 12:00Z)
- **C0** = R-ABC2's resolved recipe (695: the frozen October recipe 04be2af3, which carries no hook) + ONE export:
  `DEF_DAY_DUMP_HOOK=experiments/811-lolo-model/dumphook811c0.sh` (UNAPPLIED, the lab's: hook811 λ 0.25 on the UNSWAPPED dump).
- Same bundle (180's def-bundle), same day operator (453's 66610be7), same tree. The driver checks the frozen shas before and after
  def_day: decisions 04be2af3, def_day 11ff102e, builder e99bbd4a, emitter 0f56209d, helpers 5354e449.

## Start condition and pre-build checks (else NOT READY with the reason)
1. Readiness writes `work/GO-c0811` only after the lab's READY note: hook path + sha, the λ-0 switch, and the KA-on target.
2. The hook's sha equals the GO sha; it is committed and clean; it is neither dumphook707 (R-ABC2H) nor dumphook811 (C1).
3. **KA-off (readiness's own):** with the λ-0 switch, the hook run on 697's stack dump (2dfe327d, content-identical to R-ABC2's
   own e8401242) must return that dump's CONTENT unchanged (def_metadiff.npz_same). Tested 10-20: a hook that changes the bulk
   fails it ("array bulk differs"); a content-identical copy passes.
4. **KA-on:** the hook as def_day runs it, on the same dump, writes the lab's C0 dump byte for byte (GO `ka_on`).
5. The frozen decisions resolve to R-ABC2's recipe exactly; the hook is the only addition.
6. Build gates: no other build, scripts == HEAD == frozen, AGENTS §9 memory, disk ≥ 50 GiB.

## Gate (gate699, before READY; test699 4/4 on 695's real run)
- **The content-aware diff vs R-ABC2** (`def_metadiff`, rtol 1e-12 + atol 1e-12) is within {bulk_load_sha256, bulk_ctx_scale},
  and bulk_load_sha256 must be REAL (a no-op hook fails). The only changed export is DEF_DAY_DUMP_HOOK.
- 690g's hook lines; all of gate697/698's rows (postbuild, CLIFF OK, zero-gene 0, switches visible, same-line, novelty,
  package sha). After def_day, C0's stack dump == 697's in content, and the dump def_day used == the KA-on output in content.

## Class and tries
Transfer (811's many-line consensus bulk direction, without 707's swap). Its held-out evidence is the lab's (DEFo/DEF, C0 − B0).
This A/B/C board read is the Arc-drawn check of the swap-OFF branch. It never enters the frozen October file. Tries: 1 board read.

## The lab's prediction (quoted 10-20 02:5xZ from `experiments/811-lolo-model/PREREG-C0-ABC.md`, committed e603fe73 at 02:38:48Z, before any C0 number)
- pred (point) = 0.2739849 + 0.00269 ÷ ρ 1.1280 = **0.2763696** (range over the two class points 0.2755459 … 0.2790331).
- AS PREDICTED [0.2758696, 0.2768696]; ABOVE above it.
- **HARM: C0 ≤ 0.2758696 OR C0 < 0.2734849.** HARM closes X2 on the swap-OFF branch only.
- The lab's expectation, recorded before the read: given C1 (the folds under-read 811 by ~1.9× on top of the swap), C0 may read ABOVE.
- GO-c0811 written 02:49:41Z (6d324b66): hook 297c3832, λ-0 HOOK811_LAMBDA=0, KA-on 96204332.
