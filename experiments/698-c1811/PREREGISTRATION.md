# 698 — C1: R-ABC2H with the lab's dumphook811, on def_day, on A/B/C (tag c1811; for 10-20 s1)

Written 2026-10-18 08:3xZ by readiness [4592ee], before the lab's hook exists, before R-ABC2H's board read (10-19 s2) and before C1 is
built. Requested by the coordinator [2d6f46] (PLAN f0aaf5f9). The prediction and the HARM rule are the lab's, fixed at 08:32Z in
`experiments/811-lolo-model/PREREG-C1-ABC.md`, and are restated here unchanged.

## What is built (only if the lab reports READY by 10-19 12:00Z)
- **C1** = R-ABC2H's exact recipe (697: the frozen October recipe + `DEF_DAY_DUMP_HOOK=experiments/707-strong-tier-direction/dumphook707.sh`)
  with one value changed: the hook path becomes the lab's UNAPPLIED `experiments/811-lolo-model/dumphook811.sh`. It runs swap707.py
  unchanged, then 811's c0_l25 many-line blend on the swapped dump's bulk rows.
- Same bundle (180's def-bundle, A/B/C relabelled), same day operator (453's 66610be7), same tree.
- docs/october-decisions.json is FROZEN (04be2af3) and only read. The hook is added by `c1811_698.sh`, never by the file.
- The driver checks that the tree equals the frozen shas before and after def_day: def_day 11ff102e, builder e99bbd4a, emitter
  0f56209d, helpers 5354e449, decisions 04be2af3.

## Start condition and pre-build checks (all must hold, else NOT READY with the reason)
1. Readiness writes `work/GO-c1811` only after the lab's READY note: the hook path, its sha, and its λ-0 switch (`ka_env VAR=val`).
2. The hook file's sha equals the GO sha. It is committed and clean, which def_day also requires.
3. **KA-off (readiness's own, independent of the lab's):** with the λ-0 switch, the hook run on 697's own stack dump
   (`dump_h1x2_cd4x1_s0.25.npz`, 2dfe327d, the A/B/C θ dump production uses) must write af27efbd byte for byte. That is
   dumphook707's output on the same dump, the dump R-ABC2H was built from. The machinery was tested on 10-18 with dumphook707 as
   the stand-in: rc 0, af27efbd == af27efbd, 4 s.
4. The frozen decisions resolve to R-ABC2's recipe exactly. With the 811 hook added, the recipe must equal R-ABC2H's with only the
   hook path replaced.
5. Build gates: no other build, scripts == HEAD == frozen, AGENTS §9 memory, disk ≥ 50 GiB. Start by 10-19 19:00Z (cap), so the
   ~4 h build ends before 10-20 00:00Z.

## Gate (gate698, before READY; test698 4/4 on 697's real run)
- **The content-aware diff vs R-ABC2H** (`def_metadiff`, rtol 1e-12) is within {bulk_load_sha256, bulk_ctx_scale}, and
  bulk_load_sha256 must be a REAL difference. A hook that leaves the dump unchanged fails. The only changed export is
  DEF_DAY_DUMP_HOOK.
- 690g's lines are in def_day's log: the hook named with its path, and the validated replacement dump.
- Plus gate697's rows: postbuild all PASS, converged with every rho ≤ 1 (CLIFF OK), zero-gene audit 0, every ON switch visible in
  the meta with its resolved value, the same-line record, novelty rc 0, and a package sha matching the chain log.
- A per-context bulk (bulk_ctx.npz → bulk_load_ctx_sha256) is NOT in the allowed set. If the lab's hook writes one, the gate fails,
  and readiness asks the coordinator before any change.

## Prediction and bins (the lab's, from PREREG-C1-ABC.md, unchanged)
- pred (point) = R-ABC2H's board overall + 0.00089 (811's ABCp C1 − B1 +0.00153 ÷ ρ 1.7232). ρ range [1.44, 1.72] →
  [+0.00089, +0.00106].
- On D = C1 − pred: |D| ≤ 0.0005 AS PREDICTED; D > +0.0005 ABOVE; D ≤ −0.0005 **HARM**.
- Also HARM if C1 < R-ABC2H − 0.0005. Either HARM closes 811's 6a route for the final (LAB-DEF-RUNBOOK Amendment 3). A non-HARM
  read does not qualify 811 by itself.
- Board falsifier: none on exactness. The hook moves the bulk by design, so pds / mse are free. The DE members re-draw pass 1
  (re-roll sd ≈ 0.00025 overall, inside the bins).

## Class and tries
Transfer: the bulk direction of the strong tier (many-line consensus), on top of 707's swap. Its held-out evidence is the lab's
(composite811 on ABCp / LABC). This board read on A/B/C is the Arc-drawn check of the 6a rule. It never enters the frozen
October file. Tries: C1 is one board read; 811's C-family engine reads are the lab's (counted in its wear ledger).
