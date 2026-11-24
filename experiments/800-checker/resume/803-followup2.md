FROM THE LAB LEAD (2026-10-16 09:20Z). Your context is intact. Thank you for the sharpening read (logged, committed). One instrument
gap now blocks every bulk candidate whose ballast K moves > 5 %: 811's many-line consensus (no dossier only because of that) and
810-type composition changes. You wrote it yourself: "the engine can't yet swap the harness's bulk". Build that hook (budget $15):

1. `de803.py` arms get an optional `bulk_transform` field: a module path + sha256 + args, called on the H1 surface's bulk rows (the
   harness's own composition for its targets, before key scaling / ballast, wherever production applies a dump) and returning new
   rows of the same shape. The emitter, carriers, priors and everything else stay as E0 / the arm says.
   KA: an identity transform reproduces E0's runs bit-for-bit on both surfaces (all stats arrays, as your item-1 KA).
2. Price, 3 seeds × both surfaces, six members (DE at f 1 and f 0.35, pds/mse from emitted cells), base = `inc`:
   (a) 811's winner `c0_l25`: blend each H1 target's row norm-preservingly toward U = the equal-weight mean of the per-line unit rows
       of the D/E/F day's bank list (811's `m811.py` DAY_GROUPS) with the WHOLE H1 group held out (VCC25 H1 and Pan H1 never
       enter), λ 0.25. Reuse 811's bank loader and blend code verbatim (cite blob + sha); λ 0.10 as a second arm.
   (b) if cheap: 810's `sgu1` (the strong-tier base → pre swap) on H1's own composition parts, if the H1 harness exposes `pre`
       and the whitened base; otherwise say what blocks it.
3. Report per arm: Δ for all six members vs inc, K moves, and whether the DE members move beyond seed noise. This is the
   "DE side effects" read that rule 6 of 811 and 810's DE-H1 FIRST need. Write it to RESULTS (§8) and stop with a short message.
Rules as before (window800, EC2 only, collect only declared outputs).
