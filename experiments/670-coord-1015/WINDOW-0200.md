# Coordinator checklist for the window after dn671's END (~02:00Z 10-15)

Trigger: `experiments/671-dnmag-arm/work/dn671.log` ends with `END` and pid `$(cat experiments/671-dnmag-arm/work/dn671.pid)` has exited.
HOLD (`experiments/690-def-readiness-1015/work/HOLD`) stays in place throughout; [bf003b]'s H1 K-UP pricing and K-UP build come next.

1. dn671.log must show `READY dn671 <sha>` and `KAU PASS (== ec3a5ae6 on builder 8f8f6d5c…)`, and the live builder == 8f8f6d5c
   (pins restored). If KAU FAIL: do not commit 631; restore experiments/671-dnmag-arm/work/builder.pre631.py and remove
   src/vcc2026/dn_mag.py, tell [17adab] and [bf003b], and NOGO the DN submit (experiments/670-coord-1015/work/NOGO-1015s1) only if
   DN's own falsifiers would be suspect.
2. Commit 631: `git add -f src/vcc2026/dn_mag.py` then
   `git commit -m "631 applied: … KA1015u == ec3a5ae6 …" -- scripts/make_submission_pm_cd4.py src/vcc2026/dn_mag.py experiments/671-dnmag-arm/`.
3. Apply 690a v2 (patch690a.py sha 0af0ce21…): `--check`, apply, `scripts/test_october_recipe.py`, `experiments/644-runbook-amendments/verify644.py`
   (restore its json afterwards).
4. Apply 690c (patch690c.py sha 845aa867…): `--check`, apply; `bash -n scripts/def_day.sh`; test690c.
4b. Apply 690g (patch690g.py sha e24b09d4…, DEF_DAY_DUMP_HOOK + VCC_BULK_LOAD_CTX day-of path; needs 672 in the builder only when a hook
   emits bulk_ctx.npz): `--check`, apply; test690g; test_october_recipe stays green. Adds PENDING `dump_hook`: decide it at the 10-20
   freeze (off, or the PIE hook's path if 702 is adopted).
5. Decide (with evidence lines): `de_uprn b_k100_q20r3p10` (630arm QJRJK6OC6Zw1kjwnZqid, PARTIAL bin, "October adopts, flagged");
   `mde_sign_keygate off` (612 never read alone; idall3 BELOW). `de_dnmag d_m200` (DN wzaVLFTmlAMW3Mfw0dD8 PARTIAL, Δnmae raw −0.001646; October adopts flagged, capacity floor reported day-of; BOARD-RESULT-DN.md).
6. Commit 690a + 690c + decisions with a pathspec; message [17adab] the shas (their run690 gate needs builder == HEAD, which step 2 restores).
7. Tell [bf003b] the builder is free: 702 r1 → ntclevel417 → r2 first (their build702 is armed on dn671 exit), then rehearsals, then K-UP.
