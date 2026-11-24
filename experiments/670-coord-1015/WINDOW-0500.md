# Coordinator checklist for the window after 702 is READY (~05:00Z 10-15)

Trigger: `READY pie702w25… <sha>` in experiments/701-kup-callset/work/READY-science.log AND no bash running
experiments/702-pie-bulk-arm/build702.sh (the chooser then arms 10-15 s2 by itself; check experiments/670-coord-1015/work/choose1015s2.log).

1. Sweep 702's r1 / r2 prediction h5ads once their .vcc sha is logged (log in experiments/670-coord-1015/work/SWEEP.log). Disk ≥ 60 GiB.
2. 693p (dn_mag v2): requires [17adab]'s six H1 KA records + abc693 (v1 == DN package block, v2 == v1). Then
   `uv run python experiments/693-dnmag-speed/patch693p.py --root . --check` → apply → run_builder_tests → commit the builder +
   src/vcc2026/dn_mag2.py + def_day.sh (690a meta check moves to v2) with a pathspec; KA hashes follow in step 3's commit.
3. `nohup bash experiments/693-dnmag-speed/ka693b.sh > experiments/693-dnmag-speed/work/ka693b.nohup 2>&1 &`
   → KAU2 PASS (== ec3a5ae6) and DNV2 PASS (novelty says ALREADY BEEN SCORED for wzaVLFTmlAMW3Mfw0dD8) in work/ka693b.log, then END.
   If either FAILs: revert 693p (v1 stays on disk), tell [17adab].
4. Release [17adab]'s HOLD (rm experiments/690-def-readiness-1015/work/HOLD) after ka693b's END (or immediately if 693p waits).
5. K-UP (ku30, build701.sh) after the rehearsals, as a 10-16 candidate; tell [bf003b].
