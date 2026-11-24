FROM THE LAB LEAD (2026-10-16 10:40Z by the clock): follow-up 4. This one matters for every fresh-fold read in the lab.
Budget ≤ $6, EC2 only (the local-job HOLD is still on). Deadline: done by 10-17 12:00Z.

**Finding that prompts it.** `experiments/800-checker/panelrank800.py`, results in `experiments/800-checker/work/panelrank800.txt`:
for each held line, it ranks ABCp/<held> among L1–L6/<held> on the same arms.
- A/B/C's panel is NOT exchangeable with your six arc800 draws. ABCp is the bottom of its line's seven folds in 14/21 (arm, held)
  cells, where 3.0 are expected; its mean percentile is 0.19.
- By class:
  - key arms (495 cap6−cap3, 705): bottom in 6/6, z −1.8 to −5.5;
  - direction arms (707, sgu1, 811 c0_l25): bottom in 8/9;
  - shape arms (479p2g, 479sp2): inside the spread (ranks 2–6).
- ABC062 shows the same pattern.

**The confound to remove.** ABCp folds use A/B/C's fold rules (CD4 ×1.5 all-rows sizing, H1 ×2, A/B/C's controls and CD4
banks, the cap6+P2g(abc) incumbent). The L folds use D/E/F rules. So "the panel's targets" and "A/B/C's fold rules" are
confounded.

**Task.**
1. Build **LABC/<held>** for k562_gwps, orion_hct116 and orion_hek293t: exactly the L-fold recipe (dump801, D/E/F rules, a
   fresh day-recipe 3-donor operator fit for this panel, as L1–L6 got), with ONE change: the panel is A/B/C's 300 targets
   (`data/challenge/pert_counts.csv`). Same held truths.
   - KAs as for your L folds.
   - Record operator shas and the dump sha in the registry, status `labc` (UNAPPLIED to nothing else; it is a new fold set).
2. Run ka804's arm set on LABC (`experiments/804-known-answers/ka804.py`, unmodified, `KA804_OUT` outside the repo). The L
   arm set applies: inc = pct+P2g.
   - Collect standard-format arrays to `experiments/801-lab-folds/work/ka-labc/{bulk,meta}/` and mirror them to S3, as for
     ka-abcp.
3. **Truth-strength profile, per fold.** For every fold of the three held lines (L1–L6, ABC062, ABCp, LABC), from the
   standard arrays (`d` = per-target truth energy in the frame, on the cap3 arm), report:
   - median d, mean d, the top-30 share of Σd, the count of targets with d > 4× the fold median;
   - the Spearman correlation across folds between each profile statistic and the 495 (cap6−cap3) per-fold gain.
   Compute it on EC2, or wait for the hold to lapse (window800 prints CLEAR).

**Reading (written now, before any number):**
- LABC inside the L spread on key/direction arms → A/B/C's fold rules caused the bottom ranks.
- LABC at the bottom too → the panel's targets did. Then the truth-strength profile should name the property.

Report in RESULTS.md §14 with runs and $, then stop. I run panelrank800 with LABC added.
