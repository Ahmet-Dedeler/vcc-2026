FROM THE LAB LEAD (2026-10-16 13:22Z by the clock): follow-up 7, preparation for 10-22. Budget ≤ $1, EC2 only. Deadline 10-20 18:00Z.

**Why.** The A-only panels failed the preregistered bar (ABCp 0.26, 9/21; LAB-LOG 13:21Z). So key and direction candidates stay
decided on Arc-drawn panels, and after the release that means **D/E/F's own 300 targets on the public held lines**. The lab's
first 6a candidate (811 c0_l25) is read once there. That build must run within hours of the release window closing, without
improvisation.

**Task: a one-command "panel kit" plus a known-answer dry run.**
1. `cloud801.sh panel <TAG> <pert_counts.csv> <gene_names.csv>` (or a small driver) builds, for the three cancer held lines:
   CD4 banks, the day-recipe held-out operator fits, dump801 folds `<TAG>/<held>`, the registry entries and ka804's L arm set.
   It writes standard arrays to `work/ka-<tag>/` and mirrors them to S3. Same recipe as LABC / A-only.
   - The operator fit is fresh per panel: on D/E/F the ABCp ops cannot be reused.
2. **KA dry run:** run the kit with TAG = LABCk on `data/challenge`'s files. Every array must equal LABC's (dump708-LABC-*.npz
   shas `eba095c1…` / `e97f1fbd…` / `34ecd80e…`) to |d| < 1e-12, except where a fresh operator fit differs in the last bits.
   Record the max |d| and any such difference with its cause.
3. Write the 810 (sgu1, 707) and 811 (c0_l25, c1_l25, src_l25) job files the kit's output feeds. Use their existing A-only job
   patterns with the tag as a parameter, for the lead to launch through 810 / 811 on the day.
4. **Do not touch `data/challenge-def/`** (it does not exist yet; production owns it). The kit reads paths only.
Report in RESULTS.md §17: the command, the dry-run KA, runtime and $ per panel, the job files. Then stop.
