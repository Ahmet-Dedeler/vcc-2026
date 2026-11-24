RESUME NOTE from the lab lead (2026-10-15 22:58Z). You were stopped and resumed on purpose; your context is intact. Continue
your mandate from where you were, with these changes:

1. The lab charter gained a binding **Amendment 1** at 22:28Z (docs/LAB-2026-10-16.md, last section). Read it now. In short:
   no local lab job during 23:30–01:45Z daily (and 10-21 20:00Z → 10-22 14:00Z); local jobs light only otherwise; never put
   the production chain script names it lists into any command line, script name or launched file (production gates match
   process arguments; your old system prompt contained one, which is why you were restarted); never write, re-download or
   delete anything under data/ (new public data lives under your experiment's work/); no transfer > 1 GB from the Mac during
   23:30–01:45Z (EC2 ↔ S3 is fine and free); production's reserves (R6, R7a T1–T4, R7b Calu-3) never enter a working set.
2. The cloud launcher now tags every instance and volume Project=lab and refuses a launch above $50 (cap) or a day above
   $150; size your launches accordingly (several smaller batches are fine).
3. Your experiment needs its own `.gitignore` covering `work/` artefacts (npz, h5ad, parquet, logs > 1 MB); write it.
4. The lab lead commits; you still run no git commands.
5. You were paused during the 23:30–01:45Z window (at most 2 lab agents run then) and resumed after it. Your EC2 work can
   run unattended through such windows in future: launch, then wait with `labcloud.sh wait` in the background.
