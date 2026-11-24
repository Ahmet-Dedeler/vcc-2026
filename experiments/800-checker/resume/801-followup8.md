FROM THE LAB LEAD (2026-10-17 10:18Z by the clock): follow-up 8, preparation for 10-22. Budget ≤ $0.6, EC2 only (Mac transfers < 1 GB).

**Why.** Your kit builds folds with A/B/C's relabelled control frames. On 10-22, D/E/F's own controls exist
(`data/challenge-def/context_*.h5ad`, production's download), and the day path uses them. A D/E/F-panel fold on D/E/F's own
frames is the closer instrument.

**Task.** Add an option to `kit/kit801.sh` (and whatever it calls): `--frames <bundle dir>` (default: today's behaviour). It takes
the control frames (and the basal-expression floor, if your recipe derives one from the frames) from that bundle.
- **KA:** `--frames data/challenge` must reproduce LABCk / LABC exactly (dumps byte-identical to `dump708-LABC-*`).
- Do NOT run it on data/challenge-def (it does not exist yet).
- Note in RESULTS.md §18 which arrays the frames change, so the 10-22 read can state what differs between the two fold
  variants.
Then stop.
