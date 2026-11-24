FROM THE LAB LEAD: an instrument read, not a reopening. 801 has registered **LABC/<held>** for the three cancer held lines:
A/B/C's 300 targets built with the fresh-fold (L) recipe, D/E/F rules, inc = pct+P2g
(`experiments/801-lab-folds/work/registry.json`, status `labc`). The lab is measuring how much the fresh panels over-read
direction changes against an Arc-drawn panel (`experiments/800-checker/panelrank800.py`, `overread800.py`; LAB-LOG 10:41Z).
- Price **sgu1, sgu1all and 707 vs inc** on the three LABC folds with price810.py as frozen (record its sha; no re-tune), through
  labcloud.sh (≤ $0.3, on-demand, Project=lab).
- Write standard-format arrays to `experiments/810-strong-direction/work/labc/bulk/LABC__<held>/` (`inc.npz`, `sgu1.npz`,
  `sgu1all.npz`, `707.npz`). Use the same key names as your price-final arrays, so panelrank800 / overread800 can read them.
- Engine KAs on every fold; a failed or missing KA stops the read.
- Report in RESULTS.md as "LABC instrument read": per arm, point, CI and folds > 0. Add one wear row per arm (N +1 on
  A/B/C-panel reads). Then stop.
