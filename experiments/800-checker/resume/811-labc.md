FROM THE LAB LEAD: an instrument read, not a reopening. 801 has registered **LABC/<held>** for the three cancer held lines:
A/B/C's 300 targets built with the fresh-fold (L) recipe, D/E/F rules, inc = pct+P2g
(`experiments/801-lab-folds/work/registry.json`, status `labc`). The lab is measuring how much the fresh panels over-read
direction changes against an Arc-drawn panel (`experiments/800-checker/panelrank800.py`, `overread800.py`; LAB-LOG 10:41Z).
- Price **c0_l25, c1_l25 and src_l25 vs inc** on the three LABC folds with m811.py as frozen (`e00192e6…`; if the fold id is not
  recognised, add only a detection line for "LABC/", record the new sha and the one-line diff, and treat LABC exactly like an L
  fold), through labcloud.sh (≤ $0.5, on-demand, Project=lab).
- Write standard-format arrays to `experiments/811-lolo-model/work/run-m811labc/bulk/LABC__<held>/`.
- Engine KAs on every fold.
- Report in RESULTS.md as "LABC instrument read": per arm, point, CI and folds > 0, plus the many-line increment over
  sources-only. Add wear rows (N +1). Then stop.
