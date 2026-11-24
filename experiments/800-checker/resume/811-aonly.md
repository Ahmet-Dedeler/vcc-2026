FROM THE LAB LEAD (2026-10-16 12:44Z by the clock): an instrument read, not a reopening. 801 has registered 18 **A-only** folds,
`A1/<held>` … `A6/<held>` (`experiments/801-lab-folds/work/registry.json`, status `a-only`). They are arc800 draws with the
essential screens removed, on the L recipe, with inc = pct+P2g.
- Price **c0_l25, c1_l25 and src_l25 vs inc** on all 18 with your frozen code: m811.py plus labc811.py's wrapper pattern if
  the fold id needs it. Record the shas, and the one-line detection diff if any. Run through labcloud.sh (≤ $2, on-demand,
  Project=lab).
- Write standard arrays to `experiments/811-lolo-model/work/run-m811aonly/bulk/A<k>__<held>/`.
- Engine KAs on every fold; reproduce one L-fold known answer as you did for LABC.
- Report in RESULTS.md as "A-only instrument read": per arm, point, CI, folds > 0, and the many-line increment. Add wear rows.
  Then stop.
